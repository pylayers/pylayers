#-*- coding:Utf-8 -*-
"""

module ezone
============

.. currentmodule:: pylayers.gis.ezone

This class handles the description of an earth zone

An earth zone gathers raster and vector data
raster data comes either from srtm or aster data
vector data comes from openstreetmap

.. autosummary::
    :toctree: generated



"""
import h5py
import numpy as np
import pandas as pd
import struct
import zipfile
import pickle
import os
import pdb
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from osgeo import gdal
from scipy.interpolate import interp2d
#from imposm.parser import OSMParser
#from geomutil import *
#from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.plotutil as plu
from pylayers.util.project import *
from shapely.geometry import Polygon
from pylayers.gis.gisutil import *
import pylayers.gis.srtm as srtm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

def maxloc(f,threshold=-np.sqrt(2)):
    """ determine local maximum above a threshold

    Parameters
    ----------

    f : np.array
    threshold : float

    Returns
    -------

    g :

    """
    f_decr = np.roll(f,1,axis=1)
    f_decl = np.roll(f,-1,axis=1)
    ind = np.where((f>f_decr)&(f>f_decl)&(f>threshold))
    g = threshold*np.ones(np.shape(f))
    #
    g[ind] = f[ind]
    return(g,ind)

def enctile(lon,lat):
    """ encode tile prefix from (lon,lat)

    Parameters
    ----------

    lon : float
        longitude (degrees)
    lat : float
        latitude (degrees)

    Returns
    -------

    prefix : string
        srtm prefix filename

    Examples
    --------

        >>> from pylayers.gis.ezone import *
        >>> assert enctile(-1.5,48.5)=='N48W002'
        >>> assert enctile(0.5,48.5)=='N48E000'

    """
    if lon>0:
        slon='E'
    else:
        slon='W'

    if lat>0:
        slat='N'
    else:
        slat='S'
    # rounding and integer to string  conversion 
    if lat>0:
        clat = str(np.floor(abs(lat)).astype('int'))
    else:
        clat = str(np.ceil(abs(lat)).astype('int'))
    if lon>0:
        clon = str(np.floor(abs(lon)).astype('int'))
    else:
        clon = str(np.ceil(abs(lon)).astype('int'))
    # 2 char formating
    if len(clon)<2:
        clon ='00'+clon
    else:
        clon='0'+clon
    if len(clat)<2:
        clat ='0'+clat

    prefix = slat+clat+slon+clon
    return prefix

def dectile(prefix='N48W002'):
    """ decode tile name

    Parameters
    ----------

    _filehgt : string

    Examples
    --------

    >>> import ezone
    >>> ezone.dectile('N48W002')
    (-2.0, -1.0, 48,49)

    """
    if prefix[0]=='N':
        latmin = int(prefix[1:3])
        latmax = latmin+1
    if prefix[0]=='S':
        latmin = -int(prefix[1:3])
        latmax = latmin+1

    if prefix[3]=='W':
        st = prefix[4:]
        if st[0]=='0':
            lonmin = -int(st[1:])
        else:
            lonmin = -int(st)
        lonmax = lonmin+1

    if prefix[3]=='E':
        st = prefix[4:]
        if st[0]=='0':
            lonmin = int(st[1:])
        else:
            lonmin = int(st)

        lonmax = lonmin+1

    return (lonmin,lonmax,latmin,latmax)

def expand(A):
    """ expand numpy array

    Parameters
    ----------

    A : np.array (MxN)

    """
    M,N = A.shape
    t = np.kron(A.flatten(),np.ones(N))
    u = np.triu(np.ones((N,N))).flatten()
    v = np.kron(np.ones(M),u)
    w  = t *  v
    return(w.reshape(M,N,N).swapaxes(1,2)[:,1:,:])
    #return(w.reshape(M,N,N).swapaxes(1,2))

def conv(extent,m,mode='tocart'):
    """ convert zone to cartesian or lon lat

    Parameters
    ----------

    extent : (lonmin,lonmax,latmin,latmax)
    m  : matplotlib mapping
    mode : string
        'tocart' | 'toll'

    Returns
    -------

    out : np.array
        [xmin,xmax,ymin,ymax] if mode == 'tocart'
        [lonmin,lonmax,latmin,latmax] if mode == 'toll'

    """
    if mode=='tocart':
        pll = m(extent[0],extent[2])
        pur = m(extent[1],extent[3])
        out = np.array([pll[0],pur[0],pll[1],pur[1]])
    if mode=='toll':
        lllon,lllat = m(extent[0],extent[2],inverse=True)
        rulon,rulat = m(extent[1],extent[3],inverse=True)
        out = np.array([lllon,rulon,lllat,rulat])
    return(out)

def zone(pt,rm=1000):
    """ extract a region from a point and a radius

    Parameters
    ----------
    pt : np.array
        lon lat
    rm : float
        radius (meters)

    """

    lonc = pt[0]
    latc = pt[1]
    Rearth = 6371000.
    dphi_rad = rm/Rearth
    lonmin = lonc - dphi_rad*180/np.pi
    lonmax = lonc + dphi_rad*180/np.pi
    latmin = latc - dphi_rad*180/np.pi
    latmax = latc + dphi_rad*180/np.pi
    return (lonmin,latmin,lonmax,latmax)

class DEM(PyLayers):
    """ Class Digital Elevation Model
    """
    def __init__(self,prefix):

        self.prefix = prefix

        (lom,loM,lam,laM) = dectile(self.prefix)
        self.extent = (lom,loM,lam,laM)
        self.lon_0  = (lom+loM)/2.
        self.lat_0  = (lam+laM)/2.


        self.lL0 = np.array([lom,lam])

        self.m = Basemap(llcrnrlon = lom,
                         llcrnrlat = lam,
                         urcrnrlon = loM,
                         urcrnrlat = laM,
                         resolution = 'i',projection='cass',
                         lon_0 = self.lon_0,
                         lat_0 = self.lat_0)

    def dwlsrtm(self):
        """ download srtm tile

        Parameters
        ----------

        lat  : float
        lon  : float

        """
        downloader = srtm.SRTMDownloader()
        downloader.loadFileList()
        #ilat = int(np.floor(lat).astype('int'))
        #ilon = int(np.floor(lon).astype('int'))
        # latitude, longitude
        tile = downloader.getTile(self.lL0[1],self.lL0[0])
        self.hgts = np.array(tile.data).reshape(1201,1201)
        self.hgts[self.hgts<0]=0

    def loadsrtm(self):
        """ load hgt and lcv files from srtm directory

        """

        _filehgt = self.prefix+'.HGT'
        _filelcv = self.prefix+'.lcv'
        filehgt = pyu.getlong(_filehgt,os.path.join('gis','srtm'))
        filelcv = pyu.getlong(_filelcv,os.path.join('gis','srtm'))


        data = np.fromfile(filehgt,dtype='>i2')
        self.hgt = data.reshape(1201,1201)

        data = np.fromfile(filelcv,dtype='>i1')
        self.lcv = data.reshape(1201,1201)

    def loadaster(self,fileaster=[]):
        """ load Aster files

        """


        # construct filename from prefix
        _fileaster = 'ASTGTM2_'+self.prefix+'_dem.tif'

        if fileaster==[]:
            fileaster = pyu.getlong(_fileaster,os.path.join('gis','aster'))
        else:
            _fieleaster = pyu.getshort(fileaster)

        # zip extraction
        ext = _fileaster.split('.')
        if ext[1]=='zip':
            with zipfile.Zipfile(fileaster) as zf:
                for member in zf.infolist():
                    words = member.filename.split('/')
                    path = dest_dir
                    for word in words[:-1]:
                        drive, word = os.path.splitdrive(word)
                        head, word = os.path.split(word)
                        if word in (os.curdir, os.pardir, ''):
                            continue
                        path = os.path.join(path, word)
                    zf.extract(member, path)

        f = gdal.Open(fileaster)
        self.hgta = f.ReadAsArray()

    def show(self,**kwargs):
        """ DEM vizualisation

        Parameters
        ----------

        cmap : colormap

        """
        defaults ={'cmap': plt.cm.jet,
                   'source':'srtm',
                   'alpha':1}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        if 'fig' in kwargs:
            fig = kwargs['fig']
        else:
            fig = plt.figure()

        if 'ax' in kwargs:
            ax = kwargs['ax']
        else:
            ax =  fig.add_subplot(111)

        #im = ax.imshow(dem[ilat[0]:(ilat[-1]+1),ilon[0]:(ilon[-1]+1)],extent=(lonmin,lonmax,latmin,latmax))
        if kwargs['source']=='srtm':
            im = ax.imshow(self.hgts,extent=(self.extent[0],self.extent[1],self.extent[2],self.extent[3]),alpha=kwargs['alpha'],cmap=kwargs['cmap'])
        if kwargs['source']=='aster':
            im = ax.imshow(self.hgta,extent=(self.extent[0],self.extent[1],self.extent[2],self.extent[3]),alpha=kwargs['alpha'],cmap=kwargs['cmap'])

        # handling colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(im,cax)
        cb.set_label('Height (meters)')
        return fig,ax,divider

class Ezone(PyLayers):
    """
        An Ezone is a region of earth delimited by
        (lonmin,lonmax,latmin,latmax)

        An Ezone is stored in hdf5 format

        Attributes
        ----------

        extent : tuple (lonmin,lonmax,latmin,latmax)
            dictionnary of cities
        pll    : point lower left
        pur    : point upper right
        m      : Basemap coordinates converter


    """
    def __init__(self,prefix):
        """
        """
        self.prefix = prefix

        (lom,loM,lam,laM) = dectile(self.prefix)

        self.extent = (lom,loM,lam,laM)
        self.lon_0  = (lom+loM)/2.
        self.lat_0  = (lam+laM)/2.

        self.lL0 = np.array([lom,lam])

        self.m = Basemap(llcrnrlon = lom,
                         llcrnrlat = lam,
                         urcrnrlon = loM,
                         urcrnrlat = laM,
                         resolution = 'i',projection='cass',
                         lon_0 = self.lon_0,
                         lat_0 = self.lat_0)
        self.pll = self.m(self.extent[0],self.extent[2])
        self.pur = self.m(self.extent[1],self.extent[3])
        self.extent_c = (self.pll[0],self.pur[0],self.pll[1],self.pur[1])

    def __repr__(self):
        st = self.prefix+'\n'
        ncar = len(st)
        for c in range(ncar):
            st = st+'-'
        st = st+'\n'
        st = st+'latlon (deg) : '+str(self.extent)+'\n'
        st = st+'cartesian (meters) : ' +("[%2.3f" % self.pll[0])+' '\
                               +("%2.3f" % self.pur[0])+' '\
                               +("%2.3f" % self.pll[1])+' '\
                               +("%2.3f ] " % self.pur[1])+'\n'

        if 'dbldg' in self.__dict__:
            st = st + '\n'
            st = st + 'Buildings \n'
            st = st + '--------- \n'
            st = st + "i-longitude : "+str(self.blom)+\
            ' '+str(self.bloM)+'\n'
            st = st + "i-latitude  : "+str(self.blam)+\
            ' '+str(self.blaM)+'\n'


        return(st)

    def building(self,ltile):
        """ get building in arrays from a list of subtiles

        Parameters
        ----------

        ltile : list of string

        """
        ltile = filter(lambda x : x in self.dbldg,ltile)
        self.lpoly = []
        for it in ltile:
            h,p=self.dbldg[it]
            for pt in p:
                try:
                    poly = np.vstack((poly,np.array(pt)))
                except:
                    poly = np.array(pt)
                if (sum(pt)==0):
                    self.lpoly.append(poly[0:-1,:])
                    del poly
            try:
                th = np.hstack((th,h[:,3]))
            except:
                th = h[:,3]


        self.height = th

    def ls(self):
        files = os.listdir(os.path.join(basename,'gis','h5'))
        for f in files:
            print f

    def getdem(self):
        """ get a digital elevation model

        Parameters
        ----------

        lon : float
        lat : float

        """
        lon = self.lL0[0]
        lat = self.lL0[1]
        # Determine the srtm and aster file name
        dirsrtm  = os.environ['DIRSRTM']
        diraster = os.environ['DIRASTER']

        _filehgt  = self.prefix+'.HGT'
        _filelcv  = self.prefix+'.lcv'
        _fileaster  = 'ASTGTM2_'+self.prefix+'_dem.tif'

        filehgt = os.path.join(dirsrtm,_filehgt)
        filelcv = os.path.join(dirsrtm,_filelcv)
        fileaster = os.path.join(diraster,_fileaster)

        if (os.path.isfile(filehgt) & os.path.isfile(filelcv)):
            print "Load srtm file"
            D = DEM(self.prefix)
            D.loadsrtm()
            self.hgts = D.hgt
            self.lcv  = D.lcv
            self.m    = D.m
            self.extent = D.extent
            self.pll = self.m(self.extent[0],self.extent[2])
            self.pur = self.m(self.extent[1],self.extent[3])
            self.rebase(source='srtm')
        else:
            print "Download srtm file"
            D = DEM(self.prefix)
            D.dwlsrtm()
            self.hgts = D.hgts
        if os.path.isfile(fileaster):
            print "Load aster file"
            D = DEM()
            D.loadaster(fileaster=fileaster)
            self.hgta = D.hgta
        else:
            print "no aster file for this point"



    def loadh5old(self,_fileh5):
        """ load an Ezone from hdf5 file

        Parameters
        ----------

        _fileh5 : string

        Examples
        --------


        """

        self.fileh5 = pyu.getlong(_fileh5,os.path.join('gis','h5'))
        f = h5py.File(self.fileh5,'a')
        self.extent = f['extent'].value

        self.lon0 = (self.extent[0]+self.extent[1])/2.
        self.lat0 = (self.extent[2]+self.extent[3])/2.


        self.m = Basemap(llcrnrlon = self.extent[0],
                         llcrnrlat = self.extent[2],
                         urcrnrlon = self.extent[1],
                         urcrnrlat = self.extent[3],
                         resolution = 'i',projection='cass',
                         lon_0 = self.lon0,
                         lat_0 = self.lat0)
        self.pll = self.m(self.extent[0],self.extent[2])
        self.pur = self.m(self.extent[1],self.extent[3])

        #lcity = f['osm']
        #self.dcity= {}
        #for city in lcity:
        #    self.dcity[city] = {}
        #    self.dcity[city]['bdpt']   = f['osm'][city]['bdpt'].value
        #   self.dcity[city]['bdma']  = f['osm'][city]['bdma'].value
        #    self.dcity[city]['extent'] = f['osm'][city]['extent'].value

        #    lonm =  self.dcity[city]['bdpt'][:,0].min()
        #    lonM =  self.dcity[city]['bdpt'][:,0].max()
        #    latm =  self.dcity[city]['bdpt'][:,1].min()
        #    latM =  self.dcity[city]['bdpt'][:,1].max()
        #    self.dcity[city]['extent'] = (lonm,lonM,latm,latM)

        self.hgts = f['hgt'].value
        self.lcvs = f['lcv'].value
        f.close()

        self.rebase()
        #Nlat,Nlon = self.hgt.shape
        #self.lon = np.linspace(self.extent[0],self.extent[1],Nlon)
        #self.lat = np.linspace(self.extent[3],self.extent[2],Nlat)
        #self.lonstep = (self.extent[1]-self.extent[0])/(Nlon-1)
        #self.latstep = (self.extent[3]-self.extent[2])/(Nlat-1)
        #self.tocart(Nx=Nlon,Ny=Nlat)

    def rebase(self,source):
        """ reevaluate base

        """
        if source =='srtm':
            Nlat,Nlon = self.hgts.shape
        else:
            Nlat,Nlon = self.hgta.shape

        self.lon = np.linspace(self.extent[0],self.extent[1],Nlon)
        self.lat = np.linspace(self.extent[3],self.extent[2],Nlat)
        self.lonstep = (self.extent[1]-self.extent[0])/(Nlon-1.)
        self.latstep = (self.extent[3]-self.extent[2])/(Nlat-1.)
        self.tocart(Nx=Nlon,Ny=Nlat)

    def tocart(self,Nx=1201,Ny=1201,source='srtm'):
        """ convert to cartesian coordinates

        Parameters
        ----------

        Nx : int
            Number of points along x
        Nx : int
            Number of points along y

        """

        # x : longitude axis (axis = 1)
        x = np.linspace(self.pll[0],self.pur[0],Nx)
        # y : latitude axis (axis = 0)
        y = np.linspace(self.pll[1],self.pur[1],Ny)


        # going back in lon,lat
        #      x (axis 1)
        #   ------->
        #  |
        #  | y (axis 0)
        #  |
        lon,lat = self.m(x[None,:],y[:,None],inverse=True)

        # getting the closest integer index
        rx = np.round((lon - self.extent[0]) / self.lonstep).astype(int)
        ry = np.round((lat - self.extent[2]) / self.latstep).astype(int)

        #
        if source=='srtm':
            self.hgts_cart = self.hgts[ry,rx]
        else:
            self.hgta_cart = self.hgta[ry,rx]
        #self.lcv_cart = self.lcv[ry,rx]

        self.x = x
        # axis inversion
        self.y = y[::-1]
        self.extent_c = (self.x.min(),self.x.max(),self.y.min(),self.y.max())

    def profile(self,pa,pb,**kwargs):
        """ profile extraction between 2 points

        Parameters
        ----------

        pa : tuple
            termination point a
        pb : tuple
            termination point b
        Npt : int
            number of points
        ha : float
            antenna height a
        hb :
            antenna height b
        K : float
            K factor
        fGHz : float
            frequency in GHz
        source : string
            'aster' | 'srtm'

        """

        defaults = {'Npt':1000,
                    'ha':30,
                    'hb':1.5,
                    'K':1.3333,
                    'fGHz':.3,
                    'source':'srtm'}

        for key in defaults:
            if key not in kwargs:
                kwargs[key] = defaults[key]

        # wavelength
        lmbda = 0.3/kwargs['fGHz']

        # transmitter coordinates
        x_a,y_a = self.m(pa[0],pa[1])

        # receiver coordinates
        x_b,y_b = self.m(pb[0],pb[1])

        x = np.linspace(x_a,x_b,kwargs['Npt'])
        y = np.linspace(y_a,y_b,kwargs['Npt'])

        d = np.sqrt((x-x[0])*(x-x[0])+(y-y[0])*(y-y[0]))

        dh = d*(d[::-1])/(2*kwargs['K']*6375e3)

        #if mode=='cover':
        #    extent_c = np.array([x.min(),x.max(),y.min(),y.max()])

        lon,lat = self.m(x,y,inverse=True)

        rx = np.round((lon - self.extent[0]) / self.lonstep).astype(int)
        ry = np.round((self.extent[3]-lat) / self.latstep).astype(int)

        # add earth sphericity deviation to hgt (depends on K factor)
        if kwargs['source']=='srtm':
            height = self.hgts[ry,rx] + dh

        if kwargs['source']=='srta':
            height = self.hgta[ry,rx] + dh

        # seek for local maxima along link profile
        m,ind = maxloc(height[None,:])

        ha = height[0] + kwargs['ha']
        hb = height[-1]+ kwargs['hb']
        LOS = ha+(hb-ha)*d/d[-1]
        diff = height-LOS
        fac  = np.sqrt(2*d[-1]/(lmbda*d*d[::-1]))
        nu   = diff*fac
        num,ind  = maxloc(nu[None,:])

        #plt.plot(d,dh,'r',d,height,'b',d,m[0,:],d,LOS,'k')
        #plt.figure()
        #plt.plot(d,nu,num)

        return(height,d,dh,nu,num,m,LOS)

    def cov(self,**kwargs):
        """ coverage around a point

        Parameters
        ----------

        pc : np.array
            center point in cartesian coordinates
        Nphi : int
            Number of angular direction
        Nr : int
            Number of points along radius
        Rmax : float
            Radius maximum (meters)
        Hr : float
            Receiver height
        Ht : float
            Transmitter height
        K : float
            K factor

        """
        defaults = {'pc':(27000,12000),
                    'Nphi':'360',
                    'Nr':200,
                    'Rmax':4000,
                    'Ht':30,
                    'Hr':1.5,
                    'K':1.3333,
                    'fGHz':.3,
                    'divider':[]
                    }

        for key in defaults:
            if key not in kwargs:
                kwargs[key] = defaults[key]

        pc = kwargs['pc']
        lmbda = 0.3/kwargs['fGHz']
        phi  = np.linspace(0,2*np.pi,kwargs['Nphi'])[:,None]
        r  = np.linspace(0.02,kwargs['Rmax'],kwargs['Nr'])[None,:]

        # cartesian
        x  = pc[0] + r*np.cos(phi)
        y  = pc[1] + r*np.sin(phi)
        extent_c = np.array([x.min(),x.max(),y.min(),y.max()])

        # back to lon lat
        lon,lat = self.m(x,y,inverse=True)

        rx = np.round((lon - self.extent[0]) / self.lonstep).astype(int)
        ry = np.round((self.extent[3]-lat) / self.latstep).astype(int)

        # dem
        dem = self.hgts[ry,rx]


        # adding effect of earth equivalent curvature
        R = expand(r)
        B = r.T-R
        h_earth = (R*B)/(2*kwargs['K']*6375e3)

        # ground height + antenna height
        Ha = kwargs['Ht'] + self.hgts[ry[0,0],rx[0,0]]
        Hb = kwargs['Hr'] + dem

        # Nphi x Nr x Nr
        Hb = Hb[:,None,:]
        # LOS line
        LOS  = Ha+(Hb-Ha)*R/r.T
        diff = expand(dem)+h_earth-LOS
        fac  = np.sqrt(2*r[...,None]/(lmbda*R*B))
        nu   = diff*fac
        #num,ind  = maxloc(nu)
        numax = np.max(nu,axis=2)
        w = numax -0.1
        L = 6.9 + 20*np.log10(np.sqrt(w**2+1)-w)
        LFS = 32.4 + 20*np.log10(r)+20*np.log10(kwargs['fGHz'])
        Ltot = LFS+L

        return x,y,r,R,dem,LOS,h_earth,diff,fac,nu,numax,LFS,Ltot


    def cover(self,**kwargs):
        """ coverage around a point

        Parameters
        ----------

        pc : np.array
            center point in cartesian coordinates
        Nphi : int
            Number of angular direction
        Nr : int
            Number of points along radius
        Rmax : float
            Radius maximum (meters)
        Hr : float
            Receiver height
        Ht : float
            Transmitter height
        K : float
            K factor

        """
        defaults = {'pc':(27000,12000),
                    'Nphi':'360',
                    'Nr':200,
                    'Rmax':4000,
                    'Ht':30,
                    'Hr':1.5,
                    'K':1.3333,
                    'fGHz':.3,
                    'divider':[]
                    }

        for key in defaults:
            if key not in kwargs:
                kwargs[key] = defaults[key]

        if 'fig' not in kwargs:
            f,a = plt.subplots(1,1)
        else:
            f = kwargs['fig']
            a = kwargs['ax']

        pc = kwargs['pc']
        lmbda = 0.3/kwargs['fGHz']
        phi  = np.linspace(0,2*np.pi,kwargs['Nphi'])[:,None]
        r  = np.linspace(0.02,kwargs['Rmax'],kwargs['Nr'])[None,:]

        x  = pc[0] + r*np.cos(phi)
        y  = pc[1] + r*np.sin(phi)
        extent_c = np.array([x.min(),x.max(),y.min(),y.max()])

        # Triangulation
        triang = tri.Triangulation(x.flatten(),y.flatten())
        lon,lat = self.m(triang.x,triang.y,inverse=True)
        # back in lon,lat coordinates
        triang.x = lon
        triang.y = lat

        lon,lat = self.m(x,y,inverse=True)

        rx = np.round((lon - self.extent[0]) / self.lonstep).astype(int)
        ry = np.round((self.extent[3]-lat) / self.latstep).astype(int)

        cov = self.hgts[ry,rx]


        # adding effect of earth equivalent curvature
        R = expand(r)
        B = r.T-R
        h_earth = (R*B)/(2*kwargs['K']*6375e3)

        # ground height + antenna height
        Ha = kwargs['Ht'] + self.hgts[ry[0,0],rx[0,0]]
        Hb = kwargs['Hr'] + cov

        pdb.set_trace()
        # Nphi x Nr x Nr
        Hb = Hb[:,None,:]
        # LOS line
        LOS  = Ha+(Hb-Ha)*R/r.T
        diff = expand(cov)+h_earth-LOS
        fac  = np.sqrt(2*r[...,None]/(lmbda*R*B))
        nu   = diff*fac
        num  = maxloc(nu)
        numax = np.max(num,axis=1)
        w = numax -0.1
        L = 6.9 + 20*np.log10(np.sqrt(w**2+1)-w)
        LFS = 32.4 + 20*np.log10(r)+20*np.log10(kwargs['fGHz'])
        Ltot = -(LFS+L)

        # display coverage region
        #plt.tripcolor(triang, cov.flatten(), shading='gouraud', cmap=plt.cm.jet)
        #f,a = self.show(fig=f,ax=a,contour=False,bldg=True,height=False,coord='cartesian',extent=extent_c)
        f,a,d = self.show(fig=f,ax=a,contour=False,bldg=True,height=False,coord='lonlat',extent=self.extent)
        tc = a.tripcolor(triang, Ltot.flatten(), shading='gouraud', cmap=plt.cm.jet,vmax=-50,vmin=-130)
        #tc = a.tripcolor(triang, w.flatten(), shading='gouraud', cmap=plt.cm.jet,vmax=-50,vmin=-130)
        if kwargs['divider']==[]:
            divider = make_axes_locatable(a)
        else:
            divider=kwargs['divider']
        cax = divider.append_axes("left", size="5%", pad=0.5)
        cb = f.colorbar(tc,cax)
        cb.set_label('Loss(dB)')
        plt.axis('equal')

        return x,y,r,cov,LOS,h_earth,diff,fac,num,LFS


    def rennes(self):
        """
        Building are stored in quadTree.
        The centroid of the building is converted into an integer which
        denotes the corresponding quadtree.
        lbdg is the list of the quadtrees

        """
        Rennes = pd.read_hdf("RennesBis.h5","Building")
        # read a dictionnary of polygons
        fd = open('dpoly.pickle','rb')
        dpoly = pickle.load(fd)
        fd.close()

        keys  = np.array(dpoly.keys())
        lpoly = dpoly.values()
        #u = np.argsort(np.array(keys))

        var  = Rennes['ALT_FAITAG'].values-Rennes['ALT_SOL'].values
        pg   = np.array(map(lambda x : np.array(Polygon(x).centroid.xy).T[0],lpoly)).T
        lL0  = np.array([-2,48])
        # ent : encode in integer
        ibd  = np.array(ent(pg,lL0))
        lbdg = np.unique(ibd)

        for kbld in lbdg:
            idx = np.where(ibd==kbld)[0]
            k0 = keys[idx]
            df0 = Rennes.ix[k0]
            self.dbldg[kbld] = [df0]
            #store = pd.HDFStore('Rennesb.h5',mode="a")
            #group = 'i'+str(k)
            #store[group+'/df']=df0
            #store.close()
            lp0 = [dpoly[k] for k in k0]
            #f = h5py.File("Rennesb.h5",'a')
            z2   = np.zeros((2,))[None,:]
            # add zeros as a delimiter between polygons
            lp0z = map(lambda x:np.vstack((x,z2)),lp0)
            alp0 = reduce(lambda x,y:np.vstack((x,y)),lp0z)
            #self.dbldg[kbld].append(alp0)
            self.dbldg[kbld].append(lp0)
            #f[group+'/poly']=alp0
            #f.close()

    def show(self,**kwargs):
        """ show Ezone

        Parameters
        ----------

        title : string
        xlabel : string
        ylabel : string
        height : boolean
            display dem if True
        bldg : boolean
            display building if True
        coord : string
            'lonlat'| 'cartesian'
        source: string
            'srtm' | 'aster'
        extent : [lonmin,lomax,latmin,latmax]


        Returns
        -------

        fig,ax

        Notes
        -----

        If height is False the DEM is not displayed.
        If extent is a void list all the tile is displayed

        """
        defaults = {'title':'',
                    'xlabel':'Longitude',
                    'ylabel':'Latitude',
                    'figsize':(10,10),
                    'height':True,
                    'bldg':False,
                    'clim':(0,200),
                    'coord':'lonlat',
                    'extent':[],
                    'contour':False,
                    'source':'srtm',
                    'alpha':0.5,
                    'facecolor':'black',
                    'cmap':plt.cm.jet
                   }

        divider = []
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        if 'fig' in kwargs:
            fig = kwargs['fig']
        else:
            fig = plt.figure(figsize=kwargs['figsize'])

        if 'ax' in kwargs:
            ax = kwargs['ax']
        else:
            ax =  fig.add_subplot(111)
        # get zone limitation
        # lon,lat or cartesian

        if kwargs['extent']==[]:
            if kwargs['coord']=='lonlat':
                extent = self.extent
            if kwargs['coord']=='cartesian':
                extent_c = self.extent_c
        else:
            if kwargs['coord']=='cartesian':
                extent_c = kwargs['extent']
                extent = conv(extent_c,self.m,mode='toll')
            if kwargs['coord']=='lonlat':
                extent = kwargs['extent']


        # ploting buildings with collection of polygons
        #

        if kwargs['coord'] == 'cartesian':
            kwargs['xlabel'] = 'W-E Distance (meters)'
            kwargs['ylabel'] = 'N-S Distance (meters)'

        if kwargs['title']=='':
            kwargs['title'] = self.prefix
        ax.set_title(kwargs['title'])
        ax.set_xlabel(kwargs['xlabel'])
        ax.set_ylabel(kwargs['ylabel'])

        if (kwargs['height'] | kwargs['contour']):
            if kwargs['source']=='srtm':
                shaphgt = self.hgts.shape
            else:
                shaphgt = self.hgta.shape
            # full original x and y
            #
            if kwargs['coord']=='lonlat':
                if kwargs['source']=='srtm':
                    x = np.linspace(self.extent[0],self.extent[1],1201)
                    # warning the y axis is inversed
                    y = np.linspace(self.extent[3],self.extent[2],1201)
                    hgt = self.hgts
                else:
                    x = np.linspace(self.extent[0],self.extent[1],3601)
                    y = np.linspace(self.extent[3],self.extent[2],3601)
                    hgt = self.hgta

            if kwargs['coord']=='cartesian':
                self.tocart(source=kwargs['coord'])
                if kwargs['source']=='srtm':
                    x = np.linspace(self.extent_c[0],self.extent_c[1],1201)
                    y = np.linspace(self.extent_c[3],self.extent_c[2],1201)
                    hgt = self.hgts_cart
                else:
                    x = np.linspace(self.extent_c[0],extent_c[1],3601)
                    y = np.linspace(self.extent_c[3],extent_c[2],3601)
                    hgt = self.hgta_cart
                extent = extent_c

            # get index corresponding to the selected zone

            ix = np.where((x>=extent[0]) & (x<=extent[1]))[0]
            iy = np.where((y>=extent[2]) & (y<=extent[3]))[0]

            if kwargs['height']:
                im = ax.imshow(hgt[iy[0]:(iy[-1]+1),ix[0]:(ix[-1]+1)],
                               extent=extent,clim=kwargs['clim'],cmap=kwargs['cmap'],alpha=kwargs['alpha'])
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cb = fig.colorbar(im,cax)
                cb.set_label('Height (meters)')

            if kwargs['contour']:
                cnt = ax.contour(hgt[iy[0]:(iy[-1]+1),ix[0]:(ix[-1]+1)],N=10,extent=extent,origin='upper')

        # display buildings
        if kwargs['bldg']:
            # get subtiles corresponding to extent
            if kwargs['coord']=='cartesian':
                extent = conv(extent_c,self.m,mode='toll')
            ltiles = ext2qt(extent,self.lL0)
            # iterating over subtiles
            for ti in ltiles:
                if ti in self.dbldg.keys():
                   info = self.dbldg[ti][0]
                   poly = self.dbldg[ti][1]
                   if kwargs['coord']=='cartesian':
                       tu   = map(lambda x : self.m(x[:,0],x[:,1]),poly)
                       poly = map(lambda x : np.vstack((x[0],x[1])).T,tu)

                   if kwargs['height']:
                       fig,ax = plu.polycol(poly,
                            clim = kwargs['clim'],
                            facecolor=kwargs['facecolor'],
                            fig=fig,ax=ax)
                   else:
                       fig,ax = plu.polycol(poly,info[:,3],
                            clim = kwargs['clim'],
                            fig=fig,
                            ax=ax)


        return(fig,ax,divider)

    def loadtmp(self,_fileh5='RennesFull.h5'):
        """ load an Ezone from hdf5 file

        Parameters
        ----------

        _fileh5 :

        """

        fileh5 = pyu.getlong(_fileh5,os.path.join('gis','h5'))
        f = h5py.File(fileh5,'r',dtype=np.float32)
        self.bdpt = f['osm']['bdpt'].value
        self.bdma = f['osm']['bdma'].value
        lonm =  self.bdpt[:,0].min()
        lonM =  self.bdpt[:,0].max()
        latm =  self.bdpt[:,1].min()
        latM =  self.bdpt[:,1].max()
        self.extentc = (lonm,lonM,latm,latM)
        D = DEM()
        D.loadsrtm()
        self.extent =  (D.lonmin,D.lonmax,D.latmin,D.latmax)
        self.hgt = D.hgt
        self.lcv = D.lcv
        #vertices = np.ma.masked_array(pt, ma)

    def loadh5(self):
        """ load Ezone from an hdf5 file

        Parameters
        ----------

        prefix : string

        Notes
        -----

        Structure of the hdf5. The file has the following groups

        extent
        dem
            srtm
            aster
        bldg
            u'ia-b'
                info
                poly

        """
        _fileh5 = self.prefix+'.h5'
        fileh5 = pyu.getlong(_fileh5,os.path.join('gis','h5'))
        with h5py.File(fileh5) as fh:
            if 'extent' in fh:
                self.extent = fh['extent'][:]
            if 'dem' in fh:
                if 'srtm' in fh['dem']:
                    self.hgts = fh['dem']['srtm']['hgts'][:]
                if 'aster' in fh['dem']:
                    self.hgta = fh['dem']['aster']['hgta'][:]
            if 'bldg' in fh:
                self.dbldg={}
                for k in fh['bldg']:
                    if (('info' in fh['bldg'][k]) and
                        ('poly' in fh['bldg'][k])):
                        a = fh['bldg'][k]['info'][:]
                        b = fh['bldg'][k]['poly'][:]
                        # convert zeros separated array
                        # to a list of arrays
                        lpol = arr2lp(b)
                        self.dbldg[k] = [a,lpol]

                l1 = map(lambda x : x.replace('i',''),self.dbldg.keys())
                llon = map(lambda x: eval(x.split('-')[0]),l1)
                llat = map(lambda x: eval(x.split('-')[1]),l1)

                self.blom = min(llon)
                self.bloM = max(llon)
                self.blam = min(llat)
                self.blaM = max(llat)


    def saveh5(self):
        """ save Ezone in hdf5 format

        Parameters
        ----------

        _fileh5

        """
        _fileh5 = self.prefix+'.h5'
        fileh5 = pyu.getlong(_fileh5,os.path.join('gis','h5'))
        # if file exists open it in append mode
        f = h5py.File(fileh5,'a')

        # create missing groups
        # extent
        # dem
        if u'dem' not in f.keys():
            dem = f.create_group(u'dem')
        else:
            dem = f[u'dem']

        if hasattr(self,'extent'):
            if 'extent' not in f:
                f['extent'] = self.extent
        if u'bldg' not in f.keys():
            bldg = f.create_group(u'bldg')
        else:
            bldg=f[u'bldg']

        if hasattr(self,'hgts'):
            if u'srtm' not in dem:
                srtm  = dem.create_group(u'srtm')
            else:
                srtm = dem[u'srtm']

            if u'hgts' not in srtm:
                srtm.create_dataset('hgts',shape=self.hgts.shape,data=self.hgts)

        if hasattr(self,'lcvs'):
            if 'lcvs' not in srtm:
                srtm.create_dataset('lcvs',shape=self.lcvs.shape,data=self.lcvs)

        if hasattr(self,'hgta'):
            if u'aster' not in dem:
                aster = dem.create_group(u'aster')
            else:
                aster = dem[u'aster']

            if 'hgta' not in aster:
                aster.create_dataset('hgta',shape=self.hgta.shape,data=self.hgta)

        if hasattr(self,'dbldg'):
            # iterating on subtiles
            for k in self.dbldg:
                # if subtile does not exist create it
                if k not in bldg:
                    bldg.create_group(k)
                bldg[k]['info'] = np.array(self.dbldg[k][0])
                bldg[k]['poly'] = self.dbldg[k][1]

        f.close()

    def build(self,_fileosm,_filehgt,_filelcv):
        """ build an Ezone from heterogeneous files
        """
        pass

if __name__== "__main__":
    pass
    #ez = Ezone()
    #ez.loadh5old('N48W002.h5')
    #extent1 = (-1.7,-1.6,48.05,48.15)
    #extent1_cart = conv(extent1,ez.m,mode='tocart')
