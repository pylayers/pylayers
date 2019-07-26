# -*- coding: utf-8 -*-
"""

.. currentmodule:: pylayers.gis.ezone*

.. autosummary::

"""
import h5py
import numpy as np
import pandas as pd
import zipfile
import pickle
import os
import pdb
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.tri as tri
#from osgeo import gdal
# from imposm.parser import OSMParser
# from geomutil import *
# from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
import pylayers.util.plotutil as plu
import pylayers.antprop.loss as loss
from pylayers.util.project import *
from shapely.geometry import Polygon
from pylayers.gis.gisutil import *
import pylayers.gis.kml as gkml
import pylayers.gis.srtm as srtm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
import smopy
import copy

def maxloc(f, threshold=-np.sqrt(2)):
    """ determine local maximum above a threshold

    Parameters
    ----------

    f : np.array
    threshold : float

    Returns
    -------

    g : np.array
        values of local maximum
    ind : np.array
        index of local maximum

    Examples
    --------

        >>> import numpy as np
        >>> t = np.arange(0,6*np.pi)
        >>> f = np.sin(2*t)*cos(3*t)*sin(5*t)
        >>> g,ind = maxloc(f,threshold=0.3)

    """
    f_decr = np.roll(f,1,axis=1)
    f_decl = np.roll(f,-1,axis=1)
    ind = np.where((f>f_decr)&(f>f_decl)&(f>threshold))
    g = threshold*np.ones(np.shape(f))
    #
    g[ind] = f[ind]
    return(g, ind)

def expand(A):
    """ expand numpy array

    Parameters
    ----------

    A : np.array (MxN)

    Returns
    -------

    """

    M, N = A.shape
    t = np.kron(A.flatten(), np.ones(N))
    u = np.triu(np.ones((N, N))).flatten()
    v = np.kron(np.ones(M), u)
    w = t * v
    # return(w.reshape(M, N, N).swapaxes(1, 2)[:, 1:, :])

    return(w.reshape(M, N, N).swapaxes(1,2))

def conv(extent, m, mode='tocart'):
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
    if mode == 'tocart':
        pll = m(extent[0], extent[2])
        pur = m(extent[1], extent[3])
        out = np.array([pll[0], pur[0], pll[1], pur[1]])

    if mode == 'toll':
        lllon, lllat = m(extent[0], extent[2], inverse=True)
        rulon, rulat = m(extent[1], extent[3], inverse=True)
        out = np.array([lllon, rulon, lllat, rulat])
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
        self.hgts = data.reshape(1201,1201)

        data = np.fromfile(filelcv,dtype='>i1')
        self.lcvs = data.reshape(1201,1201)

    def loadaster(self,fileaster=[]):
        """ load Aster files

        The Aster file has the structure

        ASTGTM2_prefix_dem.tif

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

        #
        # Commented while gdal is broken in anaconda
        #f = gdal.Open(fileaster)
        #self.hgta = f.ReadAsArray()

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
    """ Earth zone

        An Ezone is a class related to a region of Earth delimited by
        (lonmin,lonmax,latmin,latmax)

        An Ezone gathers raster and vector data
            + raster data comes either from srtm or aster data
            + vector data comes from openstreetmap

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
        Parameters
        ----------

        prefix : string
            filename without extension

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
            st = st + '--------- \n'
        if 'hgta' in self.__dict__:
            st = st + 'DEM Aster (hgta) :' + str(self.hgta.shape)+'\n'
        if 'hgts' in self.__dict__:
            st = st + 'DEM srtm (hgts) :' + str(self.hgts.shape)

        return(st)

    def building(self,ltile):
        """ get building in arrays from a list of subtiles

        Parameters
        ----------

        ltile : list of strings

        """
        ltile = [ x for x in ltile if x in self.dbldg ]
        self.lpoly = []
        for it in ltile:
            h, p = self.dbldg[it]
            for pt in p:
                try:
                    poly = np.vstack((poly,np.array(pt)))
                except:
                    poly = np.array(pt)
                if (sum(pt) == 0):
                    self.lpoly.append(poly[0:-1, :])
                    del poly
            try:
                th = np.hstack((th, h[:, 3]))
            except:
                th = h[:, 3]


        self.height = th

    def ls(self):
        files = os.listdir(os.path.join(basename,'gis','h5'))
        for f in files:
            print(f)

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
            print("Load srtm file")
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
            print("Download srtm file")
            D = DEM(self.prefix)
            D.dwlsrtm()
            self.hgts = D.hgts
        if os.path.isfile(fileaster):
            print("Load aster file")
            D = DEM()
            D.loadaster(fileaster=fileaster)
            self.hgta = D.hgta
        else:
            print("no aster file for this point")



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

    def rebase(self,source='srtm'):
        """ reevaluate base

        Parameters
        ----------

        source : string
            source of data 'srtm' or 'aster'

        Notes
        -----

        This methods recalculate the longitude and latitude base based on the
        DEM source choice aster or srtm

        """
        if source =='srtm':
            Nlat,Nlon = self.hgts.shape
        elif source=='aster':
            Nlat,Nlon = self.hgta.shape

        self.lon = np.linspace(self.extent[0],self.extent[1],Nlon)
        self.lat = np.linspace(self.extent[3],self.extent[2],Nlat)
        self.lonstep = (self.extent[1]-self.extent[0])/(Nlon-1.)
        self.latstep = (self.extent[3]-self.extent[2])/(Nlat-1.)
        self.tocart(Nx=Nlon,Ny=Nlat,source=source)

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
        elif source =='aster':
            self.hgta_cart = self.hgta[ry,rx]
        else:
            print('unrecognized source')
        #self.lcv_cart = self.lcv[ry,rx]

        self.x = x
        # axis inversion
        self.y = y[::-1]
        self.extent_c = (self.x.min(), self.x.max(), self.y.min(), self.y.max())

    def profile(self, pa, pb, **kwargs):
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

        Returns
        -------

        height : np.array (,Npt)
            total heigh including eath curvature
        d : np.array(,Npt)
            horizontal distance along the link
        dh : np.array(,Npy)
            earth curvature depending on K factor
        nu : np.array(,Npt)
            Fresnel parameter
        num : np.array(1,Npt)
            Fresnel parameter above threshold (default : -sqrt(2))
        hlos : (,Npt)
            height of the line of sight
        ellFresnel : (2,N)
            Fresnel ellipsoid set of points

        """

        defaults = {'Npt': 1000,
                    'ha': 30,
                    'hb': 1.5,
                    'K': 1.3333,
                    'fGHz': .868,
                    'threshold': -np.sqrt(2),
                    'source': 'srtm'}

        for key in defaults:
            if key not in kwargs:
                kwargs[key] = defaults[key]

        # wavelength
        lmbda = 0.3/kwargs['fGHz']

        # transmitter cartesian coordinates
        x_a, y_a = self.m(pa[0], pa[1])

        # receiver cartesian coordinates
        x_b, y_b = self.m(pb[0], pb[1])

        x = np.linspace(x_a, x_b, kwargs['Npt'])
        y = np.linspace(y_a, y_b, kwargs['Npt'])

        # distance along path
        d = np.sqrt((x-x[0])*(x-x[0])+(y-y[0])*(y-y[0]))

        # equivalent earth curvature
        dh = d*(d[::-1])/(2*kwargs['K']*6375e3)

        # if mode=='cover':
        #    extent_c = np.array([x.min(),x.max(),y.min(),y.max()])

        lon, lat = self.m(x, y, inverse=True)

        Dx = (lon - self.extent[0]) / self.lonstep
        Dy = (self.extent[3]-lat) / self.latstep
        rx = np.floor(Dx).astype(int)
        ry = np.floor(Dy).astype(int)
        dx = Dx-rx
        dy = Dy-ry

        rx_old = np.round((lon - self.extent[0]) / self.lonstep).astype(int)
        ry_old = np.round((self.extent[3]-lat) / self.latstep).astype(int)

        # add earth sphericity deviation to hgt (depends on K factor)
        if kwargs['source'] == 'srtm':
            hll = self.hgts[ry, rx]
            hlr = self.hgts[ry, rx+1]
            hul = self.hgts[ry+1, rx]
            hur = self.hgts[ry+1, rx+1]
            wll = np.sqrt(dx**2+(1-dy)**2)
            wlr = np.sqrt((1-dx)**2+(1-dy)**2)
            wul = np.sqrt(dx**2+dy**2)
            wur = np.sqrt((1-dx)**2+dy**2)
            height = (wll*hll+wlr*hlr+wul*hul+wur*hur)/(wll+wlr+wul+wur)
            height_old = self.hgts[ry_old, rx_old] + dh

        if kwargs['source'] == 'aster':
            hll = self.hgta[ry, rx]
            hlr = self.hgta[ry, rx+1]
            hul = self.hgta[ry+1, rx]
            hur = self.hgta[ry+1, rx+1]
            wll = np.sqrt(dx**2+(1-dy)**2)
            wlr = np.sqrt((1-dx)**2+(1-dy)**2)
            wul = np.sqrt(dx**2+dy**2)
            wur = np.sqrt((1-dx)**2+dy**2)
            height = (wll*hll+wlr*hlr+wul*hul+wur*hur)/(wll+wlr+wul+wur)
            height_old = self.hgta[ry,rx] + dh

        # seek for local maxima along link profile

        m, ind = maxloc(height[None, :])

        ha = height[0] + kwargs['ha']
        hb = height[-1] + kwargs['hb']
        hlos = ha+(hb-ha)*d/d[-1]
        diff = height-hlos
        fac = np.sqrt(2*d[-1]/(lmbda*d*d[::-1]))
        nu = diff*fac
        imax = np.argmax(nu)
        numax = nu[imax]
        #z0 = np.zeros(np.shape(nu))
        #u1 = np.ones(np.shape(nu))
        #z0[imax:]=1
        #d1 = d - z0*d[imax]
        #h1 = hlos - u1*hlos[imax]
        #num, ind = maxloc(nu[None,:],threshold=kwargs['threshold'])
        # construction of first Fresnel ellipsoid
        pa = np.array([0, ha])
        pb = np.array([d[-1], hb])
        w = numax-0.1
        L = 6.9+20*np.log10(np.sqrt(w**2+1)+w)
        LFS = 32.4 + 20*np.log10(d[-1])+20*np.log10(kwargs['fGHz'])
        Ltot = -(LFS+L)
        ellFresnel = geu.ellipse2D(pa, pb, lmbda/2, 100 ,unit = 'meter')

        data ={}
        data= {'height':height,
               'height_old':height_old,
             'd':d,
             'dh':dh,
             'nu':nu,
             'numax':numax,
             'm':m,
             'hlos':hlos,
             'ellFresnel':ellFresnel[0],
             'LFS':LFS,
             'L':L }

        return data


    def route(self, pa, pb, **kwargs):
        """ coverage on a route

        Parameters
        ----------

        pa : np.array (1x2)
            lon,lat
        pb : np.array (Nx2)
            lon,lat
        Nr  :
        Ha  : float
            antenna a height
        Hb  : float
            antenna b height
        fGHz : float
            frequency
        K   : float
            K facteur (default 4/3)
        method :
        source :
        binterb : boolean
            interpolation (default True)

        Returns
        -------

        L : Loss
        lon
        lat

        """

        self.pa = pa
        self.pb = pb
        Nr = kwargs.pop('Nr',200)
        Ha = kwargs.pop('Ha',30)
        Hb = kwargs.pop('Hb',200)
        fGHz = kwargs.pop('fGHz',0.868)
        K = kwargs.pop('K',1.333)
        method = kwargs.pop('method','deygout')
        source = kwargs.pop('source','srtm')
        binterp = kwargs.pop('binterp',True)
        divider = kwargs.pop('divider',[])

        x_a, y_a = self.m(pa[0], pa[1])
        x_b, y_b = self.m(pb[:, 0], pb[:, 1])
        u = np.linspace(0, 1, Nr)
        x = x_a + (x_b[:,None] - x_a)*u[None,:]
        y = y_a + (y_b[:,None] - y_a)*u[None,:]

        lon, lat = self.m(x, y, inverse=True)
        # distance along path
        Dx = x - x[:, 0][:, None]
        Dy = y - y[:, 0][:, None]
        d = np.sqrt(Dx*Dx+Dy*Dy)
        # equivalent earth curvature
        dh = d*(d[:,::-1])/(2*K*6375e3)

        Dlon = (lon - self.extent[0]) / self.lonstep
        Dlat= (self.extent[3]-lat) / self.latstep
        #
        # Interpolation
        #

        if binterp:
            rlon = np.floor(Dlon).astype(int)
            rlat = np.floor(Dlat).astype(int)
            dlon = Dlon - rlon
            dlat = Dlat - rlat

            if source == 'srtm':
                hll = self.hgts[rlat, rlon]
                hlr = self.hgts[rlat, rlon+1]
                hul = self.hgts[rlat+1, rlon]
                hur = self.hgts[rlat+1, rlon+1]

            if source == 'aster':
                hll = self.hgta[rlat, rlon]
                hlr = self.hgta[rlat, rlon+1]
                hul = self.hgta[rlat+1, rlon]
                hur = self.hgta[rlat+1, rlon+1]
            #    height = self.hgta[ry, rx]

            wll = dlon**2 + (1-dlat)**2
            wlr = (1-dlon)**2 + (1-dlat)**2
            wul = dlon**2 + dlat**2
            wur = (1-dlon)**2 + dlat**2
            height = (wll*hll+wlr*hlr+wul*hul+wur*hur)/(wll+wlr+wul+wur) + dh
        else:
            rlon = np.round(Dlon).astype(int)
            rlat = np.round(Dlat).astype(int)
            if source == 'srtm':
                height = self.hgts[rlat, rlon] + dh
            if source == 'srtm':
                height = self.hgta[rlat, rlon] + dh

        self.L = loss.route(x, y, height, Ha, Hb, fGHz, K, method=method)

        return(self.L,lon,lat)

    def cover(self, **kwargs):
        """ coverage around a point

        Parameters
        ----------

        pc : np.array
            center point in lonlat coordinates
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

        Returns
        -------

        triang, L

        """
        defaults = {'pc': (-1.627449, 48.124648),
                    'Nphi': 360,
                    'Nr': 200,
                    'Rmax': 4000,
                    'Ht': 30,
                    'Hr': 1.5,
                    'K': 1.3333,
                    'fGHz': .3,
                    'source': 'srtm',
                    'divider': []
                    }

        for key in defaults:
            if key not in kwargs:
                kwargs[key] = defaults[key]

        if 'fig' not in kwargs:
            f, a = plt.subplots(1, 1)
        else:
            f = kwargs['fig']
            a = kwargs['ax']

        Ht = kwargs['Ht']
        Hr = kwargs['Hr']
        fGHz = kwargs['fGHz']
        K = kwargs['K']
        Rmax = kwargs['Rmax']
        Nr = kwargs['Nr']
        Nphi = kwargs['Nphi']

        x_c, y_c = self.m(kwargs['pc'][0], kwargs['pc'][1])
        pc = np.array([x_c, y_c])

        #lmbda = 0.3/fGHz
        # phi:  Nphi x 1
        phi = np.linspace(0, 2*np.pi, Nphi)[:, None]
        # r :   1 x Nr
        r = np.linspace(0.02, Rmax, Nr)[None, :]

        x = pc[0] + r*np.cos(phi)
        y = pc[1] + r*np.sin(phi)
        # extent_c = np.array([x.min(),x.max(),y.min(),y.max()])

        #dh = d*(d[::-1])/(2*K*6375e3)

        # Triangulation of the coverage zone
        triang = tri.Triangulation(x.flatten(), y.flatten())
        lon, lat = self.m(triang.x, triang.y, inverse=True)

        # back in lon,lat coordinates
        triang.x = lon
        triang.y = lat

        lon, lat = self.m(x, y, inverse=True)

        rx = np.floor((lon - self.extent[0]) / self.lonstep).astype(int)
        if rx.max() >self.hgts.shape[1]:
            rx = np.floor(( self.extent[1]-lon) / self.lonstep).astype(int)

        ry = np.floor((self.extent[3]-lat) / self.latstep).astype(int)
        if ry.max() >self.hgts.shape[0]:
            ry = np.floor((lat - self.extent[2]) / self.latstep).astype(int)
        # height
        #cov = self.hgts[ry, rx]
        if kwargs['source'] == 'srtm':
            height = self.hgts[ry, rx]

        if kwargs['source'] == 'aster':
            height = self.hgta[ry, rx]

        L = loss.cover(x, y, height, Ht, Hr, fGHz, K, method='deygout')
        self.triang = triang
        self.coverage  = L
        return triang, L

#        # adding effect of earth equivalent curvature
#        # r : 1 x 200
#        # R : 1 x 199 x 200
#        R = expand(r)
#        B = r.T-R
#        h_earth = (R*B)/(2*kwargs['K']*6375e3)
#
#        # ground height + antenna height
#        Ha = kwargs['Ht'] + self.hgts[ry[0, 0], rx[0, 0]]
#        Hb = kwargs['Hr'] + cov
#
#        # Nphi x Nr x Nr
#        Hb = Hb[:, None, :]
#        # LOS line
#        LOS = Ha+(Hb-Ha)*R/r.T
#        diff = expand(cov) + h_earth-LOS
#        fac = np.sqrt(2*r[...,None]/(lmbda*R*B))
#        nu = diff*fac
#        num, ind = maxloc(nu)
#        # numax : Nph x Nr
#        numax = np.max(num, axis=1)
#        w = numax -0.1
#        # L : Nph  x Nr
#        L = 6.9 + 20*np.log10(np.sqrt(w**2+1)-w)
#        # LFS : 1 x Nr
#        LFS = 32.4 + 20*np.log10(r) + 20*np.log10(kwargs['fGHz'])
#        # LFS : Nphi x Nr
#        Ltot = -(LFS+L)
#
#        return triang,LFS,L,Ltot

    def to_kmz(self, **kwargs):
        """ export to kmz file
        """
        llcrnrlon = self.extent[0]
        llcrnrlat = self.extent[2]

        urcrnrlon = self.extent[1]
        urcrnrlat = self.extent[3]

        pngsrtm = self.prefix+'_srtm.png'
        pngaster = self.prefix+'_aster.png'
        pngroute = self.prefix+'_route.png'
        pngcover = self.prefix+'_cover.png'
        kmzsrtm = self.prefix+'_srtm.kmz'
        kmzaster = self.prefix+'_aster.kmz'
        kmzroute = self.prefix+'_route.kmz'
        kmzcover = self.prefix+'_cover.kmz'

        #
        # srtm overlay
        #
        fig, ax = gkml.gearth_fig(self.extent,self.extent_c)

        #cs = ax.pcolormesh(self.hgts, cmap='jet')
        csrtm = ax.imshow(self.hgts,extent=self.extent,cmap='jet')


        fig.savefig(pngsrtm, transparent=True, format='png')
        #
        # aster overlay
        #
        fig, ax = gkml.gearth_fig(self.extent,self.extent_c)
        caster = ax.imshow(self.hgta,extent=self.extent,cmap='jet')

        fig.savefig(pngaster, transparent=True, format='png')
        #
        # route overlay
        #
        fig, ax = gkml.gearth_fig(self.extent,self.extent_c)
        sp = ax.scatter(self.pb[:,0], self.pb[:,1], c = -self.L,
                        s=30,linewidth=0,cmap='jet',vmax=-60, vmin=-120)
        fig.savefig(pngroute, transparent=True, format='png')
        #
        # cover overlay
        #
        fig, ax = gkml.gearth_fig(self.extent,self.extent_c)
        tc = ax.tripcolor(self.triang,
                          -self.coverage.flatten(),
                          shading='gouraud',
                          cmap='jet',
                          vmax=-60,
                          vmin=-120,
                          alpha = 1,
                          edgecolors='k',
                          linewidth=0.0)
        fig.savefig(pngcover, transparent=True, format='png')

        gkml.make_kml(self.extent,
                 figs = [pngroute],
                 kmzfile = kmzroute,
                 name = 'route')

        gkml.make_kml(self.extent,
                 figs = [pngcover],
                 kmzfile = kmzcover,
                 name = 'coverage')

        gkml.make_kml(self.extent,
                 figs = [pngsrtm],
                 kmzfile = kmzsrtm,
                 name = 'SRTM DSM')

        gkml.make_kml(self.extent,
                 figs = [pngaster],
                 kmzfile = kmzaster,
                 name = 'ASTER DSM')


    def showcov(self, triang, val,
                vmin=-130,
                vmax=-50,
                cmap=plt.cm.jet,
                bbuild = False,
                btile=True):
        """ Show a coverage

        Parameters
        ----------

        triang : triangulation
        val : values

        """
        lonmin = np.min(triang.x)
        lonmax = np.max(triang.x)
        latmin = np.min(triang.y)
        latmax = np.max(triang.y)
        extent = (lonmin, lonmax, latmin, latmax)
        print(extent)
        mp = smopy.Map((extent[2]+0.1, extent[0]+0.1, extent[3]-0.1,extent[1]-0.1), z=10)
        if bbuild:
            f, ax, d = self.show(fig=f,
                                 ax=ax,
                                 contour=False,
                                 btile=False,
                                 bldg=True,
                                 height=False,
                                 coord='lonlat',
                                 extent=self.extent)

        #ax = plt.gca()
        triang_ = copy.deepcopy(triang)
        if mp!=[]:
            triang_.x,triang_.y = mp.to_pixels(triang_.y, triang_.x)
        #ax = mp.show_mpl(figsize=(10,10))
        ax = plt.gca()
        tc = ax.tripcolor(triang_,
                          val.flatten(),
                          #shading='gouraud',
                          #shading='flat',
                          cmap=cmap,
                          vmax=vmax,
                          vmin=vmin,
                          alpha = 0.4,
                          edgecolors='k',
                          linewidth=0.0)
        #plt.axis('equal')
        ax = mp.show_mpl(ax=ax)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right",size="5%",pad="5%")
        cb = colorbar(tc, cax=cax)
        cb.set_label_text('Loss(dB)',fontsize=18)
        return(cb)

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
        #pg   = np.array(map(lambda x : np.array(Polygon(x).centroid.xy).T[0],lpoly)).T
        pg = np.array([ np.array(Polygon(x).centroid.xy).T[0] for x in lpoly]).T
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
            lp0 = [ dpoly[k] for k in k0 ]
            #f = h5py.File("Rennesb.h5",'a')
            z2   = np.zeros((2,))[None,:]
            # add zeros as a delimiter between polygons
            #lp0z = map(lambda x:np.vstack((x,z2)),lp0)
            #lp0z = [ np.vstack((x,z2)) for x in lp0 ]
            #for y in lp0z:
            #    x = np.vstack((x,y))
            #alp0 = x
            #alp0 = reduce(lambda x,y:np.vstack((x,y)),lp0z)
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
        btile : boolean
            display tile overlay with smopy
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
                    'btile':False,
                    'clim':(0,200),
                    'coord':'lonlat',
                    'zoom': 10,
                    'extent':[],
                    'contour':False,
                    'source':'srtm',
                    'alpha':0.5,
                    'facecolor':'black',
                    'cmap': 'gist_earth'
                   }

        divider = []
        mp = []
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        if 'fig' in kwargs:
            fig = kwargs['fig']
        else:
            fig = plt.figure(figsize=kwargs['figsize'])

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
        if 'ax' in kwargs:
            ax = kwargs['ax']
        else:
            if kwargs['btile']:
            # lon min lon max lat min lat max
            # lat min , lon min , lat max , lon max
            # 2 0 3 1
                mp = smopy.Map((extent[2], extent[0], extent[3], extent[1]), z=kwargs['zoom'])
                ax = mp.show_mpl(figsize=kwargs['figsize'],alpha=0.3)
            else:
                ax = fig.add_subplot(111)




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

            ix = np.where((x >= extent[0]) & (x <= extent[1]))[0]
            iy = np.where((y >= extent[2]) & (y <= extent[3]))[0]

            if kwargs['height']:
                im = ax.imshow(hgt[iy[0]:(iy[-1]+1),ix[0]:(ix[-1]+1)],
                               extent=extent,
                               clim = kwargs['clim'],
                               cmap = kwargs['cmap'],
                               alpha = kwargs['alpha'])
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cb = fig.colorbar(im,cax)
                cb.set_label('Height (meters)')

            if kwargs['contour']:
                cnt = ax.contour(hgt[iy[0]:(iy[-1]+1),ix[0]:(ix[-1]+1)],N=10,extent=extent,origin='upper')

        # display buildings
        # ploting buildings with collection of polygons
        #
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
                       fig, ax = plu.polycol(poly, info[:, 3],
                                             clim=kwargs['clim'],
                                             fig=fig,
                                             ax=ax)

        return fig, ax, divider, mp

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
        lonm = self.bdpt[:,0].min()
        lonM = self.bdpt[:,0].max()
        latm = self.bdpt[:,1].min()
        latM = self.bdpt[:,1].max()
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

                l1 =   [ x.replace('i','') for x in self.dbldg.keys() ]
                llon = [ eval(x.split('-')[0]) for x in l1 ]
                llat = [ eval(x.split('-')[1]) for x in l1 ]

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
