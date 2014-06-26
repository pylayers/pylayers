#-*- coding:Utf-8 -*-
"""

.. currentmodule:: pylayers.gis.ezone

This class handles the description of an earth zone

.. autosummary::
    :toctree: generated

Class Ezone
===========


"""
import h5py
import numpy as np
import pandas as pd
import struct
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
from pylayers.util.project import *
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

def maxloc(f,threshold=-0.7):
    """ determine local maximum

    Parameters
    ----------

    f : np.array

    Returns
    -------

    g :

    """
    f_decr = np.roll(f,1,axis=1)
    f_decl = np.roll(f,-1,axis=1)
    ind = np.where((f>f_decr)&(f>f_decl)&(f>threshold))
    g = np.zeros(np.shape(f))
    #
    g[ind] = f[ind]
    return(g)

def decsrtm(_filehgt='N48W002.HGT'):
    """ decode srtm file

    Parameters
    ----------

    _filehgt : string

    Examples
    --------

    >>> import ezone
    >>> ezone.dechgt('N48W002.HGT')
    (-2.0, -1.0, 48,49)

    """
    if _filehgt[0]=='N':
        latmin = eval(_filehgt[1:3])
        latmax = latmin+1

    if _filehgt[3]=='W':
        lonmin = -eval(_filehgt[5:7])
        lonmax = lonmin+1

    return (lonmin,lonmax,latmin,latmax)

def expand(A):
    M,N = A.shape
    t = np.kron(A.flatten(),np.ones(N))
    u = np.triu(np.ones((N,N))).flatten()
    v = np.kron(np.ones(M),u)
    w  = t *  v
    return(w.reshape(M,N,N).swapaxes(1,2))

def conv(extent,m,mode='tocart'):
    """ convet zone to cartesian or lon lat
    Parameters
    ----------
    extent
    m
    mode : string
        'tocart' | 'toll'

    Returns
    -------
    out : np.array
        [xmin,xmax,ymin,ymax] if mode == 'tocart'
        [lonin,lonmax,latmin,latmax] if mode == 'toll'

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



    """ Class Digital Elevation Model
    """
    def __init__(self):
        pass

    def loadsrtm(self,_filehgt='N48W002.HGT'):
        """
            load hgt and lcv files from srtm directory

        """

        (lom,loM,lam,laM) = decsrtm(_filehgt)

        self.lonmin = lom
        self.lonmax = loM
        self.latmin = lam
        self.latmax = laM
        self.lon_0  = (lom+loM)/2.
        self.lat_0  = (lam+laM)/2.

        _filelcv = _filehgt.replace('.HGT','.lcv')

        filehgt = pyu.getlong(_filehgt,'gis/hgt')
        filelcv = pyu.getlong(_filelcv,'gis/hgt')

        data = np.fromfile(filehgt,dtype='>i2')
        self.hgt = data.reshape(1201,1201)

        data = np.fromfile(filelcv,dtype='>i1')
        self.lcv = data.reshape(1201,1201)

        self.m = Basemap(llcrnrlon = self.lonmin,
                         llcrnrlat = self.latmin,
                         urcrnrlon = self.lonmax,
                         urcrnrlat = self.latmax,
                         resolution = 'i',projection='cass',
                         lon_0 = self.lon_0,
                         lat_0 = self.lat_0)

    def loadaster(self,_fileaster='ASTGTM2_N48W002_dem.tif'):
        """ Load Aster files

        Parameters
        ----------

        """
        fileaster = pyu.getlong(_fileaster,'gis/aster')
        f = gdal.Open(fileaster)
        self.hgt = f.ReadAsArray()
        self.lonmin = -2
        self.lonmax = -1
        self.latmin = 48
        self.latmax = 49



    def show(self,**kwargs):
        """ Ezone vizualisation

        Parameters
        ----------

        cmap : colormap

        """
        defaults ={'cmap': plt.cm.jet}

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
        im = ax.imshow(self.hgt,extent=(self.lonmin,self.lonmax,self.latmin,self.latmax))

        # handling colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(im,cax)
        cb.set_label('Height (meters)')

class Ezone(PyLayers):
    """
        An Ezone is a region of earth delimited by
        (lonmin,lonmax,latmin,latmax)

        An Ezone is stored in hdf5 format

        Attributes
        ----------

        extent : tuple (lonmin,lonmax,latmin,latmax)
        dcity  : dict
            dictionnary of cities
        pll    : point lower left
        pur    : point upper right
        m      : basemap coordinates change


    """
    def __init__(self):
        """
        """
        pass

    def __repr__(self):
        st = self.fileh5+'\n'
        ncar = len(st)
        for c in range(ncar):
            st = st+'-'
        st = st+'\n'
        st = st+str(self.extent)+'\n'
        st = st+'[ '+("%2.3f" % self.pll[0])+' '+("%2.3f" % self.pur[0])+' '
        st = st+("%2.3f" % self.pll[1])+' '+("%2.3f" % self.pur[1])+' ]\n'
        #st = st+str(self.pur)+'\n'
        for city in self.dcity:
            st = st+city+'\n'
            st = st+'  '+str(self.dcity[city]['extent'])+'\n'
            xll,yll = self.m(self.dcity[city]['extent'][0],self.dcity[city]['extent'][2])
            xru,yru = self.m(self.dcity[city]['extent'][1],self.dcity[city]['extent'][3])
            st = st+'  [ '+("%2.3f" % xll)+' '+("%2.3f" % yll)+' '
            st = st+("%2.3f" % xru)+' '+("%2.3f" % yru)+' ]\n'


        return(st)

    def load(self,_fileh5):
        """ load an Ezone from hdf5 file

        Parameters
        ----------

        _fileh5 : string

        Examples
        --------


        """

        self.fileh5 = pyu.getlong(_fileh5,'gis/h5')
        f = h5py.File(self.fileh5,'r')
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

        lcity = f['osm'].keys()
        self.dcity= {}
        for city in lcity:
            self.dcity[city] = {}
            self.dcity[city]['bdpt'] = f['osm'][city]['bdpt'].value
            self.dcity[city]['bdma'] = f['osm'][city]['bdma'].value
            self.dcity[city]['extent'] = f['osm'][city]['extent'].value

        #    lonm =  self.dcity[city]['bdpt'][:,0].min()
        #    lonM =  self.dcity[city]['bdpt'][:,0].max()
        #    latm =  self.dcity[city]['bdpt'][:,1].min()
        #    latM =  self.dcity[city]['bdpt'][:,1].max()
        #    self.dcity[city]['extent'] = (lonm,lonM,latm,latM)

        self.hgt = f['hgt'].value
        self.lcv = f['lcv'].value
        f.close()

        self.rebase()
        #Nlat,Nlon = self.hgt.shape
        #self.lon = np.linspace(self.extent[0],self.extent[1],Nlon)
        #self.lat = np.linspace(self.extent[3],self.extent[2],Nlat)
        #self.lonstep = (self.extent[1]-self.extent[0])/(Nlon-1)
        #self.latstep = (self.extent[3]-self.extent[2])/(Nlat-1)
        #self.tocart(Nx=Nlon,Ny=Nlat)

    def rebase(self):
        """ reevaluate base

        """
        Nlat,Nlon = self.hgt.shape
        self.lon = np.linspace(self.extent[0],self.extent[1],Nlon)
        self.lat = np.linspace(self.extent[3],self.extent[2],Nlat)
        self.lonstep = (self.extent[1]-self.extent[0])/(Nlon-1)
        self.latstep = (self.extent[3]-self.extent[2])/(Nlat-1)
        self.tocart(Nx=Nlon,Ny=Nlat)

    def tocart(self,Nx=1201,Ny=1201):
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
        lon,lat = self.m(x[np.newaxis,:],y[:,np.newaxis],inverse=True)

        # getting the closest integer index
        rx = np.round((lon - self.extent[0]) / self.lonstep).astype(int)
        ry = np.round((lat - self.extent[2]) / self.latstep).astype(int)

        #
        self.hgt_cart = self.hgt[ry,rx]
        #self.lcv_cart = self.lcv[ry,rx]

        self.x = x
        self.y = y[::-1]

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

        """

        defaults = {'Npt':1000,
                    'ha':30,
                    'hb':1.5,
                    'K':1.3333,
                    'fGHz':.3}

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
        height = self.hgt[ry,rx] + dh

        # seek for local maxima along link profile
        m = maxloc(height[np.newaxis,:])

        ha = height[0] + kwargs['ha']
        hb = height[-1]+ kwargs['hb']
        LOS = ha+(hb-ha)*d/d[-1]
        diff = height-LOS
        fac  = np.sqrt(2*d[-1]/(lmbda*d*d[::-1]))
        nu   = diff*fac
        num  = maxloc(nu[np.newaxis,:])

        #plt.plot(d,dh,'r',d,height,'b',d,m[0,:],d,LOS,'k')
        #plt.figure()
        #plt.plot(d,nu,num)

        return(height,d,dh,nu,num,m,LOS)

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
                    }

        for key in defaults:
            if key not in kwargs:
                kwargs[key] = defaults[key]

        pc = kwargs['pc']
        lmbda = 0.3/kwargs['fGHz']
        phi  = np.linspace(0,2*np.pi,kwargs['Nphi'])[:,np.newaxis]
        r  = np.linspace(0.02,kwargs['Rmax'],kwargs['Nr'])[np.newaxis,:]

        x  = pc[0] + r*np.cos(phi)
        y  = pc[1] + r*np.sin(phi)
        extent_c = np.array([x.min(),x.max(),y.min(),y.max()])
        triang = tri.Triangulation(x.flatten(),y.flatten())

        lon,lat = self.m(x,y,inverse=True)
        rx = np.round((lon - self.extent[0]) / self.lonstep).astype(int)
        ry = np.round((self.extent[3]-lat) / self.latstep).astype(int)
        cov = self.hgt[ry,rx]


        # adding effect of earth equivalent curvature
        R = expand(r)
        B = r.T-R
        hearth = (R*B)/(2*kwargs['K']*6375e3)

        # ground height + antenna height
        Ha = kwargs['Ht'] + self.hgt[ry[0,0],rx[0,0]]
        Hb = kwargs['Hr'] + cov

        Hb = Hb[:,np.newaxis,:]
        # LOS line
        LOS  = Ha+(Hb-Ha)*R/r.T
        diff = expand(cov)+hearth-LOS
        fac  = np.sqrt(2*r[...,np.newaxis]/(lmbda*R*B))
        nu   = diff*fac
        num  = maxloc(nu)
        numax = np.max(num,axis=1)

        w = numax -1
        L = 6.9 + 20*np.log10(np.sqrt(w**2+1)-w)
        LFS = 32.4 + 20*np.log10(r/1000)+20*np.log10(kwargs['fGHz']*1000)
        Ltot = -(LFS+L)

        # display coverage region
        #plt.tripcolor(triang, cov.flatten(), shading='gouraud', cmap=plt.cm.jet)
        #plt.figure(figsize=(10,10))
        print extent_c
        f,a = self.show(contour=False,bldg=True,height=False,coord='cartesian',extent=extent_c)
        plt.tripcolor(triang, Ltot.flatten(), shading='gouraud',
                      cmap=plt.cm.jet,vmax=-50,vmin=-130)
        plt.colorbar()

        plt.axis('equal')

        return x,y,r,cov,LOS,hearth,diff,fac,num,LFS


    def show(self,**kwargs):
        """ Ezone vizualization 

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

        """
        defaults = {'title':'title',
                    'xlabel':'Longitude',
                    'ylabel':'Latitude',
                    'height':True,
                    'bldg':True,
                    'coord':'lonlat',
                    'extent':[],
                    'contour':False,
                   }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        # get zone limitation
        # lon,lat or cartesian

        if kwargs['extent']==[]:
            if kwargs['coord']=='lonlat':
                extent = self.extent
            if kwargs['coord']=='cartesian':
                extent_c = (self.x.min(),self.x.max(),self.y.min(),self.y.max())
        else:
            if kwargs['coord']=='cartesian':
                extent_c = kwargs['extent']
                extent = conv(extent_c,self.m,mode='toll')
            if kwargs['coord']=='lonlat':
                extent = kwargs['extent']


        fig = plt.figure(figsize=(10,10))
        ax  = fig.add_subplot(111)

        if kwargs['bldg']:
            for city in self.dcity.keys():
                bdpt = self.dcity[city]['bdpt']
                bdma = self.dcity[city]['bdma'].astype('bool')

                # include masked point in the middle of the extent zone
                bdpt[bdma[:,0],0] = (extent[0]+extent[1])/2.
                bdpt[bdma[:,0],1] = (extent[2]+extent[3])/2.

                u = np.where((bdpt[:,0]>=extent[0]) &
                             (bdpt[:,0]<=extent[1]) &
                             (bdpt[:,1]>=extent[2]) &
                             (bdpt[:,1]<=extent[3]) )[0]
                vtx = np.ma.masked_array(bdpt[u,:],bdma[u,:])
                # if cartesian mode do the basemap conversion
                if kwargs['coord']=='cartesian':
                    # select zone
                    # convert into cartesian
                    vx,vy  = self.m(bdpt[u,0],bdpt[u,1])
                    bdmau  = bdma[u,:]
                    bdpt_c = np.hstack((vx[:,np.newaxis],vy[:,np.newaxis]))
                    vtx = np.ma.masked_array(bdpt_c,bdmau)


                ax.plot(vtx[:,0],vtx[:,1],linewidth=0.5,color='k')

        if kwargs['coord'] == 'cartesian':
            kwargs['xlabel'] = 'W-E Distance (meters)'
            kwargs['ylabel'] = 'N-S Distance (meters)'

        ax.set_title(kwargs['title'])
        ax.set_xlabel(kwargs['xlabel'])
        ax.set_ylabel(kwargs['ylabel'])

        if (kwargs['height'] | kwargs['contour']):
            shaphgt = self.hgt.shape
            # full original x and y
            if kwargs['coord']=='lonlat':
                x = self.lon
                y = self.lat
                hgt = self.hgt
            if kwargs['coord']=='cartesian':
                x = self.x
                y = self.y
                hgt = self.hgt_cart
                extent = extent_c

            # index corresponding to the selected zone
            ix = np.where((x>=extent[0]) & (x<=extent[1]))[0]
            iy = np.where((y>=extent[2]) & (y<=extent[3]))[0]

            if kwargs['height']:
                im = ax.imshow(hgt[iy[0]:(iy[-1]+1),ix[0]:(ix[-1]+1)],extent=extent)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cb = fig.colorbar(im,cax)
                cb.set_label('Height (meters)')
            if kwargs['contour']:
                cnt = ax.contour(hgt[iy[0]:(iy[-1]+1),ix[0]:(ix[-1]+1)],N=10,extent=extent,origin='upper')


        return(fig,ax)

    def loadtmp(self,_fileh5='RennesFull.h5'):
        """ load an Ezone from hdf5 file

        Parameters
        ----------

        _fileh5 :

        """

        fileh5 = pyu.getlong(_fileh5,'gis/h5')
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

    def saveh5(self,_fileh5='N48W002.h5'):
        """ save Ezone in hdf5 format

        Parameters
        ----------

        _fileh5

        """
        fileh5 = pyu.getlong(_fileh5,'gis/h5')
        #try:
        f = h5py.File(fileh5,'w')
        f.create_group('osm')
        f['extent'] = self.extent
        for city in self.dcity:
            dt = self.dcity[city]
            f['osm'].create_group(city)
            f['osm'][city].create_dataset('bdpt',shape=dt['bdpt'].shape,data=dt['bdpt'])
            f['osm'][city].create_dataset('bdma',shape=dt['bdma'].shape,data=dt['bdma'])
            f['osm'][city].create_dataset('extent',data=dt['extent'])
        #f['osm'].create_group('Rennes')
        #f['osm']['Rennes'].create_dataset('bdpt',shape=self.bdpt.shape,data=self.bdpt)
        #f['osm']['Rennes'].create_dataset('bdma',shape=self.bdma.shape,data=self.bdma)
        #f['osm']['Rennes'].create_dataset('extent',data=self.extentc)

        f.create_dataset('hgt',shape=self.hgt.shape,data=self.hgt)
        f.create_dataset('lcv',shape=self.lcv.shape,data=self.lcv)
        f.close()
        #except:
        #    f.close()
        #    raise NameError('error opening file')


    def build(self,_fileosm,_filehgt,_filelcv):
        """ build an Ezone from heterogeneous files
        """
        pass

if __name__== "__main__":
    ez = Ezone()
    ez.load('N48W002.h5')
    extent1 = (-1.7,-1.6,48.05,48.15)
    extent1_cart = conv(extent1,ez.m,mode='tocart')
