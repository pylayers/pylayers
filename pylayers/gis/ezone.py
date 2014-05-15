#-*- coding:Utf-8 -*-
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
import pyutil as pyu
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
"""

Ezone class
===========

.. autosummary::
    :toctree: generated/

"""

def maxloc(a):
    """
    Determine local maximum

    Parameters
    ----------

    a

    """
    a_decr = np.roll(a,1,axis=1)
    a_decl = np.roll(a,-1,axis=1)
    ind = np.where((a>a_decr)&(a>a_decl)&(a>-0.7))
    #pdb.set_trace()
    #ai = a[ind]
    #i_max = np.where(ai==max(ai))[0]
    #ai_max = ai[i_max]
    b = np.zeros(np.shape(a))
    #
    b[ind] = a[ind]
    return(b)

def decsrtm(_filehgt='N48W002.HGT'):
    """
        decode srtm file

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
        lonmin = -eval(_filehgt[5:8])
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
    """
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

class DEM(object):
    """
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
        """
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

class Ezone(object):
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
        """
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

    def link(self,pt,pr,**kwargs):
        """ ??
        """
        defaults = {
            'tmp':0
        }

    def extract(self,pt,pr,**kwargs):
        """
        """

        defaults = {'Npt':1000,
                    'ht':30,
                    'hr':1.5,
                    'K':1.3333,
                    'fGHz':.3}

        for key in defaults:
            if key not in kwargs:
                kwargs[key] = defaults[key]

        # wavelength
        lmbda = 0.3/kwargs['fGHz']

        # transmitter coordinates
        x_t,y_t = self.m(pt[0],pt[1])

        # receiver coordinates
        x_r,y_r = self.m(pr[0],pr[1])

        x = np.linspace(x_t,x_r,kwargs['Npt'])
        y = np.linspace(y_t,y_r,kwargs['Npt'])

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

        ha = height[0] + kwargs['ht']
        hb = height[-1]+ kwargs['hr']
        LOS = ha+(hb-ha)*d/d[-1]
        diff = height-LOS
        fac  = np.sqrt(2*d[-1]/(lmbda*d*d[::-1]))
        nu   = diff*fac
        num  = maxloc(nu[np.newaxis,:])

        plt.plot(d,dh,'r',d,height,'b',d,m[0,:],d,LOS,'k')
        #plt.figure()
        #plt.plot(d,nu,num)

        return(height,d,nu,num)

    def cover(self,**kwargs):
        """
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
        """

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
        """
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
    #pt = ez.dcity['Rennes']['bdpt']
    #ma = ez.dcity['Rennes']['bdma']
    #vertices = np.ma.masked_array(pt, ma)
    #lonmin = vertices[:,0].min()
    #lonmax = vertices[:,0].max()
    #latmin = vertices[:,1].min()
    #latmax = vertices[:,1].max()
    #ez.dcity['Rennes']['extent'] = (lonmin,lonmax,latmin,latmax)
    extent1 = (-1.7,-1.6,48.05,48.15)
    extent1_cart = conv(extent1,ez.m,mode='tocart')
    #f,a = ez.show(title='Rennes City Center',extent=extent1)
    #f,a = ez.show(title='Rennes City Center',extent=extent1_cart,coord='cartesian',bldg=True)
    #f,a = ez.show(title='Beaulieu Campus',extent=(-1.655,-1.624,48.101,48.125),coord='cartesian')
    #plt.show()

    #ez.show()
    #plt.show()
    #ez.loadtmp()
    #ez.saveh5()
#
##  Read the h5py file ( building from OSM in masked arrays)
#f = h5py.File('Rennes1.h5','r',dtype=np.float32)
#pt = f['points'].value
##u = np.where(((pt[:,0]>-1.75) & (pt[:,0]<-1.60)))
#ma = f['mask'].value
##vertices = np.ma.masked_array(pt[u[0],:], ma[u[0],:])
#vertices = np.ma.masked_array(pt, ma)
#f.close()
#
##vertices = vertices[
## getting extrema from the masked array (important)
#
#lonmin = vertices[:,0].min()
#lonmax = vertices[:,0].max()
#latmin = vertices[:,1].min()
#latmax = vertices[:,1].max()
#
#print lonmin,lonmax
#print latmin,latmax
#
#f = h5py.File(fileh5,'r',dtype=np.float32)
## construct the lon,lat grid
##lon = np.linspace(-2,-1,1201)[:,np.newaxis]
##lat = np.linspace(48,49,1201)[np.newaxis,:]
#lon = np.linspace(-2,-1,1201)
#lat = np.linspace(49,48,1201)
#x,y = m(lon,lat)
## construct the grid in cartesian coordinates
## origin is at the lower left corner
#X,Y  = np.meshgrid(x,y)
#eps = 0
#ilon = np.where((lon>=(lonmin-eps)) & (lon<=(lonmax+eps)))[0]
#ilat = np.where((lat>=(latmin-eps)) & (lat<=(latmax+eps)))[0]
#
## Relire le fichier HGT - DEM
#
#_filename = '../../../data/hgt/N48W002.HGT'
#dem = np.fromfile(_filename,dtype='>i2')
#dem = dem.reshape(1201,1201)
#fig = plt.figure(figsize=(10,20))
#ax1  = fig.add_subplot(111)
#
#
#im1 = ax1.imshow(dem[ilat[0]:(ilat[-1]+1),ilon[0]:(ilon[-1]+1)],extent=(lonmin,lonmax,latmin,latmax))
#
#divider = make_axes_locatable(ax1)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#cb = fig.colorbar(im1,cax)
#cb.set_label('Height (meters)')
#
#vertices = np.ma.masked_array(pt, ma)
#ax1.plot(vertices[:,0],vertices[:,1],linewidth=0.5,color='k')
#ax1.set_title('En longeant la Vilaine avec OSM')
#ax1.set_xlabel('Longitude')
#ax1.set_ylabel('Latitude')
#plt.show()
#
##plt.savefig('Rennes.png')
## Relire le fichier LCV
#
#path = '/home/uguen/data/lcv'
#_filename = 'N48W002.lcv'
#filename = path+ '/'+_filename
#data = np.fromfile(filename,dtype='>i1')
#sursol = data.reshape(1201,1201)
#
##plt.figure(figsize=(20,20))
##plt.imshow(sursol[ilon,ilat])
##plt.colorbar()
#
#
##def convert(lon,lat):
##    """  conversion  (lon,lat) --> (x,y)
##
##    Parameters
##    ----------
##
##    lon : np.array (,Np)
##    lat : np.array (,Np)
##
##    """
##    lonmin = min(lon)
##    lonmax = max(lon)
##    latmin = min(lat)
##    latmax = max(lat)
##
##    lon_0 = (lonmax+lonmin)/2
##    lat_0 = (latmax+latmin)/2
##
##    m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,
##            resolution='i',projection='cass',lon_0=lon_0,lat_0=lat_0)
##    x,y = m(lat,lon)
##    return x,y
##
##
##
### Each data file covers a one-degree-of-latitude by one-degree-of-longitude
### block of Earth's surface. The first seven characters indicate the southwest
### corner of the block, with N, S, E, and W referring to north, south, east, and
### west. Thus, the "N48W002.hgt" file covers latitudes 48 to 49 North and
### longitudes -2 to -1 West (this file includes Rennes). The filename extension
###" .hgt" simply stands for the word "height", meaning elevation. 
##
##
##
##lat = linspace(48,49,1201).reshape(1201)
##lon = linspace(-2,-1,1201).reshape(1201)
##
##x,y = m(lat,lon)
##city = {}
##city['Rennes']=[-1.67944,48.114722]
##city['Servon']=[-1.45944,48.12222]
##
### <codecell>
##
##_filename = 'N48W002.HGT'
##mnt = fromfile(_filename,dtype='>i2')
##mnt= flipud(mnt.reshape(1201,1201))
##        
##
### <codecell>
##
##figsize(10,10)
##
###M = mnt.max()
###print(M)
###u = where(mnt==M)
###print(u)
###plot(1169,919,'ok')
##plot(lon[1169],lat[919],'ok')
##plot(city['Rennes'][0],city['Rennes'][1],'ok')
##plot(city['Servon'][0],city['Servon'][1],'ok')
##imshow(mnt,cmap=cm.hsv,origin='lower',extent=[-2,-1,48,49])
##xlabel('Longitude (degrees)')
##ylabel('Latitude (degrees)')
##colorbar()
##
### <markdowncell>
##
### First dimension is  latitude second dimension is longitude
##
### <markdowncell>
##
### create an interpolation function
##
### <codecell>
##
##from scipy import interpolate
##import numpy as np
##
##def my_interp(X, Y, Z, x, y, spn=3):
##    xs,ys = map(np.array,(x,y))
##    z = np.zeros(xs.shape)
##    for i,(x,y) in enumerate(zip(xs,ys)):
##        # get the indices of the nearest x,y
##        xi = np.argmin(np.abs(X[0,:]-x))
##        yi = np.argmin(np.abs(Y[:,0]-y))
##        xlo = max(xi-spn, 0)
##        ylo = max(yi-spn, 0)
##        xhi = min(xi+spn, X[0,:].size)
##        yhi = min(yi+spn, Y[:,0].size)
##        # make slices of X,Y,Z that are only a few items wide
##        nX = X[xlo:xhi, ylo:yhi]
##        nY = Y[xlo:xhi, ylo:yhi]
##        nZ = Z[xlo:xhi, ylo:yhi]
##        intp = interpolate.interp2d(nX, nY, nZ)
##        z[i] = intp(x,y)[0]
##    return z
##
### <codecell>
##
##a    = Building()
##node = Nodes()
##p    = OSMParser(concurrency=4, ways_callback=a.ways)
##q    = OSMParser(concurrency=4, coords_callback=node.coords)
##
### <codecell>
##
##q.parse('map.osm')
##
### a.ways
##p.parse('map.osm')
##
##
### <codecell>
##
##latitude = []
##longitude =[]
##for b in a.building:
##    f =  a.building[b]['ref']
##    N = len(f)
##    p = np.zeros((2,N))
##    for k,id in enumerate(f):
##        longitude.append(node.latlon[id][0])
##        latitude.append(node.latlon[id][1])
##
###coord = node.latlon.values()
###longitude = map(lambda x : x[0],coord)
###latitude = map(lambda x : x[1],coord)
##lon = linspace(-2,-1,1201)
##lat = linspace(48,49,1201)
##
##lonmin = min(longitude)
##lonmax=max(longitude)
##latmin = min(latitude)
##latmax=max(latitude)
##
##steplon = (lon[-1]-lon[0])/1201
##steplat = (lat[-1]-lat[0])/1201
##
##ilonSW = int(floor((lonmin-lon[0])/steplon))
##ilatSW = int(floor((latmin-lat[0])/steplat))
##
##ilonNE = int(ceil((lonmax-lon[0])/steplon))
##ilatNE = int(ceil((latmax-lat[0])/steplat))
##
##print ilonSW,ilatSW
##print ilonNE,ilatNE
##
##print lonmin,lonmax
##print latmin,latmax   
##
##indx = range(ilonSW,ilonNE)
##indy = range(ilatSW,ilatNE)
##lon_cut = lon[indx]
##lat_cut = lat[indy]
##np.shape(mnt)
##mnt_cut = mnt[ilatSW:ilatNE,ilonSW:ilonNE]
###imshow(mnt_cut,extent=[lon_cut[0],lon_cut[-1],lat_cut[0],lat_cut[-1]])
##imshow(mnt_cut,origin='lower',interpolation='nearest')
##colorbar()
##print np.shape(mnt_cut)
##print np.shape(lon_cut)
##print np.shape(lat_cut)
##             
##
### <codecell>
##
##lon_new=linspace(lonmin,lonmax,30)
##lat_new=linspace(latmin,latmax,15)
##LAT,LON = meshgrid(lat_cut,lon_cut)
##finterpol=interpolate.interp2d(LON,LAT,mnt_cut,kind='linear',bounds_error=True)
##mntnew = finterpol(lon_new,lat_new)
##
##print np.shape(LON)
##print np.shape(LAT)
##print np.shape(mnt_cut)
##
##
##
##
### <codecell>
##
##imshow(LON)
##imshow(LAT)
##imshow(mnt_cut,origin='lower')
##
### <codecell>
##
##imshow(mntnew,origin='lower')
##
### <codecell>
##
###m = Basemap(llcrnrlon=-1.65263,llcrnrlat=48.1127,urcrnrlon=-1.62759,urcrnrlat=48.12547,
###            resolution='i',projection='cass',lon_0=-1.63,lat_0=48.115)
##m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,
##            resolution='i',projection='cass',lon_0=-1.642045,lat_0=48.119423)
##
##xT,yT=m(latT,lonT)
##
##Sb={}
##Sb['TourMath'] = {}
##Sb['TourMath']['latlon'] = (-1.642045,48.119423)
##Sb['TourMath']['pos'] = m(-1.642045,48.119423)
##Sb['Chardonnet'] = {}
##Sb['Chardonnet']['latlon'] = (-1.642218,48.107158)
##Sb['Chardonnet']['pos'] = m(-1.642218,48.107158)
##Sb['Chardonnet']['angle'] = [30,150,270]
##Sb['Chardonnet']['hsol'] =  [30,30,30]
##Sb['Ibis'] = {}
##Sb['Ibis']['Id'] = 322746
##Sb['Ibis']['latlon']=(-1.621842,48.114650)
##Sb['Ibis']['pos'] = m(-1.621842,48.114650)
##Sb['Ibis']['angle'] = [10,150,270]
##Sb['Ibis']['hsol'] =  [14,14,15]
##Sb['TourHertz'] = {}
##Sb['TourHertz']['Id'] = 357933
##Sb['TourHertz']['latlon'] = (-1.6227124,48.124767)
##Sb['TourHertz']['pos'] = m(-1.627085,48.124789)
##Sb['TourHertz']['angle'] = [30,150,270]
##Sb['TourHertz']['hsol'] =  [24,24,24]
##
### <codecell>
##
##print lon_cut[0], lon_cut[-1], lat_cut[0], lat_cut[-1]
##xl,yl=m(lon_cut[0],lat_cut[0])
##print xl,yl
##xh,yh=m(lon_cut[-1],lat_cut[-1])
##print xh,yh
##
### <codecell>
##
##fig = plt.figure(figsize=(40,40))
##ax  = fig.gca()
##
##for b in a.building:
##    f =  a.building[b]['ref']
##    N = len(f)
##    p = np.zeros((2,N))
##    for k,id in enumerate(f):
##        x,y = m(node.latlon[id][0],node.latlon[id][1])
##        p[0,k] = x
##        p[1,k] = y
##    Pol=Polygon(p)
##    fig,ax = Pol.plot(fig=fig,ax=ax)    
## 
###ax.imshow(mnt_cut,extent=[xl,xh,yl,yh],alpha=0.5,origin='lower')
###ax.plot(Sb['TourMath']['pos'][0],Sb['TourMath']['pos'][1],'or')
###ax.plot(Sb['Chardonnet']['pos'][0],Sb['Chardonnet']['pos'][1],'or')
###ax.plot(Sb['Ibis']['pos'][0],Sb['Ibis']['pos'][1],'or')
###ax.plot(Sb['TourHertz']['pos'][0],Sb['TourHertz']['pos'][1],'or')
###ax.scatter(xT,yT,c=radio.rssi,cmap=cm.hot,linewidths=0)
##
##axis('scaled')
##fig.savefig('batcampus.png')
##
### <codecell>
##
##min(xT)
##
### <codecell>
##
##max(yT)
##
### <codecell>
##
##Image(url='http://maps.googleapis.com/maps/api/streetview?size=640x400&location=48.119423,-1.642045&heading=180&pitch=40&fov=120&sensor=false',format='jpeg')
##
### <codecell>
##
##def showsite(lat=48.119423,lon=-1.642045,heading=180,pitch=40,fov=120):
##    """
##    Parameters
##    ----------
##    lat : 
##    lon : 
##    heading :
##    pitch : 
##    fov : 
##    
##    """
##    urlbase = 'http://maps.googleapis.com/maps/api/streetview?size=640x400'
##    location = str(lat)+','+str(lon)
##    heading = str(heading)
##    pitch = str(pitch)
##    fov = str(fov)
##    url = urlbase+'&location='+location+'&heading='+heading+'&pitch='+pitch+'&fov='+fov+'&sensor=false'
##    return(Image(url=url,format='jpeg'))
##
### <codecell>
##
##showsite()
##
### <codecell>
##
##showsite(lon=Sb['Chardonnet']['latlon'][0],lat=Sb['Chardonnet']['latlon'][1],heading=180)
##
### <codecell>
##
##showsite(lat=48.124571,lon=-1.627667,heading=90)
##
### <codecell>
##
##Sb['TourHertz']['latlon'][0]
##
### <codecell>
##
##Sb['TourHertz']['latlon'][1]
##
### <codecell>
##
##cartodir = '/private/staff/n/en/buguen/data/cartoradio/Rennes'
##
### <codecell>
##
##BS = pd.read_csv(cartodir+'/RennesPosition.csv')
##
### <codecell>
##
##index_orange=[]
##for j,k in enumerate(BS.Proprietaire):
##    if k=="ORANGE":
##        index_orange.append(j)
##        
##print index_orange        
##
### <codecell>
##
##site={}
###for k in BS.Latitude.keys():
##for j,k in enumerate(index_orange):
##   lat = eval(BS.Latitude[k].replace(',','.'))
##   lon = eval(BS.Longitude[k].replace(',','.'))
##   site[j]=[lat,lon]
##
### <codecell>
##
##showsite(lat=site[4][0],lon=site[4][1],heading=90,fov=100)
##
### <codecell>
##
##print lon_cut[0], lon_cut[-1], lat_cut[0], lat_cut[-1]
##xl,yl=m(lon_cut[0],lat_cut[0])
##print xl,yl
##xh,yh=m(lon_cut[-1],lat_cut[-1])
##print xh,yh
##
##fig = plt.figure(figsize=(20,20))
##ax  = fig.gca()
##
##pts = []
##for b in a.building:
##    f =  a.building[b]['ref']
##    N = len(f)
##    p = np.zeros((2,N))
##    for k,id in enumerate(f):
##        x,y = m(node.latlon[id][0],node.latlon[id][1])
##        p[0,k] = x
##        p[1,k] = y
##    Pol=Polygon(p)
##    fig,ax = Pol.plot(fig=fig,ax=ax)    
## 
##ax.imshow(mnt_cut,extent=[xl,xh,yl,yh],alpha=0.5,origin='lower')
##for k in site:
##    x,y  = m(site[k][1],site[k][0])
##    pts.append(Site(x,y))
##    ax.plot(x,y,'ob')
##axis([-2000,2000,-2000,2000])
##
### <codecell>
##
##c = Context()
##c.plot=1
##c.triangulate=1
##
### <codecell>
##
##sl = SiteList(pts) 
##
### <codecell>
##
##figsize(10,10)
##voronoi(sl,c)
##
### <codecell>
##
##c= pat.Circle((0,0),4)
##plt.show()
##
### <codecell>
##
##Image(url='http://maps.google.com/maps?q=wendys,84020&zoom=14&size=512x512&maptype=roadmap&sensor=false',format='jepg')
##
### <codecell>
##
##
##
