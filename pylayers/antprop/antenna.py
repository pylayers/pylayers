# -*- coding:Utf-8 -*-
"""

.. currentmodule:: pylayers.antprop.antenna

This module handles antennas
An antenna can be loaded from various file formats among

+ .vsh2 
+ .vsh3 
+ .sh2 
+ .sh3 
+ .mat
+ .trx

Examples
--------
    >>> import matplotlib.pyplot as plt
    >>> from pylayers.antprop.antenna import *
    >>> A = Antenna('defant.trx')
    >>> fig,ax = A.polar(fGHz=[2,3,4],phd=0)
   

Antenna Class
=============

.. autosummary::
    :toctree: generated/


Utility Functions
-----------------

.. autosummary::
    :toctree: generated/

    Antenna.__init__
    Antenna.__repr__
    Antenna.ls

    Antenna.errel
    Antenna.checkpole
    Antenna.info


    Antenna.pol2cart
    Antenna.cart2pol
    Antenna.minsh3
    Antenna.mse
    Antenna.getdelay
    Antenna.elec_delay

Synthesis Functions
-------------------

.. autosummary::
    :toctree: generated/

    Antenna.Fsynth
    Antenna.Fsynth1
    Antenna.Fsynth2s
    Antenna.Fsynth2b
    Antenna.Fsynth2
    Antenna.Fsynth3
    Antenna.Fpatt

Visualization functions
-----------------------

.. autosummary::
    :toctree: generated/

    Antenna.pattern
    Antenna.polar
    Antenna._show3
    Antenna.show3
    Antenna.plot3d
    Antenna.pol3d
    Antenna.load_trx
    Antenna.movie_vsh

Loading and Saving
------------------

.. autosummary::
    :toctree: generated/

    Antenna.loadhfss
    Antenna.loadtrx
    Antenna.loadmat
    Antenna.savevsh3
    Antenna.savesh2
    Antenna.savesh3
    Antenna.loadvsh3
    Antenna.loadsh3
    Antenna.savevsh2
    Antenna.loadsh2
    Antenna.loadvsh2

Miscellianous  functions
========================


.. autosummary::
    :toctree: generated/

    forcesympol
    compdiag
    show3D


"""

try:
    #from mayavi import mlab
    import mayavi.mlab as mlab
except:
    print 'mayavi not installed'

import doctest
import os
import glob
import re
import sys
import pdb
import Image
import numpy as np
import scipy as sp
import scipy.linalg as la
import scipy.special as special
from scipy import io
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
from pylayers.util.project import *
from pylayers.antprop.spharm import *
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import axes3d
#from scipy import sparse
from matplotlib import rc
from matplotlib import cm # colormaps
from pylayers.antprop.antssh import *
from pylayers.antprop.coeffModel import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
import pandas as pd

import matplotlib.pylab as plt
class Antenna(PyLayers):
    """ Antenna

    Attributes
    ----------

    name   : Antenna name

    Nf     : number of frequency
    Nt     : number of theta
    Np     : number of phi

    Ftheta : Normalized Ftheta    (nf,ntheta,nphi)
    Fphi   : Normalized Fphi      (nf,ntheta,nphi)
    SqG    : Square root of gain  (nf,ntheta,nphi)

    Ttheta : theta                (nf,ntheta,nphi)
    Tphi   : phi                  (nf,ntheta,nphi)

    theta  : theta base            1 x ntheta
    phi    : phi base              1 x phi

    C      : VSH Coefficients

    Methods
    -------

    info  : Display information about antenna
    vsh   : calculates Vector Spherical Harmonics
    show3 : Geomview diagram
    plot3d : 3D diagram plotting using matplotlib toolkit
    Fsynth2 : Antenna F synthesis from coeff in s2

    Antenna trx file can be stored in various order
        natural : HFSS
        bcp     : near filed chamber

    It is important when initializing an antenna object to be aware of the typ of trx file


    .trx (ASCII Vectorial antenna Pattern)

    F   Phi   Theta  Fphi  Ftheta

    """


    def __init__(self,typ='S2R2.sh3',**kwargs):
        """ class constructor

        Parameters
        ----------

        typ  : 'Omni','Gauss','WirePlate','3GPP'

        _filename : string
                    antenna file name
        directory : str
                    subdirectory of the current project where to find the antenna file
                    the file is seek in the $BASENAME/ant directory
        nf        : integer
                     number of frequency (default 104)
        ntheta    : integer
                    number of theta (default 181)
        nph       : integer
                    number of phi (default 90)
        source    : string
                source of data { 'satimo' | 'cst' | 'hfss' }
        Notes
        -----

        The supported data formats for storing antenna patterns are

        'mat': Matlab File
        'vsh2': unthresholded vector spherical coefficients
        'vsh3': thresholded vector spherical cpoefficients
        'trx' : Satimo NFC raw data
        'trx1' : Satimo NFC raw data  (deprecated)

         A = Antenna('my_antenna.mat')

        """
        defaults = {'directory': 'ant',
                    'nf':104,
                    'pol':'c',
                    'source':'satimo',
                    'ntheta':90,
                    'nphi':181,
                    'p0':0,
                    't0':np.pi/2.,
                    'p3':np.pi/6.,
                    't3':np.pi/6.,
                    'L':90,
                    'fmin':2,
                    'fmax':10,
                    'gm': 18,
                    'sllv':-18,
                    'hpbwv':6.2,
                    'hpbwh':65,
                    'fbrh':30,
                    'thtilt':0}

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        self.Nf = kwargs['nf']
        self.Nt = kwargs['ntheta']
        self.Np = kwargs['nphi']
        self.source = kwargs['source']
        self.pol = kwargs['pol']
        self.fmin = kwargs['fmin']
        self.fmax = kwargs['fmax']

        # if typ has an extension it is a file
        if isinstance(typ,str):
            AntennaName,Extension = os.path.splitext(typ)
            self.ext = Extension[1:]
            if self.ext=='':
                self.fromfile = False
            else:
                self.fromfile = True
        else:
            self.fromfile = True

        self.tau = 0
        self.evaluated = False

        if self.fromfile:
            if isinstance(typ,str):
                self._filename = typ
                if self.ext == 'vsh3':
                    self.typ='vsh'
                    self.loadvsh3()
                if self.ext == 'vsh2':
                    self.typ='vsh'
                    self.loadvsh2()
                if self.ext == 'sh3':
                    self.typ='ssh'
                    self.loadsh3()
                if self.ext == 'sh2':
                    self.typ='ssh'
                    self.loadsh2()
                if self.ext == 'trx1':
                    self.typ='trx'
                    self.load_trx(kwargs['directory'],self.Nf,self.Nt,self.Np)
                if self.ext == 'trx':
                    self.typ='trx'
                    self.loadtrx(kwargs['directory'])
                if self.ext == 'mat':
                    self.typ='mat'
                    self.loadmat(kwargs['directory'])
            elif isinstance(typ,list):
                self._filename = typ
                self.ext='hfss'
                self.loadhfss(typ, self.Nt, self.Np)

        else :
            self._filename = typ
            self.typ = typ
            if typ == 'Gauss':
                self.p0 = kwargs['p0']
                self.t0 = kwargs['t0']#np.pi/2.
                self.p3 = kwargs['p3']#np.pi/6. # 30 degrees
                self.t3 = kwargs['t3']#np.pi/6. # 30 degrees
                self.G  = 16/(self.t3*self.p3) # gain 16/\theta^{2}
                self.GdB = 10*np.log10(self.G)
                self.sqG = np.sqrt(self.G)
                self.evaluated = False
            elif typ == 'WirePlate':
                self.p0 = kwargs['p0']
                kwargs['t0'] = 5*np.pi/6.
                self.t0 =  kwargs['t0']#
                self.GdB = 5. # gain
                self.G  = pow(10.,self.GdB/10.) # gain
                self.sqG = np.sqrt(self.G)
                self.evaluated = False
            elif typ == 'Omni':
                self.GdB  = 0. # gain
                self.G  = pow(10.,self.GdB/10.) # gain
                self.sqG = np.sqrt(self.G)
                self.evaluated = False
            elif typ == '3GPP':
                self.sllv = kwargs['sllv']
                self.hpbwv = kwargs['hpbwv']
                self.hpbwh = kwargs['hpbwh']
                self.fbrh = kwargs['fbrh']
                self.thtilt = kwargs['thtilt']
                self.gm = kwargs['gm']



            elif typ == 'ssh':
                pass
            elif typ == 'vsh':
                L = kwargs['L']
                fmin = kwargs['fmin']
                fmax = kwargs['fmax']
                nf = kwargs['nf']
                self.fa = np.linspace(fmin,fmax,nf)
                # initialize a VSHCoeff with zeros

                dBr = np.zeros((nf,L+1,L),dtype='complex128')
                dBi = np.zeros((nf,L+1,L),dtype='complex128')
                dCr = np.zeros((nf,L+1,L),dtype='complex128')
                dCi = np.zeros((nf,L+1,L),dtype='complex128')

                Br = VCoeff(typ='s1',fmin=fmin,fmax=fmax,data=dBr)
                Bi = VCoeff(typ='s1',fmin=fmin,fmax=fmax,data=dBi)
                Cr = VCoeff(typ='s1',fmin=fmin,fmax=fmax,data=dCr)
                Ci = VCoeff(typ='s1',fmin=fmin,fmax=fmax,data=dCi)

                self.C = VSHCoeff(Br,Bi,Cr,Ci)

            else:
                raise NameError('antenna typ is not known')


    def __repr__(self):

        rtd = 180./np.pi
        st = ''
        if self.fromfile:
            if isinstance(self._filename,str):
                st = st + 'FileName : ' + self._filename+'\n'
                st = st + '-----------------------\n'
            else:
                for i in range(len(self._filename)):
                    st = st + 'FileName : ' + self._filename[i]+'\n'
                st = st + '-----------------------\n'
        #st = st + 'file type : ' + self.typ+'\n'
        if 'fa' in self.__dict__:
            st = st + "fmin : %4.2f" % (self.fa[0]) + "GHz\n"
            st = st + "fmax : %4.2f" % (self.fa[-1]) + "GHz\n"
            try:
                st = st + "step : %4.2f" % (1000*(self.fa[1]-self.fa[0])) + "MHz\n"
            except:
                st = st + "step : None\n"
            st = st + "Nf : %d" % (len(self.fa)) +"\n"


        if self.evaluated:
            st = st + '-----------------------\n'
            st = st + "Ntheta : %d" % (self.Nt) + "\n"
            st = st + "Nphi : %d" % (self.Np) + "\n"
            u = np.where(self.SqG==self.SqG.max())
            if len(u[0]>1):
                S = self.SqG[(u[0][0],u[1][0],u[2][0])]
                uf = u[0][0]
                ut = u[1][0]
                up = u[2][0]
            else:
                S = self.SqG[u]
                uf = u[0]
                ut = u[1]
                up = u[2]
            if self.source=='satimo':
                GdB = 20*np.log10(S)
            # see WHERE1 D4.1 sec 3.1.1.2.2
            if self.source=='cst':
                GdB = 20*np.log10(S/np.sqrt(30))
            st = st + "GmaxdB : %4.2f dB \n" % (GdB)
            st = st + "   f = %4.2f GHz \n" % (self.fa[uf])
            st = st + "   theta = %4.2f (degrees) \n" % (self.theta[ut]*rtd)
            st = st + "   phi = %4.2f  (degrees) \n" % (self.phi[up]*rtd)
        else:
            st = st + 'Not evaluated\n'


        if self.typ == 'mat':
            #st = st + self.DataFile + '\n'
            st = st + 'antenna name : '+ self.AntennaName + '\n'
            st = st + 'date : ' + self.Date +'\n'
            st = st + 'time : ' + self.StartTime +'\n'
            st = st + 'Notes : ' + self.Notes+'\n'
            st = st + 'Serie : ' + str(self.Serie)+'\n'
            st = st + 'Run : ' + str(self.Run)+'\n'
            st = st + "Nb theta (lat) : "+ str(self.Nt)+'\n'
            st = st + "Nb phi (lon) :"+ str(self.Np)+'\n'

        if self.typ == 'Gauss':
            st = st + 'Gaussian pattern' + '\n'
            st = st + 'phi0 : ' + str(self.p0) +'\n'
            st = st + 'theta0 :' + str(self.t0) + '\n'
            st = st + 'phi 3dB :' + str(self.p3) + '\n'
            st = st + 'theta 3dB :' + str(self.t3) + '\n'
            st = st + 'Gain dB :' + str(self.GdB) + '\n'
            st = st + 'Gain linear :' + str(self.G ) + '\n'
            st = st + 'sqrt G :' + str(self.sqG) + '\n'

        return(st)

    def ls(self, typ='vsh3'):
        """ list the antenna files in antenna project directory

        Parameters
        ----------

        typ : string optional
            {'mat'|'trx'|'vsh2'|'sh2'|'vsh3'|'sh3'}

        Returns
        -------

        lfile_s : list
            sorted list of all the .str file of strdir

        """

        if typ=='vsh3':
            pathname = pstruc['DIRANT'] + '/*.' + typ
        if typ=='sh3':
            pathname = pstruc['DIRANT'] + '/*.' + typ
        if typ=='mat':
            pathname = pstruc['DIRANT'] + '/*.' + typ
        if typ=='trx':
            pathname = pstruc['DIRANT'] + '/*.' + typ

        lfile_l = glob.glob(basename+'/'+pathname)
        lfile_s = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_s.append(fis)
        lfile_s.sort()

        return lfile_s


    def photo(self,directory=''):
        """ show a picture of the antenna """

        if directory == '':
            directory = os.path.join('ant','UWBAN','PhotosVideos')

        _filename = 'IMG_'+self.PhotoFile.split('-')[1]+'.JPG'
        filename = pyu.getlong(_filename,directory)
        I = Image.open(filename)
        I.show()

    def Fpatt(self,th=[],ph=[],pattern=True):
        """  generate antenna pattern

        Parameters
        ----------

        th  :
        ph  :
        pattern : boolean
            default True

        if pattern
            The pattern is generated from self.Nt and self.Np points
        else:
            phi and theta have the same length (typically ray direction)


        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.antprop.antenna  import *
            >>> A = Antenna(typ='Gauss')
            >>> A.Fpatt()
            >>> f,a = A.polar()
            >>> plt.show()

        """

        assert not self.fromfile , 'Error : this is not a pattern antenna'

        self.fa = np.linspace(self.fmin,self.fmax,self.Nf)


        if (th == []) and (ph == []):
            self.theta = np.linspace(0,np.pi,self.Nt)
            self.phi = np.linspace(0,2*np.pi,self.Np,endpoint=False)
        else :
            self.theta = th
            self.phi = ph

        if self.typ == 'Gauss':


            argth = ((self.theta-self.t0)**2)/self.t3
            e1 = np.mod(self.phi-self.p0,2*np.pi)
            e2 = np.mod(self.p0-self.phi,2*np.pi)
            e = np.array(map(lambda x: min(x[0],x[1]),zip(e1,e2)))
            argphi = (e**2)/self.p3

            if pattern :
                Fat = self.sqG * ( np.exp(-2.76*argth[:,np.newaxis]) * np.exp(-2.76*argphi[np.newaxis,:]) )
                Fap = self.sqG * ( np.exp(-2.76*argth[:,np.newaxis]) * np.exp(-2.76*argphi[np.newaxis,:]) )
                #self.theta=self.th[:,np.newaxis]
                #self.phi=self.ph[np.newaxis,:]
                self.SqG=np.ones((self.Nf,self.Nt,self.Np))
                self.SqG[:]=Fap
                self.evaluated = True
            else:
                Fat = self.sqG * ( np.exp(-2.76*argth) * np.exp(-2.76*argphi) )
                Fap = self.sqG * ( np.exp(-2.76*argth) * np.exp(-2.76*argphi) )
                Fat = np.dot(Fat[:,np.newaxis],np.ones(len(self.fa))[np.newaxis,:])
                Fap = np.dot(Fap[:,np.newaxis],np.ones(len(self.fa))[np.newaxis,:])

        if self.typ == '3GPP':
            if pattern:
                phi   = self.phi*180/np.pi-180
                theta = self.theta*180/np.pi-90
                GvdB = np.maximum(-12*((theta-self.thtilt)/self.hpbwv)**2,self.sllv)[:,None]
                GhdB = (-np.minimum(12*(phi/self.hpbwh)**2,self.fbrh)+self.gm)[None,:]
                GdB  = GhdB+GvdB
                self.sqG = np.sqrt(10**(GdB/10.))
                self.SqG = self.sqG[None,...]*np.ones(self.Nf)[:,None,None]
                if self.pol=='h':
                    Fap = self.SqG
                    Fat = np.zeros((self.Nf,len(self.theta),len(self.phi)))
                if self.pol=='v':
                    Fap = np.zeros((self.Nf,len(self.theta),len(self.phi)))
                    Fat = self.SqG
                if self.pol=='c':
                    Fap = (1./np.sqrt(2))*self.SqG
                    Fat = (1j/np.sqrt(2))*self.SqG                    
                self.evaluated = True
            else:
                phi   = self.phi*180/np.pi-180
                theta = self.theta*180/np.pi-90
                GvdB = np.maximum(-12*((theta-self.thtilt)/self.hpbwv)**2,self.sllv)
                GhdB = (-np.minimum(12*(phi/self.hpbwh)**2,self.fbrh)+self.gm)
                GdB  = GhdB+GvdB
                self.sqG = np.sqrt(10**(GdB/10.))
                if self.pol=='h':
                    Fap = self.sqG
                    Fat = np.zeros((len(self.theta),self.Nf))
                if self.pol=='v':
                    Fap = np.zeros((len(self.phi),self.Nf))
                    Fat = self.sqG
                if self.pol=='c':
                    Fap = (1./sqrt(2))*self.sqG
                    Fat = (1j/sqrt(2))*self.sqG
                   

        if self.typ == 'WirePlate':

            uth1 = np.where(self.theta < self.t0)[0]
            uth2 = np.where(self.theta >= self.t0)[0]
            p = self.t0
            q = np.pi/2.
            A=np.array(([[3*p**2,2*p,1],[p**3,p**2,p],[q**3,q**2,q]]))
            Y=np.array(([0,1,1/(1.*self.sqG)]))
            self.poly = la.solve(A,Y)

            argth1 = np.abs(self.poly[0]*self.theta[uth1]**3
                     + self.poly[1]*self.theta[uth1]**2
                             + self.poly[2]*self.theta[uth1])

            argth2 = -(1/(np.pi-self.t0)**2)*(self.theta[uth2]-self.t0)**2+1
            argth = np.hstack((argth1,argth2))[::-1]

            if pattern :
                Fat = self.sqG * (argth[:,np.newaxis])
                Fap = self.sqG * (argth[:,np.newaxis])
                #self.theta=self.theta[:,np.newaxis]
                #self.phi=self.ph[np.newaxis,:]
                self.SqG=np.ones((self.Nf,self.Nt,self.Np))
                self.SqG[:]=Fap
                self.evaluated = True
            else:
                Fat = self.sqG * argth
                Fap = self.sqG * argth
                Fat = np.dot(Fat[:,np.newaxis],np.ones(len(self.fa))[np.newaxis,:])
                Fap = np.dot(Fap[:,np.newaxis],np.ones(len(self.fa))[np.newaxis,:])

        if self.typ == 'Omni':

            if pattern :

                self.SqG = self.sqG * np.ones((self.Nf,self.Nt,self.Np))
                #self.theta = self.theta[:,np.newaxis]
                #self.phi = self.ph[np.newaxis,:]
                self.evaluated = True
                Fat = self.sqG
                Fap = self.sqG

            else:
                Fat = self.sqG * np.ones((len(self.theta),self.Nf))
                Fap = self.sqG * np.zeros((len(self.theta),self.Nf))
        
        self.Ftheta = Fat
        self.Fphi = Fap
        
        # TODO create 2 separate functions
        if not pattern:
            return (Fat,Fap)
        else:
            return (None,None)


    def loadmat(self, directory="ant"):
        """ load an antenna stored in a mat file

        Parameters
        ----------
        directory : str , optional
            default 'ant'

        Examples
        --------

            Read an Antenna file in UWBAN directory and plot a polar plot

        .. plot::
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> from pylayers.antprop.antenna import *
            >>> A = Antenna('S1R1.mat',directory='ant/UWBAN/Matfile')
            >>> f,a = A.polar(phd=0)
            >>> f,a = A.polar(thd=90,fig=f,ax=a)
            >>> txt = plt.title('S1R1 antenna : st loadmat')
            >>> plt.show()


        """
        _filemat = self._filename
        filemat = pyu.getlong(_filemat, directory)
        d = io.loadmat(filemat, squeeze_me=True, struct_as_record=False)
        ext = _filemat.replace('.mat', '')
        d = d[ext]
        #
        #
        #
        self.typ = 'mat'
        self.Date = str(d.Date)
        self.Notes = str(d.Notes)
        self.PhotoFile = str(d.PhotoFile)
        self.Serie = eval(str(d.Serie))
        self.Run = eval(str(d.Run))
        self.DataFile = str(d.DataFile)
        self.StartTime = str(d.StartTime)
        self.AntennaName = str(d.AntennaName)

        self.fa = d.freq/1.e9
        #self.theta = d.theta[:,np.newaxis]
        #self.phi = d.phi[np.newaxis,:]
        self.theta = d.theta
        self.phi = d.phi
        self.Ftheta = d.Ftheta
        self.Fphi = d.Fphi
        Gr = np.real(self.Fphi * np.conj(self.Fphi) + \
                     self.Ftheta * np.conj(self.Ftheta))
        self.SqG = np.sqrt(Gr)
        #self.Nt = np.shape(self.theta)[0]
        #self.Np = np.shape(self.phi)[1]
        self.Nt = len(self.theta)
        self.Np = len(self.phi)

        if type(self.fa) ==  float:
            self.Nf = 1
        else:
            self.Nf = len(self.fa)

        self.evaluated = True

    def load_trx(self, directory="ant", nf=104, ntheta=181, nphi=90, ncol=6):
        """ load a trx file (deprecated)

        Parameters
        ----------
        directory : str
                    directory where is located the trx file (default : ant)
        nf : float
             number of frequency points
        ntheta : float
               number of theta
        nphi : float
               number of phi


        """
        _filetrx = self._filename
        filename = pyu.getlong(_filetrx, directory)
        if ncol == 6:
            pattern = """^.*\t.*\t.*\t.*\t.*\t.*\t.*$"""
        else:
            pattern = """^.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*$"""
        fd = open(filename, 'r')
        d = fd.read().split('\r\n')
        fd.close()
        k = 0
        #while ((re.search(pattern1,d[k]) is None ) & (re.search(pattern2,d[k]) is None )):
        while re.search(pattern, d[k]) is None:
            k = k + 1

        d = d[k:]
        N = len(d)
        del d[N - 1]
        r = '\t'.join(d)
        r.replace(' ', '')
        d = np.array(r.split()).astype('float')

        #
        # TODO Parsing the header
        #
        #nf     = 104
        #nphi   = 90
        #ntheta = 181

        N = nf * nphi * ntheta
        d = d.reshape(N, 7)

        F = d[:, 0]
        PHI = d[:, 1]
        THETA = d[:, 2]
        Fphi = d[:, 3] + d[:, 4] * 1j
        Ftheta = d[:, 5] + d[:, 6] * 1j

        self.Fphi = Fphi.reshape((nf, nphi, ntheta))
        self.Ftheta = Ftheta.reshape((nf, nphi, ntheta))
        Ttheta = THETA.reshape((nf, nphi, ntheta))
        Tphi = PHI.reshape((nf, nphi, ntheta))
        Tf = F.reshape((nf, nphi, ntheta))

        self.Fphi = self.Fphi.swapaxes(1, 2)
        self.Ftheta = self.Ftheta.swapaxes(1, 2)
        Ttheta = Ttheta.swapaxes(1, 2)
        Tphi = Tphi.swapaxes(1, 2)
        Tf = Tf.swapaxes(1, 2)

        self.fa = Tf[:, 0, 0]
        self.theta = Ttheta[0, :, 0]
        #self.phi     = Tphi[0,0,:]

        #
        # Temporaire
        #
        A1 = self.Fphi[:, 90:181, :]
        A2 = self.Fphi[:, 0:91, :]
        self.Fphi = np.concatenate((A1, A2[:, ::-1, :]), axis=2)
        A1 = self.Ftheta[:, 90:181, :]
        A2 = self.Ftheta[:, 0:91, :]
        self.Ftheta = np.concatenate((A1, A2[:, ::-1, :]), axis=2)
        self.theta = np.linspace(0, np.pi, 91)
        self.phi = np.linspace(0, 2 * np.pi, 180, endpoint=False)
        self.Nt = 91
        self.Np = 180
        self.Nf = 104
        self.evaluated = True

    def pattern(self,theta=[],phi=[],typ='s3'):
        """ return multidimensionnal radiation patterns

        Parameters
        ----------

        theta   : array
            1xNt
        phi     : array
            1xNp
        typ : string
            {s1|s2|s3}

        """

        if theta == []:
            theta = np.linspace(0,np.pi,30)
        if phi == []:
            phi = np.linspace(0,2*np.pi,60)

        Nt = len(theta)
        Np = len(phi)
        Nf = len(self.fa)
        #Th = np.kron(theta, np.ones(Np))
        #Ph = np.kron(np.ones(Nt), phi)
        if typ =='s1':
            FTh, FPh = self.Fsynth1(theta, phi,pattern=True)
        if typ =='s2':
            FTh, FPh = self.Fsynth2b(theta,phi,pattern=True)
        if typ =='s3':
            FTh, FPh = self.Fsynth3(theta, phi,pattern=True )
        #FTh = Fth.reshape(Nf, Nt, Np)
        #FPh = Fph.reshape(Nf, Nt, Np)

        return(FTh,FPh)

    def coeffshow(self,**kwargs):
        """ Display antenna coefficient

            typ : string
                'ssh' | 'vsh'
            L  : maximum level
            kf : frequency index
            vmin  : float
            vmax  : float

        """
        defaults = {'typ':'vsh',
                    'L':20,
                    'kf':46,
                    'vmin':-40,
                    'vmax':0,
                    'cmap':cm.hot_r,
                    'dB':True
                   }
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        L  = kwargs['L']
        kf = kwargs['kf']

        # calculates mode energy
        # linear and log scale
        # E : f , l , m
        if kwargs['typ']=='vsh':
            E  = self.C.energy(typ='s1')
        if kwargs['typ']=='ssh':
            E  = self.S.energy(typ='s1')
        # Aem : f,l
        # calculates energy integrated over m

        Aem = np.sum(E,axis=2)
        Aem_dB = 10*np.log10(Aem)

        # Ael : f,m
        # calculates energy integrated over l

        Ael = np.sum(E,axis=1)
        Ael_dB = 10*np.log10(Ael)


        fig, ax = plt.subplots()
        fig.set_figwidth(15)
        fig.set_figheight(10)

        if kwargs['dB']:
            im = ax.imshow(10*np.log10(E[kf,:,:]),
                           vmin = kwargs['vmin'],
                           vmax = kwargs['vmax'],
                           extent =[-L,L,L,0],
                           interpolation = 'nearest',
                           cmap = kwargs['cmap'])

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        axHistx = divider.append_axes("top", 1., pad=0.5, sharex=ax)
        axHisty = divider.append_axes("left", 1., pad=0.5, sharey=ax)
        #axHistx.bar(range(-L,L),Aem)
        #axHisty.barh(range(0,L),Ael )
        axHistx.yaxis.set_ticks(np.array([0,0.2,0.4,0.6,0.8]))
        axHisty.xaxis.set_ticks(np.array([0,0.1,0.2,0.3]))
        cbar = plt.colorbar(im, cax=cax)
        fig.tight_layout()

        plt.text(-0.02,0.6 ,'levels',
             horizontalalignment='right',
             verticalalignment='top',
             transform=ax.transAxes,
             rotation =90, fontsize= 15)

        plt.text(0.6,1.1 ,'free space',
             horizontalalignment='right',
             verticalalignment='top',
             transform=ax.transAxes,
             fontsize= 15)

        plt.text(0.55,-0.1 ,'modes',
             horizontalalignment='right'
             ,verticalalignment='top', transform=ax.transAxes, fontsize= 15)

        return fig,ax


    def errel(self,kf=-1, dsf=1, typ='s3'):
        """ calculates error between antenna pattern and reference pattern

        Parameters
        ----------

        kf  : integer
            frequency index. If k=-1 integration over all frequency
        dsf : down sampling factor
        typ :

        Returns
        -------

        errelTh : float
            relative error on :math:`F_{\\theta}`
        errelPh : float
            relative error on :math:`F_{\phi}`
        errel   : float

        Notes
        -----

        .. math::

            \epsilon_r^{\\theta} =
            \\frac{|F_{\\theta}(\\theta,\phi)-\hat{F}_{\\theta}(\\theta)(\phi)|^2}
                 {|F_{\\theta}(\\theta,\phi)|^2}

            \epsilon_r^{\phi} =
            \\frac{|F_{\phi}(\\theta,\phi)-\hat{F}_{\phi}(\\theta)(\phi)|^2}
                 {|F_{\\theta}(\\theta,\phi)|^2}


        """
        #
        # Retrieve angular bases from the down sampling factor dsf
        #
        theta = self.theta[::dsf]
        phi = self.phi[::dsf]
        Nt = len(theta)
        Np = len(phi)

        #Th = np.kron(theta, np.ones(Np))
        #Ph = np.kron(np.ones(Nt), phi)

        if typ =='s1':
            FTh, FPh = self.Fsynth1(theta, phi,pattern=True)
        if typ =='s2':
            FTh, FPh = self.Fsynth2b(theta, phi,pattern=True)
        if typ =='s3':
            FTh, FPh = self.Fsynth3(theta, phi,pattern=True)

        #FTh = Fth.reshape(self.Nf, Nt, Np)
        #FPh = Fph.reshape(self.Nf, Nt, Np)
        #
        #  Jacobian
        #
        #st    = outer(sin(theta),ones(len(phi)))
        st = np.sin(theta).reshape((len(theta), 1))
        #
        # Construct difference between reference and reconstructed
        #
        if kf<>-1:
            dTh = (FTh[kf, :, :] - self.Ftheta[kf, ::dsf, ::dsf])
            dPh = (FPh[kf, :, :] - self.Fphi[kf, ::dsf, ::dsf])
            #
            # squaring  + Jacobian
            #
            dTh2 = np.real(dTh * np.conj(dTh)) * st
            dPh2 = np.real(dPh * np.conj(dPh)) * st

            vTh2 = np.real(self.Ftheta[kf, ::dsf, ::dsf] \
                 * np.conj(self.Ftheta[kf, ::dsf, ::dsf])) * st
            vPh2 = np.real(self.Fphi[kf, ::dsf, ::dsf] \
                 * np.conj(self.Fphi[kf, ::dsf, ::dsf])) * st

            mvTh2 = np.sum(vTh2)
            mvPh2 = np.sum(vPh2)

            errTh = np.sum(dTh2)
            errPh = np.sum(dPh2)
        else:
            dTh = (FTh[:, :, :] - self.Ftheta[:, ::dsf, ::dsf])
            dPh = (FPh[:, :, :] - self.Fphi[:, ::dsf, ::dsf])
            #
            # squaring  + Jacobian
            #
            dTh2 = np.real(dTh * np.conj(dTh)) * st
            dPh2 = np.real(dPh * np.conj(dPh)) * st

            vTh2 = np.real(self.Ftheta[:, ::dsf, ::dsf] \
                 * np.conj(self.Ftheta[:, ::dsf, ::dsf])) * st
            vPh2 = np.real(self.Fphi[:, ::dsf, ::dsf] \
                 * np.conj(self.Fphi[:, ::dsf, ::dsf])) * st

            mvTh2 = np.sum(vTh2)
            mvPh2 = np.sum(vPh2)

            errTh = np.sum(dTh2)
            errPh = np.sum(dPh2)

        errelTh = (errTh / mvTh2)
        errelPh = (errPh / mvPh2)
        errel =( (errTh + errPh) / (mvTh2 + mvPh2))

        return(errelTh, errelPh, errel)


    def loadhfss(self,lfa = [], Nt=72,Np=37):
        """ load antenna from HFSS file

        Parameters
        ----------

        lfa : list of antenna file
        Nt  : int
            Number of angle theta
        Np  : int
            Number of angle phi

        Notes
        -----

        One file per frequency point

        th , ph , abs_grlz,th_absdB,th_phase,ph_absdB,ph_phase_ax_ratio

        """

        # lfa : list file antenna
        self.Nf = len(lfa)
        fGHz  = []
        lacsv = []
        Fphi = np.empty((self.Nf,self.Nt,self.Np))
        Ftheta = np.empty((self.Nf,self.Nt,self.Np))
        SqG = np.empty((self.Nf,self.Nt,self.Np))

        for i in range (len(lfa)):

            fGHz.append(eval(lfa[i].split('.csv')[0][-4]))
            lacsv.append(pd.read_csv(lfa[i],
                                     header=False,
                                     sep=',',
            names=['th','ph','abs_grlz','th_absdB','th_phase','ph_absdB','ph_phase','ax_ratio'],
                                     index_col=False))

            th=lacsv[i].th.reshape(Np,Nt)*np.pi/180.
            ph=lacsv[i].ph.reshape(Np,Nt)*np.pi/180.
            Greal = lacsv[i].abs_grlz.reshape(Np,Nt)

            th_dB = lacsv[i].th_absdB.reshape(Np,Nt)
            ph_dB = lacsv[i].ph_absdB.reshape(Np,Nt)

            th_lin = pow(10,th_dB/20.)
            ph_lin = pow(10,ph_dB/20.)

            #th_phase = lacsv[i].th_phase.reshape(72,37)*np.pi/180.
            #ph_phase = lacsv[i].ph_phase.reshape(72,37)*np.pi/180.
            #axratio=lacsv[i].ax_ratio.reshape(72,37)
            Fphi[i,:,:] = ph_lin.swapaxes(1,0)
            Ftheta[i,:,:] = th_lin.swapaxes(1,0)
            SqG[i,:,:] = Greal.swapaxes(1,0)


        self.fa = np.array(fGHz)
        #self.theta = th[0,:].reshape(Nt,1)
        #self.phi = ph[:,0].reshape(1,Np)
        self.theta = th[0,:]
        self.phi = ph[:,0]
        self.Fphi=Fphi
        self.Ftheta=Ftheta
        self.SqG=SqG



    def loadtrx(self,directory):
        """ load trx file (SATIMO Near Field Chamber raw data)

        Parameters
        ----------

        directory

        self._filename: short name of the antenna file

        the file is seek in the $BASENAME/ant directory

        .. todo:
            consider using an ini file for the header

        Trx header structure

        fmin fmax Nf  phmin   phmax   Nphi    thmin    thmax    Ntheta  #EDelay
        0     1   2   3       4       5       6        7        8       9
        1     10  121 0       6.19    72      0        3.14     37      0

        """

        _filetrx = self._filename
        _headtrx = 'header_' + _filetrx
        _headtrx = _headtrx.replace('trx', 'txt')
        headtrx = pyu.getlong(_headtrx, directory)
        filename = pyu.getlong(_filetrx, directory)
    #
    # Trx header structure
    #
    # fmin fmax Nf  phmin   phmax   Nphi    thmin    thmax    Ntheta  #EDelay
    # 0     1   2   3       4       5       6        7        8       9
    # 1     10  121 0       6.19    72      0        3.14     37      0
    #
    #
        foh = open(headtrx)
        ligh = foh.read()
        foh.close()
        fmin = eval(ligh.split()[0])
        fmax = eval(ligh.split()[1])
        nf   = eval(ligh.split()[2])
        phmin = eval(ligh.split()[3])
        phmax = eval(ligh.split()[4])
        nphi = eval(ligh.split()[5])
        thmin = eval(ligh.split()[6])
        thmax = eval(ligh.split()[7])
        ntheta = eval(ligh.split()[8])
        #
        # The electrical delay in column 9 is optional
        #
        try:
            tau = eval(ligh.split()[9])  # tau : delay (ns)
        except:
            tau = 0

        #
        # Data are stored in 7 columns
        #
        # 0  1   2     3        4       5       6
        # f  phi th    ReFph   ImFphi  ReFth    ImFth
        #
        #
        fi = open(filename)
        d = np.array(fi.read().split())
        N = len(d)
        M = N / 7
        d = d.reshape(M, 7)
        d = d.astype('float')

        f = d[:, 0]
        if f[0] == 0:
            print "error : frequency cannot be zero"
        # detect frequency unit
        # if values are above 2000 its means frequency is not expressed
        # in GHz
        #
        if (f[0] > 2000):
            f = f / 1.0e9

        phi   = d[:, 1]
        theta = d[:, 2]
        #
        # type : refers to the way the angular values are stored in the file
        # Detection of file type
        #
        # Bcp
        #    f  phi theta
        #    2    1    0
        # Natural
        #    f  phi theta
        #    2    0    1
        #
        # auto detect storage mode looping
        #
        dphi = abs(phi[0] - phi[1])
        dtheta = abs(theta[0] - theta[1])

        if (dphi == 0) & (dtheta != 0):
            typ = 'bcp'
        if (dtheta == 0) & (dphi != 0):
            typ = 'natural'

        self.typ = typ
        Fphi   = d[:, 3] + d[:, 4] * 1j
        Ftheta = d[:, 5] + d[:, 6] * 1j
        #
        # Normalization
        #
        G = np.real(Fphi * np.conj(Fphi) + Ftheta * np.conj(Ftheta))
        SqG = np.sqrt(G)
        #Fphi    = Fphi/SqG
        #Ftheta  = Ftheta/SqG
        #Fphi    = Fphi
        #Ftheta  = Ftheta
        #
        # Reshaping
        #
        if typ == 'natural':
            self.Fphi = Fphi.reshape((nf, ntheta, nphi))
            self.Ftheta = Ftheta.reshape((nf, ntheta, nphi))
            self.SqG = SqG.reshape((nf, ntheta, nphi))
            Ttheta = theta.reshape((nf, ntheta, nphi))
            Tphi = phi.reshape((nf, ntheta, nphi))
            Tf = f.reshape((nf, ntheta, nphi))
        if typ == 'bcp':
            self.Fphi = Fphi.reshape((nf, nphi, ntheta))
            self.Ftheta = Ftheta.reshape((nf, nphi, ntheta))
            self.SqG = SqG.reshape((nf, nphi, ntheta))
            Ttheta = theta.reshape((nf, nphi, ntheta))
            Tphi = phi.reshape((nf, nphi, ntheta))
            Tf = f.reshape((nf, nphi, ntheta))
        #
        # Force natural order (f,theta,phi)
        # This is not the order of the satimo nfc which is  (f,phi,theta)
        #

            self.Fphi = self.Fphi.swapaxes(1, 2)
            self.Ftheta = self.Ftheta.swapaxes(1, 2)
            self.SqG = self.SqG.swapaxes(1, 2)
            Ttheta = Ttheta.swapaxes(1, 2)
            Tphi = Tphi.swapaxes(1, 2)
            Tf = Tf.swapaxes(1, 2)

        self.fa = Tf[:, 0, 0]
        self.theta = Ttheta[0, :, 0]
        self.phi = Tphi[0, 0, :]
        #
        # check header consistency
        #
        np.testing.assert_almost_equal(self.fa[0],fmin,6)
        np.testing.assert_almost_equal(self.fa[-1],fmax,6)
        np.testing.assert_almost_equal(self.theta[0],thmin,3)
        np.testing.assert_almost_equal(self.theta[-1],thmax,3)
        np.testing.assert_almost_equal(self.phi[0],phmin,3)
        np.testing.assert_almost_equal(self.phi[-1],phmax,3)

        self.Nf = nf
        self.Nt = ntheta
        self.Np = nphi
        self.tau = tau

        self.evaluated = True

    def checkpole(self, kf=0):
        """ display the reconstructed field on pole for integrity verification

        Parameters
        ----------

        kf : int
             frequency index default 0

        """
        Ft0 = self.Ftheta[kf, 0, :]
        Fp0 = self.Fphi[kf, 0, :]
        Ftp = self.Ftheta[kf, -1, :]
        Fpp = self.Fphi[kf, -1, :]
        phi = self.phi

        Ex0 = Ft0 * np.cos(phi) - Fp0 * np.sin(phi)
        Ey0 = Ft0 * np.sin(phi) + Fp0 * np.cos(phi)
        Exp = Ftp * np.cos(phi) - Fpp * np.sin(phi)
        Eyp = Ftp * np.sin(phi) + Fpp * np.cos(phi)

        plt.subplot(4, 2, 1)
        plt.plot(phi, np.real(Ex0))
        plt.subplot(4, 2, 2)
        plt.plot(phi, np.imag(Ex0))
        plt.subplot(4, 2, 3)
        plt.plot(phi, np.real(Ey0))
        plt.subplot(4, 2, 4)
        plt.plot(phi, np.imag(Ey0))
        plt.subplot(4, 2, 5)
        plt.plot(phi, np.real(Exp))
        plt.subplot(4, 2, 6)
        plt.plot(phi, np.imag(Exp))
        plt.subplot(4, 2, 7)
        plt.plot(phi, np.real(Eyp))
        plt.subplot(4, 2, 8)
        plt.plot(phi, np.imag(Eyp))

    def info(self):
        """ gives info about antenna object


        """
        print self._filename
        print "type : ", self.typ
        if self.typ == 'mat':
            print self.DataFile
            print self.AntennaName
            print self.Date
            print self.StartTime
            print self.Notes
            print self.Serie
            print self.Run
            print "Nb theta (lat) :", self.Nt
            print "Nb phi (lon) :", self.Np
        if self.typ =='bcp':
            print "--------------------------"
            print "fmin (GHz) :", self.fa[0]
            print "fmax (GHz) :", self.fa[-1]
            print "Nf   :", self.Nf
            print "thmin (rad) :", self.theta[0]
            print "thmax (rad) :", self.theta[-1]
            print "Nth  :", self.Nt
            print "phmin (rad) :", self.phi[0]
            print "phmax (rad) :", self.phi[-1]
            print "Nph  :", self.Np
        try:
            self.C.info()
        except:
            print "No vsh coefficient calculated yet"

    def polar(self,**kwargs):
        """ antenna 2D polar plot

        Parameters
        ----------

        fGHz : frequency
        phd  : phi in degrees
        thd  : theta in degrees
        GmaxdB :  max gain to be displayed

        Returns
        -------

        fig
        ax

        Examples
        --------

        .. plot::
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> from pylayers.antprop.antenna import *
            >>> A = Antenna('defant.trx')
            >>> fig,ax = A.polar(fGHz=[2,3,4],phd=0)
            >>> fig,ax = A.polar(fGHz=[2,3,4],thd=90)

        """

        if not self.evaluated:
            self.Fsynth(pattern=True)

        dtr = np.pi/180.

        defaults = {'fGHz' : np.array([4]),
                    'dyn' : 8 ,
                    'phd' : 0,
                    'legend':True,
                    'GmaxdB':5,
                    'topos':False}

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        args = {}
        for k in kwargs:
            if k not in defaults:
                args[k] = kwargs[k]

        if 'fig' not in kwargs:
            fig = plt.figure(figsize=(8, 8))
        else:
            fig = kwargs['fig']

        if 'ax' not in kwargs:
            #ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='#d5de9c')
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True )
        else:
            ax = kwargs['ax']

        rc('grid', color='#316931', linewidth=1, linestyle='-')
        rc('xtick', labelsize=15)
        rc('ytick', labelsize=15)

        DyndB = kwargs['dyn'] * 5
        GmindB = kwargs['GmaxdB'] - DyndB
        #print "DyndB",DyndB
        #print "GmindB",GmindB

        # force square figure and square axes looks better for polar, IMO

        t1 = np.arange(5, DyndB + 5, 5)
        t2 = np.arange(GmindB + 5, kwargs['GmaxdB'] + 5, 5)


        col = ['k', 'r', 'g', 'b', 'm', 'c', 'y']
        cpt = 0
       

        if len(self.fa) > 1 :
            fstep = self.fa[1]-self.fa[0]
        else :
            fstep = np.array((abs(self.fa-kwargs['fGHz'][0])+1))
        #dtheta = self.theta[1,0]-self.theta[0,0]
        #dphi = self.phi[0,1]-self.phi[0,0]
        dtheta = self.theta[1]-self.theta[0]
        dphi = self.phi[1]-self.phi[0]
        for f in kwargs['fGHz']:
            if 'col' in kwargs:
                colr  = kwargs['col']
            else:           
                colr = col[cpt]
                
            ik = np.where(abs(self.fa-f)<fstep)[0][0]
            if 'lab' in kwargs:
                chaine  = kwargs['lab']
            else:
                chaine = 'f = %3.2f GHz' %(self.fa[ik])
            # all theta
            if 'phd' in kwargs:
                itheta = np.arange(self.Nt)
                #iphi1 = np.where(abs(self.phi[0,:]-kwargs[+'phd']*dtr)<dtheta)[0][0]
                iphi1 = np.where(abs(self.phi-kwargs['phd']*dtr)<dphi)[0][0]
                Np = self.Np

                if np.mod(Np, 2) == 0:
                    iphi2 = np.mod(iphi1 + Np / 2, Np)
                else:
                    iphi2 = np.mod(iphi1 + (Np - 1) / 2, Np)

                #   0 < theta < pi/2
                #u1 = np.where((self.theta[:,0] <= np.pi / 2) &
                #              (self.theta[:,0] >= 0))[0]
                u1 = np.where((self.theta <= np.pi / 2.) & (self.theta >= 0))[0]
                #   0:Nt-1
                u2 = np.arange(self.Nt)
                #   pi/2 < theta < pi
                #u3 = np.nonzero((self.theta[:,0] <= np.pi) & ( self.theta[:,0]
                #                                              > np.pi / 2))[0]
                u3 = np.nonzero((self.theta <= np.pi) & ( self.theta > np.pi / 2))[0]
                if self.source=='satimo':
                    r1 = -GmindB + 20 * np.log10(  self.SqG[ik, u1, iphi1]+1e-12)
                if self.source=='cst':
                    r1 = -GmindB + 20 * np.log10(  self.SqG[ik, u1, iphi1]/np.sqrt(30)+1e-12)

                #r1  = self.SqG[k,u1[0],iphi1]
                negr1 = np.nonzero(r1 < 0)
                r1[negr1[0]] = 0
                if self.source=='satimo':
                    r2 = -GmindB + 20 * np.log10( self.SqG[ik, u2, iphi2]+1e-12)
                if self.source=='cst':
                    r2 = -GmindB + 20 * np.log10(  self.SqG[ik, u2, iphi2]/np.sqrt(30)+1e-12)
                #r2  = self.SqG[k,u2,iphi2]
                negr2 = np.nonzero(r2 < 0)
                r2[negr2[0]] = 0
                if self.source=='satimo':
                    r3 = -GmindB + 20 * np.log10( self.SqG[ik, u3, iphi1]+1e-12)
                if self.source=='cst':
                    r3 = -GmindB + 20 * np.log10(  self.SqG[ik, u3, iphi1]/np.sqrt(30)+1e-12)
                #r3  = self.SqG[k,u3[0],iphi1]
                negr3 = np.nonzero(r3 < 0)
                r3[negr3[0]] = 0
                r = np.hstack((r1[::-1], r2, r3[::-1], r1[-1]))

                a1 = np.arange(0, 360, 30)
                a2 = [90, 60, 30, 0, 330, 300, 270, 240, 210, 180, 150, 120]
                rline2, rtext2 = plt.thetagrids(a1, a2)

                angle = np.linspace(0, 2 * np.pi, len(r), endpoint=True)
              #  plt.title(u'V plane $\\theta$ (degrees)')

            if 'thd' in kwargs:
                iphi = np.arange(self.Np)
                #itheta = np.where(abs(self.theta[:,0]-kwargs['thd']*dtr)<dtheta)[0][0]
                #angle = self.phi[0,iphi]
                itheta = np.where(abs(self.theta-kwargs['thd']*dtr)<dtheta)[0][0]
                #angle = self.phi[0,iphi]
                angle = self.phi[iphi]
                r = -GmindB + 20 * np.log10(self.SqG[ik, itheta, iphi])
                neg = np.nonzero(r < 0)
                r[neg] = 0
               # plt.title(u'H plane - $\phi$ degrees')
                a1 = np.arange(0, 360, 30)
                a2 = [0, 30, 60, 90, 120 , 150 , 180 , 210, 240 , 300 , 330]
                rline2, rtext2 = plt.thetagrids(a1, a2)

            ax.plot(angle, r, color=colr, lw=2, label=chaine)
            rline1, rtext1 = plt.rgrids(t1, t2)
            cpt = cpt + 1
        if kwargs['legend']:
            ax.legend()
        return(fig,ax)

    #@mlab.show
    def _show3(self,newfig = True,colorbar =True,
                    name=[],title=True,**kwargs ):
        """ show3 mayavi

        fGHz : float
            frequency
        title : bool
            display title
        colorbar :
            display colorbar
        """




        if not self.evaluated:
            self.Fsynth(pattern=True)


        x, y, z, k, scalar  = self._computemesh(**kwargs)

        if newfig:
            mlab.clf()
            f=mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
        else :
            f=mlab.gcf()


        self._mayamesh = mlab.mesh(x, y, z,scalars= scalar,resolution = 1)

        if name == []:
            f.children[-1].name = 'Antenna ' + self._filename
        else :
            f.children[-1].name = name + self._filename

        if colorbar :
            mlab.colorbar()
        if title:
            mlab.title(self._filename + ' @ ' + str(self.fa[k]) + ' GHz',height=1,size=0.5)

    def _computemesh(self,**kwargs):
        """ compute mesh from theta phi

        Returns
        -------

        (x, y, z, k)
        x , y z value in carteisan axis
        k frequency point evaluated

        """
        defaults = { 'fGHz' :[],
                     'po': np.array([0,0,0]),
                     'T' : np.eye(3),
                     'minr' : 0.1,
                     'maxr' : 1 ,
                     'tag' : 'Pat',
                     'ilog' : False,
                     'title':True,

                     'ilog':False

                     }


        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        fGHz = kwargs['fGHz']
        minr = kwargs['minr']
        maxr = kwargs['maxr']
        tag  = kwargs['tag']
        ilog = kwargs['ilog']
        po = kwargs['po']
       
        # T is an unitary matrix
        T  = kwargs['T']
        if fGHz == []:
            k = len(self.fa)/2
        else :
            k = np.where(fGHz>self.fa)[0]

        r = self.SqG[k,:,:]
        th = self.theta[:,np.newaxis]
        phi = self.phi[np.newaxis,:]

        if ilog :
            r = 10*np.log10(abs(r))
        else:
            r = abs(r)
        if r.max() != r.min():
            u = (r - r.min()) /(r.max() - r.min())
        else : u = r
        r = minr + (maxr-minr) * u

        x = r * np.sin(th) * np.cos(phi) 
        y = r * np.sin(th) * np.sin(phi) 
        z = r * np.cos(th) 

        p = np.concatenate((x[...,np.newaxis],y[...,np.newaxis],z[...,np.newaxis]),axis=2)
        #
        # antenna cs -> glogal cs
        # q : Nt x Np x 3
        q = np.einsum('ij,klj->kli',T,p)
        #
        # translation
        #
        scalar=(q[...,0]**2+q[...,1]**2+q[...,2]**2)
        
        q[...,0]=q[...,0]+po[0]
        q[...,1]=q[...,1]+po[1]
        q[...,2]=q[...,2]+po[2]

        x = q[...,0]
        y = q[...,1]
        z = q[...,2]



        return x, y, z, k, scalar

    def show3(self, k=0,po=[],T=[],typ='Gain', mode='linear', silent=False):
        """ show3 geomview

        Parameters
        ----------

        k : frequency index
        po : poition of the antenna
        T  : GCS of the antenna
        typ : string
            'Gain' | 'Ftheta' | 'Fphi'
        mode : string
            'linear'| 'not implemented'
        silent : boolean
            True    | False

        Examples
        --------

            >>> from pylayers.antprop.antenna import *
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> A = Antenna('defant.trx')
            >>> #A.show3()

        """

        if not self.evaluated:
            self.Fsynth(pattern=True)

        f = self.fa[k]

        if typ == 'Gain':
            V = self.SqG[k, :, :]
        if typ == 'Ftheta':
            V = self.Ftheta[k, :, :]
        if typ == 'Fphi':
            V = self.Fphi[k, :, :]

        if po ==[]:
            po = np.array([0, 0, 0])
        if T ==[]:
            T = np.eye(3)

        _filename = 'antbody'

        geo = geu.Geomoff(_filename)
        # geo.pattern requires the following shapes
        # theta (Ntx1)
        # phi (1xNp)
        #if len(np.shape(self.theta))==1:
        #    theta = self.theta[:,np.newaxis]
        #else:
        #    theta=self.theta
        theta = self.theta
        #if len(np.shape(self.phi))==1:
        #    phi = self.phi[np.newaxis,:]
        #else:
        #    phi=self.phi
        phi = self.phi

        geo.pattern(theta,phi,V,po=po,T=T,ilog=False,minr=0.01,maxr=0.2)
        #filename = geom_pattern(self.theta, self.phi, V, k, po, minr, maxr, typ)
        #filename = geom_pattern(self.theta, self.phi, V, k, po, minr, maxr, typ)

        if not silent:
            geo.show3()

    def plot3d(self, k=0, typ='Gain', col=True):
        """ show 3D pattern in matplotlib

        Parameters
        ----------

        k : frequency index

        typ = 'Gain'
             = 'Ftheta'
             = 'Fphi'

        if col  -> color coded plot3D
        else    -> simple plot3D

        """

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)

        if typ == 'Gain':
            V = self.SqG[k, :, :]
        if typ == 'Ftheta':
            V = self.Ftheta[k, :, :]
        if typ == 'Fphi':
            V = self.Fphi[k, :, :]

        vt = np.ones(self.Nt)
        vp = np.ones(self.Np)
        Th = np.outer(self.theta, vp)
        Ph = np.outer(vt, self.phi)

        X = abs(V) * np.cos(Ph) * np.sin(Th)
        Y = abs(V) * np.sin(Ph) * np.sin(Th)
        Z = abs(V) * np.cos(Th)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        if col:
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                            cmap=cm.hot_r,shade=True)
        else:
            ax.plot3D(np.ravel(X), np.ravel(Y), np.ravel(Z))
        plt.show()

    def pol3d(self, k=0, R=50, St=1, Sp=1, silent=False):
        """ Display polarisation diagram  in 3D

           Parameters
           ----------

           k  : int
               frequency index
           R  : float
                Radius of the sphere
           St : int
               Downsampling factor along theta
           Sp : int
               Downsampling factor along phi
           silent : Boolean
                (If True the file is created and not displayed')

           The file created is named : Polar{ifreq}.list
           it is placed in the /geom directory of the project

        """
        _filename = 'Polar' + str(10000 + k)[1:] + '.list'
        filename = pyu.getlong(_filename, pstruc['DIRGEOM'])
        fd = open(filename, "w")
        fd.write("LIST\n")

        Nt = self.Nt
        Np = self.Np
        N = 10
        plth = np.arange(0, Nt, St)
        plph = np.arange(0, Np, Sp)

        for m in plph:
            for n in plth:
                #theta = self.theta[n,0]
                theta = self.theta[n]
                #print "m,theta= :",m,theta*180/np.pi
                #phi = self.phi[0,m]
                phi = self.phi[m]
                #print "n,phi=:",n,phi*180/np.pi
                B = geu.vec_sph(theta, phi)
                p = R * np.array((np.cos(phi) * np.sin(theta),
                                  np.sin(phi) * np.sin(theta),
                                  np.cos(theta)))
                fd.write('{\n')
                geu.ellipse(fd, p, B[0, :], B[1, :], self.Ftheta[k, n, m], self.Fphi[k, n, m], N)
                fd.write('}\n')
        fd.close()
        if not silent:
            chaine = "geomview " + filename + " 2>/dev/null &"
            os.system(chaine)

    def mse(self, Fth, Fph, N=0):
        """ mean square error between original and reconstructed

        Parameters
        ----------

        Fth  : np.array
        Fph  : np.array
        N    : int


        Notes
        -----

        Calculate the relative mean square error between original pattern A.Ftheta , A.Fphi and the
        pattern given as argument of the function  Fth , Fph

        The mse is evaluated on both polarization and normalized over the energy of each
        original pattern.

        The function returns the maximum between those two errors

        N is a parameter which allows to suppress value at the pole for the calculation of the error
        if N=0 all values are kept else   N < n < Nt - N

        """

        sh = np.shape(self.Ftheta)
        Nf = sh[0]
        Nt = sh[1]
        Np = sh[2]
        # plage de theta (exclusion du pole)
        pt = np.arange(N, Nt - N, 1)

        Fthr = Fth.reshape(sh)
        Fphr = Fph.reshape(sh)

        Gr = np.real(Fphr * np.conj(Fphr) + Fthr * np.conj(Fthr))
        SqGr = np.sqrt(Gr)

        Fthr = Fthr[:, pt, :].ravel()
        Fphr = Fphr[:, pt, :].ravel()
        SqGr = SqGr[:, pt, :].ravel()

        Ftho = self.Ftheta[:, pt, :].ravel()
        Fpho = self.Fphi[:, pt, :].ravel()
        SqGo = self.SqG[:, pt, :].ravel()

        Etho = np.sqrt(np.dot(np.conj(Ftho), Ftho))
        Epho = np.sqrt(np.dot(np.conj(Fpho), Fpho))
        Eo = np.sqrt(np.dot(np.conj(Ftho), Ftho) + np.dot(np.conj(Fpho), Fpho))

        errth = Ftho - Fthr
        errph = Fpho - Fphr
        Err = np.real(np.sqrt(np.dot(np.conj(errth), errth) + np.dot(np.conj(errph), errph)))

        Errth = np.real(np.sqrt(np.dot(np.conj(errth), errth)))
        Errph = np.real(np.sqrt(np.dot(np.conj(errph), errph)))

        #Errth_rel = Errth/Etho
        #Errph_rel = Errph/Epho
        Errth_rel = Errth / Eo
        Errph_rel = Errph / Eo
        Err_rel = Err / Eo

        return Err_rel, Errth_rel, Errph_rel

    def getdelay(self,delayCandidates = np.arange(-10,10,0.001)):
        """ get electrical delay

        Parameters
        ----------

        delayCandidates : ndarray
            default np.arange(-10,10,0.001)

        Returns
        -------

        electricalDelay  : float

        Author : Troels Pedersen (Aalborg University)
                 B.Uguen

        """
        if self.evaluated:
            maxPowerInd  = np.unravel_index(np.argmax(abs(self.Ftheta)),np.shape(self.Ftheta))
            elD  = delayCandidates[np.argmax(abs(
                np.dot(self.Ftheta[:,maxPowerInd[1],maxPowerInd[2]]
                       ,np.exp(2j*np.pi*self.fa[:,np.newaxis]
                               *delayCandidates[np.newaxis,:]))))]
            #electricalDelay  = delayCandidates[np.argmax(abs(
            #    np.dot(self.Ftheta[:,maxPowerInd[1],maxPowerInd[2]]
            #        ,np.exp(2j*np.pi*freq.reshape(len(freq),1)
            #          *delayCandidates.reshape(1,len(delayCandidates))))
            #        ))]
            return(elD)
        else:
            raise Warning('Antenna has not been evaluated')

    def elec_delay(self,tau):
        r""" apply an electrical delay

        Parameters
        ----------

        tau : float
            electrical delay in nanoseconds

        Notes
        -----


         This function applies an electrical delay math::`\exp{+2 j \pi f \tau)`
         on the phase of diagram math::``F_{\theta}`` and math::`F_{\phi}`

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.antprop.antenna import *
            >>> A = Antenna('S2R2.sh3')
            >>> A.Fsynth()
            >>> tau = A.getdelay()
            >>> A.elec_delay(tau)



        """

        self.tau = self.tau+tau
        if self.evaluated:
            Ftheta = self.Ftheta
            Fphi = self.Fphi
            sh = np.shape(Ftheta)
            e = np.exp(2 * np.pi * 1j * self.fa[:,np.newaxis,np.newaxis]* tau)
            #E = np.outer(e, ones(sh[1] * sh[2]))
            #Fth = Ftheta.reshape(sh[0], sh[1] * sh[2])
            #EFth = Fth * E
            #self.Ftheta = EFth.reshape(sh[0], sh[1], sh[2])
            self.Ftheta = self.Ftheta*e
            self.Fphi = self.Fphi*e
            #Fph = Fphi.reshape(sh[0], sh[1] * sh[2])
            #EFph = Fph * E
            #self.Fphi = EFph.reshape(sh[0], sh[1], sh[2])
        else:
            raise Warning('antenna has not been evaluated')


    def Fsynth(self,theta=[],phi=[],pattern=True,typ='vsh'):
        """ Perform Antenna synthesis

        Parameters
        ----------

        theta : np.array
        phi :   np.array
        pattern : boolean
            call Antenna.Fpatt or Antenna.Fsynth3

        Notes
        -----

        The antenna pattern synthesis is done either from spherical
        harmonics coefficients or from a analytical expression of the
        radiation pattern.

        """

        if ((self.fromfile) or (self.typ=='vsh') or (self.typ=='ssh')):
            Ft,Fp = self.Fsynth3(theta,phi,pattern)
        else :
            Ft,Fp=self.Fpatt(theta,phi,pattern)
        return (Ft,Fp)    


    #def Fsynth1(self, theta, phi, k=0):
    def Fsynth1(self, theta, phi,pattern=False):
        """ calculate complex antenna pattern  from VSH Coefficients (shape 1)

        Parameters
        ----------

        theta  : ndarray (1xNdir)
        phi    : ndarray (1xNdir)
        k      : int
            frequency index

        """

        Nt = len(theta)
        Np = len(phi)

        if pattern:
            theta = np.kron(theta, np.ones(Np))
            phi = np.kron(np.ones(Nt),phi)

        nray = len(theta)

        #Br = self.C.Br.s1[k, :, :]
        #Bi = self.C.Bi.s1[k, :, :]
        #Cr = self.C.Cr.s1[k, :, :]
        #Ci = self.C.Ci.s1[k, :, :]

        Br = self.C.Br.s1[:, :, :]
        Bi = self.C.Bi.s1[:, :, :]
        Cr = self.C.Cr.s1[:, :, :]
        Ci = self.C.Ci.s1[:, :, :]

        N = self.C.Br.N1
        M = self.C.Br.M1
        #print "N,M",N,M
        #
        # The - sign is necessary to get the good reconstruction
        #     deduced from observation
        #     May be it comes from a different definition of theta in SPHEREPACK
        x = -np.cos(theta)
        Pmm1n, Pmp1n = AFLegendre3(N, M, x)
        ind = index_vsh(N, M)
        n = ind[:, 0]
        m = ind[:, 1]
        #~ V, W = VW(n, m, x, phi, Pmm1n, Pmp1n)
        V, W = VW(n, m, x, phi)
        #
        # broadcasting along frequency axis
        #
        V = np.expand_dims(V,0)
        W = np.expand_dims(V,0)
        #
        #   k : frequency axis
        #   l : coeff l
        #   m
        Fth = np.eisum('klm,kilm->ki',Br,np.real(V.T)) - \
              np.eisum('klm,kilm->ki',Bi,np.imag(V.T)) + \
              np.eisum('klm,kilm->ki',Ci,np.real(W.T)) + \
              np.eisum('klm,kilm->ki',Cr,np.imag(W.T))

        Fph = -np.eisum('klm,kilm->ki',Cr,np.real(V.T)) + \
              np.eisum('klm,kilm->ki',Ci,np.imag(V.T)) + \
              np.eisum('klm,kilm->ki',Bi,np.real(W.T)) + \
              np.eisum('klm,kilm->ki',Br,np.imag(W.T))

        #Fth = np.dot(Br, np.real(V.T)) - \
        #    np.dot(Bi, np.imag(V.T)) + \
        #    np.dot(Ci, np.real(W.T)) + \
        #    np.dot(Cr, np.imag(W.T))

        #Fph = -np.dot(Cr, np.real(V.T)) + \
        #    np.dot(Ci, np.imag(V.T)) + \
        #    np.dot(Bi, np.real(W.T)) + \
        #    np.dot(Br, np.imag(W.T))

        if pattern:
            Nf = len(self.fa)
            Fth = Fth.reshape(Nf, Nt, Np)
            Fph = Fph.reshape(Nf, Nt, Np)

        return Fth, Fph



    def Fsynth2s(self,dsf=1):
        """  pattern synthesis from shape 2 vsh coefficients

        Parameters
        ----------

        phi

        Notes
        -----

        Calculate complex antenna pattern from VSH Coefficients (shape 2)
        for the specified directions (theta,phi)
        theta and phi arrays needs to have the same size

        """
        theta = self.theta[::dsf]
        phi = self.phi[::dsf]
        Nt = len(theta)
        Np = len(phi)
        theta = np.kron(theta, np.ones(Np))
        phi = np.kron(np.ones(Nt), phi)

        Ndir = len(theta)

        Br = self.C.Br.s2 # Nf x K2
        Bi = self.C.Bi.s2 # Nf x K2
        Cr = self.C.Cr.s2 # Nf x K2
        Ci = self.C.Ci.s2 # Nf x K2

        Nf = np.shape(self.C.Br.s2)[0]
        K2 = np.shape(self.C.Br.s2)[1]

        L = self.C.Br.N2 # int
        M = self.C.Br.M2 # int

        #print "N,M",N,M
        #
        # The - sign is necessary to get the good reconstruction
        #     deduced from observation
        #     May be it comes from a different definition of theta in SPHEREPACK

        x = -np.cos(theta)

        Pmm1n, Pmp1n = AFLegendre3(L, M, x)
        ind = index_vsh(L, M)

        l = ind[:, 0]
        m = ind[:, 1]

        V, W = VW2(l, m, x, phi, Pmm1n, Pmp1n)  # K2 x Ndir

        # Fth , Fph are Nf x Ndir

        tEBr = []
        tEBi = []
        tECr = []
        tECi = []

        for k in range(K2):
            BrVr = np.dot(Br[:,k].reshape(Nf,1),
                          np.real(V.T)[k,:].reshape(1,Ndir))
            BiVi = np.dot(Bi[:,k].reshape(Nf,1),
                          np.imag(V.T)[k,:].reshape(1,Ndir))
            CiWr = np.dot(Ci[:,k].reshape(Nf,1),
                          np.real(W.T)[k,:].reshape(1,Ndir))
            CrWi = np.dot(Cr[:,k].reshape(Nf,1),
                          np.imag(W.T)[k,:].reshape(1,Ndir))

            CrVr = np.dot(Cr[:,k].reshape(Nf,1),
                          np.real(V.T)[k,:].reshape(1,Ndir))
            CiVi = np.dot(Ci[:,k].reshape(Nf,1),
                          np.imag(V.T)[k,:].reshape(1,Ndir))
            BiWr = np.dot(Bi[:,k].reshape(Nf,1),
                          np.real(W.T)[k,:].reshape(1,Ndir))
            BrWi = np.dot(Br[:,k].reshape(Nf,1),
                          np.imag(W.T)[k,:].reshape(1,Ndir))

            EBr = np.sum(BrVr*np.conj(BrVr)*np.sin(theta)) + \
                  np.sum(BrWi*np.conj(BrWi)*np.sin(theta))

            EBi = np.sum(BiVi*np.conj(BiVi)*np.sin(theta)) + \
                  np.sum(BiWr*np.conj(BiWr)*np.sin(theta))

            ECr = np.sum(CrWi*np.conj(CrWi)*np.sin(theta)) + \
                + np.sum(CrVr*np.conj(CrVr)*np.sin(theta))

            ECi = np.sum(CiWr*np.conj(CiWr)*np.sin(theta)) + \
                + np.sum(CiVi*np.conj(CiVi)*np.sin(theta))

            tEBr.append(EBr)
            tEBi.append(EBi)
            tECr.append(ECr)
            tECi.append(ECi)

        #Fth = np.dot(Br, np.real(V.T)) - np.dot(Bi, np.imag(V.T)) + \
        #      np.dot(Ci, np.real(W.T)) + np.dot(Cr, np.imag(W.T))
        #Fph = -np.dot(Cr, np.real(V.T)) + np.dot(Ci, np.imag(V.T)) + \
        #      np.dot(Bi, np.real(W.T)) + np.dot(Br, np.imag(W.T))

        return np.array(tEBr),np.array(tEBi),np.array(tECr),np.array(tECi)

    def Fsynth2b(self, theta, phi,pattern=False):
        """  pattern synthesis from shape 2 vsh coefficients

        Parameters
        ----------

        theta : 1 x Nt
        phi   : 1 x Np

        Notes
        -----

        Calculate complex antenna pattern from VSH Coefficients (shape 2)
        for the specified directions (theta,phi)
        theta and phi arrays needs to have the same size

        """

        Nt = len(theta)
        Np = len(phi)

        if pattern:
            theta = np.kron(theta, np.ones(Np))
            phi = np.kron(np.ones(Nt),phi)

        Br = self.C.Br.s2 # Nf x K2
        Bi = self.C.Bi.s2 # Nf x K2
        Cr = self.C.Cr.s2 # Nf x K2
        Ci = self.C.Ci.s2 # Nf x K2

        L = self.C.Br.N2 # int
        M = self.C.Br.M2 # int

        #print "N,M",N,M
        #
        # The - sign is necessary to get the good reconstruction
        #     deduced from observation
        #     May be it comes from a different definition of theta in SPHEREPACK

        x = -np.cos(theta)

        Pmm1n, Pmp1n = AFLegendre3(L, M, x)
        ind = index_vsh(L, M)

        l = ind[:, 0]
        m = ind[:, 1]

        V, W = VW2(l, m, x, phi, Pmm1n, Pmp1n)  # K2 x Ndir

        # Fth , Fph are Nf x Ndir
        Fth = np.dot(Br, np.real(V.T)) - np.dot(Bi, np.imag(V.T)) + \
              np.dot(Ci, np.real(W.T)) + np.dot(Cr, np.imag(W.T))

        Fph = -np.dot(Cr, np.real(V.T)) + np.dot(Ci, np.imag(V.T)) + \
              np.dot(Bi, np.real(W.T)) + np.dot(Br, np.imag(W.T))

        if pattern:
            Nf = len(self.fa)
            Fth = Fth.reshape(Nf, Nt, Np)
            Fph = Fph.reshape(Nf, Nt, Np)
        return Fth, Fph


    def Fsynth2(self, theta, phi,pattern=False, typ = 'vsh'):
        """  pattern synthesis from shape 2 vsh coeff

        Parameters
        ----------

        theta : array 1 x Nt
        phi : array 1 x Np
        pattern : boolean
            default False
        typ : string
            {vsh | ssh}


        Notes
        -----

        Calculate complex antenna pattern from VSH Coefficients (shape 2)
        for the specified directions (theta,phi)
        theta and phi arrays needs to have the same size

        """

        self.Nt = len(theta)
        self.Np = len(phi)
        self.Nf = len(self.fa)

        if typ =='vsh' :

            if pattern:
                theta = np.kron(theta, np.ones(self.Np))
                phi = np.kron(np.ones(self.Nt),phi)

            Br = self.C.Br.s2
            Bi = self.C.Bi.s2
            Cr = self.C.Cr.s2
            Ci = self.C.Ci.s2

            N = self.C.Br.N2
            M = self.C.Br.M2

            #print "N,M",N,M
            #
            # The - sign is necessary to get the good reconstruction
            #     deduced from observation
            #     May be it comes from a different definition of theta in SPHEREPACK
            x = -np.cos(theta)

            Pmm1n, Pmp1n = AFLegendre3(N, M, x)
            ind = index_vsh(N, M)

            n = ind[:, 0]
            m = ind[:, 1]

            #~ V, W = VW(n, m, x, phi, Pmm1n, Pmp1n)
            V, W = VW(n, m, x, phi)


            Fth = np.dot(Br, np.real(V.T)) - np.dot(Bi, np.imag(V.T)) + \
                np.dot(Ci, np.real(W.T)) + np.dot(Cr, np.imag(W.T))
            Fph = -np.dot(Cr, np.real(V.T)) + np.dot(Ci, np.imag(V.T)) + \
                np.dot(Bi, np.real(W.T)) + np.dot(Br, np.imag(W.T))

            if pattern:
                Fth = Fth.reshape(self.Nf, self.Nt, self.Np)
                Fph = Fph.reshape(self.Nf, self.Nt, self.Np)
        else:

            cx = self.S.Cx.s2
            cy = self.S.Cy.s2
            cz = self.S.Cz.s2

            lmax = self.S.Cx.lmax
            Y ,indx = SSHFunc(lmax, theta,phi)
            Ex = np.dot(cx,Y).reshape(self.Nf,self.Nt,self.Np)
            Ey = np.dot(cy,Y).reshape(self.Nf,self.Nt,self.Np)
            Ez = np.dot(cz,Y).reshape(self.Nf,self.Nt,self.Np)

            Fth,Fph = CartToSphere (theta, phi, Ex, Ey,Ez, bfreq = True )

        self.evaluated = True
        return Fth, Fph


    def Fsynth3(self, theta = [], phi=[], pattern=True):
        r""" synthesis of a complex antenna pattern from SH coefficients
        (vsh or ssh  in shape 3)


        Ndir is the number of directions

        Parameters
        ----------

        theta : ndarray (1xNdir if not pattern)  (1xNtheta if pattern)
        phi   : ndarray (1xNdir if not pattter)  (1xNphi if pattern)

        pattern : boolean
            if True theta and phi are reorganized for building the pattern
        typ  : 'vsh' | 'ssh' | 'hfss'

        Returns
        -------

        if pattern:
            Fth   : ndarray (Ntheta x Nphi)
            Fph   : ndarray (Ntheta x Nphi)
        else:
            Fth   : ndarray (1 x Ndir)
            Fph   : ndarray (1 x Ndir)

        See Also
        --------

        pylayers.antprop.channel._vec2scalA

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.antprop.antenna import *
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> A = Antenna('defant.vsh3')
            >>> theta = np.linspace(0,np.pi,70)
            >>> phi = np.linspace(0,2*np.pi,180)
            >>> F = A.Fsynth3(theta,phi,pattern=True)

        All Br,Cr,Bi,Ci have the same (l,m) index in order to evaluate only
        once the V,W function

        If the data comes from a cst file like the antenna used in WHERE1 D4.1
        the pattern is multiplied by $\frac{4\pi}{120\pi}=\frac{1}{\sqrt{30}$

        """

        typ = self.typ
        #self._filename.split('.')[1]
        #if typ=='satimo':
        #    coeff=1.
        #if typ=='cst':
        #    coeff=1./sqrt(30)


        assert typ in ['ssh','vsh','hfss'], "Error wrong file type"

        Nf = len(self.fa)
        if theta==[]:
            theta=np.linspace(0,np.pi,45)

        if phi == []:
            phi= np.linspace(0,2*np.pi,90)

        Nt = len(theta)
        Np = len(phi)
        self.Nt = len(theta)
        self.Np = len(phi)

        if pattern:
            #self.theta = theta[:,np.newaxis]
            #self.phi = phi[np.newaxis,:]
            self.theta = theta
            self.phi = phi
            theta = np.kron(theta, np.ones(Np))
            phi = np.kron(np.ones(Nt),phi)


        if typ =='vsh':

            nray = len(theta)

            Br  = self.C.Br.s3
            lBr = self.C.Br.ind3[:, 0]
            mBr = self.C.Br.ind3[:, 1]

            Bi  = self.C.Bi.s3
            Cr  = self.C.Cr.s3
            Ci  = self.C.Ci.s3

            L = lBr.max()
            M = mBr.max()

            # vector spherical harmonics basis functions
            V, W = VW(lBr, mBr, theta, phi)

            Fth = np.dot(Br, np.real(V.T)) - \
                  np.dot(Bi, np.imag(V.T)) + \
                  np.dot(Ci, np.real(W.T)) + \
                  np.dot(Cr, np.imag(W.T))

            Fph = -np.dot(Cr, np.real(V.T)) + \
                   np.dot(Ci, np.imag(V.T)) + \
                   np.dot(Bi, np.real(W.T)) + \
                   np.dot(Br, np.imag(W.T))

            if pattern:

                Fth = Fth.reshape(Nf, Nt, Np)
                Fph = Fph.reshape(Nf, Nt, Np)


        if typ == 'ssh':
            cx = self.S.Cx.s3
            cy = self.S.Cy.s3
            cz = self.S.Cz.s3

            lmax = self.S.Cx.lmax
            Y ,indx = SSHFunc2(lmax, theta,phi)

            #k = self.S.Cx.k2[:,0]
            # same k for x y and z
            k = self.S.Cx.k2
            if pattern :
                Ex = np.dot(cx,Y[k])
                Ey = np.dot(cy,Y[k])
                Ez = np.dot(cz,Y[k])
                Fth,Fph = CartToSphere(theta, phi, Ex, Ey,Ez, bfreq = True, pattern = True )
                Fth = Fth.reshape(Nf,Nt,Np)
                Fph = Fph.reshape(Nf,Nt,Np)


            else:
                

                Ex = np.dot(cx,Y[k])
                Ey = np.dot(cy,Y[k])
                Ez = np.dot(cz,Y[k])
                Fth,Fph = CartToSphere (theta, phi, Ex, Ey,Ez, bfreq = True, pattern = False)

            self.Fphi = Fph
            self.Ftheta = Fth
            G = np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth))
            self.SqG = np.sqrt(G)


        if pattern :
            self.Fphi = Fph
            self.Ftheta = Fth
            G = np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth))
            self.SqG = np.sqrt(G)

        self.evaluated = True

        if typ == 'hfss':
            scipy.interpolate.griddata()

            Fth = self.Ftheta
            Fph = self.Fphi
        # TODO create 2 different functions for pattern and not pattern
        if not pattern:
            return Fth, Fph
        else:
            return None,None

    def gain(self,th,ph,dB=True):
        """
        Parameters
        ----------

        th : theta angle of arrival|departure in local coordinate system 
        ph : phi in 

        """
        if not isinstance(th,np.ndarray):
            th=np.array([th])
            ph=np.array([ph]) 
        Ft,Fp = self.Fsynth(theta=th,phi=ph,pattern=False)
        G   = Ft*np.conj(Ft)+Fp*np.conj(Fp)
        if dB:
            return 10*np.log10(G)
        else:
            return G 


    def movie_vsh(self, mode='linear'):
        """ animates vector spherical coeff w.r.t frequency

        Parameters
        ----------
        mode : string
            'linear' |
        """

        Brmin = abs(self.C.Br[:, 0:20, 0:20]).min()
        Brmax = abs(self.C.Br[:, 0:20, 0:20]).max()
        Bimin = abs(self.C.Bi[:, 0:20, 0:20]).min()
        Bimax = abs(self.C.Bi[:, 0:20, 0:20]).max()

        Crmin = abs(self.C.Cr[:, 0:20, 0:20]).min()
        Crmax = abs(self.C.Cr[:, 0:20, 0:20]).max()
        Cimin = abs(self.C.Ci[:, 0:20, 0:20]).min()
        Cimax = abs(self.C.Ci[:, 0:20, 0:20]).max()

        print Brmin, Brmax, Bimin, Bimax, Crmin, Crmax, Cimin, Cimax

        for k in range(self.nf):
            plt.figure()
            stf = ' f=' + str(self.fa[k]) + ' GHz'
            subplot(221)
            pcolor(abs(self.C.Br.s1[k, 0:20, 0:20]),
                   vmin=Brmin, vmax=Brmax, edgecolors='k')
            #xlabel('m',fontsize=12)
            ylabel('n', fontsize=12)
            title('$|Br_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            subplot(222)
            pcolor(abs(self.C.Bi.s1[k, 0:20, 0:20]),
                   vmin=Bimin, vmax=Bimax, edgecolors='k')
            #xlabel('m',fontsize=12)
            ylabel('n', fontsize=12)
            title('$|Bi_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            subplot(223)
            pcolor(abs(self.C.Cr.s1[k, 0:20, 0:20]),
                   vmin=Crmin, vmax=Crmax, edgecolors='k')
            xlabel('m', fontsize=12)
            #ylabel('n',fontsize=12)
            title('$|Cr_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            subplot(224)
            pcolor(abs(self.C.Ci.s1[k, 0:20, 0:20]),
                   vmin=Cimin, vmax=Cimax, edgecolors='k')
            xlabel('m', fontsize=12)
            #ylabel('n',fontsize=12)
            title('$|Ci_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            filename = str('%03d' % k) + '.png'
            savefig(filename, dpi=100)
            clf()

        command = ('mencoder',
                   'mf://*.png',
                   '-mf',
                   'type=png:w=800:h=600:fps=1',
                   '-ovc',
                   'lavc',
                   '-lavcopts',
                   'vcodec=mpeg4',
                   '-oac',
                   'copy',
                   '-o',
                   'vshcoeff.avi')
        subprocess.check_call(command)

    def minsh3(self, emax=0.05):
        """ creates vsh3 with significant coeff until given relative reconstruction error

        Parameters
        ----------

        emax : float
            error default 0.05

        Summary
        -------

        Create antenna's vsh3 file which only contains
        the significant vsh coefficients in shape 3,
        in order to obtain a reconstruction maximal error =  emax

        This function requires a reading of .trx file before being executed

        """

        #th = np.kron(self.theta, np.ones(self.Np))
        #ph = np.kron(np.ones(self.Nt), self.phi)

        Fth3, Fph3 = self.Fsynth3(self.theta, self.phi,pattern=True)
        Err = self.mse(Fth3, Fph3, 0)

        Enc = self.C.ens3()
        n = len(Enc)
        pos = 0

        while (pos < n) & (Err[0] < emax):

            Emin = Enc[pos]
            d = self.C.drag3(Emin)
            Fth3, Fph3 = self.Fsynth3(self.theta, self.phi,pattern=True)
            Err = self.mse(Fth3, Fph3, 0)

            if Err[0] >= emax:
                i = d[0][0]
                i3 = d[1][0]
                self.C.put3(i, i3)
                Fth3, Fph3 = self.Fsynth3(self.theta,self.phi,pattern=True)
                Err = self.mse(Fth3, Fph3, 0)

            pos = pos + 1

    def savevsh3(self):
        """ save antenna in vsh3 format

        Create a .vsh3 antenna file


        """

        # create vsh3 file

        _filevsh3 = os.path.splitext(self._filename)[0]+'.vsh3'
        filevsh3 = pyu.getlong(_filevsh3, pstruc['DIRANT'])

        #filevsh3 = pyu.getlong(self._filename,'ant')

        if os.path.isfile(filevsh3):
            print filevsh3, ' already exist'
        else:
            print 'create ', filevsh3, ' file'

            coeff = {}
            coeff['fmin'] = self.fa[0]
            coeff['fmax'] = self.fa[-1]
            coeff['Br.ind'] = self.C.Br.ind3
            coeff['Bi.ind'] = self.C.Bi.ind3
            coeff['Cr.ind'] = self.C.Cr.ind3
            coeff['Ci.ind'] = self.C.Ci.ind3
            coeff['Br.k'] = self.C.Br.k2
            coeff['Bi.k'] = self.C.Bi.k2
            coeff['Cr.k'] = self.C.Cr.k2
            coeff['Ci.k'] = self.C.Ci.k2
            coeff['Br.s3'] = self.C.Br.s3
            coeff['Bi.s3'] = self.C.Bi.s3
            coeff['Cr.s3'] = self.C.Cr.s3
            coeff['Ci.s3'] = self.C.Ci.s3
            io.savemat(filevsh3, coeff, appendmat=False)

    def savesh2(self):
        """ save coeff in  .sh2 antenna file

        """

        # create sh2 file
        typ = self._filename.split('.')[1]
        self.typ = typ

        _filesh2 = self._filename.replace('.'+ typ, '.sh2')
        filesh2 = pyu.getlong(_filesh2, pstruc['DIRANT'])
        if os.path.isfile(filesh2):
            print filesh2, ' already exist'
        else:
            print 'create ', filesh2, ' file'
            coeff = {}
            if type(self.fa) == float:
                fmin = self.fa
                fmax = self.fa
            else:
                fmin = self.fa[0]
                fmax = self.fa[-1]
            coeff['fmin'] = fmin
            coeff['fmax'] = fmax
            

            coeff['Cx.ind'] = self.S.Cx.ind2
            coeff['Cy.ind'] = self.S.Cy.ind2
            coeff['Cz.ind'] = self.S.Cz.ind2
            coeff['Cx.lmax']= self.S.Cx.lmax
            coeff['Cy.lmax']= self.S.Cy.lmax
            coeff['Cz.lmax']= self.S.Cz.lmax

            coeff['Cx.s2'] = self.S.Cx.s2
            coeff['Cy.s2'] = self.S.Cy.s2
            coeff['Cz.s2'] = self.S.Cz.s2
            io.savemat(filesh2, coeff, appendmat=False)

    def savesh3(self):
        """ save antenna in sh3 format

        Create a .sh3 antenna file


        """
        # create sh3 file
        typ = self._filename.split('.')[1]
        self.typ = typ
        _filesh3 = self._filename.replace('.'+ typ, '.sh3')
        filesh3 = pyu.getlong(_filesh3, pstruc['DIRANT'])
        if os.path.isfile(filesh3):
            print filesh3, ' already exist'


        else:
            print 'create ', filesh3, ' file'

            coeff = {}
            if type(self.fa) == float:
                fmin = self.fa
                fmax = self.fa
            else:
                fmin = self.fa[0]
                fmax = self.fa[-1]
            coeff['fmin'] = fmin
            coeff['fmax'] = fmax
            coeff['Cx.ind'] = self.S.Cx.ind3
            coeff['Cy.ind'] = self.S.Cy.ind3
            coeff['Cz.ind'] = self.S.Cz.ind3

            coeff['Cx.k'] = self.S.Cx.k2
            coeff['Cy.k'] = self.S.Cy.k2
            coeff['Cz.k'] = self.S.Cz.k2


            coeff['Cx.lmax']= self.S.Cx.lmax
            coeff['Cy.lmax']= self.S.Cy.lmax
            coeff['Cz.lmax']= self.S.Cz.lmax

            coeff['Cx.s3'] = self.S.Cx.s3
            coeff['Cy.s3'] = self.S.Cy.s3
            coeff['Cz.s3'] = self.S.Cz.s3


            io.savemat(filesh3, coeff, appendmat=False)

    def loadvsh3(self):
        """ Load antenna's vsh3 file

            vsh3 file contains a thesholded version of vsh coefficients in shape 3

        """

        _filevsh3 = self._filename
        filevsh3 = pyu.getlong(_filevsh3, pstruc['DIRANT'])
        self.evaluated = False

        if os.path.isfile(filevsh3):
            coeff = io.loadmat(filevsh3, appendmat=False)
            #
            # This test is to fix a problem with 2 different
            # behavior of io.loadmat
            #
            if type(coeff['fmin']) == float:
                fmin = coeff['fmin']
                fmax = coeff['fmax']
            else:
                fmin = coeff['fmin'][0][0]
                fmax = coeff['fmax'][0][0]
            # .. Warning
            # Warning modification takes only one dimension for k
            # if the .vsh3 format evolve it may not work anymore
            #
            Br = VCoeff('s3', fmin, fmax, coeff['Br.s3'],
                         coeff['Br.ind'], coeff['Br.k'][0])
            Bi = VCoeff('s3', fmin, fmax, coeff['Bi.s3'],
                         coeff['Bi.ind'], coeff['Bi.k'][0])
            Cr = VCoeff('s3', fmin, fmax, coeff['Cr.s3'],
                         coeff['Cr.ind'], coeff['Cr.k'][0])
            Ci = VCoeff('s3', fmin, fmax, coeff['Ci.s3'],
                         coeff['Ci.ind'], coeff['Ci.k'][0])
            self.C = VSHCoeff(Br, Bi, Cr, Ci)
            Nf = np.shape(Br.s3)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filevsh3, ' does not exist'

    def loadsh3(self):
        """ Load antenna's sh3 file

            sh3 file contains a thesholded version of ssh coefficients in shape 3

        """
        _filesh3 = self._filename.split('.')[0]+'.sh3'
        filesh3 = pyu.getlong(_filesh3, pstruc['DIRANT'])
        self.evaluated = False

        if os.path.isfile(filesh3):
            coeff = io.loadmat(filesh3, appendmat=False)

            #
            # This test is to fix a problem with 2 different
            # behavior of io.loadmat
            #
            if type(coeff['fmin']) == float:
                fmin = coeff['fmin']
                fmax = coeff['fmax']
            else:
                fmin = coeff['fmin'][0][0]
                fmax = coeff['fmax'][0][0]
            # .. Warning
            # Warning modification takes only one dimension for k
            # if the .sh3 format evolve it may not work anymore
            #


            if type(coeff['Cx.lmax']) == float:
                lmax = coeff['Cx.lmax']
            else:
                lmax = coeff['Cx.lmax'][0][0]

            Cx = SCoeff(typ = 's3',
                        fmin = fmin ,
                        fmax = fmax ,
                        lmax = lmax,
                        data = coeff['Cx.s3'],
                        ind =  coeff['Cx.ind'],
                        k =  np.squeeze(coeff['Cx.k']))


            Cy = SCoeff(typ= 's3',
                        fmin = fmin ,
                        fmax = fmax ,
                        lmax = lmax,
                        data = coeff['Cy.s3'],
                        ind =  coeff['Cy.ind'],
                        k =  np.squeeze(coeff['Cy.k']))



            Cz = SCoeff(typ = 's3',
                        fmin = fmin ,
                        fmax = fmax ,
                        data = coeff['Cz.s3'],
                        lmax = lmax,
                        ind =  coeff['Cz.ind'],
                        k =  np.squeeze(coeff['Cz.k']))


            if not 'S' in self.__dict__.keys():
                self.S = SSHCoeff(Cx, Cy,Cz)
            else:
                self.S.sets3(Cx,Cy,Cz)

            self.Nf = np.shape(Cx.s3)[0]
            self.fa = np.linspace(fmin, fmax, self.Nf)
        else:
            print _filesh3, ' does not exist'

    def savevsh2(self, filename = ''):
        """ save coeff in  a .vsh2 antenna file

        """

        # create vsh2 file
        if filename == '':
            _filevsh2 = self._filename.replace('.trx', '.vsh2')

        _filevsh2  = filename
        filevsh2 = pyu.getlong(_filevsh2, pstruc['DIRANT'])

        if os.path.isfile(filevsh2):
            print filevsh2, ' already exist'
        else:
            print 'create ', filevsh2, ' file'

            coeff = {}
            coeff['fmin'] = self.fa[0]
            coeff['fmax'] = self.fa[-1]
            coeff['Br.ind'] = self.C.Br.ind2
            coeff['Bi.ind'] = self.C.Bi.ind2
            coeff['Cr.ind'] = self.C.Cr.ind2
            coeff['Ci.ind'] = self.C.Ci.ind2
            coeff['Br.s2'] = self.C.Br.s2
            coeff['Bi.s2'] = self.C.Bi.s2
            coeff['Cr.s2'] = self.C.Cr.s2
            coeff['Ci.s2'] = self.C.Ci.s2

            io.savemat(filevsh2, coeff, appendmat=False)



    def loadsh2(self):
        """ load  spherical harmonics coefficient in shape  2

        """
        _filesh2 = self._filename.split('.')[0]+'.sh2'
        filesh2 = pyu.getlong(_filesh2, pstruc['DIRANT'])

        if os.path.isfile(filesh2):
            coeff = io.loadmat(filesh2, appendmat=False)
            #
            # This test is to fix a problem with 2 different
            # behavior of io.loadmat
            #
            if type(coeff['fmin']) == float:
                fmin = coeff['fmin']
                fmax = coeff['fmax']
            else:
                fmin = coeff['fmin'][0][0]
                fmax = coeff['fmax'][0][0]


            if type(coeff['Cx.lmax']) == float:
                lmax = coeff['Cx.lmax']
            else:
                lmax = coeff['Cx.lmax'][0][0]


            Cx = SCoeff(typ='s2',
                        fmin=fmin,
                        fmax=fmax,
                        lmax = lmax,
                        data=coeff['Cx.s2'],
                        ind=coeff['Cx.ind'])

            Cy = SCoeff(typ='s2',
                        fmin=fmin,
                        fmax=fmax,
                        lmax = lmax,
                        data=coeff['Cy.s2'],
                        ind=coeff['Cy.ind'])
            Cz = SCoeff(typ='s2',
                        fmin=fmin,
                        fmax=fmax,
                        lmax = lmax,
                        data=coeff['Cz.s2'],
                        ind=coeff['Cz.ind'])


            self.S = SSHCoeff(Cx, Cy,Cz)
            Nf = np.shape(Cx.s2)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filesh2, ' does not exist'



    def loadvsh2(self):
        """ load antenna from .vsh2 file format

        Load antenna's vsh2 file which only contains
        the vsh coefficients in shape 2

        """

        _filevsh2 = self._filename
        filevsh2 = pyu.getlong(_filevsh2, pstruc['DIRANT'])

        if os.path.isfile(filevsh2):
            coeff = io.loadmat(filevsh2, appendmat=False)
            #
            # This test is to fix a problem with 2 different
            # behavior of io.loadmat
            #
            if type(coeff['fmin']) == float:
                fmin = coeff['fmin']
                fmax = coeff['fmax']
            else:
                fmin = coeff['fmin'][0][0]
                fmax = coeff['fmax'][0][0]

            Br = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Br.s2'], ind=coeff['Br.ind'])
            Bi = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Bi.s2'], ind=coeff['Bi.ind'])
            Cr = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Cr.s2'], ind=coeff['Cr.ind'])
            Ci = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Ci.s2'], ind=coeff['Ci.ind'])
            self.C = VSHCoeff(Br, Bi, Cr, Ci)
            Nf = np.shape(Br.s2)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filevsh2, ' does not exist'

    def loadvsh3_old(self):
        """ Load antenna vsh coefficients in shape 3

        """

        _filevsh3 = self._filename
        filevsh3 = getlong(_filevsh3, pstruc['DIRANT'])
        fmin = 2.
        fmax = 8.
        if os.path.isfile(filevsh3):
            coeff = io.loadmat(filevsh3, appendmat=False)
            Br = VCoeff('s3', fmin, fmax, coeff['Br.s3'],
                         coeff['Br.ind'], coeff['Br.k'])
            Bi = VCoeff('s3', fmin, fmax, coeff['Bi.s3'],
                         coeff['Bi.ind'], coeff['Bi.k'])
            Cr = VCoeff('s3', fmin, fmax, coeff['Cr.s3'],
                         coeff['Cr.ind'], coeff['Cr.k'])
            Ci = VCoeff('s3', fmin, fmax, coeff['Ci.s3'],
                         coeff['Ci.ind'], coeff['Ci.k'])
            self.C = VSHCoeff(Br, Bi, Cr, Ci)
            self.fa = np.linspace(fmin, fmax, 121)
        else:
            print _filevsh3, ' does not exist'

    def pol2cart(self, ith):
        """ converts FTheta, FPhi to Fx,Fy,Fz for theta=ith

        Parameters
        ----------
        ith : theta index

        Returns
        -------

        Fx
        Fy
        Fz

        See Also
        --------

        cart2pol

        """
        Fth = self.Ftheta[:, ith, :]
        Fph = self.Fphi[:, ith, :]
        th = self.theta[ith]
        ph = self.phi

        Fx = Fth * np.cos(th) * np.cos(ph) - Fph * np.sin(ph)
        Fy = Fth * np.cos(th) * np.sin(ph) + Fph * np.cos(ph)
        Fz = (-1) * Fth * np.sin(th)

        return(Fx, Fy, Fz)

    def cart2pol(self, Fx, Fy, Fz, ith):
        """ converts Fx,Fy,Fz to Ftheta, Fphi for theta=ith

        Parameters
        ----------

        Fx  : np.array
        Fy  : np.array
        Fz  : np.array
        ith : theta index

        See Also
        --------

        pol2cart

        """
        th = self.theta[ith]
        ph = self.phi

        Fth = Fx * np.cos(th) * np.cos(ph) + Fy * np.cos(th) * np.sin(ph) - Fz * np.sin(th)
        Fph = -Fx * np.sin(ph) + Fy * np.cos(th)

        SqG = np.sqrt(np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth)))

        self.SqG[:, ith, :] = SqG
        self.Ftheta[:, ith, :] = Fth
        self.Fphi[:, ith, :] = Fph

def forcesympol(A):
    """ plot VSH transform vsh basis in 3D plot

    Parameters
    ----------

    n,m   : integer values (m<=n)
    theta : ndarray
    phi   : ndarray
    sf    : boolean
        if sf : plotted figures are saved in a *.png file
        else  : plotted figures aren't saved

    Examples
    --------

    .. plot::
        :include-source:

        >>> from pylayers.antprop.antenna import *
        >>> import  matplotlib.pyplot as plt
        >>> import numpy as np
        >>> n=5
        >>> m=3
        >>> theta = np.linspace(0,np.pi,30)
        >>> phi   = np.linspace(0,2*np.pi,60)
        >>> plotVW(n,m,theta,phi)

    """
    # calculate v and w
    if m <= n:
        theta[np.where(theta == np.pi / 2)[0]] = np.pi / 2 + \
            1e-10  # .. todo :: not clean
        x = -np.cos(theta)
        Pmm1n, Pmp1n = AFLegendre(n, m, x)

        t1 = np.sqrt((n + m) * (n - m + 1))
        t2 = np.sqrt((n - m) * (n + m + 1))
        y1 = t1 * Pmm1n[:, m, n] - t2 * Pmp1n[:, m, n]
        y2 = t1 * Pmm1n[:, m, n] + t2 * Pmp1n[:, m, n]

        Ephi = np.exp(1j * m * phi)
        cphi = np.cos(m * phi)
        if m == 0:
            sphi = 1e-10
        else:
            sphi = np.sin(m * phi)

        ny = len(y1)
        ne = len(Ephi)
        vy = np.ones(ny)
        ve = np.ones(ne)
        Y1 = np.outer(y1, ve)
        Y2 = np.outer(y2, ve)
        EPh = np.outer(vy, Ephi)

        const = (-1.0) ** n / (2 * np.sqrt(n * (n + 1)))
        V = const * Y1 * EPh
        #V[np.isinf(V)|isnan(V)]=0
        Vcos = cphi * V
        Vsin = sphi * V

        if m == 0:
            #W=np.zeros((len(theta),len(phi)))
            W = np.ones((len(theta), len(phi))) * 1e-10
        else:
            Waux = Y2 * EPh
            x1 = 1.0 / x
            W = np.outer(x1, const) * Waux
        Wcos = cphi * W
        Wsin = sphi * W

        # plot V and W
        Ntheta = np.size(theta)
        vt = np.ones(Ntheta)
        Nphi = np.size(phi)
        vp = np.ones(Nphi)
        Phi = np.outer(vt, phi)
        Theta = np.outer(theta, vp)

        #figdirV='/home/rburghel/Bureau/bases_decomposition_VW/base_V_Vsin_Vcos/'
        figdirV = './'
        ext1 = '.pdf'
        ext2 = '.eps'
        ext3 = '.png'

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(V) * np.cos(Phi) * np.sin(Theta)
        Y = abs(V) * np.sin(Phi) * np.sin(Theta)
        Z = abs(V) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirV + 'V' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Vcos) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Vcos) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Vcos) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])

        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirV + 'Vcos' + str(n) + str(m) + '.jpg'
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Vsin) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Vsin) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Vsin) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirV + 'Vsin' + str(n) + str(m) + '.jpg'
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        #figdirW='/home/rburghel/Bureau/bases_decomposition_VW/base_W_Wsin_Wcos/'
        figdirW = './'

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(W) * np.cos(Phi) * np.sin(Theta)
        Y = abs(W) * np.sin(Phi) * np.sin(Theta)
        Z = abs(W) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirW + 'W' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Wcos) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Wcos) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Wcos) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirW + 'Wcos' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Wsin) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Wsin) * np.sin(Phi) * np.sin(Theta)

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Wsin) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Wsin) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Wsin) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirW + 'Wsin' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        plt.show()

    else:
        print "Error: m>n!!!"

def compdiag(k, A, th, ph, Fthr, Fphr, typ='modulus', lang='english', fontsize=18):
    """ makes comparison between original pattern and reconstructed pattern

    Parameters
    ----------

    k : frequency index
    A : Antenna

    ph : phi base      (1 x Np)
    th : theta base    (1 x Nt)

    Fthr : Fth output of Fsynth Nf x (Ntheta*Tphi)
    Fphr : Fth output of Fsynth Nf x (Ntheta*Tphi)

    lang = 'french'
         = 'english'

    """

    Nf = np.shape(Fthr)[0]

    #Fthr = Fthr.reshape(Nf,len(th),len(ph))
    #Fphr = Fphr.reshape(Nf,len(th),len(ph))

    plt.figure()
    rc('text', usetex=True)
    Ftho = A.Ftheta
    Fpho = A.Fphi

    # limites module Fthr, Ftho, Fphr, Fpho
    maxTr = abs(Fthr[k, :, :]).max()
    maxTo = abs(Ftho[k, :, :]).max()
    MmT = max(maxTr, maxTo)

    minTr = abs(Fthr[k, :, :]).min()
    minTo = abs(Ftho[k, :, :]).min()
    mmT = min(minTr, minTo)

    maxPr = abs(Fphr[k, :, :]).max()
    maxPo = abs(Fpho[k, :, :]).max()
    MmP = max(maxPr, maxPo)

    minPr = abs(Fphr[k, :, :]).min()
    minPo = abs(Fpho[k, :, :]).min()
    mmP = min(minPr, minPo)

    # limites real Fthr, Ftho, Fphr, Fpho
    maxTrr = np.real(Fthr[k, :, :]).max()
    maxTor = np.real(Ftho[k, :, :]).max()
    MrT = max(maxTrr, maxTor)

    minTrr = np.real(Fthr[k, :, :]).min()
    minTor = np.real(Ftho[k, :, :]).min()
    mrT = min(minTrr, minTor)

    maxPrr = np.real(Fphr[k, :, :]).max()
    maxPor = np.real(Fpho[k, :, :]).max()
    MrP = max(maxPrr, maxPor)

    minPrr = np.real(Fphr[k, :, :]).min()
    minPor = np.real(Fpho[k, :, :]).min()
    mrP = min(minPrr, minPor)

    # limites real Fthr, Ftho, Fphr, Fpho
    maxTri = np.imag(Fthr[k, :, :]).max()
    maxToi = np.imag(Ftho[k, :, :]).max()
    MiT = max(maxTri, maxToi)

    minTri = np.imag(Fthr[k, :, :]).min()
    minToi = np.imag(Ftho[k, :, :]).min()
    miT = min(minTri, minToi)

    maxPri = np.imag(Fphr[k, :, :]).max()
    maxPoi = np.imag(Fpho[k, :, :]).max()
    MiP = max(maxPri, maxPoi)

    minPri = np.imag(Fphr[k, :, :]).min()
    minPoi = np.imag(Fpho[k, :, :]).min()
    miP = min(minPri, minPoi)

    # limithes arg Fth,Fph
    maxATr = np.angle(Fthr[k, :, :]).max()
    maxATo = np.angle(Ftho[k, :, :]).max()
    maT = max(maxATr, maxATo)
    minATr = np.angle(Fthr[k, :, :]).min()
    minATo = np.angle(Ftho[k, :, :]).min()
    maT0 = min(minATr, minATo)

    maxAPr = np.angle(Fphr[k, :, :]).max()
    maxAPo = np.angle(Fpho[k, :, :]).max()
    maP = max(maxAPr, maxAPo)
    minAPr = np.angle(Fphr[k, :, :]).min()
    minAPo = np.angle(Fpho[k, :, :]).min()
    maP0 = min(minAPr, minAPo)

    ax = plt.axes([0, 0, 360, 180])
    rtd = 180 / np.pi

    plt.subplot(221)
    if typ == 'modulus':
    #
    #cmap=cm.jet
        #pcolor(A.phi*rtd,A.theta*rtd,abs(Ftho[k,:,:]),vmin=0,vmax=mmT)
            #
    #cmap= gray
    #pcolor(A.phi*rtd,A.theta*rtd,abs(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
            #
    #cmap=cm.hot
        plt.pcolor(A.phi * rtd, A.theta * rtd, abs(Ftho[k, :, :]),
                   cmap=cm.hot_r, vmin=mmT, vmax=MmT)
        plt.title(r'$|F_{\theta}|$ original', fontsize=fontsize)

    if typ == 'real':
        #pcolor(A.phi*rtd,A.theta*rtd,real(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.real(Ftho[k, :, :]),
                   cmap=cm.hot_r, vmin=mrT, vmax=MrT)
        title(r'Re ($F_{\theta}$) original', fontsize=fontsize)
    if typ == 'imag':
        #pcolor(A.phi*rtd,A.theta*rtd,imag(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
        pcolor(A.phi * rtd, A.theta * rtd, np.imag(Ftho[k, :, :]),
               cmap=cm.hot_r, vmin=miT, vmax=MiT)
        title(r'Im ($F_{\theta}$) original', fontsize=fontsize)
    if typ == 'phase':
        #pcolor(A.phi*rtd,A.theta*rtd,angle(Ftho[k,:,:]),cmap=cm.gray_r,vmin=maT0,vmax=maT)
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Ftho[k, :, :]),
                   cmap=cm.hot_r, vmin=maT0, vmax=maT)
        if lang == 'french':
            plt.title(r'Arg ($F_{\theta}$) original', fontsize=fontsize)
        else:
            plt.title(r'Ang ($F_{\theta}$) original', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.ylabel(r'$\theta$ (deg)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)

    plt.subplot(222)
    if typ == 'modulus':
        plt.pcolor(A.phi * rtd, A.theta * rtd, abs(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=mmP, vmax=MmP)
        plt.title('$|F_{\phi}|$ original', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.real(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=mrP, vmax=MrP)
        plt.title('Re ($F_{\phi}$) original', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.imag(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=miP, vmax=MiP)
        plt.title('Im ($F_{\phi}$) original', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=maP0, vmax=maP)
        if lang == 'french':
            plt.title('Arg ($F_{\phi}$) original', fontsize=fontsize)
        else:
            plt.title('Ang ($F_{\phi}$) original', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)

    plt.subplot(223)
    if typ == 'modulus':
        plt.pcolor(ph * rtd, th * rtd, abs(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=mmT, vmax=MmT)
        if lang == 'french':
            plt.title(r'$|F_{\theta}|$ reconstruit', fontsize=fontsize)
        else:
            plt.title(r'$|F_{\theta}|$ reconstructed', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(ph * rtd, th * rtd, np.real(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=mrT, vmax=MrT)
        if lang == 'french':
            title(r'Re ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            title(r'Re ($F_{\theta}$) reconstructed', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(ph * rtd, th * rtd, np.imag(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=miT, vmax=MiT)
        if lang == 'french':
            plt.title(r'Im ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            plt.title(r'Im ($F_{\theta}$) reconstructed', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=maT0, vmax=maT)
        if lang == 'french':
            plt.title(r'Arg ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            plt.title(r'Ang ($F_{\theta}$) reconstructed', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.xlabel(r'$\phi$ (deg)', fontsize=fontsize)
    plt.ylabel(r'$\theta$ (deg)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)

    plt.subplot(224)
    if typ == 'modulus':
        plt.pcolor(ph * rtd, th * rtd, abs(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=mmP, vmax=MmP)
        if lang == 'french':
            plt.title('$|F_{\phi}|$ reconstruit', fontsize=fontsize)
        else:
            plt.title('$|F_{\phi}|$ reconstructed', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(ph * rtd, th * rtd, np.real(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=mrP, vmax=MrP)
        if lang == 'french':
            plt.title('Re ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Re ($F_{\phi}$) reconstructed', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(ph * rtd, th * rtd, np.imag(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=miP, vmax=MiP)
        if lang == 'french':
            plt.title('Im ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Im ($F_{\phi}$) reconstructed', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=maP0, vmax=maP)
        if lang == 'french':
            plt.title('Arg ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Ang ($F_{\phi}$) reconstructed', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.xlabel(r'$\phi$ (deg)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)

def show3D(F, theta, phi, k, col=True):
    """ show 3D matplotlib diagram

    Parameters
    ----------

    F      : ndarray (Nf,Nt,Np)
    theta  : ndarray (1xNt)
        angle
    phi    : ndarray (1xNp)
        angle
    theta  : ndarray (Nt)
    k      : int
        frequency index
    col    : boolean
        if col  -> color coded plot3D
        if col == False -> simple plot3D

    Examples
    --------

    .. plot::
        :include-source:

        >>> import matplotlib.pyplot as plt
        >>> from pylayers.antprop.antenna import *
        >>> A = Antenna('defant.vsh3')
        >>> A.Fsynth()


    Warnings
    --------

        len(theta) must be equal with shape(F)[1]
        len(phi) must be equal with shape(F)[2]

    """

    nth = len(theta)
    nph = len(phi)

    if k >= np.shape(F)[0]:
        print 'Error: frequency index k not in F defined interval'
    if nth != np.shape(F)[1]:
        print 'Error: shape mistmatch between theta and F'

    if nph != np.shape(F)[2]:
        print 'Error: shape mistmatch between phi and F'

    fig = plt.figure()
    ax = axes3d.Axes3D(fig)

    V = F[k, :, :]

    vt = np.ones(nth)
    vp = np.ones(nph)
    Th = np.outer(theta, vp)
    Ph = np.outer(vt, phi)

    X = abs(V) * np.cos(Ph) * np.sin(Th)
    Y = abs(V) * np.sin(Ph) * np.sin(Th)
    Z = abs(V) * np.cos(Th)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    if (col):
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
    else:
        ax.plot3D(np.ravel(X), np.ravel(Y), np.ravel(Z))


if (__name__ == "__main__"):
    plt.ion()
    doctest.testmod()
