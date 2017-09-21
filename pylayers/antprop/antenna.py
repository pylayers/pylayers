#-*- coding:Utf-8 -*-
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

Antenna derives from Pattern 

Examples
--------

.. plot::
    :include-source:

    >>> import matplotlib.pyplot as plt
    >>> from pylayers.antprop.antenna import *
    >>> A = Antenna()
    >>> fig,ax = A.plotG(fGHz=[2,3,4],plan='theta',angdeg=0)

Pattern Class
-------------

.. autosummary::
    :toctree: generated/

    Pattern.eval
    Pattern.gain
    Pattern.radF

Pattern Functions
=================

    Pattern.__pOmni
    Pattern.__pGauss
    Pattern.__p3gpp
    Pattern.__p3gpp

Pattern from SH coeff
=====================

    Pattern.__pvsh3
    Pattern.__psh3


Antenna Class
-------------

.. autosummary::
    :toctree: generated/


Utility Functions
=================

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
===================

.. autosummary::
    :toctree: generated/

    Antenna.Fsynth
    Antenna.Fsynth1
    Antenna.Fsynth2s
    Antenna.Fsynth2b
    Antenna.Fsynth2
    Antenna.Fsynth3

Visualization functions
=======================

.. autosummary::
    :toctree: generated/

    Antenna.pattern
    Antenna.plotG
    Antenna._show3
    Antenna.show3
    Antenna.plot3d
    Antenna.pol3d
    Antenna.load_trx
    Antenna.movie_vsh

Loading and Saving
==================

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

Miscellaneous  functions
========================


.. autosummary::
    :toctree: generated/

    forcesympol
    compdiag
    show3D


"""
#from __future__ import print_function
import doctest
import os
import glob
import re
import pdb
import sys
if sys.version_info.major==2:
    import PIL.Image as Image
    import mayavi.mlab as mlab
else:
    import image
import numpy as np
import scipy.linalg as la
from scipy import io
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
from pylayers.util.project import *
from pylayers.antprop.spharm import *
try:
    from pylayers.antprop.antvsh import vsh 
except:
    pass
from pylayers.antprop.antssh import ssh,SSHFunc2, SSHFunc, SSHCoeff, CartToSphere
from pylayers.antprop.coeffModel import *
from matplotlib import rc
from matplotlib import cm # colormaps
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from scipy.special import sici , fresnel
import pandas as pd

import matplotlib.pylab as plt

class Pattern(PyLayers):
    """ Class Pattern

    MetaClass of Antenna

    A pattern is evaluated with the 3 np.array parameters

    theta
    phi
    fGHz

    This class implements pattern methods.
    The name of a pattern method starts by p.
    Each pattern method has a unique dictionnary argument 'param'

    If self.grid dimensions are
        Nt x Np x Nf
    else:
        Ndir x Nf

    """
    def __init__(self):
        PyLayers.__init__(self)

    def __repr__(self):
        st = ''
        st = st + 'Antenna type : ' + self.typ +'\n'
        st = st+'------------------------\n'
        if 'param' in self.__dict__:
            for k in self.param:
                st = st + ' ' + k + ' : ' + str(self.param[k])+'\n'
        return (st)

    def eval(self,**kwargs):
        """  evaluate pattern functions


        Parameters
        ----------

        th: list 
            []
        ph: list 
            []
        pt : np.array (3,N)
        pr : np.array (3,N)
        azoffset : int (0) 
        Rfloor:bool
            if true add gain value to reflected ray on the floor. 
            values are append at the end of sqG.
        fGHz:list 
            []
        nth: int 
            90
        nph: int 
            181
        first: boolean 
            True if first call (to define self.param)
        grid:  boolean 
            True for pattern mode, False for Ray Tracing mode
        th0 : float 
            theta initial value
        th1 : float 
            theta finale value
        ph0 : float 
            phi initial value
        ph1 : float 
            phi final value


        Examples
        --------



        >>> from pylayers.antprop.aarray import *
        >>> A0=Antenna('Omni',param={'pol':'t','GmaxdB':0})
        >>> A1=Antenna('Gauss')
        >>> A2=Antenna('3gpp')
        >>> A3=ULArray()
        >>> A0.eval()
        >>> A1.eval()
        >>> A2.eval()
        >>> #A3.eval()

        """
        defaults = {'Rfloor':False,
                    'nth':90,
                    'nph':181,
                    'grid':True,
                    'th0':0,
                    'th1':np.pi,
                    'ph0':0,
                    'ph1':2*np.pi,
                    'azoffset':0,
                    'inplace':True
                   }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]


        if 'fGHz' not in kwargs:
            if 'fGHz' not in self.__dict__:
                self.fGHz = np.array([2.4])
        else:
            if type(kwargs['fGHz'])==np.ndarray:
                self.fGHz = kwargs['fGHz']
            else:
                self.fGHz = np.array([kwargs['fGHz']])
        self.nf = len(self.fGHz)
        self.grid = kwargs['grid']
        #
        # if th and ph are empty 
        #    if pt and pr are empty 
        #          calculates from th0,th1,nth 
        #                           ph0,phi,nph
        #    else
        #          calculates from points coordinates pt and pr
        # else
        #     take specified values
        if ('th' not in kwargs) and ('ph' not in kwargs):
            if ('pt' not in kwargs) and ('pr' not in kwargs):
                self.theta = np.linspace(kwargs['th0'],kwargs['th1'],kwargs['nth'])
                self.phi = np.linspace(kwargs['ph0'],kwargs['ph1'],kwargs['nph'],endpoint=False)
                self.grid = True
                self.full_evaluated = True
            else:
                si = kwargs['pr']-kwargs['pt']
                ssi = np.sqrt(np.sum(si*si,axis=0))
                sn = si/ssi[None,:]
                self.theta = np.arccos(sn[2,:])
                self.phi = np.mod(np.arctan2(sn[1,:],sn[0,:])+kwargs['azoffset'],2*np.pi)
                self.grid = False
                self.full_evaluated = True
                if kwargs['Rfloor']:
                    dR = np.sqrt(ssi**2 + (kwargs['pr'][2,:] + kwargs['pt'][2,:])**2)  #  reflexion length
                    thetaR = np.arccos((kwargs['pr'][2,:] + kwargs['pt'][2,:]) / dR)
                    self.theta = np.hstack([self.theta,thetaR])
                    self.phi = np.hstack([self.phi,self.phi])
                    
        else :
            assert(len(kwargs['th'])==len(kwargs['ph']))
            self.theta = kwargs['th']
            self.phi = kwargs['ph']
            self.full_evaluated = False
        if self.typ=='azel':
            self.theta=np.linspace(-np.pi,np.pi,360)
            self.phi=np.linspace(-np.pi,np.pi,360)
            self.full_evaluated = False
        self.nth = len(self.theta)
        self.nph = len(self.phi)

        #
        # evaluation of the specific Pattern__p function 
        #
        Ft,Fp = eval('self._Pattern__p'+self.typ)(param=self.param)
        if kwargs['inplace']:
            self.Ft = Ft
            self.Fp = Fp
            self.evaluated = True
            self.gain()
        else:
            return Ft,Fp
        
    def vsh(self,threshold=-1):
        if self.evaluated:
            vsh(self)
            self.C.s1tos2()
            self.C.s2tos3(threshold=threshold)

    def ssh(self,L=89,dsf=1):
        if self.evaluated:
            ssh(self,L,dsf)

    def __pOmni(self,**kwargs):
        """  omnidirectional pattern

        Parameters
        ----------

        param : dict 
            dictionnary of parameters
            + pol : string 
                't'| 'p'
            + GmaxdB : float
                0


        self.grid is used for switching between : 

          if True   angular grid : nth x nph x nf
          if False  direction    : ndir x nf

        """
        defaults = { 'param' : { 'pol' : 't', 'GmaxdB': 0 } }

        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.param = kwargs['param']
        self.GmaxdB  = self.param['GmaxdB']
        self.pol  = self.param['pol']
        G    = pow(10.,self.GmaxdB/10.) # linear gain
        if self.grid:
            # Nth x Nphx Nf
            self.sqG  = np.array(np.sqrt(G))*np.ones(len(self.fGHz))[None,None,:]
            self.evaluated = True
        else:
            # Nd x Nf
            self.sqG =  np.array(np.sqrt(G))*np.ones(len(self.fGHz))[None,:]
        Ft,Fp = self.radF()
        return Ft,Fp


    def __paperture(self,**kwargs):
        """ Aperture Pattern 

        Aperture in the (x,y) plane. Main lobe in theta=0 direction

        polar indicates the orientation of the Electric field either 'x' or 'y'
       
        See theoretical background in : 

        http://www.ece.rutgers.edu/~orfanidi/ewa/ch18.pdf

        Parameters
        ----------

        HPBW_x_deg : float 
            Half Power Beamwidth (degrees)
        HPBW_y_deg : float 
            Half Power Beamwidth (degrees)


        """
        defaults = {'param': {'HPBW_x_deg':40,
                              'HPBW_y_deg':10,
                              'Gfactor':27000,
                              'fcGHz': 27.5,
                              'polar':'x',
                              'window':'rect'
                             }}
        
        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.param = kwargs['param']

        deg_to_rad = np.pi/180.
        ld_c = 0.3/self.param['fcGHz']
        ld = 0.3/self.fGHz
        Dx = 0.886*ld_c/(self.param['HPBW_x_deg']*deg_to_rad)
        Dy = 0.886*ld_c/(self.param['HPBW_y_deg']*deg_to_rad)
        Dx_n = Dx/ld
        Dy_n = Dy/ld
        if self.grid: 
            # Nth x Nph x Nf
            theta = self.theta[:,None,None]
            phi = self.phi[None,:,None]
        else:
            # Ndir x Nf 
            theta = self.theta[:,None]
            phi = self.phi[:,None]
        
        vx = Dx_n[...,:]*np.sin(theta)*np.cos(phi) # 18.1.4
        vy = Dy_n[...,:]*np.sin(theta)*np.sin(phi) # 18.1.4

        F_nor = ((1+np.cos(theta))/2.)*np.abs(np.sinc(vx)*np.sinc(vy))
        HPBW_x = (0.886*ld/Dx)/deg_to_rad
        HPBW_y = (0.886*ld/Dy)/deg_to_rad
        Gmax = self.param['Gfactor']/(HPBW_x*HPBW_y)
        F  = np.sqrt(Gmax[...,:])*F_nor # Ndir x Nf 

        # Handling repartition on both vector components
        # enforce E.y = 0 
        if self.param['polar']=='x':
            Ft = F/np.sqrt(1+(np.cos(theta)*np.sin(phi)/np.cos(phi))**2)
            Fp = (-np.cos(theta)*np.sin(phi)/np.cos(phi))*Ft
            nan_bool = np.isnan(Fp)
            Fp[nan_bool] = F[nan_bool] 
        # enforce E.x = 0 
        if self.param['polar']=='y':
            Ft = F/np.sqrt(1+(np.cos(theta)*np.cos(phi)/np.sin(phi))**2)
            Fp = (np.cos(theta)*np.cos(phi)/np.sin(phi))*Ft
            nan_bool = np.isnan(Fp)
            Fp[nan_bool] = F[nan_bool] 
        # enforce E.x = 0 
        #
        # This is experimental 
        # How to apply the 2D windowing properly ?
        #
#        if self.param['window']!='rect':
#            Nt = self.Fp.shape[0] 
#            Np = self.Fp.shape[1] 
#            Wp = np.fft.ifftshift(np.hamming(Nt)[:,None]*np.ones(Np)[None,:])[:,:,None]
#            Wt = np.fft.ifftshift(np.ones(Nt)[:,None]*np.hamming(Np)[None,:])[:,:,None] 
#            Wu = np.fft.ifftshift(np.ones(Nt)[:,None]*np.ones(Np)[None,:])[:,:,None] 
#            Wi = np.fft.ifftshift(np.hamming(Nt)[:,None]*np.hamming(Np)[None,:])[:,:,None] 
#            W = np.fft.fftshift(np.hamming(Nt)[:,None]*np.hamming(Np)[None,:])[:,:,None] 
#            # Fp : t x p x f   ou r x f 
#            # Ft : t x p x f   ou r x f 
#
#            Kp = np.fft.ifft2(self.Fp,axes=(0,1))
#            Kt = np.fft.ifft2(self.Ft,axes=(0,1))
#            
#            self.Fp = np.fft.fft2(Kp*Wt,axes=(0,1))
#            self.Ft = np.fft.fft2(Kt*Wp,axes=(0,1))

        return Ft,Fp

    def __paperture2(self,**kwargs):
        """ Aperture Pattern 

        Aperture in the (x,y) plane. Main lobe in theta=0 direction

        polar indicates the orientation of the Electric field either 'x' or 'y'
       
        See theoretical background in : 

        http://www.ece.rutgers.edu/~orfanidi/ewa/ch18.pdf

        Parameters
        ----------

        HPBW_x_deg : float 
            Half Power Beamwidth (degrees)
        HPBW_y_deg : float 
            Half Power Beamwidth (degrees)


        """
        defaults = {'param': {'HPBW_a_deg':40,
                              'HPBW_b_deg':10,
                              'Gfactor':27000,
                              'fcGHz': 27.5,
                              'polar':'x',
                              'window':'rect'
                             }}
        
        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.param = kwargs['param']

        deg_to_rad = np.pi/180.
        ld_c = 0.3/self.param['fcGHz']
        ld = 0.3/self.fGHz
        a = 1.189*ld_c/(self.param['HPBW_a_deg']*deg_to_rad)
        b = 0.886*ld_c/(self.param['HPBW_b_deg']*deg_to_rad)
        a_n = a/ld
        b_n = b/ld
        if self.grid: 
            # Nth x Nph x Nf
            theta = self.theta[:,None,None]
            phi = self.phi[None,:,None]
        else:
            # Ndir x Nf 
            theta = self.theta[:,None]
            phi = self.phi[:,None]
        
        vx = a_n[...,:]*np.sin(theta)*np.cos(phi) # 18.1.4
        vy = b_n[...,:]*np.sin(theta)*np.sin(phi) # 18.1.4

        #F_nor = ((1+np.cos(theta))/2.)*np.abs(np.sinc(vx)*np.sinc(vy))
        F_nor = (1+np.cos(theta))/2*(np.cos(np.pi*vx)/(1-4*vx**2))*np.sinc(vy) # 18.1.3 + suppression rear radiation

        HPBW_a = (1.189*ld/a)/deg_to_rad
        HPBW_b = (0.886*ld/b)/deg_to_rad
        Gmax = self.param['Gfactor']/(HPBW_a*HPBW_b)
        F  = np.sqrt(Gmax[...,:])*F_nor # Ndir x Nf 

        # Handling repartition on both vector components
        # enforce E.y = 0 
        if self.param['polar']=='x':
            Ft = F/np.sqrt(1+(np.cos(theta)*np.sin(phi)/np.cos(phi))**2)
            Fp = (-np.cos(theta)*np.sin(phi)/np.cos(phi))*Ft
            nan_bool = np.isnan(Fp)
            Fp[nan_bool] = F[nan_bool] 
        # enforce E.x = 0 
        if self.param['polar']=='y':
            Ft = F/np.sqrt(1+(np.cos(theta)*np.cos(phi)/np.sin(phi))**2)
            Fp = (np.cos(theta)*np.cos(phi)/np.sin(phi))*Ft
            nan_bool = np.isnan(Fp)
            Fp[nan_bool] = F[nan_bool] 
        # enforce E.x = 0 
        #
        # This is experimeintal 
        # How to apply the 2D windowing properly ?
        #
#        if self.param['window']!='rect':
#            Nt = self.Fp.shape[0] 
#            Np = self.Fp.shape[1] 
#            Wp = np.fft.ifftshift(np.hamming(Nt)[:,None]*np.ones(Np)[None,:])[:,:,None]
#            Wt = np.fft.ifftshift(np.ones(Nt)[:,None]*np.hamming(Np)[None,:])[:,:,None] 
#            Wu = np.fft.ifftshift(np.ones(Nt)[:,None]*np.ones(Np)[None,:])[:,:,None] 
#            Wi = np.fft.ifftshift(np.hamming(Nt)[:,None]*np.hamming(Np)[None,:])[:,:,None] 
#            W = np.fft.fftshift(np.hamming(Nt)[:,None]*np.hamming(Np)[None,:])[:,:,None] 
#            # Fp : t x p x f   ou r x f 
#            # Ft : t x p x f   ou r x f 
#
#            Kp = np.fft.ifft2(self.Fp,axes=(0,1))
#            Kt = np.fft.ifft2(self.Ft,axes=(0,1))
#            
#            self.Fp = np.fft.fft2(Kp*Wt,axes=(0,1))
#            self.Ft = np.fft.fft2(Kt*Wp,axes=(0,1))

        return Ft,Fp

    def __phplanesectoralhorn(self,**kwargs):
        """ H plane sectoral horn 


        Parameters
        ----------
        
        rho1 : float 
            sector radius (meter)
        a1 : float
            aperture dimension along x (greatest value in meters)
        b1 : float 
            aperture dimension along y (greatest value in meters) 

        Notes
        -----

        Maximum gain in theta =0 
        Polarized along y axis (Jx=0,Jz=0)  

        """

        defaults = {'param': {'rho1':0.198,
                              'a1':0.088,  # aperture dimension along x
                              'b1':0.0126, # aperture dimension along y 
                              'fcGHz':28,
                              'GcmaxdB':19,
                              'Nx':20,
                              'Ny':20}}

        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.param = kwargs['param']
        #H-plane antenna
        rho1            = self.param['rho1']
        a1              = self.param['a1']
        b1              = self.param['b1']
        Nx              = self.param['Nx']
        Ny              = self.param['Ny']
        fcGHz           = self.param['fcGHz']
        GcmaxdB         = self.param['GcmaxdB']
        assert(a1>b1), "a1 should be greater than b1 (see fig 13.1O(a) Balanis"

        lbda   = 0.3/self.fGHz
        k      = 2*np.pi/lbda
        eta0    = np.sqrt(4*np.pi*1e-7/8.85429e-12)

        if self.grid:
            # X,Y aperture points (t,p,x,y,f)
            X = np.arange(-a1/2,a1/2,a1/(Nx-1))[None,None,:,None,None]
            Y = np.arange(-b1/2,b1/2,b1/(Ny-1))[None,None,None,:,None]
            # angular domain (theta,phi)
            Theta= self.theta[:,None,None,None,None]
            Phi = self.phi[None,:,None,None,None]
        else:
            # X,Y aperture points (r,x,y,f)
            X = np.arange(-a1/2,a1/2,a1/(Nx-1))[None,:,None,None]
            Y = np.arange(-b1/2,b1/2,b1/(Ny-1))[None,None,:,None]
            # angular domain (theta,phi)
            Theta= self.theta[:,None,None,None]
            Phi= self.phi[:,None,None,None]


        #% Aperture field Ea:
        # Ea is an approximation of the aperture field:
        # (from: C. A. Balanis, Antenna Theoy: Analysis and Design. New York
        # Wiley, 1982. ... Section 13.3.1 )

        Ea = np.cos(X*np.pi/a1)*np.exp(-.5*1j*k*((X**2)/(rho1)+(Y**2)/(rho1)))
        Jy = -Ea/eta0
        Mx = Ea

        # cosine direction
        ctsp = np.cos(Theta)*np.sin(Phi)
        cp = np.cos(Phi)
        ctcp = np.cos(Theta)*np.cos(Phi)
        sp = np.sin(Phi) 
        stcp = np.sin(Theta)*np.cos(Phi)
        stsp = np.sin(Theta)*np.sin(Phi)
        # N & L
        ejkrrp = np.exp(1j*k*( X*stcp + Y*stsp))        # exp(jk (r.r'))
        if self.grid:
            N_theta  = np.einsum('tpnmf->tpf',Jy*ctsp*ejkrrp) # 12-12 a assuming Jx,Jz=0
            N_phi    = np.einsum('tpnmf->tpf',Jy*cp*ejkrrp)   # 12-12 b ""
            L_theta  = np.einsum('tpnmf->tpf',Mx*ctcp*ejkrrp) # 12-12 c assuming My,Mz=0 
            L_phi    = np.einsum('tpnmf->tpf',-Mx*sp*ejkrrp)  # 12-12 d ""
        else:
            N_theta  = np.einsum('rnmf->rf',Jy*ctsp*ejkrrp) # 12-12 a assuming Jx,Jz=0
            N_phi    = np.einsum('rnmf->rf',Jy*cp*ejkrrp)   # 12-12 b ""
            L_theta  = np.einsum('rnmf->rf',Mx*ctcp*ejkrrp) # 12-12 c assuming My,Mz=0 
            L_phi    = np.einsum('rnmf->rf',-Mx*sp*ejkrrp)  # 12-12 d ""


        # Far-Field
        Ft  = -L_phi  - eta0*N_theta  # 12-10b p 661
        Fp  = L_theta - eta0*N_phi   # 12-10c p 661 
        G = Ft*np.conj(Ft)+Fp*np.conj(Fp)
        if self.grid:
            # Umax : ,f 
            self.Umax = G.max(axis=(0,1))
            Ft = Ft/np.sqrt(self.Umax[None,None,:])
            Fp = Fp/np.sqrt(self.Umax[None,None,:])
            # centered frequency range
            fcc = np.abs(self.fGHz-fcGHz)
            idxc = np.where(fcc==np.min(fcc))[0][0] 
            # Gain @ center frequency 
            #G = _gain(Ft[:,:,idxc],Fp[:,:,idxc])
            G = _gain(Ft,Fp)
            # effective half power beamwidth
            self.ehpbw, self.hpster  = _hpbw(G,self.theta,self.phi)
            self.Gfactor = 10**(GcmaxdB/10.)*self.ehpbw[idxc]
            Gmax = self.Gfactor/self.ehpbw
            Ft = np.sqrt(Gmax[None,None,:])*Ft
            Fp = np.sqrt(Gmax[None,None,:])*Fp
        else:
            ##
            ## Ft (r x f ) 
            ## Fp (r x f ) 
            ##
            Ft = Ft/np.sqrt(self.Umax[None,:])
            Fp = Fp/np.sqrt(self.Umax[None,:])
            Gmax = self.Gfactor/self.ehpbw
            Ft = np.sqrt(Gmax[None,:])*Ft
            Fp = np.sqrt(Gmax[None,:])*Fp

        return Ft,Fp

    def __phorn(self,**kwargs):
        """ Horn antenna 


        http://www.ece.rutgers.edu/~orfanidi/ewa/ch18.pdf (18.2) 

        Parameters
        ----------

            Half Power Beamwidth (degrees)


        """
        defaults = {'param': {'sigma_a':1.2593,
                              'sigma_b':1.0246,
                              'A_wl':16,
                              'B_wl':3,
                              'fcGHz':28.,
                              'polar':'x'
                             }}
        
        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.param = kwargs['param']

        deg_to_rad = np.pi/180.
        ld_c = 0.3/self.param['fcGHz']
        ld = 0.3/self.fGHz
        A_wl = kwargs['param']['A_wl']
        B_wl = kwargs['param']['B_wl']

        A = A_wl*ld_c
        B = B_wl*ld_c
        sigma_a = kwargs['param']['sigma_a']
        sigma_b = kwargs['param']['sigma_b']
        #b = kwargs['param']['b']
        #Ra = (A/(A-a))*RA
        #Rb = (B/(B-b))*RB
        #La = np.sqrt(Ra**2+A**2/4)
        #Lb = np.sqrt(Rb**2+B**2/4)
        #alpha = np.arctan(A/(2*Ra))
        #beta = np.arctan(B/(2*Rb))
        #Delta_a = A**2/(8*Ra)
        #Delta_b = B**2/(8*Rb)
        #sigma_a = A/np.sqrt((2*ld*Ra))
        #sigma_b = B/np.sqrt((2*ld*Rb))
        A_n = A/ld
        B_n = B/ld



        if self.grid: 
            # Nth x Nph x Nf
            theta = self.theta[:,None,None]
            phi = self.phi[None,:,None]
        else:
            # Ndir x Nf 
            theta = self.theta[:,None]
            phi = self.phi[:,None]
        
        vx = A_n[...,:]*np.sin(theta)*np.cos(phi) # 18.3.4
        vy = B_n[...,:]*np.sin(theta)*np.sin(phi) # 18.3.4

        F = ((1+np.cos(theta))/2.)*(F1(vx,sigma_a)*F0(vy,sigma_b))
        normF = np.abs(F1(0,sigma_a)*F0(0,sigma_b))**2  
        F_nor = F/np.sqrt(normF)
        efficiency = 0.125*normF # 18.4.3
        Gmax = efficiency*4*np.pi*A*B/ld**2
        F  = np.sqrt(Gmax[...,:])*F_nor # Ndir x Nf 

        # Handling repatition on both vector components
        # enforce E.y = 0 
        if self.param['polar']=='x':
            Ft = F/np.sqrt(1+(np.cos(theta)*np.sin(phi)/np.cos(phi))**2)
            Fp = (-np.cos(theta)*np.sin(phi)/np.cos(phi))*Ft
            nan_bool = np.isnan(Fp)
            Fp[nan_bool] = F[nan_bool] 
        # enforce E.x = 0 
        if self.param['polar']=='y':
            Ft = F/np.sqrt(1+(np.cos(theta)*np.cos(phi)/np.sin(phi))**2)
            Fp = (np.cos(theta)*np.cos(phi)/np.sin(phi))*Ft
            nan_bool = np.isnan(Fp)
            Fp[nan_bool] = F[nan_bool] 

        return Ft,Fp 

    def __pazel(self,**kwargs):
        """ Azimuth Elevation pattern from file

        Parameters
        ----------

        filename : ANT filename

        """



        defaults = {'param': {'filename' : '',
                              'pol':'V'}}

        f = open(kwargs['param']['filename'])
        Gthetaphi = f.readlines()
        f.close()
        Gthetaphi = np.array(Gthetaphi).astype('float')
        Gaz = Gthetaphi[360:]
        Gel = Gthetaphi[:360]

        sqGazlin = np.sqrt(pow(10,Gaz/10.))
        sqGellin = np.sqrt(pow(10,Gel/10.))

        if self.grid :
            # Nth x Nph x Nf
            if kwargs['param']['pol']=='V':
                Ft = np.ones((360,360,1))
                Fp = np.zeros((360,360,1))
                #Ft[180,:] = sqGazlin[:,None]
                #Ft[:,180] = sqGellin[:,None]
                Ft = sqGazlin[None,:,None]*sqGellin[:,None,None]
            if kwargs['param']['pol']=='H':
                Fp = np.ones((360,360,1))
                Ft = np.zeros((360,360,1))
                Fp = sqGazlin[None,:,None]*sqGellin[:,None,None]
                #self.Fp[180,:]= sqGazlin[:,None]
                #self.Fp[:,180]= sqGellin[:,None]
            if kwargs['param']['pol']=='45':
                Fp = np.ones((360,360,1))
                Ft = np.ones((360,360,1))
                # Azimuth
                Ft = (1/sqrt(2))*sqGazlin[None,:,None]*sqGellin[:,None,None]
                Fp = (1/sqrt(2))*sqGazlin[None,:,None]*sqGellin[:,None,None]
                #self.Fp[180,:]= sqGazlin[:,None]
                #self.Fp[180,:]= (1/sqrt(2))*sqGazlin[:,None]
                #Ft[180,:]= (1/sqrt(2))*sqGazlin[:,None]
                # Elevation
                #self.Fp[:,180]= (1/sqrt(2))*sqGellin[:,None]
                #Ft[:,180]= (1/sqrt(2))*sqGellin[:,None]

            #Ft = sqGthlin[:,None,None]
            #self.Fp = sqGphlin[None,:,None]
            # Ft = self.sqGmax * ( np.exp(-2.76*argth[:,None,None]) * np.exp(-2.76*argphi[None,:,None]) )
            # self.Fp = self.sqGmax * ( np.exp(-2.76*argth[:,None,None]) * np.exp(-2.76*argphi[None,:,None]) )
            self.evaluated = True
        else:
            pass
            # #
            # #  Nd x Nf
            # #
            # Ft = self.sqGmax * ( np.exp(-2.76*argth) * np.exp(-2.76*argphi) )
            # Fp = self.sqGmax * ( np.exp(-2.76*argth) * np.exp(-2.76*argphi) )
            # # add frequency axis (Ndir x Nf)
            # Ft = np.dot(Ft[:,None],np.ones(len(self.fGHz))[None,:])
            # self.Fp = np.dot(Fp[:,None],np.ones(len(self.fGHz))[None,:])
        return Ft,Fp


    def __pGauss(self,**kwargs):
        """ Gauss pattern

        Parameters
        ----------

        p0 : phi main lobe (0-2pi)
        p3 : 3dB aperture angle
        t0 : theta main lobe (0-pi)
        t3 : 3dB aperture angle

        TODO : finish implementation of polar

        """
        defaults = {'param':{'p0' : 0,
                    't0' : np.pi/2,
                    'p3' : np.pi/6,
                    't3' : np.pi/6,
                    'pol':'th'
                   }}

        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.typ='Gauss'
        self.param = kwargs['param']

        p0 = self.param['p0']
        t0 = self.param['t0']
        p3 = self.param['p3']
        t3 = self.param['t3']
        pol = self.param['pol']

        self.Gmax = 16/(t3*p3)
        self.GdB = 10*np.log10(self.Gmax)
        self.sqGmax = np.sqrt(self.Gmax)

        argth = ((self.theta-t0)**2)/t3

        e1 = np.mod(self.phi-p0,2*np.pi)
        e2 = np.mod(p0-self.phi,2*np.pi)

        e = np.array(map(lambda x: min(x[0],x[1]),zip(e1,e2)))
        argphi = (e**2)/p3
        Nf = len(self.fGHz)

        if self.grid :
            Nt = len(self.theta)
            Np = len(self.phi)
            # Nth x Nph x Nf
            # Ft = self.sqGmax * ( np.exp(-2.76*argth[:,None,None]) * np.exp(-2.76*argphi[None,:,None]) )
            # self.Fp = self.sqGmax * ( np.exp(-2.76*argth[:,None,None]) * np.exp(-2.76*argphi[None,:,None]) )
            if pol=='th':
                Ft = self.sqGmax * ( np.exp(-2.76*argth[:,None,None]) * np.exp(-2.76*argphi[None,:,None]) *np.ones(len(self.fGHz))[None,None,:])
                Fp = np.zeros((Nt,Np,Nf))
            if pol=='ph':
                Ft = np.zeros((Nt,Np,Nf))
                Fp = self.sqGmax * ( np.exp(-2.76*argth[:,None,None]) * np.exp(-2.76*argphi[None,:,None]) *np.ones(len(self.fGHz))[None,None,:])
        else:
            #
            #  Nd x Nf
            #
            Nd = len(self.theta)
            assert(len(self.phi)==Nd)
            if pol=='th':
                Ft = self.sqGmax * ( np.exp(-2.76*argth) * np.exp(-2.76*argphi) )
                Fp = np.zeros(Nd)
            if pol=='ph':
                Ft = np.zeros(Nd)
                Fp = self.sqGmax * ( np.exp(-2.76*argth) * np.exp(-2.76*argphi) )
            # add frequency axis (Ndir x Nf)
            Ft = np.dot(Ft[:,None],np.ones(len(self.fGHz))[None,:])
            Fp = np.dot(Fp[:,None],np.ones(len(self.fGHz))[None,:])
        return Ft,Fp

    def __p3gpp(self,**kwargs):
        """ 3GPP pattern

        Parameters
        ----------

        thtilt : theta tilt antenna
        hpbwv  : half power beamwidth v
        hpbwh  : half power beamwidth h
        sllv   : side lobe level
        fbrh   : front back ratio
        gm     :
        pol    : h | v | c


        if pattern
            Ft  nth x nphi x nf
            Fp  nth x nphi x nf
        else
            Ft  ndir x nf (==nth, ==nph)
            Fp  ndir x nf (==nth, ==nph)

        """
        defaults = {'param' : {'thtilt':0,  # antenna tilt
                    'hpbwv' :6.2,# half power beamwidth v
                    'hpbwh' :65, # half power beamwidth h
                    'sllv': -18, # side lobe level
                    'fbrh': 30,  # front back ratio
                    'gm': 18,    #
                    'pol':'p'    #t , p , c
                    }}


        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param'] = defaults['param']


        #if 'param' not in kwargs:
            #kwargs['param']=defaults['param']

        self.typ = "3gpp"
        self.param = kwargs['param']
        thtilt = self.param['thtilt']
        hpbwh  = self.param['hpbwh']
        hpbwv  = self.param['hpbwv']
        sllv   = self.param['sllv']
        fbrh   = self.param['fbrh']
        gm     = self.param['gm']
        pol    = self.param['pol']

        self.pol = pol
        # convert radian to degree

        phi   = self.phi*180/np.pi-180
        theta = self.theta*180/np.pi-90

        if self.grid:
            #Nth x Nph x Nf
            GvdB = np.maximum(-12*((theta-thtilt)/hpbwv)**2,sllv)[:,None,None]
            GhdB = (-np.minimum(12*(phi/hpbwh)**2,fbrh)+gm)[None,:,None]
            GdB  = GhdB+GvdB
            self.sqG = np.sqrt(10**(GdB/10.))*np.ones(self.nf)[None,None,:]
            self.evaluated = True
        else:
            #Nd x Nf
            GvdB = np.maximum(-12*((theta-thtilt)/hpbwv)**2,sllv)[:,None]
            GhdB = (-np.minimum(12*(phi/hpbwh)**2,fbrh)+gm)[:,None]
            GdB  = GhdB+GvdB
            self.sqG = np.sqrt(10**(GdB/10.))
        # radiating functions are deduced from square root of gain
        Ft,Fp = self.radF()
        return Ft,Fp

    def __pvsh1(self,**kwargs):
        """ calculate pattern from VSH Coeffs (shape 1)

        Parameters
        ----------

        theta  : ndarray (1xNdir)
        phi    : ndarray (1xNdir)
        k      : int
            frequency index

        Returns
        -------

        Ft , Fp 

        """
        assert hasattr(self,'C'),'no spherical coefficient'
        assert hasattr(self.C.Br,'s1'),'no shape 1 coeff in vsh'
        
        if self.grid:
            theta = np.kron(self.theta, np.ones(self.nph))
            phi = np.kron(np.ones(self.nth),self.phi)
        else:
            theta = self.theta
            phi = self.phi

        Nt = len(theta)
        Np = len(phi)

        if self.grid:
            theta = np.kron(theta, np.ones(Np))
            phi = np.kron(np.ones(Nt),phi)

        nray = len(theta)

        Br = self.C.Br.s1[:, :, :]
        Bi = self.C.Bi.s1[:, :, :]
        Cr = self.C.Cr.s1[:, :, :]
        Ci = self.C.Ci.s1[:, :, :]

        L = self.C.Br.L1
        M = self.C.Br.M1
        # The - sign is necessary to get the good reconstruction
        #     deduced from observation
        #     May be it comes from a different definition of theta in SPHEREPACK
        ind = index_vsh(L, M)
        l = ind[:, 0]
        m = ind[:, 1]
        #
        V, W = VW(l, m, theta, phi)
        #
        # broadcasting along frequency axis
        #
        V = np.expand_dims(V,0)
        W = np.expand_dims(V,0)
        #
        #   k : frequency axis
        #   l : axis l (theta)
        #   m : axis m (phi)
        #
        Fth = np.eisum('klm,kilm->ki',Br,np.real(V.T)) - \
              np.eisum('klm,kilm->ki',Bi,np.imag(V.T)) + \
              np.eisum('klm,kilm->ki',Ci,np.real(W.T)) + \
              np.eisum('klm,kilm->ki',Cr,np.imag(W.T))

        Fph = -np.eisum('klm,kilm->ki',Cr,np.real(V.T)) + \
              np.eisum('klm,kilm->ki',Ci,np.imag(V.T)) + \
              np.eisum('klm,kilm->ki',Bi,np.real(W.T)) + \
              np.eisum('klm,kilm->ki',Br,np.imag(W.T))

        # here Nf x Nd

        Ft = Fth.transpose()
        Fp = Fph.transpose()

        # then Nd x Nf

        if self.grid:
        # Nth x Nph x Nf
            Ft = Ft.reshape(self.nth, self.nph,self.nf)
            Fp = Fp.reshape(self.nth, self.nph,self.nf)

        # last axis should be frequency 
        assert(Ft.shape[-1]==self.nf)
        assert(Fp.shape[-1]==self.nf)
        
        return Ft, Fp

    def __pvsh3(self,**kwargs):
        """ calculate pattern from vsh3


        """
        assert hasattr(self,'C'),'no spherical coefficient'
        assert hasattr(self.C.Br,'s3'),'no shape 3 coeff in vsh'

        if self.grid:
            theta = np.kron(self.theta, np.ones(self.nph))
            phi = np.kron(np.ones(self.nth),self.phi)
        else:
            theta = self.theta
            phi = self.phi

        Br  = self.C.Br.s3
        lBr = self.C.Br.ind3[:, 0]
        mBr = self.C.Br.ind3[:, 1]

        Bi  = self.C.Bi.s3
        Cr  = self.C.Cr.s3
        Ci  = self.C.Ci.s3

        L = lBr.max()
        M = mBr.max()

        # vector spherical harmonics basis functions

        # V, W = VW(lBr, mBr, theta, phi)
        V, W = VW(lBr, mBr, theta, phi)
        Fth = np.dot(Br, np.real(V.T)) - \
              np.dot(Bi, np.imag(V.T)) + \
              np.dot(Ci, np.real(W.T)) + \
              np.dot(Cr, np.imag(W.T))

        Fph = -np.dot(Cr, np.real(V.T)) + \
               np.dot(Ci, np.imag(V.T)) + \
               np.dot(Bi, np.real(W.T)) + \
               np.dot(Br, np.imag(W.T))

        # here Nf x Nd

        Ft = Fth.transpose()
        Fp = Fph.transpose()

        # then Nd x Nf

        if self.grid:
        # Nth x Nph x Nf
            Ft = Ft.reshape(self.nth, self.nph,self.nf)
            Fp = Fp.reshape(self.nth, self.nph,self.nf)

        # last axis should be frequency 
        assert(Ft.shape[-1]==self.nf)
        assert(Fp.shape[-1]==self.nf)
        
        return Ft,Fp

    def __psh3(self,**kwargs):
        """ calculate pattern for sh3
        
        Parameters
        ----------

        """
        assert hasattr(self,'S'),'no spherical coefficient'
        assert hasattr(self.S.Cx,'s3'),'no shape 3 coeff in ssh'

        if self.grid:
            theta = np.kron(self.theta, np.ones(self.nph))
            phi = np.kron(np.ones(self.nth),self.phi)
        else:
            theta = self.theta
            phi = self.phi

        cx = self.S.Cx.s3
        cy = self.S.Cy.s3
        cz = self.S.Cz.s3

        lmax = self.S.Cx.lmax
        Y ,indx = SSHFunc2(lmax, theta,phi)

        k = self.S.Cx.k2

        if self.grid:
            Ex = np.dot(cx,Y[k])
            Ey = np.dot(cy,Y[k])
            Ez = np.dot(cz,Y[k])
            Fth,Fph = CartToSphere(theta, phi, Ex, Ey,Ez, bfreq = True, pattern = True )

            Ft = Fth.transpose()
            Fp = Fph.transpose()
            Ft = Ft.reshape(self.nth, self.nph,self.nf)
            Fp = Fp.reshape(self.nth, self.nph,self.nf)
        else:
            Ex = np.dot(cx,Y[k])
            Ey = np.dot(cy,Y[k])
            Ez = np.dot(cz,Y[k])
            Fth,Fph = CartToSphere(theta, phi, Ex, Ey,Ez, bfreq = True, pattern = False)
            Ft = Fth.transpose()
            Fp = Fph.transpose()

        assert(Ft.shape[-1]==self.nf)
        assert(Fp.shape[-1]==self.nf)
        
        return Ft,Fp

    def __pwireplate(self,**kwargs):
        """ pattern wire plate antenna

        """
        defaults = {'param':{'t0' : 5*np.pi/6,
                             'GmaxdB': 5
                   }}

        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.typ='wireplate'
        self.param = kwargs['param']
        t0 = self.param['t0']
        GmaxdB = self.param['GmaxdB']
        Gmax = pow(GmaxdB/10.,10)
        sqGmax = np.sqrt(Gmax)

        uth1 = np.where(self.theta < t0)[0]
        uth2 = np.where(self.theta >= t0)[0]
        p = t0
        q = np.pi/2.
        A = np.array(([[3*p**2,2*p,1],[p**3,p**2,p],[q**3,q**2,q]]))
        Y = np.array(([0,1,1/(1.*sqGmax)]))
        self.poly = la.solve(A,Y)

        argth1 = np.abs(self.poly[0]*self.theta[uth1]**3
                      + self.poly[1]*self.theta[uth1]**2
                      + self.poly[2]*self.theta[uth1])

        argth2 = -(1/(np.pi-t0)**2)*(self.theta[uth2]-t0)**2+1
        argth = np.hstack((argth1,argth2))[::-1]

        if self.grid:
            Ft = sqGmax * (argth[:,None])
            Fp = sqGmax * (argth[:,None])
        else:
            Fat = sqGmax * argth
            Fap = sqGmax * argth
            Ft = np.dot(Fat[:,None],np.ones(len(self.fGHz))[None,:])
            Fp = np.dot(Fap[:,None],np.ones(len(self.fGHz))[None,:])

        return Ft,Fp


    def __pcst(self,**kwargs):
        """ read antenna in text format
        """
       
        defaults = {'param':{'p' : 2,
                    'directory':'ant/FF_Results_txt_port_1_2/',
                    'fGHz':np.arange(2,6.5,0.5)}}

        if 'param' not in kwargs or kwargs['param']=={}:
            param=defaults['param']
        else:
            param=kwargs['param']
       
        self.fGHz = param['fGHz']
        self.nf = len(self.fGHz)
        
        for f in param['fGHz']:
            if ((int(f*10))%10)==0:
               _filename1 = 'E_port'+str(param['p'])+'_f'+str(int(f))+'GHz.txt'
               _filename2 = 'E_port'+str(param['p'])+'_f'+str(int(f))+'Ghz.txt'
        #    print 'toto'
            else:
                _filename1 = 'E_port'+str(param['p'])+'_f'+str(f)+'GHz.txt'
                _filename2 = 'E_port'+str(param['p'])+'_f'+str(f)+'Ghz.txt'
        
            
            filename1 = pyu.getlong(_filename1, param['directory'])
            filename2 = pyu.getlong(_filename2, param['directory'])
            
            try:
                df = pd.read_csv(filename1,sep=';')
            except:
                df = pd.read_csv(filename2,sep=';')

            columns = df.columns
            theta = (df[columns[0]]*np.pi/180).values.reshape(72,37)
            phi = (df[columns[1]]*np.pi/180).values.reshape(72,37)
            modGrlzdB = df[columns[2]]
            mFt = df[columns[3]]
            pFt = df[columns[4]]
            mFp = df[columns[5]]
            pFp = df[columns[6]]
            ratiodB = df[columns[7]]
            Ft = (10**(mFt/20)*np.exp(1j*pFt*np.pi/180)).values.reshape(72,37)
            Fp = (10**(mFp/20)*np.exp(1j*pFp*np.pi/180)).values.reshape(72,37)
            Ft = Ft.swapaxes(0,1)
            Fp = Fp.swapaxes(0,1)
            try:
                tFt=np.concatenate((tFt,Ft[...,None]),axis=2)
                tFp=np.concatenate((tFp,Fp[...,None]),axis=2)
            except:
                tFt=Ft[...,None]
                tFp=Fp[...,None]
        self.phi = phi[:,0]
        self.theta = theta[0,:]
        self.nth = len(self.theta)
        self.nph = len(self.phi)
        Ft = tFt
        Fp = tFp 
        return Ft,Fp

    def __pHertz(self,**kwargs):
        """ Hertz dipole
        """
        defaults = {'param':{'le':np.array([0,0,1])}}


        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        #k = 2*np.pi*self.fGHz[None,None,None,:]/0.3
        param=kwargs['param']

        if self.grid:
            le = param['le'][:,None,None]
            xr = np.sin(self.theta)[None,:,None]*np.cos(self.phi)[None,None,:]
            yr = np.sin(self.theta)[None,:,None]*np.sin(self.phi)[None,None,:]
            zr = np.cos(self.theta)[None,:,None]*np.ones(len(self.phi))[None,None,:]
            r = np.concatenate((xr,yr,zr),axis=0)

            xp = -np.sin(self.phi)[None,None,:]*np.ones(len(self.theta))[None,:,None]
            yp =  np.cos(self.phi)[None,None,:]*np.ones(len(self.theta))[None,:,None]
            zp = np.zeros(len(self.phi))[None,None,:]*np.ones(len(self.theta))[None,:,None]
            ph = np.concatenate((xp,yp,zp),axis=0)

            xt = np.cos(self.theta)[None,:,None]*np.cos(self.phi)[None,None,:]
            yt = np.cos(self.theta)[None,:,None]*np.sin(self.phi)[None,None,:]
            zt = -np.sin(self.theta)[None,:,None]*np.ones(len(self.phi))[None,None,:]
            th = np.concatenate((xt,yt,zt),axis=0)

            vec = le - np.einsum('kij,kij->ij',le,r)[None,...]*r
            #G = 1j*30*k*vec
            Ft = np.sqrt(3/2.)*np.einsum('kij,kij->ij',vec,th)[...,None]
            Fp = np.sqrt(3/2.)*np.einsum('kij,kij->ij',vec,ph)[...,None]
        else:
            le = param['le'][:,None]
            xr = np.sin(self.theta)*np.cos(self.phi)
            yr = np.sin(self.theta)*np.sin(self.phi)
            zr = np.cos(self.theta)
            r = np.concatenate((xr,yr,zr),axis=0)

            xp = -np.sin(self.phi)
            yp =  np.cos(self.phi)
            zp = np.zeros(len(self.phi))
            ph = np.concatenate((xp,yp,zp),axis=0)

            xt = np.cos(self.theta)*np.cos(self.phi)
            yt = np.cos(self.theta)*np.sin(self.phi)
            zt = -np.sin(self.theta)
            th = np.concatenate((xt,yt,zt),axis=0)

            vec = le - np.einsum('ki,ki->i',le,r)[None,...]*r
            #G = 1j*30*k*vec
            Ft = np.sqrt(3/2.)*np.einsum('ki,ki->i',vec,th)[...,None]
            Fp = np.sqrt(3/2.)*np.einsum('ki,ki->i',vec,ph)[...,None]

        return Ft,Fp

    def __pHuygens(self,**kwargs):
        """ Huygens source

        param : dict

        le : direction of electric current
        n  : normal to aperture
        """
        defaults = {'param':{'le':np.array([0,0,1]),
                             'n':np.array([1,0,0])}}


        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        #k = 2*np.pi*self.fGHz[None,None,None,:]/0.3
        param=kwargs['param']

        if self.grid:
            le = param['le'][:,None,None]
            n  = param['n'][:,None,None]
            xr = np.sin(self.theta)[None,:,None]*np.cos(self.phi)[None,None,:]
            yr = np.sin(self.theta)[None,:,None]*np.sin(self.phi)[None,None,:]
            zr = np.cos(self.theta)[None,:,None]*np.ones(len(self.phi))[None,None,:]
            r = np.concatenate((xr,yr,zr),axis=0)

            xp = -np.sin(self.phi)[None,None,:]*np.ones(len(self.theta))[None,:,None]
            yp =  np.cos(self.phi)[None,None,:]*np.ones(len(self.theta))[None,:,None]
            zp = np.zeros(len(self.phi))[None,None,:]*np.ones(len(self.theta))[None,:,None]
            ph = np.concatenate((xp,yp,zp),axis=0)

            xt = np.cos(self.theta)[None,:,None]*np.cos(self.phi)[None,None,:]
            yt = np.cos(self.theta)[None,:,None]*np.sin(self.phi)[None,None,:]
            zt = -np.sin(self.theta)[None,:,None]*np.ones(len(self.phi))[None,None,:]
            th = np.concatenate((xt,yt,zt),axis=0)

            vec1 = le - np.einsum('kij,kij->ij',le,r)[None,...]*r
            cro1 = np.cross(le,n,axisa=0,axisb=0,axisc=0)
            vec2 = np.cross(cro1,r,axisa=0,axisb=0,axisc=0)
            vec  = vec1-vec2

            #G = 1j*30*k*vec
            Ft = np.sqrt(3/4.)*np.einsum('kij,kij->ij',vec,th)[...,None]
            Fp = np.sqrt(3/4.)*np.einsum('kij,kij->ij',vec,ph)[...,None]
            #Ft = np.einsum('kij,kij->ij',vec,th)[...,None]
            #Fp = np.einsum('kij,kij->ij',vec,ph)[...,None]
        else:
            le = param['le'][:,None]
            xr = np.sin(self.theta)*np.cos(self.phi)
            yr = np.sin(self.theta)*np.sin(self.phi)
            zr = np.cos(self.theta)
            r = np.concatenate((xr,yr,zr),axis=0)

            xp = -np.sin(self.phi)
            yp =  np.cos(self.phi)
            zp = np.zeros(len(self.phi))
            ph = np.concatenate((xp,yp,zp),axis=0)

            xt = np.cos(self.theta)*np.cos(self.phi)
            yt = np.cos(self.theta)*np.sin(self.phi)
            zt = -np.sin(self.theta)
            th = np.concatenate((xt,yt,zt),axis=0)

            vec1 = le - np.einsum('ki,ki->i',le,r)[None,...]*r
            cro1 = np.cross(le,n,axisa=0,axisb=0,axisc=0)
            vec2 = np.cross(cro1,r,axisa=0,axisb=0,axisc=0)
            vec  = vec1-vec2
            #G = 1j*30*k*vec
            Ft = np.sqrt(3)*np.einsum('ki,ki->i',vec,th)[...,None]
            Fp = np.sqrt(3)*np.einsum('ki,ki->i',vec,ph)[...,None]

        return Ft,Fp 

    def __pArray(self,**kwargs):
        """ Array factor

        Parameters
        ----------

        Sc : np.array
            coupling S matrix

        Notes
        -----

        Nd : Number of directions
        Np : Number of points
        Nf : Number of frequency

        """

        defaults = {'param':{'Sc':[]}}

        if 'param' not in kwargs or kwargs['param']=={}:
            kwargs['param']=defaults['param']

        self.param = kwargs['param']

        lamda = (0.3/self.fGHz)
        k     = 2*np.pi/lamda

        if self.grid:
            sx = np.sin(self.theta[:,None])*np.cos(self.phi[None,:])    # Ntheta x Nphi
            sy = np.sin(self.theta[:,None])*np.sin(self.phi[None,:])    # Ntheta x Nphi
            sz = np.cos(self.theta[:,None])*np.ones(len(self.phi))[None,:]   # Ntheta x Nphi
            sx = sx.reshape(self.nth*self.nph)
            sy = sy.reshape(self.nth*self.nph)
            sz = sz.reshape(self.nth*self.nph)
        else:
            sx = np.sin(self.theta)*np.cos(self.phi)    # ,Nd
            sy = np.sin(self.theta)*np.sin(self.phi)    # ,Nd
            sz = np.cos(self.theta)                     # ,Nd

        self.s  = np.vstack((sx,sy,sz)).T         # Nd x 3
        #
        # F = exp(+jk s.p)
        #

        lshp = np.array(self.p.shape)
        if len(lshp)>2:
            Np = np.prod(lshp[1:])
            p = self.p.reshape(3,Np)
        else:
            p = self.p
        
        Np = p.shape[1]
        self.Sc = self.param['Sc']
        if self.Sc==[]:
            # Sc : Np x Np x Nf
            self.Sc = np.eye(Np)[...,None]
            #Sc2 = np.random.rand(Np,Np)[...,None]
            #pdb.set_trace()

        lshw = np.array(self.w.shape)
        if len(lshw)>2:
            Np2 = np.prod(lshw[0:-1])
            assert(Np2==Np)
            w = self.w.reshape(Np,lshw[-1])
        else:
            w = self.w
        # s : Nd x 3
        # p : 3 x Np
        #
        # sdotp : Nd x Np

        sdotp  = np.dot(self.s,p)   # s . p
        
        for a in self.la:
            if not self.grid:
                a.eval(grid=self.grid,ph=self.phi,th=self.theta)
            else:
                a.eval(grid=self.grid)
            # aFt : Nt x Np x Nf  |Nd x Nf
            # aFp : Nt x Np x Nf  |Nd x Nf
            aFt = a.Ft
            aFp = a.Fp
        #
        # Force conversion to Nd x Nf
        #

        shF = aFt.shape
        aFt = aFt.reshape(np.prod(shF[0:-1]),shF[-1])
        aFp = aFp.reshape(np.prod(shF[0:-1]),shF[-1])
        
        #
        # Same pattern on each point
        #
        aFt = aFt[:,None,:]
        aFp = aFp[:,None,:]

        #
        # Nf : frequency
        # Nd : direction
        # Np : points or array element position
        # Nu : users # Not iumplemented
        #
        # w  : Np x Nf
        # Sc : Np x Np x Nf
        #
        #
        # w' = w.Sc   Np x  Nf
        #
        # Coupling is implemented here

        # Rules : The repeated index k is the common dimension of the product
        # w    :  Np(k) x Nf(i)
        # Sc   :  Np(k) x Np(m) x Nf(i)
        # wp   :  Np(m) x Nf(i)
        wp = np.einsum('ki,kmi->mi',w,self.Sc)

        # add direction axis (=0) in w

        #if len(.w.shape)==3:
        #    self.wp   = self.wp[None,:,:,:]

        # aFT :  Nd x Np x Nf
        # E   :  Nd x Np x Nf

        E    = np.exp(1j*k[None,None,:]*sdotp[:,:,None])

        #
        # wp  : Np x Nf 
        # Fp  : Nd x Np x Nf
        # Ft  : Nd x Np x Nf
        #
        Ft = wp[None,...]*aFt*E
        Fp = wp[None,...]*aFp*E

        if self.grid:
        #
        # Integrate over the Np points (axis =1)
        # only if self.grid
        # Fp  : Nd x Nf
        # Ft  : Nd x Nf
        #
            Ft = np.sum(Ft,axis=1)
            Fp = np.sum(Fp,axis=1)
            sh = Ft.shape
            Ft = Ft.reshape(self.nth,self.nph,sh[1])
            Fp = Fp.reshape(self.nth,self.nph,sh[1])

        return Ft,Fp

    def radF(self):
        """ evaluate radiation fonction w.r.t polarization

        self.pol : 't' : theta , 'p' : phi n, 'c' : circular 

        """
        assert self.pol in ['t','p','c']
        if self.pol=='p':
            Fp = self.sqG
            if len(self.sqG.shape)==3:
                Ft = np.array([0])*np.ones(len(self.fGHz))[None,None,:]
            else:
                Ft = np.array([0])*np.ones(len(self.fGHz))[None,:]

        if self.pol=='t':
            if len(self.sqG.shape)==3:
                Fp = np.array([0])*np.ones(len(self.fGHz))[None,None,:]
            else:
                Fp = np.array([0])*np.ones(len(self.fGHz))[None,:]
            Ft = self.sqG
        if self.pol=='c':
            Fp = (1./np.sqrt(2))*self.sqG
            Ft = (1j/np.sqrt(2))*self.sqG

        return Ft,Fp

    def gain(self):
        """  calculates antenna gain
        
        Returns
        -------

        self.G  : np.array(Nt,Np,Nf) dtype:float
            linear gain 
                  or np.array(Nr,Nf)
        self.sqG : np.array(Nt,Np,Nf) dtype:float 
            linear sqare root of gain 
                  or np.array(Nr,Nf)
        self.efficiency : np.array (,Nf) dtype:float 
            efficiency 
        self.hpster : np.array (,Nf) dtype:float
            half power solid angle :  1 ~ 4pi steradian 
        self.ehpbw : np.array (,Nf) dtyp:float 
            equivalent half power beamwidth (radians)

        Notes
        -----

        .. math:: G(\theta,phi) = |F_{\\theta}|^2 + |F_{\\phi}|^2
(
        """
        self.G = np.real( self.Fp * np.conj(self.Fp)
                       +  self.Ft * np.conj(self.Ft) )


        if self.grid:
            dt = self.theta[1]-self.theta[0]
            dp = self.phi[1]-self.phi[0]
            Nt = len(self.theta)
            Np = len(self.phi)
            Gs = self.G*np.sin(self.theta)[:,None,None]*np.ones(Np)[None,:,None]
            self.efficiency = np.sum(np.sum(Gs,axis=0),axis=0)*dt*dp/(4*np.pi)

            self.sqG = np.sqrt(self.G)
            self.GdB = 10*np.log10(self.G)
            # GdBmax (,Nf)
            # Get direction of Gmax and get the polarisation state in that direction 
            # 
            self.GdBmax = np.max(np.max(self.GdB,axis=0),axis=0)
            self.umax = np.array(np.where(self.GdB==self.GdBmax))[:,0]
            self.theta_max = self.theta[self.umax[0]]
            self.phi_max = self.phi[self.umax[1]]
            M = geu.SphericalBasis(np.array([[self.theta_max,self.phi_max]]))
            self.sl = M[:,2].squeeze()
            uth = M[:,0] 
            uph = M[:,1] 
            el = self.Ft[tuple(self.umax)]*uth + self.Fp[tuple(self.umax)]*uph
            eln = el/np.linalg.norm(el)
            self.el = np.abs(eln.squeeze())
            self.hl = np.cross(self.sl,self.el)
            #assert((self.efficiency<1.0).all()),pdb.set_trace()
            self.hpster=np.zeros(len(self.fGHz))
            self.ehpbw=np.zeros(len(self.fGHz))
            for k in range(len(self.fGHz)):
                U  = np.zeros((Nt,Np))
                A = self.GdB[:,:,k]*np.ones(Nt)[:,None]*np.ones(Np)[None,:]
                u = np.where(A>(self.GdBmax[k]-3))
                U[u] = 1
                V  = U*np.sin(self.theta)[:,None]
                self.hpster[k] = np.sum(V)*dt*dp/(4*np.pi)
                self.ehpbw[k] = np.arccos(1-2*self.hpster[k])
        else:
            self.sqG = np.sqrt(self.G)
            self.GdB = 10*np.log10(self.G)





    def plotG(self,**kwargs):
        """ antenna plot gain in 2D

        Parameters
        ----------

        fGHz : frequency
        plan : 'theta' | 'phi' depending on the selected plan to be displayed
        angdeg : phi or theta in degrees, if plan=='phi' it corresponds to theta
        GmaxdB :  max gain to be displayed
        polar : boolean

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
            >>> A = Antenna('defant.vsh3')
            >>> fig,ax = A.plotG(fGHz=[2,3,4],plan='theta',angdeg=0)
            >>> fig,ax = A.plotG(fGHz=[2,3,4],plan='phi',angdeg=90)

        """

        if not self.evaluated:
            self.eval(pattern=True)

        dtr = np.pi/180.

        defaults = {'fGHz' : [],
                    'dyn' : 8 ,
                    'plan': 'phi',
                    'angdeg' : 90,
                    'legend':True,
                    'GmaxdB':20,
                    'polar':True,
                    'topos':False,
                    'source':'satimo',
                    'show':True,
                    'mode':'index',
                    'color':'black',
                    'u':0,
                    }

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
            #ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, facecolor='#d5de9c')
            if kwargs['polar']:
                ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True )
            else:
                ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']

        u = kwargs['u']
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

        #if len(self.fGHz) > 1 :
        #    fstep = self.fGHz[1]-self.fGHz[0]
        #else :
        #    fstep = np.array((abs(self.fGHz-kwargs['fGHz'][0])+1))
        #dtheta = self.theta[1,0]-self.theta[0,0]
        #dphi = self.phi[0,1]-self.phi[0,0]
        dtheta = self.theta[1]-self.theta[0]
        dphi = self.phi[1]-self.phi[0]

        if kwargs['fGHz']==[]:
            lfreq = [self.fGHz[0]]
        else:
            lfreq = kwargs['fGHz']

        for f in lfreq:
            df  = abs(self.fGHz-f)
            ik0 = np.where(df==min(df))
            ik = ik0[0][0]
            #ik=0
            chaine = 'f = %3.2f GHz' %(self.fGHz[ik])
            # all theta
            if kwargs['plan']=='theta':
                itheta = np.arange(self.nth)
                iphi1 = np.where(abs(self.phi-kwargs['angdeg']*dtr)<dphi)[0][0]
                Np = self.nph
                #   0 <= theta  <= pi/2
                u1 = np.where((self.theta <= np.pi / 2.) & (self.theta >= 0))[0]
                #   0 < theta < pi
                u2 = np.arange(self.nth)
                #   pi/2 < theta  <= pi
                u3 = np.nonzero((self.theta <= np.pi) & ( self.theta > np.pi / 2))[0]

                #
                # handle broadcasted axis =1 --> index 0
                shsqG = self.sqG.shape
                if shsqG[0]==1:
                    u1 = 0
                    u2 = 0
                    u3 = 0
                if shsqG[1]==1:
                    iphi1 = 0
                    iphi2 = 0
                if len(shsqG)==3:  # if only one frequency point
                    if shsqG[2]==1:
                        ik = 0
                else:
                    if shsqG[3]==1:
                        ik = 0

                # handle parity
                if np.mod(Np, 2) == 0:
                    iphi2 = np.mod(iphi1 + Np / 2, Np)
                else:
                    iphi2 = np.mod(iphi1 + (Np - 1) / 2, Np)

                if len(shsqG)==3:
                    arg1 = (u1,iphi1,ik)
                    arg2 = (u2,iphi2,ik)
                    arg3 = (u3,iphi1,ik)
                else:
                    if shsqG[3]==1:
                        u = 0
                    arg1 = (u1,iphi1,u,ik)
                    arg2 = (u2,iphi2,u,ik)
                    arg3 = (u3,iphi1,u,ik)

                # polar diagram
                #pdb.set_trace()
                if kwargs['polar']:
                    if kwargs['source']=='satimo':
                        r1 = -GmindB + 20 * np.log10( self.sqG[arg1]+1e-12)
                        r2 = -GmindB + 20 * np.log10( self.sqG[arg2]+1e-12)
                        r3 = -GmindB + 20 * np.log10( self.sqG[arg3]+1e-12)
                        #print max(r1)+GmindB
                        #print max(r2)+GmindB
                        #print max(r3)+GmindB
                    if kwargs['source']=='cst':
                        r1 = -GmindB + 20 * np.log10( self.sqG[arg1]/np.sqrt(30)+1e-12)
                        r2 = -GmindB + 20 * np.log10( self.sqG[arg2]/np.sqrt(30)+1e-12)
                        r3 = -GmindB + 20 * np.log10( self.sqG[arg3]/np.sqrt(30)+1e-12)

                    if type(r1)!= np.ndarray:
                        r1 = np.array([r1])*np.ones(len(self.phi))
                    if type(r2)!= np.ndarray:
                        r2 = np.array([r2])*np.ones(len(self.phi))
                    if type(r3)!= np.ndarray:
                        r3 = np.array([r3])*np.ones(len(self.phi))

                    negr1 = np.nonzero(r1 < 0)
                    negr2 = np.nonzero(r2 < 0)
                    negr3 = np.nonzero(r3 < 0)

                    r1[negr1[0]] = 0
                    r2[negr2[0]] = 0
                    r3[negr3[0]] = 0

                    r = np.hstack((r1[::-1], r2, r3[::-1], r1[-1]))

                    a1 = np.arange(0, 360, 30)
                    a2 = [90, 60, 30, 0, 330, 300, 270, 240, 210, 180, 150, 120]
                    rline2, rtext2 = plt.thetagrids(a1, a2)

                # linear diagram
                else:
                    r1 = 20 * np.log10( self.sqG[arg1]+1e-12)
                    r2 = 20 * np.log10( self.sqG[arg2]+1e-12)
                    r3 = 20 * np.log10( self.sqG[arg3]+1e-12)

                    r = np.hstack((r1[::-1], r2, r3[::-1], r1[-1]))

                # angular basis for phi
                angle = np.linspace(0, 2 * np.pi, len(r), endpoint=True)
                plt.title(u'$\\theta$ plane')

            if kwargs['plan']=='phi':
                iphi = np.arange(self.nph)
                itheta = np.where(abs(self.theta-kwargs['angdeg']*dtr)<dtheta)[0][0]
                angle = self.phi[iphi]
                if len(self.sqG.shape)==3:
                    arg = [itheta,iphi,ik]
                else:
                    arg = [itheta,iphi,u,ik]
                if kwargs['polar']:
                    if np.prod(self.sqG.shape)!=1:
                        r = -GmindB + 20 * np.log10(self.sqG[arg])
                        neg = np.nonzero(r < 0)
                        r[neg] = 0
                    else:
                        r = -GmindB+ 20*np.log10(self.sqG[0,0,0]*np.ones(np.shape(angle)))
               # plt.title(u'H plane - $\phi$ degrees')
                    a1 = np.arange(0, 360, 30)
                    a2 = [0, 30, 60, 90, 120 , 150 , 180 , 210, 240 , 300 , 330]
                    #rline2, rtext2 = plt.thetagrids(a1, a2)
                else:
                    r =  20 * np.log10(self.sqG[arg])

                plt.title(u'$\\phi$ plane ')
            # actual plotting
            if len(lfreq)>1: 
                ax.plot(angle, r, color=col[cpt], lw=2, label=chaine)
            else:
                ax.plot(angle, r, color=kwargs['color'], lw=2, label=chaine)
            cpt = cpt + 1

        if kwargs['polar']:
            rline1, rtext1 = plt.rgrids(t1, t2)
            #ax.set_rmax(t1[-1])
            #ax.set_rmin(t1[0])
        if kwargs['legend']:
            ax.legend()
        if kwargs['show']:
            plt.ion()
            plt.show()
        return(fig,ax)

class Antenna(Pattern):
    """ Antenna

    Attributes
    ----------

    name   : Antenna name

    nf     : number of frequency
    nth    : number of theta
    nph    : number of phi

    Ft     : Normalized Ftheta    (ntheta,nphi,nf)
    Fp     : Normalized Fphi      (ntheta,nphi,nf)
    sqG    : square root of gain  (ntheta,nphi,nf)


    theta  : theta base            1 x ntheta
    phi    : phi base              1 x phi

    C      : VSH Coefficients

    Methods
    -------

    info  : Display information about antenna
    vsh   : calculates Vector Spherical Harmonics
    show3 : Geomview diagram
    plot3d : 3D diagram plotting using matplotlib toolkit

    Antenna trx file can be stored in various order
        natural : HFSS
        ncp     : near filed chamber

    It is important when initializing an antenna object
    to be aware of the typ of trx file

    .trx (ASCII Vectorial antenna Pattern)

    F   Phi   Theta  Fphi  Ftheta

    """


    def __init__(self,typ='Omni',**kwargs):
        """ class constructor

        Parameters
        ----------

        typ  : 'Omni','Gauss','WirePlate','3GPP','atoll'

        _filename : string
                    antenna file name
        directory : str
                    antenna subdirectory of the current project
                    the file is seek in the $BASENAME/ant directory
        nf        : integer
                     number of frequency 
        ntheta    : integer
                    number of theta (default 181)
        nphi      : integer
                    number of phi (default 90)
        source    : string
                source of data { 'satimo' | 'cst' | 'hfss' }

        Notes
        -----

        The supported data formats for storing antenna patterns are

        'mat': Matlab File
        'vsh2': unthresholded vector spherical coefficients
        'vsh3': thresholded vector spherical cpoefficients
        'atoll': Atoll antenna file format
        'trx' : Satimo NFC raw data
        'trx1' : Satimo NFC raw data  (deprecated)

         A = Antenna('my_antenna.mat')

        """
        defaults = {'directory': 'ant',
                    'source':'satimo',
                    'ntheta':90,
                    'nphi':181,
                    'L':90, # L max
                    'param':{}
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if 'fGHz' in kwargs:
            if type(kwargs['fGHz'])==np.ndarray:
                self.fGHz=kwargs['fGHz']
            else:
                self.fGHz=np.array([kwargs['fGHz']])


        #mayavi selection
        self._is_selected=False


        self.source = kwargs['source']

        self.param = kwargs['param']

        #super(Antenna,self).__init__()
        #Pattern.__init__(self)
        #
        # if typ string has an extension it is a file
        #
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
        #determine if pattern for all theta/phi is constructed
        self.full_evaluated = False

        if self.fromfile:
            if isinstance(typ,str):
                self._filename = typ
                if self.ext == 'vsh3':
                    self.typ='vsh3'
                    self.loadvsh3()
                if self.ext == 'vsh2':
                    self.typ='vsh2'
                    self.loadvsh2()
                if self.ext == 'sh3':
                    self.typ='sh3'
                    self.loadsh3()
                if self.ext == 'sh2':
                    self.typ='sh2'
                    self.loadsh2()
                if self.ext == 'trx1':
                    self.typ='trx'
                    self.load_trx(kwargs['directory'],self.nf,self.nth,self.nph)
                if self.ext == 'trx':
                    self.typ='trx'
                    self.loadtrx(kwargs['directory'])
                if self.ext == 'mat':
                    self.typ='mat'
                    self.loadmat(kwargs['directory'])
                if self.ext == 'cst':
                    self.typ='cst'
                if self.ext == 'txt':
                    self.typ='atoll'
                    self.load_atoll(kwargs['directory'])
            elif isinstance(typ,list):
                self._filename = typ
                self.ext='hfss'
                self.loadhfss(typ, self.nth, self.nph)

        else:
            self.typ=typ
            self._filename=typ
            if self.typ=='vsh3':
                self.initvsh()
            else:
                self.eval()

    def __repr__(self):
        st = ''
        st = st + 'Antenna type : ' + self.typ +'\n'
        st = st+'------------------------\n'
        if 'param' in self.__dict__:
            for k in self.param:
                st = st + ' ' + k + ' : ' + str(self.param[k])+'\n'
        if hasattr(self,'atoll'):
            for k1 in self.atoll.keys():
                st = st + str(k1)+'\n'
                for k2 in self.atoll[k1]:
                    st = st + ' '+ str(k2)+'\n'
        st = st+'------------------------\n'
        rtd = 180./np.pi
        if self.fromfile:
            if isinstance(self._filename,str):
                st = st + 'file name : ' + self._filename+'\n'
            else:
                for i in range(len(self._filename)):
                    st = st + 'FileName : ' + self._filename[i]+'\n'
#        #st = st + 'file type : ' + self.typ+'\n'
        if 'fGHz' in self.__dict__:
            st = st + "fmin : %4.2f" % (self.fGHz[0]) + "GHz\n"
            st = st + "fmax : %4.2f" % (self.fGHz[-1]) + "GHz\n"
            try:
                st = st + "step : %4.2f" % (1000*(self.fGHz[1]-self.fGHz[0])) + "MHz\n"
            except:
                st = st + "step : None\n"
            st = st + "Nf : %d" % (len(self.fGHz)) +"\n"
#
#

        if hasattr(self,'C'):
            st = st + self.C.__repr__()

        if hasattr(self,'S'):
            st = st + self.S.__repr__()

        if self.evaluated:
            st = st + '-----------------------\n'
            st = st + '      evaluated        \n'
            st = st + '-----------------------\n'
            st = st + "Ntheta : %d" % (self.nth) + "\n"
            st = st + "Nphi : %d" % (self.nph) + "\n"
#                kwargs[k] = defaults[k]

            u = np.where(self.sqG==self.sqG.max())
            if self.grid:
                if len(u[0])>1:
                    S = self.sqG[(u[0][0],u[1][0],u[2][0])]
                    ut = u[0][0]
                    up = u[1][0]
                    uf = u[2][0]
                else:
                    S = self.sqG[u]
                    ut = u[0]
                    up = u[1]
                    uf = u[2]
            else:
                if len(u[0])>1:
                    S = self.sqG[(u[0][0],u[1][0])]
                    ud = u[0][0]
                    uf = u[1][0]
                else:
                    S = self.sqG[u]
                    ud = u[0]
                    uf = u[1]

            st = st + "GdBmax :"+str(self.GdBmax[0])+' '+str(self.GdBmax[-1])+'\n'
            st = st + "Gmax direction : .sl" + str(self.sl)+'\n'
            st = st + "Orientation of E field in Gmax direction : .el " + str(self.el)+'\n'
            st = st + "Orientation of H field in Gmax direction : .hl " + str(self.hl)+'\n'
            st = st + "effective HPBW : .ehpbw " + str(self.ehpbw[0])+' '+str(self.ehpbw[-1])+'\n'

            if self.source=='satimo':
                GdB = 20*np.log10(S)
            # see WHERE1 D4.1 sec 3.1.1.2.2
            if self.source=='cst':
                GdB = 20*np.log10(S/np.sqrt(30))
            #st = st + "GmaxdB : %4.2f dB \n" % (GdB)
            st = st + "   f = %4.2f GHz \n" % (self.fGHz[uf])
            if self.grid:
                st = st + "   theta = %4.2f (degrees) \n" % (self.theta[ut]*rtd)
                st = st + "   phi = %4.2f  (degrees) \n" % (self.phi[up]*rtd)
            else:
                st = st + " Ray n :" + str(ud)+' \n'
        else:
            st = st + 'Not evaluated\n'
#
#
#        if self.typ == 'mat':
#            #st = st + self.DataFile + '\n'
#            st = st + 'antenna name : '+ self.AntennaName + '\n'
#            st = st + 'date : ' + self.Date +'\n'
#            st = st + 'time : ' + self.StartTime +'\n'
#            st = st + 'Notes : ' + self.Notes+'\n'
#            st = st + 'Serie : ' + str(self.Serie)+'\n'
#            st = st + 'Run : ' + str(self.Run)+'\n'
#            st = st + "Nb theta (lat) : "+ str(self.nth)+'\n'
#            st = st + "Nb phi (lon) :"+ str(self.nph)+'\n'
#
#        if self.typ == 'Gauss':
#            st = st + 'Gaussian pattern' + '\n'
#            st = st + 'phi0 : ' + str(self.p0) +'\n'
#            st = st + 'theta0 :' + str(self.t0) + '\n'
#            st = st + 'phi 3dB :' + str(self.p3) + '\n'
#            st = st + 'theta 3dB :' + str(self.t3) + '\n'
#            st = st + 'Gain dB :' + str(self.GdB) + '\n'
#            st = st + 'Gain linear :' + str(self.G ) + '\n'
#            st = st + 'sqrt G :' + str(self.sqG) + '\n'

        return(st)

    def initvsh(self,lmax=45):
        """ Initialize a void vsh structure

        Parameters
        ----------

        fGHz : array
        lmax : int
            level max

        """
        nf = len(self.fGHz)
        Br = 1j * np.zeros((nf, lmax, lmax-1))
        Bi = 1j * np.zeros((nf, lmax, lmax-1))
        Cr = 1j * np.zeros((nf, lmax, lmax-1))
        Ci = 1j * np.zeros((nf, lmax, lmax-1))
        Br = VCoeff(typ='s1', fmin=self.fGHz[0], fmax=self.fGHz[-1], data=Br)
        Bi = VCoeff(typ='s1', fmin=self.fGHz[0], fmax=self.fGHz[-1], data=Bi)
        Cr = VCoeff(typ='s1', fmin=self.fGHz[0], fmax=self.fGHz[-1], data=Cr)
        Ci = VCoeff(typ='s1', fmin=self.fGHz[0], fmax=self.fGHz[-1], data=Ci)
        self.C = VSHCoeff(Br, Bi, Cr, Ci)

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
        """ show a picture of the antenna 

        Parameters
        ----------

        directory : string

        """

        if directory == '':
            directory = os.path.join('ant','UWBAN','PhotosVideos')

        _filename = 'IMG_'+self.PhotoFile.split('-')[1]+'.JPG'
        filename = pyu.getlong(_filename,directory)
        if sys.version_info.major==2:
            I = Image.open(filename)
        else:
            I = image.open(filename)

        I.show()


    def load_atoll(self,directory="ant"):
        """ load antenna from Atoll file 
        
        Atoll format provides Antenna gain given for the horizontal and vertical plane 
        for different frequencies and different tilt values 

        Parameters
        ----------

        directory : string 

        The dictionnary attol is created 

        """
        _filemat = self._filename
        fileatoll = pyu.getlong(_filemat, directory)
        fd = open(fileatoll)
        lis = fd.readlines()
        tab = []
        for li in lis:
            lispl= li.split('\t')
            if (lispl[0]!=''):
                tab.append(lispl)

        deg_to_rad = np.pi/180.
        lbs_to_kg = 0.45359237
        columns = tab[0]
        #pdb.set_trace()
        for k in np.arange(len(tab)-1):
            df = pd.DataFrame([tab[k+1]],columns=columns)
            try:
                dff=dff.append(df)
            except:
                dff= df
        self.raw = dff
        dff = dff.iloc[:,[0,8,9,10,2,5,7,14,11,16,17,13,6,12]]
        #dff = df['Name','Gain  (dBi)','FMin','FMax','FREQUENCY','Pattern','V_WIDTH','H_WIDTH','DIMENSIONS HxWxD   (INCHES)','WEIGHT (LBS)']
        dff.columns = ['Name','Fmin','Fmax','F','Gmax','G','Hpbw','H_width','V_width','HxWxD','Weight','Tilt','Etilt','Ftob']
        dff=dff.apply(lambda x :pd.to_numeric(x,errors='ignore'))
        #
        # Parse polarization in the field name
        #
        upolarp45 = ['(+45)' in x for x in dff['Name']]
        upolarm45 = ['(-45)' in x for x in dff['Name']]  
        if (sum(upolarp45)>0):
            dff.loc[upolarp45,'Polar']=45
        if (sum(upolarm45)>0):
            dff.loc[upolarm45,'Polar']=-45

        atoll = {}
        dfband = dff.groupby(['Fmin'])
        for b in dfband:
            keyband = str(b[0])+'-'+str(b[1]['Fmax'].values[0])
            atoll[keyband]={}  # band
            dfpol = b[1].groupby(['Polar'])
            for p in dfpol:
                atoll[keyband][p[0]] = {} # polar
                dftilt = p[1].groupby(['Tilt'])
                Ghor = np.empty((360,1))  # angle , tilt , frequency
                Gver = np.empty((360,1))  # angle , 
                ct = 0
                tilt = []
                for t in dftilt:
                    dffreq = t[1].groupby(['F'])
                    ct+=1
                    cf=0 
                    tilt.append(t[0])
                    freq = []
                    for f in dffreq:
                        freq.append(f[0])
                        cf+=1
                        if len(f[1])==1:
                            df = f[1]
                        else:
                            df = f[1].iloc[0:1]
                        Gmax = df['Gmax'].values
                        str1 = df.loc[:,'G'].values[0].replace('  ',' ')
                        lstr = str1.split(' ')
                        Pattern = [ eval(x) for x in lstr[0:-1]]
                        # 4 fist field / # of points
                        Nd,db,dc,Np = Pattern[0:4]
                        #print(Nd,b,c,Np)
                        tmp = np.array(Pattern[4:4+2*Np]).reshape(Np,2)
                        ah   = tmp[:,0]
                        ghor = Gmax-tmp[:,1]
                        # 4 fist field / # of points
                        da,db,dc,dd = Pattern[4+2*Np:4+2*Np+4]
                        #pdb.set_trace()
                        #print a,b,c,d
                        tmp = np.array(Pattern[4+2*Np+4:]).reshape(dc,2)
                        gver = Gmax-tmp[:,0]
                        av = tmp[:,1]
                        try:
                            Ghor = np.hstack((Ghor,ghor[:,None]))
                            Gver = np.hstack((Gver,gver[:,None]))
                        except:
                            pdb.set_trace()
                Ghor = np.delete(Ghor,0,1)
                Gver = np.delete(Gver,0,1)
                atoll[keyband][p[0]]['hor'] = Ghor.reshape(360,ct,cf)
                atoll[keyband][p[0]]['ver'] = Gver.reshape(360,ct,cf)
                atoll[keyband][p[0]]['tilt'] = np.array(tilt)
                atoll[keyband][p[0]]['freq'] = np.array(freq)
        self.atoll = atoll
        # Gmax = eval(self.df['Gain  (dBi)'].values[0])
        #fig = plt.figure()
        #ax =plt.gca(projection='polar')
        #ax =plt.gca()
        #ax.plot(H2[:,1]*deg_to_rad,Gain-H2[:,0],'r',label='vertical',linewidth=2)
        #ax.plot(H1[:,0]*deg_to_rad,Gain-H1[:,1],'b',label='horizontal',linewidth=2)
        #ax.set_rmin(-30)
        #plt.title(dir1+'/'+filename+' Gain : '+df['Gain  (dBi)'].values[0])
        #BXD-634X638XCF-EDIN.txt
        #BXD-636X638XCF-EDIN.txt        

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
            >>> f,a = A.plotG(plan='theta',angdeg=0)
            >>> f,a = A.plotG(plan='phi',angdeg=90,fig=f,ax=a)
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

        self.fGHz = d.freq/1.e9
        self.theta = d.theta
        self.phi = d.phi
        self.Ft = d.Ftheta
        self.Fp = d.Fphi

        self.Fp = self.Fp.swapaxes(0, 2)
        self.Fp = self.Fp.swapaxes(0, 1)
        self.Ft = self.Ft.swapaxes(0, 2)
        self.Ft = self.Ft.swapaxes(0, 1)
        Gr = np.real(self.Fp * np.conj(self.Fp) + \
                     self.Ft * np.conj(self.Ft))
        self.sqG = np.sqrt(Gr)
        self.nth = len(self.theta)
        self.nph = len(self.phi)

        if type(self.fGHz) ==  float:
            self.nf = 1
        else:
            self.nf = len(self.fGHz)

        self.evaluated = True
        self.grid = True

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

        TODO : DEPRECATED (Fix the Ft and Fp format with Nf as last axis)

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

        self.Fp = Fphi.reshape((nf, nphi, ntheta))
        self.Ft = Ftheta.reshape((nf, nphi, ntheta))
        Ttheta = THETA.reshape((nf, nphi, ntheta))
        Tphi = PHI.reshape((nf, nphi, ntheta))
        Tf = F.reshape((nf, nphi, ntheta))

        self.Fp = self.Fp.swapaxes(1, 2)
        self.Ft = self.Ft.swapaxes(1, 2)
        Ttheta = Ttheta.swapaxes(1, 2)
        Tphi = Tphi.swapaxes(1, 2)
        Tf = Tf.swapaxes(1, 2)

        self.fGHz = Tf[:, 0, 0]
        self.theta = Ttheta[0, :, 0]
        #self.phi     = Tphi[0,0,:]

        #
        # Temporaire
        #
        A1 = self.Fp[:, 90:181, :]
        A2 = self.Fp[:, 0:91, :]
        self.Fp = np.concatenate((A1, A2[:, ::-1, :]), axis=2)
        A1 = self.Ft[:, 90:181, :]
        A2 = self.Ft[:, 0:91, :]
        self.Ft = np.concatenate((A1, A2[:, ::-1, :]), axis=2)
        self.theta = np.linspace(0, np.pi, 91)
        self.phi = np.linspace(0, 2 * np.pi, 180, endpoint=False)
        self.nth = 91
        self.nph = 180
        self.nf = 104
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

        self.grid = True

        Nt = len(theta)
        Np = len(phi)
        Nf = len(self.fGHz)
        #Th = np.kron(theta, np.ones(Np))
        #Ph = np.kron(np.ones(Nt), phi)
        if typ =='s1':
            FTh, FPh = self.Fsynth1(theta, phi)
        if typ =='s2':
            FTh, FPh = self.Fsynth2b(theta,phi)
        if typ =='s3':
            FTh, FPh = self.Fsynth3(theta, phi)
        #FTh = Fth.reshape(Nf, Nt, Np)
        #FPh = Fph.reshape(Nf, Nt, Np)

        return(FTh,FPh)

    def coeffshow(self,**kwargs):
        """ display antenna coefficient

            typ : string
                'ssh' |'vsh'
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
            FTh, FPh = self.Fsynth1(theta, phi)
        if typ =='s2':
            FTh, FPh = self.Fsynth2b(theta, phi)
        if typ =='s3':
            FTh, FPh = self.Fsynth3(theta, phi)

        #FTh = Fth.reshape(self.nf, Nt, Np)
        #FPh = Fph.reshape(self.nf, Nt, Np)
        #
        #  Jacobian
        #
        #st    = outer(sin(theta),ones(len(phi)))
        st = np.sin(theta).reshape((len(theta), 1))
        #
        # Construct difference between reference and reconstructed
        #
        if kf!=-1:
            dTh = (FTh[kf, :, :] - self.Ft[kf, ::dsf, ::dsf])
            dPh = (FPh[kf, :, :] - self.Fp[kf, ::dsf, ::dsf])
            #
            # squaring  + Jacobian
            #
            dTh2 = np.real(dTh * np.conj(dTh)) * st
            dPh2 = np.real(dPh * np.conj(dPh)) * st

            vTh2 = np.real(self.Ft[kf, ::dsf, ::dsf] \
                 * np.conj(self.Ft[kf, ::dsf, ::dsf])) * st
            vPh2 = np.real(self.Fp[kf, ::dsf, ::dsf] \
                 * np.conj(self.Fp[kf, ::dsf, ::dsf])) * st

            mvTh2 = np.sum(vTh2)
            mvPh2 = np.sum(vPh2)

            errTh = np.sum(dTh2)
            errPh = np.sum(dPh2)
        else:
            dTh = (FTh[:, :, :] - self.Ft[:, ::dsf, ::dsf])
            dPh = (FPh[:, :, :] - self.Fp[:, ::dsf, ::dsf])
            #
            # squaring  + Jacobian
            #
            dTh2 = np.real(dTh * np.conj(dTh)) * st
            dPh2 = np.real(dPh * np.conj(dPh)) * st

            vTh2 = np.real(self.Ft[:, ::dsf, ::dsf] \
                 * np.conj(self.Ft[:, ::dsf, ::dsf])) * st
            vPh2 = np.real(self.Fp[:, ::dsf, ::dsf] \
                 * np.conj(self.Fp[:, ::dsf, ::dsf])) * st

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
        self.nf = len(lfa)
        fGHz  = []
        lacsv = []
        Fphi = np.empty((self.nf,self.nth,self.nph))
        Ftheta = np.empty((self.nf,self.nth,self.nph))
        SqG = np.empty((self.nf,self.nth,self.nph))

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


        self.fGHz = np.array(fGHz)
        #self.theta = th[0,:].reshape(Nt,1)
        #self.phi = ph[:,0].reshape(1,Np)
        self.theta = th[0,:]
        self.phi = ph[:,0]
        self.Fp=Fphi
        self.Ft=Ftheta
        self.sqG=SqG



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
            print("error : frequency cannot be zero")
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
        # nfc
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
            typ = 'nfc'
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
            self.Fp = Fphi.reshape((nf, ntheta, nphi))
            self.Ft = Ftheta.reshape((nf, ntheta, nphi))
            self.sqG = SqG.reshape((nf, ntheta, nphi))
            Ttheta = theta.reshape((nf, ntheta, nphi))
            Tphi = phi.reshape((nf, ntheta, nphi))
            Tf = f.reshape((nf, ntheta, nphi))
        if typ == 'nfc':
            self.Fp = Fphi.reshape((nf, nphi, ntheta))
            self.Ft = Ftheta.reshape((nf, nphi, ntheta))
            self.sqG = SqG.reshape((nf, nphi, ntheta))
            Ttheta = theta.reshape((nf, nphi, ntheta))
            Tphi = phi.reshape((nf, nphi, ntheta))
            Tf = f.reshape((nf, nphi, ntheta))
        #
        # Force natural order (f,theta,phi)
        # This is not the order of the satimo nfc which is  (f,phi,theta)
        #

            self.Fp = self.Fp.swapaxes(1, 2)
            self.Ft = self.Ft.swapaxes(1, 2)
            self.sqG = self.sqG.swapaxes(1, 2)
            Ttheta = Ttheta.swapaxes(1, 2)
            Tphi = Tphi.swapaxes(1, 2)
            Tf = Tf.swapaxes(1, 2)

        self.fGHz = Tf[:, 0, 0]
        self.theta = Ttheta[0, :, 0]
        self.phi = Tphi[0, 0, :]
        #
        # check header consistency
        #
        np.testing.assert_almost_equal(self.fGHz[0],fmin,6)
        np.testing.assert_almost_equal(self.fGHz[-1],fmax,6)
        np.testing.assert_almost_equal(self.theta[0],thmin,3)
        np.testing.assert_almost_equal(self.theta[-1],thmax,3)
        np.testing.assert_almost_equal(self.phi[0],phmin,3)
        np.testing.assert_almost_equal(self.phi[-1],phmax,3)

        self.nf = nf
        self.nth = ntheta
        self.nph = nphi
        self.tau = tau

        self.evaluated = True

    def checkpole(self, kf=0):
        """ display the reconstructed field on pole for integrity verification

        Parameters
        ----------

        kf : int
             frequency index default 0

        """
        Ft0 = self.Ft[kf, 0, :]
        Fp0 = self.Fp[kf, 0, :]
        Ftp = self.Ft[kf, -1, :]
        Fpp = self.Fp[kf, -1, :]
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
        print(self._filename)
        print("type : ", self.typ)
        if self.typ == 'mat':
            print(self.DataFile)
            print(self.AntennaName)
            print(self.Date)
            print(self.StartTime)
            print(self.Notes)
            print(self.Serie)
            print(self.Run)
            print("Nb theta (lat) :", self.nth)
            print("Nb phi (lon) :", self.nph)
        if self.typ =='nfc':
            print( "--------------------------")
            print( "fmin (GHz) :", self.fGHz[0])
            print( "fmax (GHz) :", self.fGHz[-1])
            print( "Nf   :", self.nf)
            print( "thmin (rad) :", self.theta[0])
            print( "thmax (rad) :", self.theta[-1])
            print( "Nth  :", self.nth)
            print( "phmin (rad) :", self.phi[0])
            print( "phmax (rad) :", self.phi[-1])
            print( "Nph  :", self.nph)
        try:
            self.C.info()
        except:
            print("No vsh coefficient calculated yet")

    #@mlab.show
    def _show3(self,bnewfig = True,
                    bcolorbar =True,
                    name=[],
                    binteract=False,
                    btitle=True,
                    bcircle=True,
                    **kwargs ):
        """ show3 mayavi

        Parameters
        ----------

        btitle : boolean
            display title
        bcolorbar : boolean
            display colorbar
        binteract : boolean 
            enable interactive mode
        newfig: boolean


        see also
        --------

        antprop.antenna._computemesh

        """

        if not self.evaluated:
            self.eval(pattern=True)

        # k is the frequency index
        if hasattr(self,'p'):
            lpshp = len(self.p.shape)
            sum_index = tuple(np.arange(1,lpshp))
            po = np.mean(self.p,axis=sum_index)
            kwargs['po']=po
        x, y, z, k, scalar  = self._computemesh(**kwargs)

        if bnewfig:
            mlab.clf()
            f=mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
        else :
            f=mlab.gcf()

        if kwargs.has_key('opacity'):
            opacity = kwargs['opacity']
        else: 
            opacity = 1

        self._mayamesh = mlab.mesh(x, y, z,
                                   scalars= scalar,
                                   resolution = 1,
                                   opacity = opacity,reset_zoom=False)

        if name == []:
            f.children[-1].name = 'Antenna ' + self._filename
        else :
            f.children[-1].name = name + self._filename

        if bcolorbar :
            mlab.colorbar()

        if btitle:
            mlab.title(self._filename + ' @ ' + str(self.fGHz[k]) + ' GHz',height=1,size=0.5)

        def circle(typ='xy',a=1.2):
            phi = np.linspace(0, 2*np.pi, 2000)
            if typ=='xy':
                return [ a*np.cos(phi) ,
                         a*np.sin(phi) ,
                         np.zeros(len(phi))
                         ]
            if typ=='yz':
                return [ np.zeros(len(phi)),
                         a*np.cos(phi) ,
                         a*np.sin(phi) 
                         ]
            if typ=='xz':
                return [ a*np.cos(phi),
                         a*np.zeros(len(phi)),
                         np.sin(phi) 
                         ]
        # draw 3D circle around pattern
        if bcircle:
            xc,yc,zc =circle('xy') # blue
            mlab.plot3d(xc,yc,zc,color=(0,0,1))
            xc,yc,zc =circle('yz') # red
            mlab.plot3d(xc,yc,zc,color=(1,0,0))
            xc,yc,zc =circle('xz') # green
            mlab.plot3d(xc,yc,zc,color=(0,1,0))

        if binteract:
            self._outline = mlab.outline(self._mayamesh, color=(.7, .7, .7))
            self._outline.visible=False
            def picker_callback(picker):
                """ Picker callback: this get called when on pick events.
                """
                if picker.actor in self._mayamesh.actor.actors:
                    self._outline.visible = not self._outline.visible
                    self._is_selected=self._outline.visible
            picker = f.on_mouse_pick(picker_callback)

        return(f)




    def _computemesh(self,**kwargs):
        """ compute mesh from theta phi

        Parameters
        ----------

        fGHz : np.array()
            default [] : takes center frequency fa[len(fa)/2]
        po   : np.array()
            location point of the antenna
        T    : np.array
            rotation matrix
        minr : float
            minimum radius in meter
        maxr : float
            maximum radius in meters
        tag : string
        ilog : boolean
        title : boolean


        Returns
        -------

        (x, y, z, k)

        x , y , z values in cartesian axis
        k frequency point evaluated

        """
        defaults = { 'fGHz' :[],
                     'po': np.array([0,0,0]),
                     'T' : np.eye(3),
                     'minr' : 0.1,
                     'maxr' : 1 ,
                     'scale':1.,
                     'tag' : 'Pat',
                     'txru' : 0,
                     'ilog' : False,
                     'title':True,

                     }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        fGHz = kwargs['fGHz']
        minr = kwargs['minr']
        maxr = kwargs['maxr']
        tag  = kwargs['tag']
        ilog = kwargs['ilog']
        txru = kwargs['txru']
        scale= kwargs['scale']

        po = kwargs['po']
        # T is an unitary matrix
        T  = kwargs['T']
        if fGHz == []:
            # self.ext == '' <=> mathematically generated => nf = 1
            if self.ext != '':
                k = len(self.fGHz)/2
            else: 
                k = 0
        else :
            if self.ext != '':
                k = np.where(self.fGHz>=fGHz)[0][0]
            else: 
                k = 0
        if len(self.Ft.shape)==3:
            r = self.sqG[:,:,k]
        else:
            r = self.sqG[:,:,txru,k]
        th = self.theta[:,None]
        phi = self.phi[None,:]

        if ilog :
            r = 10*np.log10(abs(r))
        else:
            r = abs(r)
        if r.max() != r.min():
            u = (r - r.min()) /(r.max() - r.min())
        else : u = r

        r = minr + (maxr-minr) * u
        x = scale*r * np.sin(th) * np.cos(phi)
        y = scale*r * np.sin(th) * np.sin(phi)
        z = scale*r * np.cos(th)
        if z.shape[1] != y.shape[1]:
            z = z*np.ones(y.shape[1])
        p = np.concatenate((x[...,None],
                            y[...,None],
                            z[...,None]),axis=2)
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

    def show3(self,k=0,po=[],T=[],txru=0,typ='G', mode='linear', silent=False):
        """ show3 geomview

        Parameters
        ----------

        k : frequency index
        po : poition of the antenna
        T  : GCS of the antenna
        typ : string
            'G' | 'Ft' | 'Fp'
        mode : string
            'linear'| 'not implemented'
        silent : boolean
            True    | False

        Examples
        --------

            >>> from pylayers.antprop.antenna import *
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> A = Antenna('defant.sh3')
            >>> #A.show3()

        """

        if not self.evaluated:
            self.eval(pattern=True)

        f = self.fGHz[k]

        # 3 axis : nth x nph x nf
        if len(self.Ft.shape)==3:
            if typ == 'G':
                V = self.sqG[:, :,k]
            if typ == 'Ft':
                V = self.Ft[:, :,k]
            if typ == 'Fp':
                V = self.Fp[:, :,k]
            if typ == 'Ft':
                V = self.Ft[:,:,k]

        # 4 axis : nth x nph x ntxru x nf
        if len(self.Ft.shape)==4:
            if typ == 'G':
                V = self.sqG[:, :, txru,k]
            if typ == 'Ft':
                V = self.Ft[:, : ,txru,k]
            if typ == 'Fp':
                V = self.Fp[:, :,txru,k]

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
        #    theta = self.theta[:,None]
        #else:
        #    theta=self.theta
        theta = self.theta
        #if len(np.shape(self.phi))==1:
        #    phi = self.phi[None,:]
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

        typ =  'Gain'
            = 'Ftheta'
            = 'Fphi'

        if col  -> color coded plot3D
        else    -> simple plot3D

        """

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)

        if typ == 'Gain':
            V = self.sqG[:, :,k]
        if typ == 'Ftheta':
            V = self.Ft[ :, :,k]
        if typ == 'Fphi':
            V = self.Fp[ :, :,k]

        vt = np.ones(self.nth)
        vp = np.ones(self.nph)
        Th = np.outer(self.theta, vp)
        Ph = np.outer(vt, self.phi)

        pdb.set_trace()
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

    def pol3d(self, k=0, R=50, St=4, Sp=4, silent=False):
        """ Display polarisation diagram  in 3D

           Parameters
           ----------

           k  : int
               frequency index
           R  : float
               radius of the sphere
           St : int
               downsampling factor along theta
           Sp : int
               downsampling factor along phi
           silent : Boolean
               (if True the file is created and not displayed')

           The file created is named : Polar{ifreq}.list
           it is placed in the /geom directory of the project

        """
        _filename = 'Polar' + str(10000 + k)[1:] + '.list'
        filename = pyu.getlong(_filename, pstruc['DIRGEOM'])
        fd = open(filename, "w")
        fd.write("LIST\n")

        Nt = self.nth
        Np = self.nph
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
                geu.ellipse(fd, p, B[0, :], B[1, :], self.Ft[n, m , k], self.Fp[n, m , k], N)
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

        sh = np.shape(self.Ft)
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

        Ftho = self.Ft[:, pt, :].ravel()
        Fpho = self.Fp[:, pt, :].ravel()
        SqGo = self.sqG[:, pt, :].ravel()

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

        delayCandidates : ndarray dalay in (ns)
            default np.arange(-10,10,0.001)

        Returns
        -------

        electricalDelay  : float

        Author : Troels Pedersen (Aalborg University)
                 B.Uguen

        """
        if self.evaluated:
            maxPowerInd  = np.unravel_index(np.argmax(abs(self.Ft)),np.shape(self.Ft))
            elD  = delayCandidates[np.argmax(abs(
                np.dot(self.Ft[maxPowerInd[0],maxPowerInd[1],:]
                       ,np.exp(2j*np.pi*self.fGHz[:,None]
                               *delayCandidates[None,:]))))]
            #electricalDelay  = delayCandidates[np.argmax(abs(
            #    np.dot(self.Ft[:,maxPowerInd[1],maxPowerInd[2]]
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
            >>> A.eval()
            >>> tau = A.getdelay()
            >>> A.elec_delay(tau)



        """

        self.tau = self.tau+tau
        if self.evaluated:
            Ftheta = self.Ft
            Fphi = self.Fp
            sh = np.shape(Ftheta)
            e = np.exp(2 * np.pi * 1j * self.fGHz[None,None,:]* tau)
            #E = np.outer(e, ones(sh[1] * sh[2]))
            #Fth = Ftheta.reshape(sh[0], sh[1] * sh[2])
            #EFth = Fth * E
            #self.Ft = EFth.reshape(sh[0], sh[1], sh[2])
            self.Ft = self.Ft*e
            self.Fp = self.Fp*e
            #Fph = Fphi.reshape(sh[0], sh[1] * sh[2])
            #EFph = Fph * E
            #self.Fp = EFph.reshape(sh[0], sh[1], sh[2])
        else:
            raise Warning('antenna has not been evaluated')


    def Fsynth(self,theta=[],phi=[],):
        """ Perform Antenna synthesis

        Parameters
        ----------

        theta : np.array
        phi :   np.array
            call Antenna.Fpatt or Antenna.Fsynth3

        Notes
        -----

        The antenna pattern synthesis is done either from spherical
        harmonics coefficients or from an analytical expression of the
        radiation pattern.

        """
        if ((self.fromfile) or (self.typ=='vsh') or (self.typ=='ssh')):
            Ft,Fp = self.Fsynth3(theta,phi)
            self.gain()
            self.evaluated=True
        else :
            Ft = self.Ft
            Fp = self.Fp
            self.theta = theta
            self.phi = phi
            eval('self.p'+self.typ)()
            #Ft,Fp = self.Fpatt(theta,phi,pattern)
        return (Ft,Fp)


    #def Fsynth1(self, theta, phi, k=0):

    def Fsynth1(self, theta, phi):
        """ calculate complex antenna pattern  from VSH Coefficients (shape 1)

        Parameters
        ----------

        theta  : ndarray (1xNdir)
        phi    : ndarray (1xNdir)
        k      : int
            frequency index

        Returns
        -------

        Ft , Fp 

        """

        Nt = len(theta)
        Np = len(phi)

        if self.grid:
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

        if self.grid:
            Nf = len(self.fGHz)
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

    def Fsynth2b(self, theta, phi):
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

        if self.grid:
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

        if self.grid:
            Nf = len(self.fGHz)
            Fth = Fth.reshape(Nf, Nt, Np)
            Fph = Fph.reshape(Nf, Nt, Np)

        return Fth, Fph

    def Fsynth2(self, theta, phi, typ = 'vsh'):
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

        self.nth = len(theta)
        self.nph = len(phi)
        self.nf = len(self.fGHz)

        if typ =='vsh' :

            if self.grid:
                theta = np.kron(theta, np.ones(self.nph))
                phi = np.kron(np.ones(self.nth),phi)

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

            if self.grid:
                Fth = Fth.reshape(self.nf, self.nth, self.nph)
                Fph = Fph.reshape(self.nf, self.nth, self.nph)

        if typ=='ssh':
            cx = self.S.Cx.s2
            cy = self.S.Cy.s2
            cz = self.S.Cz.s2

            lmax = self.S.Cx.lmax
            Y ,indx = SSHFunc(lmax, theta,phi)
            Ex = np.dot(cx,Y).reshape(self.nf,self.nth,self.nph)
            Ey = np.dot(cy,Y).reshape(self.nf,self.nth,self.nph)
            Ez = np.dot(cz,Y).reshape(self.nf,self.nth,self.nph)

            Fth,Fph = CartToSphere (theta, phi, Ex, Ey,Ez, bfreq = True )

        self.evaluated = True
        return Fth, Fph


    def Fsynth3(self,theta=[],phi=[],typ='vsh'):
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

        if self.grid:
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
            >>> F = A.eval(grid=True)

        All Br,Cr,Bi,Ci have the same (l,m) index in order to evaluate only
        once the V,W function

        If the data comes from a cst file like the antenna used in WHERE1 D4.1
        the pattern is multiplied by $\frac{4\pi}{120\pi}=\frac{1}{\sqrt{30}$

        """

        #typ = self.typ
        #self._filename.split('.')[1]
        #if typ=='satimo':
        #    coeff=1.
        #if typ=='cst':
        #    coeff=1./sqrt(30)


        #assert typ in ['ssh','vsh','hfss'], 
        assert (hasattr(self,'C') or hasattr(self,'S')),"No SH coeffs evaluated"

        Nf = len(self.fGHz)
        if theta==[]:
            theta=np.linspace(0,np.pi,45)

        if phi == []:
            phi= np.linspace(0,2*np.pi,90)

        Nt = len(theta)
        Np = len(phi)
        self.nth = len(theta)
        self.nph = len(phi)

        if self.grid:
            #self.theta = theta[:,None]
            #self.phi = phi[None,:]
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

            if self.grid:

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

            #self.Fp = Fph
            #self.Ft = Fth
            #G = np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth))
            #self.sqG = np.sqrt(G)


        #if self.grid:
        #    self.Fp = Fph
        #    self.Ft = Fth
        #    G = np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth))
        #    self.sqG = np.sqrt(G)

        self.evaluated = True

        #if typ == 'hfss':
        #    scipy.interpolate.griddata()

        #    Fth = self.Ft
        #    Fph = self.Fp
        # TODO create 2 different functions for pattern and not pattern
        #if not self.grid:
        return Fth, Fph
        #else:
        #    return None,None

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

        # print(Brmin, Brmax, Bimin, Bimax, Crmin, Crmax, Cimin, Cimax)

        for k in range(self.nf):
            plt.figure()
            stf = ' f=' + str(self.fGHz[k]) + ' GHz'
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

        #th = np.kron(self.theta, np.ones(self.nph))
        #ph = np.kron(np.ones(self.nth), self.phi)

        if not self.grid:
            self.grid = True
        Fth3, Fph3 = self.Fsynth3(self.theta, self.phi)
        Err = self.mse(Fth3, Fph3, 0)

        Enc = self.C.ens3()
        n = len(Enc)
        pos = 0

        while (pos < n) & (Err[0] < emax):

            Emin = Enc[pos]
            d = self.C.drag3(Emin)
            Fth3, Fph3 = self.Fsynth3(self.theta, self.phi)
            Err = self.mse(Fth3, Fph3, 0)

            if Err[0] >= emax:
                i = d[0][0]
                i3 = d[1][0]
                self.C.put3(i, i3)
                Fth3, Fph3 = self.Fsynth3(self.theta,self.phi)
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
            print( filevsh3, ' already exist')
        else:
            print( 'create ', filevsh3, ' file')

            coeff = {}
            coeff['fmin'] = self.fGHz[0]
            coeff['fmax'] = self.fGHz[-1]
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
        #typ = self._filename.split('.')[1]
        #self.typ = typ

        _filesh2 = self._filename.replace('.'+ self.typ, '.sh2')
        filesh2 = pyu.getlong(_filesh2, pstruc['DIRANT'])
        if os.path.isfile(filesh2):
            print(filesh2, ' already exist')
        else:
            print('create ', filesh2, ' file')
            coeff = {}
            coeff['fmin'] = self.fGHz[0]
            coeff['fmax'] = self.fGHz[-1]


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

        create a .sh3 antenna file

        """
        # create sh3 file
        # if self._filename has an extension 
        # it is replace by .sh3
        #typ = self._filename.split('.')[1]
        #self.typ = typ
        _filesh3 = self._filename.replace('.'+ self.typ, '.sh3')
        filesh3 = pyu.getlong(_filesh3, pstruc['DIRANT'])
        if os.path.isfile(filesh3):
            print(filesh3, ' already exist')
        else:
            print('create ', filesh3, ' file')

            coeff = {}
            coeff['fmin'] = self.fGHz[0]
            coeff['fmax'] = self.fGHz[-1]
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

            vsh3 file contains a thresholded version of vsh coefficients in shape 3

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
            self.nf = np.shape(Br.s3)[0]
            self.fGHz = np.linspace(fmin, fmax, self.nf)
        else:
            print(_filevsh3, ' does not exist')

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

            self.nf = np.shape(Cx.s3)[0]
            self.fGHz = np.linspace(fmin, fmax, self.nf)
        else:
            print(_filesh3, ' does not exist')

    def savevsh2(self, filename = ''):
        """ save coeff in  a .vsh2 antenna file

        Parameters
        ----------

        filename : string 

        """

        # create vsh2 file
        if filename == '':
            _filevsh2 = self._filename.replace('.trx', '.vsh2')

        _filevsh2  = filename
        filevsh2 = pyu.getlong(_filevsh2, pstruc['DIRANT'])

        if os.path.isfile(filevsh2):
            print(filevsh2, ' already exist')
        else:
            print('create ', filevsh2, ' file')

            coeff = {}
            coeff['fmin'] = self.fGHz[0]
            coeff['fmax'] = self.fGHz[-1]
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
            self.fGHz = np.linspace(fmin, fmax, Nf)
        else:
            print( _filesh2, ' does not exist')



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
            self.fGHz = np.linspace(fmin, fmax, Nf)
        else:
            print( _filevsh2, ' does not exist')

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
            self.fGHz = np.linspace(fmin, fmax, 121)
        else:
            print(_filevsh3, ' does not exist')

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
        Fth = self.Ft[:, ith, :]
        Fph = self.Fp[:, ith, :]
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

        self.sqG[:, ith, :] = SqG
        self.Ft[:, ith, :] = Fth
        self.Fp[:, ith, :] = Fph

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
        print("Error: m>n!!!")

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
    maxTr = abs(Fthr[:, :, k]).max()
    maxTo = abs(Ftho[:, :, k ]).max()
    MmT = max(maxTr, maxTo)

    minTr = abs(Fthr[ :, :, k ]).min()
    minTo = abs(Ftho[ :, :, k ]).min()
    mmT = min(minTr, minTo)

    maxPr = abs(Fphr[ :, :, k ]).max()
    maxPo = abs(Fpho[ :, :, k ]).max()
    MmP = max(maxPr, maxPo)

    minPr = abs(Fphr[ :, :, k ]).min()
    minPo = abs(Fpho[ :, :, k ]).min()
    mmP = min(minPr, minPo)

    # limites real Fthr, Ftho, Fphr, Fpho
    maxTrr = np.real(Fthr[ :, :, k ]).max()
    maxTor = np.real(Ftho[ :, :, k ]).max()
    MrT = max(maxTrr, maxTor)

    minTrr = np.real(Fthr[ :, :, k ]).min()
    minTor = np.real(Ftho[ :, :, k ]).min()
    mrT = min(minTrr, minTor)

    maxPrr = np.real(Fphr[ :, :, k ]).max()
    maxPor = np.real(Fpho[ :, :, k ]).max()
    MrP = max(maxPrr, maxPor)

    minPrr = np.real(Fphr[ :, :, k ]).min()
    minPor = np.real(Fpho[ :, :, k ]).min()
    mrP = min(minPrr, minPor)

    # limites real Fthr, Ftho, Fphr, Fpho
    maxTri = np.imag(Fthr[ :, :, k ]).max()
    maxToi = np.imag(Ftho[ :, :, k ]).max()
    MiT = max(maxTri, maxToi)

    minTri = np.imag(Fthr[ :, :, k ]).min()
    minToi = np.imag(Ftho[ :, :, k ]).min()
    miT = min(minTri, minToi)

    maxPri = np.imag(Fphr[ :, :, k ]).max()
    maxPoi = np.imag(Fpho[ :, :, k ]).max()
    MiP = max(maxPri, maxPoi)

    minPri = np.imag(Fphr[ :, :, k ]).min()
    minPoi = np.imag(Fpho[ :, :, k ]).min()
    miP = min(minPri, minPoi)

    # limithes arg Fth,Fph
    maxATr = np.angle(Fthr[ :, :, k ]).max()
    maxATo = np.angle(Ftho[ :, :, k ]).max()
    maT = max(maxATr, maxATo)
    minATr = np.angle(Fthr[ :, :, k ]).min()
    minATo = np.angle(Ftho[ :, :, k ]).min()
    maT0 = min(minATr, minATo)

    maxAPr = np.angle(Fphr[ :, :, k ]).max()
    maxAPo = np.angle(Fpho[ :, :, k ]).max()
    maP = max(maxAPr, maxAPo)
    minAPr = np.angle(Fphr[ :, :, k ]).min()
    minAPo = np.angle(Fpho[ :, :, k ]).min()
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
        plt.pcolor(A.phi * rtd, A.theta * rtd, abs(Ftho[ :, :, k ]),
                   cmap=cm.hot_r, vmin=mmT, vmax=MmT)
        plt.title(r'$|F_{\theta}|$ original', fontsize=fontsize)

    if typ == 'real':
        #pcolor(A.phi*rtd,A.theta*rtd,real(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.real(Ftho[ :, :, k ]),
                   cmap=cm.hot_r, vmin=mrT, vmax=MrT)
        title(r'Re ($F_{\theta}$) original', fontsize=fontsize)
    if typ == 'imag':
        #pcolor(A.phi*rtd,A.theta*rtd,imag(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
        pcolor(A.phi * rtd, A.theta * rtd, np.imag(Ftho[ :, :, k ]),
               cmap=cm.hot_r, vmin=miT, vmax=MiT)
        title(r'Im ($F_{\theta}$) original', fontsize=fontsize)
    if typ == 'phase':
        #pcolor(A.phi*rtd,A.theta*rtd,angle(Ftho[k,:,:]),cmap=cm.gray_r,vmin=maT0,vmax=maT)
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Ftho[ :, :, k ]),
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
        plt.pcolor(A.phi * rtd, A.theta * rtd, abs(Fpho[:, :, k ]),
                   cmap=cm.hot_r, vmin=mmP, vmax=MmP)
        plt.title('$|F_{\phi}|$ original', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.real(Fpho[ :, :, k ]),
                   cmap=cm.hot_r, vmin=mrP, vmax=MrP)
        plt.title('Re ($F_{\phi}$) original', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.imag(Fpho[ :, :, k ]),
                   cmap=cm.hot_r, vmin=miP, vmax=MiP)
        plt.title('Im ($F_{\phi}$) original', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Fpho[ :, :, k ]),
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
        plt.pcolor(ph * rtd, th * rtd, abs(Fthr[:, :, k ]),
                   cmap=cm.hot_r, vmin=mmT, vmax=MmT)
        if lang == 'french':
            plt.title(r'$|F_{\theta}|$ reconstruit', fontsize=fontsize)
        else:
            plt.title(r'$|F_{\theta}|$ reconstructed', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(ph * rtd, th * rtd, np.real(Fthr[:,:,k ]),
                   cmap=cm.hot_r, vmin=mrT, vmax=MrT)
        if lang == 'french':
            title(r'Re ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            title(r'Re ($F_{\theta}$) reconstructed', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(ph * rtd, th * rtd, np.imag(Fthr[ :, :, k ]),
                   cmap=cm.hot_r, vmin=miT, vmax=MiT)
        if lang == 'french':
            plt.title(r'Im ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            plt.title(r'Im ($F_{\theta}$) reconstructed', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Fthr[:,:,k]),
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
        plt.pcolor(ph * rtd, th * rtd, abs(Fphr[ :, :,k]),
                   cmap=cm.hot_r, vmin=mmP, vmax=MmP)
        if lang == 'french':
            plt.title('$|F_{\phi}|$ reconstruit', fontsize=fontsize)
        else:
            plt.title('$|F_{\phi}|$ reconstructed', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(ph * rtd, th * rtd, np.real(Fphr[ :, :,k]),
                   cmap=cm.hot_r, vmin=mrP, vmax=MrP)
        if lang == 'french':
            plt.title('Re ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Re ($F_{\phi}$) reconstructed', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(ph * rtd, th * rtd, np.imag(Fphr[ :, :,k]),
                   cmap=cm.hot_r, vmin=miP, vmax=MiP)
        if lang == 'french':
            plt.title('Im ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Im ($F_{\phi}$) reconstructed', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.angle(Fphr[ :, :,k]),
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

def BeamGauss(theta,phi,Gmax=19.77,HPBW_az=10,HPBW_el=40,Tilt=10):
    """ Beam with a Gaussian shape

    Parameters
    ----------
                
    theta : float 
          angle in degree
    phi   : float 
          angle in degree
    Gmax  : float 
    HPBW_az : float  
        Half Power Beamwidth azimuth degree
    HPBW_el : float
        Half Power Beamwidth elevation degree
    Tilt : float 
        angle in degree 

    """
    c = np.pi/180.
    az = c*(theta-(Tilt+90))*2*np.sqrt(np.log(2))
    el = c*phi*2*np.sqrt(np.log(2))
    taz = -(az/(HPBW_az*c))**2
    tel = -(el/(HPBW_el*c))**2
    gain = 10**(Gmax/10.)*np.exp(taz)*np.exp(tel)
    return(gain)

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
        >>> A.eval(grid=True)


    Warnings
    --------

        len(theta) must be equal with shape(F)[1]
        len(phi) must be equal with shape(F)[2]

    """

    nth = len(theta)
    nph = len(phi)

    if k >= np.shape(F)[0]:
        print('Error: frequency index k not in F defined interval')
    if nth != np.shape(F)[1]:
        print('Error: shape mistmatch between theta and F')

    if nph != np.shape(F)[2]:
        print('Error: shape mistmatch between phi and F')

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

class AntPosRot(Antenna):
    """ Antenna + position + Rotation
    """
    def __init__(self,name,p,T):
        Antenna.__init__(self,name)
        self.p = p 
        self.T = T 
    def _show3(self,**kwargs):
        Antenna._show3(self,newfig=False,interact=False,T=self.T,po=self.p,**kwargs)

    def field(self,p):
        """
        Parameters
        ----------

        p : np.array (N,3)

        """
        rad_to_deg = 180/np.pi
        assert p.shape[-1]==3

        if len(p.shape)==1:
            r = p[None,:]-self.p[None,:]
        else:
            r = p-self.p[None,:]
        dist = np.sqrt(np.sum(r*r,axis=-1))[:,None]
        u = r/dist
        th = np.arccos(u[:,2])
        ph = np.arctan2(u[:,1],u[:,0])
        tang = np.vstack((th,ph)).T
        #print("global",tang*rad_to_deg)
        Rt, tangl = geu.BTB_tx(tang, self.T)
        #print("local",tangl*rad_to_deg)
        self.eval(th=tangl[:,0],ph=tangl[:,1],grid=False)
        E = (self.Ft[:,None,:]*self.T[:,2][None,:,None]+self.Fp[:,None,:]*self.T[:,0][None,:,None])
        P = np.exp(-1j*2*np.pi*self.fGHz[None,None,:]*dist[...,None]/0.3)/dist[...,None]
        EP = E*P
        return(EP)
        #Rr, rangl = geu.BTB_rx(rang, self.Tr)

def _gain(Ft,Fp):
    """  calculates antenna gain
    
    Returns
    -------

    G  : np.array(Nt,Np,Nf) dtype:float
        linear gain 
              or np.array(Nr,Nf)
    sqG : np.array(Nt,Np,Nf) dtype:float 
        linear sqare root of gain 
              or np.array(Nr,Nf)
    efficiency : np.array (,Nf) dtype:float 
        efficiency 
    hpster : np.array (,Nf) dtype:float
        half power solid angle :  1 ~ 4pi steradian 
    ehpbw : np.array (,Nf) dtyp:float 
        equivalent half power beamwidth (radians)

    Notes
    -----

    .. math:: G(\theta,phi) = |F_{\\theta}|^2 + |F_{\\phi}|^2

    """
    G = np.real( Fp * np.conj(Fp)
              +  Ft * np.conj(Ft) )

    return(G)


def _hpbw(G,th,ph):
    """ half power beamwidth 

    Parameters
    ----------
    Gain : Ftheta
        Nt x Np 
    th : np.array
        ,Nt
    ph : np.array
        ,Np

    Returns
    -------

    ehpbw : effective half power beamwidth 
    hpster : half power solid angle (steradians)

    """
    #
    GdB = 10*np.log10(G)
    GdBmax = np.max(np.max(GdB,axis=0),axis=0)
    
    dt = th[1]-th[0]
    dp = ph[1]-ph[0]
    Nt = len(th)
    Np = len(ph)
    Nf = GdB.shape[2]
    hpster = np.zeros(Nf)
    ehpbw  =  np.zeros(Nf)
    for k in range(Nf):
        U  = np.zeros((Nt,Np))
        A = GdB[:,:,k]*np.ones(Nt)[:,None]*np.ones(Np)[None,:]
        u = np.where(A>(GdBmax[k]-3))
        U[u] = 1
        V  = U*np.sin(th)[:,None]
        hpster[k]  = np.sum(V)*dt*dp/(4*np.pi)
        ehpbw[k]  = np.arccos(1-2*hpster[k])

    return ehpbw,hpster 

def _efficiency(G,th,ph):
    """ determine antenna efficiency 

    Parameters
    ----------
    Gain : Ftheta
        Nt x Np 
    th : np.array
        ,Nt
    ph : np.array
        ,Np

    Returns
    -------

    oefficiency : 

    """
    #
    dt = th[1]-th[0]
    dp = ph[1]-ph[0]
    Nt = len(th)
    Np = len(ph)
    Gs = G*np.sin(th)[:,None,None]*np.ones(Np)[None,:,None]
    efficiency = np.sum(np.sum(Gs,axis=0),axis=0)*dt*dp/(4*np.pi)
    return efficiency

def _dirmax(G,th,ph):
    """ determine information in Gmax direction 

    Parameters
    ----------
    Gain : Ftheta
        Nt x Np 
    th : np.array
        ,Nt
    # GdBmax (,Nf)
    # Get direction of Gmax and get the polarisation state in that direction 
    # 
    Returns
    --------
    """

    GdB = 10*np.log10(G)
    GdBmax = np.max(np.max(GdB,axis=0),axis=0)
    umax = np.array(np.where(GdB==GdBmax))[:,0]
    theta_max = th[umax[0]]
    phi_max = ph[umax[1]]
    M = geu.SphericalBasis(np.array([[theta_max,phi_max]]))
    sl = M[:,2].squeeze()
    uth = M[:,0] 
    uph = M[:,1] 
    el = Ft[tuple(umax)]*uth + Fp[tuple(umax)]*uph
    eln = el/np.linalg.norm(el)
    el = np.abs(eln.squeeze())
    hl = np.cross(sl,el)
    return GdBmax,theta_max,phi_max,(hl,sl,el)

def F0(nu,sigma):
    """ F0 function for horn antenna pattern 
    
    Parameters
    ----------

    nu : np.array 
        (....,nf)
    sigma : np.array 
        (,nf) 

    Notes
    -----
    http://www.ece.rutgers.edu/~orfanidi/ewa/ch18.pdf

    18.3.2

    """
    nuos  = nu/sigma
    argp = nuos + sigma
    argm = nuos - sigma
    expf = np.exp(1j*(np.pi/2)*nuos**2)
    sf   = 1./sigma
    sp , cp = fresnel(argp)
    sm , cm = fresnel(argm)
    Fp = cp-1j*sp
    Fm = cm-1j*sm

    F = sf*expf*(Fp -Fm)
    return F 

def F1(nu,sigma):
    """ F1 function for horn antenna pattern 
   
    http://www.ece.rutgers.edu/~orfanidi/ewa/ch18.pdf

    18.3.3

    """
    F = 0.5*(F0(nu+0.5,sigma)+F0(nu-0.5,sigma))
    return F

if (__name__ == "__main__"):
    doctest.testmod()
