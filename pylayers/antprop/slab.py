#!/usr/bin/python
# -*- coding: latin1 -*-
"""
This module implements the Material and Slab class which are used 
in the ray tracing engine rays.py. 
It exploits heavily numpy broadcasting mechanism :


The adopted axis convention is the following

+ nf : axis = 0 frequency axis
+ nt : axis = 1 angular axis
+ p  : axis = 2 parallel polarization axis
+ o  : axis = 3 orhogonal polarization axis


.. currentmodule:: pylayers.antprop.slab

.. autosummary::
    :toctree: generated/

Interface Class
===============

.. autosummary::
    :toctree: generated/

    Interface.__init__
    Interface.RT
    Interface.pcolor
    Interface.tocolor
    Interface.loss0
    Interface.losst
    Interface.plotwrt

MatInterface Class
==================

.. autosummary::
    :toctree: generated/

     MatInterface.__init__


Mat Class
=========

.. autosummary::
    :toctree: generated/

    Mat.__init__
    Mat.eval
    Mat.info
    Mat.R


MatDB Class
===========

.. autosummary::
    :toctree: generated/

    MatDB.__init__
    MatDB.info
    MatDB.dass
    MatDB.maxindex
    MatDB.delete
    MatDB.edit
    MatDB.add
    MatDB.addgui
    MatDB.choose
    MatDB.load
    MatDB.save
    MatDB.savemat

Slab Class
==========

.. autosummary::
    :toctree: generated/

    Slab.__init__
    Slab.__repr__
    Slab.info
    Slab.conv
    Slab.ev
    Slab.filter
    Slab.excess_grdelay
    Slab.tocolor
    Slab.loss0
    Slab.losst
    Slab.editgui
    Slab.show

SlabDB Class
============

.. autosummary::
    :toctree: generated/

    SlabDB.__init__
    SlabDB.showall
    SlabDB.info
    SlabDB.dass
    SlabDB.maxindex
    SlabDB.delete
    SlabDB.edit
    SlabDB.show
    SlabDB.add
    SlabDB.addgui
    SlabDB.load
    SlabDB.loadsl
    SlabDB.save
    SlabDB.savesl

Utility Functions
==================

.. autosummary::
    :toctree: generated/

    calsig

"""
from __future__ import division, print_function, absolute_import 
import os
import sys
import string
if sys.version_info.major==2:
    import cPickle
else:
    import _pickle as cPickle
import doctest
#import objxml
import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import struct as stru
import ConfigParser
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.plotutil as plu
from pylayers.util.easygui import *
from scipy.interpolate import interp1d
import copy
import pdb
import copy


class Interface(PyLayers):
    """ Interface between 2 medium

    Attributes
    ----------

    fGHz np.array (nf,1)
    theta np.array (1,nt)
    Ip np.array (nf,nt,2,2)
    Io np.array (nf,nt,2,2)

    """
    def __init__(self, fGHz=np.array([2.4]), theta=np.array([[0.0 + 0 * 1j]]), name=''):
        """ class constructor

        Parameters
        ----------
        fGHz : np.array
            frequency in GHz (default 2.4)
        theta : np.array
        angle taken from surface normal expressed in radians

        """
        #
        # reshape theta if necessary
        # allows to use a ndim = 1 array
        #
        if theta.ndim != 2:
            theta = theta.reshape(1, len(theta))
        self.nf = len(fGHz)
        self.nt = np.shape(theta)[1]
        # f x 1
        self.fGHz = fGHz.reshape(self.nf, 1)
        # 1 x q
        self.thi = theta
        self.theta = theta
        #
        # Io : f x q x 2 x 2
        # Ip : f x q x 2 x 2
        #
        self.Ip = np.array(np.zeros([self.nf, self.nt, 2, 2]), dtype=complex)

        self.Io = np.array(np.zeros([self.nf, self.nt, 2, 2]), dtype=complex)

        self.name = name

    def RT(self, metalic=False, RT='RT'):
        r""" evaluate Reflection and Transmission matrix

        Parameters
        ----------
        metalic : boolean
        RT : string
            choose R or T

        Notes
        -----

        .. math::

            R = \left[\begin{array}[cc](R_o & 0\\0 & R_p)\end{array}\right]
            T = \left[\begin{array}[cc](T_o & 0\\0 & T_p)\end{array}\right]

        R : np.array (f , th , 2, 2)
        T : np.array (f , th , 2, 2)

        """
        sh = np.shape(self.Io)
        nf = sh[0]
        nt = sh[1]

        #
        # R and T matrices are diagonal
        #

        if 'R' in RT:
            self.R = np.array(np.zeros([nf, nt, 2, 2]), dtype=complex)
        if 'T' in RT:
            self.T = np.array(np.zeros([nf, nt, 2, 2]), dtype=complex)


        if 'R' in RT:
            self.R[:, :, 0, 0] = self.Io[:, :, 0, 1] / self.Io[:, :, 0, 0]
            self.R[:, :, 1, 1] = self.Ip[:, :, 0, 1] / self.Ip[:, :, 0, 0]

        if not metalic and 'T' in RT:
            self.T[:, :, 0, 0] = 1.0 / self.Io[:, :, 0, 0]
            self.T[:, :, 1, 1] = 1.0 / self.Ip[:, :, 0, 0]


    def pcolor(self, dB=False):
        """ display of R & T coefficients wrt frequency an angle

        Parameters
        ----------

        dB : boolean
        default False


        """
        rtd = 180 / np.pi
        #nom = self.m1.name+'|'+self.m2.name
        if dB:
            modRo = 20 * np.log10(abs(self.R[:, :, 0, 0]))
            modRp = 20 * np.log10(abs(self.R[:, :, 1, 1]))
            modTo = 20 * np.log10(abs(self.T[:, :, 0, 0]))
            modTp = 20 * np.log10(abs(self.T[:, :, 1, 1]))
        else:
            modRo = abs(self.R[:, :, 0, 0])
            modRp = abs(self.R[:, :, 1, 1])
            modTo = abs(self.T[:, :, 0, 0])
            modTp = abs(self.T[:, :, 1, 1])
        plt.hot()
        plt.subplot(221)
        plt.pcolor(self.theta[0, :] * rtd, self.fGHz[:, 0], modRo)
        plt.colorbar()
        plt.contour(self.theta[0, :] * rtd, self.fGHz[:, 0], modRo)
        plt.xlabel(r'$\theta$ (degrees) ')
        plt.ylabel('f (GHz)')
        #title('R _|_ '+nom)
        plt.title('R _|_ ')
        plt.subplot(222)
        plt.pcolor(self.theta[0, :] * rtd, self.fGHz[:, 0], modRp)
        plt.colorbar()
        plt.contour(self.theta[0, :] * rtd, self.fGHz[:, 0], modRp)
        #title('R // '+nom)
        plt.title('R // ')
        plt.xlabel(r'$\theta$ (degrees) ')
        plt.ylabel('f (GHz)')
        plt.subplot(223)
        plt.pcolor(self.theta[0, :] * rtd, self.fGHz[:, 0], modTo)
        plt.colorbar()
        plt.contour(self.theta[0, :] * rtd, self.fGHz[:, 0], modTo)
        #title('T _|_ '+nom)
        plt.title('T _|_ ')
        plt.xlabel(r'$\theta$ (degrees) ')
        plt.ylabel('f (GHz)')
        plt.subplot(224)
        plt.pcolor(self.theta[0, :] * rtd, self.fGHz[:, 0], modTp)
        plt.colorbar()
        plt.contour(self.theta[0, :] * rtd, self.fGHz[:, 0], modTp)
        #title('T // '+nom)
        plt.title('T // ')
        plt.xlabel(r'$\theta$ (degrees) ')
        plt.ylabel('f (GHz)')
        plt.show()

    def tocolor(self,fGHz):
        """ convert transmission into color

        Parameters
        ----------
        fGHz : np.array

        Returns
        -------
        col : string
            hexadecimal color


        See Also
        ---------

        pylayers.gis.layout.showGs

        """
        # nf x nt x 2 x 2
        modTo = abs(self.T[:, 0, 0, 0])
        modTp = abs(self.T[:, 0, 1, 1])
        N = len(fGHz)

        if N>3:
            M = N/3

            ared = (sum(modTo[0:M])+sum(modTp[0:M]))/(2*M)
            agreen = (sum(modTo[M:2*M])+sum(modTp[M:2*M]))/(2*M)
            ablue = (sum(modTo[2*M:])+sum(modTp[2*M:]))/(2*(N-2*M))

            # hexadsecimal convert
            vred = hex(int(np.floor(ared*255))).replace('0x','')
            vgreen = hex(int(np.floor(agreen*255))).replace('0x','')
            vblue = hex(int(np.floor(ablue*255))).replace('0x','')

            if len(vred)==1:
                vred = '0'+vred
            if len(vgreen)==1:
                vgreen = '0'+vgreen
            if len(vblue)==1:
                vblue = '0'+vblue
            col = '#'+vred+vgreen+vblue
        else:
            alpha = (sum(modTo)+sum(modTp))/(2*N)
            val = hex(int(np.floor(alpha*255))).replace('0x','')
            col = '#'+val+val+val
        return(col)


    def loss0(self, fGHz, display=False):
        """ evaluate Loss at normal incidence theta=0

        Parameters
        ----------

        fGHz : np.array (nf,1)
        display : boolean
        default (False)

        Returns
        -------

        Lo : loss in dB polarization orthogonal
        Lp : loss in dB polarization parallel

        """

        modTo = abs(self.T[:, :, 0, 0])
        modTp = abs(self.T[:, :, 1, 1])

        Lo = -20 * np.log10(modTo[:, 0])
        Lp = -20 * np.log10(modTp[:, 0])

        if display:
            plt.plot(f, Lo, 'b')
            #plot(f,Lp,'r')
            #legend(('L0 _|_','L0 //'),loc='upper right')
            plt.legend(('L0 (dB)'), loc='upper right')
            plt.xlabel('frequency (GHz)')
            plt.show()

        return(Lo, Lp)

    def losst(self, fGHz, display=False,dB=True):
        """ evaluate Loss

        Parameters
        ----------

        fGHz : np.arrray (nf)
        display : boolean
            default False
        dB : booean
            default True

        Returns
        -------

        Lo : np.array
            Loss orthogonal polarization (dB)
        Lp : np.array
            Loss parallel polarization (dB)

        Examples
        --------

        >>> from pylayers.antprop.slab import *

        See Also
        --------

        pylayers.antprop.coverage
        pylayers.antprop.loss

        """

        modTo = abs(self.T[:, :, 0, 0])
        modTp = abs(self.T[:, :, 1, 1])

        #Lo = -20 * np.log10(modTo[:, 0])
        #Lp = -20 * np.log10(modTp[:, 0])

        if dB:
            Lo = -20 * np.log10(modTo)
            Lp = -20 * np.log10(modTp)
        else:
            Lo = modTo
            Lp = modTp

        if display:
            plt.plot(fGHz, Lo, 'b')
            plt.legend(('L0 (dB)'), loc='upper right')
            plt.xlabel('frequency (GHz)')
            plt.show()

        return(Lo, Lp)

    def plotwrt(self,var='a',kv=0,**kwargs):
        """ plot R & T coefficients with respect to angle or frequency

        Parameters
        ----------

        kv : int
        variable index
        polar: string
        'po', # po | p | o (parallel+ortho | parallel | ortogonal)
        coeff: string
        'RT', # RT | R | T (Reflexion & Transmission ) | Reflexion | Transmission
        var: string
        'a', # a | f angle | frequency
        typ : string
        'm' | 'r' | 'd' | 'l20'
        mod rad deg dB

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.antprop.slab import *
            >>> import matplotlib.pylab as plt
            >>> import numpy as np
            >>> theta = np.arange(0,np.pi/2,0.01)
            >>> fGHz = np.arange(0.1,10,0.2)
            >>> sl = SlabDB('matDB.ini','slabDB.ini')
            >>> mat = sl.mat
            >>> lmatname = [mat['AIR'],mat['WOOD']]
            >>> II = MatInterface(lmatname,0,fGHz,theta)
            >>> II.RT()
            >>> fig,ax = II.plotwrt(var='a',kv=10,typ=['m'])
            >>> air = mat['AIR']
            >>> brick = mat['BRICK']
            >>> II = MatInterface([air,brick],0,fGHz,theta)
            >>> II.RT()
            >>> fig,ax = II.plotwrt(var='f',color='k',typ=['m'])
            >>> plt.ion()
            >>> plt.show()


        """
        defaults = {'typ':['l20'],
                'polar':'po', # po | p | o
                'coeff':'RT', # RT | R | T
                'att':False
               }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        #fGHz = self.fGHz[k]
        rtd = 180 / np.pi

        # filtering kwargs argument for mulplot function
        args ={}
        for k in kwargs:
            if k not in defaults.keys():
                args[k]=kwargs[k]

        if 'labels' not in kwargs.keys():
            args['labels'] = [self.name]

        args['titles'] = []
        args['typ'] = kwargs['typ']
        args['att'] = kwargs['att']

        # Reflexion
        if 'R' in kwargs['coeff']:
            if 'o' in kwargs['polar']:
                args['titles'].append(u'$R_{\perp}$')
                if var=='f': # wrt frequency
                    Ro = self.R[:, kv, 0, 0]
                    y = Ro
                if var=='a': # wrt angle
                    Ro = self.R[kv, :, 0, 0]
                    y = Ro
            if 'p' in kwargs['polar']:
                args['titles'].append(u'$R_{//}$')
                if var=='f': # wrt frequency
                    Rp = self.R[:, kv, 1, 1]
                    try:
                        y = np.vstack((y,Rp))
                    except:
                        y = Rp
                if var =='a': # wrt angle
                    Rp = self.R[kv, :, 1, 1]
                    try:
                        y = np.vstack((y,Rp))
                    except:
                        y = Rp
        # Transmission
        if 'T' in kwargs['coeff']:
            if 'o' in kwargs['polar']:
                args['titles'].append(u'$T_{\perp}$')
                if var=='f': # wrt frequency
                    To = self.T[:, kv, 0, 0]
                    try:
                        y = np.vstack((y,To))
                    except:
                        y = To
                if var =='a': # wrt angle
                    To = self.T[kv, :, 0, 0]
                    try:
                        y = np.vstack((y,To))
                    except:
                        y = To
            if 'p' in kwargs['polar']:
                args['titles'].append(u'$T_{//}$')
                if var=='f': # wrt frequency
                    Tp = self.T[:, kv, 1, 1]
                    try:
                        y = np.vstack((y,Tp))
                    except:
                        y = To
                if var =='a': # wrt angle
                    Tp = self.T[kv, :, 1, 1]
                    try:
                        y = np.vstack((y,Tp))
                    except:
                        y = To

        # setting the x axis
        if var=='f': # wrt frequency
            if len(self.fGHz)==1:
                #x = self.fGHz[np.newaxis,:]
                x = self.fGHz[:]
            else: # f x a
                x = self.fGHz[:,0]
            args['xlabels'] = ['Frequency (GHz)']
        if var=='a': # wrt angle
            if len(self.thi)==1:
                x = self.thi[0,:][:]*rtd
                #x = self.thi[0,:][np.newaxis,:]*rtd
            else: # f x a
                x = self.thi[0,:]*rtd
            args['xlabels'] = ['Angle (deg)']

        nplot = np.shape(y)[0]
        if nplot==1:
            args['ncol'] = 1
            args['nlin'] = 1
        if nplot==2:
            args['ncol'] = 1
            args['nlin'] = 2
        if nplot==4:
            args['ncol'] = 2
            args['nlin'] = 2

        fig,ax = plu.mulcplot(x,y,**args)

        return fig,ax


class MatInterface(Interface):
    r""" MatInterface : Class for Interface between two materials

    l distance from the next Interface

    Notes
    -----

    This is required for recursive utilization of this function when the
    output angle of an interface happens to be the input angle of the
    next interface. As this angle depends on materials which themselves
    depends on frequency THETA is becoming a full matrix without redundancy
    between lines.

    >>> theta = np.arange(0,np.pi/2,0.01)
    >>> fGHz = np.arange(3.1,10.6,0.2)
    >>> Nf = len(fGHz)
    >>> Nt = len(theta)
    >>> sl = SlabDB('matDB.ini','slabDB.ini')
    >>> mat = sl.mat
    >>> m1 = mat['AIR']
    >>> m2 = mat['PLASTER']
    >>> II = MatInterface([m1,m2],0,fGHz,theta)

    .. math::

    I_p = \left| \begin{array}{cc} \frac{1}{T_p} & \frac{R_p}{T_p} \\ \frac{R_p}{T_p} & \frac{1}{T_p} \end{array}\right|

    I_o = \left| \begin{array}{cc} \frac{1}{T_o} & \frac{R_o}{T_o} \\ \frac{R_o}{T_o} & \frac{1}{T_o} \end{array}\right|


    .. todo::
    MatIinterface fix the np.pi/2 NaN problem


    """
    def __init__(self, lmat, l, fGHz, theta):
        """ Fresnel reflection coefficients

        Parameters
        ----------

        lmat : [ m1 , m2] list of materials
        l : distance between interfaces
        fGHz : frequency (GHz)
        theta : angle with respect to the reflective surface normal

        """
        #if not isinstance(fGHz,np.ndarray):
        # fGHz=np.array([fGHz])
        #if not isinstance(theta,np.ndarray):
        # theta=np.array([theta])
        name = '|'.join(mat['name'] for mat in lmat)
        # Interface.__init__(self, fGHz, theta, name=name)
        super(MatInterface,self).__init__(fGHz, theta, name=name)
        self.m1 = lmat[0]
        self.m2 = lmat[1]
        # 2*np.pi* f(GHz)*eps0 = f(Ghz)/18

        epr1 = self.m1['epr'] - 1j * abs(self.m1['sigma']) * 18 / self.fGHz
        epr2 = self.m2['epr'] - 1j * abs(self.m2['sigma']) * 18 / self.fGHz

        mur1 = self.m1['mur']
        mur2 = self.m2['mur']

        #n1 = sqrte(epr1/mur1)

        n1 = np.sqrt(epr1 / mur1)
        n2 = np.sqrt(epr2 / mur2)

        ct = np.cos(self.theta)
        # // TM polarization 8.1.4 http://www.ece.rutgers.edu/~orfanidi/ewa/ch08.pdf
        nT1p = n1 / ct
        # _|_ TE polarization 8.1.4 http://www.ece.rutgers.edu/~orfanidi/ewa/ch08.pdf
        nT1o = n1 * ct

        #print np.shape(n1)
        #print np.shape(ct)
        #print "Slab cst n1 et ct ",n1[15,0],ct[0,31]
        #print "Slab cst nT1p ",nT1p[15,31]
        #print "Slab cst nT1o ",nT1o[15,31]

        #cti = pyu.sqrte(1-((n1/n2)*np.sin(self.theta))**2)
        cti = np.sqrt(1 - ((n1 / n2) * np.sin(self.theta)) ** 2)
        #CTI = np.sqrt(1-((n1/n2)*np.sin(THETA))**2)
        self.theta = np.arccos(cti)
        #print np.shape(cti)
        #print "cti ",cti[15,31]
        #print "arcos(cti) ",self.theta[15,31]
        #print '-------------------------'

        if l != 0:
            deltai = 2 * np.pi * l * n2 * cti * self.fGHz / 0.3
        else:
            deltai = 0

        nT2p = n2 / cti
        nT2o = n2 * cti

        Rp = (nT1p - nT2p) / (nT1p + nT2p)
        #Ro = (nT1o-nT2o)/(nT1o+nT2o)
        Ro = -(nT1o - nT2o) / (nT1o + nT2o) # modif Eric
        Tp = 1.0 + Rp
        To = 1.0 + Ro

        self.Ro = Ro
        self.Rp = Rp

        jdeltai = 1j*deltai
        epd = np.exp(jdeltai)
        emd = np.exp(-jdeltai)

        epdoTp = epd/Tp
        emdoTp = emd/Tp

        self.Ip[:, :, 0, 0] = epdoTp
        self.Ip[:, :, 0, 1] = Rp * emdoTp
        self.Ip[:, :, 1, 0] = Rp * epdoTp
        self.Ip[:, :, 1, 1] = emdoTp

        epdoTo = epd/To
        emdoTo = emd/To

        self.Io[:, :, 0, 0] = epdoTo
        self.Io[:, :, 0, 1] = Ro * emdoTo
        self.Io[:, :, 1, 0] = Ro * epdoTo
        self.Io[:, :, 1, 1] = emdoTo

        #print 'Slab MatInterface Ip00',self.Ip[15,31,0,0]
        #print 'Slab MatInterface Ip01',self.Ip[15,31,0,1]
        #print 'Slab MatInterface Ip10',self.Ip[15,31,1,0]
        #print 'Slab MatInterface Ip11',self.Ip[15,31,1,1]

        #print 'Slab MatInterface Io00',self.Io[15,31,0,0]
        #print 'Slab MatInterface Io01',self.Io[15,31,0,1]
        #print 'Slab MatInterface Io10',self.Io[15,31,1,0]
        #print 'Slab MatInterface Io11',self.Io[15,31,1,1]


class Mat(PyLayers,dict):
    """ Handle constitutive materials dictionnary

    Attributes
    ----------

    name : string
    name character string (default 'AIR')
    index : int
    default 1
    er : complex
    relative permittivity (w.d) (1+0j)
    mur : complex
    relative permeability (w.d) (1+0j)
    sigma : float
    conductivity (S/m) 0
    roughness : float
    (meter) 0

    """
    def __init__(self, 
                 name="AIR", 
                 index=1, 
                 epr=1 + 0.0j, 
                 mur=1 + 0.0j, 
                 sigma=0.0,
                 roughness=0.,
                 a=None,
                 b=None,
                 c=None,
                 d=None):
        """ class constructor

        Parameters
        ----------

        name : string
        index : int
        epr : complex
        mur : complex
        sigma : float
        roughness : float
        a : ITU parameter
        b : ITU parameter
        c : ITU parameter
        d : ITU parameter

        Examples
        --------

        >>> from pylayers.antprop.slab import *
        >>> M = Mat(name='Phantom',index=17,epr=2+0.15j,mur=1,sigma=4,roughness=0)

        """
        self['name'] = name
        self['index'] = index
        self['epr'] = epr
        self['mur'] = mur
        self['sigma'] = sigma
        self['roughness'] = roughness
        # Parameters a,b,c,d 
        # Table 3 Rec ITU-R P.2040.1
        #
        # epsrp  = a * fGHZ ** b 
        # sigma  = c * fGHZ ** d 
        #
        self['a']=None
        self['b']=None
        self['c']=None
        self['d']=None

    def eval(self, fGHz):
        """ evaluate Mat at given frequencies

        Parameters
        ----------

        fGHz : np.array()
            frequency (GHz)


        Notes
        -----

        w = 2*np.pi*f*1e-9
        eps0 = 8.854e-12

        100 MHz = 0.1 GHz
        10 MHz = 0.01 GHz

        sigma/(w*eps0) = sigma/(2*pi*fGHz*1e9*eps0)
        sigma/(w*eps0) = sigma/(2*pi*fGHz*1e9*8.854e-12)
        sigma/(w*eps0) = sigma/(2*pi*fGHz*1e-3*8.854)
        sigma/(w*eps0) = 18 * sigma/fGHz
        """

        self['fGHz'] = fGHz
        if self['a'] == None:
            epsc = self['epr'] - 1j * 18 * abs(self['sigma']) /  self['fGHz']
        else:
            epsr = self['a'] * self['fGHZ']**self['b']
            sigma  = self['c'] * self['fGHZ']**self['d']
            epsc = epsr - 1j * 18 * sigma /  self['fGHz']

        return(epsc)

    def info(self):
        """ display material properties
        """

        print("---------------------------------")
        for k in self:
            print(k, self[k])
        #print " "
        #print "name : ",self.name
        #print "index :",self.index
        #print "epr :",self.epr
        #print "mur : ",self.mur
        #print "sigma :",self.sigma
        #print "roughness : ",self.roughness

    def R(self, fGHz, theta):
        """ Calculate Reflection coefficient on the air|mat interface

        Parameters
        ----------

        fGHz : frequency GHz
        theta : incidence angle referenced from interface normal

        """

        air = Mat()
        lmat = [air, self]

        Nf = len(fGHz)
        Nt = np.shape(theta)[1]

        fGHz.reshape(Nf, 1)
        theta.reshape(1, Nt)
        II = MatInterface(lmat, 0, fGHz, theta)
        II.RT()
        Ro = II.Ro
        Rp = II.Rp

        return Ro, Rp


class MatDB(PyLayers,dict):
    """ MatDB Class : Material database


    Attributes
    ----------
    di : dict
    associate numeric and alphanumeric keys

    """
    def __init__(self, _fileini='matDB.ini',dm={}):
        """ class constructor

        Parameters
        ----------

        _fileini : string

        """
        self.fileini = _fileini
        self.filemat = self.fileini.replace('.ini','.mat')
        if dm=={}:
            self.load(_fileini)
        else:
            #self.update(dm)
            for matname in dm.keys():
                M=Mat(name=matname)
                ddm = dm[matname]
                for k in ddm:
                    M[k] = dm[matname][k]
                self[matname] = M


    def __repr__(self):
        st = 'List of Materials\n'
        st = st+'-------------------\n'
        for k in self:
            epsr  = "%.2f" % abs(self[k]['epr'])
            sigma = "%.2f" % abs(self[k]['sigma'])
            st = st+k+' ('+str(self[k]['index'])+')    |epsr|='+ epsr +' sigma (S/m)='+sigma+'\n'
        return(st)


    def info(self):
        """ get MatDB info

        TODO : make a __repr__
        """
        for i in self:
            S = self[i]
            S.info()

    def dass(self):
        """ create a dictionnary to associate index | name

        This dictionnary helps for quick access to a given slab

        """
        di = {}
        for name in self.keys():
            # get the integer
            index = self[name]['index']
            # associate the integer to name
            di[index] = name
        # update Class variable key association dictionnary
        self.di = di

    def maxindex(self):
        """ find the max value of the index in DB
        """

        maxi = 0
        for i in self.keys():
            v = self[i]['index']
            if (v > maxi):
                maxi = v
        return(maxi)

    def delete(self, name):
        """ Delete a material in the DB

        Parameters
        ----------
        name : string

        """
        self.__delitem__(name)

    def edit(self, name):
        """ Edit a material in the DB

        Parameters
        ----------

        name : vstring

        See Also
        --------

        pylayers.antprop.slab.dass

        """
        data = multenterbox('Material', 'Enter', ('name', 'index', 'epr', 'mur', 'sigma', 'roughness'),
                            (name, M.index, str(M.epr), str(M.mur), M.sigma, M.roughness))
        self['name'] = data[0]
        self['index'] = eval(data[1])
        self['epr'] = eval(data[2])
        self['mur'] = eval(data[3])
        self['sigma'] = eval(data[4])
        self['roughness'] = eval(data[5])
        # update keys association dictionnary
        self.dass()

    def add(self,**kwargs):
        """ add a material in the DB

        Parameters
        ----------

        name : string
            material name
        cval : float or complex
            epsilon or index
        sigma : float or complex
            conductivity
        mur : float
            relative permeability
        typ : string
            {'epsr'|'ind'|,'reim',|'THz'|'itu'}

        Notes
        -----

        Different ways to enter a material are :

        i) 'epsr' : epsr and sigma

        cval = epsr
        sigma = sigma

        ii) 'ind' : indice @ fGHz

        cval = indice

        iii) 'reim' : real(epsr) and imag(epsr) @fGHz

        iv) 'THz'



        Examples
        --------

        >>> from pylayers.antprop.slab import *
        >>> m = MatDB()
        >>> m.load('matDB.ini')
        >>> m.add(name='ConcreteJcB',cval=3.5+0*1j,alpha_cmm1=1.9,fGHz=120,typ='THz')
        >>> m.add(name='GlassJcB',cval=3.5+0*1j,alpha_cmm1=1.9,fGHz=120,typ='THz')
        >>> out = m.save('Jacob.ini')

        """
        defaults = {'name':'MAT',
                    'cval':1+0*1j,
                    'sigma':0,
                    'alpha_cmm1':1,
                    'mur':1,
                    'fGHz':1,
                    'typ':'epsr'
                   }
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]
        # get the next available index
        maxid = self.maxindex()
        M = Mat()
        M['name'] = kwargs['name']
        M['index'] = maxid + 1
        M['fGHz'] = kwargs['fGHz']
        if kwargs['typ'] == 'epsr':
            M['epr'] = kwargs['cval']
            M['sigma'] = kwargs['sigma']

        if kwargs['typ']== 'reim':
            M['epsr'] = kwargs['cval']
            M['n'] = np.sqrt(kwargs['mur']*M['epsr']) # warning check causality
            M['epr'] = np.real(M['epsr'])
            M['epr2'] = np.imag(M['epsr'])
            M['sigma'] = -M['epr2'] * M['fGHz'] / 18.

        if kwargs['typ']  == 'ind':
            M['n'] = kwargs['cval']
            M['epsr'] = kwargs['cval'] ** 2 / kwargs['mur']
            M['epr'] = np.real(M['epsr'])
            M['epr2'] = np.imag(M['epsr'])
            M['sigma'] = -M['epr2'] * M['fGHz'] / 18.
        #
        # Terahertz Dielectric Properties of Polymers Yun-Sik Jin
        # Terahertz characterization of building materials (R.Piesiewicz) El.Jou Jan 2005 Vol 41 N18
        #
        if kwargs['typ'] == 'THz':
            M['n'] = kwargs['cval']
            M['alpha_cmm1'] = kwargs['alpha_cmm1']
            M['kappa'] = 30 * M['alpha_cmm1'] / (4 * np.pi * M['fGHz'])
            M['epr'] = np.real(M['n'] ** 2 - M['kappa'] ** 2)
            M['epr2'] = np.real(2 * M['kappa'] * M['n'])
            M['sigma'] = M['epr2'] * M['fGHz'] / 18.
            M['Z'] = 1.0 / np.sqrt(M['epr'] + 1j * M['epr2'])

        if kwargs['typ']  == 'itu':
            M['a'] = a
            M['b'] = b
            M['c'] = c
            M['d'] = d

        M['mur'] = kwargs['mur']
        M['roughness'] = 0

        self[kwargs['name']] = M
        # updating dictionnary
        self.dass()

    def addgui(self, name='MAT'):
        """ Add a material in the DB

        Parameters
        ----------

        name : string
        default 'MAT'

        """
        max = self.maxindex()
        data = multenterbox('Material', 'Enter', ('name', 'index', 'epr', 'mur', 'sigma', 'roughness'),
                            (name, max + 1, '(1+0j)', '(1+0j)', 0, 0))
        M = Mat()
        M['name'] = data[0]
        M['index'] = eval(data[1])
        M['epr'] = eval(data[2])
        M['mur'] = eval(data[3])
        M['sigma'] = eval(data[4])
        M['roughness'] = eval(data[5])
        self[name] = M
        self.dass()


    def choose(self):
        """ Choose a mat from matdir
        """
        import tkFileDialog
        FD = tkFileDialog
        filename = FD.askopenfilename(filetypes=[("Mat file ", "*.ini"),
                                                 ("All", "*")],
                                      title="Please choose a Material .ini file",
                                      initialdir=matdir)
        _filename = pyu.getshort(filename)
        self.load(_filename)



    def load(self,_fileini):
        """Load a Material from a .ini file

        Parameters
        ----------

        _fileini : string
            name of the matDB file (usually matDB.ini)

        """
        fileini = pyu.getlong(_fileini, pstruc['DIRMAT'])
        materials = ConfigParser.ConfigParser()
        materials.read(fileini)
        self.di = {}
        for k,matname in enumerate(materials.sections()):
            M = Mat(name=matname)
            M['sigma'] = eval(materials.get(matname,'sigma'))
            M['roughness'] = eval(materials.get(matname,'roughness'))
            M['epr'] = eval(materials.get(matname,'epr'))
            M['index'] = k
            self.di[k]=matname
            M['mur'] = eval(materials.get(matname,'mur'))
            self[matname] = M

    def save(self,_fileini='matDB.ini'):
        """ save MatDB in an ini file

        [dict]
        id1 = name1
        [name1]
        epsr =
        mur =
        roughness =
        index =


        """
        fileini = pyu.getlong(_fileini, pstruc['DIRMAT'])
        fd = open(fileini, "w")
        config = ConfigParser.ConfigParser()
        #
        # config names
        #
        config.add_section("dict")

        for vid in self.di.keys():
            config.set("dict", str(vid), self.di[vid])

        for vid in self.di.keys():
            name = self.di[vid]
            config.add_section(name)
            try:
                config.set(name, "epr", str(self[name]['epr']))
            except:
                config.set(name, "epr", '(9+0j)')
            try:
                config.set(name, "mur", str(self[name]['mur']))
            except:
                config.set(name, "mur", '(1+0j)')
            try:
                config.set(name, "sigma", str(self[name]['sigma']))
            except:
                config.set(name, "sigma", '(0+0j)')
            try:
                config.set(name, "roughness", str(self[name]['roughness']))
            except:
                config.set(name, "roughness", '0')
            config.set(name, "index", str(self[name]['index']))

        config.write(fd)
        fd.close()

class Slab(Interface,dict):
    """ Handle a Slab

    Summary
    -------

    A Slab is a sequence of layers which each has

    - a given width
    - an integer index refering a given material in the material DB


    Attributes
    ----------

    name :
        Slab name
    nbmat :
        Number of layers
    index :
        Slab Index
    lmatname :
        list of material name
    lthick :
        list of thickness of layers
    color :
        color of slab dor display
    linewidth :
        linewidth for structure display
    mat : Associated Material Database

    evaluated : Boolean

    """
    def __init__(self, mat=[], name='NEWSLAB',ds={}):
        """ class constructor

        Parameters
        ----------

        mat : matdb
        name : string
        slab name

        """
        # if not specified choose default material database
        #super(Slab,self).__init__()
        Interface.__init__(self)
        if mat==[]:
            self.mat = MatDB()
        else:
            self.mat = mat
        self['name'] = name
        if ds=={}:
            self['index'] = 0
            # lmatname has to be set before lthick
            self['lmatname'] = ['AIR']
            self['lthick'] = [0.1]
            self['color'] = 'black'
            self['linewidth'] = 1.0
        else:
            self['index'] = ds['index']
            # lmatname has to be set before lthick
            self['lmatname'] = ds['lmatname']
            self['lthick'] = ds['lthick']
            self['color'] = ds['color']
            self['linewidth'] = ds['linewidth']
        self['evaluated'] = False
        self.conv()

    def __setitem__(self,key,value):
        """ dictionnary setter

        lmatname can be changed and lthick is imposed @ 5cm per layer
        lthick   can be changed provided the number of layer is correct

        """
        if key == "lmatname":
            nbmat = len(value)
            for na in value:
                if na not in self.mat:
                    print(self.mat.__repr__())
                    raise ValueError(na+ ' not in material Database')
            dict.__setitem__(self,"lmatname", value)
            #dict.__setitem__(self,"nbmat",nbmat)
            dict.__setitem__(self,"lthick",[0.05]*nbmat)
        
        elif key == "lthick":
            #pdb.set_trace()
            #if len(value)!=len(self['lmatname']):
            #    raise ValueError("wrong number of material layers")
            #else:
            dict.__setitem__(self,"lthick",value)
        else:        
            dict.__setitem__(self,key, value)


    def __add__(self,u):
        """ This function makes the addition between 2 Slabs

        """
        U = Slab(self.mat)
        # lmatname should be modified before lthick
        U['lmatname'] = self['lmatname']+u['lmatname']
        U['lthick']   = self['lthick']+u['lthick']
        U['name'] = self['name']+u['name']
        for i in range(len(U['lmatname'])):
            namem = U['lmatname'][i]
        U.conv()
        return(U)

    def __repr__(self):
        st = self['name']+' : '
        st = st + reduce(lambda x,y: x+' | '+y,self['lmatname'])+ ' | '
        st = st + str(self['lthick'])+'\n'
        st = st + '       ' + str(self['color'])+' '+str(self['linewidth'])+'\n'
        if self['evaluated']:
            nf = len(self.fGHz)
            nt = len(self.theta)
            if nf > 1:
                st = st + "f(GHz) : " + str((self.fGHz[0], self.fGHz[-1], nf))+'\n'
            else:

                st = st + "f(GHz) : " + str(self.fGHz[0])+'\n'

            if nt > 1:
                st= st + "theta (rad) : " + str((self.theta[0], self.theta[-1], nt))+'\n'
            else:
                st = st + "theta (rad) : " + str(self.theta[0])+'\n'
        return(st)

    def info(self):
        """ display Slab Info


        Examples
        --------

        .. plot::
            :include-source:

            >>> import numpy as np
            >>> import matplotlib.pyplot as plt
            >>> from pylayers.antprop.slab import *
            >>> sl = SlabDB('matDB.ini','slabDB.ini')
            >>> lmatname = ['PLATRE-57GHz','AIR','PLATRE-57GHz']
            >>> lthick = [0.018,0.03,0.018]
            >>> sl.add('placo',lmatname,lthick)
            >>> theta = np.arange(0,np.pi/2,0.01)
            >>> fGHz = np.array([57.5])
            >>> sl['placo'].ev(fGHz,theta)
            >>> fig,ax=sl['placo'].plotwrt(var='a',typ=['m'])
            >>> plt.ion()
            >>> plt.show()

        """
        print("------------")
        print("name : ", self['name'])
        print("nbmat : ", len(self['lmatname']))
        chaine = "[ "
        for name in self['lmatname']:
            self.mat[name].info()
            if self['evaluated']:
                epsrc = self.mat[name].eval(self.fGHz[0])
                print("epsrc : ", epsrc)
            chaine = chaine + name + ' '
            chaine = chaine + ']'
            print("index : ", self['index'])
            print("color : ", self['color'])
            print("linewidth :", self['linewidth'])
            if self['evaluated']:
                print("---------------------")
                nf = len(self.fGHz)
                nt = len(self.theta)
                if nf > 1:
                    print("f (GHz) : ", (self.fGHz[0], self.fGHz[-1], nf))
                else:
                    print("f (GHz) : ", self.fGHz[0])

                if nt > 1:
                    print("theta : ", (self.theta[0], self.theta[-1], nt))
                else:
                    print("th (rad) : ", self.theta[0])

    def conv(self):
        """ build lmat 

        """
        self['lmat'] = []

        for matname in self['lmatname']:
            mi = self.mat[matname]
            self['lmat'].append(mi)


    def ev(self, fGHz=np.array([1.0]), theta=np.linspace(0, np.pi / 2, 50),compensate=False,RT='RT'):
        """ evaluation of the Slab

        Parameters
        ----------

        fGHz : frequency GHz ( np.array([1.0]) )
        theta : np.array
            incidence angle (from normal) radians
        compensate : boolean
        RT : string 

        """

        if not isinstance(fGHz, np.ndarray):
            fGHz = np.array([fGHz])
        if not isinstance(theta, np.ndarray):
            theta = np.array([theta])
        theta_in = copy.deepcopy(theta)

        self.theta = theta
        self.fGHz = fGHz
        theta_in = copy.deepcopy(theta)

        nf = len(fGHz)
        nt = len(theta)
        #thetai = theta[0]
        #thetaf = theta[-1]
        ### WARNING thetas can be NOT sorted.
        ### thetai should be min(theta)
        ### thetaf should be max(theta)
        #th1 = np.linspace(thetai,thetaf,nt)

        metalic = False
        name1 = '|'.join(mat for mat in self['lmatname'])
        name2 = '|'.join(str(thick) for thick in self['lthick'])
        name = '(' + name1 + ')' + '(' + name2 + ')'
        #super(Slab,self).__init__(fGHz=fGHz, theta=theta, name=name)
        Interface.__init__(self, fGHz, theta, name=name)
        #self.lmat = lmat
        #self.lthick = lthick
        self.n = len(self['lmat']) + 2

        #nf = len(fGHz)
        #nt = np.shape(self.theta)[1]

        Co = np.array(np.zeros([self.nf, self.nt, 2, 2]), dtype=complex)
        Co[:, :, 0, 0] = 1
        Co[:, :, 1, 1] = 1
        # _Co= np.eye(2,dtype=complex)



        Cp = np.array(np.zeros([self.nf, self.nt, 2, 2]), dtype=complex)
        Cp[:, :, 0, 0] = 1
        Cp[:, :, 1, 1] = 1
        # _Cp = np.eye(2,dtype=complex)
        #
        # loop over the s n-1 materials 
        # lmat[0] est toujours l'air ( modifier)
        #
        for i in range(self.n - 1):
            if i == 0: # first material is AIR
                ml = Mat()
            else:
                ml = self['lmat'][i - 1]
            if i == self.n - 2:
                mr = Mat() # last material is AIR
            else:
                mr = self['lmat'][i]
            if mr['name'] == 'METAL':
                Io = np.array(np.ones([self.nf, self.nt, 2, 2]), dtype=complex)
                Io[:, :, 0, 1] = -1
                Io[:, :, 1, 0] = -1
# _Io=np.eye(2,dtype=complex)+np.eye(2)-1
                Ip = np.array(np.ones([self.nf, self.nt, 2, 2]), dtype=complex)
                Ip[:, :, 0, 1] = 1
                Ip[:, :, 1, 0] = 1
# _Ip=np.eye(2,dtype=complex)+np.eye(2)-1
            else:
                if i == self.n - 2:
                    II = MatInterface([ml, mr], 0, fGHz, theta)
                else:
                    II = MatInterface([ml, mr], self['lthick'][i], fGHz, theta)
            #
            # chains the angle, theta can be complex
            #
                theta = II.theta
            #
            # theta depends on frequency f x th
            #
            # THETA = II.THETA

                Io = II.Io
                Ip = II.Ip

            #
            # dot product Co.Io and Cp.Ip
            #
            # old version (keep it for demonstation)
            # -----------
            #for kf in range(nf):
            # for kt in range(nt):
            # T = np.dot(Co[kf,kt,:,:],Io[kf,kt,:,:])
            # Co[kf,kt,:,:] = T
            # U = np.dot(Cp[kf,kt,:,:],Ip[kf,kt,:,:])
            # Cp[kf,kt,:,:] = U
            #
            # Using Einstein summation instead of a for loop increases speed by an order of magnitude
            #

            # Co = np.einsum('ijkl,ijln->ijkn', Co, Io)
            # Cp = np.einsum('ijkl,ijln->ijkn', Cp, Ip)


            ### array broadcasing version , new increased speed in regard of einsum
            Co = np.sum(Co[...,:,:,None]*Io[...,None,:,:], axis=3)
            Cp = np.sum(Cp[...,:,:,None]*Ip[...,None,:,:], axis=3)


            if mr['name'] == 'METAL':
                metalic = True
                break

        self.Io = Co
        self.Ip = Cp

        self.RT(metalic,RT=RT)
# if compensate:
# fGHz = fGHz.reshape(nf,1,1,1)
# th1 = th1.reshape(1,nt,1,1)
# thickness = sum(self['lthick'])
# d = thickness*np.cos(th1)
# self.T = self.T*np.exp(1j*2*np.pi*fGHz*d/0.3)


        # Modification probably not compliant with coverage !!!!
        # TODO !!!
        if compensate:
            thickness = sum(self['lthick'])
            #pdb.set_trace()
            #d = thickness*np.cos(theta)
            d = thickness*np.cos(theta_in[None,:])
            self.T = self.T*np.exp(1j*2*np.pi*
                                    fGHz[:,np.newaxis,np.newaxis,np.newaxis]
                                    *d[:,:,np.newaxis,np.newaxis]
                                    /0.3)
# if 'T' in RT:
# epr = [m['epr'] for m in self['lmat']]
# epr = sum(epr)
# # theta[0] just for 1 freq
# self.costt = np.sqrt((epr-1+np.cos(theta[0])**2)/epr)
# self.sm = sum(self['lthick'])/self.costt
# self.gamma = np.cos(theta[0])/self.costt
# self.alpha = np.array(([1./epr]),dtype=complex)



        self['evaluated'] = True

    def filter(self,win,theta=0):
        """ filtering waveform

        Warning
        -------

        Not implemented yet

        Parameters
        ----------

        win : Waveform

        Returns
        -------

        wout : Waveform

        """
        # get frequency base of the waveform
        f = win.sf.x
        self.ev(f,theta)
        wout = Wafeform()
        return(wout)

    def excess_grdelay(self,fGHz=np.arange(2.4,4.0,0.1),theta=np.array([0])):
        """ calculate transmission excess delay in ns

        Parameters
        ----------

        fGHz : array
        default arange(2.4,4,0.1)
        theta : default 0

        Returns
        -------

        delayo : excess delay polarization o
        delayp : excess delay polarization p


        Examples
        --------

        #>>> from pylayers.antprop.slab import *
        #>>> from matplotlib.pylab import *
        #>>> import numpy as np
        #>>> sl = SlabDB('matDB.ini','slabDB.ini')
        #>>> s1 = sl['PARTITION']
        #>>> fGHz = np.arange(3.1,10.6,0.1)
        #>>> delayo,delayp = s1.excess_grdelay(fGHz,0)
        #>>> lineo = plt.plot(fGHz[0:-1],delayo)
        #>>> linep = plt.plot(fGHz[0:-1],delayp)
        #>>> plt.show()

        """

        assert len(fGHz)>2 , "fGHz too short needs more than one frequency point"

        df = fGHz[1]-fGHz[0]

        self.ev(fGHz,theta=theta,compensate=True)

        # f x th x p x q
        T = self.T

        To = T[:,:,0,0]
        Tp = T[:,:,1,1]

        ao = np.unwrap(np.angle(To),axis=0)
        ap = np.unwrap(np.angle(Tp),axis=0)

        delayo = -np.mean(np.diff(ao,axis=0)/(2*np.pi*df),axis=0)
        delayp = -np.mean(np.diff(ap,axis=0)/(2*np.pi*df),axis=0)

        return (delayo,delayp)

    def tocolor(self, fGHz=np.array([2.4])):
        """  convert slab properrties into a color

        Parameters
        ----------

        fGHz : np.array

        Examples
        --------

        >>> sl = SlabDB('matDB.ini','slabDB.ini')
        >>> s1 = sl['PARTITION']
        >>> col24 = s1.tocolor(np.array([2.4]))
        >>> fGHz = np.arange(0.5,8,100)
        >>> col8 = s1.tocolor(fGHz)

        """

        self.ev(fGHz, theta=np.array([0.0]),compensate=True)
        color = Interface.tocolor(self, fGHz)
        return(color)

    def loss0(self, fGHz=2.4):
        """ calculate loss for theta=0 at frequency (fGHz)

        Parameters
        ----------

        fGHz : frequency (GHz) np.array()
        default 2.4

        Returns
        -------

        Lo : np.array
        Loss at 0 deg polarization ortho
        Lp : np.array
        Loss at 0 deg polarization para



        Examples
        --------

        >>> from pylayers.antprop.slab import *
        >>> sl = SlabDB('matDB.ini','slabDB.ini')
        >>> s1 = sl['AIR']
        >>> Lo,Lp = s1.loss0(2.4)
        >>> assert (np.allclose(Lo[0],0)),'Problem with AIR slab loss'
        >>> assert (np.allclose(Lo[0],Lp[0])),'something wrong with polarization'

        """

        self.ev(fGHz, theta=np.array([0.0]),compensate=True)
        Lo, Lp = Interface.loss0(self, fGHz)
        return(Lo, Lp)

    def losst(self, fGHz, theta):
        """ Calculate loss w.r.t angle and frequency

        Parameters
        ----------

        fGHz : np.array()
        frequency (GHz)

        theta : np.array
        theta angle (radians)

        Returns
        -------

        Lo : np.array
        Loss orthogonal

        Lp : np.array
        Loss paralell

        """
        # for backward compatibility
        if type(theta)==float:
            theta = np.array([theta])

        self.ev(fGHz, theta)
        Lo, Lp = Interface.losst(self, fGHz)
        return(Lo, Lp)

#    def editgui(self):
#        """ edit a Slab in the DB
#
#        """
#        chaine1 = ""
#        chaine2 = ""
#        for i in range(self.nbmat):
#            index_mat = self.imat[i]
#            name_mat = self.mat.di[index_mat]
#            thick = str(self.thickness[i])
#            chaine1 = chaine1 + name_mat + ' '
#            chaine2 = chaine2 + thick + ' '
#
#        data = multenterbox('Slab', 'Enter',
#                            ('name', 'index', 'nbmat', 'imat',
#                             'thickness (cm)', 'color', 'linewidth'),
#                            (self.name, self.index, self.nbmat, chaine1,
#                             chaine2, self.color, str(self.linewidth)))
#
#        self.index = eval(data[1])
#        self.nbmat = eval(data[2])
#        chaine1 = data[3].split()
#        chaine2 = data[4].split()
#
#        tt = [0, 0, 0, 0, 0, 0, 0, 0]
#        th = [0, 0, 0, 0, 0, 0, 0, 0]
#        if (len(chaine1) != len(chaine2)):
#            print ('erreur edit slab')
#            return(-1)
#        for i in range(len(chaine1)):
#            nom = chaine1[i]
#            thick = chaine2[i]
#            index = self.mat[nom].index
#            tt[i] = index
#            th[i] = eval(thick)
#
#        self.imat = tt
#        self.charindex = str(tt)
#        self.thickness = th
#        self.color = data[5]
#        self.linewidth = eval(data[6])
#        self.dass()

    def show(self, fGHz=2.4, theta=np.arange(0, np.pi / 2., 0.01), dtype=np.float64, dB=False):
        """ show slab Reflection and Transmission coefficient

        Parameters
        ----------

        fGHz : float
        theta : np.array
        dtype :
        display : string
        {'modulus'}
        dB : boolean
        False

        """

        self.ev(fGHz, theta)
        if self['evaluated']:
            fig,ax=self.M.plotwrt(var='a',typ=['l20'])

        return fig,ax


class SlabDB(dict):
    """ Slab data base

    Attributes
    ----------

    DB : slab dictionnary

    """
    def __init__(self, filemat='matDB.ini', fileslab='slabDB.ini',ds={},dm={}):
        """ class constructor

        Parameters
        ----------

        filemat : string
        fileslab : string
        ds : dict 
            slab dict read from layout file. if ds == {}   load from files
        dm : dict 
            mat dict read from layout file. 

        """

        self.fileslab = fileslab

        if ds=={}:
            self.mat = MatDB()
            if (filemat != ''):
                self.mat.load(filemat)
            if (fileslab != ''):
                self.load(fileslab)
        else:
            assert(type(dm)==dict)
            assert(type(ds)==dict)
            # copier le contenu du load ici 
            #self.update(ds)
            self.mat = MatDB(_fileini=filemat,dm=dm)
            for slabname in ds:
                S = Slab(name=slabname,mat=self.mat,ds=ds[slabname])
                for k in ds[slabname]:
                    S[k]=ds[slabname][k]

                S['nmat']=len(ds[slabname]['lmatname'])
                S.conv()
                self[slabname]=S

        self.dass()

    def __repr__(self):
        st = 'List of Slabs\n'
        st = st + '-----------------------------'+'\n'+'\n'
        for i in self.keys():
            st = st + self[i].__repr__()
        #    S.info()
        #st =      "Slab file name     : " + self.fileslab+ '\n'
        #st = st + "Material file name : " +  self.mat.fileini+'\n'
        return(st)

    def showall(self):
        """ show all slabs

        """
        lsl = self.keys()
        k = len(lsl)
        nl = k / 2
        cpt = 1
        for k in lsl:
            plt.figure()
            self[k].show()

        plt.show()



    def dass(self):
        """ update conversion dictionnary

        code <--> name

        """
        di = {}
        for name in self.keys():
            index = self[name]['index']
            di[index] = name
        self.di = di

    def maxindex(self):
        """ find the max value of the index in DB
        """

        maxi = 0
        for i in self.keys():
            v = self[i]['index']
            if (v > maxi):
                maxi = v
        return(maxi)

    def delete(self, name):
        """ delete an element from the database

        Parameters
        ----------

        name : string

        """
        self.__delitem__(name)
        self.dass()

    def edit(self, name):
        """ edit a Slab in the DB

        Parameters
        ----------

        name : string

        """
        slab = self[name]
        slab.edit()

    def show(self, name='WOOD', fGHz=np.array([2.4])):
        """ evaluate and show a given slab

        Parameters
        ----------

        name : string
        fGHz : np.array

        """

        slab = self[name]
        slab.ev(fGHz=fGHz)
        fig,ax = slab.M.plotwrt(var='a')

        return fig,ax


    def add(self, name, lmatname, lthick, color='black'):
        """ add a slab from its properties


        Parameters
        ----------

        name       : string
        lmatname   : list of mat name
        lthick     : list ot float
           list of layer thickness in meters


        Examples
        --------

        Examples from the paper:

        "Reflection ant Transmission Properties of Building Materials in D-Band
        for Modeling Future mm-Wave Communication Systems "
        Martin Jacob and Thomas Kurner and Robert Geise and Radoslaw Piesiewicz
        EUCAP 2010

        .. plot::
            :include-source:


            >>> from pylayers.antprop.slab import *
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> sl = SlabDB('matDB.ini','slabDB.ini')

            >>> sl.mat.add(name='ConcreteJc',cval=3.5,alpha_cmm1=1.9,fGHz=120,typ='THz')
            >>> sl.mat.add(name='GlassJc',cval=2.55,alpha_cmm1=2.4,fGHz=120,typ='THz')
            >>> sl.add('ConcreteJc',['ConcreteJc'],[0.049])
            >>> sl.add('DoubleGlass',['GlassJc','AIR','GlassJc'],[0.0029,0.0102,0.0029])
            >>> theta = np.linspace(20,60,100)*np.pi/180
            >>> sl['ConcreteJc'].ev(120,theta)
            >>> f,a=sl['ConcreteJc'].plotwrt(var='a',typ=['l20'])
            >>> fig = plt.figure()
            >>> sl['DoubleGlass'].ev(120,theta)
            >>> f,a = sl['DoubleGlass'].plotwrt(var='a',typ=['l20'])
            >>> freq = np.linspace(110,135,50)
            >>> fig = plt.figure()
            >>> sl['DoubleGlass'].ev(freq,theta)
            >>> sl['DoubleGlass'].pcolor(dB=True)

        Exemple from paper `"[Kiani2007] Glass Characterization for Designing
        Frequency Selective Surfaces to Improve Transmission through Energy saving
        glass windows  Kiani 2007"
        <http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=4554974&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D4554974>`_
        The surface impedance is :math:`R = 4 \Omega`, the thicknesss is :math:`d = 100 nm`

        `Pilkington Spectrum OnLine applet <http://www.pilkington.com/spectrum2/default.aspx?country_code=FR>`_

        `Design of Energy Saving Windows with high Transmission at 900MHz and 1800 MHz
        <http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=530536&fileOId=624944>`_

        .. math::

            \sigma =  \\frac{1}{Rd} = 2.5 10^{6} S/m

        .. plot::
            :include-source:

            >>> from pylayers.antprop.slab import *
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> sl = SlabDB('matDB.ini','slabDB.ini')
            >>> sl.mat.add(name='CoatingPilkington',cval=1,sigma=2.5e6,typ='epsr')
            >>> sl.mat.add(name='GlassPilkington',cval = 6.9,sigma = 5e-4,typ='epsr')
            >>> sl.add('Optitherm382',['CoatingPilkington','GlassPilkington'],[100e-9,0.00382])
            >>> fGHz  = np.linspace(0.9,2.2,50)
            >>> theta = np.linspace(0,np.pi/2,100)
            >>> sl['Optitherm382'].ev(fGHz,theta)
            >>> sl['Optitherm382'].pcolor(dB=True)


        """

        U = Slab(self.mat, name)
        maxi = self.maxindex()
        U['lmatname'] = lmatname
        U['lthick'] = lthick
        U['index'] = maxi + 1
        U['color'] = color
        U['linewidth'] = 1
        U['evaluated'] = False
        U.conv()
        self[name] = U
        self.dass()

    def addgui(self, name):
        """ add a slab in the DB

        Parameters
        ----------
        name

        """

        U = Slab(self.mat, name)
        U.edit()
        self[U.name] = U
        self.dass()


    def load(self,_fileini='slabDB.ini'):
        """Load a Material from a .ini file

        Parameters
        ----------

        _fileini : string 

        """
        fileini = pyu.getlong(_fileini, pstruc['DIRMAT'])
        config = ConfigParser.ConfigParser()
        config.read(fileini)

        self.di={}
        for k,slabname in enumerate(config.sections()):
            # warning the Slab takes the whole Material Database
            S = Slab(name=slabname,mat=self.mat)
            S['lmatname']=eval(config.get(slabname,'lmatname'))
            S['nbmat']=len(S['lmatname'])
            S['color']=config.get(slabname,'color')
            S['index']= k 
            self.di[k] = slabname
            S['lthick']=eval(config.get(slabname,'lthick'))
            S['linewidth']=eval(config.get(slabname,'linewidth'))

            S.conv()
            self[slabname] = S

    def save(self,_fileini='slabDB.ini'):
        """ save SlabDB in an ini file

        Parameters
        ----------

        _fileini : string 

        """

        fileini = pyu.getlong(_fileini, pstruc['DIRSLAB'])
        fd = open(fileini, "w")
        config = ConfigParser.ConfigParser()
        #
        # config names
        #
        config.add_section("dict")
        for vid in self.di.keys():
            config.set("dict", str(vid), self.di[vid])
        for vid in self.di.keys():
            name = self.di[vid]
            config.add_section(name)
            config.set(name, 'color', str(self[name]['color']))
            config.set(name, 'linewidth', self[name]['linewidth'])
            config.set(name, 'lthick', self[name]['lthick'])
            config.set(name, 'index', self[name]['index'])
            config.set(name, 'lmatname', self[name]['lmatname'])

        config.write(fd)
        fd.close()



# class Wedge(Interface,dict):
#     """ Handle a Wedge

#     Summary
#     -------

#     A Wedge is a cone with, on its 2 faces :

#     - a Material/Slab s0
#     - a Material/slab sn


#     Attributes
#     ----------


#     mat : MatDB
#         Associated Material Database 

#     s0 : string |Mat | Slab
#         Material name |Material |Slab on face 0
#     sn : string |Mat | Slab
#         Material name |Material |Slab on face n
 

#     evaluated : Boolean

#     """
#     def __init__(self, mat=[], s0='WOOD', sn='WOOD',alpha=np.array([])):
#         """ class constructor

#         Parameters
#         ----------

#         mat :
#         name : string
#         Wedge name


#         Example
#         -------
#         >>> from pylayers.antprop.slab import *
#         >>> SDB=SlabDB()
#         >>> s0 = SDB['3D_WINDOW_GLASS']
#         >>> sn = SDB['PARTITION']
#         >>> W=Wedge(s0=s0,sn=sn)
        
#         """
#         # if not specified choose default material database

#         super(Wedge,self).__init__()

#         if mat==[]:
#             self.mat = MatDB()
#         else:
#             self.mat = mat


#         # slab 0 is a string name
#         if isinstance(s0,str):
#             self['name0']= s0
#             try:
#                 self['mat0']=self.mat[s0]
#             except:
#                 raise AttributeError('s0 is not a material from MatDB')
#         else: 
#             if isinstance(s0,Mat):
#                 self['mat0']=s0
#                 self['name0']=s0['name']
#             elif isinstance(s0,Slab):
#                 self['mat0']=s0['lmat'][0]
#                 self['name0']=s0['lmatname'][0]

#         # slab n is a string name
#         if isinstance(sn,str):
#             self['namen']= sn
#             try:
#                 self['matn']=self.mat[sn]
#             except:
#                 raise AttributeError('sn is not a material from MatDB')
#         else: 
#             if isinstance(sn,Mat):
#                 self['matn']=sn
#                 self['namen']=sn['name']
#             elif isinstance(sn,Slab):
#                 self['matn']=sn['lmat'][0]
#                 self['namen']=sn['lmatname'][0]
#         self['color'] = 'black'
#         self['linewidth'] = 1.0
#         self['evaluated'] = False
#         self['alpha']=alpha
#         self['N']=alpha/np.pi


#     # def phi0phin(self,u0,si,so):
#     #     """
#     #     Compute angle phi_0 and phi_nfrom face 0 regarding
#     #     to unit vectors u0 and un along face 0 and n respectively

#     #     Attributes
#     #     ----------

#     #     u0 : ndarray (2|3xNp)
#     #         unit vector along Np faces 0
#     #     si : ndarray (2|3x)
#     #         unit vector along incidence ray
#     #     un : ndarray (2|3xNp)
#     #         unit vector along Np faces n
#     #     """

#     #     geu.


def calsig(cval, fGHz, typ='epsr'):
    """ evaluate sigma from epsr or index at a given frequency

    Parameters
    ----------

    cval : complex value
    {epsr | epsr ^2 }
    fGHz : frequency GHz
    type :
    {'epsr' | 'ind'}


    Returns
    -------
    epr1 : ..math::
    sigma : float
    conductivity (S/m)
    delta :



    """

    if typ == 'epsr':
        epsr = cval
    if typ == 'ind':
        epsr = cval ** 2

    epr1 = np.real(epsr)
    epr2 = np.imag(epsr)
    sigma = -epr2 * fGHz / 18

    n = np.sqrt(epsr)
    n2 = -np.imag(n)
    delta = (0.3 / (2 * np.pi * fGHz * n2))

    return(epr1, sigma, delta)




if (__name__ == "__main__"):
    #plt.ion()
    doctest.testmod()
    sl = SlabDB('matDB.ini','slabDB.ini')
    s1 = sl['PILLAR']
    fGHz=np.arange(0.6,5.0,0.1)


