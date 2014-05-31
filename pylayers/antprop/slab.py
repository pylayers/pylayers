#!/usr/bin/python
# -*- coding: latin1 -*-
"""
This module exploits heavily numpy broadcasting mechanism :


Axis Convention

+ nf : axis = 0 frequency axis
+ nt : axis = 1 angular axis
+ p : axis = 2 parallel axis
+ o : axis = 2 orhogonal axis


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
    MatDB.loadmat
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
import os
import string
import cPickle
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
import pdb

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

            R = \matrix(R_o & 0\\0 & R_p)
            T = \matrix(T_o & 0\\0 & T_p)

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

        Returns
        -------

        Lo : np.array
        Loss orthogonal polarization (dB)
        Lp : np.array
        Loss parallel polarization (dB)

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
            >>> lmat = [mat['AIR'],mat['WOOD']]
            >>> II = MatInterface(lmat,0,fGHz,theta)
            >>> II.RT()
            >>> fig,ax = II.plotwrt(var='a',kv=10,typ=['m'])
            >>> air = mat['AIR']
            >>> brick = mat['BRICK']
            >>> II = MatInterface([air,brick],0,fGHz,theta)
            >>> II.RT()
            >>> fig,ax = II.plotwrt(var='f',color='k',typ=['m'])
            >>> plt.show()


        """
        defaults = {'typ':['l20'],
                'polar':'po', # po | p | o
                'coeff':'RT', # RT | R | T
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
        Interface.__init__(self, fGHz, theta, name=name)
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

        epd = np.exp(1j * deltai)
        emd = np.exp(-1j * deltai)

        self.Ip[:, :, 0, 0] = epd / Tp
        self.Ip[:, :, 0, 1] = Rp * emd / Tp
        self.Ip[:, :, 1, 0] = Rp * epd / Tp
        self.Ip[:, :, 1, 1] = emd / Tp

        self.Io[:, :, 0, 0] = epd / To
        self.Io[:, :, 0, 1] = Ro * emd / To
        self.Io[:, :, 1, 0] = Ro * epd / To
        self.Io[:, :, 1, 1] = emd / To

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
    def __init__(self, name="AIR", index=1, epr=1 + 0.0j, mur=1 + 0.0j, sigma=0.0, roughness=0.):
        """ class constructor

        Parameters
        ----------

        name : string
        index : int
        epr : complex
        mur : complex
        sigma : float
        roughness : float

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

        sigma/(w*eps0) = sigma/(2*pi*f*1e9*eps0)
        sigma/(w*eps0) = sigma/(2*pi*f*1e9*8.854e-12)
        sigma/(w*eps0) = sigma/(2*pi*f*1e-3*8.854)
        """

        self['fGHz'] = fGHz
        epsc = self['epr'] - 1j * 18 * abs(self['sigma']) /  self['fGHz']

        return(epsc)

    def info(self):
        """ display material properties
        """

        print "---------------------------------"
        for k in self:
            print k, self[k]
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
        #Ro = II.Ro
        Ro = II.Ro # modif eric plouhinec
        Rp = II.Rp

        return Ro, Rp


class MatDB(PyLayers,dict):
    """ MatDB Class : Material database


    Attributes
    ----------
    di : dict
    associate numeric and alphanumeric keys

    """
    def __init__(self, _fileini='matDB.ini'):
        """ class constructor

        Parameters
        ----------

        _fileini : string

        """
        self.fileini = _fileini
        self.filemat = self.fileini.replace('.ini','.mat')



    def info(self):
        """ get MatDB info

        TODO : make a __repr__
        """
        for i in self:
            S = self[i]
            S.info()

    def dass(self):
        """ create a dictionnary to associate index | name

        Example of created dictionnary

        {-1: 'METALIC',
        0: 'ABSORBENT',
        1: 'AIR',
        2: 'WALL',
        3: 'PARTITION',
        4: 'WINDOW',
        5: 'DOOR',
        6: 'CEIL',
        7: 'FLOOR',
        8: 'WINDOW_GLASS',
        9: 'WOOD',
        10: '3D_WINDOW_GLASS',
        11: 'WALLS',
        12: 'PILLAR',
        13: 'METAL',
        14: 'CONCRETE_15CM3D',
        15: 'CONCRETE_20CM3D',
        16: 'CONCRETE_6CM3D',
        17: 'CONCRETE_7CM3D',
        18: 'PLASTERBOARD_10CM',
        19: 'PLASTERBOARD_14CM',
        20: 'PLASTERBOARD_7CM'}

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

    def add(self, name='MAT', cval=1 + 0 * 1j, sigma=0, alpha_cmm1=1, mur=1, fGHz=1, typ='epsr'):
        """ Add a material in the DB

        Parameters
        ----------

        name : string
        material name
        cval : float or complex
        epsilon or index
        sigma : float or complex
        permeability
        mur : float
        typ : string
        {'epsr'|'ind'|,'reim',|'THz'}

        Notes
        -----

        Different ways to enter a material are :

        i) epsr : epsr and sigma
        cval = epsr
        sigma = sigma
        ii) ind : indice @ fGHz
        cval = indice
        iii) reim : real(epsr) and imag(epsr) @fGHz
        iv) THZ



        Examples
        --------

        >>> from pylayers.antprop.slab import *
        >>> m = MatDB()
        >>> m.load('matDB.ini')
        >>> m.add('ConcreteJcB',cval=3.5+0*1j,alpha_cmm1=1.9,fGHz=120,typ='THz')
        >>> m.add('GlassJcB',cval=3.5+0*1j,alpha_cmm1=1.9,fGHz=120,typ='THz')
        >>> out = m.save('Jacob.ini')

        """

        # get the next available index
        maxid = self.maxindex()
        M = Mat()
        M['name'] = name
        M['index'] = maxid + 1
        M['fGHz'] = fGHz
        if typ == 'epsr':
            M['epr'] = cval
            M['sigma'] = sigma

        if typ == 'reim':
            M['epsr'] = cval
            M['n'] = np.sqrt(mur*M['epsr']) # warning check causality
            M['epr'] = np.real(M['epsr'])
            M['epr2'] = np.imag(M['epsr'])
            M['sigma'] = -M['epr2'] * M['fGHz'] / 18

        if typ == 'ind':
            M['n'] = cval
            M['epsr'] = cval ** 2 / mur
            M['epr'] = np.real(M['epsr'])
            M['epr2'] = np.imag(M['epsr'])
            M['sigma'] = -M['epr2'] * M['fGHz'] / 18
        #
        # Terahertz Dielectric Properties of Polymers Yun-Sik Jin
        # Terahertz characterization of building materials (R.Piesiewicz) El.Jou Jan 2005 Vol 41 N°18
        #
        if typ == 'THz':
            M['n'] = cval
            M['alpha_cmm1'] = alpha_cmm1
            M['kappa'] = 30 * M['alpha_cmm1'] / (4 * np.pi * M['fGHz'])
            M['epr'] = np.real(M['n'] ** 2 - M['kappa'] ** 2)
            M['epr2'] = np.real(2 * M['kappa'] * M['n'])
            M['sigma'] = M['epr2'] * M['fGHz'] / 18
            M['Z'] = 1.0 / np.sqrt(M['epr'] + 1j * M['epr2'])

        M['mur'] = mur
        M['roughness'] = 0

        self[name] = M
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
        config = ConfigParser.ConfigParser()
        config.read(fileini)

        di = dict(config.items("dict") )

        self.di={}
        for d in di:
            self.di[eval(d)]=di[d]

        for matname in self.di.values():
            M=Mat(name=matname)
            M['sigma'] = eval(config.get(matname,'sigma'))
            M['roughness'] = eval(config.get(matname,'roughness'))
            M['epr'] = eval(config.get(matname,'epr'))
            M['index'] = eval(config.get(matname,'index'))
            M['mur'] = eval(config.get(matname,'mur'))
            self[matname] = M

        # PULSRAY compatibility : save in the old .mat format 
        self.savemat(self.filemat)

    def loadmat(self, _filemat):
        """ Load a Material from a .mat file

        Parameters
        ----------

        _filemat : string
        a short file name

        Notes 
        -----

            Deprecated this the format for PyRay

        """
        filemat = pyu.getlong(_filemat, pstruc['DIRMAT'])
        try:
            fo = open(filemat, "rb")
            data = fo.read()
            #
            # decodage des donnees lues
            #
            data_listname = data[0:1200]
            self.tname = data_listname.replace(
                "\x00", "").replace("\"", "").split()
            data_N = data[1200:1204]
            self.N = stru.unpack('i', data_N)[0]

            for i in range(self.N):
                # Creation d'un objet Mat
                M = Mat()
                delta = i * 82
                data_name = data[1204 + delta:1234 + delta]
                name = data_name.replace("\x00", "")
                M['name'] = name

                data_index = data[1234 + delta:1238 + delta]
                index = stru.unpack('i', data_index)[0]
                M['index'] = index

                data_err = data[1238 + delta:1246 + delta]
                err = stru.unpack('d', data_err)[0]
                data_eri = data[1246 + delta:1254 + delta]
                eri = stru.unpack('d', data_eri)[0]
                epr = err + 1j * eri
                M['epr'] = epr
                data_mur = data[1254 + delta:1262 + delta]
                mur = stru.unpack('d', data_mur)[0]
                data_mui = data[1262 + delta:1270 + delta]
                mui = stru.unpack('d', data_mui)[0]
                mur = mur + 1j * mui
                M['mur'] = mur
                data_sigma = data[1270 + delta:1278 + delta]
                sigma = stru.unpack('d', data_sigma)[0]
                M['sigma'] = sigma
                data_roughness = data[1278 + delta:1286 + delta]
                roughness = stru.unpack('d', data_roughness)[0]
                M['roughness'] = roughness
                self[name] = M
            fo.close()
        except:
            print "file : ", filename, "is unreachable"
        self.dass()

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

    def savemat(self, _filemat):
        """ save a .mat file (PulsRay format)

        Parameters
        ----------

        _filemat : string
        a short file name

        """
        filemat = pyu.getlong(_filemat, pstruc['DIRMAT'])
        fo = open(filemat, 'wb')
        N = len(self.di)
        data_listname = ''
        for k in range(N):
            data_listname = data_listname + "\"" + self.di[k - 1] + "\" "

        L = len(data_listname)
        if L < 1200:
            for i in range(1200 - L):
                data_listname = data_listname + "\x00"
        else:
            print " out of range in save Mat - too many materials"
        data_N = stru.pack('i', N)
        data = data_listname + data_N

        for i in range(len(self.di)):
            key = self.di[i - 1]
            M = self[key]
            data_name = key
            L = len(data_name)
            if L < 30:
                for j in range(30 - L):
                    data_name = data_name + "\x00"
            else:
                    print " Mat : name too long maximum 30 characters !"

            data_index = stru.pack('i', M['index'])
            data_err = stru.pack('d', np.real(M['epr']))
            data_eri = stru.pack('d', np.imag(M['epr']))
            data_epr = data_err + data_eri
            data_mur = stru.pack('d', np.real(M['mur']))
            data_mui = stru.pack('d', np.imag(M['mur']))
            data_mu = data_mur + data_mui
            data_sigma = stru.pack('d', M['sigma'])
            data_roughness = stru.pack('d', M['roughness'])
            data_mat = data_name + data_index + data_epr + \
                data_mu + data_sigma + data_roughness
            data = data + data_mat

        fo.write(data)
        fo.close()


class Slab(dict, Interface):
    """ Handle slab

    Summary
    -------
    A slab is a sequence of layers which has
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
    imat :
    index of layers material
    thickness :
    thickness of layers
    color :
    color of slab dor display
    linewidth :
    linewidth for structure display
    mat :
    Associated Material Database
    evaluated : Boolean

    """
    def __init__(self, mat, name='NEWSLAB'):
        """ class constructor

        Parameters
        ----------

        mat :
        name : string
        slab name

        """
        self['name'] = name
        self['index'] = 0
        self['nbmat'] = 1
        self['imat'] = (0, 0, 0, 0, 0, 0, 0, 0)
        self['thickness'] = (10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self['lthick'] = [0.1]
        self['color'] = 'black'
        self['linewidth'] = 1.0
        self.mat = mat
        self['evaluated'] = False

#    def __repr__(self):
#        if self['evaluated']:
#            st = 'fGHz : '+str(self.fGHz[0,:]) +':'+str(self.fGHz[-1,:]) +':'+str(len(self.fGHz[:,0]))+"\n"
#            st = st+ 'theta : '+str(self.theta[:,0])+':'+str(self.theta[:,-1])+':'+str(len(self.theta[0,:]))+"\n"
#            st = st + '| '
#        else:
#            st = '| '
#        for k in range(len(self['lmatname'])):
#            st = st + self['lmatname'][k]+' | '
#        st = st+'\n|'
#
#        for k in range(len(self['lmatname'])):
#            ntick = int(np.ceil(self['lthick'][k]/0.01))
#            for l in range(ntick):
#                st = st+'-'
#            st = st +'|'
#        st = st+'\n'
#        return(st)

    def info(self):
        """ Display Slab Info


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
            >>> plt.show()

        """
        print "------------"
        print "name : ", self
        print "nbmat : ", len(self['lmatname'])
        chaine = "[ "
        for name in self['lmatname']:
            self.mat[name].info()
            if self['evaluated']:
                epsrc = self.mat[name].epsc(self.fGHz[0])
                print "epsrc : ", epsrc
            chaine = chaine + name + ' '
            chaine = chaine + ']'
            print chaine
            print "index : ", self['index']
            print "imat : ", self['imat']
            print "thickness (cm) : ", self['thickness']
            print "color : ", self['color']
            print "linewidth :", self['linewidth']
            if self['evaluated']:
                print "---------------------"
                nf = len(self.fGHz)
                nt = len(self.theta)
                if nf > 1:
                    print "f (GHz) : ", (self.fGHz[0], self.fGHz[-1], nf)
                else:
                    print "f (GHz) : ", self.fGHz[0]

                if nt > 1:
                    print "theta : ", (self.theta[0], self.theta[-1], nt)
                else:
                    print "th (rad) : ", self.theta[0]

    def conv(self):
        """ build lmat and thick

        Warnings
        --------

        In .slab file thickness variable is expressed in cm

        for lthick distance are expressed in meters

        """
        #m1 = self.mat['AIR']
        self['lmat'] = []
        self['lthick'] = []
        #self.lmat.append(m1)

        for i in range(self['nbmat']):
            index_mat = self['imat'][i]
            name_mat = self.mat.di[index_mat]
            mi = self.mat[name_mat]
            self['lmat'].append(mi)
            self['lthick'].append(self['thickness'][i] * 0.01) # cm ->m

        #self.lmat.append(m1)
        #
        # ..todo::
        # fix this problem in MLAYER need thickness 0 for last medium
        #
        #self.thick.append(0.0)

    def ev(self, fGHz=np.array([1.0]), theta=np.linspace(0, np.pi / 2, 50),compensate=False,RT='RT'):
        """ evaluation of the slab

        Parameters
        ----------

        fGHz : frequency GHz ( np.array([1.0]) )
        theta : np.array
            incidence angle (from normal) radians

        """

        if not isinstance(fGHz, np.ndarray):
            fGHz = np.array([fGHz])
        if not isinstance(theta, np.ndarray):
            theta = np.array([theta])
        self.theta = theta
        self.fGHz = fGHz

        nf = len(fGHz)
        nt = len(theta)
        #thetai = theta[0]
        #thetaf = theta[-1]
        ### WARNING thetas can be NOT sorted.
        ### thetai should be min(theta)
        ### thetaf should be max(theta)
        #th1 = np.linspace(thetai,thetaf,nt)

        metalic = False
        name1 = '|'.join(mat['name'] for mat in self['lmat'])
        name2 = '|'.join(str(thick) for thick in self['lthick'])
        name = '(' + name1 + ')' + '(' + name2 + ')'
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
        # Boucle sur les n-1 matériaux
        # lmat[0] est toujours l'air (à modifier)
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
                Ip[:, :, 0, 1] = -1
                Ip[:, :, 1, 0] = -1
# _Ip=np.eye(2,dtype=complex)+np.eye(2)-1
            else:
                if i == self.n - 2:
                    II = MatInterface([ml, mr], 0, fGHz, theta)
                else:
                    II = MatInterface([ml, mr], self['lthick'][i], fGHz, theta)
            #
            # chains the angle
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
            Co=np.sum(Co[...,:,:,np.newaxis]*Io[...,np.newaxis,:,:], axis=3)
            Cp=np.sum(Cp[...,:,:,np.newaxis]*Ip[...,np.newaxis,:,:], axis=3)


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
            d = thickness*np.cos(theta)
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

        Parameters
        ----------

        win : waveform

        Returns
        -------

        wout :

        """
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
        >>> s1 = sl['PARTITION']
        >>> Lo,Lp = s1.loss0(2.4)
        >>> assert ((Lo[0]>5.54)&(Lo[0]<5.56)),'def Partition has changed'
        >>> assert (Lo[0]==Lp[0]),'something wrong with polarization'

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

    def editgui(self):
        """ edit a Slab in the DB

        """
        chaine1 = ""
        chaine2 = ""
        for i in range(self.nbmat):
            index_mat = self.imat[i]
            name_mat = self.mat.di[index_mat]
            thick = str(self.thickness[i])
            chaine1 = chaine1 + name_mat + ' '
            chaine2 = chaine2 + thick + ' '

        data = multenterbox('Slab', 'Enter',
                            ('name', 'index', 'nbmat', 'imat',
                             'thickness (cm)', 'color', 'linewidth'),
                            (self.name, self.index, self.nbmat, chaine1,
                             chaine2, self.color, str(self.linewidth)))

        self.index = eval(data[1])
        self.nbmat = eval(data[2])
        chaine1 = data[3].split()
        chaine2 = data[4].split()

        tt = [0, 0, 0, 0, 0, 0, 0, 0]
        th = [0, 0, 0, 0, 0, 0, 0, 0]
        if (len(chaine1) != len(chaine2)):
            print ('erreur edit slab')
            return(-1)
        for i in range(len(chaine1)):
            nom = chaine1[i]
            thick = chaine2[i]
            index = self.mat[nom].index
            tt[i] = index
            th[i] = eval(thick)

        self.imat = tt
        self.charindex = str(tt)
        self.thickness = th
        self.color = data[5]
        self.linewidth = eval(data[6])
        self.dass()

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
    def __init__(self, filemat='matDB.ini', fileslab='slabDB.ini'):
        """ class constructor

        Parameters
        ----------

        filemat : string
        fileslab : string

        """

        self.fileslab = fileslab
        self.fileslab = self.fileslab.replace('.ini','.slab') # WARNING !!! deprecated in new verion
        self.mat = MatDB()
        if (filemat != ''):
            self.mat.load(filemat)
        if (fileslab != ''):
            self.load(fileslab)
            self.dass()

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

    def info(self):
        """ information
        """
        print "fileslab : ", self.fileslab
        print "filemat : ", self.mat.filemat
        for i in self.keys():
            S = self[i]
            S.info()

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
        """ add a slab in dB


        Parameters
        ----------

        name       : string
        lmatname      : list of mat name
        lthick     : list ot float
            lthick  is in meters

        Warnings
        --------

        thickness is in cm in .slab

        Examples
        --------

        Example from the paper:
        "Reflection ant Transmission Properties of Building Materials in D-Band
        for Modeling Future mm-Wave Communication Systems "
        Martin Jacob and Thomas Kurner and Robert Geise and Radoslaw Piesiewicz
        EUCAP 2010

        .. plot::
            :include-source:


            from pylayers.antprop.slab import *
            import numpy as np
            import matplotlib.pylab as plt
            sl = SlabDB('matDB.ini','slabDB.ini')

            sl.mat.add('ConcreteJc',cval=3.5,alpha_cmm1=1.9,fGHz=120,typ='THz')
            sl.mat.add('GlassJc',cval=2.55,alpha_cmm1=2.4,fGHz=120,typ='THz')
            sl.add('ConcreteJc',['ConcreteJc'],[0.049])
            sl.add('DoubleGlass',['GlassJc','AIR','GlassJc'],[
                0.0029,0.0102,0.0029])
            theta = np.linspace(20,60,100)*np.pi/180
            sl['ConcreteJc'].ev(120,theta)
            sl['ConcreteJc'].plotwrt(var='a',typ=['l20'])
            fig = plt.figure()
            sl['DoubleGlass'].ev(120,theta)
            sl['DoubleGlass'].plotwrt(var='a',typ=['l20'])
            freq = np.linspace(110,135,50)
            fig = plt.figure()
            sl['DoubleGlass'].ev(freq,theta)
            sl['DoubleGlass'].pcolor(dB=True)

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

            from pylayers.antprop.slab import *
            import numpy as np
            import matplotlib.pylab as plt
            sl = SlabDB('matDB.ini','slabDB.ini')
            sl.mat.add('CoatingPilkington',cval=1,sigma=2.5e6,typ='epsr')
            sl.mat.add('GlassPilkington',cval = 6.9,sigma = 5e-4,typ='epsr')
            sl.add('Optitherm382',['CoatingPilkington',
                'GlassPilkington'],[100e-9,0.00382])
            fGHz  = np.linspace(0.9,2.2,50)
            theta = np.linspace(0,np.pi/2,100)
            sl['Optitherm382'].ev(fGHz,theta)
            sl['Optitherm382'].pcolor(dB=True)


        """

        U = Slab(self.mat, name)
        maxi = self.maxindex()
        U['lmatname'] = lmatname
        U['lthick'] = lthick
        U['index'] = maxi + 1
        U['nbmat'] = len(lmatname)
        imat = np.zeros(8).astype(int)
        thickness = np.zeros(8)
        for i in range(len(lmatname)):
            namem = lmatname[i]
            imat[i] = U.mat[namem]['index']
            thickness[i] = lthick[i] * 100 # m ->cm
        U['imat'] = tuple(imat)
        U['thickness'] = tuple(thickness) # cm
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

# Parameters
# ----------
# _filename : string

# """
# filename = pyu.getlong(_filename, pstruc['DIRSLAB'])
# fo = open(filename, 'r')
# DB = cPickle.load(fo)
# self = DB
# self.conv()
# fo.close()

# def choose(self):
# """ Choose a mat file from matdir and slab from slabdir
# """
# import tkFileDialog
# FD = tkFileDialog

# self.mat.choose()
# fileslab = FD.askopenfilename(filetypes=[("Slab file ", "*.slab"),
# ("All", "*")],
# title="Please choose a .slab file",
# initialdir=slabdir)
# _fileslab = pyu.getshort(fileslab)
# self.load(_fileslab)

    def load(self,_fileini='slabDB.ini'):
        """Load a Material from a .ini file

        """
        fileini = pyu.getlong(_fileini, pstruc['DIRMAT'])
        config = ConfigParser.ConfigParser()
        config.read(fileini)

        di = dict(config.items("dict") )
        self.di={}
        for d in di:
            self.di[eval(d)]=di[d]
        for slabname in self.di.values():
            S=Slab(name=slabname,mat=self.mat)
            S['lmatname']=eval(config.get(slabname,'lmatname'))
            S['nbmat']=len(S['lmatname'])
            S['color']=config.get(slabname,'color')
            S['index']=eval(config.get(slabname,'index'))
            S['lthick']=eval(config.get(slabname,'lthick'))
            S['linewidth']=eval(config.get(slabname,'linewidth'))
            imat=[0,0,0,0,0,0,0,0]
            for i,m in enumerate(S['lmatname']):
                imat[i]=S.mat[m]['index']
            S['imat']=tuple(imat)

            thickness=[0,0,0,0,0,0,0,0]
            for i,t in enumerate(S['lthick']):
                thickness[i]=t*100.
            S['thickness']=tuple(thickness)
            S.conv()
            self[slabname] = S
        self.savesl(self.fileslab)



    def loadsl(self, _filename):
        """ load a .slab file (PulsRay format)

        Parameters
        ---------

        _filename

        """
        filename = pyu.getlong(_filename, pstruc['DIRSLAB'])
        try:
            fo = open(filename, "rb")
        except:
            raise NameError("Exception in load : file is unreachable")

        data = fo.read()
        #
        # decodage des donnees lues
        #
        data_listname = data[0:1200]
        self.tname = data_listname.replace(
            "\x00", "").replace("\"", "").split()
        data_N = data[1200:1204]
        self.N = stru.unpack('i', data_N)[0]

        laycol = {}
        laycol['METALIC'] = 'black'
        laycol['ABSORBENT'] = 'honeydew2'
        laycol['AIR'] = 'white'
        laycol['WALL'] = 'grey20'
        laycol['WOOD'] = 'maroon'
        laycol['DOOR'] = 'maroon'
        laycol['CEIL'] = 'grey20'
        laycol['FLOOR'] = 'grey20'
        laycol['WINDOW'] = 'blue'
        laycol['WINDOW_GLASS'] = 'blue2'
        laycol['3D_WINDOW_GLASS'] = 'blue3'
        laycol['WALLS'] = 'grey20'
        laycol['PARTITION'] = 'orchid1'
        laycol['PILLAR'] = 'brown4'
        laycol['METAL'] = 'black'
        laycol['CONCRETE_15CM3D'] = 'DimGrey'
        laycol['CONCRETE_20CM3D'] = 'SlateGray'
        laycol['CONCRETE_7CM3D'] = 'LightGray'
        laycol['CONCRETE_6CM3D'] = 'LightGray'
        laycol['PLASTERBOARD_10CM'] = 'tomato'
        laycol['PLASTERBOARD_14CM'] = 'red'
        laycol['PLASTERBOARD_7CM'] = 'pink'
        laycol['TATA'] = 'yellow'
        laycol['TOTO'] = 'yellow'

        for i in range(self.N):
            # Creation d'un objet slab
            S = Slab(self.mat)

            delta = i * 168

            data_nb_mat_slab = data[1204 + delta:1208 + delta]
            #
            # !! Bug in PulsRay nb_mat is the same than nbmat
            #
            #S['nb_mat']= stru.unpack('i',data_nb_mat_slab)[0]

            data_slabname = data[1208 + delta:1238 + delta]
            name = data_slabname.replace("\x00", "")

            data_charindex = data[1238 + delta:1268 + delta]
            charindex = data_charindex.replace("\x00", "")
    # S['charindex']=charindex

            data_slab_index = data[1268 + delta:1272 + delta]
            index = stru.unpack('i', data_slab_index)[0]
            S['index'] = index
            data_nb_mat = data[1272 + delta:1276 + delta]
            nbmat = stru.unpack('i', data_nb_mat)[0]
            S['nbmat'] = nbmat
            data_index_mat = data[1276 + delta:1308 + delta]
            index_mat = stru.unpack('8i', data_index_mat)
            S['imat'] = index_mat
            data_thickness = data[1308 + delta:1372 + delta]
            thickness = stru.unpack('8d', data_thickness)
            epaisseur = thickness[0]
            # en cm
            S['thickness'] = thickness
            # insertion du slab dans la liste
            if name in laycol.keys():
                S['color'] = laycol[name]
            else:
                S['color'] = 'black'
            S['linewidth'] = int(np.ceil(epaisseur / 3.))
            S['name'] = name
            S.conv()
            self[name] = S

        fo.close()
        self.dass()

    def save(self,_fileini='slabDB.ini'):
        """ save SlabDB in an ini file
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
            lmat=[]
            for i in self[name]['imat']:
                if i !=0:
                    lmat.append(self.mat.di[i])
                else:
                    if lmat ==[]:
                        lmat=['ABSORBENT']
                        break
            config.set(name, 'lmatname', lmat)

        config.write(fd)
        fd.close()

    def savesl(self, _filename):
        """ Save a Slab database in a .slab file (PulsRay format)

        Parameters
        ----------

        _filename : string
        shortname of slabfile

        """

        filename = pyu.getlong(_filename, pstruc['DIRSLAB'])
        fo = open(filename, 'wb')
        N = len(self.di)
        tname = self.keys()
        data_listname = ''
        for k in range(N):
            data_listname = data_listname + "\"" + self.di[k - 1] + "\" "
# data_listname = string.join(tname,"\" \"")
# data_listname = "\""+data_listname+"\""
        L = len(data_listname)
        # On complete avec \x00 jusqu'a 2500 valeurs
        if L < 1200:
            for i in range(1200 - L):
                data_listname = data_listname + "\x00"
        else:
            print " out of range in save Slab - too many slabs"

        data_N = stru.pack('i', N)
        data = data_listname + data_N

        for i in range(N):
            ## Creation d'un dictionnaire slab
            key = self.di[i - 1]
            #print key
            S = self[key]
            data_nb_mat_slab = stru.pack('i', S['nbmat'])
            L = len(key)
            data_slabname = key
            if L < 30:
                for i in range(30 - L):
                    data_slabname = data_slabname + "\x00"
            else:
                print " Slab : slabname too long maximum 30 characters !"

            data_charindex = str(S['imat'])
            L = len(data_charindex)
            if L < 30:
                for i in range(30 - L):
                    data_charindex = data_charindex + "\x00"
            else:
                print " Slab : charindex too long maximum 30 characters !"

            data_slab_index = stru.pack('i', S['index'])
            data_nb_mat = stru.pack('i', S['nbmat'])

            data_index_mat = ""
            for j in range(len(S['imat'])):
                data_index_mat = data_index_mat + stru.pack('i', S['imat'][j])

            data_thickness = ""
            for j in range(len(S['thickness'])):
                data_thickness = data_thickness + \
                    stru.pack('d', S['thickness'][j])
            data_slab = data_nb_mat + data_slabname + data_charindex + data_slab_index + data_nb_mat + data_index_mat + data_thickness
            data = data + data_slab

        fo.write(data)
        fo.close()

# def savesl(self, _filename):
# """ the short filename _filename needs to have the extension .sl

# """
# ia = _filename.find('.sl')
# if (ia == -1):
# print("error : _filename needs a .sl extension")
# exit()
# else:
# # save in .sl format
# filename = pyu.getlong(_filename, pstruc['DIRSLAB'])
# fo = open(filename, 'w')
# cPickle.dump(self, fo)
# fo.close()
# # save in .slab format
# f2 = _filename.replace('.sl', pstruc['DIRSLAB'])
# self.save(f2)


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


