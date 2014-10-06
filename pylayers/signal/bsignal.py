#!/usr/bin/python
#-*- coding:Utf-8 -*-
r"""
.. currentmodule:: pylayers.signal.bsignal

.. autosummary::
    :toctree: generated

Bsignal Class
=============

.. autosummary::
    :toctree: generated/

    Bsignal.__init__
    Bsignal.__repr__
    Bsignal.extract
    Bsignal.save
    Bsignal.load
    Bsignal.setx
    Bsignal.sety
    Bsignal.stem
    Bsignal.step
    Bsignal.cformat
    Bsignal.imshow
    Bsignal.plot
    Bsignal.flatteny
    Bsignal.gating
    Bsignal.len

Usignal Class
=============

.. autosummary::
    :toctree: generated/

    Usignal.__init__
    Usignal.__repr__
    Usignal.__add__
    Usignal.__sub__
    Usignal.__mul__
    Usignal.setx
    Usignal.dx
    Usignal.width
    Usignal.expand
    Usignal.max
    Usignal.min
    Usignal.truncate
    Usignal.align
    Usignal.abs
    Usignal.eprfl
    Usignal.energy
    Usignal.fftshift
    Usignal.zright
    Usignal.zleft
    Usignal.zlr

TBsignal Class
==============

.. autosummary::
    :toctree: generated/

    TBsignal.__init__
    TBsignal.__repr__
    TBsignal.plot
    TBsignal.translate
    TBsignal.b2u

TUsignal Class
==============

.. autosummary::
    :toctree: generated/

    TUsignal.__init__
    TUsignal.__repr__
    TUsignal.diff
    TUsignal.info

Fourier Functions
-----------------

.. autosummary::
    :toctree: generated/

    TUsignal.ft
    TUsignal.fft
    TUsignal.fftsh
    TUsignal.ftshift
    TUsignal.psd
    TUsignal.Yadd_zeros2l
    TUsignal.Yadd_zeros2r
    TUsignal.esd

CIR Functions
-------------

.. autosummary::
    :toctree: generated/

    TUsignal.aggcir
    TUsignal.ecdf
    TUsignal.taumax
    TUsignal.tau_moy
    TUsignal.delays
    TUsignal.tau_rms

Visualization Functions
-----------------------

.. autosummary::
    :toctree: generated/

    TUsignal.show
    TUsignal.shift
    TUsignal.filter
    TUsignal.correlate
    TUsignal.corrgauss
    TUsignal.resample
    TUsignal.convolve

Energy Functions
----------------

.. autosummary::
    :toctree: generated/

    TUsignal.Epercent
    TUsignal.Etau0
    TUsignal.Ewin
    TUsignal.Etot
    TUsignal.Efirst
    TUsignal.Efirst_loc
    TUsignal.Efirst_corr
    TUsignal.Efirst_toath
    TUsignal.Emax

TOA Estimation Functions
------------------------

.. autosummary::
    :toctree: generated/

    TUsignal.tau_Emax
    TUsignal.toa_max2
    TUsignal.toa_new
    TUsignal.toa_win
    TUsignal.toa_max
    TUsignal.toa_th
    TUsignal.toa_cum
    TUsignal.toa_th_tmtm
    TUsignal.toa_th_tm
    TUsignal.toa_th_tmt
    TUsignal.toa_cum_tm
    TUsignal.toa_cum_tmtm
    TUsignal.toa_cum_tmt

Input Output Functions
-----------------------

.. autosummary::
    :toctree: generated/

    TUsignal.readcir
    TUsignal.readuwb

TUDsignal Class
===============

.. autosummary::
    :toctree: generated/

    TUDsignal.__init__
    TUDsignal.__repr__
    TUDsignal.fig

FBsignal Class
==============

.. autosummary::
    :toctree: generated/

    FBsignal.__init__
    FBsignal.__repr__
    FBsignal.plotri
    FBsignal.plot
    FBsignal.plotdB
    FBsignal.stem

FUsignal Class
==============

.. autosummary::
    :toctree: generated/

    FUsignal.__init__
    FUsignal.__repr__
    FUsignal.__add__
    FUsignal.__sub__
    FUsignal.__mul__
    FUsignal.__div__
    FUsignal.window
    FUsignal.get
    FUsignal.info
    FUsignal.energy
    FUsignal.applyFriis
    FUsignal.enthrsh
    FUsignal.dBthrsh
    FUsignal.zp
    FUsignal.newdf
    FUsignal.dftresamp
    FUsignal.resample
    FUsignal.symH
    FUsignal.symHz
    FUsignal.align
    FUsignal.ifft
    FUsignal.ift
    FUsignal.iftshift
    FUsignal.show
    FUsignal.decimate
    FUsignal.tap

FUDsignal Class
===============

.. autosummary::
    :toctree: generated/

    FUDsignal.__init__
    FUDsignal.__repr__
    FUDsignal.minphas
    FUDsignal.totime
    FUDsignal.iftd
    FUDsignal.ft1
    FUDsignal.ftau
    FUDsignal.cir
    FUDsignal.plot3d
    FUDsignal.ft2


FUDAsignal Class
===============

.. autosummary::
    :toctree: generated/

    FUDAsignal.__init__
    FUDAsignal.__repr__
    FUDAsignal.cut
    FUDAsignal.sort
    FUDAsignal.showtap
    FUDAsignal.tap

FHsignal Class
================

.. autosummary::
    :toctree: generated/

    FHsignal.__init__
    FHsignal.__repr__
    FHsignal.__mul__
    FHsignal.ifft
    FHsignal.unrex

Noise Class
===========

.. autosummary::
    :toctree: generated/

    Noise.__init__
    Noise.amplify
    Noise.gating

EnImpulse Class
===============

.. autosummary::
    :toctree: generated/

    EnImpulse.__init__
    EnImpulse.demo

MaskImpulse Class
=================

.. autosummary::
    :toctree: generated/

    MaskImpulse.__init__
    MaskImpulse.show

"""
import doctest
import os
import pdb
import numpy as np
import scipy as sp
import scipy.interpolate as interp
import numpy.fft as fft
import pandas as pd
import matplotlib.gridspec as gridspec
from copy import *
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from pylayers.util.pyutil import *
from pylayers.util.plotutil import *
import scipy.io as ios
from scipy.signal import cspline1d, cspline1d_eval, iirfilter, iirdesign, lfilter, firwin , correlate
#from sklearn import mixture
import scipy.stats as st


class Bsignal(PyLayers):
    r""" Signal with an embedded time base

    This class gathers a 1D signal and its axis indexation.

    The x base is not necessarily uniform

    x can have 1 or two axis

    The first axis of x and y have the same length

    By construction shape(y)[1] :=len(x), len(x) takes priority in case of observed conflict

    """

    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        r""" object constructor

        Parameters
        ----------

        x : ndarray (,Nx)
            primary axis (time or frequency)

        y : ndarray
            values  (...,Nx)
            the number of dimensions of y is arbitrary.
            the last dimension of y must be the primary axis


        """
        self.x = x
        self.y = y
        ndim = self.y.ndim
        # default naming of the axix
        if label==[]:
            self.label=[]
            for k in range(ndim):
                self.label.append('ax'+str(k))
        else:
            self.label=label
        if ndim > 1:
            shy = np.shape(self.y)
            shx = np.shape(self.x)
            lx  = shx[0]
            ly  = shy[-1]
            # last dimension of y should be equal to dimension of x
            if (ly != 1):
                if (ly != lx) :
                    print "Error in Bsignal : Dimension incompatibility "
                    print "x : ", lx
                    print "y : ", ly


    def __repr__(self):
        st = '%s :  %s  %s ' % (
                            self.__class__.__name__,
                            str(np.shape(self.x)),
                            str(np.shape(self.y)))

        for k in range(self.y.ndim):
            st = st + '\n' +self.label[k]+ ' : ' + str(self.y.shape[k])

        return(st)


    def extract(self,u):
        r""" extract a subset of signal from index

        Parameters
        ----------

        u : np.array

        Returns
        -------

        O : Usignal

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import numpy as np
            >>> x = np.arange(0,1,0.01)
            >>> y = np.sin(2*np.pi*x)
            >>> s = Bsignal(x,y)
            >>> su = s.extract(np.arange(4,20))
            >>> f,a = s.plot()
            >>> f,a = su.plot()

        """
        O = copy(self)
        O.x = O.x[u]
        try:
            O.y = O.y[u,:]
        except:
            O.y = O.y[u]
        return(O)

    def save(self, filename):
        r""" save Bsignal in Matlab File Format

        Parameters
        ----------

        filename : string

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import matplotlib.pyplot as plt
            >>> e = EnImpulse(fe=100)
            >>> fig,ax = e.plot(typ=['v'])
            >>> tit1 = plt.title('original waveform')
            >>> e.save('impulse.mat')
            >>> del e
            >>> h = TUsignal()
            >>> h.load('impulse.mat')
            >>> fig,ax = h.plot(typ=['v'])
            >>> tit2 = plt.title('retrieved waveform')

        See Also
        --------

        Bsignal.load

        """

        d = {}
        d['x'] = self.x
        d['y'] = self.y
        ios.savemat(filename, d)

    def load(self, filename):
        r""" load a Bsignal from a Matlab File

        Parameters
        ----------

        filename : string


        See Also
        --------

        Bsignal.save

        """
        d = ios.loadmat(filename)
        self.x = d['x'][0,:]
        self.y = d['y']

    def setx(self, x):
        r""" setx : set x vector

        Parameters
        ----------

        x : np.array

        Notes
        -----

        y is set to a zero vector

        Use __set__ instead

        .. todo:: check coherence of support internally

        """
        self.x = x
        Np = len(self.x)
        self.y = np.zeros(Np, dtype=float)

    def sety(self, fun):
        r""" sety : set y vector

        Parameters
        ----------

        fun : function

        """
        self.y = fun(self.x)

    def stem(self, **kwargs):
        r""" stem display

        Parameters
        ----------

        color : string
            default 'b-'

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import matplotlib.pyplot as plt
            >>> si = Bsignal()
            >>> si.x= np.arange(100)
            >>> si.y= np.arange(100)
            >>> f,a = si.stem()

        """

        ndim = self.y.ndim

        if 'fig' not in kwargs:
            fig = plt.figure()
        else:
            fig = kwargs['fig']

        if 'ax' not in kwargs:
            ax = fig.gca()
        else:
            ax = kwargs['ax']

        if ndim > 1:
            nl = len(self.y)
            for k in range(nl):
                ax.stem(self.x, self.y[k], **kwargs)
        else:
            ax.stem(self.x, self.y,**kwargs)
        return(fig,ax)

    def step(self, color='b'):
        r""" plot steps display

        Parameters
        ----------

        color : string
            default 'b-'

        """
        ndim = self.y.ndim
        if ndim > 1:
            nl = len(self.y)
            for k in range(nl):
                plt.plot(self.x, self.y[k], color, linestyle='steps')
        else:
            plt.plot(self.x, self.y, color, linestyle='steps')

    def cformat(self,**kwargs):
        r""" complex format

        Parameters
        ----------

        sax : list
            selected output axis from y default [0,1]
        typ : string
            ['l10','l20','d','r','du','ru','m','re','im']
        sel  : list of ndarray()
            data selection along selected axis, all the axis void 
            default [[],[]]
        ik : fixed axis value default (0)

        Returns
        -------

        a0 : first data axis
            this axis can be self.x or not
        a1 : second data axis
        dt : data

        ylabels : string
            label for the selected complex data format

        Notes
        -----

        This function returns 2 arrays x and y and the corresponding labels
        Convention : the last axis of y has same dimension as x
        y can have an arbitrary number of axis i.e a MIMO channel matrix
        could be x :  f and y : t x r x f

        Examples
        --------

        >>> import numpy as np
        >>> S = Bsignal()
        >>> x = np.arange(100)
        >>> y = np.arange(400).reshape(2,2,100)+1j*np.arange(400).reshape(2,2,100)
        >>> S.x = x
        >>> S.y = y
        >>> S.cformat()
        (array([[-240.        ,   43.01029996],
               [  49.03089987,   52.55272505]]), 'Magnitude (dB)')



        """
        defaults = {'sax' : [0,1],
                    'typ':'l20',
                    'sel' : [[],[]],
                    'ik' : 0
                   }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        # shape of y and number of dimensions

        shy = self.y.shape
        naxy = len(shy)

        # 2 axis selection
        ax = kwargs['sax']
        assert(len(ax)==2)

        # selected range on axis
        sel = kwargs['sel']
        ik  = kwargs['ik']

        # data selection on first axis
        if sel[0]==[]:
            # whole range
            us = 0
            uf = shy[ax[0]]
        else:
            # subset of values
            us = sel[0][0]
            uf = sel[0][-1]

        # data selection on second axis
        if sel[1]==[]:
            # whole range
            vs = 0
            vf = shy[ax[1]]
        else:
            # subset of values
            vs = sel[1][0]
            vf = sel[1][-1]


        # radians to degree conversion
        rtd = 180./np.pi

        t = kwargs['typ']

        if naxy ==1:
            dt = self.y[us:uf]
        elif naxy==2:
            dt = self.y[us:uf,vs:vf]
        elif naxy==3:
            if ((ax[0]==0) & (ax[1]==1)):
                dt = self.y[us:uf,vs:vf,ik]
            if ((ax[0]==0) & (ax[1]==2)):
                dt = self.y[us:uf,ik,vs:vf]
            if ((ax[0]==1) & (ax[1]==2)):
                dt = self.y[ik,us:uf,vs:vf]

        if t=='m':
            ylabels='Magnitude'
            dt = np.abs(dt)
        if t=='v':
            ylabels='Amplitude'
            dt = dt
        if t=='l10':
            ylabels='Magnitude (dB)'
            dt = 10*np.log10(np.abs(dt)+1e-12)
        if t=='l20':
            ylabels='Magnitude (dB)'
            dt = 20*np.log10(np.abs(dt)+1e-12)
        if t=='d':
            ylabels='Phase (deg)'
            dt = np.angle(dt)*rtd
        if t=='r':
            ylabels='Phase (rad)'
            dt = np.angle(dt)
        if t=='du':
            ylabels='Unwrapped Phase (deg)'
            dt = np.unwrap(np.angle(dt))*rtd
        if t=='ru':
            ylabels='Unwrapped Phase (rad)'
            dt = np.unwrap(np.angle(dt))
        if t=='re':
            ylabels='Real part'
            dt = np.real(dt)
        if t=='im':
            ylabels='Imaginary part'
            dt = np.imag(dt)
        if t=='gdn':
            ylabels='Group delay (ns)'
            df  = self.x[1]-self.x[0]
            dt  = -np.diff(np.unwrap(np.angle(dt)))/(2*np.pi*df)
        if t=='gdm':
            ylabels='Group distance (m)'
            df  = self.x[1]-self.x[0]
            dt = -0.3*np.diff(np.unwrap(np.angle(dt)))/(2*np.pi*df)
        if 'ylabels'  in kwargs:
            ylabels = kwargs['ylabels']

        return(dt,ylabels)

    def imshow(self,**kwargs):
        r""" imshow of y matrix

        Parameters
        ----------

        interpolation : string
            'none'|'nearest'|'bilinear'
        cmap : colormap
            plt.cm.BrBG
        aspect : string
            'auto' (default) ,'equal','scalar'
        function : string
            {'imshow'|'pcolormesh'}

        Examples
        --------


        >>> f = np.arange(100)
        >>> y = np.random.randn(50,100)+1j*np.random.randn(50,100)
        >>> F = FUsignal(f,y)

        """

        defaults = {'interpolation':'none',
                    'cmap':plt.cm.BrBG,
                    'aspect':'auto',
                    'fontsize':20,
                    'typ':'l20',
                    'function':'imshow',
                    'sax':[0,1]}

        for k in defaults.keys():
            if not kwargs.has_key(k):
                kwargs[k]=defaults[k]

        if not kwargs.has_key('fig'):
            fig = plt.figure()
        else:
            fig = kwargs['fig']

        if not kwargs.has_key('ax'):
            ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']

        # axis selection
        sax = kwargs['sax']

        if self.y.ndim>1:
            dt,ylabels = self.cformat(**kwargs)

            if 'vmin' not in kwargs:
                vmin = dt.min()
            else:
                vmin = kwargs['vmin']

            if 'vmax' not in kwargs:
                vmax = dt.max()
            else:
                vmax = kwargs['vmax']

            if kwargs['function']=='imshow':
                im = ax.imshow(dt,
                           origin = 'lower',
                           vmin = vmin,
                           vmax = vmax,
                           aspect = kwargs['aspect'],
                           #extent = (0,dt.shape[sax[0]],0,dt.shape[sax[1]]),
                           interpolation=kwargs['interpolation'],
                           cmap = kwargs['cmap'],
                           )

            ll = ax.get_xticklabels()+ax.get_yticklabels()
            for l in ll:
                l.set_fontsize(kwargs['fontsize'])


            if kwargs['function'] =='pcolormesh':
                im = ax.pcolormesh(xn,np.arange(dt.shape[0]),dt)

            cb = fig.colorbar(im)
            cb.set_label(ylabels,size=kwargs['fontsize'])

            for t in cb.ax.get_yticklabels():
                t.set_fontsize(kwargs['fontsize'])

            plt.axis('auto')
            fig.tight_layout()

            return fig,ax

    def plot(self, **kwargs):
        r""" plot signal

        Parameters
        ----------

        iy    : index of the waveform to plot (-1 = all)
        col   : string
            default  'black'
        vline : ndarray
        hline : ndarray
        unit1 : string
            default 'V'
        unit2 : string
            default 'V'
        xmin  : float
        xmax  : float
        ax    : axes instance or []
        dB    : boolean
            default False
        dist  : boolean
            default False
        display : boolean
            default True
        logx  : boolean
            defaut False
        logy  : boolean
            default False

        """

        defaults = {'iy'  :  -1,
                  'vline' : np.array([]),
                  'hline' : np.array([]),
                  'unit1' : 'V',
                  'unit2' : 'V',
                  'separated' : True,
                  'dist'  : False ,
                  'xmin'  :-1e15,
                  'xmax'  : 1e15,
                  'logx'  : False,
                  'logy'  : False,
                 }

        for key, value in defaults.items():
            if key not in kwargs:
                 kwargs[key] = value


        vline = kwargs['vline']
        hline = kwargs['hline']

        # filtering kwargs argument for plot function
        args ={}
        for k in kwargs:
            if k not in defaults.keys():
                args[k]=kwargs[k]

        conversion = 1.0
        if ((kwargs['unit1'] == 'V') & (kwargs['unit2'] == 'mV')):
            conversion = 1000

        # restriction of x support
        u = np.nonzero((self.x > kwargs['xmin']) & (self.x < kwargs['xmax']))[0]

        #
        # conver ns in meter if dist=True
        #
        if kwargs['dist']:
            x = 0.3 * self.x[u]
        else:
            x = self.x[u]

        ndim = self.y.ndim

        #
        # if ndim(y) > 1
        #
        if ndim > 1:
            yx = self.y[...,u]
            fig,ax = mulcplot(self.x[u],yx*conversion,**args)
        else:
            fig,ax = mulcplot(self.x[u],self.y[u]*conversion,**args)
        #
        # Draw vertical and horizontal lines
        #
        # To be added in mulcplot
        #

        tcolor = ['red', 'green', 'green', 'green', 'black', 'black', 'black']
        nl,nc = np.shape(ax)
        for l in range(nl):
            for c in range(nc):
                for i in range(len(vline)):
                    ax[l,c].axvline(vline[i], color=tcolor[i])
                for i in range(len(hline)):
                    ax.axhline(hline[i] * conversion, color='red')

        return(fig,ax)

    def flatteny(self,yrange=[],reversible=False):
        r""" flatten y array
        Parameters
        ----------
        yrange : array of y index values to be flattenned
        reversible : boolean
            if True the sum is place in object member yf
            else y is smashed
        """
        if self.y.ndim>1:
            if yrange==[]:
                yrange = np.arange(self.y.shape[0])
            if reversible:
                self.yf = np.sum(self.y[yrange,:],axis=0)
            else:
                self.y = np.sum(self.y[yrange,:],axis=0)

    def gating(self, xmin, xmax):
        r""" gating between xmin and xmax

        Parameters
        ----------

        xmin : float
        xmax : float

        Returns
        -------
            nothing self is modified

        Examples
        --------

        .. plot::
            :include-source:

            >>> import numpy as np
            >>> import matplotlib.pyplot as plt
            >>> from pylayers.signal.bsignal import *
            >>> x = np.linspace(-10,10,100)
            >>> y = np.sin(2*np.pi*12*x)+np.random.normal(0,0.1,len(x))
            >>> s = TUsignal(x,y)
            >>> fig,ax = s.plot(typ=['v'])
            >>> txt1 = plt.title('before gating')
            >>> plt.show()
            >>> s.gating(-3,4)
            >>> fig,ax=s.plot(typ=['v'])
            >>> txt2 = plt.title('after gating')
            >>> plt.show()

        Warnings
        --------

            When a gating is applied the removed information is lost

        """

        u = np.nonzero((self.x > xmin) & (self.x < xmax))[0]
        Nx = len(self.x)
        shy = np.shape(self.y)
        if len(shy) > 1:
            if shy[0] == Nx:
                Ny = shy[1]
                ind = 0
            else:
                Ny = shy[0]
                ind = 1

        w = np.zeros(Nx)
        w[u] = np.ones(len(u))
        if len(shy) > 1:
            if ind == 0:
                w = np.outer(w, np.ones(Ny))
            else:
                w = np.outer(np.ones(Ny), w)
        self.y = self.y * w

    def len(self):
        r""" returm length of x axis
        """
        return(len(self.x))

class Usignal(Bsignal):
    r""" Signal with an embedded uniform Base

        This class inheritate from Bsignal. The only difference
        is that the x base is supposed to be uniform

    """

    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(Usignal,self).__init__(x,y,label)
        #Bsignal.__init__(self, x, y)

    def __repr__(self):
        s = Bsignal.__repr__(self)
        return(s)

    def __add__(self, u):
        t = type(u).__name__
        if ((t == 'int') | (t == 'float')):
            U = Usignal(self.x, self.y + u)
        else:
            assert (u.y.ndim == 1) | (u.y.shape[0]==1)
            L = self.align(u)
            #u1 = L[0]
            #u2 = L[1]
            val = L.y[0:-1,:] + L.y[-1,:]
            # handle one dimensional array
            if ((val.ndim>1) & (val.shape[0]==1)):
                val = val.reshape(val.shape[1])
            U = Usignal(L.x, val)
        return(U)

    def __sub__(self, u):

        t = type(u).__name__
        if ((t == 'int') | (t == 'float')):
            U = Usignal(self.x, self.y - u)
        else:
            assert  (u.y.ndim == 1) | (u.y.shape[0]==1)
            L = self.align(u)
            #u1 = L[0]
            #u2 = L[1]
            val = L.y[0:-1,:] - L.y[-1,:]
            if ((val.ndim>1) & (val.shape[0]==1)):
                val = val.reshape(val.shape[1])
            U = Usignal(L.x,val)
            #U = Usignal(u1.x, u1.y - u2.y)
        return(U)

    def __mul__(self, u):
        t = type(u).__name__
        if ((t == 'int') | (t == 'float') | (t== 'float64') ):
            U = Usignal(self.x, self.y * u)
        else:
            assert  (u.y.ndim == 1) | (u.y.shape[0]==1)
            L = self.align(u)
            #u1 = L[0]
            #u2 = L[1]
            #U = Usignal(u1.x, u1.y * u2.y)
            val = L.y[0:-1,:] * L.y[-1,:]
            if ((val.ndim>1) & (val.shape[0]==1)):
                val = val.reshape(val.shape[1])
            U = Usignal(L.x, val)
        return(U)

    def setx(self, start, stop, dx):
        r""" set the x array of the Usignal (y=0)

        Parameters
        ----------
        start  : float
        stop   : float
        dx     : float

        Examples
        --------

        >>> u = Usignal()
        >>> u.setx(0,10,0.1)

        """
        self.x = np.arange(start, stop, dx)
        Np = len(self.x)
        self.y = 1.0 * np.zeros(Np)

    def dx(self):
        r"""  get the time step of Usignal.x

        Examples
        ---------

        >>> from pylayers.signal.bsignal import *
        >>> u = Usignal()
        >>> u.setx(0,10,0.1)
        >>> assert(u.dx()==0.1)

        """

        return(self.x[1] - self.x[0])

    def width(self):
        r""" get the extension support of the Usignal

        width is conventionnaly equal to the difference between extremities + dx

        Returns
        -------

        w : float

        """
        w = self.x[-1] - self.x[0] + self.dx()
        return(w)

    def expand(self, a):
        """ expand the time support by a scale factor a

        Parameters
        ----------
        a  : float
            expansion factor

        Returns
        -------

        Usignal  : support extended signal

        return a new Usignal with expanded factor a

        """
        dx = self.dx()
        T1 = self.width()
        T2 = a * T1
        N3 = np.ceil((T2 - T1) / (2 * dx))
        s = self
        xstart = self.x[0] - N3 * dx
        xstop = self.x[-1] + N3 * dx
        sl = Usignal()
        sr = Usignal()
        sl.setx(xstart, self.x[0], dx)
        sr.setx(self.x[-1], xstop, dx)
        som = sl + s + sr
        return(som)

    def max(self):
        """  maximum of Usignal

        Returns
        -------

        m : maximum value of self.y

        """
        y = self.y
        nd = np.ndim(y)
        if nd == 1:
            m = y.max()
        else:
            m = y.max(axis=1)

        return (m)

    def min(self):
        """  minimum of Usignal

        Returns
        -------

        m : maximum value of self.y

        """
        y = self.y
        nd = np.ndim(y)
        if nd == 1:
            m = y.min()
        else:
            m = y.min(axis=1)

        return (m)

    def truncate(self, posmin, posmax):
        """ truncate USignal in range [posmin, posmax]

        Parameters
        ----------

        posmin  : float
        posmax  : float

        Returns
        -------

        Usignal

        """
        t = type(self).__name__
        ndim = self.y.ndim

        x_new = self.x[posmin:posmax]
        if ndim > 1:
            y_new = self.y[:, posmin:posmax]
        else:
            y_new = self.y[posmin:posmax]

        if t == 'Usignal':
            U = Usignal(x_new, y_new)
        if t == 'Tchannel':
            U = Usignal(x_new, y_new)
        if t == 'TUsignal':
            U = TUsignal(x_new, y_new)
        if t == 'FUsignal':
            U = FUsignal(x_new, y_new)
        if t == 'FUDsignal':
            U = FUDsignal(x_new, y_new, self.taud)

        #if 'U' not in locals():
            #pdb.set_trace()
        return(U)

    def align(self, u2):
        """ align two Usignal on a same base

        returns a list which contains the two aligned signals

        It is assume that both signal u1 and u2 share the same dx
        This function can be improved regarding time granularity

        Parameters
        ----------

        u2 : TUsignal

        Returns
        -------

        TUsignal y extended TU signal, x bases are adjusted

        Examples
        --------

        .. plot::
            :include-source:

            >>> import matplotlib.pylab as plt
            >>> from pylayers.signal.bsignal import *
            >>> i1 = EnImpulse()
            >>> i2 = EnImpulse()
            >>> i2.translate(-10)
            >>> i3 = i1.align(i2)
            >>> fig,ax=i3.plot(typ=['v'])
            >>> plt.show()

        """
        # copy object
        u1 = self
        u1_start = u1.x[0]
        u1_stop = u1.x[-1]
        u2_start = u2.x[0]
        u2_stop = u2.x[-1]
        # get dimensions
        if u1.y.ndim>1:
            N1 = self.y.shape[0]
            M1 = self.y.shape[1]
        else:
            N1 = 1
            M1 = len(u1.y)
            u1.y = u1.y.reshape(1,M1)

        if u2.y.ndim>1:
            N2 = u2.y.shape[0]
            M2 = u2.y.shape[1]
        else:
            N2 = 1
            M2 = len(u2.y)
            u2.y = u2.y.reshape(1,M2)

        bool1 = abs(u1_start - u2_start) < 1e-10
        bool2 = abs(u1_stop - u2_stop) < 1e-10

        bool = bool1 & bool2

        if (bool):
        # same x support
            L = Usignal(u1.x, np.vstack((u1.y,u2.y)))
        else:
        # different x support
            xstart = min(u1_start, u2_start)
            xstop = max(u1_stop, u2_stop)

            b1i = abs(xstart - u1.x[0]) < 1e-15
            b2i = abs(xstart - u2.x[0]) < 1e-15
            b1f = abs(xstop - u1.x[-1]) < 1e-15
            b2f = abs(xstop - u2.x[-1]) < 1e-15

            if (b1i & b2f):
            # u1 is left u2 is right
                dx = u1.dx()
                T = xstop - xstart
                N = int(np.floor(T / dx))
                x = xstart + dx * np.arange(N)
                Y1 = np.zeros((N1,N), dtype=float)
                Y2 = np.zeros((N2,N), dtype=float)
                yleft = u1.y
                yright = u2.y
                Nleft = min(N,M1)
                Nright = min(N,M2)
                Y1[:,0:Nleft] = yleft[:,0:Nleft]
                Y2[:,-Nright:] = yright[:,0:Nright]
                U1 = Usignal(x, Y1[:,0:N])
                U2 = Usignal(x, Y2[:,0:N])

            if (b2i & b1f):
            # u2 is left u1 is right
                dx = u1.dx()
                T = xstop - xstart
                N = int(np.floor(T / dx))
                x = xstart + dx * np.arange(N)
                Y1 = np.zeros((N1,N), dtype=float)
                Y2 = np.zeros((N2,N), dtype=float)
                yleft = u2.y
                yright = u1.y
                Nleft = min(N, M2)
                Nright = min(N, M1)
                Y2[:,0:Nleft] = yleft[:,0:Nleft]
                Y1[:,-Nright:] = yright[:,0:Nright]
                U1 = Usignal(x, Y1[:,0:N])
                U2 = Usignal(x, Y2[:,0:N])

            if (b1i & b1f):
            # u2 is included in u1
                U1 = u1
                x = u1.x
                indx = np.nonzero((x >= u2_start) & (x <= u2_stop))[0]
                U2 = Usignal(x, np.zeros((N2,len(x))))
                #pdb.set_trace()
                U2.y[:,indx] = u2.y[:, 0:np.shape(indx)[0]]

            if (b2i & b2f):
            # u1 is included in u2
                U2 = u2
                x = u2.x
                indx = np.nonzero((x >= u1_start) & (x <= u1_stop))[0]
                U1 = Usignal(x, np.zeros((N1,len(x))))
                U1.y[:,indx] = u1.y

            #L = [U1, U2]
            L   = Usignal()
            L.x = U1.x
            L.y = np.vstack((U1.y,U2.y))
        return L

    def abs(self):
        """ return the absolute value of an Usignal

        Returns
        -------

        Usignal

        """
        A = Usignal(self.x, abs(self.y))
        return(A)

    def eprfl(self,axis=1):
        r""" Energy profile

        Parameters
        ----------

        axis : int

        Notes
        -----


        if axis==0

        $$\delta_x \sum_l |y(l,k)|^2$$

        if axis==1

        $$\delta_x \sum_k |y(l,k)|^2$$

        See Also
        --------

        cut

        """

        eprfl = np.real(np.sum(self.y * np.conj(self.y),axis=axis))*self.dx()
        return eprfl

    def energy(self,axis=0):
        """ calculate the energy of an Usignal

        Parameters
        ----------

        axis : int
            (default : 0)

        Returns
        -------

        energy : float

        """
        energy = self.dx() * np.sum(self.y * np.conj(self.y),axis=0)
        return(energy)

    def fftshift(self):
        self.y = fft.fftshift(self.y,axes=1)

    def zright(self, xmax):
        """ add zeros on the right until xmax

        Parameters
        ----------
        xmax : float

        """
        dx = self.dx
        cmax = max(self.x)
        assert xmax > cmax, 'Warning xmax <= cmax'
        xadd = np.linspace(cmax + dx, xmax, dx)
        yadd = np.zeros(len(xadd))
        self.x = np.hstack((self.x, xadd))
        self.y = np.hstack((self.y, yadd))

    def zleft(self, xmin):
        """ add zeros on the left  until xmin

        Parameters
        ----------
        xmin : float

        """
        dx = self.dx
        cmin = min(self.x)
        assert xmin < cmin, 'Warning xmin >= cmin'
        xadd = np.linspace(xmin, cmin - dx, dx)
        yadd = np.zeros(len(xadd))
        self.x = np.hstack((xadd, self.x))
        self.y = np.hstack((yadd, self.y))

    def zlr(self, xmin, xmax):
        """  add zeros to the left and to the right

        Parameters
        ----------

        xmin : float
            add zeros before xmin
        xmax : float
            add zeros after xmax

        Summary
        --------

        This corresponds to a gating between xmin and xmax

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> from matplotlib.pylab import *
            >>> ip = EnImpulse()
            >>> f,a = ip.plot(typ=['v'])
            >>> ip.zlr(-10,10)
            >>> f,a = ip.plot(typ=['v'])

        """
        dx = self.dx()
        cmin = min(self.x)
        cmax = max(self.x)
        if cmin > xmin:
            xaddl = np.arange(xmin, cmin - dx, dx)
            yaddl = np.zeros(len(xaddl))
            self.x = np.hstack((xaddl, self.x))
            self.y = np.hstack((yaddl, self.y))
        else:
            u = np.nonzero(self.x >= xmin)
            self.x = self.x[u[0]]
            self.y = self.y[u[0]]

        if cmax < xmax:
            xaddr = np.arange(cmax + dx, xmax, dx)
            yaddr = np.zeros(len(xaddr))
            self.x = np.hstack((self.x, xaddr))
            self.y = np.hstack((self.y, yaddr))
        else:
            u = np.nonzero(self.x <= xmax)
            self.x = self.x[u[0]]
            self.y = self.y[u[0]]

class TBsignal(Bsignal):
    """  Based signal in Time domain


    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(TBsignal,self).__init__(x,y,label)
        #Bsignal.__init__(self, x, y)
        self.label[-1] = 'time (ns)'

    def __repr__(self):
        s = Bsignal.__repr__(self)
        return(s)

    def plot(self,**kwargs):
        """ plot TBsignal

        Parameters
        ----------

        unit1 : actual unit of data
        unit2 : unit for display

        Examples
        --------

        >>> import numpy as np
        >>> from pylayers.signal.bsignal import *
        >>> from matplotlib.pylab import *
        >>> x = np.array( [ 1, 3 , 6 , 11 , 18])
        >>> y = np.array( [ 0,1 ,-5, 8 , 10])
        >>> s = Bsignal(x,y)
        >>> fi = figure()
        >>> fig,ax=s.plot(typ=['v'])
        >>> ti = title('TBsignal : plot')
        >>> show()

        """
        defaults = {'iy'  :  -1,
                  'vline' : np.array([]),
                  'hline' : np.array([]),
                  'unit1' : 'V',
                  'unit2' : 'V',
                  'dist'  : False ,
                  'xmin'  :-1e15,
                  'xmax'  : 1e15,
                  'logx'  : False,
                  'logy'  : False
                 }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

#        if tmin != -1e5:
#            xmin = tmin
#        else:
#            xmin = -1e5
#
#        if tmax != 1e5:
#            xmax = tmax
#        else:
#            xmax = tmax

        if kwargs['dist']:
            kwargs['xlabels']=['distance (m)']
        else:
            kwargs['xlabels']=['Time (ns)']


#        if kwargs['unit2'] == 'mV':
#            if kwargs['type']=='l20':
#                kwargs['ylabel']=['Voltage (dBmV)']
#            if kwargs['type']=='v':
#                kwargs['ylabel']=['Voltage (mV)']
#
#        if kwargs['unit2'] == 'V':
#            if kwargs['type']=='l20':
#                kwargs['ylabel']=['Voltage (dBmV)']
#            if kwargs['type']=='v':
#                plt.ylabel('Voltage (mV)')

        fig,ax = Bsignal.plot(self,**kwargs)
        #fig,ax = self.plot(**kwargs)

        return(fig,ax)

    def translate(self, tau):
        """  translate signal by tau

        Parameters
        ----------
        tau : float
            delay for translation

        Warnings
        --------

        Once translated original signal and translated signal might not be on the same grid

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> from matplotlib.pylab import *
            >>> ip = EnImpulse()
            >>> ip.translate(-10)
            >>> fig,ax=ip.plot(typ=['v'])
            >>> show()


        """
        self.x = self.x + tau

    def b2u(self, N):
        """ conversion into an uniform signal

        Parameters
        ----------

        N : integer
            Number of points

        Returns
        -------

        U : TUsignal


        Notes
        -----

        This function exploits linear interp1d

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import matplotlib.pyplot as plt
            >>> x = np.array( [ 1, 3 , 6 , 11 , 18])
            >>> y = np.array( [ 0,1 ,-5, 8 , 10])
            >>> sb = TBsignal(x,y)
            >>> su20 = sb.b2u(20)
            >>> su100 = sb.b2u(100)
            >>> fi = plt.figure()
            >>> st = sb.stem()
            >>> fig,ax = su20.plot(color='k')
            >>> fig,ax = su100.plot(color='r')
            >>> ti = plt.title('b2u : sb(blue) su20(black) su200(red)')

        """

        fi = interp.interp1d(self.x, self.y, kind='linear')
        xn = np.linspace(self.x[0], self.x[-1], N)
        yn = fi(xn)
        U = TUsignal(xn, yn)

        return U


class TUsignal(TBsignal, Usignal):
    """ Uniform signal in Time domain

    This class inheritates from TBsignal and Usignal

    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(TUsignal,self).__init__(x,y,label)
        #Usignal.__init__(self, x, y)

    def __repr__(self):
        s = Usignal.__repr__(self)
        return(s)

    def diff(self):
        """ numerical differentiation

        Warnings
        ---------
            the number of points is reduced by one

        Examples
        --------

        .. plot::
            :source-include:

            #>>> import numpy as np
            #>>> from pylayers.signal.bsignal import *
            #>>> from matplotlib.pyplot import *
            #>>> x  = np.arange(0,1.0,0.1)
            #>>> y  = x**3
            #>>> su = TUsignal(x,y)
            #>>> dsu   = su.diff()
            #>>> ddsu  = dsu.diff()
            #>>> dddsu = ddsu.diff()
            #>>> fi = plt.figure()
            #>>> fig,ax=su.plot(typ=['v'] ,color='k')
            #>>> fig,ax=dsu.plot(typ=['v'],color='g',fig=fig,ax=ax)
            #>>> fig,ax=ddsu.plot(typ=['v'],color='r',fig=fig,ax=ax)
            #>>> fig,ax=dddsu.plot(typ=['v'],color='b',fig=fig,ax=ax)
            #>>> ti = plt.title('TUsignal : diff')
            #>>> plt.show()

        """
        V = TUsignal()
        V.y = np.diff(self.y)
        V.x = self.x[0:-1]
        return(V)

    def info(self):
        """ display information

        """
        print 'TUsignal'
        print '--------'
        print 'shx : ', np.shape(self.x)
        print 'shy : ', np.shape(self.y)
        print 'dx :  ', self.dx()
        print 'xmin :', self.x.min()
        print 'xmax :', self.x.max()
        print 'ymin :', self.y.min()
        print 'ymax :', self.y.max()


    def awgn(self,PSDdBmpHz=-174,snr=0,seed=1,typ='psd',R=50):
        """ add a white Gaussian noise

        Parameters
        ----------

        PSDdBmpHz : float
        seed : float
        typ : string
            'psd' | 'snr'

        Returns
        -------
        n
        sn

        See Also
        --------

        bsignal.Noise

        """
        ti = self.x[0]
        tf = self.x[-1]
        tsns = self.x[1]-self.x[0]
        fsGHz = 1./tsns

        if typ=='snr':
            Ps = self.energy()/(R*(tf-ti))
            PW = Ps/10**(snr/10.)
            pWpHz = PW/(fsGHz*1e9)
            pmWpHz = pWpHz*1e3
            PSDdBmpHz = 10*np.log10(pmWpHz)

        n = Noise(ti = ti,
                   tf = tf+tsns,
                   fsGHz = fsGHz,
                   PSDdBmpHz = PSDdBmpHz,
                   R = R,
                   seed = seed)


        sn = TUsignal()
        sn.y = self.y + n.y
        sn.x = self.x

        return sn

    def fft(self, shift=False):
        """  forward fast Fourier transform


        Parameters
        ----------
        shift : boolean
            default False

        Warnings
        --------

        This fft is a scaled fft and takes into account the value of the sampling period.

        Returns
        -------

        FHsignal

        Examples
        --------

        >>> e = EnImpulse()
        >>> E = e.fft()

        """
        Np = len(self.x)
        te = self.x[1] - self.x[0]
        fe = 1.0 / te
        f = np.linspace(0, fe, Np, endpoint=False)
        y = fft.fft(self.y)
        S = FHsignal()
        S.x = f
        S.y = y
        S.y = S.y * te
        # compensation du retard 11.12.2009
        #S.y=S.y*exp(-1j*2*np.pi*f*self.x[0])
        return(S)

    def fftsh(self):
        """ return an FHsignal

        Warnings
        --------

        This fft is a scaled fft and takes into account the value of the
        sampling period.

        Returns
        -------

        FHsignal : Frequency signal with enforced Hermitian symetry

        """
        Np = len(self.x)
        te = self.x[1] - self.x[0]
        fe = 1.0 / te
        f = np.linspace(0, fe, Np, endpoint=False)
        #y = fft(self.y)
        #y = fftshift(y)
        y = fft.fftshift(self.y)
        y = fft.fft(y)
        S = FHsignal()
        S.x = f
        S.y = y
        S.y = S.y * te
        return(S)

    def filter(self, order=4, wp=0.3, ws=0.4, ftype='butter'):
        """ signal filtering

        Parameters
        ----------
        order : int
        wp    : float
        ws    : float
        ftype : string

        Returns
        -------
        O : Output filetered TUsignal


        """
        O = TUsignal()
        #ba  = iirfilter(order,[wp,ws],ftype=ftype)
        print "wp = ", wp
        print "ws = ", ws
        #ba  = iirdesign(wp,ws,1,40,ftype=ftype)
        h = firwin(1001, [wp, ws])
        O.y = lfilter(h, 1, self.y)
        O.x = self.x
        return(O)

    def ftshift(self):
        """ return the associated FUsignal

        Returns
        -------

        H : FUsignal

        See Also
        --------

        pylayers.signal.bsignal.TUsignal.fftsh
        pylayers.signal.waveform.ip_generic

        """
        A  = self.fftsh()
        AU = A.unrex()
        return(AU)

    def psd(self, Tpns=100, R=50,periodic=True):
        """ calculate power spectral density

        Parameters
        ----------

        R    : Resistance (default 50 Ohms)
            Ohms
        Tpns : real
            Signal period PRP (default 100 ns)

        .. note::

            Notice this interesting property that if time is represented in ns
            the resulting PSD is expressed in dBm/MHz because there is the
            same scale factor 1e-9 between second and nanosecond as between
            dBW/Hz and dBm/MHz

            If periodic is False the signal duration is taken as period.

        """
        P = self.esd(mode='unilateral')
        if periodic:
            P.y = P.y / (R * Tpns)
        else:
            P.y = P.y/ (R* (P.x[-1]-P.x[0]))
        return(P)

    def show(self,fig=[],ax=[],display=True,PRPns=100):
        """  show psd

        Parameters
        ----------

        fig
        ax
        display

        See Also
        ---------

        pylayers.signal.bsignal.TUsignal.psd
        pylayers.signal.bsignal.Usignal.zlr

        """
        if fig == []:
            fig = plt.gcf()

        Y = self.psd(PRPns)
        kboltzmann = 1.3806503e-23
        T = 300
        # +30 dBW/Hz -> dBm/MHz
        N0dB = 10*np.log10(kboltzmann*T)+30

        u      = np.nonzero(Y.y>N0dB)
        Pnoise = sum(Y.y[u])*Y.dx()

        ax1 = fig.add_subplot(121)
        self.plot()
        ax1.axhline(3*np.sqrt(Pnoise),color='r',linewidth=2)
        ax1.axhline(-3*np.sqrt(Pnoise),color='r',linewidth=2)
        plt.grid()
        ax2 = fig.add_subplot(122)
        #self.zlr(0,150)
        # unilateral
        plt.plot(Y.x, 10 * np.log10(Y.y), 'k')
        plt.xlabel('Frequency (GHz)')
        plt.xlim([1,10])
        ax2.axhline(N0dB,color='r',linewidth=2)
        plt.ylabel('PSD (dBm/MHz) assuming PRP='+str(PRPns)+' ns')
        plt.grid()
        if display:
            plt.show()

        return(fig,[ax1,ax2])

    def esd(self, mode='bilateral'):
        """  Calculate the energy spectral density of the U signal

        Parameters
        ----------

        mode : string
            'unilateral' | 'bilateral'

        Returns
        -------

        FHsignal : if mode == 'bilateral'
        FUsignal : if mode == 'unilateral'

        See Also
        --------

        FHsignal.unrex

        """
        te = self.dx()
        Y = self.fft()
        #YY  = (te*te)*abs(Y.y)**2
        YY = abs(Y.y) ** 2
        O = FHsignal(Y.x, YY)
        if mode == 'unilateral':
            Ou = O.unrex()
            Ou.y = 2 * Ou.y
            return(Ou)
        return(O)





    def shift(self, tau):
        """ shift the Usignal by tau (expressed in the same unit as Bsignal.x)

        This method uses fft and ifft to produce a shift which is not necessarily a
        multiple of the time step dx.

        Parameters
        ----------

        tau : float

        """
        t = self.x
        S = self.fft()
        N = len(t)
        f = S.x
        Y = S.y
        arg = (2 * np.pi * f * tau) % (2 * np.pi)
        E = np.exp(-1j * arg)
        SD = Y * E
        cc = np.array([SD[0]])
        ze = np.array([0])
        if N % 2 == 0:
            U = SD[1:N / 2 + 1]
            V = np.concatenate((cc, U, np.flipud(np.conjugate(U[0:-1]))), 1)
        else:
            U = SD[1:(N - 1) / 2 + 1]
            V = np.concatenate((cc, U, np.flipud(np.conjugate(U))), 1)
        self.y = np.real(ifft(V))
        return(self)

    def correlate(self, s, normalized = True):
        """ correlates with an other signal

        Parameters
        ----------

        s : TUsignal
        normalized : boolean
            default True

        Summary
        -------
            The time step dx should be the same

            The chosen time reference is the one of the self signal
            This means that this function can be used to detect events appearing
            in the self signal by correlating with signal ..math:`s`, maximum of
            the correlation function appears in the proper time of the self signal
            history.

            **Interpretation of the correlation levels**

            We assume that s1 is the self signal and s2 is the correlating signal

            s1 = a1 * sn(t-tau)   where sn is normalized in energy
            s2 = a2 * sn(t)             "   "   "   "  "

            the maximum of correlation appears at time tau and its value is

                    max (s1.correlate(s2,True)) = a1*a2

            In the case of the propagation simulation we have
                    s1 : the time domain UWB received waveform
                    s2 : the corresponding emitted waveform s2,

            To determine the ratio a1/a2 or a2/a1 in order to evaluate the propagation Loss L = 20 log10(a2/a1)

            This calculus assumes implicitely that the shape of the transmitted signal
            has not been modified during the propagation, which is a strong hypothesis.
            This calculus is then a minorant of the received energy.

            L = max(s1.correlate(s2,True)/s1.Energy)   with a1**2 = sum(s1**s1)*dx
    """

        x = self.x
        xcorr = np.correlate(self.y, s.y, mode='full')
        dx = self.dx()
        Ns = len(s.x)
        Nxc = len(xcorr)
        t_left = x[0] + dx * (-np.arange(Ns / 2))[-1:0:-1]
        t = np.concatenate((t_left, x))
        Nt = len(t)
        Nright = Nxc - Nt
        t_right = x[-1] + dx + dx * np.arange(Nright)
        t = np.concatenate((t, t_right))
        phi = TUsignal()
        phi.x = t
        if normalized:
            phi.y = xcorr * dx
        else:
            phi.y = xcorr
        return(phi)

    def corrgauss(self, sigma):
        """ Correlate TUsignal with a normalized gaussian of standard deviation sigma

        Parameters
        ----------
        sigma : float
            standard deviation of Gaussian

        """
        dx = self.dx()
        t = np.arange(-10 * sigma, 10 * sigma, dx)
        gn = TUsignal(t, 1. / (np.sqrt(2 * np.pi) * sigma) *
                      np.exp(-(t ** 2) / (2 * sigma ** 2)))
        pcorr = self.correlate(gn)
        fac = np.sqrt(2) * 1.15
        p = TUsignal(pcorr.x, pcorr.y * fac)
        return p

    def Efirst_loc(self, nint, E0):
        """ find the Efirst using the mean like toa_max

        Parameters
        ----------
        nint
        E0
        """
        x = self.x
        y = self.y
        M = max(self.y)
        step = M / 1e2
        thre = M - step
        while step > M / 1e5:
            u = np.nonzero(self.y > thre)[0]
            if nbint(u) < nint:
                thre = thre - step
            else:
                thre = thre + step
                step = step / 2

        w = u[1:] - u[0:-1]
        w0 = np.nonzero(w != 1)[0]
        vv = u[0:w0[0] + 1]
        ff = max(y[vv])
        Efirst = ff / E0
        Efirst = 20 * np.log10(Efirst)
        return(Efirst)

    def resample(self, x_new, kind='linear'):
        """ resample the TUsignal with a new x base

        x_new needs to be strictly included in the original x base of the Usignal.
        x is a 1D array
        y is a 2D array

        Parameters
        ----------
        x_new : ndarray
        kind  : string
            'linear' |'spline'
        """
        if kind == 'linear':
            f = interp.interp1d(self.x, self.y)
            y_new = f(x_new)

        if kind == 'spline':
            coef = splrep(self.x, self.y, s=0)
            y_new = splev(x_new, coef, der=0)

        U = TUsignal(x_new, y_new)
        return(U)

    def ft(self):
        """ return the associated FUsignal
        """
        A = self.fft()
        AU = A.unrex()
        return(AU)

    def convolve(self, u):
        """ time domain convolution
        """
        dx = u.dx()
        i0 = np.nonzero((u.x < dx) & (u.x > -dx))[0]
        ind0 = i0[0]
        N1 = len(self.x)
        N2 = len(u.x)
        t0 = self.x[0]
        t = np.arange(t0 - ind0 * dx, t0 + (N1 + N2 - 1 - ind0) * dx, dx)
        y = convolve(self.y, u.y)
        v = TUsignal(t[0:len(y)], convolve(self.y, u.y))
        return(v)

    def Yadd_zeros2l(self, N):
        """ time domain extension on the left with N zeros

        Parameters
        ----------

        N : integer
            number of additinal zero values

        See Also
        --------

        Yadd_zeros2r


        Work only for single y
        """
        te = self.dx()
        self.y = np.hstack((np.zeros(N), self.y))
        aux = np.arange(N) * te
        t0 = self.x[0]
        self.x = np.hstack((aux - (aux[-1] - t0 + te), self.x))

    def Yadd_zeros2r(self, N):
        """ time domain extension on right with N zeros

        Parameters
        ----------

        N : integer
            number of additinal zero values

        See Also
        --------

        Yadd_zeros2l

        """
        te = self.dx()
        self.y = np.hstack((self.y, np.zeros(N)))
        aux = np.arange(N) * te
        t1 = self.x[-1]
        self.x = np.hstack((self.x, aux - (aux[0] - t1 - te)))


    def Epercent(self, N=10):
        """ return N percentile delay of a cdf

        Parameters
        ----------

        N : 10

        """

        cdf, vary = self.ecdf()
        t = cdf.x
        Cdf = cdf.y
        pc = array([])
        for i in range(N - 1):
            u = np.nonzero(Cdf > (i + 1.) / N)
            tp = t[u[0]]
            pc = np.hstack((pc, tp))
        return(pc)

    def Etau0(self, tau0=0.0, Tint=1, sym=0.25, dB=True):
        """ calculate energy around delay tau0

        Parameters
        ----------
        tau0  : (ns)            (0)
        Tint  : Integration time (ns)   (1) include the system error
        sym   : symetrie factor 0.5 = symetric (0.25)
        dB    : logscale indicator (True)
        """
        #u  = nonzero((tau0 + Tint*(1-sym) > self.x) & (self.x > tau0 - Tint*sym))
        u = nonzero((tau0 + Tint > self.x) & (self.x > tau0))
        etau0 = self.dx() * sum(self.y[u] * np.conj(self.y[u]))
        if dB:
            etau0 = 10 * np.log10(etau0)
        return(etau0)

    def Ewin(self, tau, Tint=1, sym=0.25, dB=False):
        """  integrate energy around delay tau

        Parameters
        ----------
        tau   : (ns)            (0)
        Tint  : Integration time (ns)   (1) include the system error
        sym   : symetrie factor 0.5 = symetric (0.25)
        dB    : logscale indicator (True)

        """
        tstart = tau - Tint * sym
        tstop = tau + Tint * (1 - sym)
        u = np.nonzero((self.x > tstart) & (self.x < tstop))
        energy = self.dx() * sum(self.y[u] * np.conj(self.y[u]))
        if dB:
            energy = 10 * np.log10(energy)
        return(energy)

    def Etot(self, tau0=0.0, taumax=200, dB=False):
        """
        Etot  calculate the energy of the signal

        Parameters
        ----------
        tau0 : start value for integration
        dB   : (False default) if True value in dB

        usage  :

            s.Etot(tau0=10,dB=True)

        """
        u = (self.x > tau0) & (self.x < taumax)
        etot = self.dx() * sum(self.y[u] * np.conj(self.y[u]))
        if dB:
            etot = 10 * np.log10(etot)
        return(etot)

    def Efirst(self, toa, Tint=1, sym=0.25, dB=True):
        """ calculate the energy of the first path

        Parameters
        ----------
        toa  : float
            delay value
        Tint : float
            duration value (1)
        sym : float
            symmetry around delay value ( 0.25)
        dB : Boolean

        Returns
        -------

        Efirst : Energy amount in the window (in dB if dB)

        """
        u = np.nonzero((toa + Tint > self.x) & (self.x > toa))
        efirst = self.dx() * sum(self.y[u] * np.conj(self.y[u]))
        if dB:
            efirst = 10 * np.log10(efirst)
        return(efirst)

    def Efirst_corr(self, tau0, Sx, Sy, dB=True):
        """ calculate Efirst utilizing the correlation of signal emission et reponse impulsionnelle

       Parameters
       ----------
       tau0
       Sx
       Sy
       dB

        """
        te = self.dx()
        E0 = sum(Sy * Sy) * te
        n = int(np.ceil(tau0 / te))
        Correlation = np.correlate(self.y, Sy, mode='full')
        seuil = max(Correlation[len(Sx):len(Sx) + n - 200])
        v = np.nonzero(Correlation[len(Sx) + n - 200:] > seuil)[0]
        if len(v) == 0:
            ff = seuil / E0
        else:

            w = v[1:] - v[0:-1]
            w0 = np.nonzero(w != 1)[0]
            if len(w0) == 0:
                ff = max(Correlation[len(Sx) + n - 200:][v]) / E0
            else:
                vv = v[0:w0[0] + 1]
                ff = max(Correlation[len(Sx) + n - 200:][vv]) / E0

        if dB:
            Ef = 20 * np.log10(ff)

        return(Ef)

    def Efirst_toath(self, tau0, Tint=1, sym=0.25, dB=True):
        """ calculate Efirst

        Parameters
        ----------
        tau0   : Time of flight
        Tint
        sym
        dB   : if True return value in dBnJ

        """

        te = self.dx()
        n = int(np.ceil(tau0 / te))
        seuil = max(self.y[:n])
        v = np.nonzero(self.y[n:] > seuil)[0]
        if len(v) == 0:
            toa = n * te
        else:
            w = v[1:] - v[0:-1]
            w0 = np.nonzero(w != 1)[0]
            if len(w0) == 0:
                r = max(self.y[n:][v])
                toa = np.nonzero(self.y == r)[0] * te

            else:
                vv = v[0:w0[0] + 1]
                r = max(self.y[n:][vv])
                toa = np.nonzero(self.y == r)[0] * te

        u = np.nonzero((toa + Tint * (1 - sym) > self.x) & (
            self.x > toa - Tint * sym))
        efirst = te * sum(self.y[u] * np.conj(self.y[u]))

        if dB:
            efirst = 10 * np.log10(efirst)

        return(efirst)

    def taumax(self):
        r""" calculate taumax

        .. math::
            \max_{\tau} y^{2}(\tau)

        """

        y2 = (self.y) ** 2
        #
        # determine time of maximum value of ()^2
        #
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]
        tau = self.x[u]
        return(tau)

    def Emax(self, Tint=1, sym=0.5, dB=False):
        """ calculate the maximum of Energy integrated over a duration Tint

        A symetry of sym around the max value of the squared signal

        Parameters
        ----------

        Tint: float
            Integration time (ns) default 1
        sym : float
            Symmetry factor (default 0.5)
        dB  : boolean
            default False

        Notes
        -----

        W1-M1
        te     = 0.005 ns
        left  = 12
        Nright = 33
        Tint   = 45*te = 0.225 ns
        sym    = 0.25

        """
        #
        #  ( ) ^2
        #
        y2 = (self.y) ** 2
        #
        # determine time of maximum value of ()^2
        #
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]

        te = self.dx()

        Npt = int(np.ceil(Tint / te))
        Nleft = int(np.ceil(sym * Npt))
        Nright = int(np.ceil((1 - sym) * Npt))
        #
        #  Integration around the maximum value of E^2
        #  In the W1_M1 measurement
        #  te     = 0.005 ns
        #  Nleft  = 12
        #  Nright = 33
        #  Tint   = 45*te = 0.225 ns
        #  sym    = 0.25
        #

        Y = y2[u - Nleft:u + Nright]
        cumY = np.cumsum(Y)
        maxY = cumY[-1]
        Emax = maxY * te
        if dB:
            return(10 * np.log10(Emax))
        return(Emax)

    def tau_Emax(self):
        """ calculate the delay of max energy peak
        """
        y2 = (self.y) ** 2
        t = self.x
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]
        tau_Emax = t[u]
        return(tau_Emax)

    def toa_max2(self):
        """ calculate time of arrival max2 method
        """

        THRE = array([])
        V = array([])
        VL = array([])

        M = max(self.y)
        n = np.nonzero(self.y == M)[0]

        thre = M
        v = 1
        vl = 0
        THRE = np.hstack((THRE, thre))
        V = np.hstack((V, v))
        VL = np.hstack((VL, vl))

        step = M / 1e2
        thre = M - step

    #       while thre > M/1e2:
        while vl < 20:
    #       while v < 50:

            u = np.nonzero(self.y > thre)[0]
            v = nbint(u)
            h = np.nonzero(u > n)[0]
            g = delete(u, h)
            vl = nbint(g) - 1

            THRE = np.hstack((THRE, thre))
            V = np.hstack((V, v))
            VL = np.hstack((VL, vl))

            thre = thre - step

        plt.plot(1 - THRE / M, V, 'b', drawstyle='steps',
                 label='interval number')
        plt.plot(1 - THRE / M, VL, '-r', drawstyle='steps',
                 label='interval(Left) number')
        plt.xlabel('Gamma/Vmax')
        plt.legend(loc=2)
    #       ylabel('Interval Number')
        plt.show()

    def toa_new(self):
        """ calculate time of arrival (new method) 
        """
        t = self.x
        Max = max(self.y)
        nmax = np.nonzero(self.y == Max)[0]
        n = nmax
        step = Max / 1e2
        thre = Max - step

        delta = 100
        d = 0
        nint = 0
        N = np.array([])
        N = np.hstack((N, n))

        while delta > 4 * Max / 1e2:

            u = np.nonzero(self.y > thre)[0]
            hr = np.nonzero(u > n)[0]
            g = delete(u, hr)

            if nmax >= 6000:
            #set the fenetre=6000*0.005=30ns
                hl = np.nonzero(g < nmax - 6000)[0]
                u = delete(g, hl)
            else:
                u = g

            n_int = nbint(u) - 1

            if n_int == 0:
                d = d + step
            else:
                delta = d + step
                d = 0
                n = u[0]
                N = np.hstack((N, n))
                print N

            thre = thre - step
            if thre < 0:
                break
        if len(N) >= 3:
            nn = N[-3]
        else:
            nn = N[0]

        tau = t[nn]
        return tau

    def toa_win(self, w):
        """ calulate time of arrival (window method)

        Parameters
        ----------
        w : parameter between 0 and 100
        Lei takes w = 9

        """
        t = self.x
        maxbruit = max(self.y[0:1000])
        Max = max(self.y)
        nmax = np.nonzero(self.y == Max)[0]
        n = nmax
        step = Max / 1e2
        thre = Max - step

        delta = 100
        d = 0
        nint = 0
        N = np.array([])
        N = np.hstack((N, n))
        # tant delta est plus grande que w% du Max
        while delta > w * Max / 1e2:

            u = np.nonzero(self.y > thre)[0]
            hr = np.nonzero(u > n)[0]
            g = delete(u, hr)

            if nmax >= 6000:
            #set the fenetre=6000*0.005=30ns
                hl = np.nonzero(g < nmax - 6000)[0]
                u = delete(g, hl)
            else:
                u = g

            n_int = nbint(u) - 1

            if n_int == 0:
                thre = thre - step
                d = d + step
            else:
                delta = Max - maxbruit - d - step
                d = d + step
                n = u[0]
                N = np.hstack((N, n))
                thre = thre - step

            if thre < 0:
                break
        if len(N) >= 2:
            nn = N[-2]
        else:
            nn = N[0]

        tau = t[nn]
        return tau

    def toa_max(self, nint):
        """ calculate time of arrival

        descendant threshold based toa estimation

        Parameters
        ----------
        nint : integer

        """
        #
        # seek fot the maximum value of the signal
        #
        M = max(self.y)
        step = M / 1e2
    #       plot(self.x,self.y)
        thre = M - step
        while step > M / 1e5:
    #          axhline(y=thre,color='green')
            u = np.nonzero(self.y > thre)[0]
            if nbint(u) < nint:
            # down
                thre = thre - step
            else:
            # up + step reduction
                thre = thre + step
                step = step / 2

    #       plt.show()
        tau = self.x[u[0]]
        return tau

    def toa_th(self, thlos, thnlos, visibility=0):
        """ calculate time of arrival

        threshold based toa estimation using energy peak

        """
        #
        #  ( ) ^2
        #
        y2 = (self.y) ** 2
        maxy2 = max(y2)
        t = self.x

        if visibility == 'LOS':
            th = thlos * maxy2
        else:
            th = thnlos * maxy2
        #
        #In the W1-M1 measurement
        #thlos=0.05   thnlos=0.15
        #
        v = np.nonzero(y2 >= th)[0]
        toa = t[v[0]]
        return toa

    def toa_cum(self, th):
        """ calculate time of arrival

        threshold based toa estimation using cumulative energy
        """
        t = self.x
        y = self.y
        cdf, vary = self.ecdf()
        #
        #In the W1-M1 measurement th=0.15
        #
        v = np.nonzero(cdf.y >= th)[0]
        toa = t[v[0]]
        return toa

    def toa_th_tmtm(self):
        """ calculate time of arrival

        """
        y2 = (self.y) ** 2
        maxy2 = max(y2)
        t = self.x

        alpha = (np.sqrt(self.Etot()) - np.sqrt(self.Emax())) /  \
                (np.sqrt(self.Etot()) + np.sqrt(self.Emax()))
        th = alpha * maxy2

        v = np.nonzero(y2 >= th)[0]
        toa = t[v[0]]
        return toa

    def toa_th_tm(self):
        """ calculate time of arrival

        """

        y2 = (self.y) ** 2
        maxy2 = max(y2)
        t = self.x

        alpha = np.sqrt(self.Emax()) / np.sqrt(self.Etot())
        print alpha
        th = alpha * maxy2

        v = np.nonzero(y2 >= th)[0]
        toa = t[v[0]]
        return toa

    def toa_th_tmt(self):
        """ calculate time of arrival

        """
        y2 = (self.y) ** 2
        maxy2 = max(y2)
        t = self.x

        alpha = (np.sqrt(self.Etot(
        )) - np.sqrt(self.Emax())) / np.sqrt(self.Etot())
        print alpha
        th = alpha * maxy2

        v = np.nonzero(y2 >= th)[0]
        toa = t[v[0]]
        return toa

    def toa_cum_tm(self):
        """ calculate time of arrival

        """

        y2 = (self.y) ** 2
        t = self.x
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]
        cdf, vary = self.ecdf()

        alpha = np.sqrt(cdf.y[u]) / np.sqrt(cdf.y[-1])
        v = np.nonzero(cdf.y >= alpha * cdf.y[u])[0]
        toa = t[v[0]]
        return toa

    def toa_cum_tmtm(self):
        """ calculate time of arrival

        """

        y2 = (self.y) ** 2
        t = self.x
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]
        cdf, vary = self.ecdf()

        alpha = (np.sqrt(cdf.y[-1]) - np.sqrt(
            cdf.y[u])) / (np.sqrt(cdf.y[-1]) + np.sqrt(cdf.y[u]))
        v = np.nonzero(cdf.y >= alpha * cdf.y[u])[0]
        toa = t[v[0]]
        return toa

    def toa_cum_tmt(self):
        """ calculate time of arrival

        """
        y2 = (self.y) ** 2
        t = self.x
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]
        cdf, vary = self.ecdf()

        alpha = (np.sqrt(cdf.y[-1]) - np.sqrt(cdf.y[u])) / np.sqrt(cdf.y[-1])
        v = np.nonzero(cdf.y >= alpha * cdf.y[u])[0]
        toa = t[v[0]]
        return toa



    def aggcir(self,alphak,tauk):
        """ aggregation of CIR from (alphak,tauk)

        Parameters
        ----------

        alphak : ndarray
            CIR path amplitude
        tauk : ndarray
            CIR delay values

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import numpy as np
            >>> alphak = 10*np.random.rand(7)
            >>> tauk = 100*np.random.rand(7)
            >>> tau = np.arange(0,150,0.1)
            >>> y = np.zeros(len(tau))
            >>> CIR = TUsignal(tau,y)
            >>> CIR.aggcir(alphak,tauk)
            >>> f,a =CIR.plot(typ=['v'])

        """
        shy = np.shape(self.y)
        x = self.x
        eps = (x[1]-x[0])/2
        u = map(lambda t: np.where( (x>t-eps) & (x<=t+eps))[0][0],tauk)
        ynew  = np.zeros(len(x))
        ynew[u] = alphak
        if len(shy)>1:
           self.y = np.vstack((self.y,ynew))
        else:
           self.y = ynew[np.newaxis,:]


    def readcir(self,filename,outdir=[]):
        """ read channel impulse response

        Parameters
        ----------
        filename : string
            long file name if outdir is []
            short file name is outdir is <> []
        outdir : string
            output directory
        """
        if outdir <> []:
            outdir = 'output/'+outdir
            filename = getlong(filename, outdir)

        cir = ios.loadmat(filename)
        self.x = cir['t'].ravel()
        self.y = cir['cir'].ravel()


    def readuwb(self, _filename):
        """ read  Waveform from Matlab file
        """
        outdir = 'output/'+outdir
        filename = getlong(_filename, outdir)
        wfm = ios.loadmat(filename)
        d = wfm['data'][0][0]
        T0 = d.T0[0][0] / 1e-9
        Tres = d.Tres[0][0] / 1e-9
        s = d.WformOut1
        N = len(s)
        self.x = np.linspace(T0, T0 + (N - 1) * Tres, N)
        self.y = s.reshape(len(s))

    def ecdf(self, Tnoise=10, rem_noise=True, in_positivity=True, display=False, normalize=True, delay=0):
        """ calculate energy cumulative density function

        Parameters
        ----------

        Tnoise     :
            Time duration of noise only portion (default=5ns)
        rem_noise  :
            remove noise if True
        in_positivity :
            inforce positivity if True
        normalize  :
            normalize if True (Not implemented)
        display    :
            display ecdf if True
        delay      :
            give a delay for vizualization

        Returns
        -------

        ecdf , vary

        """
        #
        #  ( ) ^2
        #
        t = self.x
        y = self.y
        te = self.dx()
        y2 = y ** 2
        #
        f1 = np.cumsum(y2) * te
        # retrieve the noise only portion at the beginning of TUsignal
        #
        Nnoise = int(np.ceil(Tnoise / te))
        tn = t[0:Nnoise]
        fn = f1[0:Nnoise]
        stdy = np.std(y[0:Nnoise])
        vary = stdy * stdy
        y = t * vary
        #
        # y : linear interpolation of noise ecdf  (over whole time base)
        #
        #(ar,br)= polyfit(tn,fn,1)
        #print ar
        #y  = polyval([ar,br],t)
        if rem_noise:
            f = f1 - y
        else:
            f = f1

        #
        # inforce positivity
        #
        if in_positivity:
            pdf = np.diff(f)
            u = np.nonzero(pdf < 0)[0]
            pdf[u] = 0
            ecdf = np.cumsum(pdf)
        else:
            ecdf = f
        #
        # Normalization step
        #
        E = ecdf[-1]
        #print E

        if normalize:
            ecdf = ecdf / E
        #
        # Resizing
        #
        Nt = len(t)
        Necdf = len(ecdf)
        N = min(Nt, Necdf)
        ecdf = TUsignal(t[0:N], ecdf[0:N])
        #
        # Display
        #
        if display:
            plt.subplot(211)
            ecdf.plot()
            if normalize:
                plt.plot(t, 2 * vary * np.sqrt(2 * t) / E, 'r')
                plt.plot(t, -2 * vary * np.sqrt(2 * t) / E, 'r')
            else:
                plt.plot(t, 3 * vary * np.sqrt(2 * t), 'r')
                plt.plot(t, -3 * vary * np.sqrt(2 * t), 'r')
            plt.axvline(x=delay, color='red')
            plt.subplot(212)
            plt.plot(t, y, color='red')
            plt.plot(t, f1, color='black')
            plt.plot(t, f, color='blue')
            plt.show()

        return ecdf, vary

    def tau_moy(self, alpha=0.1, tau0=0):
        """ calculate mean excess delay starting from delay tau_0

        Parameters
        ----------

        alpha : float
        tau0 : float

        """
        t = self.x
        y = self.y

        cdf, vary = self.ecdf()
        pdf = np.diff(cdf.y)

        u = np.nonzero(cdf.y > alpha)[0]
        v = np.nonzero(cdf.y < 1 - alpha)[0]

        t = t[u[0]:v[-1]]
        pdf = pdf[u[0]:v[-1]]

        te = self.dx()
        a = np.sum(t * pdf)
        b = np.sum(pdf)
        taum = a / b

        return(taum)

    def delays(self):
        r""" calculate delay parameters and orthogonality factor from cir

        Returns
        -------

        taum :
            mean excess delay
        delayspread
            rms delay spread
        of  :
            orthogonality factor

        Neelesh Metha, Andreas Molish, Lary Greenstein "Orthogonality Factor in WCDMA Donlinks in Urban Macrocellular
        environments" 

        .. :math:

            \beta0 = 1 \frac{\sum_i=1^L}|\alpha_i|^4}{\left(\sum_i=1^L|\alpha_i|^2)^2}

        """

        self.flatteny(reversible=True)
        y2 = self.yf*self.yf
        y4 = y2*y2
        taum = sum(self.x*y2,axis=0)/sum(y2,axis=0)
        delayspread = np.sqrt(sum((self.x-taum)*(self.x-taum)*y2)/sum(y2,axis=0))
        of = 1 -sum(y4,axis=0)/sum(y2,axis=0)**2
        return taum,delayspread,of



    def tau_rms(self, alpha=0.1, tau0=0):
        r""" calculate root mean square delay spread starting from delay tau_0

        Parameters
        ----------

        alpha : float
        threshold : float
            ( delay interval is defined between :math:`\tau(\alpha)` and :math:`\tau(1 -\alpha)` )
        tau0 : float
            argument for specifying the delay start

        Notes
        -----

        .. math::

            \sqrt{\frac{\int_{\tau(\alpha)}^{\tau(1-\alpha)} (\tau-\tau_m)^{2} PDP(\tau) d\tau} {\int_{\tau(\alpha)}^{\tau(1-\alpha)} PDP(\tau) d\tau}}

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.measures.mesuwb import *
            >>> import matplotlib.pyplot as plt
            >>> M = UWBMeasure(1)
            >>> ch4 = M.tdd.ch4
            >>> f1,a1=ch4.plot(color='k')
            >>> tit0 = plt.title("WHERE1 M1 UWB Channel impulse response")
            >>> f2,a2=ch4.plot(color='k')
            >>> tit1= plt.title("WHERE1 M1 UWB Channel impulse response (Zoom 1)")
            >>> ax1=plt.axis([10,160,-90,-50])
            >>> f3,a3=ch4.plot(color='k')
            >>> tit2 = plt.title("WHERE1 M1 UWB Channel impulse response (Zoom 2)")
            >>> ax2=plt.axis([20,120,-80,-50])
            >>> plt.show()
            >>> tau_moy = ch4.tau_moy()
            >>> print "tau_moy: %2.2f ns" % tau_moy
            tau_moy: 38.09 ns
            >>> tau_rms = ch4.tau_rms()
            >>> print "tau_rms: %2.2f ns" % tau_rms
            tau_rms: 13.79 ns



        See Also
        --------

        TUsignal.ecdf
        TUsignal.tau_moy

        """

        t = self.x
        y = self.y
        cdf, vary = self.ecdf()
        pdp = np.diff(cdf.y)
        taum = self.tau_moy(tau0)

        u = np.nonzero(cdf.y > alpha)[0]
        v = np.nonzero(cdf.y < 1 - alpha)[0]

        t = t[u[0]:v[-1]]
        pdp = pdp[u[0]:v[-1]]
        te = self.dx()
        b = sum(pdp)
        m = sum(pdp * (t - taum) * (t - taum))
        taurms = np.sqrt(m / b)

        return(taurms)


class TUDsignal(TUsignal):
    """ Uniform signal in Time domain with delay

    Attributes
    ----------

    x   : ndarray
    y   : ndarray
    taud : ndarray
        direct delay
    taue : ndarray
        excess delay

    """
    def __init__(self,x=np.array([]),y=np.array([]),taud=np.array([]),taue=np.array([])):
        super(TUDsignal,self).__init__(x,y)
        #TUsignal.__init__(self, x, y)
        self.taud = taud
        self.taue = taue

    def __repr__(self):
        s = TUsignal.__repr__(self)
        return(s)

    def fig(self, N):
        """ plot a figure of the N first signals

        Parameters
        ----------

        N : int
            number of y signal to plot

        """
        x = self.x
        min = self.y.min()
        max = self.y.max()
        ec = max - min
        ecmax = ec.max()
        sh = np.shape(self.y)
        Nmax = sh[0]
        N1 = int(minimum(N, Nmax))
        y1 = self.y[0, :] + (N1 - 1) * ecmax
        yN1 = self.y[N1 - 1, :]
        for k in range(N):
            gk = str(N) + str(1) + str(k)
            plt.subplot(gk)
            plot(x, yN1[k, :])

           #r.plot(x, yN1, main='Ray response', xlab='Time (ns)', ylab='y', type='l', col='black' ,frame='False',  ylim=r.range(y1,yN1) )
           #for i in range(N1-1):
           #    yi = self.y[i+1,:] + (N1-i)*ecmax
           #    r.lines(x,yi,col='black')


class FBsignal(Bsignal):
    """
    FBsignal : Base signal in Frequency domain

    plot     : plot modulus and phase
    plotri   : plot real part and imaginary part
    plotdB   : plot modulus in dB
    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(FBsignal,self).__init__(x,y,label)
        #Bsignal.__init__(self, x, y)
        self.label[-1]='Frequency (GHz)'

    def __repr__(self):
        s = Bsignal.__repr__(self)
        return(s)

    def plotri(self, nb=-1):
        """ plot real and imaginary part

        Parameters
        ----------
        nb : int

        """
        ndim = self.y.ndim
        if ndim > 1:
            # plot all
            if (nb == -1):
                nl = len(self.y)
                for k in range(nl):
                    plt.subplot(211)
                    plt.plot(self.x, np.real(self.y[k]))
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('Real part')
                    plt.subplot(212)
                    plt.plot(self.x, np.imag(self.y[k]))
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('Imaginary part')
            # plot nb only
            else:
                plt.subplot(211)
                plt.plot(self.x, np.real(self.y[nb]))
                plt.xlabel('Frequency (GHz)')
                plt.ylabel('Real part')
                plt.subplot(212)
                plt.plot(self.x, np.imag(self.y[nb]))
                plt.xlabel('Frequency (GHz)')
                plt.ylabel('Imaginary part')
        else:
            plt.subplot(211)
            plt.stem(self.x, np.real(self.y))
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('Real part')
            plt.subplot(212)
            #plot(self.x,np.unwrap(np.angle(self.y)))
            plt.stem(self.x, np.imag(self.y))
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('Imaginary part')

    def plot(self, **kwargs):
        """ plot FBsignal

        Parameters
        ----------

        phase : boolean
            default True
        dB : boolean
            default True
        iy : index of y value to be displayed
            default [0]  only the first is displayed
        typ : string
            ['l10','l20','d','r','du','ru']
        xlabels
        ylabels

        Examples
        --------

        >>> from pylayers.signal.bsignal import *
        >>> from numpy import *
        >>> from scipy import *
        >>> S = FBsignal()
        >>> S.x = arange(100)
        >>> S.y = cos(2*pi*S.x)+1j*sin(3*pi*S.x+pi/3)
        >>> fig,ax = S.plot()
        >>> plt.show()

        See Also
        --------

        Bsignal.plot
        pylayers.util.plotutil.mulcplot

        """

        if 'typ' not in kwargs:
            kwargs['typ'] = ['l20']
            kwargs['xlabels'] = ['Frequency (GHz)']

        fig,ax = Bsignal.plot(self,**kwargs)

        return fig,ax

    def plotdB(self, mask=False, n=2, phase=True):
        """ usage : plotdB()

        Parameters
        ----------
        mask :
        n :
        phase :

        """
        ndim = self.y.ndim
        if ndim > 1:
            nl = len(self.y)
            for k in range(nl):
                if phase:
                    plt.subplot(211)
                    plt.plot(self.x, 10 * n * np.log10(
                        abs(self.y[k]) + 1.0e-15))
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('Modulus (dB)')
                    plt.subplot(212)
                    plt.plot(self.x, np.unwrap(np.angle(self.y[k])))
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('Phase')
                else:
                    plt.plot(self.x, 10 * n * np.log10(
                        abs(self.y[k]) + 1.0e-15))
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('Modulus (dB)')

        else:
            if mask:
                n = 1
                xfcc = np.array([0., 0.96, 0.96, 0.96, 0.96, 1.61, 1.61, 1.61, 1.61,
                                 1.99, 1.99, 1.99, 1.99, 3.1, 3.1, 3.1, 3.1, 10.6, 10.6, 10.6, 10.6, 15])
                xcept = np.array([0., 1.61, 1.61, 1.61, 1.61, 3.8, 3.8, 3.8, 3.8, 6.,
                                  6, 6, 6, 8.5, 8.5, 8.5, 8.5, 10.6, 10.6, 10.6, 10.6, 15])
                ycept = np.array([-90, -90, -90, -85, -85, -85, -85, -70, -70, -70, -70,
                                  -41.3, -41.3, -41.3, -41.3, -65, -65, -65, -65, -85, -85, -85])
                yfcc = np.array([-41.3, -41.3, -41.3, -75.3, -75.3, -75.3, -75.3, -53.3, -53.3,
                                 -53.3, -53.3, -51.3, -51.3, -51.3, -51.3, -41.3, -41.3, -41.3, -41.3, -51.3, -51.3, -51.3])
                xnoise = np.array([0., 15])
                ynoise = np.array([-114., -114])
                plt.step(xfcc, yfcc, 'b', linewidth=2)
                plt.step(xcept, ycept, 'r', linewidth=2)
                plt.step(xnoise, ynoise, 'k', linewidth=2)
                phase = False
                plt.axis([0, 15, -120, -40])
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('PSD (dBm/MHz)')
            if phase:
                plt.subplot(211)
            plt.plot(self.x, 10 * n * np.log10(abs(self.y) + 1.0e-15),linewidth=0.3)
            if phase:
                plt.subplot(212)
                plt.plot(self.x, np.unwrap(np.angle(self.y)))
                plt.xlabel('Frequency (GHz)')
                plt.ylabel('Phase')

    def stem(self, color='b-'):
        """ stem plot

        Parameters
        ----------
        color : string

        """

        ndim = self.y.ndim
        if ndim > 1:
            nl = len(self.y)
            for k in range(nl):

                plt.subplot(211)
                plt.stem(self.x, np.real(self.y[k]), color)
                plt.xlabel('Frequency (GHz)')
                plt.ylabel('real part)')
                plt.subplot(212)
                plt.stem(self.x, np.imag(self.y[k]), color)
                plt.xlabel('Frequency (GHz)')
                plt.ylabel('imaginary part)')
        else:

            plt.subplot(211)
            plt.stem(self.x, np.real(self.y), color)
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('real part)')
            plt.subplot(212)
            plt.stem(self.x, np.imag(self.y), color)
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('imaginary part)')


class FUsignal(FBsignal, Usignal):
    """
    FUsignal : Uniform signal in Frequency Domain

    Attributes
    ----------

    x  : nd.array((1xN))
    y  : Complex nd.array((M x N )


    Methods
    -------

    symH     : force Hermitian symetry --> FHsignal
    symHz    : force Hermitian symetry with zero padding --> FHsignal
    align    : align two FUsignal on a same frequency base
    enthrsh  : Energy thresholding thresh = 99.99 %
    dBthrsh  : dB thresholding thresh  = 40dB
    ift      : Inverse Fourier transform
    resample : resampling with a new base
    newdf    : resampling with a new df
    zp       : zero padding until len(x) = N
    plotri   : plot real part and imaginary part
    plotdB   : plot modulus in dB
    get      : get k th ray
    tap      : calculates channel taps

    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(FUsignal,self).__init__(x,y,label)
        self.isFriis = False
        #FBsignal.__init__(self, x, y)

    def __repr__(self):
        s = FBsignal.__repr__(self)
        return(s)

    def __add__(self, u):
        L = self.align(u)
        u1 = L[0]
        u2 = L[1]
        U = FUsignal(u1.x, u1.y + u2.y)
        #U =FUsignal(self.x,self.y+u.y)
        return(U)

    def __sub__(self, u):
        L = self.align(u)
        u1 = L[0]
        u2 = L[1]
        U = FUsignal(u1.x, u1.y - u2.y)
        #U =FUsignal(self.x,self.y-u.y)
        return(U)

    def __mul__(self, u):
        L = self.align(u)
        u1 = L[0]
        u2 = L[1]
        #
        # Normalisation 19/05/2009
        Df = u1.x[1] - u1.x[0]
        U = FUsignal(u1.x, u1.y * u2.y)
        #rint shape(self.x)
        #rint shape(self.y)
        #rint shape(u.y)
        # =FUsignal(self.x,self.y*u.y)
        return(U)

    def __div__(self, u):
        L = self.align(u)
        u1 = L[0]
        u2 = L[1]
        U = FUsignal(u1.x, u1.y / u2.y)
        return(U)



    def window(self, win='hamming'):
        """ windowing of FU signal

        Parameters
        ----------
        win : string
            window type ('hamming','blackman','hanning')

        Examples
        --------

        .. plot::
            :include-source:

            >>> import numpy as np
            >>> import matplotlib.pyplot as plt
            >>> from pylayers.signal.bsignal import *
            >>> x = np.arange(2,8,0.1)
            >>> y = np.ones(len(x))
            >>> U = FUsignal(x,y)
            >>> fi = plt.figure()
            >>> fig,ax = U.plot()
            >>> U.window('hamming')
            >>> fig,ax = U.plot()




        """
        Nx = len(self.x)
        shy = np.shape(self.y)
        if len(shy) > 1:
            if shy[0] == Nx:
                Ny = shy[1]
                ind = 0
            else:
                Ny = shy[0]
                ind = 1

        if win == 'hamming':
            w = np.hamming(Nx)
        if win == 'blackman':
            w = np.blackman(Nx)

        if len(shy) > 1:
            if ind == 0:
                w = np.outer(w, np.ones(Ny))
            else:
                w = np.outer(np.ones(Ny), w)

        self.y = self.y * w

    def applyFriis(self):
        r""" apply Friis factor

        The Friis factor is multiplied to y

        .. math::
            y := \frac{c}{4\pif} y

        boolean `isFriis` is set to `True`

        """
        if not self.isFriis:
            factor = -1j*0.3/(4*np.pi*self.x)
            self.y = self.y*factor[np.newaxis,:]
            self.isFriis = True


    def get(self, k):
        """
            get the kh signal

            Parameters
            ----------
            k : indes to get

            Return
            ------
            G : FUsignal

        """
        G = FUsignal()
        G.x = self.x
        G.y = self.y[k, :]
        return(G)

    def info(self):
        """
        Display Information on the FUsignal
        """
        N = len(self.x)
        sh = np.shape(self.y)
        fmin = self.x[0]
        fmax = self.x[-1]
        df = self.x[1] - self.x[0]
        T = 1.0 / df
        print 'FUsignal'
        print '--------'
        print 'N Freq    ', N
        print 'shape(y)  ', sh
        print 'Fmin (GHz) : ', fmin
        print 'Fmax (GHz) : ', fmax
        print 'Duration (ns) :', T
        print 'Frequency sampling step : ', df

    def energy(self,axis=0,Friis=False,mode='mean'):
        r""" calculate energy along given axis

        Parameters
        ----------

        axis : (default 0)
        Friis : boolean
        mode : string
            mean | center | integ | first | last

        Examples
        --------

        >>> S = FUsignal()
        >>> En1 = S.energy()
        >>> assert((En1>0.99) & (En1<1.01))

        Notes
        -----

        If Friis is true energy is multiplied by

        :math:`\frac{-j*c}{4 pi fGHz}`

        axis = 0 is ray axis

        if mode == 'mean'

        :math:`E=\frac{1}{K} \sum_k |y_k|^2`

        if mode == 'integ'

        :math:`E=\delta_x \sum_k |y_k|^2`

        if mode == 'center'

        :math:`E= |y_{K/2}|^2`

        """

        H = self.y

        if Friis:
            factor = -1.j*0.3/(4*np.pi*self.x)
            H = H*factor[np.newaxis,:]

        MH2 = abs(H * np.conjugate(H))

        if mode=='mean':
            EMH2  = MH2.mean(axis=axis)

        if mode=='integ':
            EMH2  = MH2.sum(axis=axis)*(self.x[1]-self.x[0])

        if mode=='center':
            EMH2  = MH2[:,len(self.x)/2]

        if mode=='first':
            EMH2  = MH2[:,0]

        if mode=='last':
            EMH2  = MH2[:,-1]

        return(EMH2)


    def enthrsh(self, thresh=99.99):
        """ Energy thresholding of an FUsignal

        Parameters
        ----------

        thresh : float
            threshold in percentage (default 99.99)

        Returns
        -------

        EMH2 : cumul energie H

        """
        H = self.y
        MH2 = abs(H * np.conjugate(H))
        EMH2 = MH2.sum(axis=1)
        EMH2max = EMH2.max()

        ind1 = EMH2.argsort()
        ind1rev = ind1[::-1, ]
        EMH2sorted = EMH2[ind1rev]

        EMH2cum = EMH2sorted.cumsum()
        EMH2cumnor = EMH2cum * 100 / EMH2cum[-1]
        ind2 = np.nonzero(EMH2cumnor < thresh)[0]
        indices = ind1rev[ind2]

        self.indices = indices
        self.y = H[indices, :]
        self.taud = self.taud[indices]

        return indices

    def dBthrsh(self, threshdB=40):
        """
        dBthrsh : dB thresholding of an FUsignal
        """
        H = self.y
        MH2 = abs(H * np.conjugate(H))
        EMH2 = MH2.sum(axis=1)
        EMH2dB = 10 * np.log10(EMH2)
        print EMH2dB
        EMH2dBmax = EMH2dB.max()

        ind1 = EMH2.argsort()
        print ind1

        EMH2dBsorted = EMH2dB[ind1]
        ind2 = np.nonzero(EMH2dBsorted > (EMH2dBmax - threshdB))[0]
        indices = ind1[ind2]
        print indices
        return indices

    def zp(self, N):
        """ zero padding until length N

        Parameters
        ----------
        N : int

        Returns
        -------
        FUsignal

        """
        N0 = len(self.x)
        sh = np.shape(self.y)
        ndim = np.shape(sh)[0]
        if (ndim > 1):
            N1 = sh[0]
            N2 = sh[1]
        dx = self.dx()
        if N0 < N:
            M = N - N0
            ex = np.linspace(self.x[-1] + dx, self.x[-1] + M * dx, M)
            if (ndim > 1):
                ey = np.zeros((N1, M))
            else:
                ey = np.zeros(M)
            x_new = np.hstack((self.x, ex))
            y_new = np.hstack((self.y, ey))
        U = FUsignal(x_new, y_new)
        return(U)

    def newdf(self, df):
        """ resample the FUsignal using phase and module interpolation

        Parameters
        ----------
            df
        Returns
        -------
            U : FUsignal

        """
        module = abs(self.y)
        df0 = self.x[1] - self.x[0]
        argu = np.unwrap(np.arctan2(np.imag(self.y), np.real(self.y)))
        fm = interp.interp1d(self.x, module)
        fa = interp.interp1d(self.x, argu)
        f_new = np.arange(self.x[0], self.x[-1], df)
        mod_new = fm(f_new)
        arg_new = fa(f_new)
        y_new = mod_new * (np.cos(arg_new) + 1j * np.sin(arg_new)
                           ) * np.exp(-1j * f_new * 2 * np.pi * 100)
        U = FUsignal(f_new, y_new)
        return(U)

    def dftresamp(self, df_new):
        """ non finished 

        Parameters
        ----------

        df_new :

        """
        FH = self.symH(0)
        fh = FH.ifft()


    def resample(self, x_new, kind='linear'):
        """ resample


        This function resamples the Usignal with a new x base
        which needs to be strictly included in the original x
        base of the Usignal.
        x is a 1D array
        y is a 2D array

        Parameters
        ----------

        x_new :
        kind : string

        See Also
        --------

        interp.interp1d, splrep, splev
        """

        if kind == 'linear':
            module = abs(self.y)
            argu = np.unwrap(np.arctan2(np.imag(self.y), np.real(self.y)))
            fm = interp.interp1d(self.x, module)
            fa = interp.interp1d(self.x, argu)

            mod_new = fm(x_new)
            arg_new = fa(x_new)
            y_new = mod_new * (np.cos(arg_new) + 1j * np.sin(arg_new))

        if kind == 'spline':
            module = abs(self.y)
            argu = unwrap(np.arctan2(np.imag(y), np.real(y)))
            coefm = splrep(self.x, module, s=0)
            coefa = splrep(self.x, argu, s=0)

            mod_new = splev(x_new, coefm, der=0)
            arg_new = splev(x_new, coefa, der=0)
            y_new = mod_new * (np.cos(arg_new) + 1j * np.sin(arg_new))

        U = FUsignal(x_new, y_new)
        return(U)

    def symH(self, parity):
        """ enforce Hermitian symetry

        Parameters
        ----------

        parity : integer
            0 even 1 odd

        Returns
        -------

        V  : FHsignal

        """
        f = self.x
        U = self.y
        N = len(f)
        ndim = U.ndim
        ze_x = np.array([0])
        if ndim > 1:
            ze_y = np.zeros([ndim, 1])
        else:
            ze_y = np.array([0])
        if parity == 0:
            if ndim > 1:
                Up = np.concatenate((ze_y, U, np.flipud(
                    np.conjugate(U[:, 0:-1]))), 1)
            else:
                Up = np.concatenate((ze_y, U, np.flipud(
                    np.conjugate(U[0:-1]))), 0)
            fp = np.concatenate((ze_x, f, f[0:-1] + f[-1]))
        else:
            Up = np.concatenate((ze_y, U, np.flipud(np.conjugate(U))), 1)
            fp = np.concatenate((ze_x, f, f + f[-1]))
        V = FHsignal(fp, Up)
        return V

    def symHz(self,Nz,scale='extract'):
        r""" Force Hermitian symmetry with zero padding

        Parameters
        ----------

        Nz :  int
            Number of zero above f[-1]

        Returns
        -------
        SH : FHsignal

        Warnings
        --------

        The signal is rescaled in order to conserve energy

        Let denotes the FUsignal as :math:`mathcal{X}_d`
        The Fourier matrix is :math:`\left[ \matrix{\mathbb{1}\\
                                                    \mathbb{W}_d\\
                                                    \mathbb{W}_u}\right]`


        See Also
        --------

        FHSignal

        Examples
        --------

        >>> N     = 11
        >>> x     = np.linspace(2,11,N)
        >>> y1    = sp.randn(N)+1j*sp.rand(N)
        >>> y2    = sp.randn(N)+1j*sp.rand(N)
        >>> y1[0] = 0
        >>> y2[0] = 0
        >>> y     = np.vstack((y1,y2))
        >>> S     = FUsignal(x,y)
        >>> S.symHz(10)
        FHsignal :  (47,)  (2, 47)

        """
        f = self.x
        df = self.dx()
        U = self.y[:]
        N = len(f)
        Nl = np.int(np.ceil(f[0] / df))

        ndim = U.ndim
        #
        # nline first dimension of self.y
        #
        nline = np.shape(U)[0]
        if ndim > 1:
            zl   = np.zeros([nline, Nl])
            zlm1 = np.zeros([nline, Nl-1])
        else:
            zl = np.zeros(Nl)
            zlm1 = np.zeros(Nl-1)

        #pdb.set_trace()
        if  Nz > 0:
            if ndim > 1:
                zh = np.zeros([nline, Nz])
                UZ = np.concatenate((U, zh), 1)
                UZF = np.fliplr(np.conjugate(UZ))
            else:
                zh = np.zeros(Nz)
                UZ = np.concatenate((U, zh), 1)
                UZF = np.flipud(np.conjugate(UZ))
        #
        # frequency base is divided into two parts
        #
        # fl [ 0 , (Nl-1) *df ]
        # fh [ f[-1]+df , f[-1] + Nz *df ]
        #
        fl = np.linspace(0, (Nl-1) * df, Nl)
        fh = np.linspace(f[-1] + df, f[-1] + Nz * df, Nz)
        fz = np.concatenate((f, fh), 0)
        #Up = np.concatenate((zl, UZ, UZF, zl), 1)
        #fp = np.concatenate((fl, fz, fz + fz[-1], fl + 2 * fz[-1]), 0)
        Up = np.concatenate((zl, UZ, UZF,zlm1), 1)
        fp = np.concatenate((fl, fz, fz + fz[-1], fl[0:-1] + 2 * fz[-1]), 0)

        Nfp = len(fp)
        #if scale == 'extract':
        #    scale = ((Nfp - 1) / (1.0 * N)) / 2.0
        #if scale == 'cir':
        #    scale = ((Nfp - 1) / (1.0 * N)) 
        if scale=='extract':
            scale = ((Nfp-1)/(4*N))
            #scale = Nfp/(N+1)
        if scale=='cir':
            Df = df*Nfp
            #scale = ((Nfp-1)/(2*N))/Df
            scale = ((Nfp-1)/(2*N))/Df

        #self.hermitian=True
        #self.x = fp
        #self.y = Up*scale
        V = FHsignal(fp, Up * scale)

        return V

    def align(self, u2):
        """ align 2 FUsignal

        align <=> intersection
        align : align two FUsignal on a same base
            return a list which contains the two aligned signals

        >>> i1 = EnImpulse()
        >>> i2 = EnImpulse()
        >>> i2.translate(-10)
        >>> L  = i1.align(i2)

        """

        u1 = self
        dx1 = u1.dx()
        dx2 = u2.dx()
        u1_start = u1.x[0]
        u2_start = u2.x[0]
        u1_stop = u1.x[-1]
        u2_stop = u2.x[-1]

        # it starts at the maximum of both signal
        #    stops  at the minimum of both signal
        xstart = max(u1_start, u2_start)
        xstop = min(u1_stop, u2_stop)
        dx = min(dx1, dx2)
        if tstincl(u1.x, u2.x) == 0:
            print 'Warning: Impossible to align the 2 signals'

        if (dx1 <= dx2):

            xnew = ininter(u1.x, xstart, xstop)

            dim1 = u1.len()
            pstart1 = findpos(u1.x, xnew[0])
            pstop1 = pstart1 + findpos(u1.x[pstart1:dim1], xnew[-1])
            u1 = u1.truncate(pstart1, pstop1 + 1)

            u2 = u2.resample(xnew)

        if (dx2 < dx1):

            xnew = ininter(u2.x, xstart, xstop)

            dim2 = u2.len()
            pstart2 = findpos(u2.x, xnew[0])
            pstop2 = pstart2 + findpos(u2.x[pstart2:dim2], xnew[-1])
            u2 = u2.truncate(pstart2, pstop2 + 1)

            u1 = u1.resample(xnew)

        return(u1, u2)

    def ifft(self, Npt=-1):
        r""" Inverse Fourier transform

        Parameters
        ----------

        Npt : int
            default -1 (same number as x) 

        Returns
        -------

        tc : TUsignal

        Notes
        -----

        .. math::
            x = \textrm{linspace}(0,\frac{1}{\delta f},N)

        .. math::
            y = [0,\mathbb{X}_u^T]

        """
        if Npt == -1:
            Npt = len(self.x)

        Y = self.y
        y = fft.ifft(Y, Npt)
        df = self.dx()
        x = np.linspace(0, 1 / df, Npt)
        tc = TUsignal(x, y)
        return tc

    def ift(self, Nz=1, ffts=0):
        """ Inverse Fourier Transform - returns the associated TUsignal

        Parameters
        ----------

        Nz   : Number of zeros (-1) No forcing
        ffts : 0 (no fftshift 1:fftshift)

        Returns
        -------

        uh :

        Examples
        --------

        >>> from pylayers.signal.bsignal import *
        >>> e  = EnImpulse()
        >>> E  = e.fft()
        >>> EU = E.unrex()

        Notes
        ------

        1 - Force Hermitian symmetry --> FHsignal with or without zero padding
        2 - Inverse Fourier transform with or without fftshift

        See Also
        --------

        FHsignal.ifft, FUsignal.symH, FUsignal.symHz

        """
        # enforce Hermitian Symetry
        if (Nz == -1):
            UH = self.symH(1)
        else:
            #UH = self.symHz(Nz,scale='cir')
            UH = self.symHz(Nz,scale='extract')
        # Back in time
        # UH is an FHsignal
        # uh is a TUsignal
        uh = UH.ifft(ffts)
        return(uh)

    def iftshift(self, Nz=1):
        """ Return the associated TUsignal

        Summary
        -------

        apply the inverse fftshift operator to come back in time

        """
        if (Nz == -1):
            UH = self.symH(1)
        else:
            #UH = self.symHz(Nz,scale='cir')
            UH = self.symHz(Nz,scale='extract')

        uh = UH.ifft(ffts=0)
        #uh.y=fliplr(fftshift(fliplr(uh.y)))
        uh.y = flipud(fftshift(flipud(uh.y)))
        return(uh)

    def show(self,**kwargs):
        """ pcolor visualization of Modulus and Phase
        """

        if 'fig' not in kwargs:
            fig = plt.figure()
        else:
            fig = kwargs['fig']
            kwargs.pop('fig')

        ax1 = fig.add_subplot(121)
        fig,ax1 = self.imshow(typ='l20',fig=fig,ax=ax1,**kwargs)
        ax2 = fig.add_subplot(122)
        fig,ax2= self.imshow(typ='d',fig=fig,ax=ax2,**kwargs)

        return fig,[ax1,ax2]
#       def fig(self,N):
#          """
#          """
#          x = self.x
#          min   = abs(self.y).min()
#          max   = abs(self.y).max()
#          ec    = max-min
#          ecmax = ec.max()
#          sh = shape(self.y)
#          Nmax    = sh[0]
#          N1    = int(minimum(N,Nmax))
#          y1    = abs(self.y)[0,:] + (N1-1)*ecmax
#          yN1   = abs(self.y)[N1-1,:]
#          r.par(yaxt="n")
#          r.plot(x, yN1, main='Ray Transfer function', xlab='Freq (GHz)', ylab='', type='l', col='black' ,frame='False',  ylim=r.range(y1,yN1) )
#          for i in range(N1-1):
#              yi = abs(self.y)[i+1,:] + (N1-(i+1))*ecmax
#              r.lines(x,yi,col='black')

    def decimate(self, N=2):
        """ decimate FUsignal by N
        Parameters
        ----------
        N : int
            decimation factor (default 2)

        """
        x = self.x
        y = self.y
        u = np.arange(0, len(x), N)
        xn = x[u]
        if len(np.shape(y)) == 2:
            yn = y[:, u]
        else:
            yn = y[u]

        V = FUsignal(xn, yn)

        return(V)


    def tap(self,**kwargs):
        r""" calculates channel taps

        Parameters
        ----------

        fcGHz : float
            center frequency GHz
        WGHz :  float
            bandwidth GHz
        Ntap :  number of taps
        baseband : boolean
            default : True

        Returns
        -------

        htapi


        Notes
        -----

        [Tse] David Tse, http://www.eecs.berkeley.edu/~dtse/Chapters_PDF/Fundamentals_Wireless_Communication_chapter2.pdf page 26

        """
        defaults = { 'fcGHz':4.5,
                    'WMHz':1,
                    'Ntap':100,
                    'baseband':True}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        fcGHz=kwargs['fcGHz']
        WGHz=kwargs['WMHz']
        Ntap=kwargs['Ntap']
        # yb : tau x f x 1
        if baseband:
            yb = self.y[:,:,np.newaxis]*np.exp(-2 * 1j * np.pi *self.taud[:,np.newaxis,np.newaxis] * fcGHz )
        else:
            yb = self.y[:,:,np.newaxis]
        # l : 1 x 1 x tap
        l  = np.arange(Ntap)[np.newaxis,np.newaxis,:]
        # l : tau x 1 x 1
        tau = self.taud[:,np.newaxis,np.newaxis]
        # S : tau x f x tap (form 2.34 [Tse])
        S   = np.sinc(l-tau*WMHz/1000.)
        # sum over tau : htap : f x tap
        htap = np.sum(yb*S,axis=0)
        # sum over frequency axis : htapi : tap
        # to be considered !! what about the frequency step
        htapi = np.sum(htap,axis=0)

        return htapi



class FUDsignal(FUsignal):
    """
    FUDsignal : Uniform signal in Frequency domain with delays


    Attributes
    ----------

    x    : ndarray 1xN
    y    : ndarray MxN
    taud : direct delay (Time of Flight)
    taue : excess delay

    Methods
    -------

    minphas : force minimal phase    (Not tested)
    totime  : transform to a TUD signal
    iftd    : inverse Fourier transform
    ft1     : construct CIR from ifft(RTF)
    ft2     :

    """
    def __init__(self, x=np.array([]), y=np.array([]), taud=np.array([]),label=[]):
        """ object constructor

        Parameters
        ----------

        x : np.array()
        y : np.array()
        taud : np.array(

        """
        super(FUDsignal,self).__init__(x,y,label)
        self.taud = taud
        self.taue = np.zeros(len(taud))

    def __repr__(self):
        s = FUsignal.__repr__(self)
        return(s)

    def minphas(self):
        """ construct a minimal phase FUsignal

        Notes
        -----

        - Evaluate slope of the phase
        - deduce delay
        - update delay of FUDSignal
        - Compensation of phase slope to obtain minimal phase

        This methods updates the excess delay `taue` member.

        The samplinf frequency step should be

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import numpy as np
            >>> fGHz = np.arange(2,11,0.1)
            >>> tau1 = np.array([1,2,3])[:,np.newaxis]
            >>> y = np.exp(-2*1j*np.pi*fGHz[np.newaxis,:]*tau1)/fGHz[np.newaxis,:]
            >>> H = FUDsignal(x=fGHz,y=y,taud=np.array([15,17,18]))
            >>> f,a = H.plot(typ=['ru'],xlabels=['Frequency GHz'])
            >>> t1 = plt.suptitle('Before minimal phase compensation')
            >>> H.minphas()
            >>> H.taue
            array([ 1.,  2.,  3.])
            >>> f,a = H.plot(typ=['ru'],xlabels=['Frequency GHz'])
            >>> t2 = plt.suptitle('After minimal phase compensation')

        """

        f = self.x
        phase = np.unwrap(np.angle(self.y))
        dphi = phase[:, -1] - phase[:, 0]
        df = self.x[-1] - self.x[0]
        slope = dphi / df
        #if slope >0:
        #   print 'm  inphas Warning : non causal FUSignal'
        #phi0      = +1j*slope*(f[-1]+f[0]/2)
        F, S = np.meshgrid(f, slope)
        #E   = exp(-1j*slope*f+phi0)
        E = np.exp(-1j * S * F)
        self.y = self.y * E
        self.taue = -slope / (2 * np.pi)
        # update total delay
        #self.tau = self.tau+self.taue

    def ifft(self):
        """ inverse Fourier Transform

        Examples
        --------

        >>> from pylayers.simul.link import *
        >>> L = DLink(verbose=False)
        >>> aktk = L.eval()
        >>> L.H.cut()
        >>> T1 = L.H.totime()
        >>> f,a = T1.plot(typ='v')
        >>> L.H.minphas()
        >>> T2 = L.H.totime()
        >>> f,a = T2.plot(typ='v')

        See Also
        --------

        FUsignal.ift


        """
        y = fft.ifft(self.y)
        T = 1/(self.x[1]-self.x[0])
        x = np.linspace(0,T,len(self.x))
        h = TUDsignal(x,y,self.taud,self.taue)
        return(h)


    def totime(self, Nz=1, ffts=0):
        """ transform to TUDsignal

        Parameters
        ----------

            Nz : int
                Number of zeros for zero padding
            ffts : nt
                fftshift indicator (default 0 )

        Examples
        --------

        >>> from pylayers.simul.link import *
        >>> L = DLink(verbose=False)
        >>> aktk = L.eval()
        >>> L.H.cut()
        >>> T1 = L.H.totime()
        >>> f,a = T1.plot(typ='v')
        >>> L.H.minphas()
        >>> T2 = L.H.totime()
        >>> f,a = T2.plot(typ='v')

        See Also
        --------

        FUsignal.ift


        """
        Nray = len(self.taud)
        s = self.ift(Nz, ffts)
        h = TUDsignal(s.x, fft.fftshift(s.y), self.taud,self.taue)
        return(h)


    def iftd(self, Nz=1, tstart=-10, tstop=100, ffts=0):
        """ time pasting

        Parameters
        ----------

        Nz : int
            Number of zeros
        tstart : float
        tstop  : float
        ffts   : int
            fftshift indicator


        Returns
        -------

        rf : TUsignal (1,N)


        See Also
        --------

        TUsignal.translate


        Examples
        --------


        """
        tau = self.taud+self.taue
        Nray = len(tau)
        s = self.ift(Nz, ffts)
        x = s.x
        dx = s.dx()
        x_new = np.arange(tstart, tstop, dx)
        yini = np.zeros((Nray, len(x_new)))
        rf = TUsignal(x_new, yini)
        #
        # initializes a void signal
        #
        for i in range(Nray):
            r = TUsignal(x_new, np.zeros(len(x_new)))
            si = TUsignal(x, s.y[i, :])
            si.translate(tau[i])
            r = r + si
            rf.y[i, :] = r.y
        return rf

    def ft1(self, Nz, ffts=0):
        """  construct CIR from ifft(RTF)

        Parameters
        ----------

        Nz   : number of zeros for zero padding
        ffts : fftshift indicator
            0  no fftshift
            1  apply fftshift

        Returns
        -------

        r : TUsignal


        See Also
        --------

        pylayers.signal.bsignal.


        """
        tau = self.taud+self.taue
        self.s = self.ift(Nz, ffts)
        x = self.s.x
        r = TUsignal(x, np.zeros(len(x)))

        if len(tau) == 1:
            return(self.s)
        else:
            for i in range(len(tau)):
                si = TUsignal(self.s.x, self.s.y[i, :])
                si.translate(tau[i])
                r = r + si
            return r

    def ftau(self, Nz=0, k=0, ffts=0):
        """ time superposition

        Parameters
        ----------

        Nz  : number of zeros for zero padding
        k   : starting index
        ffts = 0  no fftshift
        ffts = 1  apply fftshift

        Returns
        -------
        r : TUsignal
        """
        tau = self.taud + self.taue
        s = self.ift(Nz, ffts)
        x = s.x
        r = TUsignal(x, np.zeros(len(x)))
        si = TUsignal(s.x, s.y[k, :])
        si.translate(tau[k])
        r = r + si
        return r

    def cir(self,fGHzmin=0,fGHzmax=1000):
        """
        """
        u = (self.x>fGHzmin) & (self.y<fGHzmax)
        cir = sum(self.y)


    def plot3d(self,fig=[],ax=[]):
        """ plot in 3D

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import numpy as np
            >>> N = 20
            >>> fGHz = np.arange(1,3,1)
            >>> taud = np.sort(np.random.rand(N))
            >>> alpha = np.random.rand(N,len(fGHz))
            >>> s = FUDsignal(x=fGHz,y=alpha,taud=taud)
            >>> s.plot3d()

        """
        Ntau = np.shape(self.y)[0]
        Nf   = np.shape(self.y)[1]

        if fig==[]:
            fig = plt.figure()

        if ax == []:
            ax  = fig.add_subplot(111, projection = '3d')

        for k,f in enumerate(self.x):
            for i,j in zip(self.taud+self.taue,abs(self.y[:,k])):
                ax.plot([i,i],[f,f],[0,j],color= 'k')

        ax.set_xlabel('Delay (ns)')
        ax.set_xlim3d(0,max(self.taud+self.taue))

        ax.set_ylabel('Frequency (fGHz)')
        ax.set_ylim3d(self.x[0],self.x[-1])

        powermin = abs(self.y).min()
        powermax = abs(self.y).max()
        ax.set_zlabel('Power (linear)')
        ax.set_zlim3d(powermin,powermax)


    def ft2(self, df=0.01):
        """ build channel transfer function (frequency domain)

        Parameters
        ----------

        df : float
            frequency step (default 0.01)

        Notes
        -----

        1. get  fmin and fmax
        2. build a new base with frequency step df
        3. Initialize a FUsignal with the new frequency base
        4. build  matrix tau * f  (Nray x Nf)
        5. buildl matrix E= exp(-2 j pi f tau)
        6. resampling of FUDsignal according to f --> S
        7. apply the element wise product E .* S
        8. add all rays

        """
        fmin = self.x[0]
        fmax = self.x[-1]
        tau = self.taud+self.taue

        f = np.arange(fmin, fmax, df)

        U = FUsignal(f, np.zeros(len(f)))

        TAUF = np.outer(tau, f)
        E = np.exp(-2 * 1j * np.pi * TAUF)

        S = self.resample(f)
        ES = E * S.y
        V = sum(ES, axis=0)
        U.y = V

        return U

class FUDAsignal(FUDsignal):
    """ FUDAsignal : Uniform signal with Delays and Angles


    Attributes
    ----------

    x    : ndarray 1xN
    y    : ndarray MxN
    taud : delay
    tau1 : additional delay

    Methods
    -------

    minphas : force minimal phase    (Not tested)
    totime  : transform to a TUD signal
    iftd    : inverse Fourier transform
    ft1     : construct CIR from ifft(RTF)
    ft2     :

    """
    def __init__(self,
                 x = np.array([]),
                 y = np.array([]),
                 taud = np.array([]),
                 dod = np.array([]),
                 doa = np.array([]),
                 label = []
                 ):
        super(FUDAsignal,self).__init__(x, y,taud,label)
        # FUDsignal.__init__(self, x, y,taud)
        self.dod  = dod
        self.doa  = doa

    def __repr__(self):
        s = FUDsignal.__repr__(self)
        return(s)

    def cut(self,threshold=0.99):
        """ cut the signal at an Energy threshold level

        Parameters
        ----------

        threshold : float
            default 0.99

        """
        self.sort(typ='energy')
        E = self.eprfl()
        cumE = np.cumsum(E)/sum(E)
        v = np.where(cumE<threshold)[0]
        self.taud = self.taud[v]
        self.taue = self.taue[v]
        #self.tau = self.tau[v]
        self.doa = self.doa[v]
        self.dod = self.dod[v]
        self.y = self.y[v,:]

    def sort(self,typ='tau'):
        """ sort FUD signal

        Parameters
        ----------

        typ  : string
            which parameter to sort '
                tau : (default)
                energy

        """

        if typ == 'tau':
            u = np.argsort(self.taud+self.taue)

        if typ == 'energy':
            E = self.eprfl()
            u = np.argsort(E)[::-1]

        self.taud = self.taud[u]
        self.taue = self.taue[u]
        self.doa = self.doa[u]
        self.dod = self.dod[u]
        self.y = self.y[u,:]

    def showtap(self,**kwargs):
        """ show tap

        Parameters
        ----------

        same as tap

        See Also
        --------

        tap

        """

        # f x s  x m x tap

        htap = self.tap(**kwargs)
        # sum over time m
        Et_htap = np.sqrt(np.sum(htap*np.conj(htap),axis=2))/Nm
        # sum over s
        Er_htap = np.sum(htap,axis=1)/Ns
        corrtap = correlate(Er_htap[0,:,0],np.conj(Er_htap[0,:,0]))

    def tap(self,**kwargs):
        """ calculate channel tap

        Parameters
        ----------

        fcGHz : float
            center frequency
        WMHz : floar
            bandwidth
        Ntap : int
            number of taps (related to bandwith)
            as the bandwith increases the potential number of taps increases
        Ns : int
            number of spatial realizations
        Nm : int
            number of time samples
            the channel is sampled along a distance of half a wavelength
        Va : velocity of link termination a
        Vb : velocity of link termination b
        theta_va : float
            theta velocity termination a (in radians)
        phi_va  :
            phi  velocity termination a (in radians)
        theta_vb:
            theta velocity termination b (in radians)
        phi_vb  :
            phi velocity termination b (in radians)


        Examples
        --------

        >>> from pylayers.signal.bsignal import *

        """

        defaults = {'fcGHz':4.5,
                    'WMHz':1,
                    'Ntap':3,
                    'Ns':8,
                    'Nm':10,
                    'Va':1,  #meter/s
                    'Vb':1,  #meter/s
                    'theta_va':0,
                    'phi_va':0,
                    'theta_vb':0,
                    'phi_vb':0 }


        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        fcGHz=kwargs['fcGHz']
        WMHz=kwargs['WMHz']
        Ntap=kwargs['Ntap']
        Ns=kwargs['Ns']
        Nm=kwargs['Nm']
        Va = kwargs['Va']
        Vb = kwargs['Vb']
        # direction of link termination velocity vectors
        theta_va = kwargs['theta_va']
        theta_vb = kwargs['theta_vb']
        phi_va = kwargs['phi_va']
        phi_vb = kwargs['phi_vb']

        Nf = len(self.x)

        mmax = 0.3*WMHz*1e6/(2*fcGHz*(Va+Vb))
        lam = 0.3/fcGHz
        lamo2 = lam/2.
        fmaHz = (Va/0.3)*fcGHz
        fmbHz = (Vb/0.3)*fcGHz
        # Coherence Time
        Tca = 9/(14*np.pi*fmaHz)
        Tcb = 9/(14*np.pi*fmbHz)
        Tc  = 9/(14*np.pi*(fmaHz+fmbHz))

        # DoD DoA

        theta_a = self.dod[:,0]
        phi_a = self.dod[:,1]
        theta_b = self.doa[:,0]
        phi_b = self.doa[:,1]

        # 3 x r
        ska = np.array([np.cos(theta_a)*np.cos(phi_a),np.cos(theta_a)*np.sin(phi_a),np.sin(theta_a)])
        skb = np.array([np.cos(theta_b)*np.cos(phi_b),np.cos(theta_b)*np.sin(phi_b),np.sin(theta_b)])

        # Monte Carlo for spatial realization
        # s x m x tap
        ua0 = (np.cos(theta_va)+1)/2
        va0 =  phi_va/(2*np.pi)
        ub0 = (np.cos(theta_vb)+1)/2
        vb0 =  phi_vb/(2*np.pi)
        # standard deviation of  velocity vector orientation is inversely
        # proportional to velocity magnitude
        ua = (((1/(Va+0.1))*np.random.rand(Ns)+ua0)%1)[:,np.newaxis,np.newaxis]
        va = (((1/(Va+0.1))*np.random.rand(Ns)+va0)%1)[:,np.newaxis,np.newaxis]
        ub = (((1/(Vb+0.1))*np.random.rand(Ns)+ub0)%1)[:,np.newaxis,np.newaxis]
        vb = (((1/(Vb+0.1))*np.random.rand(Ns)+vb0)%1)[:,np.newaxis,np.newaxis]

        # uniform sampling over the sphere
        tha = np.arccos(2*va-1)
        pha = 2*np.pi*ua
        thb = np.arccos(2*vb-1)
        phb = 2*np.pi*ub

        vax = np.cos(tha)*np.cos(pha)
        vay = np.cos(tha)*np.sin(pha)
        vaz = np.sin(tha)*np.cos(pha*0)

        vaxy = np.concatenate([vax[np.newaxis,np.newaxis,np.newaxis,...],vay[np.newaxis,np.newaxis,np.newaxis,...]])
        va = np.concatenate([vaxy,vaz[np.newaxis,np.newaxis,np.newaxis,...]])

        vbx = np.cos(thb)*np.cos(phb)
        vby = np.cos(thb)*np.sin(phb)
        vbz = np.sin(thb)*np.cos(phb*0)

        vbxy = np.concatenate([vbx[np.newaxis,np.newaxis,np.newaxis,...],vby[np.newaxis,np.newaxis,np.newaxis,...]])

        # 3 x r x f x s x m x tap
        vb = np.concatenate([vbxy,vbz[np.newaxis,np.newaxis,np.newaxis,...]])

        # beta : r x f x s x m x tap
        betaa = np.sum(ska[:,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]*va,axis=0)
        betab = np.sum(skb[:,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]*vb,axis=0)


        # m discrete time axis
        # r x f x s x m x tap
        m = np.linspace(0,mmax,Nm)[np.newaxis,np.newaxis,np.newaxis,:,np.newaxis]
        # r x f x s x m x tap
        l  = np.arange(Ntap)[np.newaxis,np.newaxis,np.newaxis,np.newaxis,:]
        # l : r x f x s x m x tap
        tau = self.taud[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]+ \
              self.taue[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]

        ba  = betaa*Va*m/(0.3*WMHz*1e6)
        bb  = betab*Vb*m/(0.3*WMHz*1e6)
        tau2 = tau + ba + bb
        # S : r x f x s x m x tap (form 2.34 [D. Tse])
        S   = np.sinc(l-tau2*WMHz/1000.)
        # sum over r :  f x s  x m x tap
        htap = np.sum(S*self.y[...,np.newaxis,np.newaxis,np.newaxis]*np.exp(-2*1j*np.pi*fcGHz*tau2),axis=0)

        # f x s  x m x tap
        htap  = htap.reshape(Nf,Ns,Nm,Ntap)
        Et_htap = np.sqrt(np.sum(htap*np.conj(htap),axis=2))/Nm
        Er_htap = np.sum(htap,axis=1)/Ns
        corrtap = correlate(Er_htap[0,:,0],np.conj(Er_htap[0,:,0]))
        return(htap,Et_htap,Er_htap,corrtap)

class FHsignal(FUsignal):
    """
    FHsignal : Hermitian uniform signal in Frequency domain

    ifft  : inverse Fourier transform --> TUsignal
    unrex : unredundant extraction    --> FUsignal

    """
    def __init__(self, x=np.array([]), y=np.array([])):
        FUsignal.__init__(self, x, y)

    def __repr__(self):
        s = FUsignal.__repr__(self)
        return(s)

    def __mul__(self, u):
        x = self.x
        # rescaling 19/05/2009
        #Df= x[1]-x[0]
        U = FHsignal(x, self.y * u.y)
        return(U)

    def ifft(self, ffts=0, tt='centered'):
        """ Inverse Fourier Transform

        Parameters
        ----------

        ffts : int
            0 no fftshift (default)
            1 apply fftshift
        tt : string
            'centered'  default

        Returns
        -------

        a real TUsignal

        Examples
        --------

        >>> e  = EnImpulse(fe=200)
        >>> E  = e.fft()
        >>> ee = E.ifft()
        >>> assert(abs(sum(e.y-ee.y))<1e-13)

        """
        Np = len(self.x)
        df = self.x[1] - self.x[0]
        Df = Np * df
        te = 1.0 / Df
        Tww = 1.0 / df
        if (tt == 'centered'):
            t = np.linspace(-0.5 * Tww + te / 2, 0.5 *
                            Tww + te / 2, Np, endpoint=False)
        else:
            t = np.linspace(0, Tww, Np, endpoint=False)

        y = fft.ifft(self.y)
        s = TUsignal()
        s.x = t
        if (ffts == 0):
            s.y = np.real(y)
        else:
            nd = np.ndim(y)
            #print "Nombre de dimensions : ",nd
            #print "shape (y)  : ",shape(y)
            if (nd > 1):
                s.y = np.real(fft.fftshift(y, axes=[1]))
            else:
                s.y = np.real(fft.fftshift(y))

        s.y = s.y * Df
        return(s)

    def unrex(self):
        r""" extraction of the non redundant part of a real signal
        expressed in the frequency domain

        Examples
        --------

        >>> x = np.arange(0,6,1)
        >>> y = np.arange(0,6,1)
        >>> s = TUsignal(x,y)
        >>> S = s.fft()
        >>> U = S.unrex()

        """
        N = len(self.x)
        if np.mod(N, 2) == 0:
            xu = self.x[1:(N + 2) / 2]
            yu = self.y[1:(N + 2) / 2]
        else:
            xu = self.x[1:(N + 1) / 2]
            yu = self.y[1:(N + 1) / 2]

        O = FUsignal(xu, yu)

        return(O)

#class RefPulse154(TUsignal):
#       """
#       Reference pulse of the IEEE 802.15.4a standard
#       """
#   def __init__(self,x=np.array([]),numchanxxfc=4,band=3,thresh=10,fe=20):


class Noise(TUsignal):
    """ Create noise
    """
    def __init__(self, ti=0,tf = 100, fsGHz = 50, PSDdBmpHz=-174, NF=0, R=50, seed=1):
        """ object constructor

        Parameters
        ----------

        ti        : float
            time start (ns)
        tf        : float
            time stop (ns)
        fsGHZ     : float
            sampling frequency
        PSDdBmpHz : float
            Power Spectral Density Noise (dBm/Hz)
        R         : float
            50 Ohms
        NF        : 0
        seed      : []

        """
        TUsignal.__init__(self)
        self._tsns  = 1./fsGHz
        self._ti = ti
        self._tf = tf
        self._fsGHz = fsGHz
        self._PSDdBmpHz = PSDdBmpHz
        self._NF = NF
        self._R = R
        self._seed=seed
        self.eval()

    @property
    def ti(self):
        return(self._ti)

    @property
    def tf(self):
        return(self._tf)

    @property
    def tsns(self):
        return(self._tsns)

    @property
    def R(self):
        return(self._R)

    @property
    def seed(self):
        return(self._seed)

    @property
    def NF(self):
        return(self._NF)

    @property
    def PSDdBmpHz(self):
        return(self._PSDdBmpHz)

    @property
    def fsGHz(self):
        return(self._fsGHz)

    #-------

    @ti.setter
    def ti(self,ti):
        self._ti = ti
        self.eval()

    @tf.setter
    def tf(self,tf):
        self._tf = tf
        self.eval()

    @tsns.setter
    def tsns(self,tsns):
        self._tsns = tsns
        self._fsGHz = 1./tsns
        self.eval()

    @R.setter
    def R(self,R):
        self._R = R
        self.eval()

    @fsGHz.setter
    def fsGHz(self,fsGHz):
        self._fsGHz=fsGHz
        self._tsns=1./fsGHz
        self.eval()

    @NF.setter
    def NF(self,NF):
        self._NF = NF
        self.eval()

    @seed.setter
    def seed(self,seed):
        self._seed = seed
        np.random.seed(seed)
        self.eval()

    @PSDdBmpHz.setter
    def PSDdBmpHz(self,PSDdBmpHz):
        self._PSDdBmpHz = PSDdBmpHz
        #self.eval()


    def eval(self):
        """ noise evaluation
        """
        p = self._PSDdBmpHz + self._NF
        pmW = 10 ** (p / 10.)  # DSP : dBm/Hz -> mW/Hz
        pW = pmW / 1e3         # DSP : mw/Hz  -> W/Hz
        self.PW = pW * (self._fsGHz * 1e9)   # Power : p * Bandwith Hz
        self.vrms = np.sqrt(self._R*self.PW)
        self.x = np.arange(self.ti, self.tf, self.tsns)
        N = len(self.x)
        n = self.vrms * np.random.randn(N)
        self.y   = n
        self.var = np.var(n)
        self.Pr  = self.var/self._R

    def __repr__(self):
        st = ''
        st = st+ 'Sampling frequency : '+ str(self.fsGHz)+' GHz\n'
        st = st+ 'ti  : '+ str(self.ti)+'ns \n'
        st = st+ 'tf  : '+ str(self.tf)+'ns \n'
        st = st+ 'ts  : '+ str(self.tsns)+'ns \n'
        st = st + '-------------\n'
        st = st+ 'DSP : ' + str(self.PSDdBmpHz)+ ' dBm/Hz\n'
        st = st+ 'NF : ' + str(self.NF)+ ' dB\n'
        st = st+ 'Vrms : '+ str(self.vrms)+ ' Volts\n'
        st = st+ 'Variance : '+ str(self.var)+ ' V^2\n'
        st = st+ 'Power /'+str(self.R)+' Ohms : '+ str(10*np.log10(self.PW)-60)+ ' dBm\n'
        st = st+ 'Power realized /'+str(self.R)+' Ohms : '+ str(10*np.log10(self.Pr)-60)+ ' dBm\n'
        return(st)

    def psd(self,mask=True):
        """
        Parameters
        ----------

        mask : boolean
            True
        """
        w2 = TUsignal.psd(self,periodic=False)
        w2.plotdB(mask=mask)

    def amplify(self, GdB, NF):
        sel

    def fgating(self, fcGHz, BGHz, window='rect'):
        """ apply a frequency gating

        Parameters
        ----------

        fcGHz : float
        BGHz  : float
        window : string
            'rect'

        """
        N = self.fft()
        if len(self.x) % 2 == 0:
            parity = 0
        else:
            parity = 1
        U = N.unrex()
        f = U.x
        f1 = fcGHz - BGHz / 2.
        f2 = fcGHz + BGHz / 2.
        u = np.nonzero((f > f1) & (f < f2))[0]
        gate = np.zeros(len(f))
        if window=='rect':
            gate[u] = np.ones(len(u))

        G = FUsignal(f, gate)
        V = G * U
        NF = V.symH(parity)
        nf = NF.ifft()
        return(nf)
        #fe = 1./(self.x[1]-self.x[0])
        #fN = fe/2
        #print "fN : ",fN
        #wp = (fcGHz-BGHz/2.)/fN
        #ws = (fcGHz+BGHz/2.)/fN
        #print "fN : ",wp
        #print "fN : ",ws
        #o  = self.iirfilter(order=4,wp=wp,ws=ws)
    #    self.x = o.x
    #    self.y = o.y


class EnImpulse(TUsignal):
    """
    Create an Energy normalized Gaussian impulse (Usignal)

    EnImpulse(x,fc,band,thresh,fe)

          fc     (GHz)   (def = 4GHz)
          band   (GHz)   (def = 3GHz)
          thresh (dB)    (def = 10dB)
          fe     (GHz)   (def = 100GHz)


    """
    def __init__(self, x=np.array([]), fc=4, band=3, thresh=10, fe=20):
        TUsignal.__init__(self)
        Tp = (2 / (band * np.pi)) * np.sqrt(abs(thresh) * np.log(10) /20.)
        coeff = np.sqrt(2 * np.sqrt(2)/ (Tp * np.sqrt(np.pi)))

        if len(x) == 0:
            te = 1.0 / fe
            Tww = 10 * Tp
            Ni = round(Tww / (2 * te))
            # Tww/2 multiple de te
            Tww = 2 * te * Ni
            x = np.linspace(-0.5 * Tww, 0.5 * Tww, 2 * Ni + 1)

        y = coeff * np.exp(-(x / Tp) ** 2) * np.cos(2 * np.pi * fc * x)
        self.x = x
        self.y = y
        self.Tp = Tp
        self.fc = fc

    def demo():
        """ small demo in the docsting

        Examples
        --------

        >>> from pylayers.signal.bsignal import *
        >>> ip    = EnImpulse(fc=4,band=3,thresh=10,fe=100)
        >>> Eip1  = ip.energy()
        >>> ESDu  = ip.esd(mode='unilateral')
        >>> ESDb  = ip.esd(mode='bilateral')
        >>> df    = ESDu.dx()
        >>> Eipu  = sum(ESDu.y)*df
        >>> Eipb  = sum(ESDb.y)*df
        >>> erru  = Eip1-Eipu
        >>> errb  = Eip1-Eipb

        """
        pass

class MaskImpulse(TUsignal):
    """
    MaskImpulse : Create an Energy normalized Gaussian impulse (Usignal)

          fc     (GHz)   (def = 4GHz)
          band   (GHz)   (def = 3GHz)
          thresh (dB)    (def = 10dB)
          fe     (GHz)   (def = 100GHz)

    Examples
    --------

        >>> ip    = EnImpulse(fc=4,band=3,thresh=10,fe=20)
        >>> Eip1  = ip.energy()
        >>> ESDip = ip.esd()
        >>> df    = ESDip.dx()
        >>> Eip2  = sum(ESDip.y)*df
        >>> err   = Eip1-Eip2

    """
    def __init__(self, x=np.array([]), fc=4, band=3, thresh=10, Tp=100, Pm=-41.3, R=50, fe=100):
        """ object constructor

        Parameters
        ----------
        fc     : center frequency (GHz)
        band   : bandwidth (GHz)
        Tp     : Pulse repetion rate
        R      : Resistance
        Pm     : PSD max
        fe     : sampling freqeuncy
        thresh : definition of band at Pm - thresh (dB)

        """
        self.fc = fc
        self.band = band
        self.thresh = thresh
        self.Tp = Tp
        self.Pm = Pm
        self.R = R
        self.fe = fe

        Usignal.__init__(self)
        #alpha  = 1./(2*np.sqrt(abs(thresh)*np.log(10)/20))
        alpha = 1. / (2 * np.sqrt(abs(thresh) * np.log(10) / 10))
        tau = 1 / (alpha * band * np.pi * np.sqrt(2))
        A = np.sqrt(2 * R * Tp * 10 ** (Pm / 10)) / (tau * np.sqrt(np.pi))
        if len(x) == 0:
            te = 1.0 / fe
            Tw = 10. / band
            Ni = round(Tw / (2 * te))
            # Tww/2 multiple de te
            Tww = 2 * te * Ni
            x = np.linspace(-0.5 * Tww, 0.5 * Tww, 2 * Ni + 1)

        y = A * np.exp(-(x / tau) ** 2) * np.cos(2 * np.pi * fc * x)
        self.x = x
        self.y = y

    def show(self):
        plt.subplot(211)
        self.plot()
        plt.subplot(212)
        P = self.psd(self.Tp, self.R)
        P.plotdB(mask=True)


def test():
    dx1 = 0.01
    x1 = np.arange(-5, 5, dx1)
    s1 = EnImpulse(x1, 4, 2, 10)

    S1 = s1.fft()

    is1 = S1.ifft()

    S12 = S1 * S1
    s12 = S12.ifft(1)

# Nouveau signal moins echantillonne

    dx2 = 0.01
    x2 = np.arange(-3, 3, dx2)
    s2 = EnImpulse(x2, 4, 2, 10)
    S2 = s2.fft()

    U1 = S1.unrex()
    U2 = S2.unrex()

    U3 = U1 * U2

    H3 = U3.symH(1)
    # Comparer H3 et S12

    s3 = H3.ifft(1)
    plt.figure()
    plt.subplot(221)
    plt.title('s1')
    s1.plot()
    plt.subplot(222)
    s12.plot()
    plt.title('s12 = s1*s1')
    plt.subplot(223)
    s2.plot()
    plt.title('s2')
    plt.subplot(224)
    s3.plot()
    plt.title('s3= s1*s2')

    plt.figure()
    S1.plot()
    plt.title('S1')

    plt.figure()
    S12.plot()
    plt.title('S12 = S1*S1')

    plt.figure()
    S2.plot()
    plt.title('S2')

    plt.figure()
    H3.plot()
    plt.title('H3 ')

    plt.show()


    ##
#
#

if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
    #ip1 = EnImpulse(fc=4.493,band=0.499,thresh=3,fe=40)
    #ip2 = EnImpulse(fc=4.493,band=0.499,thresh=3,fe=40)
    #ip2.translate(1.123)
    #ip3 = EnImpulse(fc=4.493,band=0.499,thresh=3,fe=40)
    #ip3.translate(2.067)
    #ip4 = EnImpulse(fc=4.493,band=0.499,thresh=3,fe=40)
    #ip4.translate(3.45677)
    #p = ip1+ip2+ip3+ip4
#       ip.zlr(-10,10)
#    s=TUsignal()
#    s.load('fitest.mat')
     #tau=s.toa()
     #print tau
#   print "> "
#   print "> generate a normalized energy impulse"
#   print "> "
#   print "i1=EnImpulse()"
#   i1=EnImpulse()
#   print "i1.plot()"
#   i1.plot()
#   print "> "
#   print "> translation of the signal acts only on x base "
#   print "> "
#   figure()
#   i1.translate(3)
#   print"i1.translate(3)"
#   i1.plot()
#   print "> "
#   print "> I1 is the fft of i1 "
#   print "> "
#   I1=i1.fft()
#   print "I1=i1.fft()"
#   print "> "
#   print "> Then we extract the unredundant part of I1 "
#   print "> "
#   U =I1.unrex()
#   print "U =I1.unrex()"
#   V =U.symH(1)
#       show()

#   W =U.symHz(4)

#       x = np.arange(-5,5,0.01)
#       s1 = EnImpulse(x,4,2,10)
#       s1.plot()
#       s2 = EnImpulse(x,5,3,10)
#       s2.translate(5)
#
#       u  = s1.convolve(s2)
#
#       S1 = s1.ft()
#       S2 = s2.ft()
#
#       S1S2  = S1*S2
#       S1S2H = S1S2.symHz(100)
#       s1s2h = S1S2H.ifft(1)
#       v = s1s2h
