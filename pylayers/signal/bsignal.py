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
    Bsignal.mean
    Bsignal.min
    Bsignal.max
    Bsignal.len
    Bsignal.append
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

Usignal Class
=============

.. autosummary::
    :toctree: generated/

    Usignal.__init__
    Usignal.__repr__
    Usignal.setx
    Usignal.dx
    Usignal.width
    Usignal.expand
    Usignal.max
    Usignal.min
    Usignal.truncate
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
    TBsignal.b2tu
    TBsignal.b2fu

TUsignal Class
==============

.. autosummary::
    :toctree: generated/

    TUsignal.__init__
    TUsignal.__repr__
    TUsignal.__add__
    TUsignal.diff
    TUsignal.info
    TUsignal.align
    TUsignal.filter
    TUsignal.ftshift
    TUsignal.show
    TUsignal.esd
    TUsignal.shift



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

Input Output Functions
-----------------------

.. autosummary::
    :toctree: generated/

    TUsignal.readcir
    TUsignal.readuwb

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

    x has 1 axis

    x and the last axes of y have the same length

    By construction shape(y)[-1] :=len(x), len(x) takes priority in case of observed conflict

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
        label : list of labels

        """
        self.x = x
        self.y = y
        ndim = self.y.ndim
        if ndim==1:
            self.y=self.y.reshape((1,len(self.y)))
            ndim = 2
        assert((len(x)==self.y.shape[-1])or (self.y.shape[-1]==1) )," Dimension problem"
        # default naming of the axis
        if label==[]:
            self.label=[]
            for k in range(ndim):
                self.label.append('ax'+str(k))
        else:
            self.label=label

        shx = np.shape(self.x)
        shy = np.shape(self.y)
        self.N  = shx[0]
        ly  = shy[-1]
        # multi axes indexation
        self.uax = np.hstack((np.ones(ndim-1),np.r_[self.N]))
        # last dimension of y should be equal to first dimension of x
        #if (ly != self.N) :
        #    print "Error in Bsignal : Dimension incompatibility "
        #    print "x : ", self.N
        #    print "y : ", ly


    def __repr__(self):
        st = '%s :  %s  %s ' % (
                            self.__class__.__name__,
                            str(np.shape(self.x)),
                            str(np.shape(self.y)))

        #for k in range(self.y.ndim):
        #    st = st + '\n' +self.label[k]+ ' : ' + str(self.y.shape[k])

        return(st)

    def mean(self):
        """ mean value of the signal
        """
        S = type(self)()
        S.__dict__=self.__dict__
        S.x = self.x
        S.y = np.mean(self.y,axis=0)[None,:]
        return(S)

    def min(self):
        """ min value of the signal (module if complex)
        """
        S = type(self)()
        S.__dict__=self.__dict__
        S.x = self.x
        if np.iscomplex(self.y).any():
            S.y = np.min(np.abs(self.y),axis=0)[None,:]
        else:
            S.y = np.min(self.y,axis=0)[None,:]
        return(S)

    def max(self):
        """ max value of the signal (module if complex)
        """
        S = type(self)()
        S.__dict__=self.__dict__
        S.x = self.x
        if np.iscomplex(self.y).any():
            S.y = np.max(np.abs(self.y),axis=0)[None,:]
        else:
            S.y = np.max(self.y,axis=0)[None,:]

        return(S)

    def append(self,bs):
        """ append bs to Bsignal
        """
        assert((self.x==bs.x).all())
        self.y = np.vstack((self.y,bs.y))

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
        O.y = O.y[...,u]
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
            >>> e = TUsignal()
            >>> e.EnImpulse(feGHz=100)
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
        if self.y.ndim==1:
            self.y = self.y[None,:]
        self.N = len(self.x)

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
            >>> si.y= np.arange(100)[None,:]
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
                           extent = (self.x[0],self.x[-1],0,self.y.shape[0]),
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

            return fig, ax

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

        Returns
        -------

        fig : figure
        ax : np.array of axes

        See Also
        --------

        pylayers.util.plotutil.mulcplot


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
                  'idx'   :[0,0,0,0,0,0,0]
                 }

        for key, value in defaults.items():
            if key not in kwargs:
                 kwargs[key] = value


        vline = kwargs['vline']
        hline = kwargs['hline']

        idx = kwargs['idx']
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
        if ndim == 4:
            yx = self.y[idx[0],idx[1],idx[2],u]
            fig,ax = mulcplot(self.x[u],yx*conversion,**args)
        if ndim == 3:
            yx = self.y[idx[0],idx[1],u]
            fig,ax = mulcplot(self.x[u],yx*conversion,**args)
        if ndim == 2:
            yx = self.y[idx[0],u]
            fig,ax = mulcplot(self.x[u],yx*conversion,**args)
        if ndim==1:
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

    def __repr__(self):
        s = Bsignal.__repr__(self)
        return(s)

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
        if (len(self.x)>1):
            return(self.x[1] - self.x[0])
        else:
            return(0)

    def resample(self, x_new, kind='linear'):
        """ resample the Usignal with a new x base

        x_new needs to be strictly included in the original x base of the Usignal.

        x is a 1D array
        y is a 2D array

        if y is complex the interpolation is done on module and unwrapped
        phase separately

        Parameters
        ----------

        x_new : ndarray
        kind  : string
            'linear' |'spline'
        """
        if kind == 'linear':
            if np.iscomplex(self.y).any():
                module = abs(self.y)
                argu = np.unwrap(np.arctan2(np.imag(self.y), np.real(self.y)))
                fm = interp.interp1d(self.x, module)
                fa = interp.interp1d(self.x, argu)

                mod_new = fm(x_new)
                arg_new = fa(x_new)
                y_new = mod_new * (np.cos(arg_new) + 1j * np.sin(arg_new))
            else:
                f = interp.interp1d(self.x, self.y)
                y_new = f(x_new)

        if kind == 'spline':
            if np.iscomplex(self.y).any():
                module = abs(self.y)
                argu = unwrap(np.arctan2(np.imag(y), np.real(y)))
                coefm = splrep(self.x, module, s=0)
                coefa = splrep(self.x, argu, s=0)

                mod_new = splev(x_new, coefm, der=0)
                arg_new = splev(x_new, coefa, der=0)
                y_new = mod_new * (np.cos(arg_new) + 1j * np.sin(arg_new))
            else:
                coef = splrep(self.x, self.y, s=0)
                y_new = splev(x_new, coef, der=0)

        U = type(self)(x_new, y_new)
        return(U)

    def alignc(self, u2):
        """ align 2 Usignal

        alignc <=> intersection
        alignc : align two Usignal on a same base
            return a list which contains the two aligned signals

        Returns
        -------

        L : Usignal
            concatenated signal L1.y and L2.y


        """
        u1 = self
        # nothing to align
        if len(u1.x)==len(u2.x):
            if (u1.x==u2.x).all():
                return u1.x,u1.y,u2.y

        sh1 = u1.y.shape
        sh2 = u2.y.shape
        naxis1 = len(sh1)
        naxis2 = len(sh2)
        assert(naxis1==naxis2),"Problem signal haven't the same number of axis"

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

        assert(xstop>=xstart), "error the 2 signal have disjoint support"

        dx = min(dx1, dx2)
        if tstincl(u1.x, u2.x) == 0:
            print 'Warning: impossible to align the 2 signals'

        if (dx1 <= dx2):

            xnew = ininter(u1.x, xstart, xstop)

            dim1 = u1.len()
            #pstart1 = findpos(u1.x, xnew[0])
            pstart1 = np.where(u1.x==xnew[0])[0]
            #pstop1 = pstart1 + findpos(u1.x[pstart1:dim1], xnew[-1])
            pstop1 = pstart1 + np.where(u1.x[pstart1:dim1]==xnew[-1])[0]
            u1 = u1.truncate(pstart1, pstop1 + 1)

            u2 = u2.resample(xnew)

        if (dx2 < dx1):

            xnew = ininter(u2.x, xstart, xstop)

            dim2 = u2.len()
            #pstart2 = findpos(u2.x, xnew[0])
            pstart2 = np.where(u2.x==xnew[0])[0]
            #pstop2 = pstart2 + findpos(u2.x[pstart2:dim2], xnew[-1])
            pstop2 = pstart2 + np.where(u2.x[pstart2:dim2]==xnew[-1])[0]

            u2 = u2.truncate(pstart2, pstop2 + 1)
            u1 = u1.resample(xnew)
        #L   = Usignal()
        #L.x = u1.x
        #L.y = np.vstack((u1.y,u2.y))
        #u1.y=u1.y[...,None]
        #u2.y=u2.y[...,None]
        #pdb.set_trace()
        #L.y = np.concatenate((u1.y,u2.y),axis=naxis1)

        return(u1.x,u1.y,u2.y)

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
            y_new = self.y[..., posmin:posmax]
        else:
            y_new = self.y[posmin:posmax]

        U = type(self)(x_new, y_new)

        return(U)


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
        shy = self.y.shape
        nay = len(shy)
        energy = self.dx() * np.sum(self.y * np.conj(self.y),axis=nay-1)
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
            >>> ip = TUsignal()
            >>> ip.EnImpulse()
            >>> f,a = ip.plot(typ=['v'])
            >>> ip.zlr(-10,10)
            >>> f,a = ip.plot(typ=['v'])

        """
        dx = self.dx()
        cmin = min(self.x)
        cmax = max(self.x)
        shy = self.y.shape
        shyl = shy[0:-1]
        if cmin > xmin:
            xaddl = np.arange(xmin, cmin - dx, dx)
            nshy = list(shyl)+[len(xaddl)]
            yaddl = np.zeros(nshy)
            self.x = np.hstack((xaddl, self.x))
            self.y = np.concatenate((yaddl, self.y))
        else:
            u = np.nonzero(self.x >= xmin)
            self.x = self.x[u[0]]
            self.y = self.y[...,u[0]]

        if cmax < xmax:
            xaddr = np.arange(cmax + dx, xmax, dx)
            nshy = list(shyl)+[len(xaddr)]
            yaddl = np.zeros(nshy)
            self.x = np.hstack((self.x, xaddr))
            self.y = np.concatenate((self.y, yaddr))
        else:
            u = np.nonzero(self.x <= xmax)
            self.x = self.x[u[0]]
            self.y = self.y[...,u[0]]

    def window(self, win='hamming'):
        """ windowing Usignal
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
            >>> U = Usignal(x,y)
            >>> fig,ax = U.plot()
            >>> U.window('hamming')
            >>> fig,ax = U.plot()


        """
        if win=='rect':
            exit
        if win == 'hamming':
            w = np.hamming(len(self.x))[None,:]
        if win == 'blackman':
            w = np.blackman(len(self.x))[None,:]
        if win == 'hanning':
            w = np.hanning(len(self.x))[None,:]
        self.windowed = True
        self.y = self.y * w

class TBsignal(Bsignal):
    """  Based signal in Time domain


    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        Bsignal.__init__(self, x, y)
        self.label = 'time (ns)'

    def __repr__(self):
        s = Bsignal.__repr__(self)
        return(s)

    def plot(self,**kwargs):
        """ plot TBsignal

        Parameters
        ----------

        unit1 : actual unit of data
        unit2 : unit for display
        dist : boolean
        xmin : float
        xmax : float
        logx : boolean
        logy :boolean

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

    def integ(self,Tns,Tsns=50):
        """ integtation of alphak tauk

        used energy detector for IEEE 802.15.6 standard


        """
        u1 = np.where(self.x<(self.x[0]+Tns))[0]
        u2 = np.where((self.x<(self.x[0]+Tsns+Tns)) &
                      (self.x>=(self.x[0]+Tsns)))[0]
        Hp = np.sum((self.y[u1])**2)
        Hi = np.sum((self.y[u2])**2)

        return(Hp,Hi)

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
            >>> ip = TUsignal()
            >>> ip.EnImpulse()
            >>> ip.translate(-10)
            >>> fig,ax=ip.plot(typ=['v'])
            >>> show()


        """
        self.x = self.x + tau

    def b2fud(self, N=300):
        r""" conversion into a FUDsignal


        Parameters
        ----------

        N :

        Notes
        -----

        This method is assuming that each element of TBsignal is a delta function.

        $$ h = \sum__k y \delta(x-x_k)$$

        $$ H = \sum__k y \exp(2j\pi f x_k)$$

        """
        # difference of times 
        dtau = self.x[1:]-self.x[0:-1]
        # determine the minimum value
        mindtau = np.min(dtau)
        # fix maximum frequency as the inverse of the minimum delay between
        # delta functions
        fmax  = 1./mindtau
        # create an uniform frequency base
        f = np.linspace(fmax/(1.0*N),fmax,N)
        z = np.sum(self.y[:,None]*np.exp(-2*1j*f[None,:]*np.pi*self.x[:,None]),axis=0)
        H = FUDsignal(f,z,taud=self.x)
        return H

    def b2tu(self, N):
        """ conversion into a TUsignal

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
            >>> su20 = sb.b2tu(20)
            >>> su100 = sb.b2tu(100)
            >>> fi = plt.figure()
            >>> st = sb.stem()
            >>> fig,ax = su20.plot(color='k')
            >>> fig,ax = su100.plot(color='r')
            >>> ti = plt.title('b2tu : sb(blue) su20(black) su200(red)')

        """

        fi = interp.interp1d(self.x, self.y, kind='linear')
        xn = np.linspace(self.x[0], self.x[-1], N)
        yn = fi(xn)
        U = TUsignal(xn, yn)

        return U


class TUsignal(TBsignal, Usignal):
    """ Uniform signal in time domain

    This class inheritates from TBsignal and Usignal

    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(TUsignal,self).__init__(x,y,label)

    def __add__(self, u):
        t = type(u).__name__
        if ((t == 'int') | (t == 'float')):
            U = type(self)(self.x, self.y + u)
        if type(u)==type(self):
            x,y1,y2 = self.align(u)
            U = type(self)()
            U.x = x
            U.y = y1 + y2
        return(U)

    def __sub__(self, u):
        t = type(u).__name__
        if ((t == 'int') | (t == 'float')):
            U = type(self)(self.x, self.y - u)
        if type(u)==type(self):
            x,y1,y2 = self.align(u)
            U = type(self)()
            U.x = x
            U.y = y1 - y2
        return(U)

    def __mul__(self, u):
        t = type(u).__name__
        if ((t == 'int') | (t == 'float') | (t== 'float64') ):
            U = type(self)()
            U.x = self.x
            U.y = self.y*u
        if issubclass(type(self),type(u)):
            x,y1,y2 = self.align(u)
            U = type(self)()
            U.x = x
            U.y = y1 * y2

        return(U)


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

        >>> from pylayers.signal.bsignal import *
        >>> e = TUsignal()
        >>> e.EnImpulse()
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
            >>> i1 = TUsignal()
            >>> i2 = TUsignal()
            >>> i1.EnImpulse()
            >>> i2.EnImpulse()
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
                U2.y[:,indx] = u2.y[:, 0:np.shape(indx)[0]]

            if (b2i & b2f):
            # u1 is included in u2
                U2 = u2
                x = u2.x
                indx = np.nonzero((x >= u1_start) & (x <= u1_stop))[0]
                U1 = Usignal(x, np.zeros((N1,len(x))))
                U1.y[:,indx] = u1.y

            #L = [U1, U2]
            #L   = Usignal()
            #L.x = U1.x
            #L.y = np.vstack((U1.y,U2.y))
        return U1.x,U1.y,U2.y
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

        O : Output filtered TUsignal


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
        Y  = self.fft()
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


    def EnImpulse(self,**kwargs):
        """
        Create an Energy normalized Gaussian impulse (Usignal)

        Parameters
        ----------

        fcGHz : float
        WGHz : float
        threshdB : float
        feGHz : float

        """
        defaults = {'fcGHz' : 4,
                    'WGHz' : 3,
                    'threshdB' : 10,
                    'feGHz' : 20
                   }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        WGHz = kwargs.pop('WGHz')
        fcGHz = kwargs.pop('fcGHz')
        feGHz = kwargs.pop('feGHz')
        threshdB = kwargs.pop('threshdB')

        #TUsignal.__init__(self)
        Tp = (2 / (WGHz * np.pi)) * np.sqrt(abs(threshdB) * np.log(10) /20.)
        coeff = np.sqrt(2 * np.sqrt(2)/ (Tp * np.sqrt(np.pi)))

        te = 1.0 / feGHz
        Tww = 10 * Tp
        Ni = round(Tww / (2 * te))
        # Tww/2 multiple de te
        Tww = 2 * te * Ni
        x = np.linspace(-0.5 * Tww, 0.5 * Tww, 2 * Ni + 1)

        y = coeff * np.exp(-(x / Tp) ** 2) * np.cos(2 * np.pi * fcGHz * x)
        self.x = x
        self.y = y[None,:]
        self.Tp = Tp
        self.fcGHz = fcGHz

    #    def demo():
    #        """ small demo in the docsting
    #
    #        Examples
    #        --------
    #
    #        >>> from pylayers.signal.bsignal import *
    #        >>> ip    = EnImpulse(fc=4,band=3,thresh=10,fe=100)
    #        >>> Eip1  = ip.energy()
    #        >>> ESDu  = ip.esd(mode='unilateral')
    #        >>> ESDb  = ip.esd(mode='bilateral')
    #        >>> df    = ESDu.dx()
    #        >>> Eipu  = sum(ESDu.y)*df
    #        >>> Eipb  = sum(ESDb.y)*df
    #        >>> erru  = Eip1-Eipu
    #        >>> errb  = Eip1-Eipb
    #
    #        """
    #        pass


    def MaskImpulse(self,**kwargs):
        """
        MaskImpulse : Create an Energy normalized Gaussian impulse (Usignal)

        def __init__(self, x=np.array([]), fc=4, band=3, thresh=10, Tp=100, Pm=-41.3, R=50, fe=100):

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

    #    def show(self):
    #        plt.subplot(211)
    #        self.plot()
    #        plt.subplot(212)
    #        P = self.psd(self.Tp, self.R)
    #        P.plotdB(mask=True)


class FBsignal(Bsignal):
    """
    FBsignal : Base signal in Frequency domain

    plot     : plot modulus and phase
    plotri   : plot real part and imaginary part
    plotdB   : plot modulus in dB
    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        Bsignal.__init__(self, x, y)
        self.label='Frequency (GHz)'

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


class FUsignal(FBsignal,Usignal):
    """
    FUsignal : Uniform signal in Frequency Domain

    Attributes
    ----------

    x  : nd.array((1xN))
    y  : Complex nd.array((... x N )


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
    window   :

    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(FUsignal,self).__init__(x,y,label)
        self.isFriis = False

    def __repr__(self):
        s = FBsignal.__repr__(self)
        return(s)

    def __add__(self, u):
        x,y1,y2 = self.alignc(u)
        U = FUsignal(x, y1 + y2)
        return(U)

    def __sub__(self, u):
        x,y1,y2 = self.alignc(u)
        U = FUsignal(x, y1 - y2)
        return(U)

    def __mul__(self, u):
        x,y1,y2 = self.alignc(u)
        U = FUsignal(x, y1 * y2)
        return(U)

    def __div__(self, u):
        x,y1,y2 = self.alignc(u)
        U = FUsignal(x, y1 / y2)
        return(U)


    def applyFriis(self):
        r""" apply Friis factor

        The Friis factor is multiplied to y

        .. math::
            y := \( \frac{-j c}{4 \pi f} \) y

            x is frequency in GHz

        boolean `isFriis` is set to `True`

        """

        if not self.isFriis:
            factor = -1j*0.3/(4*np.pi*self.x)
            factor = factor.reshape(self.uax)
            self.y = self.y*factor
            self.isFriis = True



    def get(self, k):
        """
        get the kh signal

        Parameters
        ----------
        k : indes to get

        Returns
        -------

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

    def energy(self,axis=1,Friis=False,mode='mean'):
        r""" calculate energy along given axis

        Parameters
        ----------

        axis : (default 0)
        Friis : boolean
        mode : string
            mean | center | integ | first | last

        Examples
        --------

        >>> ip = TUsignal()
        >>> ip.EnImpulse()
        >>> En1 = ip.energy()
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



    def symH(self, parity=0):
        """ enforce Hermitian symetry

        Parameters
        ----------

        parity : integer
            0 even 1 odd

        Returns
        -------

        V  : FHsignal

        """
        #assert self.x[0]!=0
        f = self.x
        U = self.y
        N = len(f)

        ys = U.shape
        ze_x = np.array([0])
        ze_y = np.zeros((ys[0],1))

        if parity == 0:
            Up = np.concatenate((ze_y, U, np.flipud(np.conjugate(U[:, 0:-1]))), 1)
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
        scale : string
            default : 'extract' | 'cir'

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
        >>> FH =  S.symHz(10)

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
        if scale=='extract':
            scale = np.sqrt(Nfp/(2.*N))

        if scale=='cir':
            Df = df*Nfp
            #scale = ((Nfp-1)/(2*N))/Df
            scale = ((Nfp-1)/(2.*N))/Df

        #self.hermitian=True
        #self.x = fp
        #self.y = Up*scale
        V = FHsignal(fp, Up * scale)

        return V


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
        >>> e = TUsignal()
        >>> e.EnImpulse()
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
        if 'vmin' in kwargs: 
            del kwargs['vmin']
        if 'vmax' in kwargs: 
            del kwargs['vmax']
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
        defaults = {'fcGHz':4.5,
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



class FHsignal(FUsignal):
    """
    FHsignal : Hermitian uniform signal in Frequency domain

    Methods
    -------

    ifft  : inverse Fourier transform --> TUsignal
    unrex : unredundant extraction    --> FUsignal

    """
    #def __init__(self, x=np.array([]), y=np.array([]),label=[]):
    #    FUsignal.__init__(self, x, y,label)

    def __repr__(self):
        s = FUsignal.__repr__(self)
        return(s)

    def __mul__(self, u):
        x = self.x
        # rescaling 19/05/2009
        #Df= x[1]-x[0]
        U = FHsignal(x, self.y * u.y)
        return(U)

    def ifft(self, ffts=0, centered=True):
        """ Inverse Fourier Transform

        Parameters
        ----------

        ffts : int
            0 no fftshift (default)
            1 apply fftshift
        centered:  boolean

        Returns
        -------

        a real TUsignal

        Examples
        --------

        >>> e  = TUsignal()
        >>> e.EnImpulse(feGHz=200)
        >>> E  = e.fft()
        >>> ee = E.ifft()
        >>> assert(abs(sum(e.y-ee.y).all())<1e-13)


        """
        Np = len(self.x)
        df = self.x[1] - self.x[0]
        Df = Np * df
        te = 1.0 / Df
        Tww = 1.0 / df
        if centered:
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
        # even case
        if np.mod(N, 2) == 0:
            xu = self.x[1:(N + 2) / 2]
            yu = self.y[:,1:(N + 2) / 2]
        # odd case
        else:
            xu = self.x[1:(N + 1) / 2]
            yu = self.y[:,1:(N + 1) / 2]

        O = FUsignal(x=xu, y=yu)

        return(O)

#class RefPulse154(TUsignal):
#       """
#       Reference pulse of the IEEE 802.15.4a standard
#       """
#   def __init__(self,x=np.array([]),numchanxxfc=4,band=3,thresh=10,fe=20):


class Noise(TUsignal):
    """ Create noise
    """
    def __init__(self,
                 ti=0,
                 tf = 100,
                 fsGHz = 50,
                 PSDdBmpHz = -174,
                 NF = 0,
                 R = 50, seed=1):
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
        NF        : float
            (Default Value : 0)
        seed      : float

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
        r""" noise evaluation

        $$N_0 = 4 \times 10^{-21}$$

        """
        e  = self._PSDdBmpHz + self._NF
        emJ = 10 ** (e / 10.)  # DSP : dBm/Hz -> mW/Hz
        eJ = emJ / 1e3         # DSP : mw/Hz  -> W/Hz
        self.PW = eJ * (self._fsGHz * 1e9)   # Power : p * Bandwith Hz
        self.vrms = np.sqrt(self._R*self.PW)
        self.x = np.arange(self.ti, self.tf, self.tsns)
        N = len(self.x)
        n = self.vrms * np.random.randn(N)
        self.y   = n
        self.var = np.var(n)
        self.Pr  = self.var/self._R
        self.Er  = self.Pr/(self._fsGHz*1e9)

    def __repr__(self):
        st = ''
        st = st+ 'Sampling frequency : '+ str(self.fsGHz)+' GHz\n'
        st = st+ 'ti  : '+ str(self.ti)+'ns \n'
        st = st+ 'tf  : '+ str(self.tf)+'ns \n'
        st = st+ 'ts  : '+ str(self.tsns)+'ns \n'
        st = st+ 'N   : '+ str(len(self.x))+'\n'
        st = st + '-------------\n'
        st = st+ 'DSP : ' + str(self.PSDdBmpHz)+ ' dBm/Hz\n'
        st = st+ '    : ' + str(10**(self.PSDdBmpHz/10.)*1e-3)+ ' Joules\n'
        st = st + '-------------\n'
        st = st+ 'Noise Figure : ' + str(self.NF)+ ' dB\n'
        st = st+ 'Vrms : '+ str(self.vrms)+ ' Volts\n'
        st = st+ 'Variance : '+ str(self.var)+ ' V^2\n'
        st = st+ 'Power (dBm) /'+str(self.R)+' Ohms : '+ str(10*np.log10(self.PW)-60)+ ' dBm\n'
        st = st+ 'Power realized /'+str(self.R)+' Ohms : '+ str(10*np.log10(self.Pr)-60)+ ' dBm\n'
        return(st)

    def ppsd(self,mask=True):
        """ plot Power Spectral Density

        Parameters
        ----------

        mask : boolean
            True
        """
        W = TUsignal.psd(self,periodic=False)
        W.plotdB(mask=mask)

    def amplify(self, GdB, NF):
        sel

    def fgating(self, fcGHz, WGHz, window='rect'):
        """ apply a frequency gating

        Parameters
        ----------

        fcGHz : float
        WGHz  : float
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
        f1 = fcGHz - WGHz / 2.
        f2 = fcGHz + WGHz / 2.
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
        #wp = (fcGHz-WGHz/2.)/fN
        #ws = (fcGHz+WGHz/2.)/fN
        #print "fN : ",wp
        #print "fN : ",ws
        #o  = self.iirfilter(order=4,wp=wp,ws=ws)
    #    self.x = o.x
    #    self.y = o.y



def test():
    dx1 = 0.01
    x1 = np.arange(-5, 5, dx1)
    s1 = TUsignal()
    s1 = s1.EnImpulse(x1, fcGHz=4, WGHz=2, feGHz=10)

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
    pass
    #plt.ion()
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
