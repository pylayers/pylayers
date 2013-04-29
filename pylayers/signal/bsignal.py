#!/usr/bin/python
#-*- coding:Utf-8 -*-
import doctest
import os
import pdb
import numpy as np
import scipy as sp
import scipy.interpolate as interp
import numpy.fft as fft
from copy import *
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from pylayers.util.pyutil import *
import scipy.io as ios
from scipy.signal import cspline1d, cspline1d_eval, iirfilter, iirdesign, lfilter, firwin


class Bsignal(object):
    """ Signal with an embedded time base

    This class gathers a 1D signal and its axis indexation.

    This class offers a transparent back and forth mechanism between time and frequency domain.

    The base is not necessarily uniform

    x can have 1 or two axis

    The first axis of x and y have the same length

    By construction len(y):=len(x), len(x) has priority in case of conflict

    """

    def __init__(self, x=np.array([]), y=np.array([])):
        """
        Parameters
        ----------
        x : ndarray
        y : ndarray

        """
        self.x = x
        self.y = y
        ndim = self.y.ndim
        if ndim > 1:
            shy = np.shape(self.y)
            lx = max(np.shape(self.x))
            if (shy[1] != lx):
                print "Error in Bsignal : Dimension incompatibility "
                print "x : ", lx
                print "y : ", shy

    def __repr__(self):
        return '%s :  %s  %s' % (
                            self.__class__.__name__,
                            str(np.shape(self.x)),
                            str(np.shape(self.y)))




    def save(self, filename):
        """
        Save Bsignal in Matlab File Format

        Parameters
        ----------
            filename

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.signal.bsignal import *
            >>> import matplotlib.pyplot as plt
            >>> e = EnImpulse()
            >>> e.plot()
            >>> e.save('impulse.mat')
            >>> del e
            >>> h = TUsignal()
            >>> h.load('impulse.mat')
            >>> h.plot()
        """

        d = {}
        d['x'] = self.x
        d['y'] = self.y
        ios.savemat(filename, d)

    def load(self, filename):
        """ load a Bsignal saved in a Matlab File

        Parameters
        ----------
        filename : string

        """
        d = ios.loadmat(filename)
        self.x = d['x'][:, 0]
        self.y = d['y'][:, 0]

    def setx(self, x):
        """ setx : set x vector

        Parameters
        ----------

        Notes
        -----
        y is set to the corresponding zero vector
        """
        self.x = x
        Np = len(self.x)
        self.y = np.zeros(Np, dtype=float)

    def sety(self, function):
        """ sety : set y vector

        Parameters
        ----------
        function 

        """
        self.y = function(self.x)

    def stem(self, color='b-'):
        """ stem display

        Parameters
        ----------
        color : string 
            default 'b-'

        """
        ndim = self.y.ndim
        if ndim > 1:
            nl = len(self.y)
            for k in range(nl):
                plt.stem(self.x, self.y[k], color)
        else:
            plt.stem(self.x, self.y, color)

    def step(self, color='b'):
        """ plot steps display

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
    
    def imshow(self,interpolation=None,cmap=plt.cm.BrBG,aspect='auto',dB=False):
        """ imshow of y matrix
        """
        if self.y.ndim>1:
            if not dB:
                vmin = abs(self.y.min())
                vmax = abs(self.y.max())
                vm   = max(vmin,vmax)
                vmin = -vm
                vmax = +vm
                val  = self.y
                cmap = plt.cm.BrBG
            else:
                vmin = -100
                vmax = 10*np.log10(abs(self.y.max()))
                val  = 10*np.log10(abs(self.y)+1e-10)
                cmap = plt.cm.hot
            plt.imshow(val,
                       origin = 'lower',
                       vmin = vmin,
                       vmax = vmax,
                       aspect = aspect,
                       extent = (self.x[0],self.x[-1],0,self.y.shape[0]),
                       interpolation=interpolation,
                      # cmap=plt.cm.PiYG)
                       cmap=cmap)
            plt.colorbar()
            plt.axis('auto')

    def plot(self, 
             iy = 0, 
             col = 'black', 
             vline = np.array([]),
             hline = np.array([]),
             unit1 = 'V',
             unit2 = 'V', 
             xmin = -1e5,
             xmax = 1e5,
             ax =[],
             dB = False, 
             dist = False, 
             display = True,
             logx = False, 
             logy = False,
            ):
        """ plot signal

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
        ndim = self.y.ndim
        conversion = 1.0
        if ((unit1 == 'V') & (unit2 == 'mV')):
            conversion = 1000

        u = np.nonzero((self.x > xmin) & (self.x < xmax))[0]
        #
        if dist:
            x = 0.3 * self.x[u]
        else:
            x = self.x[u]
        #
        # More than 1 y
        #
        if ndim > 1:
            nl = len(self.y)
            if (iy == -1):
                sety = range(nl)
            else:
                sety = [iy]
            for k in sety:
                if dB:
                    y = 20 * np.log10(abs(self.y[k, u] * conversion) + 1e-12)
                else:
                    y = self.y[k, u] * conversion

                if ax != []:
                    if logx & logy:
                        ax.loglog(x, abs(y),
                                  color=col)
                    elif logx:
                        ax.semilogx(x, y,
                                    color=col)
                    elif logy:
                        ax.semilogy(x, abs(y),
                                    color=col)
                    else:
                        ax.plot(x, y, color=col)
                else:
                    if logx & logy:
                        loglog(x, abs(y),
                               color=col)
                    elif logx:
                        semilogx(x, y, color=col)
                    elif logy:
                        semilogy(x, abs(y),
                                 color=col)
                    else:
                        plt.plot(x, y, color=col)
        #
        # Only one y
        #
        else:
            if dB:
                y = 20 * np.log10(abs(self.y[u] * conversion) + 1e-12)
            else:
                y = self.y[u] * conversion
            if ax != []:
                if logx & logy:
                    ax.loglog(x, abs(y), color=col)
                elif logx:
                    ax.semilogx(x, y, color=col)
                elif logy:
                    ax.semilogy(x, abs(y),
                                color=col)
                else:
                    ax.plot(x, y, color=col)
            else:
                if logx & logy:
                    plt.loglog(x, abs(y),
                               color=col)
                elif logx:
                    plt.semilogx(x, y, color=col)
                elif logy:
                    plt.semilogy(x, abs(y),
                                 color=col)
                else:
                    plt.plot(x, y, color=col)
        #
        # Draw vertical and horizontal lines
        #
        tcolor = ['red', 'green', 'green', 'green', 'black', 'black', 'black']
        for i in range(len(vline)):
            if ax != []:
                ax.axvline(vline[i], color=tcolor[i])
            else:
                axvline(vline[i], color=tcolor[i])

        for i in range(len(hline)):
            if ax != []:
                ax.axhline(hline[i] * conversion, color='red')
            else:
                axhline(hline[i] * conversion, color='red')

#    def plotdB(self):
#        """
#        plotdB()    obsolete use plot with dB=True option instead
#        """
#        ndim = self.y.ndim
#        if ndim > 1 :
#            nl = len(self.y)
#            for k in range(nl):
#                plot(self.x,20*np.log10(self.y[k]+1e-12))
#        else:
#            plot(self.x,20*np.log10(abs(self.y)+1e-12))
#            xlabel('Time (ns)')
#            ylabel('dB')
    def flatteny(self,yrange=[],reversible=False):
        """ flatten y array
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
        """ gating beween xmin and xmax

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
            >>> s.plot()
            >>> txt1 = plt.title('before gating')
            >>> plt.show()
            >>> s.gating(-3,4)
            >>> s.plot()
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
        """ returm length of x axis
        """
        return(len(self.x))


class Usignal(Bsignal):
    """ Signal with an embedded uniform Base

        This class inheritate from Bsignal. The only difference
        is that the x base is supposed to be uniform

    """

    def __init__(self, x=np.array([]), y=np.array([])):
        Bsignal.__init__(self, x, y)

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
        if ((t == 'int') | (t == 'float')):
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
        """ set the x array of the Usignal (y=0)

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
        """  get the time step of Usignal.x

        Examples
        ---------

        >>> from pylayers.signal.bsignal import *
        >>> u = Usignal()
        >>> u.setx(0,10,0.1)
        >>> assert(u.dx()==0.1)

        """
        return(self.x[1] - self.x[0])

    def width(self):
        """ get the extension support of the Usignal

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
        ------
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
            U = FUDsignal(x_new, y_new, self.tau0)

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
            >>> i3.plot()
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
            L = Usignal(u1.x, np.vstack(u1.y,u2.y))
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
                indx = plt.find((x >= u2_start) & (x <= u2_stop))
                U2 = Usignal(x, np.zeros((N2,len(x))))
                #pdb.set_trace()
                U2.y[:,indx] = u2.y[:, 0:np.shape(indx)[0]]

            if (b2i & b2f):
            # u1 is included in u2
                U2 = u2
                x = u2.x
                indx = plt.find((x >= u1_start) & (x <= u1_stop))
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

    def energy(self):
        """ calculate the energy of an  Usignal

        Returns
        -------
        energy : float

        """
        energy = self.dx() * sum(self.y * np.conj(self.y))
        return(energy)

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
        xmax : float

        Summary
        --------
        This is a gating between xmin and xmax

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.signal import *
            >>> from matplotlib.pylab import *
            >>> ip = EnImpulse()
            >>> ip.plot()
            >>> ip.zlr(-10,10)
            >>> ip.plot()
            >>> show()

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
    def __init__(self, x=np.array([]), y=np.array([])):
        Bsignal.__init__(self, x, y)

    def __repr__(self):
        s = Bsignal.__repr__(self)
        return(s)

    def plot(self,
             iy=0,
             col='black',
             vline=np.array([]),
             hline=np.array([]),
             showlabel=[True, True],
             unit1 = 'V',
             unit2 = 'V',
             ax=[],
             tmin = -1e5,
             tmax = +1e5,
             dB = False,
             dist=False,
             logx=False,
             logy=False,
             ):
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
        >>> s.plot()
        >>> ti = title('TBsignal : plot')
        >>> show()

        """
        if tmin != -1e5:
            xmin = tmin
        else:
            xmin = -1e5
        if tmax != 1e5:
            xmax = tmax
        else:
            xmax = tmax
        Bsignal.plot(self, iy=iy, col=col, vline=vline, hline=hline, unit1=unit1, unit2=unit2,
                     xmin=xmin, xmax=xmax, ax=ax, dB=dB, dist=dist, logx=logx,
                     logy=logy)
        if showlabel[0]:
            if dist:
                plt.xlabel('distance (m)')
            else:
                plt.xlabel('Time (ns)')
        if showlabel[1]:
            if unit2 == 'mV':
                if dB:
                    plt.ylabel('Voltage (dBmV)')
                else:
                    plt.ylabel('Voltage (mV)')
            if unit2 == 'V':
                if dB:
                    plt.ylabel('Voltage (dBV)')
                else:
                    plt.ylabel('Voltage (V)')

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
            >>> ip.plot()
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
            >>> sb.stem()
            >>> su20.plot(col='k')
            >>> su100.plot(col='r')
            >>> ti = plt.title('b2u : sb(blue) su20(black) su200(red)')
            >>> plt.show()

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
    def __init__(self, x=np.array([]), y=np.array([])):
        Usignal.__init__(self, x, y)

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

            >>> import numpy as np
            >>> from pylayers.signal.bsignal import *
            >>> from matplotlib.pyplot import *
            >>> x  = np.arange(0,1.0,0.1)
            >>> y  = x**3
            >>> su = TUsignal(x,y)
            >>> dsu   = su.diff()
            >>> ddsu  = dsu.diff()
            >>> dddsu = ddsu.diff()
            >>> fi = plt.figure()
            >>> su.plot(col='k')
            >>> dsu.plot(col='g')
            >>> ddsu.plot(col='r')
            >>> dddsu.plot(col='b')
            >>> ti = plt.title('TUsignal : diff')
            >>> plt.show()

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
        FHsignal : Frequency signal with hermitian symmetry

        """
        Np = len(self.x)
        te = self.x[1] - self.x[0]
        fe = 1.0 / te
        f = np.linspace(0, fe, Np, endpoint=False)
        #y=fft(self.y)
        #y     = fftshift(y)
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

        Examples
        --------

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

        Examples
        --------

            >>> from pylayers.signal.bsignal import *

        """
        A  = self.fftsh()
        AU = A.unrex()
        return(AU)

    def psd(self, Tpns=100, R=50):
        """ calculate power spectral density

        Parameters
        ----------
        R    : Resistance (default 50 Ohms)
        Tpns : real 
            PRP (default 100 ns)

        .. note::
            If time is in ns the resulting PSD is expressed in dBm/MHz (~10-9)

        """
        P = self.esd(mode='unilateral')
        P.y = P.y / (R * Tpns)
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

        psd
        zlr
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

#       def waterfall(self,N,typ='l'):
#          """
#          """
#          min   = self.y.min()
#          max   = self.y.max()
#          ec    = max-min
#          ecmax = ec.max()
#          sh = np.shape(self.y)
#          Nmax    = sh[0]
#          N1    = int(minimum(N,Nmax))
#
#          x     = self.x
#          r.postscript('waterfall.eps')
#          y1    = self.y[0,:] + (N1-1)*ecmax
#          yN1   = self.y[N1-1,:]
#          r.plot(x, yN1,  main='Waterfall', xlab='Time (ns)', ylab='y', type=typ, col='black' ,  ylim=r.range(y1,yN1) )
#          for i in range(N1-1):
#              yi = self.y[i+1,:] + (N1-i)*ecmax
#          #       stem(self.x,self.y[i+1,:]+(i+1)*ecmax)
#              if (typ=='o'):
#                  r.points(x,yi,col='black')
#              r.lines(x,yi,col='black')
#          #show()
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
#print shape(t_left),shape(x)
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
            u = plt.find(self.y > thre)
            if nbint(u) < nint:
                thre = thre - step
            else:
                thre = thre + step
                step = step / 2

        w = u[1:] - u[0:-1]
        w0 = plt.find(w != 1)
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
        i0 = plt.find((u.x < dx) & (u.x > -dx))
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
        """
        te = self.dx()
        self.y = np.hstack((self.y, np.zeros(N)))
        aux = np.arange(N) * te
        t1 = self.x[-1]
        self.x = np.hstack((self.x, aux - (aux[0] - t1 - te)))


#-------------------------------------------------------------
# Energy content 
#
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
#       u  = nonzero((tau0 + Tint*(1-sym) > self.x) & (self.x > tau0 - Tint*sym))
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
        v = plt.find(Correlation[len(Sx) + n - 200:] > seuil)
        if len(v) == 0:
            ff = seuil / E0
        else:

            w = v[1:] - v[0:-1]
            w0 = plt.find(w != 1)
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
        v = plt.find(self.y[n:] > seuil)
        if len(v) == 0:
            toa = n * te
        else:
            w = v[1:] - v[0:-1]
            w0 = plt.find(w != 1)
            if len(w0) == 0:
                r = max(self.y[n:][v])
                toa = plt.find(self.y == r) * te

            else:
                vv = v[0:w0[0] + 1]
                r = max(self.y[n:][vv])
                toa = plt.find(self.y == r) * te

        u = np.nonzero((toa + Tint * (1 - sym) > self.x) & (
            self.x > toa - Tint * sym))
        efirst = te * sum(self.y[u] * np.conj(self.y[u]))

        if dB:
            efirst = 10 * np.log10(efirst)

        return(efirst)

    def taumax(self):
        """ calculate taumax

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
        u = plt.find(y2 == maxy2)

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
        u = plt.find(y2 == maxy2)
        tau_Emax = t[u]
        return(tau_Emax)

    def toa_max2(self):
        """ calculate time of arrival max2 method
        """

        THRE = array([])
        V = array([])
        VL = array([])

        M = max(self.y)
        n = plt.find(self.y == M)

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

            u = plt.find(self.y > thre)
            v = nbint(u)
            h = plt.find(u > n)
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
        nmax = plt.find(self.y == Max)
        n = nmax
        step = Max / 1e2
        thre = Max - step

        delta = 100
        d = 0
        nint = 0
        N = np.array([])
        N = np.hstack((N, n))

        while delta > 4 * Max / 1e2:

            u = plt.find(self.y > thre)
            hr = plt.find(u > n)
            g = delete(u, hr)

            if nmax >= 6000:
            #set the fenetre=6000*0.005=30ns
                hl = plt.find(g < nmax - 6000)
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
        nmax = plt.find(self.y == Max)
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

            u = plt.find(self.y > thre)
            hr = plt.find(u > n)
            g = delete(u, hr)

            if nmax >= 6000:
            #set the fenetre=6000*0.005=30ns
                hl = plt.find(g < nmax - 6000)
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
            u = plt.find(self.y > thre)
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
        v = plt.find(y2 >= th)
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
        v = plt.find(cdf.y >= th)
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

        v = plt.find(y2 >= th)
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

        v = plt.find(y2 >= th)
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

        v = plt.find(y2 >= th)
        toa = t[v[0]]
        return toa

    def toa_cum_tm(self):
        """ calculate time of arrival 

        """

        y2 = (self.y) ** 2
        t = self.x
        maxy2 = max(y2)
        u = plt.find(y2 == maxy2)
        cdf, vary = self.ecdf()

        alpha = np.sqrt(cdf.y[u]) / np.sqrt(cdf.y[-1])
        v = plt.find(cdf.y >= alpha * cdf.y[u])
        toa = t[v[0]]
        return toa

    def toa_cum_tmtm(self):
        """ calculate time of arrival 

        """

        y2 = (self.y) ** 2
        t = self.x
        maxy2 = max(y2)
        u = plt.find(y2 == maxy2)
        cdf, vary = self.ecdf()

        alpha = (np.sqrt(cdf.y[-1]) - np.sqrt(
            cdf.y[u])) / (np.sqrt(cdf.y[-1]) + np.sqrt(cdf.y[u]))
        v = plt.find(cdf.y >= alpha * cdf.y[u])
        toa = t[v[0]]
        return toa

    def toa_cum_tmt(self):
        """ calculate time of arrival 

        """
        y2 = (self.y) ** 2
        t = self.x
        maxy2 = max(y2)
        u = plt.find(y2 == maxy2)
        cdf, vary = self.ecdf()

        alpha = (np.sqrt(cdf.y[-1]) - np.sqrt(cdf.y[u])) / np.sqrt(cdf.y[-1])
        v = plt.find(cdf.y >= alpha * cdf.y[u])
        toa = t[v[0]]
        return toa

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
        f1 = cumsum(y2) * te
        # retrieve the noise only portion at the beginning of TUsignal
        #
        Nnoise = int(np.ceil(Tnoise / te))
        tn = t[0:Nnoise]
        fn = f1[0:Nnoise]
        stdy = std(y[0:Nnoise])
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
            pdf = diff(f)
            u = plt.find(pdf < 0)
            pdf[u] = 0
            ecdf = cumsum(pdf)
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
        pdf = diff(cdf.y)

        u = plt.find(cdf.y > alpha)
        v = plt.find(cdf.y < 1 - alpha)

        t = t[u[0]:v[-1]]
        pdf = pdf[u[0]:v[-1]]

        te = self.dx()
        a = sum(t * pdf)
        b = sum(pdf)
        taum = a / b

        return(taum)

    def delays(self):
        """ calculate delay parameters and othogonality factor from cir
        Returns
        -------
        taum :
            mean excess delay
        delayspread
            rms delay spread
        of  :
            orthogonality factor

        Orthogonality Factor in WCDMA Donlinks in Urban Macrocellular
        environments : Neelesh Metha, Andreas Molish, Lary Greenstein

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
        """ calculate root mean square delay spread starting from delay tau_0

        Parameters
        ----------
        alpha : 
        tau0 :

        """
        t = self.x
        y = self.y
        cdf, vary = self.ecdf()
        pdf = diff(cdf.y)
        taum = self.tau_moy(tau0)

        u = plt.find(cdf.y > alpha)
        v = plt.find(cdf.y < 1 - alpha)

        t = t[u[0]:v[-1]]
        pdf = pdf[u[0]:v[-1]]
        te = self.dx()
        b = sum(pdf)
        m = sum(pdf * (t - taum) * (t - taum))
        taurms = np.sqrt(m / b)

        return(taurms)


class TUDsignal(TUsignal):
    """ Uniform signal in Time domain with delay

    Attributes
    ----------
    x   : ndarray
    y   : ndarray
    tau : float or ndarray


    .. todo::
        build a tunfold time restitution of delays  

    """
    def __init__(self, x=np.array([]), y=np.array([]), tau=np.array([])):
        TUsignal.__init__(self, x, y)
        self.tau = tau

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
    def __init__(self, x=np.array([]), y=np.array([])):
        Bsignal.__init__(self, x, y)
    
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

    def plot(self, phase=True, dB=True,
             iy=np.array([0]),
             ax=[],fig=[],name=''):
        """ plot

        Parameters
        ----------
        phase : boolean 
            default True 
        dB : boolean
            default True 
        iy : index of y value to be displayed
            default [0]  only first the line is displayed

        Examples  
        --------

        >>> from pylayers.signal.bsignal import *
        >>> from numpy import *
        >>> from scipy import *
        >>> S = FBsignal()
        >>> S.x = arange(100)
        >>> S.y = cos(2*pi*S.x)+1j*sin(3*pi*S.x+pi/3)
        >>> S.plot()
        >>> plt.show()

        """
        
        if phase :
            nrow = 2
        else :
            nrow = 1


        if fig==[]:
            fig,axs=plt.subplots(nrows=nrow,ncols=1,sharex=True,num=name)
        elif ax== []:
            axs=[]
            axs.append(fig.add_subplot(2,1,1))
            axs.append(fig.add_subplot(2,1,2))
#            ff,axs=plt.subplots(nrows=nrow,ncols=1,sharex=True,num=name)

        ndim = self.y.ndim
        if ndim > 1:
            for k in iy:
                if phase:
                    if dB:
                        axs[0].plot(self.x, 20 *
                                 np.log10(abs(self.y[k])))
                    else:
                        axs[0].plot(self.x,
                                 abs(self.y[k]))
                    axs[0].set_ylabel('Modulus')
                    axs[1].plot(self.x,np.unwrap(np.angle(self.y[k])))
                    axs[1].set_xlabel('Frequency (GHz)')
                    axs[1].set_ylabel('Phase (rad)')
                else:
                    if dB:
                        axs.plot(self.x, 20 * np.log10(abs(self.y[k])))
                    else:
                        axs.plot(self.x, abs(self.y[k]))
                    axs.set_xlabel('Frequency (GHz)')
                    axs.set_ylabel('Modulus')
        else:
            if phase:
                if dB:
                    axs[0].plot(self.x, 20 * np.log10(abs(self.y)))
                else:
                    axs[0].plot(self.x, abs(self.y))
                axs[0].set_ylabel('Modulus')
                #plot(self.x,np.unwrap(angle(self.y)))
                axs[1].plot(self.x,np.unwrap(np.angle(self.y)))
                axs[1].set_xlabel('Frequency (GHz)')
                axs[1].set_ylabel('Phase (rad)')
            else:
                axs.plot(self.x, abs(self.y))
                axs.set_xlabel('Frequency (GHz)')
                axs.set_ylabel('Modulus')
        
        return (fig,axs)

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
            plt.plot(self.x, 10 * n * np.log10(abs(self.y) + 1.0e-15))
            if phase:
                plt.subplot(212)
                plt.plot(self.x, np.unwrap(np.angle(self.y)))
                plt.xlabel('Frequency (GHz)')
                plt.ylabel('Phase')

    def stem(self, color='b-'):
        """ stem(self,color='b- ')
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
                ylabel('imaginary part)')
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
    FUsignal : Uniform signal in Frequency domain

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
    plot     : plot modulus and phase
    plotri   : plot real part and imaginary part
    plotdB   : plot modulus in dB
    get      : get k th ray

    """
    def __init__(self, x=np.array([]), y=np.array([])):
        FBsignal.__init__(self, x, y)

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
            >>> U.plot()
            >>> U.window('hamming')
            >>> U.plot()




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

    def energy(self, axis=0):
        """ calculate energy along given axis

        Parameters
        ----------
        axis : (default 0)

        Examples
        --------
            >>> e   = EnImpulse()
            >>> En1 = e.energy()
            >>> E   = e.esd()
            >>> En2 = E.energy()
            >>> assert((En1>0.99) & (En1<1.01))

        """
        H = self.y
        MH2 = abs(H * np.conjugate(H))
        EMH2 = MH2.sum(axis=axis)

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
        ind2 = plt.find(EMH2cumnor < thresh)
        indices = ind1rev[ind2]

        self.indices = indices
        self.y = H[indices, :]
        self.tau0 = self.tau0[indices]

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
        ind2 = plt.find(EMH2dBsorted > (EMH2dBmax - threshdB))
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
        """

        """
        FH = self.symH(0)
        fh = FH.ifft()

    def resample(self, x_new, kind='linear'):
        """
        resample(self,x_new,kind='linear'):

        This function resample the Usignal with a new x base
        which needs to be strictly included in the original x
        base of the Usignal.
        x is a 1D array
        y is a 2D array
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
        FHsignal
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

    def symHz(self, Nz):
        """ Force Hermitian symmetry with zero padding

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
        >>> SH    = S.symHz(10)

        """
        f = self.x
        df = self.dx()
        U = self.y
        N = len(f)
        Nl = np.ceil(f[0] / df)

        ndim = U.ndim
        #
        # nline first dimension of self.y
        #
        nline = np.shape(U)[0]
        if ndim > 1:
            zl = np.zeros([nline, Nl])
        else:
            zl = np.zeros(Nl)

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
        fl = np.linspace(0, (Nl - 1) * df, Nl)
        fh = np.linspace(f[-1] + df, f[-1] + Nz * df, Nz)
        fz = np.concatenate((f, fh), 0)
        Up = np.concatenate((zl, UZ, UZF, zl), 1)
        fp = np.concatenate((fl, fz, fz + fz[-1], fl + 2 * fz[-1]), 0)

        Nfp = len(fp)
        scale = ((Nfp - 1) / (1.0 * N)) / 2.0

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

    def ifft(self, Npt=0):
        """
        """
        if Npt == 0:
            Npt = len(self.x)
        Y = self.y
        y = ifft(Y, Npt)
        df = self.dx()
        x = np.linspace(0, 1 / df, Npt)
        tc = TUsignal(x, y)
        return tc

    def ift(self, Nz=1, ffts=0):
        """
            return the associated TUsignal

        Algorithm
        ---------
            1) Force Hermitian symmetry --> FHsignal
                with or without zero padding
            2) Inverse Fourier transform
                with or without fftshift

        Parameters
        ----------
            Nz   : Number of zeros (-1) No forcing
            ffts : 0 (no fftshift 1:fftshift)

        >>> e  = EnImpulse()
        >>> E  = e.fft()
        >>> EU = E.unrex()
        """
        # forcage de symetrie Hermitienne
        if (Nz == -1):
            UH = self.symH(1)
        else:
            UH = self.symHz(Nz)
        # retour en temps
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
            UH = self.symHz(Nz)

        uh = UH.ifft(ffts=0)
        #uh.y=fliplr(fftshift(fliplr(uh.y)))
        uh.y = flipud(fftshift(flipud(uh.y)))
        return(uh)

    def show(self):
        """ pcolor visualization of Modulus and Phase
        """
        N = np.shape(self.y)[0]
        tt = np.arange(N)
        plt.subplot(121)
        plt.pcolor(self.x, tt, np.abs(self.y))
        plt.title('modulus')
        plt.colorbar()
        plt.subplot(122)
        plt.pcolor(self.x, tt, np.angle(self.y))
        plt.title('Phase (rd)')
        plt. colorbar()
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



class FUDsignal(FUsignal):
    """
    FUDsignal : Uniform signal in Frequency domain with delays


    Attributes
    ----------
        x    : ndarray 1xN
        y    : ndarray MxN
        tau0 : delay
        tau1 : additional delay

    Methods
    -------

    minphas : force minimal phase    (Not tested)
    totud   : transform to a TUD signal
    iftd    : inverse Fourier transform
    ft1     : construct CIR from ifft(RTF)
    ft2     :
    """
    def __init__(self, x=np.array([]), y=np.array([]), tau0=np.array([])):
        FUsignal.__init__(self, x, y)
        self.tau0 = tau0
        self.tau1 = 0.0

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

        """
        f = self.x
        phase = np.unwrap(np.angle(self.y))
        dphi = phase[:, -1] - phase[:, 0]
        df = self.x[-1] - self.x[0]
        slope = dphi / df
        #if slope >0:
        #   print 'minphas Warning : non causal FUSignal'
        #phi0      = +1j*slope*(f[-1]+f[0]/2)
        F, S = np.meshgrid(f, slope)
        #E   = exp(-1j*slope*f+phi0)
        E = np.exp(-1j * S * F)
        self.y = self.y * E
        self.tau1 = -slope / (2 * np.pi)

    def totud(self, Nz=1, ffts=0):
        """ transform to TUDsignal

        Parameters
        ----------
            Nz     : int
                Number of zeros for zero padding
            ffts   : nt
                fftshift indicator (default 0 )
        """
        tau = self.tau0 + self.tau1
        Nray = len(tau)
        s = self.ift(Nz, ffts)
        tud = TUDsignal(s.x, s.y, tau)
        return(tud)

    def iftd(self, Nz=1, tstart=-10, tstop=100, ffts=0):
        """ time pasting

        Parameters
        ----------

        Nz : int
        tstart : float
        tstop  : float
        ffts   : int
            fftshift indicator

        """
        tau = self.tau0
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

        """
        tau = self.tau0
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
        tau = self.tau0
        s = self.ift(Nz, ffts)
        x = s.x
        r = TUsignal(x, np.zeros(len(x)))
        si = TUsignal(s.x, s.y[k, :])
        si.translate(tau[k])
        r = r + si
        return r

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
        tau = self.tau0

        f = np.arange(fmin, fmax, df)

        U = FUsignal(f, np.zeros(len(f)))

        TAUF = np.outer(tau, f)
        E = np.exp(-2 * 1j * np.pi * TAUF)

        S = self.resample(f)
        ES = E * S.y
        V = sum(ES, axis=0)
        U.y = V

        return U

class FUDAsignal(FUsignal):
    """
    FUDAsignal : Uniform signal in Frequency domain with delays and angles


    Attributes
    ----------
        x    : ndarray 1xN
        y    : ndarray MxN
        tau0 : delay
        tau1 : additional delay

    Methods
    -------

    minphas : force minimal phase    (Not tested)
    totud   : transform to a TUD signal
    iftd    : inverse Fourier transform
    ft1     : construct CIR from ifft(RTF)
    ft2     :
    """
    def __init__(self, 
                 x = np.array([]), 
                 y = np.array([]),
                 tau0 = np.array([]),
                 dod = np.array([]),
                 doa = np.array([])):

        FUsignal.__init__(self, x, y)
        self.tau0 = tau0
        self.dod  = dod
        self.doa  = doa
        self.tau1 = 0.0

    def __repr__(self):
        s = FUDsignal.__repr__(self)
        return(s)

    def minphas(self):
        """ construct a minimal phase FUsignal

        Notes
        -----

        - Evaluate slope of the phase
        - deduce delay
        - update delay of FUDSignal
        - Compensation of phase slope to obtain minimal phase

        """
        f = self.x
        phase = np.unwrap(np.angle(self.y))
        dphi = phase[:, -1] - phase[:, 0]
        df = self.x[-1] - self.x[0]
        slope = dphi / df
        #if slope >0:
        #   print 'minphas Warning : non causal FUSignal'
        #phi0      = +1j*slope*(f[-1]+f[0]/2)
        F, S = np.meshgrid(f, slope)
        #E   = exp(-1j*slope*f+phi0)
        E = np.exp(-1j * S * F)
        self.y = self.y * E
        self.tau1 = -slope / (2 * np.pi)

    def totud(self, Nz=1, ffts=0):
        """ transform to TUDsignal

        Parameters
        ----------
            Nz     : int
                Number of zeros for zero padding
            ffts   : nt
                fftshift indicator (default 0 )
        """
        tau = self.tau0 + self.tau1
        Nray = len(tau)
        s = self.ift(Nz, ffts)
        tud = TUDsignal(s.x, s.y, tau)
        return(tud)

    def iftd(self, Nz=1, tstart=-10, tstop=100, ffts=0):
        """ time pasting

        Parameters
        ----------

        Nz : int
        tstart : float
        tstop  : float
        ffts   : int
            fftshift indicator

        """
        tau = self.tau0
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

        """
        tau = self.tau0
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
        tau = self.tau0
        s = self.ift(Nz, ffts)
        x = s.x
        r = TUsignal(x, np.zeros(len(x)))
        si = TUsignal(s.x, s.y[k, :])
        si.translate(tau[k])
        r = r + si
        return r

    def ft2(self, df=0.01):
        """ build channel transfer function (frequency domain)

        Parameters
        ----------
        df : float 
            frequency step (dafault 0.01)
        
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
        tau = self.tau0

        f = np.arange(fmin, fmax, df)

        U = FUsignal(f, np.zeros(len(f)))

        TAUF = np.outer(tau, f)
        E = np.exp(-2 * 1j * np.pi * TAUF)

        S = self.resample(f)
        ES = E * S.y
        V = sum(ES, axis=0)
        U.y = V

        return U

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
        """
        Inverse Fourier Transform
        Parameters
        ----------
            ffts = 0 : no fftshift (default)
            ffts = 1 : apply fftshift
            tt = 'centered'  default

        Returns
        -------

            return a real TUsignal

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
        """
        extraction of the non redundant part of a real signal
        expressed in the frequency domain

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
    """
    Create noise
    """
    def __init__(self, Tobs=100, fe=50, DSPdBmpHz=-174, NF=0, R=50, seed=[]):
        """
        Parameters
        ----------
        Tobs      : Time duration
        fe        : sampling frequency
        DSPdBmpHz : Power Spectral Density Noise (dBm/Hz)
        R         : 50 Ohms
        NF        : 0
        R         : 50
        seed      : []
        """
        TUsignal.__init__(self)
        P = DSPdBmpHz + NF
        pmW = 10 ** (P / 10.)  # mW/Hz
        pW = pmW / 1e3  # W/Hz
        #PW    = pW*2*(fe*1e9)        #   W
        PW = pW * (fe * 1e9)  # W
        std = np.sqrt(R * PW)
        if seed != []:
            np.random.seed(seed)
        self.x = np.arange(0, Tobs, 1. / fe)
        N = len(self.x)
        n = std * np.random.randn(N)
        self.y = n

    def amplify(self, GdB, NF):
        pass

    def gating(self, fcGHz, BGHz, window='rect'):
        """
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

    Example
    -------

        """
    def __init__(self, x=np.array([]), fc=4, band=3, thresh=10, fe=20):
        TUsignal.__init__(self)
        Tp = (2 / (band * np.pi)) * np.sqrt(abs(thresh) * np.log(10) / 20)
        coeff = np.sqrt(2 * np.sqrt(2) / (Tp * np.sqrt(np.pi)))
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
        """
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

    >>> ip    = EnImpulse(fc=4,band=3,thresh=10,fe=20)
    >>> Eip1  = ip.energy()
    >>> ESDip = ip.esd()
    >>> df    = ESDip.dx()
    >>> Eip2  = sum(ESDip.y)*df
    >>> err   = Eip1-Eip2

    """
    def __init__(self, x=np.array([]), fc=4, band=3, thresh=10, Tp=100, Pm=-41.3, R=50, fe=100):
        """
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
        A = np.sqrt(2 * R * Tp * 10 ** (Pm / 10)) / (tau * np.sqrt(pi))
        if len(x) == 0:
            te = 1.0 / fe
            Tw = 10. / band
            Ni = round(Tw / (2 * te))
            # Tww/2 multiple de te
            Tww = 2 * te * Ni
            x = np.linspace(-0.5 * Tww, 0.5 * Tww, 2 * Ni + 1)

        y = A * np.exp(-(x / tau) ** 2) * np.cos(2 * pi * fc * x)
        self.x = x
        self.y = y

    def show(self):
        """
        """
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
