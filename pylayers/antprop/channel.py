# -*t coding:Utf-8 -*-
r"""
.. currentmodule:: pylayers.antprop.channel

.. autosummary::
    :toctree: generated

TBchannel class
===============

.. autosummary::
    :toctree: generated/

    TBchannel.tau_Emax
    TBchannel.tau_moy
    TBchannel.tau_rms
    TBchannel.delays
    TBchannel.SalehValenzuela

TUchannel class
===============

.. autosummary::
    :toctree: generated/

    TUchannel.psd
    TUchannel.awgn
    TUchannel.Etau0
    TUchannel.Ewin
    TUchannel.Etot
    TUchannel.Emax


TUDchannel
==========

.. autosummary::
    :toctree: generated/

Tchannel
=========

Members
-------

    filcal : calibration file
    win : string type of window {'rect'| }


.. autosummary::
    :toctree: generated/

    Tchannel.saveh5
    Tchannel.loadh5
    Tchannel.apply
    Tchannel.applywavC
    Tchannel.chantap
    Tchannel.applywavB
    Tchannel.applywavA
    Tchannel.plotd
    Tchannel.plotad
    Tchannel.doadod
    Tchannel.energy
    Tchannel.wavefig
    Tchannel.rayfig
    Tchannel.rssi
    Tchannel.cut
    Tchannel.sort
    Tchannel.showtap
    Tchannel.tap
    Tchannel.minphas
    Tchannel.ifft
    Tchannel.totime
    Tchannel.iftd
    Tchannel.ft1
    Tchannel.ftau
    Tchannel.cir
    Tchannel.plot3d
    Tchannel.ft2
    Tchannel.frombuf
    Tchannel.capacity
    Tchannel.calibrate
    Tchannel.pdp

"""
import doctest
import pdb
import numpy as np
import scipy as sp
import pylab as plt
import struct as stru
import scipy.stats as st
import numpy.fft as fft
import pylayers.util.pyutil as pyu
import pylayers.signal.bsignal as bs
import pylayers.util.geomutil as geu
import pylayers.antprop.antenna as ant
from pylayers.antprop.raysc import GrRay3D
from pylayers.util.project import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
try:
    import h5py
except:
    print 'h5py is not installed: Ctilde(object cannot be saved)'

class TBchannel(bs.TBsignal):
    """ radio channel in non uniform delay domain
    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        #super(TUsignal,self).__init__(x,y,label)
        bs.TBsignal.__init__(self,x,y,label)


    def tau_Emax(self):
        """ calculate the delay of max energy pildeeak

        .. math::
            \max_{\tau} y^{2}(\tau)
        """
        y2 = (self.y) ** 2
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]
        tau_Emax = self.x[u]
        return(tau_Emax)

    def tau_moy(self, alpha=0.1, tau0=0):
        """ calculate mean excess delay starting from delay tau_0

        Parameters
        ----------

        alpha : float
        tau0 : float

        """
        t = self.x
        y = self.y

        #cdf, vary = self.ecdf()

        cdf = np.cumsum(self.y,axis=1)
        cdf = cdf/cdf[:,-1][:,None]

        pdf = np.diff(cdf.y)

        u = np.nonzero(cdf.y > alpha)[0]
        v = np.nonzero(cdf.y < 1 - alpha)[0]

        t = t[u[0]:v[-1]]
        pdf = pdf[u[0]:v[-1]]

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



        See Also
        --------

        TUsignal.ecdf
        TUsignal.tau_moy

        """

        t = self.x
        y = self.y
        #cdf, vary = self.ecdf()
        #pdp = np.diff(cdf.y)

        cdf = np.cumsum(self.y,axis=1)
        cdf = cdf/cdf[:,-1][:,None]

        taum = self.tau_moy(tau0)

        u = np.nonzero(cdf.y > alpha)[0]
        v = np.nonzero(cdf.y < 1 - alpha)[0]

        t = t[u[0]:v[-1]]
        pdp = pdp[u[0]:v[-1]]
        b = sum(pdp)
        m = sum(pdp * (t - taum) * (t - taum))
        taurms = np.sqrt(m / b)

        return(taurms)

    def toFD(self,fGHz=np.linspace(2,5,256)):
        """ Transform to Frequency domain

        Returns
        -------

        H : Tchannel

        """

        z = np.sum(self.y[:,None]*np.exp(-2*1j*fGHz[None,:]*np.pi*self.x[:,None]),axis=0)
        H = Tchannel(x=fGHz,y=z,tau=self.x)
        return H

    def SalehValenzuela(self,**kwargs):
        """ generic Saleh and Valenzuela Model

        Parameters
        ----------

        Lam : clusters Poisson Process parameter (ns)
        lam : rays Poisson Process parameter (ns)
        Gam : clusters exponential decay factor
        gam : rays exponential decay factor
        T   : observation duration

        Examples
        --------

        >>> from pylayers.antprop.channel import *
        >>> C=TBchannel()
        >>> C.SalehValenzuela()
        >>> C.stem()

        """
        defaults = { 'Lam' : .1,
                     'lam' : .5,
                     'Gam' : 30,
                     'gam' : 5 ,
                     'T'   : 100}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        Lam = kwargs['Lam']
        lam = kwargs['lam']
        Gam = kwargs['Gam']
        gam = kwargs['gam']
        T   = kwargs['T']
        Nr  = 1.2*T/Lam
        Nc  = 1.2*T/lam
        e1 = st.expon(1./Lam)
        e2 = st.expon(1./lam)

        # cluster time of arrival
        tc   = np.cumsum(e1.rvs(Nr))
        tc   = tc[np.where(tc<T)]
        Nc   = len(tc)
        tauc = np.kron(tc,np.ones((1,Nr)))[0,:]

        # rays time of arrival
        taur = np.cumsum(e2.rvs((Nr,Nc)),axis=0).ravel()

        # exponential decays of cluster and rays
        etc = np.exp(-tauc/(1.0*Gam))
        etr = np.exp(-taur/(1.0*gam))
        et = etc*etr
        tau = tauc+taur

        # filtering < T and reordering in delay domain
        tau = tau[np.where(tau<T)]
        et = et[np.where(tau<T)]
        u = np.argsort(tau)
        taus = tau[u]
        ets  = et[u]*np.sign(np.random.rand(len(u))-0.5)

        # delays and amplitudes
        self.x = taus
        self.y = ets

class TUchannel(TBchannel,bs.TUsignal):
    """ Uniform channel in delay domain
    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        super(TBchannel,self).__init__(x,y,label)

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
            g = np.delete(u, h)
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
        """ estimate time of arrival (new method)

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
            g = np.delete(u, hr)

            if nmax >= 6000:
            #set the fenetre=6000*0.005=30ns
                hl = np.nonzero(g < nmax - 6000)[0]
                u = np.delete(g, hl)
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
                #print N

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
            g = np.delete(u, hr)

            if nmax >= 6000:
            #set the fenetre=6000*0.005=30ns
                hl = np.nonzero(g < nmax - 6000)[0]
                u = np.delete(g, hl)
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
            number of intervals

        """
        #
        # seek fot the maximum value of the signal
        #
        M = self.y.max()
        step = M / 1e2
    #       plot(self.x,self.y)
        thre = M - step
        while step > M / 1e5:
    #          axhline(y=thre,color='green')
            u = np.where(self.y > thre)[0]
            # nbint : number of contiguous intervals
            if pyu.nbint(u) < nint:
            # down
                thre = thre - step
            else:
            # up + step reduction
                thre = thre + step
                step = step / 2.

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

    def awgn(self,PSDdBmpHz=-174,snr=0,seed=1,typ='psd',R=50):
        """ add a white Gaussian noise

        Parameters
        ----------

        PSDdBmpHz : float
        snr : float
        seed : float
        typ : string
            'psd' | 'snr'
        R : float

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

        sn.y = self.y + n.y[0:len(self.x)]
        sn.x = self.x

        return sn,n

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
        """ Etot  calculate the energy of the signal

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
        ecdf = bs.TUsignal(t[0:N], ecdf[0:N])
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

class TUDchannel(TUchannel):
    """ Uniform channel in Time domain with delay

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
        super(TUDchannel,self).__init__(x,y)
        #TUsignal.__init__(self, x, y)
        self.taud = taud
        self.taue = taue

    def __repr__(self):
        s = TUchannel.__repr__(self)
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

class Tchannel(bs.FUsignal):
    """ Handle the transmission channel

    The transmission channel TChannel is obtained through combination of the propagation
    channel and the antenna transfer functions from both transmitter and receiver.

    Members
    -------

        ray transfer functions  (nray,nfreq)
    dod  :
        direction of depature (rad) [theta_t,phi_t]  nray x 2
    doa  :
        direction of arrival (rad)  [theta_r,phi_r]  nray x 2
    tauk :
        delay ray k in ns

    Methods
    -------

    imshow()
    apply(W)
    applywavB(Wgam)
    applywavC(Wgam)
    chantap(fcGHz,WGHz,Ntap)
    doddoa()
    wavefig(w,Nray)
    rayfig(w,Nray)
    rssi(ufreq)

    See Also
    --------

    pylayers.antprop.Ctilde.prop2tran


    """
    def __init__(self,
                x = np.arange(0,2,1),
                y = np.arange(0,2,1),
                tau  = np.array(([],)),
                dod  = np.array(([[],[]])).T,
                doa  = np.array(([[],[]])).T,
                label = ''):
        """ class constructor

        Parameters
        ----------

        x  :  , nfreq
            frequency GHz
        y  :  nray x nfreq
            path amplitude
        tau   :  1 x nray
            path delay (ns)
        dod   :  direction of departure (nray x 2)
        doa   :  direction of arrival   (nray x 2)

        """
        self.taud = tau
        self.taue = np.zeros(len(tau))
        # FUDsignal.__init__(self, x, y,taud)
        self.dod  = dod
        self.doa  = doa
        # , Nf
        # Nd x Nf x Np x Nu
        self.label = label
        self.win = 'rect'
        self.isFriis = False
        self.windowed = False
        self.calibrated = False
        self.filcal="calibration.mat"
        bs.FUsignal.__init__(self,x=x,y=y,label='Channel')


    def __repr__(self):
        st = 'Tchannel : Ray transfer function (Nray x Nr x Nt x Nf)\n'
        st = st+'-----------------------------------------------------\n'
        st = st + 'freq : '+str(self.x[0])+' '+str(self.x[-1])+' '+str(len(self.x))+"\n"
        st = st + 'shape  : '+str(np.shape(self.y))+"\n"
        st = st + 'tau (min, max) : '+str(min(self.taud))+' '+str(max(self.taud))+"\n"
        st = st + 'dist (min,max) : '+str(min(0.3*self.taud))+' '+str(max(0.3*self.taud))+"\n"
        if self.isFriis:
            st = st + 'Friis factor -j c/(4 pi f) has been applied'

        if self.calibrated:
            st = st+'\n calibrated : Yes\n'
        else:
            st = st+'\n calibrated : No\n'

        if self.windowed:
            st = st+' windowed : Yes\n'
            st = st+self.win+'\n'
        else:
            st = st+' windowed : No\n'

        return(st)

        return(st)


    def saveh5(self,Lfilename,idx,a,b,Ta,Tb):
        """ save Ctilde object in hdf5 format

        Parameters
        ----------

        Lfilename  : string
            Layout filename
        Tilde
            file identifier number
        a : np.ndarray
            postion of point a (transmitter)
        b : np.ndarray
            postion of point b (receiver)
        Ta : np.ndarray
            rotation matrice of antenna a
        Tb : np.ndarray
            rotation matrice of antenna b

        """
        _Lfilename=Lfilename.split('.')[0]
        filename= _Lfilename +'_' + str(idx).zfill(5) + '.h5'
        filenameh5=pyu.getlong(filename,pstruc['DIRH'])

        f=h5py.File(filenameh5,'w')

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f.attrs['a']=a
            f.attrs['b']=b
            f.attrs['Ta']=Ta
            f.attrs['Tb']=Tb
            # keys not saved as attribute of h5py file
            for k,va in self.__dict__.items():
                f.create_dataset(k,shape = np.shape(va),data=va)
            f.close()
        except:
            f.close()
            raise NameError('Channel Tchannel: issue when writting h5py file')

    def loadh5(self,Lfilename,idx, output = True):
        """ load Ctilde object in hdf5 format

        Parameters
        ----------

        Lfilename  : string
            Layout filename
        idx : int
            file identifier number
        output : bool
            return an output precised in return

        Returns
        -------

        if output:
        (a,b,Ta,Tb)

        with
            a = np.ndarray
                postion of point a (transmitter)
            b = np.ndarray
                postion of point b (receiver)
            Ta = np.ndarray
                rotation matrice of antenna a
            Tb = np.ndarray
                rotation matrice of antenna b


        """
        filename = Lfilename.split('.')[0] +'_' + str(idx).zfill(5) + '.h5'
        filenameh5 = pyu.getlong(filename,pstruc['DIRH'])

        f=h5py.File(filenameh5, 'r')
        try:
            # keys not saved as attribute of h5py file
            for k,va in f.items():
                # if k != 'tau1':
                #     setattr(self,str(k),va[:])
                # else :
                setattr(self,str(k),va)

            a = f.attrs['a']
            b = f.attrs['b']
            Ta = f.attrs['Ta']
            Tb = f.attrs['Tb']
            f.close()

            self.__init__(self.x, self.y, self.taud, self.dod, self.doa)

            if output :
                return a,b,Ta,Tb

        except:
            f.close()
            raise NameError('Channel Tchannel: issue when reading h5py file')

    def _saveh5(self,filenameh5,grpname):
        """ save Tchannel object in hdf5 format compliant with Link Class

        Parameters
        ----------

        filenameh5  : str
            file name of h5py file Link format
        grpname  : int
            groupname in filenameh5

        """


        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:

            fh5=h5py.File(filename,'a')
            if not grpname in fh5['H'].keys():
                fh5['H'].create_group(grpname)
            else :
                print 'Warning : H/'+grpname +'already exists in '+filenameh5
            f=fh5['H/'+grpname]

            for k,va in self.__dict__.items():
                #print k,va
                f.create_dataset(k,shape = np.shape(va),data=va)
            fh5.close()
        except:
            fh5.close()
            raise NameError('Channel Tchannel: issue when writting h5py file')

    def _loadh5(self,filenameh5,grpname):
        """ Load H object in hdf5 format compliant with Link Class

        Parameters
        ----------

        filenameh5  : str
            file name of h5py file Link format
        grpname  : int
            groupname in filenameh5

        """
        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        try:
            fh5=h5py.File(filename,'r')

            f = fh5['H/'+grpname]

            # keys not saved as attribute of h5py file
            for k,va in f.items():
                if k !='isFriis':
                    try:
                        setattr(self,str(k),va[:])
                    except:
                        setattr(self,str(k),va)
                else :
                    setattr(self,str(k),va)
            fh5.close()
            self.__init__(self.x, self.y, self.taud, self.dod, self.doa)

        except:
            fh5.close()
            raise NameError('Channel Tchannel: issue when reading h5py file')


    def apply(self, W):
        """ apply FUsignal W to the Tchannel

        Parameters
        ----------

        W :  Bsignal.FUsignal

        It exploits multigrid convolution from Bsignal.


        Returns
        -------

        mport ipdb
        V : FUDAsignal

        Notes
        -----

        Returns :math:`W(f) H_k(f)`

            + W may have a more important number of points and a smaller frequency band.
            + If the frequency band of the waveform exceeds the one of the
            Transmission Channel, a warning is sent.
            + W is a FUsignal whose shape doesn't need to be homogeneous with FUChannel H

        """

        U = W * self 
        V = Tchannel(x= U.x, y = U.y, tau = self.taud, dod = self.dod, doa= self.doa)

        return(V)


    def applywav(self, Wgam):
        """ apply waveform (time domain ) to obtain the
            rays impulses response

            this is the 2015 vectorized method for applying 
            wav on Tchannel

        Parameters
        ----------

        Wgam : waveform

        Returns
        -------

        rir  : array, 
            impulse response for each ray separately
            the size of the array is (nb_rays, support_length)
            support_length is calculated in regard of the 
            delays of the channel

            

        Notes
        ------

            The overall received signal is built in time domain

            Wgam is applied on each Ray Transfer function

        See Also
        --------

        pylayers.signal.channel.rir

        """

        # product in frequency domain between Channel (self) and waveform
        Y = self.apply(Wgam)
        # back in time domain
        rir = Y.rir(Nz=500,ffts=1)
        return rir


    def get_cir(self,Wgam):
        """ get Channel impulse response of the channel 
            for a given waveform

        Parameters
        ----------

        Wgam : waveform

        Returns
        -------

        ri  : TUsignal

            impulse response for each ray separately


        """

        rir = self.applywav(Wgam)
        cir = np.sum(rir.y,axis=0)
        return bs.TUsignal(rir.x, cir)
        

    def applywavC(self, w, dxw):
        """ apply waveform method C
        DEPRECATED

        Parameters
        ----------
        w :
            waveform
        dxw :

        Notes
        -----

        The overall received signal is built in time domain
        w is apply on the overall CIR


        """
        print DeprecationWarning(
            'WARNING : Tchannel.applywavC is going to be replaced by Tchannel.applywav')
        H = self.H
        h = H.ft1(500, 1)
        dxh = h.dx()
        if (abs(dxh - dxw) > 1e-10):
            if (dxh < dxw):
                # reinterpolate w
                f = interp1d(w.x, w.y)
                x_new = arange(w.x[0], w.x[-1], dxh)[0:-1]
                y_new = f(x_new)
                w = bs.TUsignal(x_new, y_new)
            else:
                # reinterpolate h
                f = interp1d(h.x, h.y)
                x_new = arange(h.x[0], h.x[-1], dxw)[0:-1]
                y_new = f(x_new)
                h = bs.TUsignal(x_new, y_new)

        ri = h.convolve(w)
        return(ri)

    def chantap(self,**kwargs):
        """ channel tap

        Parameters
        ----------

        fcGHz :
            WGHz  :
        Ntap  : int

        """

        defaults = {'fcGHz':4.5,
                    'WGHz':1,
                    'Ntap':100}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        h = bs.Tchannel(x=self.x, y=self.y, tau=self.taud)
        htap = h.chantap(**kwargs)
        return htap

    def applywavB(self, Wgam):
        """ apply waveform method B (time domain )

        DEPRECATED

        Parameters
        ----------

        Wgam : waveform

        Returns
        -------

        ri  : TUDsignal

            impulse response for each ray separately

        Notes
        ------

            The overall received signal is built in time domain

            Wgam is applied on each Ray Transfer function

        See Also
        --------

        pylayers.signal.bsignal.TUDsignal.ft1

        """

        print DeprecationWarning(
            'WARNING : Tchannel.applywavB is going to be replaced by Tchannel.applywav')

        # product in frequency domain between Channel (self) and waveform
        Y = self.apply(Wgam)
        # back in time domain
        ri = Y.ft1(Nz=500,ffts=1)

        return(ri)

    def applywavA(self, Wgam, Tw):
        """ apply waveform method A

        DEPRECATED

        Parameters
        ----------

        Wgam :
        Tw   :

        The overall received signal is built in frequency domain

        See Also
        --------

        pylayers.signal.bsignal

        """
        print DeprecationWarning(
            'WARNING : Tchannel.applywavA is going to be replaced by Tchannel.applywav')
        Hab = self.H.ft2(0.001)
        HabW = Hab * Wgam
        RI = HabW.symHz(10000)
        ri = RI.ifft(0,'natural')
        ri.translate(-Tw)
        return(ri)


    def plotd (self, d='doa', **kwargs):
        """plot direction of arrival and departure

        Parameters
        ----------

        d: 'doa' | 'dod'
            display direction of departure | arrival
        fig : plt.figure
        ax : plt.axis
        phi: tuple (-180, 180)
            phi angle
        normalize: bool
            energy normalized
        reverse : bool
            inverse theta and phi represenation
        polar : bool
            polar representation
        cmap: matplotlib.cmap
        mode: 'center' | 'mean' | 'in'
            see bsignal.energy
        s : float
            scatter dot size
        fontsize: float
        edgecolors: bool
        colorbar: bool
        title : bool

        """
        defaults = {
                    'fig': [],
                    'ax': [],
                    'phi':(-180, 180),
                    'normalize':False,
                    'reverse' : True,
                    'cmap':plt.cm.hot_r,
                    'mode':'center',
                    's':30,
                    'fontsize':12,
                    'edgecolors':'none',
                    'polar':False,
                    'colorbar':False,
                    'title':False,
                    'xa':[],
                    'xb':[]
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value


        di = getattr(self, d, 'doa')


        # remove non plt.scatter kwargs
        phi = kwargs.pop('phi')
        the = (0,180)
        fontsize = kwargs.pop('fontsize')
        polar = kwargs.pop('polar')
        fig = kwargs.pop('fig')
        ax = kwargs.pop('ax')
        colorbar = kwargs.pop('colorbar')
        reverse = kwargs.pop('reverse')
        normalize = kwargs.pop('normalize')
        mode =kwargs.pop('mode')
        title =kwargs.pop('title')
        xa = kwargs.pop('xa')
        xb = kwargs.pop('xb')
        if fig == []:
            fig = plt.figure()


        Etot = self.energy(mode=mode) + 1e-15


        if normalize:
            Emax = max(Etot)
            Etot = Etot / Emax
        #
        # col  = 1 - (10*log10(Etot)-Emin)/(Emax-Emin)
        # WARNING polar plot require radian angles
        #
        if polar :
            al = 1.
            alb = 180. / np.pi

            phi=np.array(phi)
            the=np.array(the)

            if reverse :
                phi[0] = phi[0]*np.pi/180
                phi[1] = phi[1]*np.pi/180
                the[0] = the[0]
                the[1] = the[1]
            else :
                phi[0] = phi[0]
                phi[1] = phi[1]
                the[0] = the[0]*np.pi/180
                the[1] = the[1]*np.pi/180
        else :
            al = 180. / np.pi
            alb = 180. / np.pi

        col = 10 * np.log10(Etot)
        kwargs['c'] = col
        if len(col) != len(di):
            print "len(col):", len(col)
            print "len(di):", len(dir)
        if ax == []:
            ax = fig.add_subplot(111, polar=polar)
        if reverse :
            scat = ax.scatter(di[:, 1] * al, di[:, 0] * alb, **kwargs)
            ax.axis((phi[0], phi[1], the[0], the[1]))
            ax.set_xlabel('$\phi(^{\circ})$', fontsize=fontsize)
            ax.set_ylabel("$\\theta_t(^{\circ})$", fontsize=fontsize)
        else:
            scat = ax.scatter(di[:, 0] * al, di[:, 1] * alb, **kwargs)
            ax.axis((the[0], the[1], phi[0], phi[1]))
            ax.set_xlabel("$\\theta_t(^{\circ})$", fontsize=fontsize)
            ax.set_ylabel('$\phi(^{\circ})$', fontsize=fontsize)

        if title:
            ax.set_title(d, fontsize=fontsize+2)
        if colorbar:
            b = plt.colorbar(scat,cax=ax)
            if normalize:
                b.set_label('dB')
            else:
                b.set_label('Path Loss (dB)')

            for t in b.ax.get_yticklabels():
                t.set_fontsize(fontsize)

        return (fig, ax)


    def plotad(self,a='phi', **kwargs):
        """plot angular delays

         Parameters
        ----------

        d: 'doa' | 'dod'
            display direction of departure | arrival
        typ : 'ns' | 'm'
            display delays in nano seconds ( ns) or meter (m)
        fig : plt.figure
        ax : plt.axis
        a : str
            angle 'theta' | 'phi'
        normalize: bool
            energy normalized
        reverse : bool
            inverse theta and phi represenation
        polar : bool
            polar representation
        cmap: matplotlib.cmap
        mode: 'center' | 'mean' | 'in'
            see bsignal.energy
        s : float
            scatter dot size
        fontsize: float
        edgecolors: bool
        colorbar: bool
        titel : bool
        'clipval': float
            remove values below clipval in dB
        """
        defaults = { 'fig': [],
                    'ax': [],
                    'normalize':False,
                    'cmap':plt.cm.hot_r,
                    'mode':'center',
                    's':30,
                    'fontsize':12,
                    'edgecolors':'none',
                    'polar':False,
                    'colorbar':False,
                    'taumin':[],
                    'taumax':[],
                    'typ':'m',
                    'title':False,
                    'clipval': -2500,
                    'd':'doa'
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value


        # remove non plt.scatter kwargs

        fontsize = kwargs.pop('fontsize')
        polar = kwargs.pop('polar')
        fig = kwargs.pop('fig')
        ax = kwargs.pop('ax')
        colorbar = kwargs.pop('colorbar')
        normalize = kwargs.pop('normalize')
        mode =kwargs.pop('mode')
        dmin = kwargs.pop('taumin')
        dmax = kwargs.pop('taumax')
        title = kwargs.pop('title')
        typ = kwargs.pop('typ')
        clipval = kwargs.pop('clipval')
        do = kwargs.pop('d')
        if fig == []:
            fig = plt.figure()

        if do=='doa':
            di = self.doa
        elif do=='dod':
            di = self.dod

        if a == 'theta':
            ang = np.array((0,180))
        else :
            ang = np.array((-180,180))

        delay = self.taud
        if typ =='m':
            delay = delay*0.3

        if dmin == []:
            dmin = 0.#min(delay)
        if dmax == []:
            dmax= max(delay)



        if self.isFriis :
            Etot = self.energy(mode=mode,Friis=False) + 1e-15
        else :
            Etot = self.energy(mode=mode,Friis=True) + 1e-15


        if normalize:
            Emax = max(Etot)
            Etot = Etot / Emax
        #
        #
        #
        # col  = 1 - (10*log10(Etot)-Emin)/(Emax-Emin)
        # WARNING polar plot require radian angles
        #
        #
        if polar :
            al = 1.

        else :
            al = 180. / np.pi


        col = 10 * np.log10(Etot)
        cv = np.where(col >= clipval)[0]
        kwargs['c'] = col[cv]
        if len(col) != len(di):
            print "len(col):", len(col)
            print "len(di):", len(dir)
        if ax == []:
            ax = fig.add_subplot(111, polar=polar)
        if a == 'phi':
            scat = ax.scatter(di[cv, 1] * al, delay[cv], **kwargs)
            ax.axis((ang[0], ang[1], dmin, dmax))

            ax.set_xlabel(r"$\phi(^{\circ})$", fontsize=fontsize)
            if typ == 'm' :
                ax.set_ylabel("distance (m)", fontsize=fontsize-2)
            else :
                ax.set_ylabel(r"$\phi(^{\circ})$", fontsize=fontsize-2)
        elif a == 'theta':
            scat = ax.scatter(di[cv, 0] * al, delay[cv], **kwargs)
            ax.axis((ang[0], ang[1], dmin,dmax))
            ax.set_xlabel(r"$\\theta_t(^{\circ})$", fontsize=fontsize)
            if typ == 'm' :
                ax.set_ylabel("distance (m)", fontsize=fontsize-2)
            else :
                ax.set_ylabel(r"$\phi(^{\circ})$", fontsize=fontsize-2)
            if title :
                ax.set_title('DoA vs delay (ns)', fontsize=fontsize+2)
        if colorbar:
            b=fig.colorbar(scat)
            if normalize:
                b.set_label('dB')
            else:
                b.set_label('Path Loss (dB)')

        return (fig, ax)


    def doadod(self, **kwargs):
        """ doadod scatter plot

        Parameters
        ----------

        phi: tuple (-180, 180)
            phi angle
        normalize: bool
            energy normalized
        reverse : bool
            inverse theta and phi represenation
        polar : bool
            polar representation
        cmap: matplotlib.cmap
        mode: 'center' | 'mean' | 'in'
            see bsignal.energy
        s : float
            scatter dot size
        fontsize: float
        edgecolors: bool
        colorbar bool

        Summary
        --------

        scatter plot of the DoA-DoD channel structure
        the energy is colorcoded over all couples of DoA-DoD


        """
        defaults = {
                    'phi':(-180, 180),
                    'normalize':False,
                    'reverse' : True,
                    'cmap':plt.cm.hot_r,
                    'mode':'center',
                    's':30,
                    'fontsize':12,
                    'edgecolors':'none',
                    'polar':False,
                    'mode':'mean',
                    'xa':0,
                    'xb':0
                    }


        fig = plt.figure()
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        ax1  = fig.add_subplot(121,polar=kwargs['polar'])
        ax2  = fig.add_subplot(122,polar=kwargs['polar'])


        if kwargs['xa']<kwargs['xb']:
            fig,ax = self.plotd(d='dod',fig=fig,ax=ax1,**kwargs)
            fig,ax = self.plotd(d='doa',fig=fig,ax=ax2,**kwargs)
        else:
            fig,ax = self.plotd(d='doa',fig=fig,ax=ax1,**kwargs)
            fig,ax = self.plotd(d='dod',fig=fig,ax=ax2,**kwargs)
        return fig,ax

    def field(self):

        tau  = self.tk[:,None,None,None]
        fGHz = self.x[None,None,None,:]
        E = np.exp(-2*1j*tau*fGHz)
        F = self.y*E
        return np.sum(F,axis=0)
        #f = bs.FUsignal(x=self.x,y=np.sum(F,axis=0))
        #return(f)


    def energy(self,mode='mean',Friis=True,sumray=False):
        """ calculates channel energy including antennas spatial filtering

        Parameters
        ----------

        mode : string
            center | mean | integ    (different manner to get the value)
        Friis : boolean
            apply the Frris coeff(2/(4p pi f)
        sumray: boolean
            ray energy cummulation indicator

        """
        #
        #  r x f
        #  axis 1 : ray
        #  axis 1 : frequency
        #
        Etot = bs.FUsignal.energy(self,axis=1,mode=mode,Friis=Friis)
        if sumray:
            Etot = np.sum(Etot,axis=0)
        return Etot

    def wavefig(self, w, Nray=5):
        """ display

        Parameters
        ----------

        w      :  waveform
        Nray   :  int
            number of rays to be displayed

        """
        # Construire W
        W = w.ft()
        # Appliquer W
        Y = self.apply(W)
        # r.require('graphics')
        # r.postscript('fig.eps')
        # r('par(mfrow=c(2,2))')
        # Y.fig(Nray)
        y = Y.iftd(100, 0, 50, 0)
        y.fig(Nray)
        # r.dev_off()
        # os.system("gv fig.eps ")
        # y.fidec()
        # Sur le FUsignal retourn
        # A gauche afficher le signal sur chaque rayon
        # A droite le meme signal decal
        # En bas a droite le signal resultant

    def rayfig(self, k, W, col='red'):
        """ build a figure with rays

        Parameters
        ----------

        k : ray index
        W : waveform    (FUsignal)

        Notes
        -----

        W is apply on k-th ray and the received signal is built in time domain

        """
        # get the kth Ray  Transfer function
        Hk = bs.FUDsignal(self.H.x, self.H.y[k,:])

        dxh = Hk.dx()
        dxw = W.dx()
        w0 = W.x[0]    # fmin W
        hk0 = Hk.x[0]   # fmin Hk

        # on s'arrange pour que hk0 soit egal a w0 (ou hk0 soit legerement inferieur a w0)
        if w0 < hk0:
            np = ceil((hk0 - w0) / dxh)
            hk0_new = hk0 - np * dxh
            x = arange(hk0_new, hk0 + dxh, dxh)[0:-1]
            Hk.x = hstack((x, Hk.x))
            Hk.y = hstack((zeros(np), Hk.y))

        if (abs(dxh - dxw) > 1e-10):
            if (dxh < dxw):
                # reinterpolate w
                print " resampling w"
                x_new = arange(W.x[0], W.x[-1] + dxh, dxh)[0:-1]
                Wk = W.resample(x_new)
                dx = dxh
            else:
                # reinterpolate h
                print " resampling h"
                x_new = arange(Hk.x[0], Hk.x[-1] + dxw, dxw)[0:-1]
                Hk = Hk.resample(x_new)
                dx = dxw
                Wk = W

        # qHk.x[0]==Wk.x[0]

    def rssi(self,ufreq=0) :
        """ Compute RSSI value for a frequency index

        Parameters
        ----------

        ufreq : int
            index in the frequency range


        Returns
        -------

        RSSI: float
        RSSI value in dB

        Notes
        -----

        This function will be deprecated by energy function

        """

        Ak   = self.y[:, ufreq]
        Pr   = np.sum(Ak*np.conj(Ak))
        akp   = Ak*np.exp(-2*1j*np.pi*self.x[ufreq]*self.tk)
        Prp   = np.abs(np.sum(akp))**2
        PrdB  = 10*np.log10(Pr)
        PrpdB = 10*np.log10(Prp)

        return PrdB,PrpdB

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
            import ipdb
            ipdb.set_trace()
            u = np.argsort(E,axis=0)[::-1]
            u = u[:,0,0]

        self.taud = self.taud[u]
        self.taue = self.taue[u]
        self.doa = self.doa[u]
        self.dod = self.dod[u]
        self.y = self.y[u,...]

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
        WMHz : float
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
            >>> H = Tchannel(x=fGHz,y=y,taud=np.array([15,17,18]))
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


        """
        y = fft.ifft(self.y)
        T = 1/(self.x[1]-self.x[0])
        x = np.linspace(0,T,len(self.x))
        h = TUDchannel(x,y,self.taud,self.taue)
        return(h)


    def totime(self, Nz=1, ffts=0):
        """ transform to TUDchannel

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
        h = TUDchannel(s.x, fft.fftshift(s.y), self.taud,self.taue)
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
        rf = bs.TUsignal(x_new, yini)
        #
        # initializes a void signal
        #
        for i in range(Nray):
            r = bs.TUsignal(x_new, np.zeros(len(x_new)))
            si = bs.TUsignal(x, s.y[i, :])
            si.translate(tau[i])
            r = r + si
            rf.y[i, :] = r.y
        return rf

    def rir(self, Nz, ffts=0):
        """  construct ray impulse response

        Parameters
        ----------

        Nz   : number of zeros for zero padding
        ffts : fftshift indicator
            0  no fftshift
            1  apply fftshift

        Returns
        -------


        rir : TUsignal


        See Also
        --------

        pylayers.signal.bsignal.


        """

        tau = self.taud + self.taue
        taumin = min(tau) 
        taumax = max(tau)
        dtau = (taumax-taumin)
        self.s = self.ift(Nz, ffts)


        shy = self.s.y.shape
        dx = self.s.x[1]-self.s.x[0]
        N  = np.ceil(dtau/dx)+shy[-1]
        itau = np.floor((tau-taumin)/dx).astype(int)
        
        U = np.ones((shy[0],shy[-1]),dtype=int)
        CU = np.cumsum(U,axis=1)-1 #-1 to start @ value 0 

        rir  = np.zeros((shy[0],N))
        col1 = np.repeat(np.arange(shy[0],dtype=int),shy[-1])
        col2 = (CU+itau[:,None]).ravel()
        index = np.vstack((col1,col2)).T
    
        rir[index[:,0],index[:,1]] = self.s.y.ravel()
        t = np.linspace(taumin,taumax,N)
        return bs.TUsignal(x=t, y=rir)


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
        tau = self.taud + self.taue
        self.s = self.ift(Nz, ffts)
        x = self.s.x
        r = bs.TUsignal(x=x, y=np.zeros(self.s.y.shape[1:]))


        if len(tau) == 1:
            return(self.s)
        else:
            for i in range(len(tau)):
                si = bs.TUsignal(self.s.x, self.s.y[i, :])
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
        r  = bs.TUsignal(x, np.zeros(len(x)))
        si = bs.TUsignal(s.x, s.y[k, :])
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
            >>> s = Tchannel(x=fGHz,y=alpha,taud=taud)
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

        U = bs.FUsignal(f, np.zeros(len(f)))

        TAUF = np.outer(tau, f)
        E = np.exp(-2 * 1j * np.pi * TAUF)

        S = self.resample(f)
        ES = E * S.y
        V = sum(ES, axis=0)
        U.y = V

        return U

    def frombuf(self,S,sign=-1):
        """ load a buffer from vna

        Parameters
        ----------

        S : buffer
        sign : int (+1 |-1)  for complex reconstruction

        """
        N = len(self.x)
        u = np.arange(0,N)*2
        v = np.arange(0,N)*2+1
        S21 = (S[u]+sign*1j*S[v]).reshape((1,N))
        self.y = S21

    def capacity(self,Pt,T=290,mode='blast'):
            """  calculates channel Shannon capacity (no csi)

            Parameters
            ----------

            Pt : Power transmitted
            T : Temperrature Kelvin
            mode : string

            Returns
            -------

            C : Channel capacity (bit/s)

            """
            kB = 1.3806488e-23
            N0 = kB*T
            dfGHz = self.x[1]-self.x[0]
            BGHz  = self.x[-1]-self.x[0]
            Pb = N0*BGHz*1e9
            H2 = self.y*np.conj(self.y)
            snr = Pt[:,None]*H2[None,:]/Pb
            c = np.log(1+snr)/np.log(2)
            C = np.sum(c,axis=1)*dfGHz
            SNR = np.sum(snr,axis=1)*dfGHz

            return(C,SNR)

    def calibrate(self,filecal='calibration.mat',conjugate=False):
        """ calibrate

        Parameters
        ----------

        filecal : string
            calibration file name  "calibration.mat"
        conjugate : boolean
            default False


        """
        self.filecal = filecal
        Hcal = Tchannel()
        Hcal.load(filecal)
        assert (len(self.x) == len(Hcal.x)),"calibration file has not the same number of points"
        if not self.calibrated:
            if not(conjugate):
                self.y = self.y/Hcal.y
            else:
                self.y = self.y/np.conj(Hcal.y)
            self.calibrated = not self.calibrated
        else:
            if not(conjugate):
                self.y = self.y*Hcal.y
            else:
                self.y = self.y*np.conj(Hcal.y)
            self.calibrated = not self.calibrated

    def pdp(self,win='hamming',calibrate=True):
        """ calculates power delay profile

        Parameters
        ----------

        win : string
            window name
        """
        self.win = win
        if calibrate and not self.calibrated:
            self.calibrate()

        if not self.windowed:
            self.window(win=win)

        # inverse Fourier transform

        pdp = self.ift(ffts=1)
        return pdp

class Ctilde(PyLayers):
    """ container for the 4 components of the polarimetric ray channel

    Attributes
    ----------

    Ctt : bsignal.FUsignal
    Ctp : bsignal.FUsignal
    Cpt : bsignal.FUsignal
    Cpp : bsignal.FUsignal

    tauk : ndarray delays
    tang : ndarray angles of departure
    rang : ndarray angles of arrival
    tangl : ndarray angles of departure (local)
    rangl : ndarray angles of arrival (local)

    fGHz : np.array
        frequency array
    nfreq : int
        number of frequency point
    nray  : int
        number of rays

    Methods
    -------

    choose
    load
    mobility
    doadod
    show
    energy
    sort
    prop2tran

    """
    def __init__(self):
        """ class constructor

        Notes
        -----

        transpose == False   (r,f)
        transpose == True    (f,r)

        A Ctilde object can be :
            + returned from eval method of a Rays object.
            + generated from a statistical model of the propagation channel

        """
        self.fail = False
        self.islocal = False
        self.Tt = np.eye(3)
        self.Tr = np.eye(3)

    def __repr__(self):
        s = 'Ctilde : Ray Propagation Channel Matrices'+'\n---------\n'
        if hasattr(self, 'Cpp'):
            s = s + str(np.shape(self.Cpp.y))+'\n'
        if hasattr(self, 'nray'):
            s = s + 'Nray : ' + str(self.nray)+'\n'
            s = s + 'fmin(GHz) : ' + str(self.Cpp.x[0])+'\n'
            s = s + 'fmax(GHz): ' + str(self.Cpp.x[-1])+'\n'
            s = s + 'Nfreq : ' + str(self.nfreq)+'\n'
        return(s)

    


    def saveh5(self,Lfilename,idx,a,b):
        """ save Ctilde object in hdf5 format

        Parameters
        ----------

        Lfilename  : string
            Layout filename
        idx : int
            file identifier number
        a : np.ndarray
            postion of point a (transmitter)
        b : np.ndarray
            postion of point b (receiver)


        """

        Lfilename=Lfilename.split('.')[0]
        _filename= Lfilename +'_' + str(idx).zfill(5) + '.hdf5'

        filename=pyu.getlong(_filename,pstruc['DIRCT'])

        # save channel in global basis
        if self.islocal:
            self.locbas(b2g=True)

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f=h5py.File(filename,'w')
            f.create_dataset('Tt',shape=np.shape(self.Tt),data=self.Tt)
            f.create_dataset('Tr',shape=np.shape(self.Tr),data=self.Tr)
            f.create_dataset('tang',shape=np.shape(self.tang),data=self.tang)
            f.create_dataset('rang',shape=np.shape(self.rang),data=self.rang)
            f.create_dataset('tauk',shape=np.shape(self.tauk),data=self.tauk)

            f.create_dataset('fGHz',shape=np.shape(self.fGHz),data=self.fGHz)


            f.create_dataset('Ctt_y',shape=np.shape(self.Ctt.y),data=self.Ctt.y)
            f.create_dataset('Cpp_y',shape=np.shape(self.Cpp.y),data=self.Cpp.y)
            f.create_dataset('Cpt_y',shape=np.shape(self.Cpt.y),data=self.Cpt.y)
            f.create_dataset('Ctp_y',shape=np.shape(self.Ctp.y),data=self.Ctp.y)

            f.create_dataset('Tx',shape=np.shape(a),data=a)
            f.create_dataset('Rx',shape=np.shape(b),data=b)

            f.close()
        except:
            f.close()
            raise NameError('Channel.Ctilde: issue when writting h5py file')




    def loadh5(self,Lfilename,idx,output=True):
        """ load Ctilde object in hdf5 format

        Parameters
        ----------

        Lfilename  : string
            Layout filename
        idx : int
            file identifier number
        output : bool
            return an output precised in return

        Returns
        -------

        if output:
        (Layout filename , Tx position, Rx position)

        """

        _Lfilename=Lfilename.split('.')[0]
        _filename= _Lfilename +'_' + str(idx).zfill(5) + '.hdf5'
        filename=pyu.getlong(_filename,pstruc['DIRCT'])

        try:
            f=h5py.File(filename,'r')
            self.fGHz = f['fGHz'][:]
            self.tang = f['tang'][:]
            self.rang = f['rang'][:]
            self.tauk = f['tauk'][:]

            self.Tt = f['Tt'][:]
            self.Tr = f['Tr'][:]

            Ctt = f['Ctt_y'][:]
            Cpp = f['Cpp_y'][:]
            Ctp = f['Ctp_y'][:]
            Cpt = f['Cpt_y'][:]

            self.Ctt = bs.FUsignal(self.fGHz, Ctt)
            self.Ctp = bs.FUsignal(self.fGHz, Ctp)
            self.Cpt = bs.FUsignal(self.fGHz, Cpt)
            self.Cpp = bs.FUsignal(self.fGHz, Cpp)
            tx = f['Tx'][:]
            rx = f['Rx'][:]


            self.nfreq = len(self.fGHz)
            self.nray = np.shape(self.Cpp.y)[0]

            f.close()
        except:
            f.close()
            raise NameError('Channel.Ctilde: issue when reading h5py file')

        if output :
            return (Lfilename ,tx,rx)


    def _saveh5(self,filenameh5,grpname):
        """ save Ctilde object in hdf5 format compliant with Link Class

        Parameters
        ----------

        filenameh5  : str
            file name of h5py file Link format
        grpname  : int
            groupname in filenameh5


        """

        if self.islocal:
            self.locbas(b2g=True)


        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:

            fh5=h5py.File(filename,'a')
            if not grpname in fh5['Ct'].keys():
                fh5['Ct'].create_group(grpname)
            else :
                print 'Warning : Ct/'+grpname +'already exists in '+filenameh5
            f=fh5['Ct/'+grpname]

            # save channel in global basis
            f.create_dataset('Tt',shape=np.shape(self.Tt),data=self.Tt)
            f.create_dataset('Tr',shape=np.shape(self.Tr),data=self.Tr)
            f.create_dataset('tang',shape=np.shape(self.tang),data=self.tang)
            f.create_dataset('rang',shape=np.shape(self.rang),data=self.rang)
            f.create_dataset('tauk',shape=np.shape(self.tauk),data=self.tauk)

            f.create_dataset('fGHz',shape=np.shape(self.fGHz),data=self.fGHz)


            f.create_dataset('Ctt_y',shape=np.shape(self.Ctt.y),data=self.Ctt.y)
            f.create_dataset('Cpp_y',shape=np.shape(self.Cpp.y),data=self.Cpp.y)
            f.create_dataset('Cpt_y',shape=np.shape(self.Cpt.y),data=self.Cpt.y)
            f.create_dataset('Ctp_y',shape=np.shape(self.Ctp.y),data=self.Ctp.y)

            fh5.close()
        except:
            fh5.close()
            raise NameError('Channel.Ctilde: issue when writting h5py file')

    def los(self,pa=np.r_[0,0,0],pb=np.r_[3,0,0],fGHz=np.r_[2.4],tang=np.r_[[1.57,3.14]],rang=np.r_[[1.57,0]]):
        """ Line of site channel

        Parameters
        ----------

        d(m)
        fGHz (,Nf)
        tang (1x2)
        rang (1x2)

        """
        self.pa = pa
        self.pb = pb
        self.fGHz = fGHz
        self.nfreq = len(fGHz)
        self.nray = 1
        si = pb-pa
        d = np.r_[np.sqrt(np.sum(si*si))]
        self.tauk = d/0.3
        tht =  np.arccos(si[2])
        pht =  np.arctan2(si[1],si[0])
        thr =  np.arccos(-si[2])
        phr =  np.arctan2(-si[1],-si[0])
        self.tang = np.array([tht,pht]).reshape((1,2))
        self.rang = np.array([thr,phr]).reshape((1,2))
        U = np.ones(len(self.fGHz),dtype=complex)/d[0]
        Z = np.zeros(len(self.fGHz),dtype=complex)
        self.Ctt = bs.FUsignal(self.fGHz, U)
        self.Ctp = bs.FUsignal(self.fGHz, Z)
        self.Cpt = bs.FUsignal(self.fGHz, Z)
        self.Cpp = bs.FUsignal(self.fGHz, U)

    def _loadh5(self,filenameh5,grpname):
        """ load Ctilde object in hdf5 format

        Parameters
        ----------

        filenameh5  : str
            file name of h5py file Link format
        grpname  : int
            groupname in filenameh5


        """

        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        try:
            fh5=h5py.File(filename,'r')
            f = fh5['Ct/'+grpname]

            self.fGHz = f['fGHz'][:]
            self.tang = f['tang'][:]
            self.rang = f['rang'][:]
            self.tauk = f['tauk'][:]

            self.Tt = f['Tt'][:]
            self.Tr = f['Tr'][:]
            Ctt = f['Ctt_y'][:]
            Cpp = f['Cpp_y'][:]
            Ctp = f['Ctp_y'][:]
            Cpt = f['Cpt_y'][:]

            self.Ctt = bs.FUsignal(self.fGHz, Ctt)
            self.Ctp = bs.FUsignal(self.fGHz, Ctp)
            self.Cpt = bs.FUsignal(self.fGHz, Cpt)
            self.Cpp = bs.FUsignal(self.fGHz, Cpp)


            self.nfreq = len(self.fGHz)
            self.nray = np.shape(self.Cpp.y)[0]

            fh5.close()
        except:
            fh5.close()
            raise NameError('Channel.Ctilde: issue when reading h5py file')

    def load(self, filefield, transpose=False):
        """ load a Ctilde from a .field file

        Load the three files .tauk .tang .rang which contain respectively
        delay , angle of departure , angle of arrival.
maicher
        Parameters
        ----------

        filefield  : string
        transpose  : boolean
            default False

        Examples
        --------

        >>> from pylayers.antprop.channel import *
        >>> from pylayers.simul.simulem import *
        >>> S = Simul()
        >>> S.load('default.ini')
        >>> C = Ctilde()
        >>> out = C.load(pyu.getlong(S.dtud[1][1],'output'))

        """
        filetauk = filefield.replace('.field', '.tauk')
        filetang = filefield.replace('.field', '.tang')
        filerang = filefield.replace('.field', '.rang')
        try:
            fo = open(filefield, "rb")
        except:
            raise NameError( "file " + filefield + " is unreachable")

        # decode filename (*.field file obtained from evalfield simulation)
        nray = stru.unpack('i', fo.read(4))[0]
        nfreq = stru.unpack('i', fo.read(4))[0]
        if nfreq == 0:
            print " Error : incorrect number of frequency points in .field"
            self.fail = True
            return

        n = nray * nfreq
        buf = fo.read()
        fo.close()

        CMat = np.ndarray(shape=(n, 8), buffer=buf)
        c11 = CMat[:, 0] + CMat[:, 1]*1j
        c12 = CMat[:, 2] + CMat[:, 3]*1j
        c21 = CMat[:, 4] + CMat[:, 5]*1j
        c22 = CMat[:, 6] + CMat[:, 7]*1j

        c11 = c11.reshape(nray, nfreq)
        c12 = c12.reshape(nray, nfreq)
        c21 = c21.reshape(nray, nfreq)
        c22 = c22.reshape(nray, nfreq)

        if transpose:
            c11 = c11.transpose()
            c12 = c12.transpose()
            c21 = c21.transpose()
            c22 = c22.transpose()

        #
        # Temporary freq --> read filefreq
        #
        fGHz = np.linspace(2, 11, nfreq)

        self.Ctt = bs.FUsignal(fGHz, c11)
        self.Ctp = bs.FUsignal(fGHz, c12)
        self.Cpt = bs.FUsignal(fGHz, c21)
        self.Cpp = bs.FUsignal(fGHz, c22)
        self.nfreq = nfreq
        self.nray = nray

        try:
            fo = open(filetauk, "rb")
        except:
            self.fail = True
            print "file ", filetauk, " is unreachable"
        # decode filetauk
        if not self.fail:
            nray_tauk = stru.unpack('i', fo.read(4))[0]
            #print "nb rays in .tauk file: ", nray_tauk
            buf = fo.read()
            fo.close()
            nray = len(buf) / 8
            #print "nb rays 2: ", nray
            self.tauk = np.ndarray(shape=nray, buffer=buf)
            #if nray_tauk != nray:
            #    print nray_tauk - nray
        self.tauk = self.tauk

        # decode the angular files (.tang and .rang)
        try:
            fo = open(filetang, "rb")
        except:
            self.fail = True
            print "file ", filetang, " is unreachable"
        if not self.fail:
            nray_tang = stru.unpack('i', fo.read(4))[0]
            buf = fo.read()
            fo.close()
            # coorectif Bug evalfield
            tmp = np.ndarray(shape=(nray_tang, 2), buffer=buf)
            self.tang = tmp[0:nray,:]
        try:
            fo = open(filerang, "rb")
        except:
            self.fail = True
            print "file ", filerang, " is unreachable"

        if not self.fail:
            nray_rang = stru.unpack('i', fo.read(4))[0]
            buf = fo.read()
            fo.close()
            # corectif Bug evalfield
            tmp = np.ndarray(shape=(nray_rang, 2), buffer=buf)
            self.rang = tmp[0:nray,:]

    def mobility(self, v, dt):
        """ modify channel for uniform mobility

        Parameters
        ----------

        v  : float
            velocity (m/s)
        dt : float
            delta t (s)

        Notes
        -----

        Calculate a channel field from Ctilde and v(terminal vitese)
        and dt(time of deplacement)

        dt en s  (observation time between 2 Rx position)
        v en m/s (vitesse de changement de Rx)

        Returns
        -------

        tau : modified Ctilde

        """

        c = 0.3  # m/ns celerity of light
        tauk = self.tauk
        tang = self.tang
        rang = self.rang

        rk = tauk * c
        rk_mod = abs(rk)
        sk_ch = rk / rk_mod

        # cos_alph =dot(v/abs(v),sk_ch)

        cos_alph = (v * sk_ch) / abs(v)
        self.cos_alph = cos_alph
        rk_ch = rk_mod * cos_alph * abs(v) * dt
        sk_ch_ch = (rk + v * dt) / (rk_ch + cos_alph * abs(v) * dt)
        tauk_ch = (abs(rk_ch) * sk_ch_ch) / c

        return(tauk_ch)

    def plotd (self, d='doa', **kwargs):
        """ plot direction of arrival/departure

        Parameters
        ----------

        d: string
            'doa' | 'dod'
            display direction of departure | arrival
        fig : plt.figure
        ax : plt.axis
        phi: tuple (-180, 180)
            phi angle
        normalize: bool
            energy normalized
        reverse : bool
            inverse theta and phi represenation
        polar : bool
            polar representation
        cmap: matplotlib.cmap
        mode: 'center' | 'mean' | 'in'
            see bsignal.energy
        s : float
            scatter dot size
        fontsize: float
        edgecolors: bool
        colorbar: bool
        title : bool

        """
        defaults = {
                    'fig': [],
                    'ax': [],
                    'phi':(-180, 180),
                    'normalize':False,
                    'reverse' : True,
                    'cmap':plt.cm.hot_r,
                    'mode':'center',
                    's':30,
                    'fontsize':12,
                    'edgecolors':'none',
                    'polar':False,
                    'colorbar':False,
                    'title' : False
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value


        if d =='dod':
            tit = 'DOD : A'
            di = getattr(self, 'tang')
        elif d == 'doa':
            tit = 'DOA : B'
            di = getattr(self, 'rang')
        else :
            raise AttributeError('d attribute can only be doa or dod')


        # remove non plt.scatter kwargs
        phi = kwargs.pop('phi')
        the = (0,180)
        fontsize = kwargs.pop('fontsize')
        polar = kwargs.pop('polar')
        fig = kwargs.pop('fig')
        ax = kwargs.pop('ax')
        colorbar = kwargs.pop('colorbar')
        reverse = kwargs.pop('reverse')
        normalize = kwargs.pop('normalize')
        mode = kwargs.pop('mode')
        title = kwargs.pop('title')

        if fig == []:
            fig = plt.figure()


        Ett, Epp, Etp, Ept = self.energy(mode=mode,Friis=True)
        Etot = Ett+Epp+Etp+Ept + 1e-15

        if normalize:
            Emax = max(Etot)
            Etot = Etot / Emax
        #
        #
        #
        # col  = 1 - (10*log10(Etot)-Emin)/(Emax-Emin)
        # WARNING polar plot require radian angles
        if polar :
            al = 1.
            alb = 180. / np.pi

            phi=np.array(phi)
            the=np.array(the)

            if reverse :
                phi[0] = phi[0]*np.pi/180
                phi[1] = phi[1]*np.pi/180
                the[0] = the[0]
                the[1] = the[1]
            else :
                phi[0] = phi[0]
                phi[1] = phi[1]
                the[0] = the[0]*np.pi/180
                the[1] = the[1]*np.pi/180
        else :
            al  = 180. / np.pi
            alb = 180. / np.pi

        col = 10 * np.log10(Etot)
        kwargs['c'] = col
        if len(col) != len(di):
            print "len(col):", len(col)
            print "len(di):", len(dir)
        if ax == []:
            ax = fig.add_subplot(111, polar=polar)

        if reverse :
            scat = ax.scatter(di[:, 1] * al, di[:, 0] * alb, **kwargs)
            ax.axis((phi[0], phi[1], the[0], the[1]))
            ax.set_xlabel('$\phi(^{\circ})$', fontsize=fontsize)
            ax.set_ylabel("$\\theta_t(^{\circ})$", fontsize=fontsize)

        else:
            scat = ax.scatter(di[:, 0] * al, di[:, 1] * alb, **kwargs)
            ax.axis((the[0], the[1], phi[0], phi[1]))
            ax.set_xlabel("$\\theta_t(^{\circ})$", fontsize=fontsize)
            ax.set_ylabel('$\phi(^{\circ})$', fontsize=fontsize)

        if title:
            ax.set_title(tit, fontsize=fontsize+2)

        ll = ax.get_xticklabels()+ax.get_yticklabels()
        for l in ll:
            l.set_fontsize(fontsize)

        if colorbar:
            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("right",size="5%",pad=0.05)
            clb = plt.colorbar(scat,ax=ax)
            if normalize:
                clb.set_label('dB',size=fontsize)
            else:
                clb.set_label('Path Loss (dB)',size=fontsize)

            for t in clb.ax.get_yticklabels():
                t.set_fontsize(fontsize)

        return (fig, ax)

    def doadod(self, **kwargs):
        """ doadod scatter plot

        Parameters
        ----------

        phi : tuple (-180, 180)
            phi angle
        normalize : bool
            energy normalized
        reverse : bool
            inverse theta and phi representation
        polar : bool
            polar representation
        cmap : matplotlib.cmap
        mode : string
            'center' | 'mean' | 'in'
        s : float
            scatter dot size
        fontsize : float
        edgecolors : bool
        colorbar : bool
        xa :
        xb :

        Summary
        --------

        scatter plot of the DoA-DoD channel structure
        the energy is color coded over all couples of DoA-DoD

        Examples
        --------

        >>> from pylayers.antprop.channel import *


        See Also
        --------

        pylayers.signal.bsignal.energy

        """
        defaults = {
                    'phi':(-180, 180),
                    'normalize':False,
                    'reverse' : True,
                    'cmap':plt.cm.hot_r,
                    'mode':'center',
                    's':30,
                    'fontsize':12,
                    'edgecolors':'none',
                    'polar':False,
                    'mode':'mean',
                    'colorbar':False,
                    'xa':0,
                    'xb':1
                    }


        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        xa = kwargs.pop('xa')
        xb = kwargs.pop('xb')

        if 'fig' not in kwargs:
            fig = plt.gcf()
            kwargs['fig']=fig
        else:
            fig = kwargs['fig']


        ax1  = fig.add_subplot(121,polar=kwargs['polar'])
        ax2  = fig.add_subplot(122,polar=kwargs['polar'])

        if xa<xb:
            fig,ax1 = self.plotd(d='dod',ax=ax1,**kwargs)
            fig,ax2 = self.plotd(d='doa',ax=ax2,**kwargs)
        else:
            fig,ax1 = self.plotd(d='doa',ax=ax1,**kwargs)
            fig,ax2 = self.plotd(d='dod',ax=ax2,**kwargs)

        return fig,[ax1,ax2]


    def locbas(self, Tt=[], Tr=[],b2g=False):
        """ global reference frame to local reference frame

        If Tt and Tr are [] the global channel is  retrieved

        Parameters
        ----------

        Tt  : Tx rotation matrix 3x3
            default []
        Tr  : Rx rotation matrix 3x3
            default []
        b2g: bool
            back to global reference frame

        Returns
        -------

        Cl : Ctilde local/global
            depends on self.islocal boolean value

        Examples
        --------

        """
        # get Ctilde frequency axes

        fGHz = self.fGHz
        # if rot matrices are passed
        if (Tt != []) & (Tr != []):
            if self.islocal:
                if (hasattr(self,'Tt')) & (hasattr(self,'Tr')):
                    # run locbas to return to global basis
                    self.locbas(b2g=True)
                else:
                    raise NameError('Channel has no self.Tt or self.Tr')
            self.Tt = Tt
            self.Tr = Tr
            self.islocal = True
        # if a return to global is requested
        elif b2g:
            if self.islocal :
                if (hasattr(self,'Tt')) & (hasattr(self,'Tr')):
                    self.Tt = self.Tt.transpose()
                    self.Tr = self.Tr.transpose()
                    self.islocal = False
                else:
                    raise NameError ('self.Tt and self.Tr should exist')
            else:
                print "nothing to do to return in global basis"
                return self
        # if Tt and Tr == []
        else:
            return self

        # get angular axes
        # Rt (2x2)
        # Rr (2x2)
        #
        # tang : r x 2
        # rang : r x 2
        #
        # Rt : 2 x 2 x r
        # Rr : 2 x 2 x r
        #
        # tangl : r x 2
        # rangl : r x 2
        #


        Rt, tangl = geu.BTB_tx(self.tang, self.Tt)
        Rr, rangl = geu.BTB_rx(self.rang, self.Tr)
        #
        # update direction of departure and arrival
        #

        self.tangl = tangl
        self.rangl = rangl

        #uf = np.ones(self.nfreq)

        #
        # r0 : r x 1(f)
        #

        #r0 = np.outer(Rr[0, 0,:], uf)
        r0 = Rr[0,0,:][:, np.newaxis]
        #r1 = np.outer(Rr[0, 1,:], uf)
        r1 = Rr[0,1,:][:, np.newaxis]

        t00 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t01 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rr[1, 0,:], uf)
        r0 = Rr[1, 0,:][:, np.newaxis]
        #r1 = np.outer(Rr[1, 1,:], uf)
        r1 = Rr[1, 1,:][:, np.newaxis]

        t10 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t11 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rt[0, 0,:], uf)
        r0 = Rt[0, 0, :][:, np.newaxis]
        #r1 = np.outer(Rt[1, 0,:], uf)
        r1 = Rt[1, 0, :][:, np.newaxis]

        Cttl = t00 * r0 + t01 * r1
        Cptl = t10 * r0 + t11 * r1

        #r0 = np.outer(Rt[0, 1,:], uf)
        r0 = Rt[0, 1, :][:, np.newaxis]
        #r1 = np.outer(Rt[1, 1,:], uf)
        r1 = Rt[1, 1, :][:, np.newaxis]

        Ctpl = t00 * r0 + t01 * r1
        Cppl = t10 * r0 + t11 * r1

        self.Ctt = bs.FUsignal(fGHz, Cttl)
        self.Ctp = bs.FUsignal(fGHz, Ctpl)
        self.Cpt = bs.FUsignal(fGHz, Cptl)
        self.Cpp = bs.FUsignal(fGHz, Cppl)

        return self

    def Cg2Cl(self, Tt=[], Tr=[]):
        """ global reference frame to local reference frame

        If Tt and Tr are [] the global channel is  retrieved

        Parameters
        ----------

        Tt  : Tx rotation matrix 3x3
            default []
        Tr  : Rx rotation matrix 3x3
            default []

        Returns
        -------

        Cl : Ctilde local

        Examples
        --------

        """
        # get frequency axes

        fGHz = self.fGHz

        if (Tt !=[]) & (Tr!=[]):
            self.Tt = Tt
            self.Tr = Tr
        else:
            if (hasattr(self,'Tt')) & (hasattr(self, 'Tr')):
                self.Tt = self.Tt.transpose()
                self.Tr = self.Tr.transpose()
            else:
                return

        # get angular axes
        # Rt (2x2)
        # Rr (2x2)
        #
        # tang : r x 2
        # rang : r x 2
        #
        # Rt : 2 x 2 x r
        # Rr : 2 x 2 x r
        #
        # tangl : r x 2
        # rangl : r x 2
        #

        Rt, tangl = geu.BTB_tx(self.tang, self.Tt)
        Rr, rangl = geu.BTB_rx(self.rang, self.Tr)

        #
        # update direction of departure and arrival
        #

        self.tang = tangl
        self.rang = rangl

        #uf = np.ones(self.nfreq)

        #
        # r0 : r x 1(f)
        #

        #r0 = np.outer(Rr[0, 0,:], uf)
        r0 = Rr[0,0,:][:,np.newaxis]
        #r1 = np.outer(Rr[0, 1,:], uf)
        r1 = Rr[0,1,:][:,np.newaxis]

        t00 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t01 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rr[1, 0,:], uf)
        r0 = Rr[1, 0,:][:,np.newaxis]
        #r1 = np.outer(Rr[1, 1,:], uf)
        r1 = Rr[1, 1,:][:,np.newaxis]

        t10 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t11 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rt[0, 0,:], uf)
        r0 = Rt[0,0,:][:,np.newaxis]
        #r1 = np.outer(Rt[1, 0,:], uf)
        r1 = Rt[1,0,:][:,np.newaxis]

        Cttl = t00 * r0 + t01 * r1
        Cptl = t10 * r0 + t11 * r1

        #r0 = np.outer(Rt[0, 1,:], uf)
        r0 = Rt[0,1,:][:,np.newaxis]
        #r1 = np.outer(Rt[1, 1,:], uf)
        r1 = Rt[1,1,:][:,np.newaxis]

        Ctpl = t00 * r0 + t01 * r1
        Cppl = t10 * r0 + t11 * r1


        self.Ctt = bs.FUsignal(fGHz, Cttl)
        self.Ctp = bs.FUsignal(fGHz, Ctpl)
        self.Cpt = bs.FUsignal(fGHz, Cptl)
        self.Cpp = bs.FUsignal(fGHz, Cppl)



        return self


    def show(self, **kwargs):
        """ show the propagation channel

        Parameters
        ----------

        typ   : 'm', 'l20' , 'r'
        cmap  : colormap
            default hot
        fontsize : int
            default 14

        """

        defaults = {'typ': 'm',
                    'cmap': plt.cm.hot,
                    'fontsize':14}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if 'fig' not in kwargs:
            kwargs['fig'] = plt.figure()


        ax1 = kwargs['fig'].add_subplot(221)
        fig, ax1 = self.Ctt.imshow(ax=ax1,**kwargs)
        ax1.set_xlabel('Frequency (GHz)',fontsize=kwargs['fontsize'])
        ax1.set_title(u'$C_{\\theta\\theta}$',fontsize=kwargs['fontsize'])

        ax2 = kwargs['fig'].add_subplot(222)
        fig, ax2 = self.Ctp.imshow(ax=ax2,**kwargs)
        ax2.set_xlabel('Frequency (GHz)',fontsize=kwargs['fontsize'])
        ax2.set_title(u'$C_{\\theta\phi}$',fontsize=kwargs['fontsize'])

        ax3 = kwargs['fig'].add_subplot(223)
        fig, ax3 = self.Cpt.imshow(ax=ax3,**kwargs)
        ax3.set_xlabel('Frequency (GHz)',fontsize=kwargs['fontsize'])
        ax3.set_title(u'$C_{\phi\\theta}$',fontsize=kwargs['fontsize'])

        ax4 = kwargs['fig'].add_subplot(224)
        fig, ax4 = self.Cpp.imshow(ax=ax4,**kwargs)
        ax4.set_xlabel('Frequency (GHz)',fontsize=kwargs['fontsize'])
        ax4.set_title(u'$C_{\phi\phi}$',fontsize=kwargs['fontsize'])

        return fig, (ax1, ax2, ax3, ax4)

    def check_reciprocity(self, C):
        """ check channel reciprocity

        Parameters
        ----------

        C : Ctilde


        Notes
        -----

        This is not properly implemented

        """
        issue=[]
        assert np.allclose(self.tauk, C.tauk)
        for r in range(self.nray):
            if not np.allclose(self.Ctt.y[r,:], C.Ctt.y[r,:]):
                issue.append(r)
        if len(issue) == 0:
            print "Channel is reciprocal"
        else: 
            print "WARNING Reciprocity issue WARNING"
            print len(issue),'/',self.nray, 'rays are not reciprocal,'
            print "rays number with an issue :",issue

        # assert np.allclose(self.tang,C.rang)
        # assert np.allclose(self.rang,C.tang)


    def energy(self,mode='mean',Friis=True,sumray=False):
        """ calculates energy on each channel

        Parameters
        ----------

        mode : string
            'mean'
        Friis: boolean
            True
        sumray: boolean
            False


        Returns
        -------

        ECtt  : Energy on co channel    tt
        ECpp  : Energy on co channel    pp
        ECtp  : Energy on co channel    tp
        ECpt  : Energy on co channel    pt


        See Also
        --------

        pylayers.signal.bsignal.FUsignal.energy

        Notes
        -----

        r x f+
          axis 0 : ray
          axis 1 : frequency


        """

        #
        #  r x f
        #  axis 0 : ray
        #  axis 1 : frequency
        #

        ECtt = self.Ctt.energy(axis=1,Friis=Friis,mode=mode)
        ECtp = self.Ctp.energy(axis=1,Friis=Friis,mode=mode)
        ECpt = self.Cpt.energy(axis=1,Friis=Friis,mode=mode)
        ECpp = self.Cpp.energy(axis=1,Friis=Friis,mode=mode)

        if sumray:
            ECtt = np.sum(ECtt,axis=0)
            ECtp = np.sum(ECtp,axis=0)
            ECpt = np.sum(ECpt,axis=0)
            ECpp = np.sum(ECpp,axis=0)

        return ECtt, ECpp, ECtp, ECpt

    def cut(self,threshold=0.99):
        """ cut rays from a energy threshold

        Parameters
        ----------

        threshold : float
            default 0.99

        """
        Ett, Epp, Etp, Ept = self.energy()
        Etot = Ett+Epp+Etp+Ept
        u = np.argsort(Etot)
        cumE = np.cumsum(Etot)/sum(Etot)
        v = np.where(cumE<threshold)[0]

        self.tauk = self.tauk[v]
        self.tang = self.tang[v,:]
        self.rang = self.rang[v,:]
        self.Ctt.y = self.Ctt.y[v,:]
        self.Cpp.y = self.Cpp.y[v,:]
        self.Ctp.y = self.Ctp.y[v,:]
        self.Cpt.y = self.Cpt.y[v,:]

    def sort(self,typ='tauk'):
        """ sort Ctilde with respect to typ (default tauk)

        Parameters
        ----------

        typ  : string
            sort w.r.t
                'tauk'   : delay (default)
                'att'    : theta Tx
                'atp'    : phi Tx
                'art'    : theta Rx
                'arp'    : phi Rx
                'energy' : energy

        """
        if typ == 'tauk':
            u = np.argsort(self.tauk)
        if typ == 'att':
            u = np.argsort(self.tang[:, 0])
        if typ == 'atp':
            u = np.argsort(self.tang[:, 1])
        if typ == 'art':
            u = np.argsort(self.rang[:, 0])
        if typ == 'arp':
            u = np.argsort(self.rang[:, 1])
        if typ == 'energy':
            Ett, Epp, Etp, Ept = self.energy()
            Etot = Ett+Epp+Etp+Ept
            u = np.argsort(Etot)
        self.tauk = self.tauk[u]
        self.tang = self.tang[u,:]
        self.rang = self.rang[u,:]

        self.Ctt.y = self.Ctt.y[u,:]
        self.Cpp.y = self.Cpp.y[u,:]
        self.Ctp.y = self.Ctp.y[u,:]
        self.Cpt.y = self.Cpt.y[u,:]

    def prop2tran(self,a=[],b=[],Friis=True,debug=True):
        r""" transform propagation channel into transmission channel

        Parameters
        ----------

        a : antenna or array a
        b : antenna or array b

        Ta : np.array(3x3)
           unitary matrice for antenna orientation
        Tb : np.array(3x3)
           unitary matrice for antenna orientation
        Friis : boolean
            if True scale with :math:`-j\frac{\lambda}{f}`
        debug : boolean
            if True the antenna gain for each ray is stored

        Returns
        -------

        H : Tchannel(bs.FUsignal)


        """
        freq  = self.fGHz
        nfreq = self.nfreq
        nray  = self.nray
        sh = np.shape(self.Ctt.y)

        # select default antennas
        # omni polar theta 't' <=> vertical polarization
        #
        if a ==[]:
            a = ant.Antenna('Omni',param={'pol':'t','GmaxdB':0},fGHz=self.fGHz)
        if b ==[]:
            b = ant.Antenna('Omni',param={'pol':'t','GmaxdB':0},fGHz=self.fGHz)

        a.eval(th=self.tangl[:, 0], ph=self.tangl[:, 1], grid=False)
        Fat = bs.FUsignal(a.fGHz, a.Ft)
        Fap = bs.FUsignal(a.fGHz, a.Fp)
        b.eval(th=self.rangl[:, 0], ph=self.rangl[:, 1], grid=False)
        Fbt = bs.FUsignal(b.fGHz, b.Ft)
        Fbp = bs.FUsignal(b.fGHz, b.Fp)

        # Cg2cl should be applied here
        #

        #
        #  C  = 2 x 2 x r x f
        #  Ctt : r x f
        #  Fa = 2 x r x f
        #  Fb = 2 x r x f
        #
        #  (r x f ) (r x Nt x f )
        t1 = self.Ctt * Fat + self.Ctp * Fap
        t2 = self.Cpt * Fat + self.Cpp * Fap

        # depending on siso or mimo case
        if len(t1.y.shape)==3:
            T1 = t1.y[:,None,:,:]
            T2 = t2.y[:,None,:,:]
        else:
            T1 = t1.y[:,None,None,:]
            T2 = t2.y[:,None,None,:]
        if len(Fbt.y.shape)==3:
            FBt = Fbt.y[:,:,None,:]
            FBp = Fbp.y[:,:,None,:]
        else:
            FBt = Fbt.y[:,None,None,:]
            FBp = Fbp.y[:,None,None,:]

        # determine the common interval on frequency axis
        if np.sum(t1.x!=Fbt.x)>0:
            t1x_int = (np.round(t1.x*100)).astype(int)
            Fbtx_int = (np.round(Fbt.x*100)).astype(int)
            inter = np.intersect1d(t1x_int,Fbtx_int)
            ut = np.in1d(t1x_int,inter)
            uf = np.in1d(Fbtx_int,inter)
        else:
            ut = np.arange(len(t1.x))
            uf = np.arange(len(Fbt.x))
        assert(len(t1.x[ut])==len(Fbt.x[uf])),"problem in common index plage calculation"

        alpha1 = np.einsum('ljkm,lkim->ljim',FBt[...,uf],T1[...,ut])
        alpha2 = np.einsum('ljkm,lkim->ljim',FBp[...,uf],T2[...,ut])

        #alpha = t1 * Fbt + t2 * Fbp
        # Nd x Nr x Nt x Nf
        alpha = alpha1 + alpha2

        self.fGHz = t1.x[ut]

        H = Tchannel(x = self.fGHz,
                     y = alpha,
                     tau = self.tauk,
                     dod = self.tang,
                     doa = self.rang)

        if debug :
            H.alpha=alpha
            H.Fat=Fat.y
            H.Fap=Fap.y
            H.Fbt=Fbt.y
            H.Fbp=Fbp.y
            H.Gat=10*np.log10(np.sum(Fat.y*np.conj(Fat.y),axis=1)/len(Fat.x))
            H.Gap=10*np.log10(np.sum(Fap.y*np.conj(Fap.y),axis=1)/len(Fap.x))
            H.Gbt=10*np.log10(np.sum(Fbt.y*np.conj(Fbt.y),axis=1)/len(Fbt.x))
            H.Gbp=10*np.log10(np.sum(Fbp.y*np.conj(Fbp.y),axis=1)/len(Fbp.x))

        if Friis:
            H.applyFriis()

        # average w.r.t frequency
        Nf   = H.y.shape[-1]
        H.ak = np.real(np.sqrt(np.sum(H.y * np.conj(H.y)/Nf, axis=1)))
        H.tk = H.taud

        return(H)

if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
