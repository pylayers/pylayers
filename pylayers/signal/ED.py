#-*- coding:Utf-8 -*-
r""""
.. currentmodule:: pylayers.signal.ED

.. autosummary::
    :toctree: generated

Energy Detector
===============

.. autosummary::
    :toctree: generated/

    ED.__init__
    ED.__repr__
    ED.apply
    ED.moment
    ED.pdf
    ED.errprob
    ED.roc
    ED.plot

References
----------

    Cite the relevant literature, e.g. [1]_.  You may also cite these
    references in the notes section above.

.. [1] O. McNoleg, "The integration of GIS, remote sensing,
       expert systems and adaptive co-kriging for environmental habitat
       modelling of the Highland Haggis using object-oriented, fuzzy-logic
       and neural-network techniques," Computers & Geosciences, vol. 22,
       pp. 585-588, 1996.

"""
from scipy.signal import *
import scipy.special as spe
import pandas as pd
import scipy.stats as st
from pylab import *
from pylayers.signal.DF import *
import pylayers.signal.bsignal as bs
from pylayers.util.project import *
import matplotlib.pyplot as plt
import doctest

class ED(PyLayers):
    """
    Energy Detector Class

    This class implements an energy detector receiver


    """
    def __init__(self,**kwargs):
        """ Constructor

        Parameters
        ----------

        fs : sampling frequency
        fc : center frequency
        beta : scale factor :math:`$\beta$`
        Tns  : time integration in nanoseconds
        BGHz : Bandwidth in GHz
        pfa  : false alarm probability
        wt   : filter B ripple
        R    : Resistance for Tension/Power conversion 
        NF   : Noise factor in dB (default 0dB)
        gpass: gain in the pass band
        gstop : gain in the stop band

        Notes
        -----

        An energy detector is determined by 2 filter a bandpass elliptic filter 
        :math:`H_B(f)` and a lowpass filter which corresponds to the time 
        integration over duration :math:`T` :math:`H_T(f)`

        """

        defaults = {'fsGHz':50,
                    'fcGHz':4.4928,
                    'Tns':64,
                    'pfa':0.01,
                    'beta':1,
                    'wt':0.01,
                    'gpass':0.5,
                    'gstop':30,
                    'NF':0,
                    'R':50}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.fsGHz = kwargs['fsGHz']
        self.fcGHz = kwargs['fcGHz']
        self.R = kwargs['R']
        BGHz = kwargs['BGHz']
        Tns = kwargs['Tns']
        wt  = kwargs['wt']
        #
        # Filter B : Input Band Pass Filter
        #
        # Elliptic filter
        #
        
        self.ts  = 1./self.fsGHz
        #fN  = self.fsGHz/2.
        #w1  = (self.fcGHz-BGHz/2.)/fN
        #w2  = (self.fcGHz+BGHz/2.)/fN
        #wp  = [w1,w2]
        #ws  = [w1-wt,w2+wt]

        #self.filterB  = DF()
        #self.filterB.ellip_bp(wp,ws,kwargs['gpass'],kwargs['gstop'])
        #self.filterB.butter(order=6,w=wp,typ='bandpass')

        #
        # Filter T : Time Averaging Filter (Integrator)  (Tns => N samples)
        #

        self.N  = np.ceil(Tns/self.ts)
        #b  = (Tns*1e-9/self.N)*ones(self.N)
        #self.NfilterT = len(b)
        #a  = array([1])
        #self.filterT = DF(b,a)

        #
        #
        #

        self.BGHz  = BGHz
        self.Tns   = Tns
        self.beta  = kwargs['beta']
        self.pfa   = kwargs['pfa']
        self.NF    = kwargs['NF']

        self.update()

    def update(self):
        """ updating 
        
        """
        #
        # calculates ED moments, order and scale
        #
        # self.moment(typ='struc')
        # kB = 1,3806488  10-23
        kB = 1.3806488e-23
        # TB : Noise temperature
        TB  = 290*10**(self.NF/10.)
        self.N0  = kB*TB
        self.order = 2*self.BGHz*self.Tns
        self.muy   = self.order*self.N0
        self.vary  = 4*self.BGHz*self.Tns*self.N0**2
        self.scale = np.sqrt(self.vary/(2*self.order))

    def __repr__(self):
        self.update()
        st = 'Energy Detector'+'\n'
        st = st+'---------------'+'\n\n'
        st = st+'Center frequency : '+str(self.fcGHz)+' GHz\n'
        st = st+'Bandwidth : '+str(self.BGHz)+' GHz \n'
        st = st+'Integration Time : '+str(self.Tns)+' ns \n'
        st = st+'Resistance : '+str(self.R)+' Ohms\n'
        st = st+'Mean (y) : '+str(self.muy/1e-9)+' nJ \n'
        st = st+'Std (y) : '+str(np.sqrt(self.vary)/1e-9)+' nJ \n'
        st = st+'Order (df) : '+str(self.order)+'  \n'
        st = st+'Scale : '+str(self.scale)+'  \n'
        st = st+'2BT : '+str(2*self.BGHz*self.Tns)+' \n'

        return(st)

    def apply(self,s,typ='ideal'):
        r""" apply signal 

        Parameters
        ----------

        s : input signal
        typ : 'real' | 'ideal'

        Returns
        -------

        y : :math:`F_2\{\beta^2 F_1\{x\}^2\}`

            :math:`F{x}` means a filtering of x with filter :math:`F`

        """
        if typ=='real':
            # B filtering
            self.xf1 = self.filterB.filter(s.y)

        if typ=='ideal':
            # B filtering  --> FH signal
            S = s.fft()
            self.gate = np.zeros(len(S.x))
            self.gate[0]=1
            u1 = np.where((S.x>(self.fcGHz-self.BGHz/2.)) &
                          (S.x<(self.fcGHz+self.BGHz/2.)))[0]
            u2 = np.where((S.x>(self.fsGHz-self.fcGHz-self.BGHz/2.)) &
                          (S.x<(self.fsGHz-self.fcGHz+self.BGHz/2.)))[0]
            self.gate[u1] = 1
            self.gate[u2] = 1
            S.y  = S.y*self.gate
            self.xf1 = S.ifft().y

        # quadrator and scaling
        self.xq  = ((self.beta)**2)*(self.xf1*self.xf1)/self.R
        # T filtering
        #xfq = self.filterT.filter(self.xq)
        xfq = pd.rolling_mean(self.xq,self.N)*self.N*1e-9*self.ts
        self.y = bs.TUsignal(s.x[self.N:],xfq[self.N:])

        return(self.y)

    def moment(self,typ='struc'):
        r""" evaluate mean and variance of output signal y

        Parameters
        ----------

        typ : string
          'struc' structural evaluation of moments
          'emp' empirical evaluation of moments

        Notes
        -----

        :math:`R_y(\tau)=2R_x^2(\tau)+R_x^2(0)`

        This function updates the order and scale or the Energy Detector

        """
        self.update()
        # filterB bandpass filter
        self.filterB.freqz(display=False)
        # filterT time integration (=low pass filter)
        self.filterT.freqz(display=False)
        w = bs.Noise(fsGHz=self.fsGHz)

        if typ=='struc':
            # ACF of filter B
            self.RB = self.filterB.wnxcorr(var=w.var)
            # F.de Coulom T&T signal eq 8.148
            # ACF of xq
            self.Rx = self.RB*self.RB*2
            # Power of xq
            RB0= self.RB.y[0]
            #
            self.Rx.y = self.Rx.y+RB0**2

            Phix = np.real(fft.fft(self.Rx.y))
            HH   = self.filterT.H.symH(0)
            HH.y = HH.y/(self.Tns*1e-9)
            Phiy = Phix * np.real(HH.y * np.conj(HH.y))
            self.muy = self.RB.max()[0]*self.beta*self.beta
            self.vary = mean(Phiy)

        if typ=='emp':
            self.muy  = self.y.y[self.NfilterT:].mean()
            self.vary = self.y.y[self.NfilterT:].var()

        self.order = (2*(self.muy)**2)/self.vary
        self.scale = np.sqrt(self.vary/(2*self.order))

    def pdf(self,E,Ei=0):
        """ calculates decision problem densities under H0 and H1

        Parameters
        ----------

        E  : float
            Legitimate Signal Energy
        Ei : float
            Interferer Signal Energy

        See Also
        --------

        scipy.stats.chi2
        scipy.stats.ncx2

        """

        self.E = E
        self.Ei = Ei

        # H0 : Noise only (no interferer)
        # central chi square
        if Ei==0:
            self.pH0 = st.chi2(self.order,
                               scale=self.scale)
        else:
        # H0 : Noise + interference
        # non central chi square
            self.pH0 = st.ncx2(self.order,
                               scale=self.scale,
                               nc=Ei/(self.scale))

        # H1 : Signal + interference + Noise
        # non central chi square
        self.pH1 = st.ncx2(self.order,
                           scale=self.scale,
                           nc=(E+Ei)/(self.scale))

    def errprob(self,E,Ei=0,thresh=2e-9,typ='thresh'):
        """ calculates error probability

        Parameters
        ----------

        E  : float
            Legitimate Signal Energy
        Ei : float
            Interferer Signal Energy
        thresh : float
        typ : string
            'thresh'
            'scale': for IEEE.802.11.6 standard
            'both' : both kind of decision

        Returns
        -------

        errp : error probability

        Notes
        -----

        See Also
        --------

        scipy.stats.ncf
        ED.pdf

        """

        self.pdf(E,Ei=Ei)

        if typ=='thresh':
            self.p10 = self.pH0.sf(thresh)
            self.p01 = 1-self.pH1.sf(thresh)
            errp = 0.5*(self.p10+self.p01)

        if typ=='scale':
            ncf = st.ncf(dfn = self.order,
                         dfd = self.order,
                         scale = 1,
                         nc = self.E/self.scale)
            errp = 1 - ncf.sf(1)

        if typ=='both':
            self.p10 = self.pH0.sf(thresh)
            self.p01 = 1-self.pH1.sf(thresh)
            pe1 = 0.5*(self.p10+self.p01)
            ncf = st.ncf(dfn=self.order,
                         dfd=self.order,
                         scale=1,
                         nc=self.E/self.scale)
            pe2= 1 - ncf.sf(1)
            errp =(pe1,pe2)

        return(errp)

    def plot(self,E):
        """ plot decision problem densities

        Parameters
        ----------

        E : float
            value of signal energy

        """

        self.pdf(E)
        e  = np.linspace(0,max(self.muy*5,2*self.E),200)
        p0 = self.pH0.pdf(e)
        p1 = self.pH1.pdf(e)
        #plt.plot(10*np.log10(E),p0,'b',linewidth=4)
        #plt.plot(10*np.log10(E),p1,'r',linewidth=4)
        plt.plot(e,p0,'b',linewidth=4)
        plt.plot(e,p1,'r',linewidth=4)

    def roc(self,E,typ='thresh'):
        """ determines receiver operating curve

        Parameters
        ----------

        E   : float
            Energy (J)
        typ : string
            'thresh' | 'both' | 'scale'


        Returns
        -------

        terr : error probability
        tpfa : false alarm probability
        tpd  : detection probability


        Notes
        -----

        Iterate on energy from 0 to 10*E

        Examples
        --------

        ... plot:
            include-source:

            >>> from pylayers.signal.ED import *
            >>> EnDet = ED(BGHz=0.5,Tns=20)


        """

        thresh = np.linspace(10*E,0,300)
        terr = []
        tpfa = []
        tpd =[]

        for t in thresh:
            perr = self.errprob(E,thresh=t,typ=typ)
            terr.append(perr)
            tpfa.append(self.p10)
            tpd.append(1-self.p01)

        terr = np.array(terr)
        tpfa = np.array(tpfa)
        tpd = np.array(tpd)

        return(terr,tpfa,tpd)



if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
