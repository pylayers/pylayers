from scipy.signal import *
import scipy.special as spe
import scipy.stats as st
from pylab import *
from pylayers.signal.DF import *
import pylayers.signal.bsignal as bs
from pylayers.util.project import *
import matplotlib.pyplot as plt
import doctest

r"""
.. currentmodule:: pylayers.signal.ED

.. autosummary::
    :toctree: generated

"""
class ED(PyLayers):
    """
    Energy Detector Class

    This class implements an Energy Detector receiver

    """
    def __init__(self,**kwargs):
        """

        fs : sampling frequency
        fc : center frequency
        beta : scale factor
        Tns  : time integration
        BGHz : Bandwidth
        pfa : false alarm probability
        wt :
        gpass:
        gstop :

        """

        defaults = {'fsGHz':100,
                    'fcGHz':4,
                    'BGHz':3,
                    'Tns':10,
                    'pfa':0.01,
                    'beta':1,
                    'wt':0.01,
                    'gpass':0.5,
                    'gstop':30}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        fs = kwargs['fsGHz']
        fc = kwargs['fcGHz']
        BGHz = kwargs['BGHz']
        Tns = kwargs['Tns']
        wt  = kwargs['wt']
        #
        # Filter B : Input Band Pass Filter
        #
        ts  = 1./fs
        fN  = fs/2.
        w1  = (fc-BGHz/2.)/fN
        w2  = (fc+BGHz/2.)/fN
        wp  = [w1,w2]
        ws  = [w1-wt,w2+wt]


        self.filterB  = DF()
        self.filterB.ellip_bp(wp,ws,kwargs['gpass'],kwargs['gstop'])
        #self.filterB.butter(order=6,w=wp,typ='bandpass')

        #
        # Filter T : Time Averaging Filter (Integrator)  (Tns => N samples)
        #

        N  = np.ceil(Tns/ts)
        b  = (1./N)*ones(N)
        a  = array([1])
        self.filterT = DF(b,a)

        #
        #
        #
        self.beta  = kwargs['beta']
        self.pfa   = kwargs['pfa']

        #
        # calculates ED moments, order and scale
        #

        self.moment(typ='struc')

    def apply(self,x):
        r"""

        Parameters
        ----------

        x : input signal

        Returns
        -------

        y : $F2{\beta^2 F1{x}^2}$

        F{x} means a filtering of x with filter F

        """
        self.xf1 = self.filterB.filter(x.y)
        self.xq  = ((self.beta)**2)*(self.xf1*self.xf1)
        y        = bs.TUsignal(x.x,self.filterT.filter(self.xq))
        return(y)

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

        This functions updates the order and scale or the ED

        """
        self.filterB.freqz(display=False)
        self.filterT.freqz(display=False)
        w = bs.Noise()

        if typ=='struc':
            self.RB = self.filterB.wnxcorr(var=w.var)
            # F.de Coulom T&T signal eq 8.148
            self.Rx = self.RB*self.RB*2
            RB0= self.RB.y[0]

            self.Rx.y = self.Rx.y+RB0**2

            Phix = np.real(fft.fft(self.Rx.y))
            HH   = self.filterT.H.symH(0)
            Phiy = Phix * np.real(HH.y * np.conj(HH.y))
            self.muy = self.RB.max()[0]*self.beta*self.beta
            self.vary = mean(Phiy)

        if typ=='emp':
            self.muy = 1
            self.vary = 1

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

        if Ei==0:
            self.pH0 = st.chi2(self.order,
                               scale=self.scale)
        else:
            self.pH0 = st.ncx2(self.order,
                               scale=self.scale,
                               nc=Ei/(self.scale))

        self.pH1 = st.ncx2(self.order,
                           scale=self.scale,
                           nc=(E+Ei)/(self.scale))

    def errprob(self,E,Ei=0,thresh=0,typ='thresh'):
        """ calculates error probability

        Parameters
        ----------

        E  : float
        Ei : float
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
            value of energy

        """

        self.pdf(E)
        e = np.linspace(0,max(self.muy*5,2*self.E),200)
        p0 = self.pH0.pdf(e)
        p1 = self.pH1.pdf(e)
        #plt.plot(10*np.log10(E),p0,'b',linewidth=4)
        #plt.plot(10*np.log10(E),p1,'r',linewidth=4)
        plt.plot(e,p0,'b',linewidth=4)
        plt.plot(e,p1,'r',linewidth=4)

    def roc(self,E,typ='thresh'):
        """ determine receiver operating curve

        Parameters
        ----------

        E   : float
            Energy
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
