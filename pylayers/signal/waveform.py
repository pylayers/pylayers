# -*- coding:Utf-8 -*-
import doctest
import os
import logging
import pdb
import numpy as np
import scipy as sp
import scipy.io as io
import scipy.signal as si
import scipy.linalg as la
import ConfigParser
import matplotlib.pylab as plt
import pylayers.signal.bsignal as bs
from pylayers.util     import easygui
from pylayers.measures import mesuwb

class Waveform(dict):
    """

    Attributes
    ----------

    st  : time domain
    sf  : frequency domain
    sfg : frequency domain integrated

    Methods
    -------

    eval
    show2
    ip_generic
    fromfile
    fromfile2
    read
    gui
    show

    """
    def __init__(self,**kwargs):
        """

        Parameters
        ----------

        'typ' : string
            'generic',
         'bandGHz': float
            0.499
         'fcGHz': float
            4.493
         'fsGHz': float
            100,
         'threshdB':
              3,
         'twns': float
            30

        typ  :  'generic','W1compensate','W1offset'

        """
        defaults = {'typ':'generic',
                'bandGHz': 0.499,
                'fcGHz': 4.493,
                'feGHz': 100,
                'threshdB': 3,
                'twns': 30}

        for key, value in defaults.items():
            if key not in kwargs:
                self[key] = value
            else:
                self[key] = kwargs[key]
        self.eval()

    def eval(self):
        u""" evaluate waveform

        The :math:`\lambda/4*\pi` factor which is necessary to get the proper budget
        link ( from the Friis formula) is introduced in this function.

        """

        if self['typ']  == 'generic':
            [st,sf]=self.ip_generic()
        #elif self['typ']  == 'mbofdm':
        #    [st,sf]=self.mbofdm()
        elif self['typ'] == 'W1compensate':
            [st,sf]=self.fromfile()
        elif self['typ'] == 'W1offset':
            [st,sf]=self.fromfile2()
        else:
            logging.critical('waveform typ not recognized, check your config \
                             file')

        self.st       = st
        self.sf       = sf
        self.f        = self.sf.x

        ygamma        = -1j*0.3/(4*np.pi*self.f)
        self.gamm     = bs.FUsignal(x=self.f,y=ygamma)
        self.sfg      = self.sf*self.gamm
        self.sfgh     = self.sfg.symH(0)
        self.stgh     = self.sfgh.ifft(1)

    def info(self):
        """  display information about waveform

        Results
        -------

        >>> from pylayers.signal.waveform import *
        >>> w = Waveform(typ='generic',bandGHz=0.499,fcGHz=4.49,feGHz=100,threshdB=3,twns=30)
        >>> w.show()
        >>> plt.show()


        """
        if self['typ']=='generic':
            for k in self.keys():
                print k , " : ",self[k]
        else:
            print "typ:",self['typ']

    def show2(self,Tpns=1000):
        """ show2

        Parameters
        ----------

        Tpns : float

        """
        plt.subplot(211)
        self.st.plot()
        plt.subplot(212)
        psd = self.st.psd(Tp,50)
        plt.title('Tp = '+str(Tp))
        psd.plotdB(mask=True)

    def ip_generic(self):
        """   Create an Energy normalized Gaussian impulse (Usignal)

        ip_generic(self,parameters)


        """
        Tw = self['twns']
        fcGHz = self['fcGHz']
        band = self['bandGHz']
        thresh = self['threshdB']
        feGHz = self['feGHz']
        te = 1.0/feGHz

        self['te'] = te
        Np     = feGHz*Tw
        self['Np']=Np
        #x      = np.linspace(-0.5*Tw+te/2,0.5*Tw+te/2,Np,endpoint=False)
        #x     = arange(-Tw,Tw,te)

        w = bs.EnImpulse(fcGHz=fcGHz,WGHz=band,threshdB=thresh,feGHz=feGHz)
        #W = w.ft()
        W = w.ftshift()
        return (w,W)

    def ref156(self):
        """ reference pulse of IEEE 802.15.6 UWB standard
        """
        Tw     = self['twns']
        fc     = self['fcGHz']
        band   = self['bandGHz']
        thresh = self['threshdB']
        fe     = self['feGHz']
        te     = 1./fe
        beta = 0.5
        Tns  = 1./0.4992
        x    =  np.linspace(-0.5*Tw+te/2,0.5*Tw+te/2,Np,endpoint=False)
        z    = x/T
        t1  = np.sin(np.pi*(1-beta)*z)
        t2  = np.cos(np.pi*(1+beta)*z)
        t3  = (np.pi*z)*(1-(4*beta*z)**2)
        y   = (t1+4*beta*z*t2)/t3


    def fromfile(self):
        """
        get the measurement waveform from WHERE1 measurement campaign

        This function is not yet generic

        >>> from pylayers.signal.waveform import *
        >>> wav = Waveform(typ='W1compensate')
        >>> wav.show()

        """

        M = mesuwb.UWBMeasure(1,h=1)
        w = bs.TUsignal()

        ts = M.RAW_DATA.timetx[0]
        tns = ts*1e9
        te = tns[1]-tns[0]

        y  = M.RAW_DATA.tx[0]

        # find peak position  u is the index of the peak
        # yap :after peak
        # ybp : before peak
        # yzp : zero padding
        maxy = max(y)
        u = np.where(y ==maxy)[0][0]
        yap = y[u:]
        ybp = y[0:u]

        yzp = np.zeros(len(yap)-len(ybp)-1)

        tnsp = np.arange(0,tns[-1]-tns[u]+0.5*te,te)
        tnsm = np.arange(-(tns[-1]-tns[u]),0,te)

        y = np.hstack((yzp,np.hstack((ybp,yap))))
        tns = np.hstack((tnsm,tnsp))

        #
        # Warning (check if 1/sqrt(30) is not applied elsewhere
        #
        w.x = tns
        w.y = y*(1/np.sqrt(30))

        #  w : TUsignal
        #  W : FUsignal (Hermitian redundancy removed)

        W   = w.ftshift()
        return (w,W)

    def fromfile2(self):
        """
        get the measurement waveform from WHERE1 measurement campaign

        This function is not yet generic

        >>> from pylayers.signal.waveform import *
        >>> wav = Waveform(typ='W1offset')
        >>> wav.show()

        """
        M = mesuwb.UWBMeasure(1,1)
        w = bs.TUsignal()

        ts = M.RAW_DATA.timetx[0]
        tns = ts*1e9
        te = tns[1]-tns[0]

        y  = M.RAW_DATA.tx[0]

        # find peak position  u is the index of the peak
        # yap :after peak
        # ybp : before peak
        # yzp : zero padding
#        maxy = max(y)
#        u = np.where(y ==maxy)[0][0]
#        yap = y[u:]
#        ybp = y[0:u]

        yzp = np.zeros(len(y)-1)

#        tnsp = np.arange(0,tns[-1]-tns[u]+0.5*te,te)
#        tnsm = np.arange(-(tns[-1]-tns[u]),0,te)
        N=len(ts)-1
        tnsm = np.linspace(-tns[-1],-te,N)
        y = np.hstack((yzp,y))
        tns = np.hstack((tnsm,tns))

        #
        # Warning (check if 1/sqrt(30) is not applied elsewhere
        #
        w.x = tns
        w.y = y*(1/np.sqrt(30))

        #  w : TUsignal
        #  W : FUsignal (Hermitian redundancy removed)

        W   = w.ftshift()
        return (w,W)


    def read(self,config):
        """
        Parameters
        ----------

            config : ConfigParser object

        Returns
        -------
            w      : waveform
        """

        par = config.items("waveform")
        for k in range(len(par)):
            key = par[k][0]
            val = par[k][1]
            if key == "bandGHz":
                self[key] = float(val)
            if key == "fcGHz":
                self[key] = float(val)
            if key == "feGHz":
                self[key] = float(val)
            if key == "threshdB":
                self[key] = float(val)
            if key == "twns":
                self[key] = float(val)
            if key == "typ":
                self[key] = val

        self.eval()


    def gui(self):
        """
        Get the Waveform parameter
        """
        if self['typ'] == 'generic':
            self.st.plot()
            show()
            wavegui = multenterbox('','Waveform Parameter',
            ('Tw (ns) integer value',
             'fc (GHz)',
             'band (GHz)',
             'thresh (dB)',
             'fe (GHz) integer value'),
            ( self['twns'] ,
            self['fcGHz'] ,
            self['bandGHz'] ,
            self['threshdB'],
            self['feGHz']
            ))

            self.parameters['Twns']    = eval(wavegui[0])
            self.parameters['fcGHz']    = eval(wavegui[1])
            self.parameters['bandGHz']  = eval(wavegui[2])
            self.parameters['threshdB'] = eval(wavegui[3])
            self.parameters['feGHz']    = eval(wavegui[4])

            [st,sf]       = self.ip_generic()
            self.st       = st
            self.sf       = sf
            st.plot()
            show()

    def show(self,fig=[]):
        """
        """
        # title construction
        if fig ==[]:
            fig = plt.figure()
        title =''
        for pk in self.keys():
            val   = self[pk]
            title = title + pk + ': '
            if type(val) != 'str':
                title = title + str(val) + ' '
        #plt.title(title)
        ax1 = fig.add_subplot(2,1,1)
        ax1.plot(self.st.x,self.st.y)
        plt.xlabel('time (ns)')
        plt.ylabel('level in linear scale')
        ax2 = fig.add_subplot(2,1,2)
        ax2.plot(self.sf.x,abs(self.sf.y))
        plt.xlabel('frequency (GHz)')
        plt.ylabel('level in linear scale')
        fig.suptitle(title)


if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
