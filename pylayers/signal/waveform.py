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
from pylayers.signal   import bsignal as bs
from pylayers.util     import easygui
from pylayers.measures import mesuwb

class Waveform:
    """
        Waveform Class
    """
    def __init__(self,parameters=[]):
        """
           wavetype   =   'generic' =   'mbofdm' =   'fromfile'
           st  : time domain
           sf  : frequency domain
           sfg : frequency domain integrated
        """
        param= {'type' : 'generic',
                'band': 0.499,
                'fc': 4.493,
                'fe': 100,
                'thresh': 3,
                'tw': 30}
        if parameters==[]:
            self.parameters = param
        else:
            self.parameters = parameters
        self.eval()

    def eval(self):
        """ evaluate waveform

            The lambda/4*pi factor which is necessary to get the proper budget 
            link ( from the Friis formula) is introduced in this function.
        """

        if self.parameters['type']  == 'generic':
            [st,sf]=self.ip_generic()
        elif self.parameters['type']  == 'mbofdm':
            [st,sf]=self.mbofdm()
        elif self.parameters['type'] == 'W1compensate':
            [st,sf]=self.fromfile()
        elif self.parameters['type'] == 'W1offset':
            [st,sf]=self.fromfile2()
        else:
            logging.critical('waveform type not recognized, check your config \
                             file')

        self.st       = st
        self.sf       = sf
        self.f        = self.sf.x
        ygamma        = -1j*0.3/(4*np.pi*self.f)
        self.gamm     = bs.FUsignal(self.f,ygamma)
        self.sfg      = self.sf*self.gamm
        self.sfgh     = self.sfg.symH(0)
        self.stgh     = self.sfgh.ifft(1)

    def info(self):
        """  display information about waveform

        Results
        -------

        >>> from pylayers.signal.waveform import *
        >>> param= {'type' : 'generic',\
            'band': 0.499,\
            'fc': 4.493,\
            'fe': 100,\
            'thresh': 3,\
            'tw': 30}
        >>> w = Waveform(param)
        >>> # W = {}
        >>> # W['t']=w.st.x
        >>> # W['p']=w.st.y
        >>> w.show()
        >>> plt.show()         


        """
        if self.parameters['type']=='generic':
            for k in self.parameters.keys():
                print k , " : ",self.parameters[k]
        else:
            print "type:",self.parameters['type']

    def show2(self,Tp=1000):
        """ show2

        Parameters
        ----------
        Tp : float
        
        """
        plt.subplot(211)
        self.st.plot()
        plt.subplot(212)
        psd = self.st.psd(Tp,50)
        plt.title('Tp = '+str(Tp))
        psd.plotdB(mask=True)

    def ip_generic(self):
        """    Create an Energy normalized Gaussian impulse (Usignal)
            ip_generic(self,parameters)
        

        """
        Tw     = self.parameters['tw']
        fc     = self.parameters['fc']
        band   = self.parameters['band']
        thresh = self.parameters['thresh']
        fe     = self.parameters['fe']
        te     = 1.0/fe

        self.parameters['te'] = te
        Np     = fe*Tw
        self.parameters['Np']=Np
        x      = np.linspace(-0.5*Tw+te/2,0.5*Tw+te/2,Np,endpoint=False)
        #x     = arange(-Tw,Tw,te)

        w   = bs.EnImpulse(x,fc,band,thresh,fe)
        #W = w.ft()
        W   = w.ftshift()
        return (w,W)

    def fromfile(self):
        """
        get the measurement waveform from WHERE1 measurement campaign

        This function is not yet generic

        >>> from pylayers.signal.waveform import *
        >>> wav = Waveform({'type':'file'})
        >>> wav.show()

        """
        M = mesuwb.UWBMesure(1,1)
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
        >>> wav = Waveform({'type':'file'})
        >>> wav.show()

        """
        M = mesuwb.UWBMesure(1,1)
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
        wparam =  {}
        for k in range(len(par)):
            key = par[k][0]
            val = par[k][1]
            if key == "band":
                wparam[key] = float(val)
            if key == "fc":
                wparam[key] = float(val)
            if key == "fe":
                wparam[key] = float(val)
            if key == "thresh":
                wparam[key] = float(val)
            if key == "tw":
                wparam[key] = float(val)
            if key == "type":
                wparam[key] = val
 
        self.parameters = wparam
        self.eval()


    def gui(self):
        """
        Get the Waveform parameter
        """
        if self.parameters['type'] == 'generic':
            self.st.plot()
            show()
            wavegui = multenterbox('','Waveform Parameter',
            ('Tw (ns) integer value',
             'fc (GHz)',
             'band (GHz)',
             'thresh (dB)',
             'fe (GHz) integer value'),
            ( self.parameters['Tw'] ,
            self.parameters['fc'] ,
            self.parameters['band'] ,
            self.parameters['tresh'],
            self.parameters['fe']
            ))

            self.parameters['Tw']    = eval(wavegui[0])
            self.parameters['fc']    = eval(wavegui[1])
            self.parameters['band']  = eval(wavegui[2])
            self.parameters['tresh'] = eval(wavegui[3])
            self.parameters['fe']    = eval(wavegui[4])

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
        for pk in self.parameters.keys():
            val   = self.parameters[pk]
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
    doctest.testmod()
