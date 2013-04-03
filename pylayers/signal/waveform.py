# -*- coding:Utf-8 -*-
import doctest
import os
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
            link ( from the fris formula) is introduced in this function.
        """

        if self.parameters['type']  == 'generic':
            [st,sf]=self.ip_generic()
        if self.parameters['type']  == 'mbofdm':
            [st,sf]=self.mbofdm()
        if self.parameters['type'] == 'file':
            [st,sf]=self.fromfile()

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
        for k in self.parameters.keys():
            print k , " : ",self.parameters[k]

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
        get the measurement waveform
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
    
        M = mesuwb.UWBMesure(1,1)
        w   = bs.TUsignal()
        timeTx = M.RAW_DATA.timetx[0]*1e9
        timeI = timeTx[1]-timeTx[0]
        Tx = M.RAW_DATA.tx[0]
        u=np.where(Tx==max(Tx))[0][0]
    
        Tx3=Tx[u:]
        Tx2=Tx[0:u]
        Tx1=np.zeros(len(Tx3)-len(Tx2))
        
        t=np.arange(0,timeTx[-1]-timeTx[u]+0.5*timeI,timeI)
        tm=np.arange(-(timeTx[-1]-timeTx[u]),0,timeI)
    
        rTx=np.hstack((Tx1,np.hstack((Tx2,Tx3))))
        rt=np.hstack((tm,t))
        w.x=rt
        w.y=rTx*(1/np.sqrt(30))
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

    def show(self):
        """
        """
        title =''
        for st in self.parameters.keys():
            val   = self.parameters[st]
            title = title + st + ': '
            if type(val) != 'str':
                title = title + str(val) + ' '
        fig = plt.figure()
        plt.title(title)
        fig.add_subplot(211)
        self.st.plot()
        fig.add_subplot(212)
        plt.plot(self.sf.x,abs(self.sf.y))



if __name__ == "__main__":
    doctest.testmod()
