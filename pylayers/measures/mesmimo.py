from pylayers.signal.bsignal import *
from pylayers.util.project import *
from pylayers.gis.readvrml import *
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import scipy.linalg as la
#
# This class handles the data coming from the MIMO Channel Sounder IETR lab
#
class MIMO(object):
    def __init__(self,
                 _filename='',
                 rep = '',
                 Nf = 1601,
                 fminGHz = 1.8,
                 fmaxGHz =2.2,
                 calibration=True,
                 time=True,
                 Nz = 100,
                 Nt = 4,
                 Nr = 8):
        """

        Parameters
        ----------

        filename : string
        rep : string
        fminGHz : float
        fmaxGHz : float
        Nf
        calibration : Boolean
        Nz : int
            Number of Zeros
        nT : int
            (default = 1)

        Notes
        -----

        Data are placed in the directory mesdir + rep directory

        """

        self.freq = np.linspace(fminGHz,fmaxGHz,Nf)

        self.Nt = Nt
        self.Nr = Nr

        self.Nf = Nf

        if _filename <> '':
            self.filename = mesdir+rep+_filename
            # load file
            self.loadraw(self.Nt)
            if calibration:
                self.calibration()
                if time:
                    # reshaping for using ift (todo update ift for MDA !!)
                    Hcal = FUsignal(self.Hcal.x,np.reshape(self.Hcal.y,(Nt*Nr,Nf)))
                    hcal = Hcal.ift(Nz=Nz,ffts=1)
                    shh = hcal.y.shape
                    self.hcal = TUsignal(hcal.x,np.reshape(hcal.y,(Nr,Nt,shh[-1])))


    def __repr__(self):
        str = 'MIMO Object'+'\n'
        str = 'Nr : '+str(self.Nr)+ '\n'
        str = 'Nt : '+str(self.Nr)+ '\n'
        str = 'Nf : '+str(self.Nf)+ '\n'
        return(str)

    def __sub__(self,m):
        N = MIMO()
        N.freq = self.freq
        N.Nt = self.Nt
        N.Nr = self.Nr
        N.Hcal = self.Hcal - m.Hcal
        return(N)

    def loadraw(self,nT=1):
        """ load a MIMO NrxNt raw data sounder file

            The sounder output file is a 2 columns ASCII csv file
            Module (dB) ;  Angle (Degree)

        """

        fd  = open(self.filename)
        lis = fd.readlines()
        fd.close()
        module  = []
        phasedeg   = []

        for l in lis:
            l.replace('\r\n','')
            g = l.split(';')
            module.append(float(g[0]))
            phasedeg.append(float(g[1]))

        m   = np.array(module)
        phi = np.array(phasedeg)*np.pi/180.
        m   = m.reshape(self.Nr*self.Nt,self.Nf)
        phi = phi.reshape(self.Nr*self.Nt,self.Nf)
        y   = 10**(m/20)*np.exp(1j*phi)
        #
        # Nr x Nt x Nf    (8x4x1601)
        #
        y   = y.reshape(self.Nr,self.Nt,self.Nf)

        self.H = FUsignal(self.freq,y)

    def calibration(self):
        """ Apply calibration files

        """
        for iR in range(self.Nr):
            for iT in range(self.Nt):
                filename = 'Calib'+str(iT+1)+'x'+str(iR+1)+'.txt'
                C = MIMO(filename,'/calibration/',calibration=False,Nt=self.Nt)
                try:
                    #tc = np.vstack((tc,C.H.y[iR*4+iT,:]))
                    tc = np.vstack((tc,C.H.y[iR,iT,:]))
                except:
                    #tc = C.H.y[iR*4+iT,:]
                    tc = C.H.y[iR,iT,:]

        tc = tc.reshape(self.Nr,self.Nt,self.Nf)

        self.C = FUsignal(C.freq,tc)

        self.Hcal = self.H/self.C

        del self.H
        del self.C

    def capacity(self,trho=np.arange(-40,120,1)):
        """ calculates channel Shannon capacity

        thrho : np.array (SNR)

            log_2(det(I+(Et/N0Nt)HH^{H})
        """

        tC = np.zeros((len(trho),len(self.freq)))

        for ir in range(len(trho)):
            rho = 10**(trho[ir]/10)
            for ik in range(len(self.freq)):
                #pdb.set_trace()
                H   = self.Hcal.y[:,:,ik]
                HHH = np.dot(H,np.conj(H).T)
                IR  = np.eye(self.Nr)
                tC[ir,ik] = np.log(la.det(IR+(rho/self.Nt)*HHH))/np.log(2)

        self.C=tC


    def mulcplot(self,mode,**kwargs):
        """
        """

        defaults = { 'types' : ['m'],
                   'titles' : np.array([['11','12','13','14'],
                 ['21','22','23','34'],
                 ['31','32','33','34'],
                 ['41','42','43','44'],
                 ['51','52','53','54'],
                 ['61','62','63','64'],
                 ['71','72','73','74'],
                 ['81','82','83','84']]),
                   'ylabels':np.array([['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','','']]),
                   'xlabels':np.array([['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['fGHz','fGHz','fGHz','fGHz']]),
                   'labels':np.array([['calibrated','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','','']])
                 }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if mode=='f':
            fig,ax = self.Hcal.plot(**kwargs)
        else:
            kwargs['xlabels'] = np.array([['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['','','',''],
                 ['t(ns)','t(ns)','t(ns)','t(ns)']]),
            fig,ax = self.hcal.plot(**kwargs)

        return(fig,ax)

        return fig,ax

    def plot(self,frequency=True,phase=False,dB=True,cal=True,fig=[],ax=[],color='k'):
        """

        """
        if fig==[]:
            fig,ax=plt.subplots(8,self.Nt,sharex=True,sharey=True)
        if cal:
            H = self.Hcal
        else:
            H = self.H
        for iR in range(self.Nr):
            for iT in range(self.Nt):
                k = iR*4+iT
                if frequency:
                    if not phase:
                        if dB:
                            #ax[iR,iT].plot(H.x,20*np.log10(abs(H.y[k,:])),color=color) 
                            ax[iR,iT].plot(H.x,20*np.log10(abs(H.y[iR,iT,:])),color=color) 
                        else:
                            #ax[iR,iT].plot(H.x,abs(H.y[k,:]),color='k') 
                            ax[iR,iT].plot(H.x,abs(H.y[iR,iT,:]),color='k') 
                    else:
                        #ax[iR,iT].plot(H.x,np.unwrap(np.angle(H.y[k,:])),color=color) 
                        ax[iR,iT].plot(H.x,np.unwrap(np.angle(H.y[iR,iT,:])),color=color) 
                else:
                        ax[iR,iT].plot(self.h.x,abs(self.h.y[iR,iT,:]),color=color) 
                if (iR==7):         
                    ax[iR,iT].set_xlabel('f (GHz)') 
                ax[iR,iT].set_title(str(iR+1)+'x'+str(iT+1)) 
        return(fig,ax)

