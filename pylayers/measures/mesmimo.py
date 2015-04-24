from pylayers.signal.bsignal import *
from pylayers.util.project import *
from pylayers.gis.readvrml import *
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import numpy.linalg as la
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
        st = 'MIMO Object'+'\n'
        st = st + 'Nr : '+str(self.Nr)+ '\n'
        st = st + 'Nt : '+str(self.Nt)+ '\n'
        st = st + 'Nf : '+str(self.Nf)+ '\n'
        return(st)

    def __sub__(self,m):
        N = MIMO()
        N.freq = self.freq
        N.Nt = self.Nt
        N.Nr = self.Nr
        N.Hcal = self.Hcal - m.Hcal
        return(N)

    def loadraw(self):
        """ load a MIMO Nr x Nt raw data sounder file

    
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


    def normalize(self):
        """ Normalization of H

        """

        HdH,U,S,V = self.transfer()
        self.rg = np.sqrt(trace(HdH)/(self.Nt*self.Nr))
        self.Hcal = self.Hcal/self.rg
        self.normalize=True



    def transfer(self):
        """ calculate transfer matrix

        Returns
        -------

        HdH : Hermitian transfer matrix
        U   : Unitary tensor
        L   : Singular values
        V   : = Ud (in that case because HdH Hermitian)

        HdH = U L U^{\dagger}

        """

        H   = self.Hcal.y
        Hd  = np.conj(self.Hcal.y.swapaxes(0,1))
        HdH = np.einsum('ijk,jlk->ilk',Hd,H)
        HdH  = HdH.swapaxes(0,2)
        U,S,V  = la.svd(HdH)

        return (HdH,U,S,V)

    def Bcapacity(self,Pt=1e-3,Tp=273):
        """ calculates BLAST deterministic MIMO channel capacity

        Parameters
        ----------

        Pt : float
            the total power is assumed uniformaly distributed over the whole bandwidth
        Tp : Receiver Temperature (K)

        Returns
        -------

        C : capacity (bit/s)

            log_2(det(I+(Et/(N0Nt))HH^{H})

        """

        fGHz  = self.Hcal.x
        Nf    = len(fGHz)
        BGHz  = fGHz[-1]-fGHz[0]
        dfGHz = fGHz[1]-fGHz[0]

        # White Noise definition
        #
        # Boltzman constantf    = len(fGHz)

        kB = 1.03806488e-23

        # N0 ~ J ~ W/Hz ~ W.s

        N0 = kB*Tp


        # Evaluation of the transfer tensor

        HdH,U,S,V = self.transfer()

        It  = np.eye(self.Nt)
        Ir  = np.eye(self.Nr)

        #Ps = (Pt/Nf)/(self.Nt)
        Ps = Pt/(self.Nt)
        Pb = N0*BGHz*1e9
        #coeff = Ps/Pb
        #M     = It[None,...] + coeff*HdH
        #detM  = la.det(M)
        #logdetM = np.real(np.log(detM)/np.log(2))
        #C1  = dfGHz*logdetM
        CB  = dfGHz*np.sum(np.log(1+(Ps/Pb)*S)/np.log(2),axis=1)
        #return(M,detM,logdetM,C1,C2,S)
        return(CB)

    def WFcapacity(self,Pt=1e-3,Tp=273):
        """ calculates deterministic MIMO channel capacity

        Parameters
        ----------

        Pt :  the total power to be distributed over the different spatial
            channels using water filling
        Tp : Receiver Temperature (K)

        Returns
        -------

        C : capacity (bit/s)

            log_2(det(It + HH^{H})

        """

        fGHz  = self.Hcal.x
        Nf    = len(fGHz)
        BGHz  = fGHz[-1]-fGHz[0]
        dfGHz = fGHz[1]-fGHz[0]

        Hp  = self.Hcal.y.swapaxes(1,2)
        H   = Hp.swapaxes(0,1)
        Hd  = np.conj(H.swapaxes(1,2))

        # White Noise definition
        #
        # Boltzman constant

        kB = 1.03806488e-23

        # N0 ~ J ~ W/Hz ~ W.s

        N0 = kB*Tp


        # Evaluation of the transfer tensor

        HdH,U,ld,V = self.transfer()

        It  = np.eye(self.Nt)

        Ir  = np.eye(self.Nr)

        pb = N0*dfGHz*1e9*np.ones((self.Nf,self.Nt))
        pt = Pt/((self.Nf-1)*self.Nt)
        mu = pt
        Q0 = np.maximum(0,mu-pb/ld)
        u  = np.where(Q0>0)[0]
        Nnz1  = len(u)
        
        Peff = np.sum(Q0)
        deltamu = pt/100.
        while np.abs(Peff-Pt)>1e-16:
            mu = mu + deltamu
            Q = np.maximum(0,mu-pb/ld)
            Peff = np.sum(Q)
            #print "mu , Peff : ",mu,Peff
            if Peff>Pt:
                mu = mu - deltamu
                deltamu = deltamu/2.
        #print Peff
        v  = np.where(Q>0)[0]
        Nnz2  = len(v)
        #print self.Nf*self.Nt-Nnz1,self.Nf*self.Nt-Nnz2
        #Qoptn = Qopt/pb
        Qn   = Q/pb
        #Qn_e = Qn[:,:,None]*It[None,:,:]
        #pb_e = pb[:,:,None]*It[None,:,:]
        #print "Q : ",Qn_e.shape
        #print "Hd : ",Hd.shape
        #print "H : ",H.shape
        # k :frequency (i.e 1601)
        #QHd  = np.einsum('kij,kil->kil',Qn_e,Hd)
        #QHdH = np.einsum('kil,klj->kij',QHd,H)
        #print "QHdH : ",QHdH.shape
        #M = It[None,...] + QHdH
        # Pt / df ~ J
        # print 10*np.log10(Ps*np.real(HH[0,0,0])/Pb)
        #detM  = la.det(M)
        #logdetM = np.real(np.log(detM)/np.log(2))
        #C1  = dfGHz*logdetM
        Cwf  = dfGHz*np.sum(np.log(1+(Qn)*ld)/np.log(2),axis=1)
        #C   = dfGHz*np.log(la.det(IR[None,...]+(Pt/self.Nt)*HH/(N0*dfGHz)))/np.log(2)

        return(Cwf,Q)

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

