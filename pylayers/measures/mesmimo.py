#!/usr/bin/python
#-*- coding:Utf-8 -*-
from pylayers.signal.bsignal import *
from pylayers.antprop.channel import *
from pylayers.antprop.aarray import *
from pylayers.util.project import *
from pylayers.antprop.channel import *
from pylayers.gis.readvrml import *
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import numpy.linalg as la


class MIMO(object):
    """ This class handles the data coming from the MIMO Channel Sounder IETR lab

    Parameters
    ----------
    H    : raw channel matrix in frequencey domain
    Hcal : calibrated channel matrix in frequency domain
    hcal : channel matrix in time domain

    """
    def __init__(self,**kwargs):
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

        defaults = { '_filename':'',
                    'rep':'',
                    'Nf':1601,
                    'fminGHz' : 1.8,
                    'fmaxGHz' :2.2,
                    'calibration':True,
                    'time':True,
                    'Nz' : 100,
                    'Nt' : 4,
                    'Nr' : 8,
                    'Aat': [],
                    'Aar': []
                  }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        _filename = kwargs.pop('_filename')
        rep = kwargs.pop('rep')
        Nf =  kwargs.pop('Nf')
        fminGHz = kwargs.pop('fminGHz')
        fmaxGHz = kwargs.pop('fmaxGHz')
        calibration = kwargs.pop('calibration')
        time = kwargs.pop('time')
        Nz = kwargs.pop('Nz')
        Nt = kwargs.pop('Nt')
        Nr = kwargs.pop('Nr')

        self.Aat = kwargs.pop('Aat')
        self.Aar = kwargs.pop('Aar')
        if self.Aar == []:
            self.Aar = AntArray(N=[8,1,1])
        if self.Aat == []:
            self.Aat = AntArray(N=[4,1,1])

        self.freq = np.linspace(fminGHz,fmaxGHz,Nf)

        self.Nf  = Nf
        self.rep = rep
        self.Nt = Nt
        self.Nr = Nr

        #pdb.set_trace()
        if _filename <> '':
            self.filename = mesdir + rep + _filename
            # load file
            self.loadraw()
            if calibration:
                self.calibration()
                if time:
                    # reshaping for using ift (todo update ift for MDA !!)
                    #Hcal = TChannel(x=self.Hcal.x,y=np.reshape(self.Hcal.y,(Nt*Nr,Nf)))
                    Hcal = Tchannel(self.Hcal.x,np.reshape(self.Hcal.y,(Nt*Nr,Nf)))
                    hcal = Hcal.ift(Nz=Nz,ffts=1)
                    shh = hcal.y.shape
                    self.hcal = TUsignal(hcal.x,np.reshape(hcal.y,(Nr,Nt,shh[-1])))


    def __repr__(self):
        st = 'MIMO Object'+'\n'
        st = st + 'axe 0  Nr : '+str(self.Nr)+ '\n'
        st = st + 'axe 1  Nt : '+str(self.Nt)+ '\n'
        st = st + 'axe 2  Nf : '+str(self.Nf)+ '\n'
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

        self.H = Tchannel(x=self.freq,y=y)

    def calibration(self):
        """ Apply calibration files

        """
        for iR in range(self.Nr):
            for iT in range(self.Nt):
                _filename = 'Calib'+str(iT+1)+'x'+str(iR+1)+'.txt'
                C = MIMO(_filename=_filename,rep='/calibration/',calibration=False,Nt=self.Nt)
                try:
                    #tc = np.vstack((tc,C.H.y[iR*4+iT,:]))
                    tc = np.vstack((tc,C.H.y[iR,iT,:]))
                except:
                    #tc = C.H.y[iR*4+iT,:]
                    tc = C.H.y[iR,iT,:]

        # Nr x Nt x Nf
        tc = tc.reshape(self.Nr,self.Nt,self.Nf)

        # C.freq , Nf

        self.C = Tchannel(x=C.freq,y=tc)


        self.Hcal = self.H/self.C

        del self.H
        del self.C


    def calHa(self,**kwargs):
        """ calculate the Ha function (angular domain representation)

        fcGHz : float
        duR   : grid step in uR
        duT   : grid step in uT
        time  : boolean
        taumin  : float 0
        taumax  : float
        Nz   : int (20000)

        See : David Tse (7.70 pp 373)

        """
        defaults = { 'fcGHz':2,
                    'duR':0.05,
                    'duT':0.05,
                    'time':False,
                    'taumin':0,
                    'taumax':80,
                    'Nz':20000
                    }
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        fcGHz = kwargs.pop('fcGHz')
        duR = kwargs.pop('duR')
        duT = kwargs.pop('duT')
        time = kwargs.pop('time')
        taumin = kwargs.pop('taumin')
        taumax = kwargs.pop('taumax')
        Nz = kwargs.pop('Nz')

        # f : m x n x uR x f
        fGHz = self.freq[None,None,None,:]
        # m : m x n x  uR x f
        m = np.arange(self.Nr)[:,None,None,None]
        # uR : m x n x uR x f
        uR = np.arange(-1,1,duR)[None,None,:,None]
        # eR : m x n x uR x f
        eR = np.exp(-1j*np.pi*m*uR*fGHz/fcGHz)
        # S : m x n x uR x f
        S  = self.Hcal.y[:,:,None,:] * eR
        #  SR : n x uR x uT x f
        SR = np.sum(S,axis=0)[:,:,None,:]
        # n : n x uR x uT x f
        n = np.arange(self.Nt)[:,None,None,None]
        # uT : n x uR x uT x f
        uT = np.arange(-1,1,duT)[None,None,:,None]
        # eT : n x uR x uT x f
        eT = np.exp(-1j*np.pi*n*uT*fGHz/fcGHz)
        # summation along axix m and n
        self.Ha = np.sum(SR*eT,axis=0)

        self.uR = np.arange(-1,1,duR)
        self.uT = np.arange(-1,1,duT)

        NuR = len(self.uR)
        NuT = len(self.uT)
        Nf  = len(self.freq)

        if time:
            #T = fft.ifft(self.h,axis=2)
            #self.h = abs(fft.fftshift(T,axes=2))
            Ha = FUsignal(self.freq,np.reshape(self.Ha,(NuR*NuT,Nf)))
            ha = Ha.ift(Nz=Nz,ffts=1)
            ut = np.where((h.x>taumin) & (h.x<taumax))[0]
            xlim = ha.x[ut]
            ylim = ha.y[...,ut]
            npts = len(ut)
            self.ha = TUsignal(xlim,np.reshape(ylim,(NuR,NuT,npts)))

    def normalize(self):
        """ Normalization of H

        """

        HdH,U,S,V = self.transfer()
        HdH = HdH.swapaxes(0,2)
        self.rg = np.real(np.sqrt(np.trace(HdH)/(self.Nt*self.Nr)))
        self.Hcal.y = self.Hcal.y/self.rg
        self.normalize=True


    def svd(self):
        """ singular value decomposition of matrix H

        Parameters
        ----------

        The native H matrix is currently (nr x nt x nf ). For applying a
        broadcasted svd a reshaping in (nf x nr x nt ) is required.
        In the future,  it would be a good thing to define the MIMO matrix as

        nf x na x nb structure from the begining

        or

        ns x nf x na x nb

        Returns
        -------

        U  : nf x nr x nr
        D  : nf x min(nr,nt)
        Vh : nf x nt x nt

        """
        # H  : nr x nt x nf
        H  = self.Hcal.y
        # H1  : nf x nt x nr
        H1 = H.swapaxes(0,2)
        # H2  : nf x nr x nt
        H2 = H1.swapaxes(1,2)
        U,D,Vh = la.svd(H2)
        return(U,D,Vh)


    def transfer(self):
        """ calculate transfer matrix.
            it involves H and Hd against svd() which acts only over H.

        Returns
        -------

        HdH : Hermitian transfer matrix  (nf x nt x nt )
        U   : Unitary tensor  (nf x nt x nt )
        S   : Singular values (nf x nt)
        V   : = Ud (in that case because HdH Hermitian)  (nf x nt x nt)

        HdH = U L U^{\dagger}

        """

        # H  : nr x nt x nf
        H   = self.Hcal.y
        # Hd : nt x nr x nf
        Hd  = np.conj(self.Hcal.y.swapaxes(0,1))
        #HdH : nt x nt x nf
        HdH = np.einsum('ijk,jlk->ilk',Hd,H)
        # HdH : nf x nt x nt
        HdH  = HdH.swapaxes(0,2)
        #U   : nf x nt x nt
        #S   : nf x nt
        #V   : nf x nt x nt
        U,S,V  = la.svd(HdH)

        return (HdH,U,S,V)

    def Bcapacity(self,Pt=np.array([1e-3]),Tp=273):
        """ calculates BLAST deterministic MIMO channel capacity

        Parameters
        ----------

        Pt : np.array (,NPt)
            the total power is assumed uniformaly distributed over the whole bandwidth
        Tp : Receiver Temperature (K)

        Returns
        -------

        C   : spectral efficiency (bit/s/Hz)
            np.array (Nf,NPt)
        rho : SNR
            np.array (Nf,Nt,NPt)

            log_2(det(I+(Et/(N0Nt))HH^{H})

        """

        fGHz  = self.Hcal.x
        Nf    = len(fGHz)
        BGHz  = fGHz[-1]-fGHz[0]
        dfGHz = fGHz[1]-fGHz[0]

        if type(Pt)==float:
            Pt=np.array([Pt])

        # White Noise definition
        #
        # Boltzman constantf    = len(fGHz)

        kB = 1.03806488e-23

        # N0 ~ J ~ W/Hz ~ W.s

        N0 = kB*Tp


        # Evaluation of the transfer tensor
        #
        # HdH :

        HdH,U,S,V = self.transfer()

        #singular value decomposition of channel tensor (broadcasted along frequency axis)


        Us,D,Vsh = self.svd()


        # Vsh : nf x nt x nt

        It  = np.eye(self.Nt)
        Ir  = np.eye(self.Nr)

        #Ps = (Pt/Nf)/(self.Nt)
        Ps  = Pt/(self.Nt)
        #Ps1 = Pt/(self.Nt*self.Nf)

        # equi amplitude vector (nf,nt,1)
        #wu  = np.sqrt(Ps[None,None,None,:]*np.ones((self.Nf,self.Nt))[:,:,None,None]/self.Nf)
        # spatial subchanel weights (nf,nt,1)
        #Vshwu = np.einsum('kijp,kjlp->kilp',Vsh[:,:,:,None],wu)
        # nf x nt x 1 x power
        # Ps2   = Vshwu*np.conj(Vshwu)

        Pb  = N0*BGHz*1e9

        #Pb2 = N0*dfGHz*1e9*np.ones((self.Nf,self.Nt))


        # rho : nf x nt x power
        #S2 = np.real(D[:,:,None]*np.conj(D[:,:,None]))
        #
        rho  = (Ps[None,None,:]/Pb)*S[:,:,None]
        #rho1 = (Ps1[None,None,:]/Pb2[:,:,None])*S[:,:,None]
        #rho2 = (Ps2[:,:,0,:]/Pb2[:,:,None])*S2
        #pdb.set_trace()
        #coeff = Ps/Pb
        #M     = It[None,...] + coeff*HdH
        #detM  = la.det(M)
        #logdetM = np.real(np.log(detM)/np.log(2))
        #C1  = dfGHz*logdetM

        #CB  = dfGHz*np.sum(np.log(1+rho)/np.log(2),axis=1)
        #CB  = dfGHz*np.sum(np.log(1+rho)/np.log(2))
        CB   = dfGHz*np.sum(np.log(1+rho)/np.log(2),axis=1)
        #CB1  = dfGHz*np.sum(np.log(1+rho1)/np.log(2),axis=1)
        #CB2  = dfGHz*np.sum(np.log(1+rho2)/np.log(2),axis=1)
        #return(M,detM,logdetM,C1,C2,S)
        return(rho,CB)

    def Scapacity(self,Pt=1e-3,Tp=273):
        """ equivalent SISO capacity
        """

        pass

    def BFcapacity(self,Pt=np.array([1e-3]),Tp=273):
        """ calculates the capacity in putting all the power on the more important mode

        Parameters
        ----------

        Pt : np.array  (,NPt)
            Transmitted power
        Tp : float
            Noise Temperature

        """
        fGHz  = self.Hcal.x
        Nf    = len(fGHz)
        BGHz  = fGHz[-1]-fGHz[0] # bandwidth
        dfGHz = fGHz[1]-fGHz[0]  # frequency step

        #
        # swaping axes
        #   self.Hcal.y  (Nr,Nt,Nf)
        #   Hp           (Nr,Nf,Nt)
        #   H            (Nf,Nr,Nt)
        #   Hd           (Nf,Nt,Nr)
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

        # pb : Nf x Nt
        pb = N0*dfGHz*1e9*np.ones((self.Nf,self.Nt))
        pt = Pt/((self.Nf-1))*np.array([1,0,0,0])[None,:]
        #print pt.shape
        Qn   = pt/pb
        rho  = Qn*ld
        #print Qn
        #print Qn.shape

        Cbf  = dfGHz*np.sum(np.log(1+rho)/np.log(2),axis=1)
        #C   = dfGHz*np.log(la.det(IR[None,...]+(Pt/self.Nt)*HH/(N0*dfGHz)))/np.log(2)
        return(Cbf,Qn)



    def WFcapacity(self,Pt=np.array([1e-3]),Tp=273):
        """ calculates deterministic MIMO channel capacity

        Parameters
        ----------

        Pt :  the total power to be distributed over the different spatial
            channels using water filling
        Tp : Receiver Noise Temperature (K)

        Returns
        -------

        C : capacity (bit/s)
        rho : SNR (in linear scale)

            log_2(det(It + HH^{H})

        """

        fGHz  = self.Hcal.x
        Nf    = len(fGHz)
        BGHz  = fGHz[-1]-fGHz[0]
        dfGHz = fGHz[1]-fGHz[0]

        # White Noise definition
        #
        # Boltzman constant

        kB = 1.03806488e-23

        # N0 ~ J ~ W/Hz ~ W.s

        N0 = kB*Tp


        # Evaluation of the transfer tensor

        HdH,U,ld,V = self.transfer()
        It = np.eye(self.Nt)
        Ir = np.eye(self.Nr)

        #
        # Iterative implementation of Water Filling algorithm
        #

        pb = N0*dfGHz*1e9*np.ones((self.Nf,self.Nt))
        pt = Pt[None,None,:]/((self.Nf-1)*self.Nt)
        mu = pt
        Q0 = np.maximum(0,mu-pb[:,:,None]/ld[:,:,None])
        u  = np.where(Q0>0)[0]

        Peff = np.sum(np.sum(Q0,axis=0),axis=0)
        deltamu = pt
        while (np.abs(Peff-Pt)>1e-16).any():
            mu = mu + deltamu
            Q = np.maximum(0,mu-pb[:,:,None]/ld[:,:,None])
            Peff = np.sum(np.sum(Q,axis=0),axis=0)
            #print "mu , Peff : ",mu,Peff
            usup = np.where(Peff>Pt)[0]
            mu[:,:,usup] = mu[:,:,usup]- deltamu[:,:,usup]
            deltamu[:,:,usup] = deltamu[:,:,usup]/2.
        Qn   = Q/pb[:,:,None]
        rho  = Qn*ld[:,:,None]

        Cwf  = dfGHz*np.sum(np.log(1+rho)/np.log(2),axis=1)


        return(rho,Cwf)


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

    def plot(self,**kwargs):
        """ plot channel

        Pramaters
        ---------

        frequency:True
        phase:True
        dB:True
        cal:True
        fig:[]
        ax:[]
        color':'k'

        """
        defaults = {'frequency':True,
                    'phase':False,
                    'dB':True,
                    'cal':True,
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        frequency = kwargs.pop('frequency')
        phase = kwargs.pop('phase')
        dB = kwargs.pop('dB')
        cal = kwargs.pop('cal')

        fig,ax=plt.subplots(8,self.Nt,sharex=True,sharey=True,**kwargs)

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
                            #ax[iR,iT].plot(H.x,20*np.log10(abs(H.y[iR,iT,:])),color='k')
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
                    ax[iR,iT].plot(H.x,np.unwrap(np.angle(H.y[iR,iT,:])),color='k')
                else:
                        ax[iR,iT].plot(self.hcal.x,abs(self.hcal.y[iR,iT,:]),color='k')
                if (iR==7):
                    ax[iR,iT].set_xlabel('Frequency (GHz)')
                    ax[iR,iT].set_title(str(iR+1)+'x'+str(iT+1))
        return(fig,ax)

