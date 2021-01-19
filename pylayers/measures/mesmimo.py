#!/usr/bin/python
#-*- coding:Utf-8 -*-
from __future__ import print_function
from pylayers.signal.bsignal import *
from pylayers.antprop.aarray import *
from pylayers.util.project import *
from pylayers.antprop.channel import *
from pylayers.gis.readvrml import *
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import scipy as sp
import scipy.special as ss
import numpy.linalg as la
from time import sleep
import math as mt
from pylayers.measures.vna.E5072A import *
"""
.. curentmodule:: pylayers.antprop.mesmimo

.. autosummary::
    :members:

"""

class MIMO(object):
    """ This class handles the data coming from a MIMO Channel Sounder

    Parameters
    ----------
    H    : raw channel matrix in frequency domain
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
                    'Aar': [],
                    'snrdB': np.linspace(0,25,100)
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
        self.snrdB = kwargs.pop('snrdB')

        self.Aat = kwargs.pop('Aat')
        self.Aar = kwargs.pop('Aar')
        if self.Aar == []:
            #self.Aar = AntArray(N=[8,1,1])
            self.Aar = UArray(N=[8,1,1])
        if self.Aat == []:
            #self.Aat = AntArray(N=[4,1,1])
            self.Aat = UArray(N=[4,1,1])


        self.Nf  = Nf
        self.freq = np.linspace(fminGHz,fmaxGHz,Nf)
        self.rep = rep
        self.Nt = Nt
        self.Nr = Nr

        #pdb.set_trace()
        if _filename != '':
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
        module   = []
        phasedeg = []

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

        #MIMO
        # Nr x Nt x Nf
        tc = tc.reshape(self.Nr,self.Nt,self.Nf)
        # C.freq , Nf
        self.C = Tchannel(x=C.freq,y=tc)
        self.Hcal = self.H/self.C
        del self.H
        del self.C

    def  calHa(self,**kwargs):
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
        defaults = {'fcGHz':2,
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

        Transfered to Mchannel
        DONE
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

        Transfered to Mchannel

        Done

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

        Pb  = N0*BGHz*1e9   # Watt

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
        # Bandwidth
        BGHz  = fGHz[-1]-fGHz[0]
        # Frequency step
        dfGHz = fGHz[1]-fGHz[0]

        # White Noise definition
        #
        # Boltzman constant

        kB = 1.03806488e-23

        # N0 ~ J ~ W/Hz ~ W.s

        N0 = kB*Tp


        # Evaluation of the transfer tensor

        HdH,U,ld,V = self.transfer()

        # Identity matrices

        It = np.eye(self.Nt)
        Ir = np.eye(self.Nr)

        #
        # Iterative implementation of Water Filling algorithm
        #

        # pb : (nf,nt)   noise power (Watt)
        pb = N0*dfGHz*1e9*np.ones((self.Nf,self.Nt))
        # pt : (nf,nt,power)  Total power uniformly spread over (nt*nf-1)
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

    # def ber(self,cmd='QPSK',m = 4, snrdB = np.linspace(0,25,100)):
    #     """computation of bit error rate

    #     Parameters
    #     ----------
    #     cmd : 'QPSK' or M-QAM

    #     M : number of bit (int) (2 or 4 or 8)
    #     """
        
    #     snr = 10**(snrdB/10.)
    #     M = 2**m
    #     eta = np.log(M, 2)

    #     if cmd == 'QPSK':
    #         berqpsk = 0.5 * ss.erfc(sqrt(snr))

    #     if cmd == 'M-PSK':
    #         bermpsk =  1 / eta * ss.erfc(sqrt(snr * eta) * np.sin(np.pi / M))

    #     if cmd == 'M-QAM':
    #         bermqam = 2 / eta * (1 - 1 / sqrt(M)) * ss.erfc(sqrt(3 * snr * eta/(2 * (M - 1)))
        
    #     return(berqpsk,bermpsk,bermqam)

    # def berplot(self):
    #     """plot BER functions
    #     """

    #     berqpsk,bermpsk,bermqam = self.ber(cmd='',m = 4, snrdB = np.linspace(0,25,100))

    #     if cmd == 'QPSK':
    #         plt.semilogy(snrdB,berqpsk,label='QPSK')

    #     if cmd == 'M-PSK':
    #         plt.semilogy(snrdB,bermqpsk,label='QPSK')

    #     if cmd == 'M-QAM':
    #         plt.semilogy(snrdB,bermqam,label='4-PSK')

    #     sns.set_style("darkgrid")
    #     plt.ylim([10**-9, 0.5])
    #     plt.figure(figsize=(20,20))
    #     plt.xlabel('SNR(dB)',fontsize=15)
    #     plt.ylabel('Bit Error Rate',fontsize=15)
    #     plt.legend(loc='best')
    #     plt.title("Digital Modulation Bit Error Rate")
    #     plt.show()

    def linear_ZF(self,cmd='QPSK',m = 4, snrdB = np.linspace(0,25,100)):
        """linear Zero Forcing precoding
        Parameters
        ----------
        
        """

        # H  : nr x nt x nf
        H = self.Hcal.y
        # Hd : nt x nr x nf
        Hd  = np.conj(self.Hcal.y.swapaxes(0,1))
        H_inv = np.linalg.inv(H)
        H_inv_d = np.transpose(H_inv)
        tr_mat = np.matrix.trace(H_inv*H_inv_d)
        beta = sqrt(self.Nt/(tr_mat))
        W_zf = np.dot(beta,H_inv)



    def linear_MMSE(self,cmd='QPSK',m = 4, snrdB = np.linspace(0,25,100)):
        """linear MMSE precoding

        Parameters
        ----------
        
        """
        # H  : nr x nt x nf
        H = self.Hcal.y
        # Hd : nt x nr x nf
        Hd  = np.conj(self.Hcal.y.swapaxes(0,1))
        HHd =np.einsum('ijk,jlk->ilk',H,Hd)
        Hh = np.transpose(H)
        H_inv = np.linalg.inv(H)
        H_inv_d = np.transpose(H_inv)
        tr_mat = np.matrix.trace(H_inv*H_inv_d)
        beta = sqrt(self.Nt/(tr_mat))
        
        Pt = np.logspace(-3,1,100)
        kB = 1.3806488e-23
        N0 = kB*273
        B  = 400e6
        Pb = N0*B
        A = np.linalg.inv(HHd + snr)
        B = np.dot(Hh,A)
        W_mmse = beta * B



    # def meas(self):
    #     """ Allows meas from VNA and Scanner
    #     """

    #     defaults = { 'lavrg':'['1','999']',
    #                  'lif':'['1000','300000','500000']',
    #                  'lpoints' : '[201,401,601,801,1601]',
    #                  'Nf':1601,
    #                  'fminGHz' : 1.8,
    #                  'fmaxGHz' :2.2,
    #                  'calibration':True,
    #                  'time':True,
    #                  'Nmeas' : 100,
    #                  'Nt' : 4,
    #                  'Nr' : 8,
    #                  'Aat': [],
    #                  'Aar': []
    #               }

    #     for k in defaults:
    #         if k not in kwargs:
    #             kwargs[k]=defaults[k]

    #     fminGHz = kwargs.pop('fminGHz')
    #     fmaxGHz = kwargs.pop('fmaxGHz')
    #     lavrg   =  kwargs.pop('lavrg')
    #     lif     = kwargs.pop('lif')
    #     lpoints = kwargs.pop('lpoints')
    #     Nmeas = kwargs.pop('Nmeas')


    #     ##################
    #     ### VNA
    #     #################


    #     # FROM MAIN OF E5072A.py
    #     vna = SCPI("129.20.33.201",verbose=False)
    #     ident = vna.getIdent()
    #     print "Talking to : ",ident
    #     vna.write("FORM:DATA REAL")
    #     #vna.write("SENS:AVER:ON")
    #     vna.select(param='S21',chan=1)
    #     print "channel "+str(chan)+ " selected"
    #     vna.setf(startGHz=1.8,stopGHz=2.2)
    #     print "fstart (GHz) : ",startGHz
    #     print "fstop (fGHz) : ",stopGHz


    #     ######
    #     vna.setf(fminGHz,fmaxGHz)
    #     prefix = 'cal_'
    #     S = []
    #     lt = []

    #     tic = time.time()

    #     for i in lif:
    #         vna.write(":SENS1:BAND " + str(i))
    #         for n in lpoints:
    #             fGHz = np.linspace(startGHz,stopGHz,n)
    #             vna.setnpoint(n)
    #             com = ":CALC1:DATA:SDAT?\n"
    #             npts = vna.getnpoints()
    #             print "Nbrs of points : ",npts
    #             S = vna.getdata(n)
    #             lt.append(time.time())
    #             try:
    #                 S21.append(S)
    #             except:
    #                 S21=S
    #             S.save(prefix+str(n))
    #             #for k in range(Nmeas):
    #                 #S = vna.getdata(Npoints=Npoints)
    #                 #lt.append(time.time())
    #                 #try:
    #                     #S21.append(S)
    #                 #except:
    #                     #S21=S
    #     toc = time.time()
    #     print toc-tic
    #     #lt.append(toc-tic)
    #     #lS.append(S21)
    #     #del S21
    #     #vna.close()
    #     #S21.save('calibration.mat')


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


    def grid(self,M,
             OR=np.array([3.4,0.73]),
             OT=np.array([5.29,6.65]),
             cT=np.array([-0.07,0]),
             cR=np.array([0.07,0])):
        """ Evaluate the data on a grid in the plane

            Parameters
            ----------

            M : np.array() (Nx x Ny)
            OR : np.array (,2)
                Origin of receiver [3.4,0.73]
            OT : np.array (,2)
                Origin of transmitter [5.29,6.65]
            cR : np.array (,2)
                array receiving vector [0.07,0]
            cT : np.array (,2)
                array transmitting vector [-0.07,0]

            Notes
            -----

           Updated object members

            self.grid : M  (Nx x Ny x 2)
            self.gloc : TUsignal (x (,ntau) y (Nx x Ny,ntau) )

        """

        aR = cR[0]/np.sqrt(cR[0]**2+cR[1]**2)
        bR = cR[1]/np.sqrt(cR[0]**2+cR[1]**2)

        aT = cT[0]/np.sqrt(cT[0]**2+cT[1]**2)
        bT = cT[1]/np.sqrt(cT[0]**2+cT[1]**2)
        # mapping
        uT = (aT*(M[...,0]-OT[0])+bT*(M[...,1]-OT[1]))/np.sqrt((M[...,0]-OT[0])**2+(M[...,1]-OT[1])**2)
        uR = (aR*(M[...,0]-OR[0])+bR*(M[...,1]-OR[1]))/np.sqrt((M[...,0]-OR[0])**2+(M[...,1]-OR[1])**2)
        # sampling in uR and uT
        uuR = self.uR
        uuT = self.uT
        # index in uR and uT
        iUr=np.array(map(lambda x : np.where(abs(uuR-x)==(abs(uuR-x)).min())[0][0], np.ravel(uR)))
        iUt=np.array(map(lambda x : np.where(abs(uuT-x)==(abs(uuT-x)).min())[0][0], np.ravel(uT)))

        self.grid = M
        shM = M.shape
        self.gloc = TUsignal(self.h.x,self.h.y[iUr,iUt,:])
        #self.gloc = self.h[iUr,iUt,:]
        #shL =  gloc.shape
        #assert(shL[0]==shM[0]*shM[1])
        #self.gloc = np.reshape(gloc,(shM[0],shM[1],shL[1]))


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

    def showgrid(self,**kwargs):
        """ show the data on a spatial grid

        Parameters
        ----------

        layout:[],
        s:50,
        vmin : 0,
        vmax: 0.5,
        linewidth:0,
        fig:[],
        ax:[],
        save:True,
        filename:'showgrid1',
        title:'',
        save:True,
        dB : False,
        OR : np.array([3.4,0.73]),
        OT : np.array([5.29,6.65]),
        cR : np.array([0.07,0]),
        cT : np.array([-0.07,0]),
        target : np.array([]),
        gating : False,
        dynamic : 30


        Notes
        -----

        This function accepts a Layout as input and allows to display
        a projection of the spatio-delay volume on a 2D grid.


        """
        defaults = { 'layout':[],
                    's':50,
                    'vmin' : 0,
                    'vmax': 0.5,
                    'linewidth':0,
                    'fig':[],
                    'ax':[],
                    'save':True,
                    'filename':'showgrid1',
                    'title':'',
                    'save':True,
                    'dB':False,
                    'OR' : np.array([3.4,0.73]),
                    'OT' : np.array([5.29,6.65]),
                    'cR' : np.array([0.07,0]),
                    'cT' : np.array([-0.07,0]),
                    'target' : np.array([]),
                    'gating':False
                   }


        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        OR = kwargs['OR']
        OT = kwargs['OT']
        cR = kwargs['cR']
        cT = kwargs['cT']

        ULAR = OR+np.arange(8)[:,np.newaxis]*cR-3.5*cR
        ULAT = OT+np.arange(4)[:,np.newaxis][::-1]*cT-1.5*cT

        if kwargs['gating']:
            dTM = np.sqrt((self.grid[...,0]-OT[0])**2+(self.grid[...,1]-OT[1])**2)
            dRM = np.sqrt((self.grid[...,0]-OR[0])**2+(self.grid[...,1]-OR[1])**2)
            # dM : Nx,Ny
            dM  = dTM+dRM
            # dM : ,Nx x Ny
            dM = np.ravel(dM)
            # 6 sigma = 1/400MHz
            # 6 sigma = 2.5ns
            # sigma = (2.5/6)
            # alpha = 1/(2 sigma^2) = 2*(2.5)**2/36 = 0.347
            #
            alpha = 0.347
            # Gaussian gate
            # Laplacian gate
            # Nx x Ny x Ntau
            self.gate = np.exp(-alpha*(dM[:,np.newaxis]/0.3-self.gloc.x[np.newaxis,:])**2)
            data = self.gloc.y*self.gate
            data = np.sum(abs(data),axis=1)
        else:
            data = np.sum(abs(self.gloc.y),axis=1)

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=(10,10))
            ax  = fig.add_subplot(111)
        else:
            fig=kwargs['fig']
            ax = kwargs['ax']

        if kwargs['dB']:
            data = 20*np.log10(data)
            vmax = data.max()
            # clipping @ vmax - dynamic
            vmin = vmax-kwargs['dynamic']
        else:
            vmin = data.min()
            vmax = data.max()

        scat = ax.scatter(self.grid[...,0],
                               self.grid[...,1],
                               c= data,
                               s=kwargs['s'],
                               vmin=vmin,
                               vmax=vmax,
                               linewidth=kwargs['linewidth'])

        cb = plt.colorbar(scat)
        if kwargs['dB']:
            cb.set_label('Level (dB)')
        else:
            cb.set_label('Linear Level')


        # plot ULAs

        ax.plot(ULAR[:,0],ULAR[:,1],'+b')
        ax.plot(ULAT[:,0],ULAT[:,1],'+g')
        plt.axis('off')

        # plot target

        if kwargs['target']!=[]:
            target = ax.scatter(kwargs['target'][0],kwargs['target'][1],c='black',s=100)

        # display layout
        if kwargs['layout'] != []:
            L = kwargs['layout']
            #fig,ax = L.showG('s',fig=fig,ax=ax,nodes=False)
            L.display['ednodes']=False
            L.display['nodes']=False
            L.display['title']=kwargs['title']
            fig,ax = L.showG('s',fig=fig,ax=ax,nodes=False)

        if kwargs['save']:
            fig.savefig(kwargs['filename']+'.pdf')
            fig.savefig(kwargs['filename']+'.png')

        return fig,ax

    def animgrid(self,**kwargs):
        """
        """

        defaults = { 'layout':[],
                    's':100,
                    'vmin' : 0,
                    'vmax': 0.5,
                    'linewidth':0,
                    'fig':[],
                    'ax':[],
                    'filename':'animgrid1',
                    'save':True,
                    'abs':True,
                    'title':'',
                   }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value


        if kwargs['fig']==[]:
            fig = plt.figure(figsize=(20,20))
            ax  = fig.add_subplot(111)

        if kwargs['layout'] != []:
            L = kwargs['layout']
            fig,ax = L.showG('s',fig=fig,ax=ax,nodes=False)

        Nframe = self.gloc.y.shape[1]
        if kwargs['abs']:
            scat = ax.scatter(self.grid[...,0],
                               self.grid[...,1],
                               c=abs(self.gloc.y[:,0]),
                               s=kwargs['s'],
                               vmin=kwargs['vmin'],
                               vmax=kwargs['vmax'],
                               linewidth=kwargs['linewidth'])
        else:
            scat = ax.scatter(self.grid[...,0],
                               self.grid[...,1],
                               c=self.gloc.y[:,0],
                               s=kwargs['s'],
                               vmin=kwargs['vmin'],
                               vmax=kwargs['vmax'],
                               linewidth=kwargs['linewidth'])

        title  = ax.text(0.1,0.9,kwargs['title'],transform=ax.transAxes,fontsize=18)
        cb   = plt.colorbar(scat)
        delay_template = '%d : tau = %5.2f (ns) d= %5.2f (m)'
        delay_text  = ax.text(0.1,0.9,'',transform=ax.transAxes,fontsize=18)

        def init():
            delay_text.set_text('')
            if kwargs['abs']:
                scat.set_array(abs(self.gloc.y[:,0]))
            else:
                scat.set_array(self.gloc.y[:,0])
            return scat,delay_text

        def animate(i):
            delay_text.set_text(delay_template%(i,self.gloc.x[i],self.gloc.x[i]*0.3))
            if kwargs['abs']:
                scat.set_array(abs(self.gloc.y[:,i]))
            else:
                scat.set_array(abs(self.gloc.y[:,i]))
            return scat,delay_text

        anim = animation.FuncAnimation(fig,
                                       animate,
                                       init_func=init,
                                       frames=Nframe,
                                       interval=1,
                                       blit=True)
        if kwargs['save']:
            anim.save(kwargs['filename']+'.mp4', fps=5)
        return fig,ax,anim

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

