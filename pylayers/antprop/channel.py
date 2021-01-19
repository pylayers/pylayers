# -*t coding:Utf-8 -*-
"""
.. currentmodule:: pylayers.antprop.channel

.. autosummary::
    :members:

"""
from __future__ import print_function
import doctest
import pdb
import numpy as np
import numpy.ma as ma
import numpy.linalg as la
import scipy as sp
import scipy.signal as si
import pylab as plt
import struct as stru
import scipy.stats as st
import scipy.optimize as optimize
import numpy.fft as fft
from scipy.io import loadmat
import pylayers.util.pyutil as pyu
import pylayers.signal.bsignal as bs
import pylayers.util.geomutil as geu
import pylayers.antprop.antenna as ant
from pylayers.util.project import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fmin
import copy

try:
    import h5py
except:
    print('h5py is not installed: Ctilde(object cannot be saved)')

class AFPchannel(bs.FUsignal):
    """ Angular Frequency Profile channel

    Attributes
    ----------

    x : np.array
        frequency ,Nf
    y : np.array
        Amplitude Na,Nf
    tx : np.array
        tx coordinate (,3)
    rx : np.array
        rx coordinates (,3)
    az : np.array (,Na)
        AFP azimutal range in radians
    theta : link elevation angle
    phi : link (txrx) azimuth angle (with offset)
    tau : link delay (ns)

    offset : float angle in radians
        azimuth offset w.r.t global frame

    Methods
    -------

    norm2
    construct
    electrical_delay
    loadmes
    toadp
    estimate
    peak
    specular_model

    """
    def __init__(self,x = np.array([]),
                      y = np.array([]),
                      tx = np.array([]),
                      rx = np.array([]),
                      az = np.array([]),
                      label = '',
                      _filename = '',
                      refinement = False,
                      ang_offset = 0
                      ):
        bs.FUsignal.__init__(self,x=x,y=y,label='AFP')
        if len(self.x)!=0:
            self.fcGHz = self.x[int(len(self.x)/2)]
        self.tx = tx
        self.rx = rx
        self.label = label
        self.ang_offset = ang_offset
        self.refinement = refinement
        self.az = az
        self._filename = _filename

    def __add__(self,other):
        assert(self.y.shape == other.y.shape)
        assert((self.x == other.x).all())
        assert((self.az == other.az).all())
        A = AFPchannel()
        A.x = self.x
        A.fcGHz = A.x[int(len(A.x)/2)]
        A.az = self.az
        A.y = self.y + other.y
        return A

    def __sub__(self,other):
        assert(self.y.shape == other.y.shape)
        assert((self.x == other.x).all())
        assert((self.az == other.az).all())
        A = AFPchannel()
        A.x = self.x
        A.fcGHz = A.x[int(len(A.x)/2)]
        A.az = self.az
        Amod = np.abs(self.y)-np.abs(other.y)
        A.y = Amod*np.exp(1j*np.angle(self.y))
        return A

    def __repr__(self):
        cv = 180/np.pi
        s = 'Angular Frequency Profile: '+self.label+'\n'
        s = s + 'Tx : '+str(self.tx)+'\n'
        s = s + 'Rx : '+str(self.rx)+'\n'
        return s

    def norm2(self):
        return np.real(np.sum(self.y*np.conj(self.y)))

    def construct(self,tx):
        S = AFPchannel(x=self.x,y=np.zeros(self.y.shape),az=self.az)
        for k in range(tx.shape[0]):
            C = AFPchannel()
            C.specular_model(tx[k],self.x,self.az)
            S = S + C
        return S

    def loadmes(self,_filename,_filecal,fcGHz=32.6,BW=1.6,win='rect',ang_offset=0.37,ext='txt',dirmeas='meas/Espoo',refinement=False):
        """ Load measurement file

        Measurement files and the associated back to back calibration files
        are placed in the mes directory of the project.

        Parameters
        ----------

        _filename : string
            data matfile name
        _filecal : string
            calibration matfile name
        fcGHz : float
            center frequency
        BW : float
            measurement bandwidth
        win : string
            window type in ['rect','hamming','blackman']
        ang_offset : float
            angle in radian
        ext : string
            file extension 'txt' | '.mat'
        diremeas : string

        Notes
        -----

        This function updates :
        + self.x (frequency GHz)
        + self.y
        + self.az azimuth radians

        The calibration file _filecal (.mat file) should be added in the data directory
        In practice for Espoo B2B.mat

        See Also
        --------

        pylayers.util.pyutil.getlong


        """
        self._filename = _filename
        self.BW = BW
        self.fcGHz = fcGHz
        self.fmin = fcGHz-BW/2.
        self.fmax = fcGHz+BW/2.
        self.win = win
        self.refinement = refinement
        # read calibration file (Matlab file) in the same directory as measurments (convention)
        filecal = pyu.getlong(_filecal,dirmeas)
        U = loadmat(filecal)
        cal_trf = U['cal_trf'][:,0]
        # read measurement file (.txt or Mat file)
        filename = pyu.getlong(_filename,dirmeas)
        if ext=='txt':
            D = np.loadtxt(filename,skiprows=2)# load Back 2 Back calibration file
            amp = D[:,2::2]
            ang = D[:,3::2]
        else:
            D = loadmat(filename)
            amp = D['amp']
            ang = D['ang']
            rotationangle = D['rotationangle'].squeeze()

        # load Back 2 Back calibration file

        #
        # Transfer function reconstruction
        #

        self.Na  = amp.shape[0]
        self.Nf  = amp.shape[1]
        #
        # select apodisation window
        #
        if win=='hamming':
            window = np.hamming(self.Nf)
        elif win=='blackman':
            window = np.blackman(self.Nf)
        else:
            window = np.ones(self.Nf)
        #
        # complex transfer function
        #

        self.x = np.linspace(self.fmin,self.fmax,self.Nf)
        self.fcGHz = self.x[int(len(self.x)/2)]
        self.y = amp*np.exp(1j*ang*np.pi/180.)*cal_trf[None,:]*window
        #
        # if extension is txt file comes from ESPOO measurement
        #
        # self.az : 5.86 -> 1.94
        if ext=='txt':
            self.azmes = (360-D[:,0])*np.pi/180.
            self.az = self.azmes + ang_offset - 2*np.pi
            u = np.where(self.az<0)
            self.az[u] = self.az[u] + 2*np.pi
        else:
            self.azmes = rotationangle*np.pi/180.
            self.az = ang_offset - self.azmes
            u = np.where(self.az<0)
            self.az[u] = self.az[u] + 2*np.pi

    def  electrical_delay(self,tauns=0):
        """ electrical delay

        Parameters
        ----------
        tauns : float

        """
        self.y = self.y * np.exp(-2*1j*np.pi*self.x*tauns)

    def toadp(self,imax=-1):
        """ convert afp into adp (frequency->delay)

        Notes
        -----

        tx and rx need to be defined

        """
        # x : delay (starting at 0 ns)
        # y : ifft axis 1 (frequency)
        x = np.linspace(0,(len(self.x)-1)/(self.x[-1]-self.x[0]),len(self.x))
        y = np.fft.ifft(self.y,axis=1)
        if imax!=-1:
            y = y[:,0:imax]
            x = x[0:imax]
        adp = ADPchannel(x=x,
                         y=y,
                         az=self.az,
                         tx=self.tx,
                         rx=self.rx,
                         fcGHz=self.fcGHz,
                         _filename=self._filename,
                         refinement=self.refinement,
                         ang_offset = self.ang_offset)
        return adp

    def get_path(self):
        E = float(self.norm2())
        tEk = [E/200.,E]
        A = copy.copy(self)
        for k in range(10):
        #while ((tEk[-1]>tEk[0]) and (tEk[-1]>E/100)):
            xe,C,A = A.estimate()
            Energy = A.norm2()
            print(xe)
            tEk.append(Energy)
            try:
                txe = np.vstack((txe,xe))
            except:
                txe = xe[None,:]
        return(txe)


    def estimate(self,taumax=200,phimax=2*np.pi):
        """ estimate specular model parameters

        Parameters
        ----------

        taumax : float
        phimax : float

        See Also
        --------

        specular_model

        """
        def cost(xk,f,phi):
            B = AFPchannel()
            #B.specular_model2(xk,f,phi)
            B.specular_model(xk,f,phi)
            C = self - B
            return C.norm2()
        x_0 = self.peak()
        x_est = optimize.fmin_l_bfgs_b(cost,
                                       x_0,
                                       args=(self.x,self.az),
                                       disp=0,
                                       approx_grad=1,
                                       bounds=((0,2*x_0[0]),
                                               (0,taumax),
                                               (0,phimax)))[0]
        #x_est = optimize.fmin_l_bfgs_b(cost,
        #                               x_0,
        #                               args=(self.x,self.az),
        #                               disp=0,
        #                               approx_grad=1,
        #                               bounds=((0,2*x_0[0]),
        #                                       (0,2),
        #                                       (0,2),
        #                                       (0,taumax),
        #                                       (0,phimax)))[0]
        #x_est = optimize.fmin(cost,x_0,args=(self.x,self.az))
        Ck = AFPchannel()
        Ck.specular_model(x_est,self.x,self.az)
        #Ck.specular_model2(x_est,self.x,self.az)
        D = self - Ck
        return x_est,Ck, D

    def peak(self):
        adpself = self.toadp()
        ak,tauk,phik = adpself.peak()
        x = np.array([ak,tauk,phik])
        #x = np.array([ak,1,0,tauk,phik])
        return(x)

    def specular_model(self,x,fGHz,phi,wH=[],HPBW=10*np.pi/180,GmaxdB=21):
        """ Creates an AFP from a discrete specular model

        Parameters
        ----------
        x : [a0,a1,..,aK,tau0,tau1,...,tauk,phi0,...,phiK]
        fGHz :
        phi :
        wH : windowing on frequency axis
        HPBW : Half Power Beamwidth

        Examples
        --------
        >>> import numpy as np
        >>> rs = np.random.seed(1)
        >>> E = st.expon(0.5)
        >>> K = 5
        >>> tauk = 250*np.random.rand(K)
        >>> alphak = E.rvs(K)
        >>> phik  = 2*np.pi*np.random.rand(K)
        >>> xk   = np.hstack((alphak,tauk,phik))
        >>> A = AFPchannel()
        >>> fGHz = np.linspace(27,29,2001)
        >>> wH = np.ones(len(fGHz))
        >>> phi = np.linspace(0,2*np.pi,73)
        >>> A.specular_model(xk,fGHz,phi,wH)


        """
        K = int(len(x)/3)
        assert(len(x)==3*K)
        ak = x[0:K][:,None,None]
        tk = x[K:2*K][:,None,None]
        pk = x[2*K:3*K][:,None,None]
        # tf : paths (0) , freq (1), angle (2)
        if wH ==[]:
            wH = np.ones(len(fGHz))
        tf  = ak*np.exp(-2*1j*np.pi*fGHz[None,:,None]*tk)*wH[None,:,None]
        dphi = pk - phi[None,None,:]
        Gmax = 10**(GmaxdB/10.)
        g = np.exp(-(2*np.sqrt(np.log(2))*dphi/HPBW)**2)
        tfg = tf*g
        self.x = fGHz
        self.fcGHz = self.x[int(len(self.x)/2)]
        # self.y : angle(0) fGHz(1)
        self.y = np.sum(tfg,axis=0).T
        self.az = phi
        #h    = np.fft.ifft(H)

    def specular_model2(self,x,fGHz,phi,wH=[],HPBW=10*np.pi/180,GmaxdB=21):
        """ Creates an AFP from a discrete specular model

        Parameters
        ----------
        x : [a0,a1,..,aK,tau0,tau1,...,tauk,phi0,...,phiK]
        fGHz :
        phi :
        wH : windowing on frequency axis
        HPBW : Half Power Beamwidth

        Examples
        --------
        >>> import numpy as np
        >>> rs = np.random.seed(1)
        >>> E = st.expon(0.5)
        >>> K = 5
        >>> tauk = 250*np.random.rand(K)
        >>> alphak = E.rvs(K)
        >>> phik  = 2*np.pi*np.random.rand(K)
        >>> xk   = np.hstack((alphak,tauk,phik))
        >>> A = AFPchannel()
        >>> fGHz = np.linspace(27,29,2001)
        >>> wH = np.ones(len(fGHz))
        >>> phi = np.linspace(0,2*np.pi,73)
        >>> A.specular_model(xk,fGHz,phi,wH)


        """
        K = int(len(x)/3)
        Nf = len(fGHz)
        assert(len(x)==5*K)
        ak = x[0:1*K][:,None,None]
        bk = x[1*K:2*K][:,None,None]
        ck = x[2*K:3*K][:,None,None]
        tk = x[3*K:4*K][:,None,None]
        pk = x[4*K:5*K][:,None,None]
        # tf : paths (0) , freq (1), angle (2)
        if wH ==[]:
            wH = np.ones(len(fGHz))

        a = ak * bk*(fGHz[None,:,None]-fGHz[int(Nf/2)])**ck
        tf = a*np.exp(-2*1j*np.pi*fGHz[None,:,None]*tk)*wH[None,:,None]
        dphi = pk - phi[None,None,:]
        Gmax = 10**(GmaxdB/10.)

        g = np.exp(-(2*np.sqrt(np.log(2))*dphi/HPBW)**2)
        tfg = tf*g
        self.x = fGHz
        self.fcGHz = self.x[int(len(self.x)/2)]
        # self.y : angle(0) fGHz(1)
        self.y = np.sum(tfg,axis=0).T
        self.az = phi
        #h    = np.fft.ifft(H)


class ADPchannel(bs.TUsignal):
    """ Angular Delay Profile channel

    Attributes
    ----------

    az : array
        azimuth in radian
    ang_offset :
    theta : float
    phi : float
    tau : float
    _filename : string
        short filename for saving


    """
    def __init__(self,
            x = np.array([]),
            y = np.array([]),
            az = np.array([]),
            tx = np.array([]),
            rx = np.array([]),
            fcGHz=28,
            _filename='',
            refinement = False,
            ang_offset = 0,
            ):
        """
        Parameters
        ----------

        x : np.array
            delay
        y : np.array
            angle x delay
        az : np.array
            azimuth angle
        tx : np.array
            tx coordinates
        rx : np.array
            rx coordinates
        _filename :
        refinement : boolean
            False
        offset  :

        """
        bs.TUsignal.__init__(self, x=x, y=y,label='ADP')
        self.az = az
        self.tx = tx
        self.rx = rx
        self._filename = _filename
        self.fcGHz = fcGHz
        self.refinement = refinement
        self.ang_offset = ang_offset
        if ((len(self.tx) !=0 ) and (len(self.rx)!= 0)):
            v = self.tx - self.rx
            distLOS = np.linalg.norm(v)
            self.taulos_geo = distLOS/0.3
            self.anglos_geo = np.arctan2(v[1],v[0])*180/np.pi
            if self.anglos_geo<0:
                self.anglos_geo += 360
            LFS = -(32.4 + 20*np.log10(fcGHz) + 20*np.log10(distLOS))
            self.alphalos_geo = 10**(LFS/10.)

            if self.anglos_geo<0:
                self.anglos_geo = 2*np.pi+self.anglos_geo

        alphapeak,taupeak,angpeak = self.peak(refinement=refinement)
        self.angpeak_est = angpeak*180/np.pi
        self.taupeak_est = taupeak
        self.alphapeak_est = alphapeak


        self._filename = _filename

    def __repr__(self):
        cv = 180/np.pi
        s = 'Angular Delay Profile object \n'
        stx = "%.2f, %.2f, %.2f" % (self.tx[0],self.tx[1],self.tx[2])
        srx = "%.2f, %.2f, %.2f" % (self.rx[0],self.rx[1],self.rx[2])
        s = s + 'Tx : '+ stx + '\n'
        s = s + 'Rx : '+ srx + '\n'
        s = s + 'Angular offset (degrees) : '+str(cv*self.ang_offset)+'\n'
        s = s + 'agstart : '+str(self.az[0]*cv)+' agstop : '+str(self.az[-1]*cv)+'\n'
        if hasattr(self,'alphalos_geo'):
            alphalosdB = 10*np.log10(self.alphalos_geo)
            alphapeak_estdB = 20*np.log10(self.alphapeak_est)
            s = s + 'alpha (geo): '+ '%.2f dB' % alphalosdB+' (est): '+ '%.2f dB' % alphapeak_estdB + ' GdB : '+str(10*np.log10(self.alphapeak_est**2/self.alphalos_geo))+' \n'
            s = s + 'tau (geo): '+ str(self.taulos_geo)+' (est): '+str(self.taupeak_est)+' \n'
            s = s + 'ang (geo): '+ str(self.anglos_geo)+' (est): '+str(self.angpeak_est)+'\n'
        return(s)

    def peak(self, refinement=False):
        """ evaluate peak of PADP

        Parameters
        ----------
        refinment : boolean
            provide a refined version of angular estimation

        Returns
        -------

        alphapeak, taupeak , phipeak

        """
        alphapeak = np.max(np.abs(self.y))
        iphi, itau = np.where(np.abs(self.y)==alphapeak)
        taupeak = self.x[itau][0]
        if refinement:
            pr = np.abs(self.y)[iphi-1:iphi+2,itau].squeeze()
            azr = self.az[iphi-1:iphi+2]
            Id = np.sum(pr)
            In = np.sum(pr*azr)
            phipeak = In/Id
        else:
            phipeak = self.az[iphi]

        return alphapeak, taupeak, phipeak[0]

    def cut(self,imin=0,imax=1000):
        self.y = self.y[:,imin:imax]
        self.x = self.x[imin:imax]

    def correlate(self,adp,thresholddB=-105):
        """ correlate ADP with an other ADP

        Parameters
        ----------

        adp : ADPchannel

        Returns
        -------

        rhoE  : energy ratio of padp Eadp/Eself
        rhoEc : energy ratio of centered padp Ecadp/Ecself
        rho   : normalized intercorrelation  :  <self-mean(self),adp-mean(adp)>/Eself
        rhon  : intercorrelation of normalized padp   <self_normalized,adp_normalized>

        Notes
        -----
        This can be used to compare a measured PADP with a Ray tracing PADP


        """
        #import ipdb
        #ipdb.set_trace()
        #
        # apply the min dB level thresholding
        #
        tmp_self = np.abs(self.y)
        tmp_adp  = np.abs(adp.y)
        u1 = np.where(20*np.log10(tmp_self)>thresholddB)
        u2 = np.where(20*np.log10(tmp_adp)>thresholddB)
        padp_self = np.zeros(tmp_self.shape)
        padp_adp = np.zeros(tmp_adp.shape)

        padp_self[u1] = tmp_self[u1]
        padp_adp[u2] = tmp_adp[u2]

        padpc_self = padp_self-np.mean(padp_self)
        padpc_adp  = padp_adp-np.mean(padp_adp)

        Eself = np.max(si.correlate2d(padp_self,padp_self,mode='same'))
        Ecself = np.max(si.correlate2d(padpc_self,padpc_self,mode='same'))
        Eadp = np.max(si.correlate2d(padp_adp,padp_adp,mode='same'))
        Ecadp = np.max(si.correlate2d(padpc_adp,padpc_adp,mode='same'))
        #Eself =  np.sum(padp_self*padp_self)
        #Ecself = np.sum(padpc_self*padpc_self)
        #Eadp = np.sum(padp_adp*padp_adp)
        #Ecadp = np.sum(padpc_adp*padpc_adp)

        padpcn_self = padpc_self/np.sqrt(Ecself)
        padpcn_adp = padpc_adp/np.sqrt(Ecadp)

        rhoE = Eadp/Eself
        rhoEc = Ecadp/Ecself
        #rho  = np.sum(padpc_self*padpc_adp)/Eself
        #rhoc = np.sum(padpc_self*padpc_adp)/Ecself
        #rhon = np.sum(padpcn_self*padpcn_adp)
        rho  = np.max(si.correlate2d(padpc_self,padpc_adp,mode='same'))/Eself
        rhoc = np.max(si.correlate2d(padpc_self,padpc_adp,mode='same'))/Ecself
        rhon = np.max(si.correlate2d(padpcn_self,padpcn_adp,mode='same'))

        return rhoE,rhoEc,rho,rhoc,rhon

    def svd(self):
        """ perform singular value decomposition of the PADP

        Notes
        -----
        It creates a dictionnay
            {'sv':sv,'b':b}

        """
        [U,S,V]=la.svd(self.y)
        self.d = {}
        for k,sv in enumerate(S):
            b = sv*np.dot(U[:,k][:,None],V[k,:][None,:])
            self.d[k] = {'sv':sv,'b':b}

    def imshow(self,**kwargs):
        """ show Angular Delay Profile

        Parameters
        ----------
        origin: string
            'lower'
        vmax : -65,
        vmin : -120,
        interpolation : string 
             'nearest',
        alpha:1,
        imin = 0
        imax = -1
        dB   = True
        fig  = []
        ax   = []
        fonts = 18
        label = ''
        blos = True
        orientation = -1
        bcolorbar = False
        ang_offset = 450

        """
        defaults = {'origin':'lower',
                    'vmax' : -65,
                    'vmin' : -120,
                    'interpolation' : 'nearest',
                    'alpha':1,
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        imin = kwargs.pop('imin',0)
        imax = kwargs.pop('imax',-1)
        dB   = kwargs.pop('dB',True)
        fig  = kwargs.pop('fig',[])
        ax   = kwargs.pop('ax',[])
        fonts = kwargs.pop('fonts',18)
        label = kwargs.pop('label','')
        blos = kwargs.pop('blos',True)
        orientation = kwargs.pop('orientation',-1)
        bcolorbar = kwargs.pop('colorbar',False)
        ang_offset = kwargs.pop('ang_offset',450)

        if fig==[]:
            fig = plt.figure()
        if ax==[]:
            ax = fig.add_subplot(111)

        #rd2deg = 180/np.pi
        #extent = (self.az[-1]*rd2deg+agoffset,
        #          self.az[0]*rd2deg+agoffset,
        #          self.x[imin],self.x[imax])
        #extent = (self.az[0]*rd2deg,
        #          self.az[-1]*rd2deg,
        #          self.x[imin],self.x[imax])

        agmin = self.az.min()*180/np.pi
        agmax = self.az.max()*180/np.pi
        extent = (agmin,agmax,self.x[imin],self.x[imax])
        if orientation==-1:
            padp = np.abs(self.y)[::-1,imin:imax].T
        else:
            padp = np.abs(self.y)[:,imin:imax].T
        if dB:
            padp  = 20*np.log10(padp)
            pdb.set_trace()
            im = ax.imshow(padp,extent=extent,aspect='auto',**kwargs)
            #plt.axis('equal')

        if blos:
            a1 = ang_offset + self.angpeak_est
            ax.scatter(a1,self.taupeak_est,marker='*',s=70,color='r')
            if hasattr(self,'anglos_geo'):
                a2 = ang_offset + self.anglos_geo
                ax.scatter(a2,self.taulos_geo,marker='D',s=70,color='g')

        if bcolorbar:
            cbar = plt.colorbar(im)

            if dB:
                cbar.set_label(label+' dB',fontsize=fonts)
            else:
                cbar.set_label(label+' linear',fontsize=fonts)

            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fonts)

        ax.set_ylabel('Propagation delay [ns]',fontsize=fonts)
        ax.set_xlabel('Angle[deg]',fontsize=fonts)

        #ax.title('PADP',fontsize=fonts)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fonts) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fonts) 

        return fig,ax

    def clean(self,threshold_dB=20):
        """  clean ADP

        Parameters
        ----------

        threshold_dB : float

        Notes
        -----

        All values below Max -threshold are set to zero

        """
        Na = self.y.shape[0]
        P = np.real(self.y*np.conj(self.y))
        MaxdB = 10*np.log10(np.max(P))
        u = np.where(10*np.log10(P) < MaxdB-threshold_dB)

        self.y[u] = 0+0j

    def pap(self,
            fcGHz=28,
            fontsize=18,
            figsize=(10,10),
            Gmax=22.68,
            Gmin=19,
            threshdB=-95,
            label='',
            color='k',
            fig=[],
            ax=[],
            xlabel=True,
            ylabel=True,
            legend=True):
        """ Calculate Power Angular Profile

        Parameters
        ----------

        fcGHz : float
        fontsize : int
        figsize : tuple
        fig :
        ax  :
        xlabel : boolean
        ylabel : boolean
        legen : boolean

        Returns
        -------

        fig,ax

        """

        Na = self.y.shape[0]
        # integration over frequency
        # adp (angle)
        Gtyp = (Gmax+Gmin)/2.
        Py   = np.real(self.y*np.conj(self.y))
        pdp0 = np.sum(Py,axis=0)
        pdp0dB = 10*np.log10(pdp0)
        u = pdp0dB > threshdB
        adp  = np.sum(Py[:,u],axis=1)
        #mPya = np.median(Py,axis=0)
        #mPya = np.mean(Py,axis=0)
        #sPy = Py-mPya[None,:]
        #adp = np.sum(Pyz,axis=1)
        u = np.where(adp==max(adp))[0]
        if fig==[]:
            fig = plt.figure(figsize=figsize)
        else:
            fig = fig
        if ax == []:
            ax  = fig.add_subplot(111)
        else:
            ax = ax
        #ax.plot(self.az*180/np.pi,10*np.log10(adp),color='r',label=r'$10\log_{10}(\sum_{\tau} PADP(\phi,\tau))$',linewidth=1.5)
        #ag = np.linspace(45,260,len(adp))
        ag = self.az*180/np.pi
        ax.plot(ag, #360self.az*180/np.pi,
                10*np.log10(adp)-Gtyp,
                color=color,
                label=label,
                linewidth=1.5)
        ax.vlines(self.anglos_geo,ymin=-130,ymax=-40,linestyles='dashed',color='red')
        ax.hlines(-120,xmin=ag[0],xmax=ag[-1],linestyles='dashed',color='black')
        #ax.set_ylim(-80,-60)
        if xlabel:
            ax.set_xlabel('Angle [deg]',fontsize=fontsize)
        if ylabel:
            ax.set_ylabel('level (dB)',fontsize=fontsize)

        #ax.set_title(self._filename,fontsize=fontsize)
        if legend:
            plt.legend(loc='best')
        return fig,ax

    def app(self,**kwargs):
        """ Calculate Angular Power Profile
        """
        Na = self.y.shape[0]
        app = np.real(np.sum(self.y*np.conj(self.y),axis=1))

    def pltcir(self,phideg,Gain=21):
        """ plot Channel Impulse Response

        Parameters
        ----------

        phideg : f

        Returns
        -------

        fig,ax
        u

        """


        phi = phideg*np.pi/180.
        dang = np.abs(self.az - phi)
        u = np.where(dang==np.min(dang))[0][0]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        FS = -(32.4+20*np.log10(self.x*0.3)+20*np.log10(self.fcGHz))
        plt.semilogx(self.x,20*np.log10(np.abs(self.y[u,:]))-Gain)
        plt.semilogx(self.x,FS,'k',linewidth=2)
        plt.show()
        return fig,ax,u

    def pdp_v(self,**kwargs):
        """ Calculate and plot Power Delay Profile

        Parameters
        ----------

        fcGHz : float

        """

        defaults = { 'figsize':(10,10),
                     'fontsize':18,
                     'fig' : [],
                     'ax': [],
                     'xlabel': True,
                     'ylabel': True,
                     'legend': True,
                     'losdelay': True,
                     'freespace': True,
                     'desembeded': False,
                     'noisefloor': False,
                     'typic':True,
                     'semilogx':True,
                     'bcir':False,
                     'raw': False,
                     'Gmax':22.68,
                     'Gmin':19,
                     'threshdB':75,
                     'imax':-1,
                     'Tilt':10,
                     'HPBW':10,
                     'dphi':5,
                     'marker':'*',
                     'color':'k',
                     'label':'',
                     'linewidth':1
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        Gmax = kwargs.pop('Gmax')
        Gmin = kwargs.pop('Gmin')
        imax = kwargs.pop('imax')
        threshdB = kwargs.pop('threshdB')
        Gtyp = (Gmax+Gmin)/2.
        # get peak value of the PADP
        alpha,tau,phi = self.peak()
        Na = self.y.shape[0]
        # pdp : power delay profie 
        Py = np.real(self.y*np.conj(self.y))
        pap0 = np.sum(Py,axis=1)
        pap0dB = 10*np.log10(pap0)
        u = pap0dB>np.percentile(pap0dB,threshdB)
        pdp = np.sum(Py[u,:],axis=0)
        pdp = pdp[0:imax]
        x = self.x[0:imax]
        # spdp : square root of power delay profie 
        spdp = TUchannel(x=x,y=np.sqrt(pdp))
        u  = np.where(pdp==max(pdp))[0]
        FS = -(32.4+20*np.log10(x*0.3)+20*np.log10(self.fcGHz))
        AttmaxdB = 20*np.log10(alpha)
        #Gmax = AttmaxdB-FS[u] 
        #Gmax_r = np.round(Gmax[0]*100)/100.
        #
        # The -3dB is specific to the Aalto measurement and desembeding (1/2)
        #
        pdp_min = 10*np.log10(pdp)-Gmax-1
        pdp_max = 10*np.log10(pdp)-Gmin-1
        pdp_typ = 10*np.log10(pdp)-Gtyp-1
        uflashing = np.where(pdp_typ>FS)

        umin = np.where(pdp_min>-118)
        pdp_min_thr = pdp_min[umin]
        umax = np.where(pdp_max>-118)
        pdp_max_thr = pdp_max[umax]
        PL = -10*np.log10(np.sum(10**(pdp_min_thr/10.)))

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else:
            fig = kwargs['fig']
        if kwargs['ax'] == []:
            ax  = fig.add_subplot(111)
        else:
            ax = kwargs['ax']

        if kwargs['semilogx']:
            if kwargs['raw']:
                ax.semilogy(10*np.log10(pdp),x,color='r',label=r'$10\log_{10}(\sum_{\phi} PADP(\phi))$',linewidth=0.5)
            #ax.semilogx(np.array([tau]),np.array([AttmaxdB]),color='k')

            if kwargs['desembeded']:
                ax.semilogy(pdp_min,x,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmax),color='green')
                ax.semilogy(pdp_max,x,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmin),color='red')

            if kwargs['typic']:
                ax.semilogy(pdp_typ,x,label=kwargs['label'],color=kwargs['color'],linewidth=kwargs['linewidth'])
                ax.semilogy(pdp_typ[uflashing],x[uflashing],label=kwargs['label'],color='red',linewidth=kwargs['linewidth'])

            if kwargs['freespace']:
                if kwargs['typic']:
                    ax.semilogy(FS,x,color=kwargs['color'],linewidth=kwargs['linewidth']+1,label='Free Space path profile')
                else:
                    ax.semilogy(FS,x,color='k',linewidth=2,label='Free Space path profile')

            if kwargs['losdelay']:
                ax.hlines(self.taupeak_est,xmin=-130,xmax=-40,linestyles='dashed',color='blue')
                ax.hlines(self.taulos_geo,xmin=-130,xmax=-40,linestyles='dashed',color='red')

            if kwargs['noisefloor']:
                ax.vlines(-130,ymin=0,ymax=x[-1],linestyles='dashed',color='black')

            #ax.set_xlim(10,1000)
            if kwargs['xlabel']:
                ax.set_ylabel('Delay (ns) log scale',fontsize=kwargs['fontsize']) 

            if kwargs['bcir']:
                phi = self.angpeak_est*np.pi/180.
                dang = np.abs(self.az - phi)
                u = np.where(dang==np.min(dang))[0][0]
                ax.semilogx(20*np.log10(np.abs(self.y[u,:]))-Gmax,x,color='r')
                ax.semilogx(20*np.log10(np.abs(self.y[u,:]))-Gmin,x,color='g')
        else:
            if kwargs['raw']:
                ax.plot(10*np.log10(pdp),x,color='r',label=r'$10\log_{10}(\sum_{\phi} PADP(\phi))$',linewidth=0.5)
            ax.plot(np.array([AttmaxdB]),np.array([tau]),color='k')

            if kwargs['desembeded']:
                ax.plot(pdp_min,x,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmax))
                ax.plot(pdp_max,x,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmin))

            if kwargs['typic']:
                ax.plot(pdp_typ,x,label=kwargs['label'],color=kwargs['color'])
                ax.scatter(pdp_typ[uflashing],x[uflashing],s=80,c='red')

            if kwargs['freespace']:
                if kwargs['typic']:
                    ax.plot(FS,x,color=kwargs['color'],linewidth=kwargs['linewidth']+1,label='Free Space path profile')
                    ax.plot(FS-(Gmax-Gmin)/2,x,color='blue',linewidth=0.5,label='Free Space path profile')
                    ax.plot(FS+(Gmax-Gmin)/2,x,color='blue',linewidth=0.5,label='Free Space path profile')
                else:
                    ax.plot(FS,x,color='k',linewidth=2,label='Free Space path profile')

            if kwargs['losdelay']:
                ax.hlines(self.taupeak_est,xmin=-130,xmax=-40,linestyles='dashed',color='blue')
                ax.hlines(self.taulos_geo,xmin=-130,xmax=-40,linestyles='dashed',color='red')

            if kwargs['noisefloor']:
                ax.vlines(-130,ymin=0,ymax=x[-1],linestyles='dashed',color='red')
            #ax.set_xlim(0,1000)
            if kwargs['xlabel']:
                ax.set_ylabel('Delay (ns)',fontsize=kwargs['fontsize']) 

            if kwargs['bcir']:
                phi = self.angpeak_est*np.pi/180.
                dang = np.abs(self.az - phi)
                u = np.where(dang==np.min(dang))[0][0]
                ax.plot(20*np.log10(np.abs(self.y[u,:]))-Gmax,x,'r')
                ax.plot(20*np.log10(np.abs(self.y[u,:]))-Gmin,x,'g')

        if kwargs['ylabel']:
            ax.set_xlabel('level (dB)',fontsize=kwargs['fontsize']) 
        #ax.set_title(self._filename+' '+str(PL))
        if kwargs['legend']:
            plt.legend(loc='best') 

        ax.set_ylim(0,x[-1])
        return fig,ax

    def pdp(self,**kwargs):
        """ Calculate the Power Delay Profile

        Parameters
        ----------

        fcGHz : float
        figsize':(1010)
        fontsize':18
        fig' : []
        ax': []
        xlabel': True
        ylabel': True
        legend': True
        losdelay': True
        freespace': True
        desembeded': False
        typic':True
        semilogx':True
        bcir':False
        raw': False
        Gmax':22.68
        Gmin':19
        Tilt':10
        HPBW':10

        Returns
        -------

        tau
        pdp

        """

        defaults = { 'figsize':(10,10),
                     'fontsize':18,
                     'fig' : [],
                     'ax': [],
                     'xlabel': True,
                     'ylabel': True,
                     'legend': True,
                     'losdelay': True,
                     'freespace': True,
                     'desembeded': False,
                     'typic':True,
                     'semilogx':True,
                     'bcir':False,
                     'raw': False,
                     'bplot':True,
                     'Gmax':22.68,
                     'Gmin':19,
                     'Tilt':10,
                     'HPBW':10,
                     'dphi':5,
                     'marker':'*',
                     'color':'k',
                     'label':'',
                     'linewidth':1
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]
        # get antenna gain extremum
        # typical value is chosen as the mean value
        Gmax = kwargs.pop('Gmax')
        Gmin = kwargs.pop('Gmin')
        Gtyp = (Gmax+Gmin)/2.
        # get peak value of the PADP
        # it is assume that this retreave the LOS component
        alpha, tau, phi = self.peak()
        # Na : number of angular steps
        Na = self.y.shape[0]
        # pdp : power delay profile
        pdp = np.real(np.sum(self.y*np.conj(self.y),axis=0))
        # delay index of pdp maximum
        u = np.where(pdp==max(pdp))[0]
        # omnidirectional free space path loss
        # spdp : square root of power delay profile
        spdp = TUchannel(x=self.x,y=np.sqrt(pdp))
        u = np.where(pdp==max(pdp))[0]
        FS = -(32.4+20*np.log10(self.x*0.3)+20*np.log10(self.fcGHz))
        AttmaxdB = 20*np.log10(alpha)
        #Gmax = AttmaxdB-FS[u]
        #Gmax_r = np.round(Gmax[0]*100)/100.
        #
        # The -3dB is specific to the Aalto measurement and desembeding (1/2)
        #
        pdpdB = 10*np.log10(pdp)
        pdp_min = pdpdB-Gmax-1
        pdp_max = pdpdB-Gmin-1
        pdp_typ = pdpdB-Gtyp-1

        umin = np.where(pdp_min>-118)
        pdp_min_thr = pdp_min[umin]
        umax = np.where(pdp_max>-118)
        pdp_max_thr = pdp_max[umax]

        PL = -10*np.log10(np.sum(10**(pdp_min_thr/10.)))

        if kwargs['bplot']:
            if kwargs['fig']==[]:
                fig = plt.figure(figsize=kwargs['figsize'])
            else:
                fig = kwargs['fig']
            if kwargs['ax'] == []:
                ax  = fig.add_subplot(111)
            else:
                ax = kwargs['ax']

            if kwargs['semilogx']:
                if kwargs['raw']:
                    ax.semilogx(self.x,10*np.log10(pdp),color='r',label=r'$10\log_{10}(\sum_{\phi} PADP(\phi))$',linewidth=0.5)
                #ax.semilogx(np.array([tau]),np.array([AttmaxdB]),color='k')

                if kwargs['desembeded']:
                    ax.semilogx(self.x,pdp_min,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmax),color='green')
                    ax.semilogx(self.x,pdp_max,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmin),color='red')

                if kwargs['typic']:
                    ax.semilogx(self.x,pdp_typ,label=kwargs['label'],color=kwargs['color'],linewidth=kwargs['linewidth'])

                if kwargs['freespace']:
                    if kwargs['typic']:
                        ax.semilogx(self.x,FS,color=kwargs['color'],linewidth=kwargs['linewidth']+1,label='Free Space path profile')
                    else:
                        ax.semilogx(self.x,FS,color='k',linewidth=2,label='Free Space path profile')

                if kwargs['losdelay']:
                    ax.vlines(self.taupeak_est,ymin=-130,ymax=-40,linestyles='dashed',color='blue')
                    ax.vlines(self.taulos_geo,ymin=-130,ymax=-40,linestyles='dashed',color='red')

                #ax.set_xlim(10,1000)
                if kwargs['xlabel']:
                    ax.set_xlabel('Delay (ns) log scale',fontsize=kwargs['fontsize'])

                if kwargs['bcir']:
                    phi = self.angpeak_est*np.pi/180.
                    dang = np.abs(self.az - phi)
                    u = np.where(dang==np.min(dang))[0][0]
                    ax.semilogx(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmax,color='r')
                    ax.semilogx(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmin,color='g')
            else:
                if kwargs['raw']:
                    ax.plot(self.x,10*np.log10(pdp),color='r',label=r'$10\log_{10}(\sum_{\phi} PADP(\phi))$',linewidth=0.5)
                ax.plot(np.array([tau]),np.array([AttmaxdB]),color='k')

                if kwargs['desembeded']:
                    ax.plot(self.x,pdp_min,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmax))
                    ax.plot(self.x,pdp_max,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmin))

                if kwargs['typic']:
                    ax.plot(self.x,pdp_typ,label=kwargs['label'],color=kwargs['color'])

                if kwargs['freespace']:
                    if kwargs['typic']:
                        ax.plot(self.x,FS,color=kwargs['color'],linewidth=kwargs['linewidth']+1,label='Free Space path profile')
                    else:
                        ax.plot(self.x,FS,color='k',linewidth=2,label='Free Space path profile')

                if kwargs['losdelay']:
                    ax.vlines(self.taupeak_est,ymin=-130,ymax=-40,linestyles='dashed',color='blue')
                    ax.vlines(self.taulos_geo,ymin=-130,ymax=-40,linestyles='dashed',color='red')

                #ax.set_xlim(0,1000)
                if kwargs['xlabel']:
                    ax.set_xlabel('Delay (ns)',fontsize=kwargs['fontsize']) 

                if kwargs['bcir']:
                    phi = self.angpeak_est*np.pi/180.
                    dang = np.abs(self.az - phi)
                    u = np.where(dang==np.min(dang))[0][0]
                    ax.plot(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmax,'r')
                    ax.plot(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmin,'g')

            if kwargs['ylabel']:
                ax.set_ylabel('level (dB)',fontsize=kwargs['fontsize']) 
            ax.set_title(self._filename+' '+str(PL))
            if kwargs['legend']:
                plt.legend(loc='best')

            return fig,ax
        else:
            return (self.x,pdp)

#        PL = -10*np.log10(np.sum(10**(pdp_min_thr/10.)))
#        return self.x,pdp

#        if kwargs['fig']==[]:
#            fig = plt.figure(figsize=kwargs['figsize'])
#        else:
#            fig = kwargs['fig']
#        if kwargs['ax'] == []:
#            ax  = fig.add_subplot(111)
#        else:
#            ax = kwargs['ax']
#
#        if kwargs['semilogx']:
#            if kwargs['raw']:
#                ax.semilogx(self.x,10*np.log10(pdp),color='r',label=r'$10\log_{10}(\sum_{\phi} PADP(\phi))$',linewidth=0.5)
#            #ax.semilogx(np.array([tau]),np.array([AttmaxdB]),color='k')
#
#            if kwargs['desembeded']:
#                ax.semilogx(self.x,pdp_min,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmax),color='green')
#                ax.semilogx(self.x,pdp_max,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmin),color='red')
#
#            if kwargs['typic']:
#                ax.semilogx(self.x,pdp_typ,label=kwargs['label'],color=kwargs['color'],linewidth=kwargs['linewidth'])
#
#            if kwargs['freespace']:
#                if kwargs['typic']:
#                    ax.semilogx(self.x,FS,color=kwargs['color'],linewidth=kwargs['linewidth']+1,label='Free Space path profile')
#                else:
#                    ax.semilogx(self.x,FS,color='k',linewidth=2,label='Free Space path profile')
#
#            if kwargs['losdelay']:
#                ax.vlines(self.taupeak_est,ymin=-130,ymax=-40,linestyles='dashed',color='blue')
#                ax.vlines(self.taulos_geo,ymin=-130,ymax=-40,linestyles='dashed',color='red')
#
#            #ax.set_xlim(10,1000)
#            if kwargs['xlabel']:
#                ax.set_xlabel('Delay (ns) log scale',fontsize=kwargs['fontsize']) 
#
#            if kwargs['bcir']:
#                phi = self.angpeak_est*np.pi/180.
#                dang = np.abs(self.az - phi)
#                u = np.where(dang==np.min(dang))[0][0]
#                ax.semilogx(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmax,color='r')
#                ax.semilogx(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmin,color='g')
#        else:
#            if kwargs['raw']:
#                ax.plot(self.x,10*np.log10(pdp),color='r',label=r'$10\log_{10}(\sum_{\phi} PADP(\phi))$',linewidth=0.5)
#            ax.plot(np.array([tau]),np.array([AttmaxdB]),color='k')
#
#            if kwargs['desembeded']:
#                ax.plot(self.x,pdp_min,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmax))
#                ax.plot(self.x,pdp_max,label=r'$10\log_{10}(\sum_{\phi} PADP(\phi)) - $'+str(Gmin))
#
#            if kwargs['typic']:
#                ax.plot(self.x,pdp_typ,label=kwargs['label'],color=kwargs['color'])
#
#            if kwargs['freespace']:
#                if kwargs['typic']:
#                    ax.plot(self.x,FS,color=kwargs['color'],linewidth=kwargs['linewidth']+1,label='Free Space path profile')
#                else:
#                    ax.plot(self.x,FS,color='k',linewidth=2,label='Free Space path profile')
#
#            if kwargs['losdelay']:
#                ax.vlines(self.taupeak_est,ymin=-130,ymax=-40,linestyles='dashed',color='blue')
#                ax.vlines(self.taulos_geo,ymin=-130,ymax=-40,linestyles='dashed',color='red')
#
#            #ax.set_xlim(0,1000)
#            if kwargs['xlabel']:
#                ax.set_xlabel('Delay (ns)',fontsize=kwargs['fontsize']) 
#
#            if kwargs['bcir']:
#                phi = self.angpeak_est*np.pi/180.
#                dang = np.abs(self.az - phi)
#                u = np.where(dang==np.min(dang))[0][0]
#                ax.plot(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmax,'r')
#                ax.plot(self.x,20*np.log10(np.abs(self.y[u,:]))-Gmin,'g')
#
#        if kwargs['ylabel']:
#            ax.set_ylabel('level (dB)',fontsize=kwargs['fontsize']) 
#        ax.set_title(self._filename+' '+str(PL))
#        if kwargs['legend']:
#            plt.legend(loc='best')
#
#        return fig,ax

    def tomap(self,L,**kwargs):
        """ surimpose PADP on the Layout

        Parameters
        ----------

        L : Layout
        xmin : 10
        xmax : 400
        ymin : 10
        ymax : 400,
        Nx  :3000,
        Ny  :3000,
        'cmap':'jet',
        'mode':'image',
        'excess':'los',
        'figsize':(20,20),
        'thmindB':-110,
        'thmaxdB':-108,
        'vmindB':-110,
        'vmaxdB':-60,
        'offset':0,
        'display':True,
        'compensated':True,
        'tauns_excess':0


        """
        xmin = kwargs.pop('xmin',0)
        ymin = kwargs.pop('ymin',0)
        xmax = kwargs.pop('xmax',20)
        ymax = kwargs.pop('ymax',20)
        mode = kwargs.pop('mode','sbounce')
        vmindB = kwargs.pop('vmindB',-110)
        vmaxdB = kwargs.pop('vmaxdB',-60)
        thmindB = kwargs.pop('thmindB',-110)
        thmaxdB = kwargs.pop('thmaxdB',-108)
        Nx = kwargs.pop('Nx',3000)
        Ny = kwargs.pop('Ny',3000)
        cmap = kwargs.pop('cmap','jet')
        offset = kwargs.pop('offset',0)
        excess = kwargs.pop('excess','los')
        display = kwargs.pop('display',True)
        compensated = kwargs.pop('compensated',False)
        tauns_excess = kwargs.pop('tauns_excess',0)
        figsize = kwargs.pop('figsize',(20,20))

        if 'fig' not in kwargs:
            fig = plt.figure(figsize=figsize)
        else:
            fig = kwargs['fig']

        if 'ax' not in kwargs:
            ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']

        #
        # Prepare the array for spatial information in horizontal plane x,y
        # Nx and Ny should be large enough
        #
        Z  = np.zeros((Nx,Ny),dtype=complex)

        #
        # spatial indexation in x and  y
        #
        xr = np.linspace(xmin,xmax,Nx)
        yr = np.linspace(xmin,xmax,Ny)

        # distance Tx Rx in the horizontal plane (2D)
        dtx_rx_2D = np.sqrt((self.tx[0]-self.rx[0])**2+(self.tx[1]-self.rx[1])**2)
        # distance Tx Rx in the horizontal plane (3D)
        dtx_rx = np.sqrt((self.tx[0]-self.rx[0])**2+(self.tx[1]-self.rx[1])**2+(self.tx[2]-self.rx[2])**2)
        # distance  Tx ground Rx (3D)
        dtx_gr_rx = np.sqrt(dtx_rx_2D**2+(self.tx[2]+self.rx[2])**2)
        assert(dtx_gr_rx > dtx_rx)
        # difference of heights beween Tx and Rx
        deltah = np.abs(self.tx[2]-self.rx[2])

        #
        # Dt = vec(P,Tx)
        # Dr = vec(Rx,P)
        #
        dxt =(self.tx[0]-xr)[:,None]
        dyt =(self.tx[1]-yr)[None,:]
        #
        # nwt : distance between Tx and each point of the plane
        # nwr : distance between Rx and each point of the plane
        #
        nwt = np.sqrt(dxt*dxt+dyt*dyt)

        dxr =(xr-self.rx[0])[:,None]
        dyr =(yr-self.rx[1])[None,:]
        nwr = np.sqrt(dxr*dxr+dyr*dyr)

        # dsbounce : elliposidal distance (single bounce hypothesis)
        dsbounce = nwt+nwr

        # maximal ellipsoidal distance on the Z selected region
        dmax = dsbounce.max()
        taumax = dmax/0.3

        # determine index of maximal distance
        if self.x.max()>taumax:
            itaumax = np.where(self.x>taumax)[0][0]
        else:
            itaumax=len(self.x)-1
        # convert maximal distance into maximal delay (self.x is delay)
        taumax = self.x[itaumax]
        # determine coefficient  between delay and index ( ns --> integer)
        tau2idx = taumax/itaumax

        # Determine the angle of arrival
        # direction of arrival normalization of the vector

        dxrn = dxr/nwr
        dyrn = dyr/nwr

        # angle of arrival in [-pi,pi]
        phi = np.arctan2(dyrn,dxrn)-offset*np.pi/180
        # back in [0-2pi]
        phi = (1-np.sign(phi))*np.pi+phi

        #iphi=((315-phi*180/np.pi)/5).astype(int)
        iphi=((360-phi*180/np.pi)/5).astype(int)

        drpt = np.sqrt(dxr*dxr+dyr*dyr+dxt*dxt+dyt*dyt)
        dpr = np.sqrt(dxr*dxr+dyr*dyr)

        if mode=='sbounce':
            iid = np.round((np.sqrt(dxt*dxt+dyt*dyt)+np.sqrt(dxr*dxr+dyr*dyr))/(0.3*tau2idx)).astype('int')
        else:
            #d = np.round(np.sqrt(dxr*dxr+dyr*dyr)/(0.3*0.625)).astype('int')
            #d = np.round(np.sqrt(dxr*dxr+dyr*dyr)/(0.3*0.625)).astype('int')
            alpha = np.arctan(deltah/drpt)
            dv = dpr/np.cos(alpha)
            iid = np.round(dv/(0.3*tau2idx)).astype('int')
            #pdb.set_trace()

        #
        # create indexation for spatial region Z
        #
        ix = np.arange(Nx)[:,None]
        iy = np.arange(Ny)[None,:]

        # ird : index for delays (d for delays)
        ird = iid[ix,iy].ravel()
        # irp : index for directio of arrival (p for phi)
        irp = iphi[ix,iy].ravel()
        #
        #  (d < dmax ) and (d>dlos+tauns_excess)
        #  iphi >= 0 and iphi < Nphimax
        ilos    = np.round((dtx_rx/(0.3*tau2idx))).astype(int)
        iground = np.round((dtx_gr_rx/(0.3*tau2idx))).astype(int)
        iexcess = np.round(tauns_excess/tau2idx).astype(int)

        if excess=='los':
            ud = np.where((ird<itaumax) & (ird>ilos+iexcess))

        if excess=='ground':
            ud = np.where((ird<itaumax) & (ird>iground+iexcess))

        up = np.where((irp>=0) & (irp<len(self.az)))

        # determine the index of points in a corona wich satisfy jointly the
        # condition on delays and angles
        #
        u = np.intersect1d(ud,up)

        # ravelize Z (2D -> 1D)
        rz = Z.ravel()
        # filling rz with self.y nphi,Ntau
        rz[u] = self.y[irp[u],ird[u]]

        #
        # back to matrix form
        #

        Z = rz.reshape(Nx,Ny)
        lmbda = 0.3/self.fcGHz
        sqG = 10

        Z_compensated = Z*(4*np.pi*dtx_rx)/(sqG*lmbda)

        if compensated:
            ZdB = 20*np.log10(np.abs(Z_compensated.T))
        else:
            ZdB = 20*np.log10(np.abs(Z.T))

        mask = ((ZdB.all()>thmindB) and (ZdB.all()<thmaxdB))
        #mzdB = ma.masked_array(ZdB,mask)

        ZdBmax = ZdB.max()
        ZdBmin = ZdB.min()
        #
        # constructing figure
        #
        if display:
            #fig=plt.figure(figsize=figsize)
            fig,ax = L.showG('s', fig=fig, ax=ax, labels=0)
            #plt.axis('on')
            ax.imshow(ZdB, extent=(xr[0],xr[-1],yr[0],yr[-1]),
                    cmap = cmap,
                    origin = 'lower',
                    alpha = 0.9,
                    vmin = ZdBmax - 60,
                    vmax = ZdBmax, interpolation = 'nearest')
            #plt.imshow(mzdB,alpha=0.9,origin='lower')
            ax.plot(self.tx[0],self.tx[1],'og')
            ax.plot(self.rx[0],self.rx[1],'ob')
            #plt.colorbar()
            ax.set_title(self._filename)
            #plt.savefig(self._filename+'.png')

        #return Z,np.linspace(xr[0],xr[-1],Nx),np.linspace(yr[0],yr[-1],Ny)
        return fig,ax

    def polarplot(self,**kwargs):
        """  polar plot of PADP

        Parameters
        -----------

        fig
        ax
        figsize
        typ : string
        Ndec : int
            decimation factor (1)
        imax : int
            max value 150
        vmin : float
            -120
        vmax : float
            -50
        cmap : colormap
        title : PADP

        Returns
        -------
        fig , ax , pc (colormash) 


        """
        defaults = { 'fig':[],
                     'ax':[],
                     'figsize':(10,10),
                     'typ':'l20',
                     'Ndec':1,
                     'vmin':-120,
                     'vmax':-50,
                     'imax':150,
                     'alpha':1.,
                     'bcolorbar':True,
                     'cmap': plt.cm.jet,
                     'title':'PADP'
                   }

        cvel = 0.3
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if kwargs['fig'] == []:
            fig = plt.figure(figsize=kwargs['figsize'])
        else:
            fig = kwargs.pop('fig')
        if kwargs['ax'] == []:
            ax = fig.add_subplot(111,polar=True)
        else:
            ax = kwargs.pop('ax')

        imax = kwargs.pop('imax')
        Ndec = kwargs.pop('Ndec')
        vmin = kwargs.pop('vmin')
        vmax = kwargs.pop('vmax')
        cmap = kwargs.pop('cmap')
        alpha = kwargs.pop('alpha')
        title = kwargs.pop('title')

        rho,theta = np.meshgrid(self.x*cvel,self.az)
        # convert y data in desired format
        dt,ylabels = self.cformat(**kwargs)
        val = dt[:,0::Ndec][:,0:int(imax/Ndec)]
        th  = theta[:,0::Ndec][:,0:int(imax/Ndec)]
        rh  = rho[:,0::Ndec][:,0:int(imax/Ndec)]
        #vmin = np.min(val)
        #vmax = np.max(val)
        #Dynamic = max_val-vmin
        pc  = ax.pcolormesh(th,rh,val,cmap=cmap,vmin=vmin, vmax=vmax, alpha=alpha)
        #ptx = ax.plot(self.az,self.x*cvel,'or')
        if kwargs['bcolorbar']:
            fig.colorbar(pc,orientation='horizontal')
        ax.set_title(title)
        return fig,ax,pc
        #ax.axis('equal')
        #ax.axis('equal')


    def toafp(self,fmin):
        """ angular delay profile -> angular frequency profile
        """
        x = np.linspace(0,(len(self.x)-1)/(self.x[-1]-self.x[0]),len(self.x))+fmin
        y = np.fft.fft(self.y,axis=1)
        afp = AFPchannel(x=x,
                         y=y,
                         az=self.az,
                         tx=self.tx,
                         rx=self.rx,
                         _filename = self._filename,
                         refinement = self.refinement,
                         ang_offset = self.ang_offset)

        return afp

class TBchannel(bs.TBsignal):
    """ radio channel in non uniform delay domain
    """
    def __init__(self, x=np.array([]), y=np.array([]),label=[]):
        #super(TUsignal,self).__init__(x,y,label)
        bs.TBsignal.__init__(self,x,y,label)


    def tau_Emax(self):
        """ calculate the delay of max energy peak

        .. math::
            \max_{\tau} y^{2}(\tau)
        """
        y2 = (self.y) ** 2
        maxy2 = max(y2)
        u = np.nonzero(y2 == maxy2)[0]
        tau_Emax = self.x[u]
        return(tau_Emax)

    def tau_moy(self, alpha=0.1, threshold_dB = 20, tau0=0):
        """ calculate mean excess delay starting from delay tau0

        Parameters
        ----------

        alpha : float
        tau0 : float

        """

        u = np.max(self.y*self.y)
        v = 10**(np.log10(u)-threshold_dB/10.)

        uf = np.where(self.y*self.y > v)

        num = np.sum(self.y[uf]*self.y[uf]*self.x[uf[-1]])
        den = np.sum(self.y[uf]*self.y[uf])
        taum = num/den


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
        of = 1 - sum(y4,axis=0)/sum(y2,axis=0)**2

        return taum,delayspread,of


    def Kfactor(self,threshold_dB=20,dB=True):
        """ determine Ricean K factor 

        Parameters
        -----------

        Threshold_dB : float 
            Only the energy above threshold is taken into account
        dB : boolean 
            if True value in dB is returned

        """
        t = self.x
        y = self.y
        
        u  = np.max(self.y*self.y)
        v  = 10**(np.log10(u)-threshold_dB/10.)
        vmax = np.where(self.y*self.y==u)
        Pmax = self.y[vmax]*self.y[vmax]
        uf   = np.where(self.y*self.y>v)
        Ptot = np.sum(self.y[uf]*self.y[uf])
        K = Pmax/(Ptot-Pmax)
        if dB:
            K=10*np.log10(K)
        return K[0]

    def tau_rms(self, alpha=0.1,threshold_dB=20, tau0=0):
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

        See Also
        --------

        TUsignal.ecdf
        TUsignal.tau_moy

        """

        t = self.x
        y = self.y
        #cdf, vary = self.ecdf()
        #pdp = np.diff(cdf.y)

        u = np.max(self.y*self.y)
        v = 10**(np.log10(u)-threshold_dB/10.)

        uf = np.where(self.y*self.y>v)

        taum = self.tau_moy(tau0,threshold_dB=threshold_dB)

        num = np.sum(self.y[uf]*self.y[uf]*(self.x[uf[-1]]-taum)**2)
        den = np.sum(self.y[uf]*self.y[uf])

        taurms = np.sqrt(num/den)

        return taurms

    def toFD(self,fGHz=np.linspace(2,5,256)):
        """ Transform to Frequency domain

        Parameters
        ----------
        fGHz : ,Nf
            frequency in GHz

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
        >>> f,a = C.stem()

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
                #print(N)

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
        print(alpha)
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
        print(alpha)
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
        return P

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
            >>> # CIR = TUsignal(tau,y)
            >>> # CIR.aggcir(alphak,tauk)
            >>> # f,a =CIR.plot(typ=['v'])

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
           self.y = ynew[None,:]
        self.y = np.delete(self.y,0,0)

    def readcir(self,filename,outdir=[]):
        """ read channel impulse response

        Parameters
        ----------
        filename : string
            long file name if outdir is []
            short file name is outdir is != []
        outdir : string
            output directory
        """
        if outdir != []:
            outdir = 'output/'+outdir
            filename = getlong(filename, outdir)

        cir = ios.loadmat(filename)
        self.x = cir['t'].ravel()
        self.y = cir['cir'].ravel()


    def readuwb(self, _filename):
        """ read  Waveform from Matlab file

        Parameters
        ----------

        _filename : file name with extension (.mat)

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
        s1 = "Time domain channel with delay \n"
        s = TUchannel.__repr__(self)
        s = s1+s
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

class Mchannel(bs.FUsignal):
    """ Handle the measured channel

    """
    def __init__(self,
                x ,
                y ,
                **kwargs):
        """ class constructor

        Parameters
        ----------

        x  :  , nfreq
            frequency GHz
        y  :  Nm x Nr x Nt x Nf
            measured channel

        """
        defaults = {
        'Aat': [],
        'Aar': [],
        'calibrated':True,
        'label' :'',
        'filename':'',
        'mes':''
        }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.calibrated = kwargs.pop('calibrated')
        self.label = kwargs.pop('label')
        self.filename = kwargs.pop('filename')
        self.mes = kwargs.pop('mes')
        self.Aat = kwargs.pop('Aat')
        self.Aar = kwargs.pop('Aar')
        sh = y.shape
        self.Nm = sh[0]
        self.Nr = sh[1]
        self.Nt = sh[2]
        self.Nf = sh[3]

        bs.FUsignal.__init__(self,x=x,y=y,label='Mchannel')


    def __repr__(self):
        st = bs.FUsignal.__repr__(self)
        if self.calibrated:
            st = st + 'Calibrated'
        else:
             st = st + 'Not calibrated'
        return(st)



    def eig(self,HdH=False):
        """ calculate eigen values of the transfer matrix.
            it involves H and Hd against svd() which acts only over H.

        Returns
        -------

        HdH : Hermitian transfer matrix  (nf x nt x nt )
        U   : Unitary tensor  (nf x nt x nt )
        S   : Singular values (nf x nt)
        V   : = Ud (in that case because HdH Hermitian)  (nf x nt x nt)

        HdH = U L U^{\dagger}

        """

        # H  : nm x nr x nt x nf
        H   = self.y
        # Hd : nm x nt x nr x nf
        Hd  = np.conj(self.y.swapaxes(1,2))

        if HdH:
            #T : nm x nt x nt x nf
            T = np.einsum('uijk,ujlk->uilk',Hd,H)
        else:
            #T : nm x nr x nr x nf
            T = np.einsum('uijk,ujlk->uilk',H,Hd)
        # HdH : nm x nf x nr x nr
        T  = T.swapaxes(1,3)
        #U   : nm x nf x (nr|nt) x (nr|nt)
        #S   : nm x nf x (nr|nt)
        #V   : nm x nf x (nr|nt) x (nr|nt)
        U,S,V  = la.svd(T)


        return (U,S,V)

    def Bcapacity(self,Pt=np.array([1e-3]),Tp=273):
        """ calculates BLAST deterministic MIMO channel capacity

        Parameters
        ----------

        Pt : np.array (,NPt)
            the total power is assumed uniformaly distributed over the whole bandwidth
        Tp : Receiver Temperature (K)

        Returns
        -------

        C   : sum rate or spectral efficiency (bit/s)
            np.array (Nf,NPt)
        rho : SNR
            np.array (Nf,Nt,NPt)

            log_2(det(I+(Et/(N0Nt))HH^{H})

        Notes
        -----

        The returned value is homogeneous to bit/s the aggregated capacity is
        obtrained by a simple summation
        of the returned quantity. To obtain the sum rate or the spectral
        efficiency in (bit/s/Hz ) the returned value should be divided by the
        frequency step dfGHz

        """

        fGHz  = self.x
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

        U,S,V = self.eig(HdH=True)


        Ps  = Pt/(self.Nt)


        Pb  = N0*BGHz*1e9   # Watt


        # S : nm x nf x nr
        # rho : nm x nf x nr x power
        #
        rho  = (Ps[None,None,None,:]/Pb)*S[:,:,:,None]

        CB   = dfGHz*np.sum(np.log(1+rho)/np.log(2),axis=2)

        return(rho,CB)

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

        fGHz  = self.x
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

        # Evaluation of the transfer HHd tensor
        
        U,ld,V = self.eig(HdH=True)
       
        #
        # Iterative implementation of Water Filling algorithm
        #

        # pb : (nm,nf,nt)   noise power (Watt)
        pb = N0*dfGHz*1e9*np.ones((self.Nm,self.Nf,self.Nt))
        # pt : (nm,nf,nt,power)  Total power uniformly spread over (nt*nf-1)
        pt = Pt[None,None,None,:]/((self.Nf-1)*self.Nt)
        mu = pt
        Q0 = np.maximum(0,mu-pb[:,:,:,None]/ld[:,:,:,None])
        u  = np.where(Q0>0)[0]

        Peff = np.sum(np.sum(Q0,axis=1),axis=1)
        deltamu = pt
        while (np.abs(Peff-Pt)>1e-16).any():
            mu = mu + deltamu
            Q = np.maximum(0,mu-pb[:,:,:,None]/ld[:,:,:,None])
            Peff = np.sum(np.sum(Q,axis=1),axis=1)
            #print "mu , Peff : ",mu,Peff
            usup = np.where(Peff>Pt)[0]
            mu[:,:,:,usup] = mu[:,:,:,usup]- deltamu[:,:,:,usup]
            deltamu[:,:,:,usup] = deltamu[:,:,:,usup]/2.
        Qn   = Q/pb[:,:,:,None]
        rho  = Qn*ld[:,:,:,None]

        Cwf  = dfGHz*np.sum(np.log(1+rho)/np.log(2),axis=2)


        return(rho,Cwf)

    def plot2(self,fig=[],ax=[],mode='time'):

        if fig ==[]:
            fig = plt.gcf()
        if ax ==[]:
            ax = plt.gca()

        if mode=='time':
            cir = self.ift(ffts=1)
            y = cir.y
            x = cir.x
        else:
            y = self.y
            x = self.x

        my = np.mean(np.abs(y),axis=0)
        yc = np.abs(y)-my[None,...] # TD centered      ; TDc.shape : (85, 4, 8, 801)

        yc2 = np.abs(yc)**2 # square of TD centered     ; TDc2.shape : (85, 4, 8, 801)
        vary = np.mean(yc2,axis=0) #variance of TD  ; varTD.shape : (4, 8, 801)

        cpt = 0
        for r in range(self.Nr):
            for t in range(self.Nt):
                #cpt=cpt+1
                #ax = plt.subplot(self.Nr,self.Nt,cpt)
                #l1, = ax.plot(self.x,np.sqrt(vary[r,t,:]),color='k',linewidth=1,alpha=1)
                #l1, = ax.plot(self.x,np.sqrt(vary[r,t,:]),linewidth=1,alpha=1)
                #l2, = ax.plot(self.x,my[r,t,:],color='r',linewidth=1,alpha=1)
                l2, = ax.plot(x,my[r,t,:],linewidth=1,alpha=1)
                ticksx = ax.axes.get_xticklabels()
                ticksy = ax.axes.get_yticklabels()
                plt.setp(ticksx, visible=True)
                plt.setp(ticksy, visible=True)

                if (r == 0) & (t==1):
                    #l1, = ax.plot(self.x,np.sqrt(vary[r,t,:]),color='k',label='sd',linewidth=1,alpha=1)
                    l2, = ax.plot(x,np.abs(my[r,t,:]),color='r',label='mean',linewidth=1,alpha=1)
                if (r == 3) & (t==0):
                    plt.setp(ticksx, visible=True)
                    ax.axes.set_xticks(np.arange(x[0],x[-1],0.2))
                    plt.setp(ticksy, visible=True)
                if (r == 0) & (t==3):
                    plt.title(r'Evolution of the mean and the standard deviation of $\mathbf{H}(f)$',fontsize=12)
                if (r == 1) & (t==0):
                    ax.axes.set_ylabel('Amplitude (linear scale $\in [0,1]$)',fontsize=15)
                if (r == 3) & (t == 3):
                    ax.axes.set_xlabel('Frequency (GHz)',fontsize=15)




class Tchannel(bs.FUsignal):
    """ Handle the transmission channel

    The transmission channel TChannel is obtained through combination of the propagation
    channel and the antenna transfer functions from both transmitter and receiver.
    This channel contains all the spatial information for each individual ray.
    Warning : This is a frequency domain channel deriving from bs.FUsignal

    Attributes
    ----------

        ray transfer functions  (nray,nfreq)
    dod  :
        direction of depature (rad) [theta_t,phi_t]  nray x 2
    doa  :
        direction of arrival (rad)  [theta_r,phi_r]  nray x 2
    tau  :
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
        self.filcal = "calibration.mat"
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
                position of point a (transmitter)
            b = np.ndarray
                position of point b (receiver)
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
                print('Warning : H/'+grpname +'already exists in '+filenameh5)
            f=fh5['H/'+grpname]

            for k,va in self.__dict__.items():
                #print(k,va)
                f.create_dataset(k,shape = np.shape(va),data=va)
            fh5.close()
        except:
            fh5.close()
            raise NameError('Channel Tchannel: issue when writting h5py file')

    def _loadh5(self,filenameh5,grpname,**kwargs):
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


    def apply(self, W=[]):
        """ apply FUsignal W to the Tchannel

        Parameters
        ----------

        W :  Bsignal.FUsignal

        It exploits multigrid convolution from Bsignal.


        Returns
        -------

        V : FUDAsignal

        Notes
        -----

        Returns :math:`W(f) H_k(f)`

            + W may have a more important number of points and a smaller frequency band.
            + If the frequency band of the waveform exceeds the one of the
            transmission channel, a warning is sent.
            + W is a FUsignal whose shape doesn't need to be homogeneous with FUChannel H

        """
        if W!=[]:
            U = W * self
        else:
            U = self

        V = Tchannel(x= U.x, y = U.y, tau = self.taud, dod = self.dod, doa= self.doa)

        return(V)


    def applywav(self, Wgam=[]):
        """ apply waveform (time domain ) to obtain the rays impulses response

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


    def getcir(self,BWGHz=1,Nf=40000,fftshift=False):
        """ get the channel impulse response

        Parameters
        ----------

        BWGHz : Bandwidth
        Nf    : Number of frequency points
        fftshift : boolean

        See Also
        --------

        pylayers.simul.link.DLink.plt_cir

        """
        fGHz  = np.linspace(0,BWGHz,Nf)
        dfGHz = fGHz[1]-fGHz[0]
        tauns = np.linspace(0,1/dfGHz,Nf)
        # E : r x nr x nt x f
        E    = np.exp(-2*1j*np.pi*self.taud[:,None,None,None]*fGHz[None,None,None,:])
        # self.y : r x nr x nt x f
        if self.y.shape[3]==E.shape[3]:
            H    = np.sum(E*self.y,axis=0)
        else:
            if self.y.shape[3]==1:
                H    = np.sum(E*self.y,axis=0)
            else:
                H    = np.sum(E*self.y[:,:,:,0][:,:,:,None],axis=0)
        # back in time - last axis is frequency (axis=2)
        cir  = np.fft.ifft(H,axis=2)
        if fftshift:
            cir = np.fft.fftshift(cir,axes=2)
            tauns = np.linspace(-Nf/(2*BWGHz),Nf/(2*BWGHz)-1/BWGHz,Nf)

        cir = bs.TUsignal(x=tauns,y=cir)

        return(cir)

    def get_cir(self,Wgam=[]):
        """ get Channel impulse response of the channel
            for a given waveform

        Parameters
        ----------

        Wgam : waveform

        Returns
        -------

        ri  : TUsignal

            impulse response for each ray separately

        See Also
        --------

        pylayers.antprop.channel.rir

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
        print(DeprecationWarning(
            'WARNING : Tchannel.applywavC is going to be replaced by Tchannel.applywav'))
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

    def baseband(self,**kwargs):
        """ Channel transfer function in baseband

        Parameters
        ----------

        fcGHz : center frequency
        WMHz  : bandwidth in MHz
        Nf    : Number of frequency points

        """

        defaults = {'fcGHz':4.5,
                    'WMHz':20,
                    'Nf':100}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        fcGHz = kwargs['fcGHz']
        WMHz  = kwargs['WMHz']
        Nf = kwargs['Nf']

        # self.y : Nray x Nr x Nt x Nf
        # self.taud : (,Nray)

        # complex amplitude in baseband
        # Nray x Nr x Nt x Nf1
        abb  = self.y*np.exp(-2 * 1j * np.pi *self.taud[:,None,None,None] * fcGHz )
        fMHz = np.linspace(-WMHz/2.,WMHz/2,Nf)
        E = np.exp(-2*1j*np.pi*fMHz[None,None,None,:]*1e-3*self.taud[:,None,None,None])
        y = np.sum(abb*E,axis=0)
        #y = np.sum(E,axis=0)
        H = bs.FUsignal(x=fMHz,y=y)

        return(H)

    def chantap(self,**kwargs):
        """ channel tap

        Parameters
        ----------

        fcGHz : center frequency
        WGHz  : bandwidth
        Ntap  : int

        """

        defaults = {'fcGHz':4.5,
                    'WGHz':1,
                    'Ntap':100}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        fcGHz=kwargs['fcGHz']
        WGHz=kwargs['WGHz']
        Ntap=kwargs['Ntap']
        # yb : tau x f x 1
        yb = self.y[:,:,None]*np.exp(-2 * 1j * np.pi *self.taud[:,None,None] * fcGHz )
        # l : 1 x 1 x tap
        l  = np.arange(Ntap)[None,None,:]
        # l : tau x 1 x 1
        tau = self.tau0[:,None,None]
        # S : tau x f x tap
        S   = np.sinc(l-tau*WGHz)
        # htap : f x tap
        htap  = np.sum(yb*S,axis=0)
        htapi = np.sum(htap,axis=0)

        return htapi

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

        print(DeprecationWarning(
            'WARNING : Tchannel.applywavB is going to be replaced by Tchannel.applywav'))

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
        print(DeprecationWarning(
            'WARNING : Tchannel.applywavA is going to be replaced by Tchannel.applywav'))
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
            inverse theta and phi representation
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
                    'phi': (-180, 180),
                    'reverse' : True,
                    'cmap': plt.cm.hot_r,
                    'vmin': [],
                    'vmax': [],
                    'mode': 'center',
                    's': 30,
                    'fontsize':12,
                    'edgecolors':'none',
                    'b3d':False,
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
#        b3d = kwargs.pop('b3d')
        the = (0,180)
        fontsize = kwargs.pop('fontsize')
        polar = kwargs.pop('polar')
        fig = kwargs.pop('fig')
        ax = kwargs.pop('ax')
        colorbar = kwargs.pop('colorbar')
        reverse = kwargs.pop('reverse')
        mode = kwargs.pop('mode')
        title =kwargs.pop('title')
        xa = kwargs.pop('xa')
        xb = kwargs.pop('xb')
        b3d = kwargs.pop('b3d')
        vmin = kwargs.pop('vmin')
        vmax = kwargs.pop('vmax')

        if fig == []:
            fig = plt.figure()


        Etot = self.energy(mode=mode) + 1e-15
        EtotdB = 10*np.log10(Etot)
        if vmax == []:
            vmax = EtotdB.max()
        if vmin == []:
            vmin = EtotdB.min()
        EtotdB = np.minimum(EtotdB,vmax)
        EtotdB = np.maximum(EtotdB,vmin)


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
        col = ((EtotdB - vmin)/(vmax-vmin)).squeeze()
        kwargs['c'] = col
        kwargs['s'] = 200*col
        kwargs['vmin'] = 0.
        kwargs['vmax'] = 1.

        if len(col) != len(di):
            print("len(col):", len(col))
            print("len(di):", len(di))
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



        Etot = self.energy(mode=mode) + 1e-15


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
            print("len(col):", len(col))
            print("len(di):", len(dir))
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

        Notes
        -----

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
                    'b3d':False,
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

        tau  = self.tau[:,None,None,None]
        fGHz = self.x[None,None,None,:]
        E = np.exp(-2*1j*tau*fGHz)
        F = self.y*E
        return np.sum(F,axis=0)
        #f = bs.FUsignal(x=self.x,y=np.sum(F,axis=0))
        #return(f)


    def energy(self,mode='mean',sumray=False):
        """ calculates channel energy including antennas spatial filtering

        Parameters
        ----------

        mode : string
            center | mean | integ  (different manner to get the value)
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

        if self.isFriis:
            Etot = bs.FUsignal.energy(self,axis=1,mode=mode,Friis=False)
        else:
            Etot = bs.FUsignal.energy(self,axis=1,mode=mode,Friis=True)

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
                print(" resampling w")
                x_new = arange(W.x[0], W.x[-1] + dxh, dxh)[0:-1]
                Wk = W.resample(x_new)
                dx = dxh
            else:
                # reinterpolate h
                print(" resampling h")
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

        PrdB: float
            RSSI value in dB
            $$10\log_10 |a_k|^2$$

        PrpdB : float
            RSSI in a tap over a frequency band ufreq
            $$10\log_10 |a_k e^{-2j \pi f \tau|^2$$


        Notes
        -----

            This function will be deprecated by energy function

        """
        # Amplitude
        Ak    = self.y[:, ufreq]
        # Power
        Pr    = np.sum(Ak*np.conj(Ak))
        # Complex amplitude
        akp   = Ak*np.exp(-2*1j*np.pi*self.x[ufreq]*self.taud)
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
        v = np.where(cumE[0,:]<threshold)[0]
        self.taud = self.taud[v]
        self.taue = self.taue[v]
        #self.tau = self.tau[v]
        self.doa = self.doa[v,:]
        self.dod = self.dod[v,:]
        self.y = self.y[v,...]

    def sort(self,typ='tau'):
        """ sort FUD signal

        Parameters
        ----------

        typ  : string
            which parameter to sort '
                'tau' : (default)
                'energy'

        """

        if typ == 'tau':
            u = np.argsort(self.taud+self.taue)

        if typ == 'energy':
            E = self.eprfl()
            u = np.argsort(E,axis=0)[::-1]
            u = u[:,0,0]

        self.taud = self.taud[u]
        self.taue = self.taue[u]
        self.doa = self.doa[u]
        self.dod = self.dod[u]
        self.y = self.y[u,...]
        return(u)

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
        Et_htap = np.sqrt(np.sum(htap*np.conj(htap),axis=i-1))/Nm
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
        ua = (((1/(Va+0.1))*np.random.rand(Ns)+ua0)%1)[:,None,None]
        va = (((1/(Va+0.1))*np.random.rand(Ns)+va0)%1)[:,None,None]
        ub = (((1/(Vb+0.1))*np.random.rand(Ns)+ub0)%1)[:,None,None]
        vb = (((1/(Vb+0.1))*np.random.rand(Ns)+vb0)%1)[:,None,None]

        # uniform sampling over the sphere
        tha = np.arccos(2*va-1)
        pha = 2*np.pi*ua
        thb = np.arccos(2*vb-1)
        phb = 2*np.pi*ub

        vax = np.cos(tha)*np.cos(pha)
        vay = np.cos(tha)*np.sin(pha)
        vaz = np.sin(tha)*np.cos(pha*0)

        vaxy = np.concatenate([vax[None,None,None,...],vay[None,None,None,...]])
        va = np.concatenate([vaxy,vaz[None,None,None,...]])

        vbx = np.cos(thb)*np.cos(phb)
        vby = np.cos(thb)*np.sin(phb)
        vbz = np.sin(thb)*np.cos(phb*0)

        vbxy = np.concatenate([vbx[None,None,None,...],vby[None,None,None,...]])

        # 3 x r x f x s x m x tap
        vb = np.concatenate([vbxy,vbz[None,None,None,...]])

        # beta : r x f x s x m x tap
        betaa = np.sum(ska[:,:,None,None,None,None]*va,axis=0)
        betab = np.sum(skb[:,:,None,None,None,None]*vb,axis=0)


        # m discrete time axis
        # r x f x s x m x tap
        m = np.linspace(0,mmax,Nm)[None,None,None,:,None]
        # r x f x s x m x tap
        l  = np.arange(Ntap)[None,None,None,None,:]
        # l : r x f x s x m x tap
        tau = self.taud[:,None,None,None,None]+ \
              self.taue[:,None,None,None,None]

        ba  = betaa*Va*m/(0.3*WMHz*1e6)
        bb  = betab*Vb*m/(0.3*WMHz*1e6)
        tau2 = tau + ba + bb
        # S : r x f x s x m x tap (form 2.34 [D. Tse])
        S   = np.sinc(l-tau2*WMHz/1000.)
        # sum over r :  f x s  x m x tap
        htap = np.sum(S*self.y[...,None,None,None]*np.exp(-2*1j*np.pi*fcGHz*tau2),axis=0)

        # f x s  x m x tap
        htap  = htap.reshape(Nf,Ns,Nm,Ntap)
        Et_htap = np.sqrt(np.sum(htap*np.conj(htap),axis=2))/Nm
        Er_htap = np.sum(htap,axis=1)/Ns
        corrtap = correlate(Er_htap[0,:,0],np.conj(Er_htap[0,:,0]))
        return(htap,Et_htap,Er_htap,corrtap)

    # def minphas(self):
    #     """ construct a minimal phase FUsignal

    #     - Evaluate slope of the phase
    #     - deduce delay
    #     - update delay of FUDSignal
    #     - Compensation of phase slope to obtain minimal phase

    #     This methods updates the excess delay taue member.

    #     The samplinf frequency step should be

    #     # Examples
    #     # --------

    #     # .. plot::
    #     #     :include-source:

    #     #     >>> from pylayers.signal.bsignal import *
    #     #     >>> import numpy as np
    #     #     >>> fGHz = np.arange(2,11,0.1)
    #     #     >>> tau1 = np.array([1,2,3])[:,None]
    #     #     >>> y = np.exp(-2*1j*np.pi*fGHz[None,:]*tau1)/fGHz[None,:]
    #     #     >>> H = Tchannel(x=fGHz,y=y,tau=np.array([15,17,18]))
    #     #     >>> f,a = H.plot(typ=['ru'],xlabels=['Frequency GHz'])
    #     #     >>> t1 = plt.suptitle('Before minimal phase compensation')
    #     #     >>> H.minphas()
    #     #     >>> H.taue
    #     #     array([ 1.,  2.,  3.])
    #     #     >>> f,a = H.plot(typ=['ru'],xlabels=['Frequency GHz'])
    #     #     >>> t2 = plt.suptitle('After minimal phase compensation')

    #     """

    #     f = self.x
    #     phase = np.unwrap(np.angle(self.y))
    #     dphi = phase[:, -1] - phase[:, 0]
    #     df = self.x[-1] - self.x[0]
    #     slope = dphi / df
    #     #if slope >0:
    #     #   print 'm  inphas Warning : non causal FUSignal'
    #     #phi0      = +1j*slope*(f[-1]+f[0]/2)
    #     F, S = np.meshgrid(f, slope)
    #     #E   = exp(-1j*slope*f+phi0)
    #     E = np.exp(-1j * S * F)
    #     self.y = self.y * E
    #     self.taue = -slope / (2 * np.pi)
    #     # update total delay
    #     #self.tau = self.tau+self.taue

    def sumf(self):
        tau  = self.taud + self.taue
        fGHz = self.x
        E    = np.exp(-2*1j*np.pi*fGHz[None,None,None,:]*tau[:,None,None,None])
        Y    = self.y*E
        sy   = np.sum(Y,axis=0)
        U    = bs.FUsignal(x=fGHz,y=np.squeeze(sy))
        return U


    def ifft(self):
        """ inverse Fourier Transform

        Examples
        --------

        >>> from pylayers.simul.link import *
        >>> L = DLink(verbose=False)
        >>> aktk = L.eval(force=True)
        >>> L.H.cut()
        >>> #T1 = L.H.totime()
        >>> #f,a = T1.plot(typ='v')
        >>> #L.H.minphas()
        >>> #T2 = L.H.totime()
        >>> #f,a = T2.plot(typ='v')


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

        >>> #from pylayers.simul.link import *
        >>> #L = DLink(verbose=False)
        >>> #aktk = L.eval()
        >>> #L.H.cut()
        >>> #T1 = L.H.totime()
        >>> #f,a = T1.plot(typ='v')
        >>> #L.H.minphas()
        >>> #T2 = L.H.totime()
        >>> #f,a = T2.plot(typ='v')

        See Also
        --------

        FUsignal.ift


        """
        Nray = len(self.taud)
        s = self.ift(Nz, ffts)
        sy_shifted = fft.fftshift(s.y,axes=-1)
        h = TUDchannel(s.x, sy_shifted, self.taud,self.taue)
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

        t0 = self.s.x[0]
        te = self.s.x[-1]

        shy = self.s.y.shape
        dx = self.s.x[1]-self.s.x[0]
        # Delta Tau + Npoints
        N  = np.ceil(dtau/dx)+shy[-1]

        # convert tau in an integer offset
        # taumin ray is not shifted
        itau = np.floor((tau-taumin)/dx).astype(int)

        U = np.ones((shy[0],shy[-1]),dtype=int)
        CU = np.cumsum(U,axis=1)-1 #-1 to start @ value 0

        rir  = np.zeros((shy[0],N))
        col1 = np.repeat(np.arange(shy[0],dtype=int),shy[-1])
        col2 = (CU+itau[:,None]).ravel()
        index = np.vstack((col1,col2)).T

        rir[index[:,0],index[:,1]] = self.s.y.ravel()
        t = np.linspace(t0+taumin,te+taumax,N)
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
            >>> #s = Tchannel(x=fGHz,y=alpha,tau=taud)
            >>> #s.plot3d()

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
        V = np.sum(ES, axis=0)
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
            T : Temperature (Kelvin)
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
        """ calibrate data 

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


    def scatterers(self,pMS,pBS,mode='sbounce'):
        """

        Parameters
        ----------

        mode: Boolean
              - image | sbounce (single bounce)

        pMS : np.array (,3)

        pBS : np.array (,3)
        ang_offset : float  (degrees)

        Returns
        -------

        xs: estimated x scatterer coordinate
        ys: estimated y scatterer coordinate

        """
        def fun(p_s,pBS,PMs,tau,phi):
            """ function to be minimized

            Parameters
            ----------

            p_s: np.array (,2)
                 2D scatterer coordinates containing the xs and ys coordinates

            tau   : float (in ns)
                    estimated delay
            phi   : float (in deg.)
                    estimated angle

            pMS : np.array (,2)
                  containing the x and y MS coordinates
            pBS : np.array (,2)
                  containing the x and y and z BS coordinates

            """

            d0  = np.sqrt((pBS[0] - p_s[0])**2 +
                          (pBS[1] - p_s[1])**2 +
                          (pBS[2] - p_s[2])**2 )  # distance between the BS and the estimated scatterer
            d1  = np.sqrt((pMS[0] - p_s[0])**2 +
                          (pMS[1] - p_s[1])**2 +
                          (pMS[2] - p_s[2])**2) # distance between the estimated scatterer and the MS
            xs = p_s[0]
            ys = p_s[1]
            zs = p_s[2]

            # Equations to be minimized

            r_fres = 0.3*tau
            eq1 = r_fres - (d0 + d1)
            eq2 = xs - pBS[0] - d0 * np.cos(phi)
            eq3 = ys - pBS[1] - d0 * np.sin(phi)

            return np.abs(eq1) + np.abs(eq2) + np.abs(eq3)

        Nscat = np.shape(self.y)[0]
        phi = self.doa[:,1]
        tau = self.taud

        if mode == 'image': # image principal mode
            xs = pBS[0] + tau*0.3*np.cos(phi)
            ys = pBS[1] + tau*0.3*np.sin(phi)

        if mode=='sbounce': # single bounce mode
            xs = np.array([])
            ys = np.array([])
            zs = np.array([])

            for k in range(Nscat): # loop over the all detected MPCs
                #zguess = (pMS[0:2] + pBS[0:2])/2. # Initial guess
                zguess = (pMS + pBS)/2. # Initial guess

                #z  = fmin(fun,zguess,(pBS[0:2],pMS[0:2],tau[k],phi[k]),disp=False) # minimizing Eq.
                z  = fmin(fun,zguess,(pBS,pMS,tau[k],phi[k]),disp=False) # minimizing Eq.
                xs = np.append(xs,z[0])
                ys = np.append(ys,z[1])
                zs = np.append(zs,z[2])

        return(xs,ys,zs)


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
        # by default C is expressed between the global frames
        self.islocal = False
        # by default antenna rotation matrices are identity 
        self.Ta = np.eye(3)
        self.Tb = np.eye(3)
        self.fGHz = np.array([2.4])
        # a single ray
        self.nray = 1
        self.Ctt = bs.FUsignal(x=self.fGHz,y=np.array([[1]]))
        self.Ctp = bs.FUsignal(x=self.fGHz,y=np.array([[0]]))
        self.Cpt = bs.FUsignal(x=self.fGHz,y=np.array([[0]]))
        self.Cpp = bs.FUsignal(x=self.fGHz,y=np.array([[1]]))
        #
        self.tang = np.array([[np.pi/2,np.pi/2]])
        self.rang = np.array([[np.pi/2,3*np.pi/2]])
        #
        self.tangl = np.array([[np.pi/2,np.pi/2]])
        self.rangl = np.array([[np.pi/2,3*np.pi/2]])

    def __repr__(self):
        s = 'Ctilde : Ray Propagation Channel Tensor (2x2xrxf)'+'\n---------\n'
        if self.islocal:
            s = s + 'between antennas local frames\n'
        else:
            s = s + 'between termination global frames\n'

        s = s + 'Nray : ' + str(self.nray)+'\n'
        if self.Cpp.x[0]!=self.Cpp.x[-1]:
            s = s + 'Nfreq : ' + str(len(self.Cpp.x))+'\n'
            s = s + 'fmin(GHz) : ' + str(self.Cpp.x[0])+'\n'
            s = s + 'fmax(GHz): ' + str(self.Cpp.x[-1])+'\n'
        else:
            s = s + 'fGHz : ' + str(self.Cpp.x[0])+'\n'
        s = s + '---1st ray 1st freq ---\n'
        s = s + 'global angles (th,ph) degrees : \n'
        s = s + str(np.round(self.tang[0,:]*1800/np.pi)/10.)+'\n'
        s = s + str(np.round(self.rang[0,:]*1800/np.pi)/10.)+'\n'
        s = s + 'local angles (th,ph) degrees : \n'
        s = s + str(np.round(self.tangl[0,:]*1800/np.pi)/10.)+'\n'
        s = s + str(np.round(self.rangl[0,:]*1800/np.pi)/10.)+'\n'
        s = s + '    | '+ str(self.Ctt.y[0,0])+'   '+str(self.Ctp.y[0,0])+' |\n'
        s = s + '    | '+ str(self.Cpt.y[0,0])+'   '+str(self.Cpp.y[0,0])+' |\n'
        return(s)

    def inforay(self,iray,ifreq=0):
        """ provide information about a specific ray

        """
        dray   = self.tauk[iray]*0.3
        draydB = 20*np.log10(1./dray)

        Ctt = self.Ctt.y[iray,ifreq]
        Ctp = self.Ctp.y[iray,ifreq]
        Cpt = self.Cpt.y[iray,ifreq]
        Cpp = self.Cpp.y[iray,ifreq]


        Cttc = Ctt*dray
        Ctpc = Ctp*dray
        Cppc = Cpp*dray
        Cptc = Cpt*dray
        if self.islocal:
            print("between local frames")
            print("--------------------")
        else:
            print("between global frames")
            print("--------------------")
        print('distance losses {:.2f}'.format(draydB))
        if (np.abs(Cttc)!=0):
            CttdB = 20*np.log10(np.abs(Ctt))
            CttcdB = 20*np.log10(np.abs(Cttc))
        else:
            CttdB = -np.inf
            CttcdB = -np.inf
        if (np.abs(Cppc)!=0):
            CppdB = 20*np.log10(np.abs(Cpp))
            CppcdB = 20*np.log10(np.abs(Cppc))
        else:
            CppdB = -np.inf
            CppcdB = -np.inf
        if (np.abs(Ctpc)!=0):
            CtpdB = 20*np.log10(np.abs(Ctp))
            CtpcdB =20*np.log10(np.abs(Ctpc))
        else:
            CtpdB = -np.inf
            CtpcdB = -np.inf
        if (np.abs(Cptc)!=0):
            CptdB = 20*np.log10(np.abs(Cpt))
            CptcdB = 20*np.log10(np.abs(Cptc))
        else:
            CptdB = -np.inf
            CptcdB = -np.inf
        print('Without distance losses (Interactions only)')
        print("-----------------------------------------------")
        print('co-pol (tt,pp)    dB : {:.2f} {:.2f}'.format(CttcdB,CppcdB))
        print('cross-pol (tt,pp) dB : {:.2f} {:.2f}'.format(CtpcdB,CptcdB))

        print('With distance losses (Interactions + distance)')
        print("-----------------------------------------------")
        print('co-pol (tt,pp)    dB : {:.2f} {:.2f}'.format(CttdB,CppdB))
        print('cross-pol (tp,pt) dB : {:.2f} {:.2f}'.format(CtpdB,CptdB))


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
        # new call to locbas
        if self.islocal:
            self.locbas()

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f=h5py.File(filename,'w')
            f.create_dataset('Ta',shape=np.shape(self.Ta),data=self.Ta)
            f.create_dataset('Tb',shape=np.shape(self.Tb),data=self.Tb)
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

            self.Ta = f['Ta'][:]
            self.Tb = f['Tb'][:]

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
        # back to global frame
        if self.islocal:
            self.locbas()


        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:

            fh5=h5py.File(filename,'a')
            if not grpname in fh5['Ct'].keys():
                fh5['Ct'].create_group(grpname)
            else :
                print('Warning : Ct/'+grpname +'already exists in '+filenameh5)
            f=fh5['Ct/'+grpname]

            # save channel in global basis
            f.create_dataset('Ta',shape=np.shape(self.Ta),data=self.Ta)
            f.create_dataset('Tb',shape=np.shape(self.Tb),data=self.Tb)
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

    def los(self,**kwargs):
        """ Line of site channel

        Parameters
        ----------

        d(m)
        fGHz (,Nf)
        tang (1x2)
        rang (1x2)

        """
        defaults = {'pa':np.r_[197,189.8,1.65]
                   ,'pb': np.r_[220,185,6]
                   ,'fGHz':np.r_[32.6] 
                   ,'Ta':np.eye(3)
                   ,'Tb':np.array([[0.28378894, -0.8972627,  -0.33820628],
                        [-0.57674955, -0.44149706,  0.68734293],
                         [-0.76604425,  0.,         -0.64278784]])
                   }

        for k in defaults: 
            if k not in kwargs:
                kwargs[k]=defaults[k]
        self.pa = kwargs['pa']
        self.pb = kwargs['pb'] 
        self.fGHz = kwargs['fGHz']
        self.Ta = kwargs['Ta']
        self.Tb = kwargs['Tb']
        self.nray = 1

        si = self.pb-self.pa

        d = np.r_[np.sqrt(np.sum(si*si))]
        si = si/d
        self.tauk = d/0.3
        #
        # ka = - kb for LOS 
        #
        tha =  np.arccos(si[2])
        pha =  np.arctan2(si[1],si[0])
        thb =  np.arccos(-si[2])
        phb =  np.arctan2(-si[1],-si[0])

        self.tang = np.array([tha,pha]).reshape((1,2))
        self.rang = np.array([thb,phb]).reshape((1,2))

        U = np.ones(len(self.fGHz),dtype=complex)/d[0]
        Z = np.zeros(len(self.fGHz),dtype=complex)
         
        self.Ctt = bs.FUsignal(self.fGHz, U)
        self.Ctp = bs.FUsignal(self.fGHz, Z)
        self.Cpt = bs.FUsignal(self.fGHz, Z)
        self.Cpp = bs.FUsignal(self.fGHz, U)
        self.locbas()

    def _loadh5(self,filenameh5,grpname,**kwargs):
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

            self.Ta = f['Ta'][:]
            self.Tb = f['Tb'][:]

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
                    'fontsize':22,
                    'edgecolors':'none',
                    'b3d':False,
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
        b3d = kwargs.pop('b3d')
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
            print("len(col):", len(col))
            print("len(di):", len(dir))
        if b3d:
            ax = fig.add_subplot(111,projection='3d')
            ax.scatter(1.05*array(xa),1.05*array(ya),1.05*array(za),'b')
            ax.scatter(1.05*array(xb),1.05*array(yb),1.05*array(zb),'r')
        else:
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

        Notes
        -----

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
                    'b3d':False,
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


    def locbas(self,**kwargs):
        """ global reference frame to local reference frame

        If Tt and Tr are [] the global channel is  retrieved

        Parameters
        ----------

        Ta  : rotation matrix 3x3  side a 
            default []
        Tb  : rotation matrix 3x3  side b
            default []

        Returns
        -------
        
        This method affects the boolean islocal 
        This method update the ray propagation channel in either local or global frame
        self.Ta and self.Tb are updated with input parameters Ta an Tb 

        C : ray propagation channel (2x2xrxf) complex 
            either local or global depends on self.islocal boolean value

        Examples
        --------

        >>> C = Ctilde()
        >>> Ta = MEulerAngle(np.pi/2,np.pi/2,np.pi/2.)
        >>> Tb = MEulerAngle(np.pi/3,np.pi/3,np.pi/3.)
        >>> C.locbas(Ta=Ta,Tb=Tb)

        """
        
        # get Ctilde frequency axes

        fGHz = self.fGHz
#       if rotation matrices are passed in argument 
#       back to global if local 
        if ('Ta' in kwargs) & ('Tb' in kwargs):
            if self.islocal:
                self.locbas()
                self.islocal=False
            self.Tb = kwargs['Tb']
            self.Ta = kwargs['Ta']
        # angular axes
        #
        # tang : r x 2
        # rang : r x 2
        #
        # Ra : 2 x 2 x r
        # Rb : 2 x 2 x r
        #
        # tangl : r x 2
        # rangl : r x 2
        #
        tangl,Ra = geu.BTB(self.tang, self.Ta)
        rangl,Rb = geu.BTB(self.rang, self.Tb)

        if self.islocal:
            Ra = Ra.transpose((1,0,2))
            self.islocal=False
        else:
            Rb = Rb.transpose((1,0,2))
            self.islocal=True

        #
        # update direction of departure and arrival
        #

        self.tangl = tangl
        self.rangl = rangl

        #uf = np.ones(self.nfreq)

        #
        # r0 : r x 1(f)
        #

        #r0 = rb00
        r0 = Rb[0,0,:][:, None]
        #r1 = rb01
        r1 = Rb[0,1,:][:, None]

        t00 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t01 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = rb10 
        r0 = Rb[1, 0,:][:, None]
        #r1 = rb11 
        r1 = Rb[1, 1,:][:, None]

        t10 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t11 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = ra00 
        r0 = Ra[0, 0, :][:, None]
        #r1 = ra10 
        r1 = Ra[1, 0, :][:, None]

        Cttl = t00 * r0 + t01 * r1
        Cptl = t10 * r0 + t11 * r1

        #r0 = ra01
        r0 = Ra[0, 1, :][:, None]
        #r1 = ra11
        r1 = Ra[1, 1, :][:, None]

        Ctpl = t00 * r0 + t01 * r1
        Cppl = t10 * r0 + t11 * r1

        self.Ctt = bs.FUsignal(fGHz, Cttl)
        self.Ctp = bs.FUsignal(fGHz, Ctpl)
        self.Cpt = bs.FUsignal(fGHz, Cptl)
        self.Cpp = bs.FUsignal(fGHz, Cppl)

        #return self

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
            self.Ta = Tt
            self.Tb = Tr
        else:
            if (hasattr(self,'Ta')) & (hasattr(self, 'Tb')):
                self.Ta = self.Ta.transpose()
                self.Tb = self.Tb.transpose()
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

        tangl , Ra = geu.BTB(self.tang, self.Ta)
        rangl , Rb = geu.BTB(self.rang, self.Tb)
        Rb = Rb.transpose((1,0,2))
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
        r0 = Rr[0,0,:][:,None]
        #r1 = np.outer(Rr[0, 1,:], uf)
        r1 = Rr[0,1,:][:,None]

        t00 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t01 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rr[1, 0,:], uf)
        r0 = Rr[1, 0,:][:,None]
        #r1 = np.outer(Rr[1, 1,:], uf)
        r1 = Rr[1, 1,:][:,None]

        t10 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t11 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rt[0, 0,:], uf)
        r0 = Rt[0,0,:][:,None]
        #r1 = np.outer(Rt[1, 0,:], uf)
        r1 = Rt[1,0,:][:,None]

        Cttl = t00 * r0 + t01 * r1
        Cptl = t10 * r0 + t11 * r1

        #r0 = np.outer(Rt[0, 1,:], uf)
        r0 = Rt[0,1,:][:,None]
        #r1 = np.outer(Rt[1, 1,:], uf)
        r1 = Rt[1,1,:][:,None]

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
            print("Channel is reciprocal")
        else: 
            print("WARNING Reciprocity issue WARNING")
            print(len(issue),'/',self.nray, 'rays are not reciprocal,')
            print("rays number with an issue :",issue)

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

    def cut(self,threshold_dB=50):
        """ cut rays from a energy threshold

        Parameters
        ----------

        threshold : float
            default 0.99

        """
        Ett, Epp, Etp, Ept = self.energy()
        Etot = Ett+Epp+Etp+Ept
        u = np.argsort(Etot)[::-1]
        #cumE = np.cumsum(Etot[u])/sum(Etot)
        profdB = 10*np.log10(Etot[u]/np.max(Etot))
        #v1 = np.where(cumE<threshold)[0]
        v = np.where(profdB>-threshold_dB)[0]
        w = u[v]
        self.selected = w
        self.Eselected = Etot[w]
        self.tauk = self.tauk[w]
        self.tang = self.tang[w,:]
        self.rang = self.rang[w,:]
        self.Ctt.y = self.Ctt.y[w,:]
        self.Cpp.y = self.Cpp.y[w,:]
        self.Ctp.y = self.Ctp.y[w,:]
        self.Cpt.y = self.Cpt.y[w,:]

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

    def prop2tran(self,a=[],b=[],Friis=True,debug=False):
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

        a.eval(th = self.tangl[:, 0], ph = self.tangl[:, 1])
        Fat = bs.FUsignal(a.fGHz, a.Ft)
        Fap = bs.FUsignal(a.fGHz, a.Fp)
        #b.eval(th=self.rangl[:, 0], ph=self.rangl[:, 1], grid=False)
        b.eval(th = self.rangl[:, 0], ph = self.rangl[:, 1])

        Fbt = bs.FUsignal(b.fGHz, b.Ft)
        Fbp = bs.FUsignal(b.fGHz, b.Fp)

        #
        #  C  :  2 x 2 x r x f
        #
        #  Ctt : r x f     (complex FUsignal)
        #  Cpp : r x f     (complex FUsignal)
        #  Ctp : r x f     (complex FUsignal)
        #  Cpt : r x f     (complex FUsignal)
        #
        #  a.Ft = r x (Na) x f  (complex ndarray)
        #  a.Fp = r x (Na) x f  (complex ndarray)
        #  b.Ft = r x (Nb) x f  (complex ndarray)
        #  b.Fp = r x (Nb) x f  (complex ndarray)
        #
        #  (r x f ) (r x Nt x f )
        #
        # This exploit * overloading in FUsignal 

        t1 = self.Ctt * Fat + self.Ctp * Fap
        t2 = self.Cpt * Fat + self.Cpp * Fap

        # depending on SISO or MIMO case
        # the shape of the received fields T1 and T2 
        # 
        # In MIMO case
        #   a.Ft.y.shape == (r x Na x f) 
        #   a.Fp.y.shape == (r x Na x f)
        # In SISO case 
        #   a.Ft.y.shape == (r x f) 
        #   a.Fp.y.shape == (r x f)
        #
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


        return H

if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
