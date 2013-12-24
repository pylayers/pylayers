# -*- coding:Utf-8 -*-
import doctest
import pdb
import numpy as np
import scipy as sp
import pylab as plt
import struct as stru
import pylayers.util.pyutil as pyu
import pylayers.signal.bsignal as bs
import pylayers.util.geomutil as geu
from pylayers.antprop.raysc import GrRay3D
from pylayers.util.project import *


class Ctilde(object):
    """ container for the 4 components of the polarimetric ray channel

    Attributes
    ----------

    Ctt : FUsignal
    Ctp : FUsignal
    Cpt : FUsignal
    Cpp : FUsignal
    fGHz : np.array
        frequency array
    nfreq : int
        number of frequency point
    nray  : int
        number of rays

    """
    def __init__(self):
        """
            transpose == False   (r,f)
            transpose == True    (f,r)

            A Ctilde object is the output of eval method
            of a Rays object.

        """
        self.fail = False

    def __repr__(self):
        s = 'Ctilde'+'\n---------\n'
        if hasattr(self,'Cpp'):
            s = s +str(np.shape(self.Cpp.y))+'\n'
        if hasattr(self,'nray'):
            s = s+ 'Nray : '+ str(self.nray)+'\n'
            s = s+ 'fmin(GHz) : '+ str(self.Cpp.x[0])+'\n'
            s = s+ 'fmax(GHz): '+ str(self.Cpp.x[-1])+'\n'
            s = s+ 'Nfreq : '+ str(self.nfreq)+'\n'
        s = s+ '\n methods :'+'\n---------\n'
        s = s+ 'prop2tran(a=theta,b=phi)\n'
        s = s+ 'energy()\n'
        s = s+ 'doadod(cmap=plt.cm.hot_r,s=30,fontsize=12,phi=(0,360))\n'
        s = s+ 'mobility(v,dt)\n'
        s = s+ 'show(mode=linear)\n'
        s = s+ 'sort()\n'
            
        return(s)


    def choose(self):
        """ Choose a field file in tud directory
        """
        import tkFileDialog as FD
        filefield = FD.askopenfilename(filetypes=[("Files field  ", "*.field"),
                                                  ("All", "*")],
                                       title="Please choose a .field file",
                                       initialdir=tuddir)
        self.load(filefield, transpose=False)

    def load(self, filefield, transpose=False):
        """ load a Ctilde from a .field file

        Load the three files .tauk .tang .rang which contain  respectively
        delay , angle of departure , angle of arrival.

        Parameters
        ----------

        filefield  : string
        transpose  : boolean
            default False

        Examples
        --------

        >>> from pylayers.antprop.channel import *
        >>> from pylayers.simul.simulem import *
        >>> S = Simul()
        >>> S.load('default.ini')
        >>> C = Ctilde()
        >>> C.load(pyu.getlong(S.dtud[1][1],'output'))

        """
        filetauk = filefield.replace('.field', '.tauk')
        filetang = filefield.replace('.field', '.tang')
        filerang = filefield.replace('.field', '.rang')
        try:
            fo = open(filefield, "rb")
        except:
            raise NameError( "file "+ filefield+ " is unreachable")

        # decode filename (*.field file obtained from evalfield simulation)
        nray = stru.unpack('i', fo.read(4))[0]
        print "nray : ", nray
        nfreq = stru.unpack('i', fo.read(4))[0]
        print "nfreq : ", nfreq
        if nfreq == 0:
            print " Error : incorrect number of frequency points in .field"
            self.fail = True
            return

        n = nray * nfreq
        buf = fo.read()
        fo.close()

        CMat = np.ndarray(shape=(n,8),buffer=buf)
        c11 = CMat[:, 0] + CMat[:, 1]*1j
        c12 = CMat[:, 2] + CMat[:, 3]*1j
        c21 = CMat[:, 4] + CMat[:, 5]*1j
        c22 = CMat[:, 6] + CMat[:, 7]*1j

        c11 = c11.reshape(nray, nfreq)
        c12 = c12.reshape(nray, nfreq)
        c21 = c21.reshape(nray, nfreq)
        c22 = c22.reshape(nray, nfreq)

        if transpose:
            c11 = c11.transpose()
            c12 = c12.transpose()
            c21 = c21.transpose()
            c22 = c22.transpose()

        #
        # Temporary freq --> read filefreq
        #
        fGHz = np.linspace(2, 11, nfreq)

        self.Ctt = bs.FUsignal(fGHz, c11)
        self.Ctp = bs.FUsignal(fGHz, c12)
        self.Cpt = bs.FUsignal(fGHz, c21)
        self.Cpp = bs.FUsignal(fGHz, c22)
        self.nfreq = nfreq
        self.nray = nray

        try:
            fo = open(filetauk, "rb")
        except:
            self.fail = True
            print "file ", filetauk, " is unreachable"
        # decode filetauk
        if not self.fail:
            nray_tauk = stru.unpack('i', fo.read(4))[0]
            print "nb rays in .tauk file: ", nray_tauk
            buf = fo.read()
            fo.close()
            nray = len(buf) / 8
            print "nb rays 2: ", nray
            self.tauk = np.ndarray(shape=nray, buffer=buf)
            if nray_tauk != nray:
                print nray_tauk - nray
        self.tauk = self.tauk

        # decode the angular files (.tang and .rang)
        try:
            fo = open(filetang, "rb")
        except:
            self.fail = True
            print "file ", filetang, " is unreachable"
        if not self.fail:
            nray_tang = stru.unpack('i', fo.read(4))[0]
            buf = fo.read()
            fo.close()
            # coorectif Bug evalfield
            tmp = np.ndarray(shape=(nray_tang, 2), buffer=buf)
            self.tang = tmp[0:nray, :]
        try:
            fo = open(filerang, "rb")
        except:
            self.fail = True
            print "file ", filerang, " is unreachable"

        if not self.fail:
            nray_rang = stru.unpack('i', fo.read(4))[0]
            buf = fo.read()
            fo.close()
            # corectif Bug evalfield
            tmp = np.ndarray(shape=(nray_rang, 2), buffer=buf)
            self.rang = tmp[0:nray, :]

    def mobility(self, v, dt):
        """ Modify channel for uniform mobility

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

        #cos_alph =dot(v/abs(v),sk_ch)

        cos_alph = (v * sk_ch) / abs(v)
        self.cos_alph = cos_alph
        rk_ch = rk_mod * cos_alph * abs(v) * dt
        sk_ch_ch = (rk + v * dt) / (rk_ch + cos_alph * abs(v) * dt)
        tauk_ch = (abs(rk_ch) * sk_ch_ch) / c

        return(tauk_ch)

    def doadod(self, cmap=plt.cm.hot_r, s=30,fontsize = 12,phi=(0,360)):
        """ doadod scatter plot

        Parameters
        -----------

        cmap : color map
        s    : float
            size (default 30)
        fontsize : integer
            default 12

        Summary
        --------

        scatter plot of the DoA-DoD channel structure
        the energy is colorcoded over all couples of DoA-DoD

        """
        dod = self.tang
        doa = self.rang
        # determine Energy in each channel 
        Ett,Epp,Etp,Ept = self.energy()
        Etot = Ett+Epp+Etp+Ept
        Emax = max(Etot)
        Etot = Etot / Emax + 1e-7
        Emax = max(10 * np.log10(Etot))
        Emin = min(10 * np.log10(Etot))
        #
        #
        #
        #col  = 1 - (10*log10(Etot)-Emin)/(Emax-Emin)
        al = 180. / np.pi
        col = 10 * np.log10(Etot)
        if len(col) != len(dod):
            print "len(col):", len(col)
            print "len(dod):", len(dod)
        plt.subplot(121)
        plt.scatter(dod[:, 0] * al, dod[:, 1] * al, s=s, c=col,
                    cmap=cmap, edgecolors='none')
        #scatter(dod[:,0]*al,dod[:,1]*al,s=s)
        plt.axis((0, 180, phi[0], phi[1]))
        #plt.xticks(fontsize=20)
        #plt.yticks(fontsize=20)
        #a = plt.colorbar()
        #for t in a.ax.get_yticklabels():
        #    t.set_fontsize(18)
        #a.set_label('dB')
        plt.xlabel("$\\theta_t(\degree)$", fontsize=fontsize)
        plt.ylabel('$\phi(\degree)$', fontsize=fontsize)
        #ylabel('$\phi_t(\degree)$',fontsize=18)
        plt.title('DoD',fontsize=fontsize+2)
        plt.subplot(122)
        plt.scatter(doa[:, 0] * al, doa[:, 1] * al, s=30, c=col,
                    cmap=plt.cm.hot_r, edgecolors='none')
        plt.axis((0, 180, phi[0], phi[1]))
        #plt.xticks(fontsize=20)
        #plt.yticks(fontsize=20)
        b = plt.colorbar()
        b.set_label('dB')
        #for t in b.ax.get_yticklabels():
        #    t.set_fontsize(20)
        plt.xlabel("$\\theta_r(\degree)$", fontsize=fontsize)
        plt.title('DoA', fontsize=fontsize+2)
        plt.ylabel("$\phi_r (\degree)$", fontsize=fontsize)
        plt.axis

    def show(self, display=False, mode='linear'):
        """ show the propagation channel 
        
        Parameters
        ----------
        
        display : True or False
        mode    : 'linear', 'dB'
        
        """

        f = abs(self.Ctt.x)
        u = np.argsort(self.tauk)
        tt = self.tauk[u]
        utt = abs(self.Ctt.y[u, :])
        utp = abs(self.Ctp.y[u, :])
        upt = abs(self.Ctp.y[u, :])
        upp = abs(self.Cpp.y[u, :])

        uttmax = utt.max()
        utpmax = utp.max()
        uptmax = upt.max()
        uppmax = upp.max()

        #vmax=max(uttmax,utpmax,uptmax,uppmax)
        if mode == 'linear':
            plt.figure()
            plt.subplot(221)
            plt.pcolor(f, tt, utt)
            plt.colorbar()
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.title('Ctt')
            plt.subplot(222)
            plt.pcolor(f, tt, utp)
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.colorbar()
            plt.title('Ctp')
            plt.subplot(223)
            plt.pcolor(f, tt, upt)
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.colorbar()
            plt.title('Cpt')
            plt.subplot(224)
            plt.pcolor(f, tt, upp)
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.colorbar()
            plt.title('Cpp')
        else:
            plt.figure()
            plt.subplot(221)
            plt.pcolor(f, tt, 20 * np.log10(utt + 1e-5))
            plt.colorbar()
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.title('Ctt')
            plt.subplot(222)
            plt.pcolor(f, tt, 20 * np.log10(utp + 1e-5))
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.colorbar()
            plt.title('Ctp')
            plt.subplot(223)
            plt.pcolor(f, tt, 20 * np.log10(upt + 1e-5))
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.colorbar()
            plt.title('Cpt')
            plt.subplot(224)
            plt.pcolor(f, tt, 20 * np.log10(upp + 1e-5))
            plt.xlabel('f (GHz)')
            plt.ylabel('delay (ns)')
            plt.colorbar()
            plt.title('Cpp')
        if display:
            plt.show()

    def energy(self):
        """ Calculates energy on each channel

        Returns
        -------

        ECtt  : Energy on co channel    tt 
        ECpp  : Energy on co channel    pp 
        ECtp  : Energy on co channel    tp 
        ECpt  : Energy on co channel    pt 

        
        See Also
        --------

        pylayers.signal.bsignal.FUsignal.energy

    

        """

        ECtt = self.Ctt.energy(1)
        ECtp = self.Ctp.energy(1)
        ECpt = self.Cpt.energy(1)
        ECpp = self.Cpp.energy(1)
        #Eco = ECtt + ECpp
        #Ecross = ECtp + ECpt

        return ECtt,ECpp,ECtp,ECpt

    def sort(self,typ='tauk'):
        """ sort Ctilde with respect to tauk

        Parameters
        ----------

        typ  : string 
            which parameter to sort '
                tauk : (default)
                att  : theta Tx
                atp  : phi Tx
                art  : theta Rx 
                arp  : phi Rx 
                energy 

        """
        if typ=='tauk':
            u = np.argsort(self.tauk)
        if typ=='att':
            u = np.argsort(self.tang[:,0])
        if typ=='atp':
            u = np.argsort(self.tang[:,1])
        if typ=='art':
            u = np.argsort(self.rang[:,0])
        if typ=='arp':
            u = np.argsort(self.rang[:,1])
        if typ=='energy':
            Ett,Epp,Etp,Ept=self.energy()
            Etot = Ett+Epp+Etp+Ept 
            u = np.argsort(Etot)

        self.tauk = self.tauk[u]
        self.tang = self.tang[u,:]
        self.rang = self.rang[u,:]

        self.Ctt.y = self.Ctt.y[u, :]
        self.Cpp.y = self.Cpp.y[u, :]
        self.Ctp.y = self.Ctp.y[u, :]
        self.Cpt.y = self.Cpt.y[u, :]

    def prop2tran(self,a='theta',b='theta'):
        """ transform propagation channel into transmission channel

        Parameters
        ----------

        a : string or antenna array
            polarization antenna a ( 'theta' | 'phi' | 'ant' ) 
        b : string or antenna array
            polarization antenna b ( 'theta' | 'phi' | 'ant' ) 

            0 : theta
            1 : phi

        Returns
        -------

        H : Tchannel(bs.FUDAsignal)


        """
        freq = self.fGHz
        nfreq = self.nfreq
        nray  = self.nray
        sh = np.shape(self.Ctt.y)

        if type(a)==str:
            Fat = np.zeros(nray*nfreq).reshape((nray,nfreq))
            Fap = np.zeros(nray*nfreq).reshape((nray,nfreq))
            if a=='theta':
                Fat = np.ones((nray,nfreq))
            if a=='phi':
                Fap = np.ones((nray,nfreq))
            Fat = bs.FUsignal(self.fGHz,Fat)
            Fap = bs.FUsignal(self.fGHz,Fap)
        else:
            Fat , Fap = a.Fsynth3(self.rang[:, 0],self.rang[:,1],pattern=False)
            Fat = Fat.transpose()
            Fap = Fap.transpose()
            Fat = bs.FUsignal(a.fa,Fat)
            Fap = bs.FUsignal(a.fa,Fap)

        if type(b)==str:
            Fbt = np.zeros(nray*nfreq).reshape((nray,nfreq))
            Fbp = np.zeros(nray*nfreq).reshape((nray,nfreq))
            if b=='theta':
                Fbt = np.ones((nray,nfreq))
            if b=='phi':
                Fbp = np.ones((nray,nfreq))
            Fbt = bs.FUsignal(self.fGHz,Fbt)
            Fbp = bs.FUsignal(self.fGHz,Fbp)
        else:
            Fbt , Fbp = b.Fsynth3(self.rang[:, 0],self.rang[:,1],pattern=False)

            Fbt = Fbt.transpose()
            Fbp = Fbp.transpose()
            Fbt = bs.FUsignal(b.fa,Fbt)
            Fbp = bs.FUsignal(b.fa,Fbp)
        t1 = self.Ctt * Fat + self.Cpt * Fap
        t2 = self.Ctp * Fat + self.Cpp * Fap
        alpha = t1 * Fbt + t2 * Fbp

        H = Tchannel(alpha.x,alpha.y,self.tauk,self.tang,self.rang)
        return(H)

    def vec2scal(self):
        """ calculate scalChannel from VectChannel and antenna

        Returns
        -------

        slach : ScalChannel

        Note 
        ----

        deprecated this is now replaced by the function prop2tran

        """
        fGHz = self.fGHz
        sh = np.shape(self.Ctt.y)
        Ftt = bs.FUsignal(fGHz, np.ones(sh))
        Ftp = bs.FUsignal(fGHz, np.zeros(sh))
        Frt = bs.FUsignal(fGHz, np.ones(sh))
        Frp = bs.FUsignal(fGHz, np.zeros(sh))
        scalch = ScalChannel(self, Ftt, Ftp, Frt, Frp)
        return(scalch)

    # Inclusion of realistic antenna behaviour
    # Specify transmitter antenna and receiver antenna

    def vec2scalA(self, At, Ar, alpha=1.0):
        """
        
        Parameters
        ----------

        At : transmitter antenna
        Ar : receiver antenna
        alpha : normalization factor
    
        Notes
        -----

        Calculate ScalChannel by combining the propagation channel VectChannel
        with realistic antennas transfer function


        """

        Ftt, Ftp = At.Fsynth3(self.rang[:, 0], self.rang[:, 1],pattern=False)

        Frt, Frp = Ar.Fsynth3(self.rang[:, 0], self.rang[:, 1],pattern=False)

        Ftt = Ftt.transpose()
        Ftp = Ftp.transpose()
        Frt = Frt.transpose()
        Frp = Frp.transpose()

        Ftt = bs.FUsignal(At.fa, Ftt * alpha)
        Ftp = bs.FUsignal(At.fa, Ftp * alpha)
        Frt = bs.FUsignal(Ar.fa, Frt * alpha)
        Frp = bs.FUsignal(Ar.fa, Frp * alpha)

        scalch = ScalChannel(self, Ftt, Ftp, Frt, Frp)
        return(scalch)

    def info(self):
        """ Info (Nf,Nray,shape(y))
        """

        print "Nfreq  :", self.nfreq
        print "Nray :", self.nray
        print "shape Ctt :", np.shape(self.Ctt.y)
        print "shape Ctp :", np.shape(self.Ctp.y)
        print "shape Cpt :", np.shape(self.Cpt.y)
        print "shape Cpp :", np.shape(self.Cpp.y)


def Cg2Cl(Cg, Tt, Tr):
    """ global reference frame to local reference frame

    Parameters
    ----------

    Cg  : Ctilde global
    Tt  : Tx rotation matrix
    Tr  : Rx rotation matrix

    Returns
    -------

    Cl : Ctilde local 

    Examples
    --------

    """
    import copy
    
    # don't loose the global channel
    Cl = copy.deepcopy(Cg)
    
    # get frequency axes    
    fGHz = Cl.fGHz
    
    # get angular axes

    Rt, tangl = BTB_tx(Cg.tang, Tt)
    Rr, rangl = BTB_rx(Cg.rang, Tr)

    VCl.tang = tangl
    VCl.rang = rangl

    uf = np.ones(VCg.nfreq)
    r0 = np.outer(Rr[0, 0, :], uf)
    r1 = np.outer(Rr[0, 1, :], uf)

    # print "shape r0 = ",np.shape(r0)
    # print "shape VCg.Ctt.y = ",np.shape(VCg.Ctt.y)
    # print "shape r1 = ",np.shape(r1)
    # print "shape VCg.Cpt.y = ",np.shape(VCg.Cpt.y)

    t00 = r0 * VCg.Ctt.y + r1 * VCg.Cpt.y
    t01 = r0 * VCg.Ctp.y + r1 * VCg.Cpp.y

    r0 = np.outer(Rr[1, 0, :], uf)
    r1 = np.pouter(Rr[1, 1, :], uf)
    t10 = r0 * VCg.Ctt.y + r1 * VCg.Cpt.y
    t11 = r0 * VCg.Ctp.y + r1 * VCg.Cpp.y

    r0 = np.outer(Rt[0, 0, :], uf)
    r1 = np.outer(Rt[1, 0, :], uf)

    Cttl = t00 * r0 + t01 * r1
    Cptl = t10 * r0 + t11 * r1

    r0 = np.outer(Rt[0, 1, :], uf)
    r1 = np.outer(Rt[1, 1, :], uf)
    Ctpl = t00 * r0 + t01 * r1
    Cppl = t10 * r0 + t11 * r1

    VCl.Ctt = bs.FUsignal(fGHz, Cttl)
    VCl.Ctp = bs.FUsignal(fGHz, Ctpl)
    VCl.Cpt = bs.FUsignal(fGHz, Cptl)
    VCl.Cpp = bs.FUsignal(fGHz, Cppl)

    return VCl


class Tchannel(bs.FUDAsignal):
    """ Handle the transmission channel 

    The transmission channel TChannel is obtained through combination of the propagation
    channel and the antenna transfer functions from both transmitter and
    receiver.

    Members
    -------

        ray transfer functions  (nray,nfreq)
    dod  :
        direction of depature (rad) [theta_t,phi_t]  nray x 2
    doa  :
        direction of arrival (rad)  [theta_r,phi_r]  nray x 2
    tauk :
        delay ray k in ns


    """
    def __init__(self,fGHz,alpha,tau,dod,doa):
        """

        Parameters
        ----------

        fGHz  :  1 x nfreq
        alpha :  channel amplitude (nray x nfreq)
        tau   :  delay (1 x nray)
        dod   :  direction of departure (nray x 2)
        doa   :  direction of arrival  (nray x 2)

        """

        bs.FUDAsignal.__init__(self,fGHz,alpha,tau,dod,doa)
    
    def __repr__(self):
        st =''
        st = st + 'freq :'+str(self.x[0])+' '+str(self.x[-1])+' '+str(len(self.x))+"\n"
        st = st + 'shape  :'+str(np.shape(self.y))+"\n"
        st = st + 'tau :'+str(min(self.tau0))+' '+str(max(self.tau0))+"\n"
        st = st + 'dist :'+str(min(0.3*self.tau0))+' '+str(max(0.3*self.tau0))+"\n"
        return(st)

    def info(self):
        """ display information

        """
        print 'Ftt,Ftp,Frt,Frp'
        print 'dod,doa,tau'
        print 'H - FUDsignal '
        print 'tau min , tau max :', min(self.tau), max(self.tau)
        self.H.info()

    def imshow(self):
        """ imshow vizualization of H

        """
        self.H
        sh = np.shape(self.H.y)
        itau = np.arange(len(self.tau))
        plt.imshow(abs(self.H.y))
        plt.show()

    def apply(self, W):
        """ Apply a FUsignal W to the ScalChannel.

        Parameters
        ----------
        W :  Bsignal.FUsignal

        It exploits multigrid convolution from Bsignal.

        Notes
        -----
            + W may have a more important number of points and a smaller frequency band.
            + If the frequency band of the waveform exceeds the one of the ScalChannei, a warning is sent.
            + W is a FUsignal whose shape doesn't need to be homogeneous with FUDsignal H

        """

        #H = self.H
        U = self * W
        V = bs.FUDAsignal(U.x, U.y, self.tau0,self.dod,self.doa)

        return(V)

    def applywavC(self, w, dxw):
        """ apply waveform method C

        Parameters
        ----------
        w     :
            waveform
        dxw

        Notes
        -----

        The overall received signal is built in time domain
        w is apply on the overall CIR

        """

        H = self.H
        h = H.ft1(500, 1)
        dxh = h.dx()
        if (abs(dxh - dxw) > 1e-10):
            if (dxh < dxw):
                # reinterpolate w
                f = interp1d(w.x, w.y)
                x_new = arange(w.x[0], w.x[-1], dxh)[0:-1]
                y_new = f(x_new)
                w = TUsignal(x_new, y_new)
            else:
                # reinterpolate h
                f = interp1d(h.x, h.y)
                x_new = arange(h.x[0], h.x[-1], dxw)[0:-1]
                y_new = f(x_new)
                h = bs.TUsignal(x_new, y_new)

        ri = h.convolve(w)
        return(ri)

    def chantap(self,**kwargs):


        defaults = {
                    'fcGHz':4.5,
                    'WGHz':1,
                    'Ntap':100
        }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value



        h = bs.FUDsignal(self.x,self.y,self.tau0)
        htap = h.chantap(**kwargs)
        return htap
 
        



    def applywavB(self, Wgam):
        """ apply waveform method B (time domain )

        Parameters
        ----------

        Wgam : waveform including gamma factor

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
        #
        # return a FUDsignal
        #
        Y = self.apply(Wgam)
        ri = Y.ft1(500, 1)

        return(ri)

    def applywavA(self, Wgam, Tw):
        """ apply waveform method A

        Parameters
        ----------

        Wgam :
        Tw   :

        The overall received signal is built in frequency domain

        """
        Hab = self.H.ft2(0.001)
        HabW = Hab * Wgam
        RI = HabW.symHz(10000)
        ri = RI.ifft(0, 'natural')
        ri.translate(-Tw)
        return(ri)

    def doddoa(self):
        """ doddoa() : DoD / DoA diagram

        """
        dod = self.dod
        doa = self.doa
        #
        #col  = 1 - (10*np.log10(Etot)-Emin)/(Emax-Emin)
        Etot = self.H.energy()
        Etot = Etot / max(Etot)
        al = 180 / np.pi
        col = 10 * np.log10(Etot)
        print len(dod[:, 0]), len(dod[:, 1]), len(col[:])
        plt.subplot(121)
        plt.scatter(dod[:, 0] * al, dod[:, 1] * al, s=15, c=col,
                    cmap=plt.cm.gray_r, edgecolors='none')
        a = colorbar()
        #a.set_label('dB')
        plt.xlabel("$\\theta_t(\degree)$", fontsize=18)
        plt.ylabel('$\phi_t(\degree)$', fontsize=18)
        title('DoD')
        plt.subplot(122)
        plt.scatter(doa[:, 0] * al, doa[:, 1] * al, s=15, c=col,
                    cmap=plt.cm.gray_r, edgecolors='none')
        b = colorbar()
        b.set_label('dB')
        plt.title('DoA')
        plt.xlabel("$\\theta_r(\degree)$", fontsize=18)
        plt.ylabel("$\phi_r (\degree)$", fontsize=18)
        plt.show()

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
        #r.require('graphics')
        #r.postscript('fig.eps')
        #r('par(mfrow=c(2,2))')
        #Y.fig(Nray)
        y = Y.iftd(100, 0, 50, 0)
        y.fig(Nray)
        #r.dev_off()
        #os.system("gv fig.eps ")
        #y.fidec()
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
        Hk = bs.FUDsignal(self.H.x, self.H.y[k, :])

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
                print " resampling w"
                x_new = arange(W.x[0], W.x[-1] + dxh, dxh)[0:-1]
                Wk = W.resample(x_new)
                dx = dxh
            else:
                # reinterpolate h
                print " resampling h"
                x_new = arange(Hk.x[0], Hk.x[-1] + dxw, dxw)[0:-1]
                Hk = Hk.resample(x_new)
                dx = dxw
                Wk = W

        # qHk.x[0]==Wk.x[0]

    def RSSI(self,ufreq=0) :
        """ Compute RSSI value from a 
        specific frequency of the transmission channel
    
        Parameters
        ----------
        
        ufreq : int
            index in the frequency range


        Returns
        -------
            RSSI: float
                RSSI value

        """
    
        Tk = np.real(self.y[:,ufreq])
        return(20*np.log(np.sum(Tk**2)))

if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
