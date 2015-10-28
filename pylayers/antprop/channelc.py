# -*- coding:Utf-8 -*-
"""
.. currentmodule:: pylayers.antprop.channelc

DEPRECATED

VectChannel Class
=================


.. autosummary::
   :toctree: generated/

   VectChannel.__init__
   VectChannel.show3_old
   VectChannel.show3

ScalChannel Class
=================

.. autosummary::
   :toctree: generated/

   ScalChannel.__init__
   ScalChannel.info
   ScalChannel.imshow
   ScalChannel.apply
   ScalChannel.applywavC
   ScalChannel.applywavB
   ScalChannel.applywavA
   ScalChannel.doddoa
   ScalChannel.wavefig
   ScalChannel.rayfig

VectLOS Class
=============

.. autosummary::
   :toctree: generated/

   VectLOS.__init__
   VectLOS.cir

"""

import doctest
import pdb
import numpy as np
import scipy as sp
import pylab as plt
import struct as stru
from pylayers.antprop.channel import *
import pylayers.util.pyutil as pyu
import pylayers.signal.bsignal as bs
import pylayers.util.geomutil as geu
from pylayers.antprop.raysc import GrRay3D
from pylayers.util.project import *

class VectChannel(Ctilde):
    """ container for a vector representation of the propagation channel

    Attributes
    -----------

    Ctt   FUsignal (Nray x Nf canal )
    Cpp
    Cpt
    Ctp

    built in vec2scal1

    Frt   Fusignal  (Nray x Nf antenna )
    Frp
    Ftt
    Ftp

    fGHz  : frequency
    tauk  : delay
    tang  : dod
    rang  : doa

    Methods
    -------

    init(S,itx,irx)
        S is a simulation object, itx and irx are index of tx and rx
    show(display=False,mode='linear')
        display vect channel
    doadod()
        scatter plot DoA - DoD
    vec2scal(fGHz)
        build scal channel without antenna
    vec2scal1(fGHz)
        build scal channel with antenna

    """

    def __init__(self, S, itx, irx, transpose=False):
        """ class constructor

        Parameters
        ----------
        S
             Simulation
        itx
            tx number
        irx
            rx number
        transpose
            antenna transposition indicator
        """
        # .. todo::
        #
        #   a verifier -ravel-
        self.fail = False
        _filefield = S.dfield[itx][irx]
        filefield = pyu.getlong(_filefield,pstruc['DIRTUD'])
        _filetauk = S.dtauk[itx][irx]
        filetauk = pyu.getlong(_filetauk,pstruc['DIRTUD'])
        _filetang = S.dtang[itx][irx]
        filetang = pyu.getlong(_filetang,pstruc['DIRTUD'])
        _filerang = S.drang[itx][irx]
        filerang = pyu.getlong(_filerang,pstruc['DIRTUD'])

        """
        .. todo::

              Revoir Freq
        """
        # old version
        #freq      = S.freq()
        #self.freq = freq
        self.fGHz = S.fGHz

        #
        # pour show3 de gr on a besoin de filetra et indoor
        # pas beau
        #
        self.filetra = S.dtra[itx][irx]

        self.L = S.L

        #try:
        #     fo = open(filetauk, "rb")
        #except:
        #     self.fail=True
        #     print "file ",filetauk, " is unreachable"
        # decode filetauk
        #if not self.fail:
        #      nray_tauk = unpack('i',fo.read(4))[0]
        #     print "nb rayons dans .tauk : ",nray_tauk
        #     buf = fo.read()
        #      fo.close()
        #      nray = len(buf)/8
        #     print "nb rayons 2: ",nray
        #      self.tauk = ndarray(shape=nray,buffer=buf)
        #      if nray_tauk != nray:
        #          print itx , irx
        #          print nray_tauk - nray
        #self.tauk = self.tauk

        Ctilde.__init__(self)
        self.load(filefield, transpose)

        # decode the angular files (.tang and .rang)
#         #try:
#                 fo = open(filetang, "rb")
#         except:
#             self.fail=True
#             print "file ",filetang, " is unreachable"
#          if not self.fail:
#              nray_tang = unpack('i',fo.read(4))[0]
#                 buf = fo.read()
#              fo.close()
#              # coorectif Bug evalfield
#              tmp  = ndarray(shape=(nray_tang,2),buffer=buf)
#              self.tang = tmp[0:nray,:]
#          try:
#             fo = open(filerang, "rb")
#         except:
#             self.fail=True
#             print "file ",filerang, " is unreachable"
#
#          if not self.fail:
#              nray_rang = stru.unpack('i',fo.read(4))[0]
#             buf  = fo.read()
#              fo.close()
#              # correctif Bug evalfield
#              tmp = ndarray(shape=(nray_rang,2),buffer=buf)
#              self.rang = tmp[0:nray,:]
        #sh  = shape(self.Ctt.y)

        """
         .. todo::

             Express Ftt and Ftp in global frame from Tt and ant_tx
             Express Frt and Frp in global frame from Tt and ant_tx

        """

        #self.Ftt = FUsignal(fGHz,np.ones(sh))
        #self.Ftp = FUsignal(fGHz,np.zeros(sh))
        #self.Frt = FUsignal(fGHz,np.ones(sh))
        #self.Frp = FUsignal(fGHz,np.zeros(sh))

    def show3_old(self, id=0):
        """ geomview visualization old version

        This function provides a complete ray tracing vsualization
        of the channel structure. The rays are color coded as a fonction
        of their energy.

        Parameters
        ----------

        id : int
            index of filetra

        """
        E = self.Ctt.energy() + self.Ctp.energy() + \
            self.Cpt.energy() + self.Cpp.energy()
        u = argsort(E)
        v = u[-1::-1]
        Es = E[v]

        gr = GrRay3D()
        gr.load(self.filetra, self.L)

        filename = pyu.getlong("grRay" + str(id) + "_col.list",pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        fo.write("{<strucTxRx.off}\n")

        Emax = Es[0]

        rayset = len(Emax)
        for i in range(rayset):
            j = v[i]
            r = gr.ray3d[j]
            col = np.array([1, 0, 0])    # red
            fileray = r.show3(False, False, col, j)
            fo.write("{< " + fileray + " }\n")

        k = i + 1
        rayset = len(where((Es >= 0.1 * Emax) & (Es < 0.5 * Emax))[0])
        for i in range(rayset):
            j = v[i + k]
            r = gr.ray3d[j]
            col = np.array([0, 0, 1])    # blue
            fileray = r.show3(False, False, col, j)
            fo.write("{< " + fileray + " }\n")

        k = i + 1
        rayset = len(where((Es >= 0.01 * Emax) & (Es < 0.1 * Emax))[0])
        for i in range(rayset):
            j = v[i + k]
            r = gr.ray3d[j]
            col = np.array([0, 1, 1])    # cyan
            fileray = r.show3(False, False, col, j)
            fo.write("{< " + fileray + " }\n")

        k = i + 1
        rayset = len(where((Es >= 0.001 * Emax) & (Es < 0.01 * Emax))[0])
        for i in range(rayset):
            j = v[i + k]
            r = gr.ray3d[j]
            col = np.array([0, 1, 0])    # green
            fileray = r.show3(False, False, col, j)
            fo.write("{< " + fileray + " }\n")

        k = i + 1
        rayset = len(where(Es < 0.001 * Emax)[0])
        for i in range(rayset):
            j = v[i + k]
            r = gr.ray3d[j]
            col = np.array([1, 1, 0])    # yellow
            fileray = r.show3(False, False, col, j)
            fo.write("{< " + fileray + " }\n")

        fo.close()

        chaine = "geomview " + filename + " 2>/dev/null &"
        os.system(chaine)

    def show3(self, seuildb=100):
        """ geomview vizualization

        This function provides a complete ray tracing visualization 
        of the radio channel. Rays are color coded as a fonction of 
        their energy.

        Parameters
        ----------

        seuildb : float
             default 100
        """
        E = self.Ctt.energy() + self.Ctp.energy() + \
            self.Cpt.energy() + self.Cpp.energy()
        u = argsort(E)
        v = u[-1::-1]
        Es = E[v]

        gr = GrRay3D()
        gr.load(self.filetra, self.L)

        filename = pyu.getlong("grRay" + str(seuildb) + "_col.list", pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        fo.write("{<strucTxRx.off}\n")

        Emax = Es[0]

        rayset = len(v)
        db = 20 * np.log10(Es)
        c = 1 - (db > -seuildb) * (db + seuildb) / seuildb
        app = round(np.log10(Es / Emax))
        lw = app - min(app)

        for i in range(rayset):
            j = v[i]
            r = gr.ray3d[j]
            col = np.array([c[i], c[i], c[i]])
            l = int(lw[i])
            fileray = r.show3(False, False, col, j, l)
            #fileray =r.show3(False,False,col,j)
            fo.write("{< " + fileray + " }\n")

        fo.close()

        chaine = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
        os.system(chaine)

    def _vec2scal(self):
        """ calculate scalChannel from VectChannel and antenna

        DEPRECATED
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

    def _vec2scalA(self, At, Ar, alpha=1.0):
        """ calculate scalChannel from Vectchannel

        DEPRECATED

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

        Ftt, Ftp = At.Fsynth3(self.rang[:, 0], self.rang[:, 1], pattern=False)

        Frt, Frp = Ar.Fsynth3(self.rang[:, 0], self.rang[:, 1], pattern=False)

        Ftt = Ftt.transpose()
        Ftp = Ftp.transpose()
        Frt = Frt.transpose()
        Frp = Frp.transpose()

        Ftt = bs.FUsignal(At.fGHz, Ftt * alpha)
        Ftp = bs.FUsignal(At.fGHz, Ftp * alpha)
        Frt = bs.FUsignal(Ar.fGHz, Frt * alpha)
        Frp = bs.FUsignal(Ar.fGHz, Frp * alpha)

        scalch = ScalChannel(self, Ftt, Ftp, Frt, Frp)
        return(scalch)


class ScalChannel(object):
    """
    DEPRECATED
    ScalChannel Class :

    The ScalChannel is obtained from combination of the propagation
    channel and the antenna transfer function from both transmitting
    and receiving antennas

    Members
    -------

    H    : FUDSignal
        ray transfer functions  (nray,nfreq)
    dod  :
        direction of depature (rad) [theta_t,phi_t]  nray x 2
    doa  :
        direction of arrival (rad)  [theta_r,phi_r]  nray x 2
    tauk :
        delay ray k in ns


    """
    def __init__(self, VC, Ftt, Ftp, Frt, Frp):

        self.Ftt = Ftt
        self.Ftp = Ftp
        self.Frt = Frt
        self.Frp = Frp

        t1 = VC.Ctt * Frt + VC.Cpt * Frp
        t2 = VC.Ctp * Frt + VC.Cpp * Frp
        t3 = t1 * Ftt + t2 * Ftp

        self.dod = VC.tang
        self.doa = VC.rang
        self.tau = VC.tauk

        self.H = bs.FUDsignal(t3.x, t3.y, VC.tauk)

        # thresholding of rays
        if (VC.nray > 1):
            indices = self.H.enthrsh()
            self.dod = self.dod[indices, :]
            self.doa = self.doa[indices, :]
            self.tau = self.tau[indices, :]

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

        H = self.H
        U = H * W
        V = bs.FUDsignal(U.x, U.y, H.tau0)

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
                h = TUsignal(x_new, y_new)

        ri = h.convolve(w)
        return(ri)

    def applywavB(self, Wgam):
        """ apply waveform method B (time domain )

        Parameters
        ----------
        Wgam :
            waveform including gamma factor

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
        #ri      = Y.ft1(500,0)
        # Le fftshift est activ
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
        """ doddoa() : DoD | DoA diagram

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

        # on s'arrange que Hk.x[0]==Wk.x[0]
#         if Wk.x[0]!=Hk.x[0]:
#             x=arange(Wk.x[0],Hk.x[0],dx)
#             if Hk.x[0]!=x[0]:
#                 Hk.x=hstack((x,Hk.x[1:]))
#                 nz=len(x)
#                 Hk.y=hstack((zeros(nz),Hk.y))
#             else:
#                 Hk.x=hstack((x,Hk.x[0:]))
#                 nz=len(x)
#                 Hk.y=hstack((zeros(nz),Hk.y))
#

        self.Hk = Hk
        self.Wk = Wk
        Rk = Hk * Wk
        self.Rk = Rk
        rk = Rk.iftshift()

        plot(rk.x, rk.y, col)

        return(rk)


class VectLOS(Ctilde):
    def __init__(self, d, fmin=2, fmax=11, Nf=180):
        self.tauk = np.array([d / 0.3])
        fGHz = np.linspace(fmin, fmax, Nf)
        c1 = 1.0 / d * np.ones(len(fGHz))
        c2 = zeros(len(fGHz))
        c1.reshape(1, Nf)
        c2.reshape(1, Nf)
        self.freq = freq
        self.Ctt = bs.FUsignal(fGHz, c1)
        self.Ctp = bs.FUsignal(fGHz, c2)
        self.Cpt = bs.FUsignal(fGHz, c2)
        self.Cpp = bs.FUsignal(fGHz, c1)
        self.tang = array([0])
        self.rang = array([0])
        self.nray = 1

    def  cir(self, wav):
        """ Channel Impulse Response

        Parameters
        ----------
        wav :

        """
        SCO = self.vec2scal()
        ciro = SCO.applywavB(wav.sfg)
        return(ciro)

