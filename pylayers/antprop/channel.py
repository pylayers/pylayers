# -*- coding:Utf-8 -*-
import doctest
import numpy as np
import scipy as sp
import pylab as plt
import struct as stru
import pylayers.util.pyutil as pyu
import pylayers.signal.bsignal as bs
import pylayers.util.geomutil as geu
from pylayers.antprop.rays import GrRay3D
from pylayers.util.project import *


class Ctilde(object):
    """ container for the 4 components of the polarimetric 1 ray channel

    Attributes
    ----------

    Ctt : FUsignal
    Ctp : FUsignal
    Cpt : FUsignal
    Cpp : FUsignal
    nfreq : int
        number of frequency point
    nray  : int
        number of rays

    """
    def __init__(self):
        """
            transpose == False   (nray,nfreq)
            transpose == True    (nfreq,nray)
        """
        self.fail = False

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

        Parameters
        ----------

        filefield  : string
        transpose  : boolean
            default False

        """
        filetauk = filefield.replace('.field', '.tauk')
        filetang = filefield.replace('.field', '.tang')
        filerang = filefield.replace('.field', '.rang')
        try:
            fo = open(filefield, "rb")
        except:
            print "file ", filefield, " is unreachable"
            exit()

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
        freq = np.linspace(2, 11, nfreq)
        self.Ctt = bs.FUsignal(freq, c11)
        self.Ctp = bs.FUsignal(freq, c12)
        self.Cpt = bs.FUsignal(freq, c21)
        self.Cpp = bs.FUsignal(freq, c22)
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
            print "nb rayons dans .tauk : ", nray_tauk
            buf = fo.read()
            fo.close()
            nray = len(buf) / 8
            print "nb rayons 2: ", nray
            self.tauk = np.ndarray(shape=nray, buffer=buf)
            if nray_tauk != nray:
                print itx, irx
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
            # coorectif Bug evalfield
            tmp = np.ndarray(shape=(nray_rang, 2), buffer=buf)
            self.rang = tmp[0:nray, :]

    def doadod(self, cmap=plt.cm.hot_r, s=30):
        """ doadod scatter plot

        Parameters
        -----------
        cmap : color map
        s    : float
            size (default 30)

        Summary
        --------

        scatter plot of the DoA-DoD channel structure
        the energy is colorcoded over all couples of DoA-DoD

        """
        dod = self.tang
        doa = self.rang
        # determine Energy in copol and crospol
        Eco, Ecross = self.energy()
        Etot = Eco + Ecross
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
        plt.axis((0, 180, 0, 360))
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        #a=colorbar()
        #for t in a.ax.get_yticklabels():
        #    t.set_fontsize(18)
        #a.set_label('dB')
        plt.xlabel("$\\theta_t(\degree)$", fontsize=14)
        plt.ylabel('$\phi(\degree)$', fontsize=14)
        #ylabel('$\phi_t(\degree)$',fontsize=18)
        plt.title('DoD', fontsize=18)
        plt.subplot(122)
        plt.scatter(doa[:, 0] * al, doa[:, 1] * al, s=30, c=col,
                    cmap=plt.cm.hot_r, edgecolors='none')
        plt.axis((0, 180, 0, 360))
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        #b=colorbar()
        #b.set_label('dB',fontsize=18)
        #for t in b.ax.get_yticklabels():
        #    t.set_fontsize(20)
        plt.xlabel("$\\theta_r(\degree)$", fontsize=14)
        plt.title('DoA', fontsize=18)
        plt.ylabel("$\phi_r (\degree)$", fontsize=14)

    def show(self, display=False, mode='linear'):
        """ show
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
        """ Calculate energy on co and cross channel

            Returns
            -------
            Eco    : Energy on co channel    tt + pp
            Ecross : Energy on cross channel tp + pt
        """
        ECtt = self.Ctt.energy(1)
        ECtp = self.Ctp.energy(1)
        ECpt = self.Cpt.energy(1)
        ECpp = self.Cpp.energy(1)
        Eco = ECtt + ECpp
        Ecross = ECtp + ECpt
        return Eco, Ecross

    def sort(self, tauk):
        """ sort Ctilde with respect of tauk

        Parameters
        ----------
        tauk :

        """
        u = argsort(tauk)

        self.Ctt = self.Ctt[u, :]
        self.Cpp = self.Cpp[u, :]
        self.Ctp = self.Ctp[u, :]
        self.Cpt = self.Cpt[u, :]

    def vec2scal(self):
        """ calculate scalChannel from VectChannel and antenna

        Returns
        -------

        slach : ScalChannel

        """
        freq = self.freq
        sh = np.shape(self.Ctt.y)
        Ftt = bs.FUsignal(freq, np.ones(sh))
        Ftp = bs.FUsignal(freq, np.zeros(sh))
        Frt = bs.FUsignal(freq, np.ones(sh))
        Frp = bs.FUsignal(freq, np.zeros(sh))
        scalch = ScalChannel(self, Ftt, Ftp, Frt, Frp)
        return(scalch)

    # Inclusion of realistic antenna behaviour
    # Specify transmitter antenna and receiver antenna

    def vec2scalA(self, At, Ar, alpha=1.0):
        """
        vec2scalA(self,At,Ar,alpha=1.0):

        At = transmitter antenna
        Ar = receiver antenna

        Calculate ScalChannel by combining the propagation channel VectChannel
        with realistic antennas transfer function

        alpha : normalization factor

        """

        Ftt, Ftp = At.Fsynth3(self.rang[:, 0], self.rang[:, 1])

        Frt, Frp = Ar.Fsynth3(self.rang[:, 0], self.rang[:, 1])

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
        """
        Info (Nf,Nray,shape(y))
        """
        print "Nf   :", len(self.freq)
        print "Nray :", len(self.tauk)
        print "shape y :", shape(self.Ctt.y)

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

    freq  : frequency
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
    vec2scal(freq)
        build scal channel without antenna
    vec2scal1(freq)
        build scal channel with antenna

    """

    def __init__(self, S, itx, irx, transpose=False):
        """
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
        self.freq = S.freq

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

        #self.Ftt = FUsignal(freq,np.ones(sh))
        #self.Ftp = FUsignal(freq,np.zeros(sh))
        #self.Frt = FUsignal(freq,np.ones(sh))
        #self.Frp = FUsignal(freq,np.zeros(sh))

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

    def mobility(self, v, dt):
        """ Modify channel for uniform mobility

        Parameters
        ----------
        v  : float
            velocity (m/s)
        dt : float
            delta t (s)

        Notes
        ------
        Calculate new channel field from the old one and v(terminal vitese) 
        and dt(time of deplacement)

        dt en s  (le temps d'observation entre deux position de Rx)
        v en m/s (vitesse de changement de Rx)

        Returns
        -------

        VC : modified VectChannel

        """

        c = 0.3  # m/ns vitesse de lumiere
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


def VCg2VCl(VCg, Tt, Tr):
    """ global reference frame to local reference frame

        Parameters
        ----------
        VCg : VectChannel global
        Tt  : Tx rotation matrix
        Tr  : Rx rotation matrix

        Returns
        -------
        VCl : VectChannel (local)

    """
    import copy

    VCl = copy.deepcopy(VCg)
    freq = VCl.freq
    Rt, tangl = BTB_tx(VCg.tang, Tt)
    Rr, rangl = BTB_rx(VCg.rang, Tr)

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

    VCl.Ctt = bs.FUsignal(freq, Cttl)
    VCl.Ctp = bs.FUsignal(freq, Ctpl)
    VCl.Cpt = bs.FUsignal(freq, Cptl)
    VCl.Cpp = bs.FUsignal(freq, Cppl)

    return VCl


class ScalChannel(object):
    """
    ScalChannel Class :

    The ScalChannel is obtained from combination of the propagation
    channel and the antenna transfer function from both transmitting
    and receiving antennas

    Members
    -------

    H    : FUDSignal
        Rays Transfer function  (nray,nfreq)
    dod  :
        Direction of depature (rad) [theta_t,phi_t]  nrayx2
    doa  :
        Direction of arrival (rad)  [theta_r,phi_r]  nrayx2
    tauk : Delay ray k


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
        sh = shape(self.H.y)
        itau = arange(len(self.tau))
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
            + If the frequency band of the waveform exceed the one of the ScalChannel a warning is sent.
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
            waveform including gamme factor

        Returns
        -------
        ri  : TUDsignal
            impulse response for each ray separately

        Notes
        ------

            The overall received signal is built in time domain

            Wgam is applied on each Ray Transfer function

        """
        #
        # return a FUDsignal
        #
        Y = self.apply(Wgam)
        #ri      = Y.ft1(500,0)
        # Le fftshift est activé
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
        freq = np.linspace(fmin, fmax, Nf)
        c1 = 1.0 / d * np.ones(len(freq))
        c2 = zeros(len(freq))
        c1.reshape(1, Nf)
        c2.reshape(1, Nf)
        self.freq = freq
        self.Ctt = bs.FUsignal(freq, c1)
        self.Ctp = bs.FUsignal(freq, c2)
        self.Cpt = bs.FUsignal(freq, c2)
        self.Cpp = bs.FUsignal(freq, c1)
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

if __name__ == "__main__":
    #import tkFileDialog
    #FD = tkFileDialog
    doctest.testmod()
############################
## test Channel, Simulation, Antenna
############################
#
#     #################################
#     # initialisation de la simulation
#     #################################
#
#     #
#     # choix des antennes
#     #
#     fileantt = FD.askopenfilename(filetypes = [("Fichiers vsh3","*.vsh3"),
#     ("All", "*")],
#      title="Please choose an antenna file",
#      initialdir=antdir)
#
#     fileantr = FD.askopenfilename(filetypes = [("Fichiers vsh3","*.vsh3"),
#     ("All", "*")],
#      title="Please choose an antenna file",
#      initialdir=antdir)
#
#     #
#     # creation d'objet Simulation
#     #
#     S=Simulation(fileantTx=_fileantt,fileantRx=_fileantr)
#
#     #
#     # choix du fichier field
#     #
#     filefield = FD.askopenfilename(filetypes = [("Fichiers field","*.field"),
#      ("All", "*")],
#      title="Please choose a field file",
#      initialdir=tuddir)
#
#     filetauk = filefield.replace('.field','.tauk')
#     filetang = filefield.replace('.field','.tang')
#     filerang = filefield.replace('.field','.rang')
#     filefreq = filefield.replace('.field','.freq')
#     # il faudra recuperer aussi filetra
#
#      S.filefield.append(getshort(filefield))
#     S.filetauk.append(getshort(filetauk))
#     S.filetang.append(getshort(filetang))
#     S.filerang.append(getshort(filerang))
#     #S.filetra.append(getshort(filetra))

#     fichsimulin = FD.askopenfilename(filetypes = [("Fichiers simul ","*.simul"),
#          ("All", "*")],
#          title="Please choose a simulation file",
#          initialdir=simuldir)
##          parent=root)
#     _fichsimulin = getshort(fichsimulin)
#     S = Simulation(_fichsimulin)
#
#     #
#     # Frequency range for Siradel measurements
#     #
#     f     = arange(2,11.05,0.05)
#
#     #################################
#     # definition de la forme d'onde a l'emission
#     #################################
#     Tw = 10
#     x  = arange(-Tw,Tw,0.005)
#     w  = EnImpulse(x,5,3,10)
#
#     W = w.ft()
#     f = W.x
#     ygamma = -1j*0.3/(4*np.pi*f)
#     gamm  = FUsignal(f,ygamma)
#     Wgam  = W*gamm
#     wgam  = Wgam.ift()
#
###       #
###       # plot waveform
###       #
###       figure()
###       plot(w.x,w.y,linewidth=3)
###       ylabel('(V)',fontsize=28)
###       xlabel('retard(ns)',fontsize=28)
###       xticks(fontsize=26)
###       yticks(fontsize=26)
###       title('w(t)',fontsize=30)
###       grid(linewidth=2)
###       show()
##
##      ################################
##      # recuperation du canal de propagation simule
##      ################################
##
##      VC   = VectChannel(S,0,False)
##      #
##      # If file from Siradel
##      #
##      #VC   = VectChannel(S,0,True)
##
##      ################################
##      # calcul du canal de transmission
##      ################################
##
##      #
##      # avec antennes omnis
##      #
##      SC      = VC.vec2scal()
##
##      #
##      # avec antennes simulees
##      #
##      VCA     = VC
##      SCA     = VCA.vec2scalA(S.tx.A,S.rx.A)
##
##      ################################
##      # calcul du signal recu
##      ################################
##
##      #
##      # avec antennes omnis
##      #
##      r  = SC.applywavB(Wgam)
##
##      #
##      # avec antennes simulees
##      #
##      rA = SCA.applywavB(Wgam)
##
##      ################################
##      # visualisation de resultats
##      ################################
##
##      figure()
##      # pour monocone on a un gain max de 8.59**2
##      #plot(r.x,r.y*70,'b-')
##      # pour dipole on a un gain max de 1.68**2
##      #plot(r.x,r.y*2.8,'b-')
##
##      #
##      # calcul de la valeur maximale du gain d'antenne
##      #
##
##      gain_At=np.sqrt((abs(SCA.Ftt.y))**2+(abs(SCA.Ftp.y))**2)
##      gain_Ar=np.sqrt((abs(SCA.Frt.y))**2+(abs(SCA.Frp.y))**2)
##
##      G=gain_At.max()*gain_Ar.max()
##
##      figure()
##      plot(r.x,r.y*G,'b-')
##      plot(rA.x,rA.y,'r-')
##      ax=gca()
##      ax.ticklabel_format(style='sci',axis='y')
##      ax.yaxis.major.formatter.set_powerlimits((0,0))
##      xlabel('retard (ns)',fontsize=28)
##      ylabel('(V)',fontsize=28)
##      xticks(fontsize=24)
##      yticks(fontsize=24)
##      title('r(t)',fontsize=30)
##      show()
##
###       #legend(('antennes omnis','antennes simulees'))
###       show()
##
################################
### test changement d'antennes
################################
##      #S.tx.gantenna()
##      #S.rx.gantenna()
##
##      #SCA   = VCA.vec2scalA(S.tx.A,S.rx.A)
##      #rA = SCA.applywavB(Wgam)
##
##      #gain_At=np.sqrt((abs(SCA.Ftt.y))**2+(abs(SCA.Ftp.y))**2)
##      #gain_Ar=np.sqrt((abs(SCA.Frt.y))**2+(abs(SCA.Frp.y))**2)
##
##      #G=gain_At.max()*gain_Ar.max()
##
##      #figure()
##      #plot(r.x,r.y*G,'b-')
##      #plot(rA.x,rA.y,'r-')
##      #xlabel('retard (ns)',fontsize=28)
##      #ylabel('(V)',fontsize=28)
##      #xticks(fontsize=24)
##      #yticks(fontsize=24)
##      #title('r(t)',fontsize=30)
##
##      #legend(('ideale','dipole'))
##      #show()
##
##
#
###############################
## test passage Rg -> Rl
###############################
#     itx=0
#     irx=0
#     VCg = VectChannel(S,itx,irx,False)
#     alpha = np.pi*rand()
#     beta  = np.pi*rand()
#     gamma = np.pi*rand()
#     Ta = MEulerAngle(alpha,beta,gamma)
#     alpha = np.pi*rand()
#     beta  = np.pi*rand()
#     gamma = np.pi*rand()
#     Tb = MEulerAngle(alpha,beta,gamma)
#
#     VCl=VCg2VCl(VCg,Ta,Tb)
#
#     SCl  = VCl.vec2scal()
#
#     VClA  = VCl
#     SClA  = VClA.vec2scalA(S.tx.A,S.rx.A)
#
#     rl = SCl.applywavB(Wgam)
#     rlA = SClA.applywavB(Wgam)
#
#     gain_At=np.sqrt((abs(SClA.Ftt.y))**2+(abs(SClA.Ftp.y))**2)
#     gain_Ar=np.sqrt((abs(SClA.Frt.y))**2+(abs(SClA.Frp.y))**2)
#     G=gain_At.max()*gain_Ar.max()
#
#     figure()
#     plot(rl.x,rl.y*G,'b-')
#     plot(rlA.x,rlA.y,'r-')
#     ax=gca()
#     ax.ticklabel_format(style='sci',axis='y')
#     ax.yaxis.major.formatter.set_powerlimits((0,0))
#     xlabel('retard (ns)',fontsize=28)
#     ylabel('(V)',fontsize=28)
#     xticks(fontsize=24)
#     yticks(fontsize=24)
#     title('r(t)',fontsize=30)
#     show()
#
