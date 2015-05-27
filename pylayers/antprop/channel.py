# -*- coding:Utf-8 -*-
"""
.. currentmodule:: pylayers.antprop.channel

.. automodule::
    :members:

.. autoclass:: Ctilde
    :members:

.. autoclass:: Tchannel
    :members:

"""
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
try:
    import h5py
except:
    print 'h5py is not installed: Ctilde(object cannot be saved)'

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
        self.fail = False
        self.islocal = False
        self.Tt = np.eye(3)
        self.Tr = np.eye(3)

    def __repr__(self):
        s = 'Ctilde'+'\n---------\n'
        if hasattr(self, 'Cpp'):
            s = s + str(np.shape(self.Cpp.y))+'\n'
        if hasattr(self, 'nray'):
            s = s + 'Nray : ' + str(self.nray)+'\n'
            s = s + 'fmin(GHz) : ' + str(self.Cpp.x[0])+'\n'
            s = s + 'fmax(GHz): ' + str(self.Cpp.x[-1])+'\n'
            s = s + 'Nfreq : ' + str(self.nfreq)+'\n'
        return(s)

    def choose(self):
        """ Choose a field file in tud directory

        DEPRECATED 
        """
        import tkFileDialog as FD
        filefield = FD.askopenfilename(filetypes=[("Files field  ", "*.field"),
                                                  ("All", "*")],
                                       title="Please choose a .field file",
                                       initialdir=tuddir)
        self.load(filefield, transpose=False)


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
        if self.islocal:
            self.locbas(b2g=True)

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f=h5py.File(filename,'w')
            f.create_dataset('Tt',shape=np.shape(self.Tt),data=self.Tt)
            f.create_dataset('Tr',shape=np.shape(self.Tr),data=self.Tr)
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

            self.Tt = f['Tt'][:]
            self.Tr = f['Tr'][:]

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

        if self.islocal:
            self.locbas(b2g=True)


        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        # try/except to avoid loosing the h5 file if 
        # read/write error
        try:

            fh5=h5py.File(filename,'a')
            if not grpname in fh5['Ct'].keys(): 
                fh5['Ct'].create_group(grpname)
            else :
                print 'Warning : Ct/'+grpname +'already exists in '+filenameh5
            f=fh5['Ct/'+grpname]

            # save channel in global basis
            f.create_dataset('Tt',shape=np.shape(self.Tt),data=self.Tt)
            f.create_dataset('Tr',shape=np.shape(self.Tr),data=self.Tr)
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

    def los(self,pa=np.r_[0,0,0],pb=np.r_[3,0,0],fGHz=np.r_[2.4],tang=np.r_[[1.57,3.14]],rang=np.r_[[1.57,0]]):
        """ Line of site channel

        Parameters
        ----------

        d (m) 
        fGHz (,Nf)
        tang (1x2)
        rang (1x2)

        """
        self.pa = pa
        self.pb = pb
        self.fGHz = fGHz
        self.nfreq = len(fGHz)
        self.nray = 1
        si = pb-pa
        d = np.r_[np.sqrt(np.sum(si*si))]
        self.tauk = d/0.3
        tht =  np.arccos(si[2])
        pht =  np.arctan2(si[1],si[0])
        thr =  np.arccos(-si[2])
        phr =  np.arctan2(-si[1],-si[0])
        self.tang = np.array([tht,pht]).reshape((1,2))
        self.rang = np.array([thr,phr]).reshape((1,2))
        U = np.ones(len(self.fGHz),dtype=complex)/d[0]
        Z = np.zeros(len(self.fGHz),dtype=complex)
        self.Ctt = bs.FUsignal(self.fGHz, U)
        self.Ctp = bs.FUsignal(self.fGHz, Z)
        self.Cpt = bs.FUsignal(self.fGHz, Z)
        self.Cpp = bs.FUsignal(self.fGHz, U)

    def _loadh5(self,filenameh5,grpname):
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

            self.Tt = f['Tt'][:]
            self.Tr = f['Tr'][:]
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

    def load(self, filefield, transpose=False):
        """ load a Ctilde from a .field file

        Load the three files .tauk .tang .rang which contain respectively
        delay , angle of departure , angle of arrival.
maicher
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
        >>> out = C.load(pyu.getlong(S.dtud[1][1],'output'))

        """
        filetauk = filefield.replace('.field', '.tauk')
        filetang = filefield.replace('.field', '.tang')
        filerang = filefield.replace('.field', '.rang')
        try:
            fo = open(filefield, "rb")
        except:
            raise NameError( "file " + filefield + " is unreachable")

        # decode filename (*.field file obtained from evalfield simulation)
        nray = stru.unpack('i', fo.read(4))[0]
        nfreq = stru.unpack('i', fo.read(4))[0]
        if nfreq == 0:
            print " Error : incorrect number of frequency points in .field"
            self.fail = True
            return

        n = nray * nfreq
        buf = fo.read()
        fo.close()

        CMat = np.ndarray(shape=(n, 8), buffer=buf)
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
            #print "nb rays in .tauk file: ", nray_tauk
            buf = fo.read()
            fo.close()
            nray = len(buf) / 8
            #print "nb rays 2: ", nray
            self.tauk = np.ndarray(shape=nray, buffer=buf)
            #if nray_tauk != nray:
            #    print nray_tauk - nray
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
            self.tang = tmp[0:nray,:]
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
            self.rang = tmp[0:nray,:]

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
                    'fontsize':12,
                    'edgecolors':'none',
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
            print "len(col):", len(col)
            print "len(di):", len(dir)
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

        Summary
        --------

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


    def locbas(self, Tt=[], Tr=[],b2g=False):
        """ global reference frame to local reference frame

        If Tt and Tr are [] the global channel is  retrieved

        Parameters
        ----------

        Tt  : Tx rotation matrix 3x3
            default []
        Tr  : Rx rotation matrix 3x3
            default []
        b2g: bool
            back to global reference frame

        Returns
        -------

        Cl : Ctilde local/global
            depends on self.islocal boolean value

        Examples
        --------

        """
        # get frequency axes

        fGHz = self.fGHz
        # if rot matrices are passed
        if (Tt != []) & (Tr != []):
            if self.islocal:
                if (hasattr(self,'Tt')) & (hasattr(self,'Tr')):
                    # run locbas to return to global basis
                    self.locbas(b2g=True)
                else:
                    raise NameError('Channel has no self.Tt or self.Tr')
            self.Tt = Tt
            self.Tr = Tr
            self.islocal = True
        # if a return to global is requested
        elif b2g:
            if self.islocal :
                if (hasattr(self,'Tt')) & (hasattr(self,'Tr')):
                    self.Tt = self.Tt.transpose()
                    self.Tr = self.Tr.transpose()
                    self.islocal = False
                else:
                    raise NameError ('self.Tt and self.Tr should exist')
            else:
                print "nothing to do to return in global basis"
                return self
        # if Tt and Tr == []
        else:
            return self

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


        Rt, tangl = geu.BTB_tx(self.tang, self.Tt)
        Rr, rangl = geu.BTB_rx(self.rang, self.Tr)
        #
        # update direction of departure and arrival
        #

        self.tangl = tangl
        self.rangl = rangl

        #uf = np.ones(self.nfreq)

        #
        # r0 : r x 1(f)
        #

        #r0 = np.outer(Rr[0, 0,:], uf)
        r0 = Rr[0,0,:][:, np.newaxis]
        #r1 = np.outer(Rr[0, 1,:], uf)
        r1 = Rr[0,1,:][:, np.newaxis]

        t00 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t01 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rr[1, 0,:], uf)
        r0 = Rr[1, 0,:][:, np.newaxis]
        #r1 = np.outer(Rr[1, 1,:], uf)
        r1 = Rr[1, 1,:][:, np.newaxis]

        t10 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t11 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rt[0, 0,:], uf)
        r0 = Rt[0, 0, :][:, np.newaxis]
        #r1 = np.outer(Rt[1, 0,:], uf)
        r1 = Rt[1, 0, :][:, np.newaxis]

        Cttl = t00 * r0 + t01 * r1
        Cptl = t10 * r0 + t11 * r1

        #r0 = np.outer(Rt[0, 1,:], uf)
        r0 = Rt[0, 1, :][:, np.newaxis]
        #r1 = np.outer(Rt[1, 1,:], uf)
        r1 = Rt[1, 1, :][:, np.newaxis]

        Ctpl = t00 * r0 + t01 * r1
        Cppl = t10 * r0 + t11 * r1

        self.Ctt = bs.FUsignal(fGHz, Cttl)
        self.Ctp = bs.FUsignal(fGHz, Ctpl)
        self.Cpt = bs.FUsignal(fGHz, Cptl)
        self.Cpp = bs.FUsignal(fGHz, Cppl)

        return self

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
            self.Tt = Tt
            self.Tr = Tr
        else:
            if (hasattr(self,'Tt')) & (hasattr(self, 'Tr')):
                self.Tt = self.Tt.transpose()
                self.Tr = self.Tr.transpose()
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

        Rt, tangl = geu.BTB_tx(self.tang, self.Tt)
        Rr, rangl = geu.BTB_rx(self.rang, self.Tr)

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
        r0 = Rr[0,0,:][:,np.newaxis]
        #r1 = np.outer(Rr[0, 1,:], uf)
        r1 = Rr[0,1,:][:,np.newaxis]

        t00 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t01 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rr[1, 0,:], uf)
        r0 = Rr[1, 0,:][:,np.newaxis]
        #r1 = np.outer(Rr[1, 1,:], uf)
        r1 = Rr[1, 1,:][:,np.newaxis]

        t10 = r0 * self.Ctt.y + r1 * self.Cpt.y
        t11 = r0 * self.Ctp.y + r1 * self.Cpp.y

        #r0 = np.outer(Rt[0, 0,:], uf)
        r0 = Rt[0,0,:][:,np.newaxis]
        #r1 = np.outer(Rt[1, 0,:], uf)
        r1 = Rt[1,0,:][:,np.newaxis]

        Cttl = t00 * r0 + t01 * r1
        Cptl = t10 * r0 + t11 * r1

        #r0 = np.outer(Rt[0, 1,:], uf)
        r0 = Rt[0,1,:][:,np.newaxis]
        #r1 = np.outer(Rt[1, 1,:], uf)
        r1 = Rt[1,1,:][:,np.newaxis]

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
        ax1.set_xlabel('f (GHz)',fontsize=kwargs['fontsize'])
        ax1.set_title(u'$C_{\\theta\\theta}$',fontsize=kwargs['fontsize'])

        ax2 = kwargs['fig'].add_subplot(222)
        fig, ax2 = self.Ctp.imshow(ax=ax2,**kwargs)
        ax2.set_xlabel('f (GHz)',fontsize=kwargs['fontsize'])
        ax2.set_title(u'$C_{\\theta\phi}$',fontsize=kwargs['fontsize'])

        ax3 = kwargs['fig'].add_subplot(223)
        fig, ax3 = self.Cpt.imshow(ax=ax3,**kwargs)
        ax3.set_xlabel('f (GHz)',fontsize=kwargs['fontsize'])
        ax3.set_title(u'$C_{\phi\\theta}$',fontsize=kwargs['fontsize'])

        ax4 = kwargs['fig'].add_subplot(224)
        fig, ax4 = self.Cpp.imshow(ax=ax4,**kwargs)
        ax4.set_xlabel('f (GHz)',fontsize=kwargs['fontsize'])
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
        assert np.allclose(self.tauk, C.tauk)
        for r in range(self.nray):
            if np.allclose(self.Ctt.y[r,:], C.Ctt.y[r,:]):
                print r

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

    def cut(self,threshold=0.99):
        """ cut rays from a energy threshold

        Parameters
        ----------

        threshold : float
            default 0.99

        """
        Ett, Epp, Etp, Ept = self.energy()
        Etot = Ett+Epp+Etp+Ept
        u = np.argsort(Etot)
        cumE = np.cumsum(Etot)/sum(Etot)
        v = np.where(cumE<threshold)[0]

        self.tauk = self.tauk[v]
        self.tang = self.tang[v,:]
        self.rang = self.rang[v,:]
        self.Ctt.y = self.Ctt.y[v,:]
        self.Cpp.y = self.Cpp.y[v,:]
        self.Ctp.y = self.Ctp.y[v,:]
        self.Cpt.y = self.Cpt.y[v,:]

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

    def prop2tran(self,a='theta',b='theta',Friis=True,debug=True):
        r""" transform propagation channel into transmission channel

        Parameters
        ----------

        a : string or antenna array
            polarization antenna a ( 'theta' | 'phi' | 'ant' )
            default (theta)
        b : string or antenna array
            polarization antenna b ( 'theta' | 'phi' | 'ant' )
            default (theta)

            0 : theta
            1 : phi

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

        H : Tchannel(bs.FUDAsignal)


        """
        freq = self.fGHz
        nfreq = self.nfreq
        nray  = self.nray
        sh = np.shape(self.Ctt.y)
        if type(a) == str:

            if a == 'theta':
                Fat = np.ones((nray, nfreq))
                Fap = np.zeros(nray*nfreq).reshape((nray, nfreq))

            if a == 'phi':
                Fat = np.zeros(nray*nfreq).reshape((nray, nfreq))
                Fap = np.ones((nray, nfreq))

            Fat = bs.FUsignal(self.fGHz, Fat)
            Fap = bs.FUsignal(self.fGHz, Fap)

        else:
            if a.fromfile :
                Fat, Fap = a.Fsynth3(self.tangl[:, 0], self.tangl[:, 1], pattern=False)
                Fat = Fat.transpose()
                Fap = Fap.transpose()
                Fat = bs.FUsignal(a.fa, Fat)
                Fap = bs.FUsignal(a.fa, Fap)
            else:
                Fat, Fap = a.Fpatt(self.tangl[:, 0], self.tangl[:, 1], pattern=False)
                Fat = bs.FUsignal(a.fa, Fat)
                Fap = bs.FUsignal(a.fa, Fap)


        if type(b) == str:

            if b == 'theta':
                Fbt = np.ones((nray, nfreq))
                Fbp = np.zeros(nray*nfreq).reshape((nray, nfreq))
            if b == 'phi':
                Fbp = np.ones((nray, nfreq))
                Fbt = np.zeros(nray*nfreq).reshape((nray, nfreq))

            Fbt = bs.FUsignal(self.fGHz, Fbt)
            Fbp = bs.FUsignal(self.fGHz, Fbp)
        else:

            if b.fromfile :
                Fbt, Fbp = b.Fsynth3(self.rangl[:, 0], self.rangl[:, 1], pattern=False)
                Fbt = Fbt.transpose()
                Fbp = Fbp.transpose()
                Fbt = bs.FUsignal(b.fa, Fbt)
                Fbp = bs.FUsignal(b.fa, Fbp)
            else:
                Fbt, Fbp = b.Fpatt(self.rangl[:, 0], self.rangl[:, 1], pattern=False)
                Fbt = bs.FUsignal(b.fa, Fbt)
                Fbp = bs.FUsignal(b.fa, Fbp)
        # Ctt : r x f
        # Cg2cl should be applied here
        #

        #
        #  C  = 2 x 2 x r x f
        #  Fa = 2 x r x f
        #  Fb = 2 x r x f
        #t1 = self.Ctt * Fat + self.Cpt * Fap
        #t2 = self.Ctp * Fat + self.Cpp * Fap
        t1 = self.Ctt * Fat + self.Ctp * Fap
        t2 = self.Cpt * Fat + self.Cpp * Fap
        alpha = t1 * Fbt + t2 * Fbp



        H = Tchannel(alpha.x, alpha.y, self.tauk, self.tang, self.rang)

        if debug :
            H.alpha=alpha
            H.Fat=Fat
            H.Fap=Fap
            H.Fbt=Fbt
            H.Fbp=Fbp
            H.Gat=10*np.log10(np.sum(Fat.y*np.conj(Fat.y),axis=1)/len(Fat.x))
            H.Gap=10*np.log10(np.sum(Fap.y*np.conj(Fap.y),axis=1)/len(Fap.x))
            H.Gbt=10*np.log10(np.sum(Fbt.y*np.conj(Fbt.y),axis=1)/len(Fbt.x))
            H.Gbp=10*np.log10(np.sum(Fbp.y*np.conj(Fbp.y),axis=1)/len(Fbp.x))

        if Friis:
            H.applyFriis()

        # average w.r.t frequency
        Nf   = H.y.shape[1]
        H.ak = np.real(np.sqrt(np.sum(H.y * np.conj(H.y)/Nf, axis=1)))
        H.tk = H.taud

        return(H)

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



class Tchannel(bs.FUDAsignal):
    """ Handle the transmission channel

    The transmission channel TChannel is obtained through combination of the propagation
    channel and the antenna transfer functions from both transmitter and receiver.

    Members
    -------

        ray transfer functions  (nray,nfreq)
    dod  :
        direction of depature (rad) [theta_t,phi_t]  nray x 2
    doa  :
        direction of arrival (rad)  [theta_r,phi_r]  nray x 2
    tauk :
        delay ray k in ns

    Methods
    -------

    imshow()
    apply(W)
    applywavB(Wgam,Tw)
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
                fGHz = np.arange(0,2,1),
                alpha= np.arange(0,2,1),
                tau  = np.array(([],)),
                dod  = np.array(([[],[]])).T,
                doa  = np.array(([[],[]])).T):
        """ class constructor

        Parameters
        ----------

        fGHz  :  1 x nfreq
            frequency GHz
        alpha :  nray x nfreq
            path amplitude
        tau   :  1 x nray
            path delay (ns)
        dod   :  direction of departure (nray x 2)
        doa   :  direction of arrival   (nray x 2)

        """
        super(Tchannel,self).__init__(fGHz, alpha, tau, dod, doa)

    def __repr__(self):
        st = ''
        st = st + 'freq : '+str(self.x[0])+' '+str(self.x[-1])+' '+str(len(self.x))+"\n"
        st = st + 'shape  : '+str(np.shape(self.y))+"\n"
        st = st + 'tau (min, max) : '+str(min(self.taud))+' '+str(max(self.taud))+"\n"
        st = st + 'dist :'+str(min(0.3*self.taud))+' '+str(max(0.3*self.taud))+"\n"
        if self.isFriis:
            st = st + 'Friis factor -j c/(4 pi f) has been applied'
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
                postion of point a (transmitter)
            b = np.ndarray
                postion of point b (receiver)
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
                print 'Warning : H/'+grpname +'already exists in '+filenameh5
            f=fh5['H/'+grpname]

            for k,va in self.__dict__.items():
                f.create_dataset(k,shape = np.shape(va),data=va)
            fh5.close()
        except:
            fh5.close()
            raise NameError('Channel Tchannel: issue when writting h5py file')

    def _loadh5(self,filenameh5,grpname):
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


    def apply(self, W):
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
            Transmission Channel, a warning is sent.
            + W is a FUsignal whose shape doesn't need to be homogeneous with FUDsignal H

        """

        U = self * W
        V = bs.FUDAsignal(U.x, U.y, self.taud, self.dod, self.doa)

        return(V)

    def applywavC(self, w, dxw):
        """ apply waveform method C

        Parameters
        ----------
        w :
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
        """ channel tap

        Parameters
        ----------

        fcGHz :
            WGHz  :
        Ntap  :

        """

        defaults = {'fcGHz':4.5,
                    'WGHz':1,
                    'Ntap':100}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value



        h = bs.FUDsignal(self.x, self.y, self.taud)
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
        # return a TUsignal
        #
        #import ipdb
        #ipdb.set_trace()
        Y = self.apply(Wgam)
        ri = Y.ft1(Nz=500,ffts=1)

        return(ri)

    def applywavA(self, Wgam, Tw):
        """ apply waveform method A

        Parameters
        ----------

        Wgam :
        Tw   :

        The overall received signal is built in frequency domain

        See Also
        --------

        pylayers.signal.bsignal

        """
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
                    'fontsize':12,
                    'edgecolors':'none',
                    'polar':False,
                    'colorbar':False,
                    'title':False
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value


        di = getattr(self, d, 'doa')


        # remove non plt.scatter kwargs
        phi = kwargs.pop('phi')
        the = (0,180)
        fontsize = kwargs.pop('fontsize')
        polar = kwargs.pop('polar')
        fig = kwargs.pop('fig')
        ax = kwargs.pop('ax')
        colorbar = kwargs.pop('colorbar')
        reverse = kwargs.pop('reverse')
        normalize = kwargs.pop('normalize')
        mode =kwargs.pop('mode')
        title =kwargs.pop('title')
        xa = kwargs.pop('xa')
        xb = kwargs.pop('xb')
        if fig == []:
            fig = plt.figure()


        Etot = self.energy(mode=mode) + 1e-15


        if normalize:
            Emax = max(Etot)
            Etot = Etot / Emax
        #
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

        col = 10 * np.log10(Etot)
        kwargs['c'] = col
        if len(col) != len(di):
            print "len(col):", len(col)
            print "len(di):", len(dir)
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
            if normalize:
                b.set_label('dB')
            else:
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



        if self.isFriis :
            Etot = self.energy(mode=mode,Friis=False) + 1e-15
        else :
            Etot = self.energy(mode=mode,Friis=True) + 1e-15


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
            print "len(col):", len(col)
            print "len(di):", len(dir)
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

        Summary
        --------

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

    def energy(self,mode='mean',Friis=True,sumray=False):
        """ calculates channel energy including antennas spatial filtering

        Parameters
        ----------

        mode : string
            center | mean | integ    (different manner to get the value)
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
        Etot = bs.FUsignal.energy(self,axis=1,mode=mode,Friis=Friis)
        if sumray:
            Etot = np.sum(Etot,axis=0)
        return Etot




    #def doadod(self, cmap=plt.cm.hot_r, s=30,fontsize = 12,phi=(0, 360),norm=False,polar=False):
    # def doadod(self,**kwargs):
    #     """ doadod scatter plot

    #     Parameters
    #     -----------

    #     cmap : color map
    #     s    : float
    #         size (default 30)
    #     fontsize : integer
    #         default 12

    #     Summary
    #     --------

    #     scatter plot of the DoA-DoD channel structure
    #     the energy is colorcoded over all couples of DoA-DoD

    #     """

    #     defaults = {'cmap' : plt.cm.hot_r,
    #                 's': 30,
    #                 'fontsize' : 12,
    #                 'phi':(-180,180),
    #                 'normalize':False,
    #                 'polar':False,
    #                 'mode':'center'}

    #     for k in defaults:
    #         if k not in kwargs:
    #             kwargs[k] = defaults[k]

    #     args = {}
    #     for k in kwargs:
    #         if k not in defaults:
    #             args[k] = kwargs[k]
                
    #     the = (0,180)
    #     dod = self.dod
    #     doa = self.doa
    #     #
    #     # determine Energy in each channel
    #     #
    #     Etot = self.energy(mode=kwargs['mode']) +1e-15

    #     # normalization
    #     if kwargs['normalize']:
    #         Emax = max(Etot)
    #         Etot = Etot / Emax

    #     Emax = max(10 * np.log10(Etot))
    #     Emin = min(10 * np.log10(Etot))

    #     #
    #     #
    #     #
    #     # col  = 1 - (10*log10(Etot)-Emin)/(Emax-Emin)
    #     #

    #     if kwargs['polar'] :
    #         al = 1.
    #         phi=np.array(phi)
    #         the=np.array(the)            
    #         phi[0] = phi[0]*np.pi/180
    #         phi[1] = phi[1]*np.pi/180
    #         the[0] = the[0]*np.pi/180
    #         the[1] = the[1]*np.pi/180
    #     else :
    #         al = 180./np.pi

    #     col = 10 * np.log10(Etot)

    #     if len(col) != len(dod):
    #         print "len(col):", len(col)
    #         print "len(dod):", len(dod)

    #     plt.subplot(121, polar=kwargs['polar'])
    #     if kwargs['reverse']:
    #         plt.scatter(dod[:, 1] * al, dod[:, 0] * al,
    #                     s=kwargs['s'], c=col,
    #                     cmap=kwargs['cmap'],
    #                     edgecolors='none')
    #         plt.axis((kwargs['phi'][0], kwargs['phi'][1],the[0],the[1]))
    #         plt.xlabel('$\phi(^{\circ})$', fontsize=kwargs['fontsize'])
    #         plt.ylabel("$\\theta_t(^{\circ})$", fontsize=kwargs['fontsize'])
    #     else:
    #         plt.scatter(dod[:, 0] * al, dod[:, 1] * al,
    #                     s=kwargs['s'], c=col,
    #                     cmap=kwargs['cmap'],
    #                     edgecolors='none')
    #         plt.axis((the[0], the[1], kwargs['phi'][0], kwargs['phi'][1]))
    #         plt.xlabel("$\\theta_t(^{\circ})$", fontsize=kwargs['fontsize'])
    #         plt.ylabel('$\phi(^{\circ})$', fontsize=kwargs['fontsize'])
    #     # ylabel('$\phi_t(^{\circ})$',fontsize=18)
    #     plt.title('DoD', fontsize=kwargs['fontsize']+2)


    #     plt.subplot(122, polar=kwargs['polar'])
    #     if kwargs['reverse']:
    #         plt.scatter(doa[:, 1] * al, doa[:, 0] * al, s=30, c=col,
    #                     cmap=plt.cm.hot_r, edgecolors='none')
    #         plt.axis((kwargs['phi'][0], kwargs['phi'][1],the[0],the[1]))
    #         plt.xlabel("$\phi_r (^{\circ})$", fontsize=kwargs['fontsize'])
    #         plt.ylabel("$\\theta_r(^{\circ})$", fontsize=kwargs['fontsize'])
            
    #     else :
    #         plt.scatter(doa[:, 0] * al, doa[:, 1] * al, s=30, c=col,
    #                     cmap=plt.cm.hot_r, edgecolors='none')
    #         plt.axis((the[0], the[1], kwargs['phi'][0], kwargs['phi'][1]))
    #         plt.xlabel("$\\theta_r(^{\circ})$", fontsize=kwargs['fontsize'])
    #         plt.ylabel("$\phi_r (^{\circ})$", fontsize=kwargs['fontsize'])

    #     plt.title('DoA', fontsize=kwargs['fontsize']+2)

    #     # plt.xticks(fontsize=20)
    #     # plt.yticks(fontsize=20)
    #     b = plt.colorbar()
    #     if kwargs['normalize']:
    #         b.set_label('dB')
    #     else:
    #         b.set_label('Path Loss (dB)')
    #     # for t in b.ax.get_yticklabels():
    #     #    t.set_fontsize(20)
    #     plt.xlabel("$\\theta_r(^{\circ})$", fontsize=kwargs['fontsize'])
    #     plt.title('DoA', fontsize=kwargs['fontsize']+2)
    #     plt.ylabel("$\phi_r (^{\circ})$", fontsize=kwargs['fontsize'])



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

    def rssi(self,ufreq=0) :
        """ Compute RSSI value for a frequency index

        Parameters
        ----------

        ufreq : int
            index in the frequency range


        Returns
        -------

        RSSI: float
        RSSI value in dB

        Notes
        -----

        This function will be deprecated by energy function

        """

        Ak   = self.y[:, ufreq]
        #tauk = np.abs(self.x[:, ufreq])
        Pr   = np.sum(Ak*np.conj(Ak))
        Prp  = np.sum(Ak)*np.conj(np.sum(Ak))
        #rssiF   = np.sum(Ak*exp(-2*1j*,)
        PrdB  = 10*np.log10(Pr)
        PrpdB = 10*np.log10(Prp)

        return PrdB,PrpdB


if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
