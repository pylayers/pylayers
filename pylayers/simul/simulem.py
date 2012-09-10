#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""
     This module run the Ray Tracing simulation, it is a connection with PulsRay
"""
import doctest
import os
import getopt
import sys
import time
import ConfigParser
import pdb
import cPickle
import numpy as np
import scipy as sp
import scipy.io as spio
import matplotlib.pylab as plt
import struct
from pylayers.util import easygui
from pylayers.antprop.slab import Slab, SlabDB, Mat, MatDB
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
import pylayers.antprop.waveform as wvf
import pylayers.signal.bsignal as bs
from pylayers.simul.radionode import RadioNode
# Handle Layout
from pylayers.gis.layout import Layout
# Handle Rays
from pylayers.antprop.rays import GrRay3D, GrRayTud
# Handle VectChannel and ScalChannel
from pylayers.antprop import channel
#from   Channel import *
# Handle directory hierarchy
from pylayers.util.project import *
# Handle UWB measurements
from pylayers.measures import mesuwb as muwb


def spafile(_filename, point, sdir):
    """
        create a .spa file for  Ray Tracing

        Parameters
        ----------
            _filename
                shortname of the file .spa
            point
                3d coordinates string
            sdir
                save directory relative to $BASENAME

    """
    filespaTx = _filename
    filename = pyu.getlong(_filename, sdir)
    fspa = open(filename, "w")
    fspa.write("0\n")
    fspa.write("1\n")
    chaine = str(point).replace('[', '').replace(']', '')
    fspa.write(chaine + "\n")
    fspa.close()


class Palch(object):
    """ Launching parameters class 

    Methods
    -------

    info
    load
        load from Project launch directory
    save
        save to Project launch directory
    gui

    """
    #def __init__(self,angTx=8,ISBang=80,ethreshold=0.1,maxdeep=16,typalgo=0):
    def __init__(self, filename):
        self.filename = filename
        self.load()

    def info(self):
        """ display information
        """
        print "----------------------------------------------"
        print "            Launching Parameter               "
        print "----------------------------------------------"

        print "angTx      : Tx angular step ( degrees)     : ", self.angTx
        print "ISBang     : ISB angular sector ( degrees ) : ", self.ISBang
        print "ethreshold : Exploration Threshold (linear) : ", self.ethreshold
        print "maxdeep    : Tree deep max (integer value)  : ", self.maxdeep
        print "typalgo    : Type of algo (default 0)       : ", self.typalgo

    def info(self):
        for i, j in enumerate(self.__dict__.keys()):
            print j, ':', self.__dict__.values()[i]

    def load(self):
        filepalch = pyu.getlong(self.filename, 'launch')
        fi = open(filepalch)
        l = fi.read()
        u = l.split()
        self.angTx = eval(u[0])
        self.ISBang = eval(u[1])
        self.ethreshold = eval(u[2])
        self.maxdeep = eval(u[3])
        self.typalgo = eval(u[4])
        fi.close()

    def save(self):
        filepalch = pyu.getlong(self.filename, 'launch')
        fi = open(filepalch, 'w')
        fi.write(str(self.angTx) + '\n')
        fi.write(str(self.ISBang) + '\n')
        fi.write(str(self.ethreshold) + '\n')
        fi.write(str(self.maxdeep) + '\n')
        fi.write(str(self.typalgo) + '\n')
        fi.close()

    def gui(self):
        """
        Get the Launching parameter .palch
        """
        palchgui = multenterbox('', 'Launching Parameter',
                                ('Tx angular step (degrees)',
                                 'ISB angular sector (degrees)',
                                 'Exploration threshold (linear)',
                                 'Tree deep max (integer value)',
                                 'type of algo (default 0)'),
                                (self.angTx, self.ISBang, self.ethreshold,
                                 self.maxdeep, self.typalgo))

        if palchgui is not None:
            self.angTx = eval(palchgui[0])
            self.ISBang = eval(palchgui[1])
            self.ethreshold = eval(palchgui[2])
            self.maxdeep = eval(palchgui[3])
            self.typalgo = eval(palchgui[4])
            self.save()


class Patra(object):
    """
    Tracing parameters class
    """
    #def __init__(self,maxdeep=2,distdiff=2.0,var2D3D=1):
    def __init__(self, filename):
        self.filename = filename
        self.load()

    def info(self):
        print "----------------------------------------------"
        print "            Tracing  Parameter                "
        print "----------------------------------------------"
        print "Max deep     : ", self.maxdeep
        print "distdiff     : ", self.distdiff
        print "var2D3D  0=2D 1=3D    : ", self.var2D3D
        for i, j in enumerate(self.__dict__.keys()):
            print j, ':', self.__dict__.values()[i]

    def load(self):
        filepatra = pyu.getlong(self.filename, 'trace')
        fi = open(filepatra)
        l = fi.read()
        u = l.split()
        self.maxdeep = eval(u[0])
        self.distdiff = eval(u[1])
        self.var2D3D = eval(u[2])

    def save(self):
        """ save
        """
        filepatra = pyu.getlong(self.filename, 'trace')
        fi = open(filepatra, 'w')
        fi.write(str(self.maxdeep) + '\n')
        fi.write(str(self.distdiff) + '\n')
        fi.write(str(self.var2D3D) + '\n')
        fi.close()

    def gui(self):
        """ get the Launching parameter .palch
        """
        patragui = multenterbox('', 'Launching Parameter',
                                ('Max Deep ',
                                 'DistDiff',
                                 '2D3D '),
                                (self.maxdeep, self.distdiff, self.var2D3D))
        if patragui is not None:
            self.maxdeep = eval(patragui[0])
            self.distdiff = eval(patragui[1])
            self.var2D3D = eval(patragui[2])
            self.save()


class Pafreq(object):
    """ frequency setting
    """
    #def __init__(self,maxdeep=2,distdiff=2.0,var2D3D=1):
    def __init__(self,filename):
        self.filename = filename
        self.load()

    def info(self):
        """
           display frequency base information
        """
        print "----------------------------------------------"
        print "    Propagation Channel Frequency Ramp        "
        print "----------------------------------------------"
        print "fGHz min : ", self.fghzmin
        print "fGHz max : ", self.fghzmax
        print "Nf : ", self.Nf
        for i, j in enumerate(self.__dict__.keys()):
            print j, ':', self.__dict__.values()[i]

    def load(self):
        filefreq = pyu.getlong(self.filename, 'tud')
        fi = open(filefreq)
        l = fi.read()
        u = l.split()
        self.fghzmin = eval(u[0])
        self.fghzmax = eval(u[1])
        self.Nf = eval(u[2])

    def save(self):
        filefreq = pyu.getlong(self.filename, 'tud')
        fi = open(filefreq, 'w')
        fi.write(str(self.fghzmin) + ' ')
        fi.write(str(self.fghzmax) + ' ')
        fi.write(str(self.Nfp) + '\n')
        fi.close()

    def gui(self):
        """
        Get the Launching parameter .palch
        """
        pafreqgui = multenterbox('', 'Propagation Channel frequency ',
                                 ('fp_min (GHz) ',
                                  'fp_max (GHz) ',
                                  'Nfp  '),
                                 (self.fp_minGHz, self.fp_maxGHz, self.Nfp))
        if pafreqgui is not None:
            self.fp_minGHz = eval(pafreqgui[0])
            self.fp_maxGHz = eval(pafreqgui[1])
            self.Nfp = eval(pafreqgui[2])
            self.save()


class Patud(object):
    """ tratotud parameters
    """
    def __init__(self, purc=100, num=-1, nrmax=500):
        self.purc = purc
        self.num = num
        self.nrmax = nrmax

    def info(self):
        """ info 

        Examples
        --------
        >>> from pylayers.simul.simulem import *
        >>> p=Patud()
        >>> p.info()

        """
        print "----------------------------------------------"
        print "            tratotud  parameters              "
        print "----------------------------------------------"
        print "num (-1 all rays) : ", self.num
        print "nrmax : ", self.nrmax
        print "purc : ", self.purc
        for i, j in enumerate(self.__dict__.keys()):
            print j, ':', self.__dict__.values()[i]

    def gui(self):
        """ gui for tratotud parameters
        """
        tudgui = multenterbox('', 'Launching Parameter',
                              ('num',
                               'nrmax',
                               'purc'),
                              (self.num, self.nrmax, self.purc))
        self.num = eval(tudgui[0])
        self.nrmax = eval(tudgui[1])
        self.purc = eval(tudgui[2])


class Launch(object):
    """
    Class to receive the data from .lch files

    Attributes
    ----------

    Tx        array 1x3
    Ray_exist
    nstr
    deep
    x
    y
    node_phii
    node_phir
    edge_length
    edge_type
    tail
    head


    Methods
    ------
    load      : load a .lch file
    show      : view a .lch file
    info      : info about .lch
    launching :

    """
    def info(self):
        """ get __dict__ info
        """
        print len(self.x)
        print len(self.tail)
        for i, j in enumerate(self.__dict__.keys()):
            print j, ':', self.__dict__.values()[i]

    def choose(self):
        """ Choose a Launching  file in launchdir

        """
        import tkFileDialog
        FD = tkFileDialog
        filelch = FD.askopenfilename(filetypes=[("Fichiers Launching ", "*.lch"),
                                                ("All", "*")],
                                     title="Please choose a Launching file",
                                     initialdir=lchdir)
        _filelch = pyu.getshort(filelch)
        self.load(_filelch)

    def load(self, filelch):
        """
            load a .lch file
        """

        filelch = pyu.getlong(filelch, 'launch')
        fd = open(filelch, "rb")
        data = fd.read()
        fd.close()

        start = 0
        stop = start + 1024
        dt = data[start:stop]
        filestr = dt.replace("\x00", "")
        self.filestr = pyu.getshort(filestr)

        start = stop
        stop = start + 1024
        dt = data[start:stop]
        fileslab = dt.replace("\x00", "")
        self.fileslab = pyu.getshort(fileslab)

        start = stop
        stop = start + 1024
        dt = data[start:stop]
        filepalch = dt.replace("\x00", "")
        self.filepalch = pyu.getshort(filepalch)

        start = stop
        stop = start + 1024
        dt = data[start:stop]
        filespa = dt.replace("\x00", "")
        self.filespa = pyu.getshort(filespa)

        self.Tx = np.array([0.0, 0.0, 0.0])
        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Tx[0] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Tx[1] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Tx[2] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 4
        dt = data[start:stop]
        self.Ray_exist = stru.unpack('i', dt)[0]

        start = stop
        stop = start + 4
        dt = data[start:stop]
        self.node_num = stru.unpack('i', dt)[0]
        node_num = self.node_num

        self.tail = zeros(node_num - 1, dtype='int')
        self.head = zeros(node_num - 1, dtype='int')
        self.nstr = zeros(node_num, dtype='int')
        self.deep = zeros(node_num, dtype='int')
        self.x = zeros(node_num, dtype='float')
        self.y = zeros(node_num, dtype='float')
        self.node_phii = zeros(node_num, dtype='float')
        self.node_phid = zeros(node_num, dtype='float')
        self.edge_length = zeros(node_num - 1, dtype='float')
        self.edge_type = zeros(node_num - 1, dtype='float')

        for k in range(node_num - 1):
            start = stop
            stop = start + 4
            dt = data[start:stop]
            self.tail[k] = stru.unpack('i', dt)[0]

        for k in range(node_num - 1):
            start = stop
            stop = start + 4
            dt = data[start:stop]
            self.head[k] = stru.unpack('i', dt)[0]

        for k in range(node_num):
            start = stop
            stop = start + 4
            dt = data[start:stop]
            self.nstr[k] = stru.unpack('i', dt)[0]

        for k in range(node_num):
            start = stop
            stop = start + 4
            dt = data[start:stop]
            self.deep[k] = stru.unpack('i', dt)[0]

        for k in range(node_num):
            start = stop
            stop = start + 8
            dt = data[start:stop]
            self.x[k] = stru.unpack('d', dt)[0]

        for k in range(node_num):
            start = stop
            stop = start + 8
            dt = data[start:stop]
            self.y[k] = stru.unpack('d', dt)[0]

        for k in range(node_num):
            start = stop
            stop = start + 8
            dt = data[start:stop]
            self.node_phii[k] = stru.unpack('d', dt)[0]

        for k in range(node_num):
            start = stop
            stop = start + 8
            dt = data[start:stop]
            self.node_phid[k] = stru.unpack('d', dt)[0]

        for k in range(node_num - 1):
            start = stop
            stop = start + 8
            dt = data[start:stop]
            self.edge_length[k] = stru.unpack('d', dt)[0]

        for k in range(node_num - 1):
            start = stop
            stop = start + 4
            dt = data[start:stop]
            self.edge_type[k] = stru.unpack('i', dt)[0]

    def show(self, L, deepmax=1):
        """
        Launching display

        Parameters
        ----------

        L        : Layout
        deepmax  : display until deepmax (def=1)

        """
        sl = SlabDB()
        sl.mat = MatDB()
        sl.mat.load(self.fileslab.replace('.slab', '.mat'))
        sl.load(self.fileslab)
        #L      = Layout()
        #L.sl   = sl
        #G      = Graph(sl,filename=self.filestr)
        #fig    = figure(facecolor='white')
        #sp     = fig.add_subplot(111)
        L.showGs()
        #indoor.display['Visu']=False
        #indoor.display['alpha']=1.0
        #indoor.display['Node']=True
        #indoor.display['Edge']=True
        #indoor.show(fig,sp)

        Nseg = self.node_num - 1
        ita = self.tail - 1
        ihe = self.head - 1
        sdeep = self.deep[1::]
        usdeep = np.unique(sdeep)
        Mdeep = max(usdeep)
        plt.axis('scaled')
        plt.axis('off')
        for k in usdeep:
            if k <= deepmax:
                u = np.nonzero(sdeep == k)
                itak = ita[u[0]]
                ihek = ihe[u[0]]
                pt = np.vstack((self.x[itak], self.y[itak]))
                ph = np.vstack((self.x[ihek], self.y[ihek]))
                geu.displot(pt, ph, str(k / (1.0 * Mdeep)))

        """
        pz   =  empty((2,))
        pn   = zeros((2,))
        for i in range(Nseg):
            pz = vstack((pz,pt[:,i],ph[:,i],pn))
            m1   = np.array([0,0,1])
                mask = np.kron(ones((2,Nseg)),m1)
                pzz  = pz[1:,:].T
                vertices = np.ma.masked_np.array(pzz,mask)
                plot(vertices[0,:],vertices[1,:],color='black')
                show()


        for i in range(self.node_num-1):
            ita = self.tail[i]-1
            ihe = self.head[i]-1
            xt = self.x[ita]
            yt = self.y[ita]
            xh = self.x[ihe]
            yh = self.y[ihe]
            plot([xt,xh],[yt,yh],color='black',linewidth=1)

        show()
        """


class Simul(object):
    """ Simulation Class 

    Methods
    -------
        gui()
            graphical user interface
        choose()
            choose a simulation file in simuldir
        deepsave
            save simulation parameter
        info()
            info about the simulation
        showray(itx,irx,iray)
            show a single ray
        help()
            help on using Simulation object
        structure2()
            get .str2 file (ASCII description structure)
        freq()
            return the frequency base
        getlaunch(k)
            return the launching tree kth Transmittter
        save()
            save Simulation file
        layout
            load a Layout file 
        load()
            load Simulation
        show3()
            geomview vizualization
        show3l(itx,irx)
            geomview vizualization of link itx irx
        show()
            2D visualization of simulation with furniture
        launching()
            ray launching
        tracing(k)
            ray tracing
        tratotud(k,l)
            convert ray for tud
        field(l)
            evaluate field
        run(runlevel)
            run simulation
        CIR(itx,irx,store_level=0,alpha=1.0)
            Channel Impulse Response calculation


    Attributes
    ----------
        fileconf
        filestr
        filemat
        fileslab
        filepalch
        filepatra
        filefreq
        filespaTx
        filespaRx
        fileantTx
        fileantRx


    filefield
    filetauk
    filetang
    filerang

    tx
    rx

    palch
    patra
    patud
    progress

    indoor
    mat
    sl

    Notes
    ------

    This class group together all the parametrization files of a simulation

    Directory list file :

        fileconf

    Constitutive parameters file :

        fileslab
        filemat

    Ray Tracing parameters files :

        filepalch
        filetra

    Frequency list file :

        filefreq

    Position of Tx and Rx files :

        filespaTx
        filespaRx

    .. todo::
        migrer vers importations propres
    """
    def __init__(self, _filesimul='SimulDefault.ini'):

        self.filesimul = _filesimul
        self.filemat = []
        self.fileslab = []
        self.filespaTx = []
        self.filespaRx = []
        self.fileantTx = []
        self.fileantRx = []
        self.filepatra = []
        self.filepalch = []

        self.patud=Patud()

        self.freq = np.linspace(2.4,2.4,1,endpoint=True)

        self.progress = -1  # simulation not loaded
        self.filelch = []
        self.filestr = []
        self.filetra = []
        self.filetud = []
        self.filetang = []
        self.filerang = []
        self.filetauk = []
        self.filefield = []

        self.claunching = []
        self.ctracing = []
        self.ctratotud = []
        self.fileconf = []
        self.cfield = []
        #
        self.dout = {}
        self.dlch = {}
        self.dtra = {}
        self.dtud = {}
        self.dtang = {}
        self.drang = {}
        self.dtauk = {}
        self.dfield = {}
        self.dcir = {}

        self.freq = np.linspace(2, 11, 181)
        self.config = ConfigParser.ConfigParser()
        try:
            self.load()
        except:
            print('simulation file does not exist')

    def gui(self):
        """ gui to modify the simulation file
        """
        simulgui = multenterbox('', 'Simulation file',
                                ('filesimul',
                                 'filestr',
                                 'filefreq',
                                 'filespaTx',
                                 'filespaRx',
                                 'fileantTx'
                                 'fileantRx'),
                                (self.filesimul,
                                 self.filestr,
                                 self.filefreq,
                                 self.filespaTx,
                                 self.filespaRx,
                                 self.fileantTx,
                                 self.fileantRx))
        if simulgui is not None:
            self.filesimul = simulgui[0]
            self.filestr = simulgui[1]
            self.filefreq = simulgui[6]
            self.filespaTx = simulgui[7]
            self.filespaRx = simulgui[8]
            self.fileantTx = simulgui[9]
            self.fileantRx = simulgui[10]

    def clean(self, level=1):
        """ clean

       Notes
       -----
       obsolete

        Parameters
        ----------
            level  = 1


        """
        if level > 0:
            for itx in range(self.tx.N):
                filename = pyu.getlong(self.filelch[itx], 'launch')
                print filename
        if level > 1:
            for itx in range(self.tx.N):
                for irx in range(self.rx.N):
                    filename = pyu.getlong(self.filetra[itx][irx], 'trace')
                    print filename
        if level > 2:
            for itx in range(self.tx.N):
                for irx in range(self.rx.N):
                    filename = pyu.getlong(self.filetud[itx][irx], 'tud')
                    print filename
                    filename = pyu.getlong(self.filetauk[itx][irx], 'tud')
                    print filename
        if level > 3:
            for itx in range(self.tx.N):
                for irx in range(self.rx.N):
                    filename = pyu.getlong(self.filetang[itx][irx], 'tud')
                    print filename
                    filename = pyu.getlong(self.filerang[itx][irx], 'tud')
                    print filename
                    filename = pyu.getlong(self.filefield[itx][irx], 'tud')
                    print filename

    def deepsave(self):
        """ deepsave
        """
        sigui = multenterbox('', 'Simulation file',
                             ('filesimul', 'filestr', 'fileantTx',
                                 'fileantRx'),
                             (self.filesimul, self.filestr, self.fileantTx, self.fileantRx))
        self.filesimul = sigui[0]
        self.filestr = sigui[1]
        self.fileantTx = sigui[2]
        self.fileantRx = sigui[3]

        # filepalch
        fipalch = self.filesimul.replace('.simul', '.palch')
        self.palch.filename = fipalch
        self.palch.save()
        self.filepalch = fipalch

        # filepatra
        fipatra = self.filesimul.replace('.simul', '.patra')
        self.patra.filename = fipatra
        self.patra.save()
        self.filepatra = fipatra

        # filefreq
        fipafreq = self.filesimul.replace('.simul', '.freq')
        self.pafreq.filefreq = fipafreq
        self.pafreq.save()
        self.filefreq = fipafreq

        # filemat
        fimat = self.filesimul.replace('.simul', '.mat')
        self.mat.save(fimat)
        self.filemat = fimat

        # fileslab
        fislab = self.filesimul.replace('.simul', '.slab')
        self.sl.save(fislab)
        self.fileslab = fislab

        # filespaTx
        fispaTx = self.filesimul.replace('.simul', 'Tx.spa')
        self.tx.filespa = fispaTx
        self.tx.save()
        self.filespaTx = fispaTx

        # filespaRx
        fispaRx = self.filesimul.replace('.simul', 'Rx.spa')
        self.rx.filespa = fispaRx
        self.rx.save()
        self.filespaRx = fispaRx

        #
        # save simulation file
        #
        self.save()

    def save(self):
        """ save simulation file

        """
        filesimul = pyu.getlong(self.filesimul,"ini")
        fd = open(filesimul, "w")
        config = ConfigParser.ConfigParser()

        #
        # files section
        #
        config.add_section("files")
        config.set("files", "conf", self.fileconf)
        config.set("files", "struc", self.filestruc)
        config.set("files", "slab", self.fileslab)
        config.set("files", "mat", self.filemat)
        config.set("files", "tx", self.filespaTx)
        config.set("files", "rx", self.filespaRx)
        config.set("files", "txant", self.fileantTx)
        config.set("files", "rxant", self.fileantRx)
        config.set("files", "patra", self.filepatra)
        config.set("files", "palch", self.filepalch)

        #
        # tud section
        #

        config.add_section("tud")

        config.set("tud", "purc", self.patud.purc)
        config.set("tud", "nrmax", self.patud.nrmax)
        config.set("tud", "num", self.patud.num)

        #
        # frequency section
        #

        config.add_section("frequency")

        config.set("frequency", "fghzmin", self.freq[0])
        config.set("frequency", "fghzmax", self.freq[-1])
        config.set("frequency", "Nf", len(self.freq))

        #
        # output section 
        #

        #filelch exists
        if self.progress > 0:
            config.add_section("output")
            for k in range(len(self.filelch)):
                _fileout = "out" + "???"
                filename = self.filelch[k]
                fi.write(filename + '\n')
                config.set("launch", str(k + 1), filename)

        # filetra exists
        for k in range(len(self.filelch)):
            if self.progress > 1:
                config.add_section("trace")
                for l in arange(len(self.filetra[k])):
                    filename = self.filetra[k][l]
                    fi.write(filename + '\n')
                    config.set("trace", "rx" + str(l + 1), filename)

        # .tang exists
        # .rang exists
        # .tud exists
            if self.progress > 2:
                config.add_section("tud")
                config.add_section("tang")
                config.add_section("rang")

                for l in arange(len(self.filetud[k])):
                    ftud = self.filetud[k][l]
                    config.set("tud", "rx" + str(l + 1), ftud)

                for l in arange(len(self.filetang[k])):
                    ftang = self.filetang[k][l]
                    config.set("tang", "rx" + str(l + 1), ftang)

                for l in arange(len(self.filerang[k])):
                    frang = self.filerang[k][l]
                    config.set("rang", "rx" + str(l + 1), frang)

        # .field exist
        # .tauk exist
            if self.progress > 3:
                config.add_section("tauk")
                config.add_section("field")
                for l in arange(len(self.filetud[k])):
                    ftauk = self.filetud[k][l]
                    config.set("tauk", "rx" + str(l + 1), ftauk)

                for l in arange(len(self.filefield[k])):
                    ffield = self.filefield[k][l]
                    config.set("field", "rx" + str(l + 1), ffield)

        config.write(fd)
        fd.close()

    def choose(self):
        """
            Choose a simulation file in simuldir

        """
        import tkFileDialog as FD
        fichsimul = FD.askopenfilename(filetypes=[("Fichiers simul ", 
                                                   "*.simul"),
                                                  ("All", "*")],
                                       title="Please choose a simulation file", 
                                       initialdir=simuldir)
        self.filesimul = pyu.getshort(fichsimul)
        self.load()

    def load(self, _filesimul):
        """ load a simulation configuration file

         each transmiter simulation results in the creation of an .ini file
         with the following sections
         related to the results obtained for different receivers

        Parameters
        ----------
        _filesimul   : file in the simul directory of the Project

        """
        self.filesimul = _filesimul
        filesimul = pyu.getlong(self.filesimul, "ini")

        self.config.read(filesimul)

        sections = self.config.sections()

        filetx = self.config.get("files", "tx")
        filerx = self.config.get("files", "rx")

        fileanttx = self.config.get("files", "txant")
        fileantrx = self.config.get("files", "rxant")

        self.tx = RadioNode(filetx, fileanttx, 1)
        self.rx = RadioNode(filerx, fileantrx, 2)
#
# Launching and Tracing parameters
#
#        self.palch     = Palch(self.config.get("files","palch"))
        self.patra     = Patra(self.config.get("files","patra"))
        #_outfilename = "out"+self.ntx+".ini"
        #self.outfilename = pyu.getlong(_outfilename,"simul")
#  Load Simulation Mat File
#
        self.filemat = self.config.get("files", "mat")
        self.mat = MatDB()
        self.mat.load(self.filemat)
#
#  Load Simulation Slab File
#
        self.fileslab= self.config.get("files", "slab")
        self.sl = SlabDB()
        self.sl.mat = self.mat
        self.sl.load(self.fileslab)
#
# Load layout from .str or .str2 file
#
        self.filestr = self.config.get("files", "struc")
        self.L = Layout(self.filemat,self.fileslab)
        self.L.load(self.filestr)
#
# Frequency base
#
        if "frequency" in sections:
            self.freq = np.linspace(float(self.config.get("frequency", "fghzmin")),
                                    float(self.config.get("frequency", "fghzmax")),
                                    int(self.config.get("frequency", "nf")),
                                    endpoint=True)

            # update .freq file in tud directory 

            self.filefreq = pyu.getlong("def.freq", "tud")
            fd = open(self.filefreq, "w")
            chaine = self.config.get("frequency", "fghzmin") + ' ' + \
                self.config.get("frequency", "fghzmax") + ' ' + \
                self.config.get("frequency", "nf")
            fd.write(chaine)
            fd.close
#
# Simulation Progress
#
        self.output = {}
        if "output" in sections:
            for itx in self.config.options("output"):
                _filename = self.config.get("output", itx)
                filename = pyu.getlong(_filename, "simul")
                self.dout[int(itx)] = filename
                output = ConfigParser.ConfigParser()
                output.read(filename)
                secout = output.sections()
                self.dtra[int(itx)] = {}
                self.dtud[int(itx)] = {}
                self.dtang[int(itx)] = {}
                self.drang[int(itx)] = {}
                self.dtauk[int(itx)] = {}
                self.dfield[int(itx)] = {}
                self.dcir[int(itx)] = {}
                if "launch" in secout:
                    self.progress = 1
                    keys_launch = output.options("launch")
                    for kl in keys_launch:
                        self.dlch[int(kl)] = output.get("launch", kl)
                if "trace" in secout:
                    self.progress = 2
                    keys_tra = output.options("trace")
                    for kt in keys_tra:
                        self.dtra[int(itx)][int(kt)] = output.get("trace", kt)

                if "tang" in secout:
                    self.progress = 3
                    keys_tang = output.options("tang")
                    for kt in keys_tang:
                        self.dtang[int(itx)][int(kt)] = output.get("tang", kt)
                        self.drang[int(itx)][int(kt)] = output.get("rang", kt)
                        self.dtud[int(itx)][int(kt)] = output.get("tud", kt)
                if "field" in secout:
                    self.progress = 4
                    keys_field = output.options("field")
                    for kt in keys_field:
                        self.dfield[int(itx)][int(kt)] = output.get(
                            "field", kt)
                        self.dtauk[int(itx)][int(kt)] = output.get("tauk", kt)
                if "cir" in secout:
                    self.progress = 5
                    keys_cir = output.options("cir")
                    for kt in keys_cir:
                        self.dcir[int(itx)][int(kt)] = output.get("cir", kt)

                self.output[int(itx)] = output
        #
        # Waveform section
        #
        self.wav = wvf.Waveform()
        self.wav.read(self.config)

    def layout(self, _filestruc, _filemat='simul8.mat', _fileslab='simul8.slab'):
        """ load a layout in the simulation oject

        Parameters
        ----------
        _filestruc : string
            short file name of the Layout object
        _filemat   : string
            short file name of the Mat object  (default simul8.mat)
        _fileslab  : string
            short file name of the Slab object (default simul8.slab)

        Examples
        --------

        >>> from pylayers.simul.simulem import *
        >>> S = Simul()
        >>> S.layout('sircut.str')

        """

        self.filestr = _filestruc
        self.filemat = _filemat
        self.fileslab = _fileslab

        self.L = Layout(_filemat,_fileslab)
        self.L.load(_filestruc)

    def help(self):
        """
            use >>> instead
        """
        print "S.help()  : display help"
        print "S.info()  : display simulation info"
        print "S.palch.gui() : gui launching parameters "
        print "S.patra.gui() : gui tracing parameters"
        print "S.tx.gpoint() : gui Tx point "
        print "S.rx.gpoint() : gui Rx point"
        print "S.tx.position = hstack((Tx[1,:],1.2)).reshape(3,1)"
        print "S.rx.position = hstack((Rx[40,:],1.2)).reshape(3,1)"
        print "S.run(runlevel) : run launching tracing tratotud evalfield "
        print "S.launching() : exec ray launching "
        print "L = S.getlaunch(0) : get ray launching 0"
        print "S.tracing(0)  : exec ray tracing with tx 0"
        print "gr=GrRay3D() "
        print "gr.load(S.filetra[0][0],L.Gs)"
        print "gr.ray3d[0].show3() "
        print "S.tratotud(0) : exec tratotud with tx 0"
        print "S.field() : exec evalfield"
        print "VC = VectChannel(S,0,0,False)"
        print "SC = VC.vec2scalA(S.tx.A,S.rx.A)"

    def show(self, itx=[-1], irx=[-1], furniture=True, s=8, c='b', traj=False, num=False):
        """ show simulation

            Parameters
            -----------
            itx        : list of tx indexes
            irx        : list of rx indexes
            furniture  : boolean for METALIC furniture display
            s          : scale fir scatter plot  (default 8)
            c          : color for scatter plot  (default 'b')
            traj       : boolean  (def False)
            num        : boolean  (def False)


            Examples
            --------
            >>> import matplotlib.pyplot as plt
            >>> from pylayers.simul.simulem import *
            >>> S = Simul()
            >>> S.load('where2.ini')
            >>> S.L.loadfur('FurSiradel.ini')
            >>> S.show()
            >>> plt.show()



        """
        fig = plt.gcf()
        ax = fig.gca()
        #self.L.display['scaled']=False
        self.L.showGs(ax, show=False)
        #
        if furniture:
            if 'lfur' in self.L.__dict__:
                for fur in self.L.lfur:
                    if fur.Matname == 'METAL':
                        fur.show(fig, ax)
            else:
                print "Warning : no furniture file loaded"

        if irx[0] == -1:
            ax.scatter(self.rx.position[:, 0],
                       self.rx.position[:, 1], c='b', s=s, alpha=0.5)
            #ax.scatter(self.rx.position[0,0],self.rx.position[0,1],c='k',s=s,linewidth=0)
            #ax.scatter(self.rx.position[1,0],self.rx.position[1,1],c='b',s=s,linewidth=0)
            #ax.scatter(self.rx.position[2,0],self.rx.position[2,1],c='g',s=s,linewidth=0)
            #ax.scatter(self.rx.position[3,0],self.rx.position[3,1],c='c',s=s,linewidth=0)
        else:
            for k in irx:
                #ax.plot(self.rx.position[k-1,0],self.rx.position[k-1,1],'ob')
                ax.scatter(self.rx.position[k - 1, 0], self.rx.position[
                    k - 1, 1], c='b', s=s, alpha=0.5)
                if num:
                    ax.text(self.rx.position[k - 1, 0], self.rx.position[
                        k - 1, 1], str(k), color='blue')

        if itx[0] == -1:
            ax.scatter(self.tx.position[:, 0],
                       self.tx.position[:, 1], c='r', s=s)
            if num:
                for k in range(302):
                    ax.text(self.tx.position[k - 1, 0], self.tx.position[
                        k - 1, 1], str(k), color='black')
        else:
            #ax.plot(self.tx.position[itx-1,0],self.tx.position[itx-1,1],'k-')
            if traj:
                cpt = 1
            for k in itx:
                ax.scatter(self.tx.position[k - 1, 0], self.tx.position[
                    k - 1, 1], c=c, s=s, linewidth=0)
                if num:
                    if traj:
                        ax.text(self.tx.position[k - 1, 0], self.tx.position[k - 1, 1], str(cpt), color='black')
                        cpt = cpt + 1
                    else:
                        ax.text(self.tx.position[k - 1, 0], self.tx.position[k - 1, 1], str(k), color='black')
        #for k in range(self.tx.N):
        #    ax.text(self.tx.position[0,k],self.tx.position[1,k],str(k+1),color='black')
        #    ax.scatter(self.tx.position[0,:],self.tx.position[0,:],

        #for k in range(self.rx.N):
        #   ax.text(self.rx.position[0,k],self.rx.position[1,k],str(k),color='black')

    def PL(self, itx):
        """ plot Path Loss 

        itx 
        """
        td = []
        tEa = []
        tEo = []
        for irx in self.dcir[itx].keys():
            d = self.delay(itx, irx) * 0.3
            cira, ciro = self.loadcir(itx, irx)
            Ea = cira.energy()[0]
            Eo = ciro.energy()[0]
            td.append(d)
            tEa.append(Ea)
            tEo.append(Eo)

        plt.semilogx(td, 10 * np.log10(tEa), 'xr')
        plt.semilogx(td, 10 * np.log10(tEo), 'xb')
        plt.show()
        return td, tEa, tEo

    def loadcir(self, itx, irx):
        """
        Parameters
        ----------
        itx : Tx index
        irx : Rx index

        Returns
        -------

        cir(itx,irx)
        """
        _filecir = self.dcir[itx][irx] + '.mat'
        ext = str(itx)
        if len(ext) == 1:
            ext = '00' + ext
        if len(ext) == 2:
            ext = '0' + ext

        filecir = pyu.getlong(_filecir, 'cir/Tx' + ext)
        D = spio.loadmat(filecir)

        kxa = 'ta' + str(irx)
        kya = 'cira' + str(irx)

        kxo = 'to' + str(irx)
        kyo = 'ciro' + str(irx)

        cira = bs.TUsignal(D[kxa], D[kya][:, 0])
        ciro = bs.TUsignal(D[kxo], D[kyo][:, 0])

        return(cira, ciro)

    def pltcir(self, itx=1, irx=1, mode='linear', noise=False, color='b'):
        """ plot Channel Impulse Reresponse

        Parameters
        ----------
        itx : Tx index
        irx : Rx index
        mode : str
            {'linear' 'dB'}
            noise : boolean
        color : string
            default b

        >>> S = Simul()
        >>> S.load('where2.ini')
        >>> S.run2(1,[1])
        >>> S.pltcir(1,1,mode='linear',noise=False,color='k')

        """
        _filecir = self.dcir[itx][irx] + '.mat'
        filecir = pyu.getlong(_filecir, 'cir/Tx' + str('%0.3d' % itx))
        D = spio.loadmat(filecir)
        plt.subplot('211')
        self.show([itx], [irx])
        plt.subplot('212')
        kxa = 'ta' + str(irx)
        kya = 'cira' + str(irx)
        kxo = 'to' + str(irx)
        kyo = 'ciro' + str(irx)
        ta = D[kxa]
        to = D[kxo]

        Tobs = ta[-1] - ta[0]
        te = ta[1] - ta[0]
        if noise:
            na = bs.Noise(Tobs + te, 1. / te)
            naf = na.gating(4.493, 0.5)

        Tobs = to[-1] - to[0]
        te = to[1] - to[0]
        if noise:
            no = bs.Noise(Tobs + te, 1. / te)
            nof = no.gating(4.493, 0.5)

        cira = bs.TUsignal(ta, D[kya][:, 0])
        ciro = bs.TUsignal(to, D[kyo][:, 0])

        if mode == 'linear':
            #plt.plot(ta,naf.y,color='k',label='Noise')
            plt.plot(ta, D[kya], label='Rx ' + str(irx), color=color)
            '''if noise:
                naf.plot(col='k')
            cira.plot(col=color)'''
        else:
            '''if noise:
                naf.plotdB(col='k')
            cira.plotdB()'''
            plt.plot(ta, 20 * np.log10(
                abs(D[kya])), label='Rx ' + str(irx), color=color)
        plt.legend()
        plt.show()
        #plt.savefig('Tx'+str(itx),format=pdf,dpi=300)

    def scatter(self, itx, irx, values,
                cmap=plt.cm.gray,
                s=30,
                spl=221,
                title='',
                vaxis=((-30, 10, 2, 18)),
                vmin=0,
                vmax=1,
                colbool=False,
                cblabel='dB'):
        """
            Parameters
            ----------
            itx
            irx
            values
            cmap
            s
            spl
            title
            vaxis
            vmin
            vmax
            colbool
            clabel

        """
        fig = plt.gcf()
        ax = fig.add_subplot(spl)
        xtx = self.tx.position[itx, 0]
        ytx = self.tx.position[itx, 1]
        xrx = self.rx.position[irx, 0]
        yrx = self.rx.position[irx, 1]
        self.L.display['title'] = title
        self.L.showGs(ax)
        for furk in siradel.siradel_furniture.keys():
            fur = siradel.siradel_furniture[furk]
            if fur.Matname == 'METAL':
                fur.show(fig, ax)

        #self.show(furniture=True)
        plt.axis(vaxis)
        b1 = ax.scatter(xtx, ytx, s=s, c=values, cmap=cmap,
                        linewidths=0, vmin=vmin, vmax=vmax)
        ax.scatter(xrx, yrx, s=30, c='b', linewidths=0)
        if colbool:
            cb = colorbar(b1)
            cb.set_label(cblabel, fontsize=14)
        return(b1)

    def info(self, itx=[], irx=[]):
        """ display simulation information

         Parameters
         ----------
         itx : Tx index
         irx : Rx index

        """
        print self.filesimul
        print '------------------------------------------'
        try:
            print "Layout Info : \n", self.L.info()
        except:
            print "provide a Layout in the simulation : S.L "
            print ">>> S.layout(filename.str) "
            print "or "
            print ">>> S.layout(filename.str2 "
            print "or "
            print ">>> S.layout(filename.str,filemat,filemat) "
            print "default files exists for filemat and fileslab "

            return
        try:
            print "Tx Info :\n", self.tx.info()
        except:
            print "provide a tx in the simulation : S.tx "
            return
        try:
            print "Rx Info :\n", self.rx.info()
        except:
            print "provide a rx in the simulation : S.rx "
            return

        #    print "Tx : ",self.tx.points[itx]
            print "Rx : ", self.rx.points[itx]
            print "Delay (ns) :", self.delay(itx, irx)
#            print "Distance (m) :",0.3/self.delay(itx,irx)
#            print ""
#            if itx in self.dlch.keys():
#                print "-----"
#                print "Launching "
#                print "-----"
#                print " ",self.dlch[itx]
#            if irx in self.dtra[itx].keys():
#                print "-----"
#                print "Tracing "
#                print "-----"
#                print " ",self.dtra[itx][irx]
#                gr = GrRay3D()
#                gr.load(self.dtra[itx][irx],self.L)
#                gr.info()
#            if irx in self.dtud[itx].keys():
#                print "-----"
#                print "Tud parameters "
#                print "-----"
#                print " ",self.dtud[itx][irx]
#                print " ",self.dtang[itx][irx]
#                print " ",self.drang[itx][irx]
#                gt = GrRay3D.GrRayTud()
#                gt.load(self.dtud[itx][irx], self.dtang[itx][irx], self.drang[itx][irx],self.sl)
#            if irx in self.dtauk[itx].keys():
#                print self.dtauk[itx][irx]
#                print self.dfield[itx][irx]
#                VC = self.VC(itx,irx)
#            if irx in self.dcir[itx].keys():
#                print self.dcir[itx][irx]

    def info2(self):
        for i, j in enumerate(self.__dict__.keys()):
            print j, ':', self.__dict__.values()[i]

    def filtray(self, itx, irx, tau0, tau1, col='b'):
        """ filter rays 

        Parameters
        ----------
        itx : 
        irx :
        tau0 :
        tau1 :
        col :

        Display ray and nstr
        """
        gr = GrRay3D()
        gr.load(self.dtra[itx][irx], self.L)
        self.L.display['Thin'] = True
        self.L.display['Node'] = False
        self.L.display['NodeNum'] = False
        self.L.display['EdgeNum'] = False
        plt.axis('scaled')
        #self.L.show(fig,ax,nodelist,seglist)
        self.L.showGs()
        delays = gr.delay()
        rayset = np.nonzero((delays >= tau0) & (delays <= tau1))[0]
        fig = plt.gcf()
        ax = fig.get_axes()[0]
        gr.show(ax, rayset, col=col, node=False)
        plt.title('Tx' + str(itx) + '-Rx' + str(irx) + ' : ' + str(
            tau0) + ' < tau < ' + str(tau1))

    def showray(self, itx, irx, iray=np.array([]), fig=[], ax=[]):
        """ show layout and rays for a radio link

        Parameters
        ----------
        itx  : tx index
        irx  : rx index
        iray : list of rays to be displayed ndarray


        """

        #if ax==[]:
        #    fig = plt.figure()
        #    ax  = fig.add_subplot('111')

        gr = GrRay3D()
        gr.load(self.dtra[itx][irx], self.L)
        if len(iray == 1):
            ray = gr.ray3d[iray[0]]
            nstr = ray.nstr[1:-1]
            uneg = np.nonzero(nstr < 0)
            upos = np.nonzero((nstr > 0) & (nstr <= self.L.Ne))
            uceil = np.nonzero(nstr == self.L.Ne + 1)
            ufloor = np.nonzero(nstr == self.L.Ne + 2)
        #seglist  = nstr[upos[0]]-1
        #nodelist = -nstr[uneg[0]]-1
        #seglist2 = S.L.segpt(nodelist)
        self.L.display['Thin'] = True
        self.L.display['Node'] = False
        self.L.display['NodeNum'] = False
        self.L.display['EdgeNum'] = False
        #self.L.show(fig,ax)
        #print ray.nn
        #print len(ray.nstr)
        #print ray.nstr
        #print nodelist
        #print seglist
        #print seglist2
        #seglist = hstack((seglist,seglist2))
        self.L.display['Node'] = False
        self.L.display['Thin'] = False
        self.L.display['NodeNum'] = True
        self.L.display['EdgeNum'] = True
        plt.axis('scaled')
        #self.L.show(fig,ax,nodelist,seglist)
        fig, ax = self.L.showGs(show=False)
        gr.show(ax, iray, col='b', node=False)

        if len(iray) == 1:
            plt.title(str(nstr))
        else:
            plt.title('Tx' + str(itx) + '-Rx' + str(irx) +
                      ' ' + str(min(iray)) + ' ' + str(max(iray)))

    def show3l(self, itx, irx):
        """ geomview display of a specific link

        g = S.show3l(itx,irx)

        Parameters
        ----------
        itx
            transmitter index
        irx
            receiver index

        """
        filetra = self.dtra[itx][irx]
        gr = GrRay3D()
        gr.load(filetra, self.L)
        gr.show3()

        return(gr)

    def show3(self):
        """ geomview display of the simulation configuration
        """
        try:
            self.tx.save()
        except:
            print('tx set is no defined')
        try:
            self.rx.save()
        except:
            print('rx set is no defined')
        _filename = self.filesimul.replace('.ini','.off')
        filename = pyu.getlong(_filename, "geom")
        fo = open(filename, "w")
        fo.write("LIST\n")
        try:
            sttx = "{<"+self.tx.filegeom+"}\n"
        except:
            sttx = "\n"
        try:
            strx = "{<"+self.rx.filegeom+"}\n"
        except:
            strx = "\n"
        try:
            stst = "{<"+self.L.filegeom+"}\n"
        except:
            stst = "\n"
        fo.write(sttx)
        fo.write(strx)
        fo.write(stst)
        fo.write("{</usr/share/geomview/geom/xyz.vect}\n")
        fo.close()
        command = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
        os.system(command)

    def freq(self, GUI=False):
        """
            return the frequency base from the content of filefreq
        Parameters
        ----------
        GUI
            Boolean


        """
        filefreq = self.filefreq
        filefreq = pyu.getlong(filefreq, 'tud')
        fo = open(filefreq)
        l = fo.readline().split()
        fmin = eval(l[0])
        fmax = eval(l[1])
        Nf = eval(l[2])
        fo.close()
        if GUI == True:
            val = multenterbox('Enter frequency ramp', '',
                               ('start', 'stop', 'N'), 
                               (str(fmin), str(fmax), str(Nf)))
            fmin = eval(val[0])
            fmax = eval(val[1])
            Nf = eval(val[2])
            fo = open(filefreq, 'w')
            data = val[0] + ' ' + val[1] + ' ' + val[2] + '\n'
            fo.write(data)
            fo.close()

        freq = np.linspace(fmin, fmax, Nf, endpoint=True)
        return(freq)

    def getlaunch(self, k=1):
        """
         get the kth launch

        Parameters
        ----------
        k :  launching index
        """
        if str(k) in self.dlch.keys():
            filename = self.dlch[str(k)]
            L = Launch()
            L.load(filename)
            return(L)
        else:
            print("Tx not available")

        def gettra(self, itx, irx):
            """
             Gtra = S.gettra(itx,irx)
            """
            if (itx < self.tx.N) & (irx < self.rx.N):
                Gtra = GrRay3D()
                Gtra.load(self.filetra[itx][irx], self.indoor)
                return Gtra
            else:
                print "gettra warning : wrong tx or rx index"

        def gettud(self, itx, irx):
            """
             Gtud = S.gettud(itx,irx)
            """
            if (itx < self.tx.N) & (irx < self.rx.N):
                Gtud = GrRayTud()
                Gtud.load(self.filetud[itx][irx], self.sl)
                return Gtud
            else:
                print "gettud warning : wrong tx or rx index"
        #def gettud(self,k,l):
        #def getfield(self,k,l):

    def launching(self, itx=1):
        """
           start the launching program and get the results files

        Parameters
        ----------
        itx : int 
            transmiter index

        """
        chaine = "launching -str  " + self.config.get("files", "struc") + \
            " -slab " + self.config.get("files", "slab") + \
            " -palch " + self.config.get("files", "palch") + \
            " -spa " + self.filespaTx + \
            " -conf " + self.config.get("files", "conf")

        print chaine

        self.claunching.append(chaine)
        os.system(chaine)
        aux = os.popen(chaine, "r")
        recup = aux.read()
        aux.close()

        print recup

        self.recup = recup
        aux = recup.splitlines()
        len_aux = recup.count("\n")

        filelch = []
        for i in range(len_aux):
            if aux[i].find("filelchout") != -1:
                aux[i] = aux[i].replace('filelchout : ', '')
                fshort = pyu.getshort(aux[i])
                filelch.append(fshort)
                self.dlch[itx] = fshort

        self.filelch = filelch
        self.progress = 1
        self.nTx = len(self.filelch)

    def tracing(self, itx, irx):
        """
        tracing(ntx)  : exec tracing

        Parameters
        ----------
        itx  : tx index
        irx  : rx index

        """
        if (self.progress >= 1):
            #chaine = "tracing -lch " + self.filelch[l] + " -patra " + self.filepatra + " -spa " + self.filespaRx + " -conf "  + self.fileconf
            chaine = "tracing -lch " + self.dlch[itx] + \
                " -patra " + self.config.get("files", "patra") + \
                "  -spa " + self.filespaRx + \
                " -conf " + self.config.get("files", "conf")
            print chaine
            self.ctracing.append(chaine)
            aux = os.popen(chaine, "r")
            recup = aux.read()
            print recup
            aux.close()
            self.recup = recup
            aux = recup.splitlines()
            len_aux = recup.count("\n")
            filetra = []
            for i in range(len_aux):
                if aux[i].find("filetraout") != -1:
                    aux[i] = aux[i].replace('filetraout : ', '')
                    filetra.append(pyu.getshort(aux[i]))
                    self.dtra[itx][irx] = filetra[0]
            print filetra
            #self.filetra.insert(ntx,filetra)
            self.progress = 2
        else:
            print "No launching available"

    def tratotud2(self, itx, irx):
        """ convert tracing in .tud

        Parameters
        ----------

         itx : integer 
         irx : integer 

        """
        nrmin = self.config.get("tud", "nrmax")
        num = self.config.get("tud", "num")
        purc = self.config.get("tud", "purc")
        chaine = "tratotud -tra " + self.dtra[itx][irx] + \
            " -min " + nrmin + \
            " -purc " + purc + \
            " -num " + num +  \
            " -conf " + self.config.get("files", "conf")
        self.ctratotud.append(chaine)
        print chaine
        aux = os.popen(chaine, "r")
        os.system('echo $?')
        recup = aux.read()
        aux.close()
    #   self.recup=recup
        aux = recup.splitlines()
        print aux
        len_aux = recup.count("\n")
        for i in range(len_aux):
            if aux[i].find("filetudout") != -1:
                aux[i] = aux[i].replace('filetudout : ', '')
                filename = pyu.getshort(aux[i])
                print filename
                self.dtud[itx][irx] = filename
            elif aux[i].find("filetangout") != -1:
                aux[i] = aux[i].replace('filetangout : ', '')
                filename = pyu.getshort(aux[i])
                print filename
                self.dtang[itx][irx] = filename
            elif aux[i].find("filerangout") != -1:
                aux[i] = aux[i].replace('filerangout : ', '')
                filename = pyu.getshort(aux[i])
                print filename
                self.drang[itx][irx] = filename

    def tratotud(self, itx, irx):
        """
         tratotud(itx,irx)  convert tracing in .tud

         l : launching index

         loop on all file filetra[l][k]

         .. todo::
            check if the loop cannot be done in C

        """
        print " t2t ", itx, irx
        if (self.progress >= 2):
            nrmin = self.config.get("tud", "nrmax")
            num = self.config.get("tud", "num")
            purc = self.config.get("tud", "purc")
            filetud = []
            filetang = []
            filerang = []
            print nrmin
            print purc
            print num
            print self.dtra[itx][irx]
            chaine = "tratotud -tra " + self.dtra[itx][irx] + \
                " -min " + nrmin + \
                            " -purc " + purc + \
                " -num " + num +  \
                            " -conf " + self.config.get("files", "conf")
            self.ctratotud.append(chaine)
            print chaine
            aux = os.popen(chaine, "r")
            os.system('echo $?')
            recup = aux.read()
            aux.close()
        #   self.recup=recup
            aux = recup.splitlines()
            len_aux = recup.count("\n")
        #   self.aux=aux
        #   .. todo::
        #       remove the loop
            for i in range(len_aux):
                if aux[i].find("filetudout") != -1:
                    aux[i] = aux[i].replace('filetudout : ', '')
                    filename = pyu.getshort(aux[i])
                    #filetud.append(filename)
                    self.dtud[itx][irx] = filename
                elif aux[i].find("filetangout") != -1:
                    aux[i] = aux[i].replace('filetangout : ', '')
                    filename = pyu.getshort(aux[i])
                    #filetang.append(filename)
                    self.dtang[itx][irx] = filename
                elif aux[i].find("filerangout") != -1:
                    aux[i] = aux[i].replace('filerangout : ', '')
                    filename = pyu.getshort(aux[i])
                    #filerang.append(filename)
                    self.drang[itx][irx] = filename

            #self.filetud.insert(l,filetud)
            #self.filetang.insert(l,filetang)
            #self.filerang.insert(l,filerang)
            self.progress = 3
        else:
            print "No tracing available"

    def field2(self, itx, irx):
        """
        field calculation for Tx Rx

        Parameters
        ----------
            ntx
                 launching index
            nrx
                 tracing index
        """
        chaine = "evalfield -tud " + self.dtud[itx][irx] + \
                 " -slab " + self.config.get("files", "slab") + \
                 " -mat " + self.config.get("files", "mat") + \
                 " -freq " + pyu.getshort(self.filefreq) + \
                 " -conf " + self.config.get("files", "conf")

        self.cfield.append(chaine)
        print chaine
        os.system(chaine)
        aux = os.popen(chaine, "r")
        recup = aux.read()
        aux.close()
       # self.recup=recup
        aux = recup.splitlines()
        print aux
        len_aux = recup.count("\n")
        for i in range(len_aux):
            print 'dbg ": ', aux[i]
            if aux[i].find("filefieldout") != -1:
                aux[i] = aux[i].replace('filefieldout : ', '')
                filename = pyu.getshort(aux[i]).replace(' ', '')
                self.dfield[itx][irx] = filename
            elif aux[i].find("filetaukout") != -1:
                aux[i] = aux[i].replace('filetaukout : ', '')
                filename = pyu.getshort(aux[i]).replace(' ', '')
                self.dtauk[itx][irx] = filename

    def field(self, itx, irx):
        """
        field calculation for Tx Rx

        Parameters
        ----------
            ntx
                 launching index
            nrx
                 tracing index
        """
        if (self.progress >= 3):
            filefield = []
            filetauk = []

            chaine = "evalfield -tud " + self.dtud[itx][irx] + \
                     " -slab " + self.config.get("files", "slab") + \
                     " -mat " + self.config.get("files", "mat") + \
                     " -freq " + pyu.getshort(self.filefreq) + \
                     " -conf " + self.config.get("files", "conf")

            self.cfield.append(chaine)
            print chaine
            os.system(chaine)
            aux = os.popen(chaine, "r")
            recup = aux.read()
            aux.close()
           # self.recup=recup
            aux = recup.splitlines()
            len_aux = recup.count("\n")
            for i in range(len_aux):
                print 'dbg ": ', aux[i]
                if aux[i].find("filefieldout") != -1:
                    aux[i] = aux[i].replace('filefieldout : ', '')
                    filename = pyu.getshort(aux[i]).replace(' ', '')
                    self.dfield[itx][irx] = filename
                    #filefield=pyu.getshort(aux[i])
                    #filefield=filefield.replace(' ','')
                    #filefield.append(pyu.getshort(aux[i]).replace(' ',''))
                elif aux[i].find("filetaukout") != -1:
                    aux[i] = aux[i].replace('filetaukout : ', '')
                    filename = pyu.getshort(aux[i]).replace(' ', '')
                    self.dtauk[itx][irx] = filename
                    #filetauk=pyu.getshort(aux[i])
                    #filetauk=filetauk.replace(' ','')
                    #filetauk.append(pyu.getshort(aux[i]).replace(' ',''))

            #self.filefield.insert(l,filefield)
            #self.filetauk.insert(l,filetauk)

            self.progress = 4
        else:
            print "no .tud file available, try tratotud"

    def saveconfig(self):
        """
            save current config
            To be  remove .simul is obsolete
        """
        filesimul = pyu.getlong(self.filesimul, "simul")
        fileini = filesimul.replace(".simul", ".ini")
        f2 = open(fileini, "w")
        self.config.write(f2)
        f2.close()

    def run(self, runlevel=5):
        """
        run the simulation until a given runlevel

        runlevel 1  : launching
        runlevel 2  : tracing
        runlevel 3  : tratotud
        runlevel 4  : field
        runlevel 5  : cir

        """
        t0 = time.clock()
        if self.progress < 1:
            #
            # Analyse le numero de launching a lancer
            # et construit le fichier spa associé
            #
            ltx = self.tx.points.keys()
            etx = self.dlch.keys()
            for k in ltx:
                if str(k) not in etx:
                    self.ntx = str(k)
                    break
            _filename = 'spaTx.spa'
            self.filespaTx = _filename
            filename = pyu.getlong(_filename, "launch")
            fspa = open(filename, "w")
            fspa.write("0\n")
            fspa.write("1\n")
            chaine = str(self.tx.points[k]).replace('[', '').replace(']', '')
            fspa.write(chaine + "\n")
            fspa.close()
            print "---------------"
            print "Start Launching Tx : ", self.ntx
            print "---------------"
            self.launching(self.ntx)
            t1 = time.clock()
            print "elapsed ", t1 - t0
            if "output" not in  self.config.sections():
                self.config.add_section("output")
            _outfilename = "out" + self.ntx + ".ini"
            self.outfilename = pyu.getlong(_outfilename, "simul")
            self.config.set("output", self.ntx, _outfilename)
            self.output.add_section("launch")
            self.output.set("launch", self.ntx, self.dlch[self.ntx])
            fd = open(self.outfilename, "w")
            self.output.write(fd)
            fd.close()
            self.saveconfig()
            self.dtra[self.ntx] = {}
            self.dtud[self.ntx] = {}
            self.dtang[self.ntx] = {}
            self.drang[self.ntx] = {}
            self.dtauk[self.ntx] = {}
            self.dfield[self.ntx] = {}
            self.dcir[self.ntx] = {}

        if (self.progress < 2) & (runlevel > 1):
            #
            # Analyse les numero de rx a choisir
            # et construit le fichier spa associé
            #
            lrx = self.rx.points.keys()
            erx = self.dtra[self.ntx].keys()
            for k in lrx:
                print "k= ", k
                print "erx", erx
                if str(k) not in erx:
                    self.nrx = str(k)
                    break
            _filename = 'spaRx.spa'
            self.filespaRx = _filename
            filename = pyu.getlong(_filename, "trace")
            fspa = open(filename, "w")
            fspa.write("0\n")
            fspa.write("1\n")
            chaine = str(self.tx.points[k]).replace(
                '[', '').replace(']', '')
            fspa.write(chaine + "\n")
            fspa.close()
            print "     --------------------"
            print "     Start Tracing  Rx : ", self.nrx
            print "     --------------------"
            #self.filetra=[]
            #for l in range(len(self.filelch)):
            #    print "launching N° : ",l
            t2 = time.clock()
            self.tracing(self.ntx, self.nrx)

            if "trace" not in  self.output.sections():
                self.output.add_section("trace")
            self.output.set("trace", self.nrx, self.dtra[self.ntx][self.nrx])
            fd = open(self.outfilename, "w")
            self.output.write(fd)
            fd.close()

            #print "elapsed time",time.clock()-t2
            #gr = GrRay3D()
            #for nf in self.filetra[l]:
            #    gr.load(nf,self.L)
            #    gr.save(nf)

        if (self.progress < 3) & (runlevel > 2):
            print "---------------"
            print "Start Tratotud"
            print "---------------"
            #self.filetud=[]
            #self.filetang=[]
            #self.filerang=[]
            #for l in range(len(self.filelch)):
            print "launching N° : ", self.ntx
            t3 = time.clock()
            self.tratotud(self.ntx, self.nrx)
            if "tud" not in  self.output.sections():
                self.output.add_section("tud")
            if "tang" not in  self.output.sections():
                self.output.add_section("tang")
            if "rang" not in  self.output.sections():
                self.output.add_section("rang")
            self.output.set("tud", self.ntx, self.dtud[self.ntx][self.nrx])
            self.output.set("tang", self.ntx, self.dtang[self.ntx][self.nrx])
            self.output.set("rang", self.ntx, self.drang[self.ntx][self.nrx])
            fd = open(self.outfilename, "w")
            self.output.write(fd)
            fd.close()
            print "elapsed time", time.clock() - t3
            print "---------------"
            print "Stop Tratotud"
            print "---------------"
        if (self.progress < 4) & (runlevel > 3):
            print "---------------"
            print "Start Field"
            print "---------------"
            #self.filetauk=[]
            #self.filefield=[]
            #for l in range(len(self.filelch)):
            print "launching N° : ", self.ntx
            t3 = time.clock()
            self.field(self.ntx, self.nrx)
            if "field" not in  self.output.sections():
                self.output.add_section("field")
            if "tauk" not in  self.output.sections():
                self.output.add_section("tauk")
            self.output.set("field", self.ntx, self.dfield[self.ntx][self.nrx])
            self.output.set("tauk", self.ntx, self.dtauk[self.ntx][self.nrx])
            fd = open(self.outfilename, "w")
            self.output.write(fd)
            fd.close()
            print "elapsed time", time.clock() - t3
            print "---------------"
            print "Stop Field "
            print "---------------"
        if (self.progress < 5) & (runlevel > 4):
            if "waveform" in self.config.sections():
                par = self.config.items("waveform")
                self.wparam = {}
                for k in range(len(par)):
                    key = par[k][0]
                    val = par[k][1]
                    if key == "band":
                        self.wparam[key] = float(val)
                    if key == "fc":
                        self.wparam[key] = float(val)
                    if key == "fe":
                        self.wparam[key] = float(val)
                    if key == "thresh":
                        self.wparam[key] = float(val)
                    if key == "tw":
                        self.wparam[key] = float(val)
                    if key == "type":
                        self.wparam[key] = val

                self.wav = wvf.Waveform(self.wparam)
                Ntx = self.tx.N
                Nrx = self.rx.N
                alpha = np.sqrt(1. / 30.0)
                # S.CIR(arange(Ntx),arange(Nrx),self.wav,store_level=16+8+4+2+1,alpha=alpha)
                #stx      = S.dlch.keys()
                #srx      = S.dtra[stx[0]].keys()
                #S.CIR(stx,srx,store_level=16+8+4+2+1,alpha=alpha)
                self.CIR([self.ntx], [self.nrx], store_level=16 +
                         8 + 4 + 2 + 1, alpha=alpha)
                if "cir" not in  self.output.sections():
                    self.output.add_section("cir")
                #self.output.set("cir",self.ntx,self.dcir[self.ntx][self.nrx])
            else:
                raise("Error no waveform in the config file ")

    def run2(self, itx, srx=[], cirforce=True):
        """
            run the simulation for 1 tx and a set of rx

            index point start in 1

            Parameters
            ----------

            itx      : tx index
            srx      : list of rx index
            cirforce : boolean


            Example
            -------

            >>> from pylayers.simul.simulem import *
            >>> itx = 1
            >>> srx = [1,2,3]
            >>> S   = Simul()
            >>> S.load('where2.ini')
            >>> out = S.run2(itx,srx)


        """
        #t0 = time.clock()
        if srx == []:
            srx = self.rx.points.keys()

        if itx not in self.dlch.keys():
            # launching itx does not exist
            self.filespaTx = 'tx' + str(itx) + '.spa'
            point = self.tx.points[itx]
            spafile(self.filespaTx, point, "launch")
            print "---------------"
            print "Start Launching Tx : " + str(itx)
            print "---------------"
            self.launching(itx)
            if "output" not in  self.config.sections():
                self.config.add_section("output")
            self.dtra[itx] = {}
            self.dtud[itx] = {}
            self.dtang[itx] = {}
            self.drang[itx] = {}
            self.dtauk[itx] = {}
            self.dfield[itx] = {}
            self.dcir[itx] = {}
            self.output[itx] = ConfigParser.ConfigParser()
            self.output[itx].add_section("launch")
            self.output[itx].set("launch", str(itx), self.dlch[itx])
        else:
            print "launching exists"

        _outfilename = "out" + str(itx) + ".ini"
        outfilename = pyu.getlong(_outfilename, "simul")
        self.config.set("output", str(itx), _outfilename)
        #
        # save output files
        #
        fd = open(outfilename, "w")
        self.output[itx].write(fd)
        fd.close()
        self.saveconfig()
        #
        # Loop over a set of rx
        #
        for irx in srx:
            tracingabort = False
            if irx not in self.dtra[itx].keys():
                self.filespaRx = 'rx' + str(irx) + '.spa'
                point = self.rx.points[irx]
                spafile(self.filespaRx, point, "trace")
                print "--------------------"
                print "Start Tracing  Rx : " + str(irx)
                print "--------------------"
                self.tracing(itx, irx)
                if irx in self.dtra[itx].keys():
                    if "trace" not in  self.output[itx].sections():
                        self.output[itx].add_section("trace")
                    self.output[itx].set("trace", str(
                        irx), self.dtra[itx][irx])
                    fd = open(outfilename, "w")
                    self.output[itx].write(fd)
                    fd.close()
                else:
                    print 'Erreur dans tracing ', irx
                    tracingabort = True

            if not tracingabort:
                if irx not in self.dtud[itx].keys():
                    print "---------------"
                    print "Start Tratotud ", irx
                    print "---------------"
                    self.tratotud2(itx, irx)
                    if "tud" not in  self.output[itx].sections():
                        self.output[itx].add_section("tud")
                    if "tang" not in  self.output[itx].sections():
                        self.output[itx].add_section("tang")
                    if "rang" not in  self.output[itx].sections():
                        self.output[itx].add_section("rang")
                    if irx in self.dtud[itx].keys():
                        self.output[itx].set("tud", str(
                            irx), self.dtud[itx][irx])
                    if irx in self.dtang[itx].keys():
                        self.output[itx].set("tang", str(irx),
                                             self.dtang[itx][irx])
                    if irx in self.drang[itx].keys():
                        self.output[itx].set("rang", str(irx),
                                             self.drang[itx][irx])
                    gt = self.gt(itx, irx)
                    gt.save(self.dtud[itx][irx], 
                            self.dtang[itx][irx],
                            self.drang[itx][irx])
                    fd = open(outfilename, "w")
                    self.output[itx].write(fd)
                    fd.close()
                if irx not in self.dfield[itx].keys():
                    print "---------------"
                    print "Start Field  ", irx
                    print "---------------"
                    self.field2(itx, irx)
                    if "field" not in  self.output[itx].sections():
                        self.output[itx].add_section("field")
                    if "tauk" not in  self.output[itx].sections():
                        self.output[itx].add_section("tauk")
                    if irx in self.dfield[itx].keys():
                        self.output[itx].set("field", str(irx),
                                             self.dfield[itx][irx])
                    if irx in self.dtauk[itx].keys():
                        self.output[itx].set("tauk", str(irx),
                                             self.dtauk[itx][irx])
                    fd = open(outfilename, "w")
                    self.output[itx].write(fd)
                    fd.close()

                if ((irx not in self.dcir[itx].keys()) | (cirforce == True)):
                    print "---------------"
                    print "Start CIR      ", irx
                    print "---------------"
                    if "waveform" in self.config.sections():
                        par = self.config.items("waveform")
                        self.wparam = {}
                        for k in range(len(par)):
                            key = par[k][0]
                            val = par[k][1]
                            if key == "band":
                                self.wparam[key] = float(val)
                            if key == "fc":
                                self.wparam[key] = float(val)
                            if key == "fe":
                                self.wparam[key] = float(val)
                            if key == "thresh":
                                self.wparam[key] = float(val)
                            if key == "tw":
                                self.wparam[key] = float(val)
                            if key == "type":
                                self.wparam[key] = val

                        self.wav = wvf.Waveform(self.wparam)
                        alpha = np.sqrt(1. / 30.0)
                        self.CIR([itx], [irx], store_level=16 + 8 +
                                 4 + 2 + 1, alpha=alpha)
                        if "cir" not in  self.output[itx].sections():
                            self.output[itx].add_section("cir")
                        if irx in self.dcir[itx].keys():
                            self.output[itx].set("cir", str(irx), self.dcir[
                                itx][irx])
                        fd = open(outfilename, "w")
                        self.output[itx].write(fd)
                        fd.close()
                    else:
                        raise("Error no waveform in the config file ")
                        return(False)

        return(True)

    def gt(self, itx, irx):
        """ gtud
        """
        gt = GrRayTud()
        filetud = self.dtud[itx][irx]
        filetang = self.dtang[itx][irx]
        filerang = self.drang[itx][irx]
        gt.load(filetud, filetang, filerang, self.sl)
        return(gt)

    def delay(self, itx, irx):
        """
            calculate LOS link delay

            Parameters
            ----------

            itx
            irx

        """
        tx = self.tx.points[itx]
        rx = self.rx.points[irx]
        df = tx - rx
        dist = np.sqrt(np.dot(df, df))
        return(dist / 0.3)

    def gr(self, itx, irx):
        """
           return a  cluster or fays from link itx-irx
        """
        gr = GrRay3D()
        gr.load(self.dtra[itx][irx], self.L)
        return(gr)

    def VC(self, itx, irx):
        """
            return Vect Channel for link itx irx
        """

        VCl = channel.VectChannel(self, itx, irx, False)
        return(VCl)

    def CIR(self, itx, irx, store_level=0, alpha=1.0, ext='', rep='cir'):
        """
        Calculate a set of channel impulse responses

        Parameters
        ----------
            itx : transmitter index iterable set
            irx : receiver index iterable set
            wav : applied waveform (st,sf,sfg)
                this waveform includes the gamma factor
            store_level : binary mask
                             bit  0 : save CVC
                             bit  1 : save CVCO
                             bit  2 : save CSCA
                             bit  3 : save CIRo
                             bit  4 : save CIRa
            alpha : normalization factor

         Notes
         -----
         A factor ::math`\\sqrt{1}{30}` is necessary when applying the antenna

         Examples
         --------

        """

        if (store_level) & 1 == 1:
            self.CVC = []
        if (store_level) & 2 == 2:
            self.CSCO = []
        if (store_level) & 4 == 4:
            self.CSCA = []
        if (store_level) & 8 == 8:
            self.CIRo = []
        if (store_level) & 16 == 16:
            self.CIRa = []
            racine = 'cir-'
        for l in itx:
            CVC = []
            CSCO = []
            CSCA = []
            CIRo = []
            CIRa = []
            for k in irx:
                D = {}
                #D['Tx'] = S.tx.position[:,l]
                D['Tx'] = self.tx.points[l]
                if ext == '':
                    txrx = 'tx' + str('%0.3d' % l) + '-rx' + str('%0.3d' % k)
                else:
                    txrx = ext

                _filename = racine + txrx
                self.dcir[l][k] = _filename
                rep = rep + '/Tx' + str('%0.3d' % l)
                filename = pyu.getlong(_filename, rep)
                VCl = channel.VectChannel(self, l, k, False)
                CVC.append(VCl)
                if not VCl.fail:
                    #pdb.set_trace()
                    SCO = VCl.vec2scal()
                    CSCO.append(SCO)
                    SCA = VCl.vec2scalA(self.tx.A, self.rx.A, alpha=alpha)
                    CSCA.append(SCA)
                    ciro = SCO.applywavB(self.wav.sfg)
                    CIRo.append(ciro)
                    cira = SCA.applywavB(self.wav.sfg)
                    CIRa.append(cira)
                    #D['Rx'+k]    = self.rx.position[:,k]
                    D['Rx' + str(k)] = self.rx.points[int(k)]
                    D['to' + str(k)] = ciro.x
                    D['ta' + str(k)] = cira.x
                    D['ciro' + str(k)] = ciro.y
                    D['cira' + str(k)] = cira.y
                    spio.savemat(filename, D)
                else:
                    CSCO.append([])
                    CSCA.append([])
                    CIRo.append([])
                    CIRa.append([])

            if (store_level) & 1 == 1:
                self.CVC.append(CVC)
                self.CSCO.append(CSCO)
            if (store_level) & 4 == 4:
                self.CSCA.append(CSCA)
            if (store_level) & 8 == 8:
                self.CIRo.append(CIRo)
            if (store_level) & 16 == 16:
                self.CIRa.append(CIRa)

if (__name__ == "__main__"):
    doctest.testmod()
#    S = Simul()
#    narg = len(sys.argv)
#    if narg >=2:
#        S.load('where2.ini')
#        itx = int(sys.argv[1])
#        if narg == 4:
#            rx = range(int(sys.argv[2]),int(sys.argv[3]),1)
#        else:
#            srx = range(1,len(S.rx.points.keys())+1,1)
#
#        S.run2(itx,srx)
#    else:
#        S.load('where2.ini')
#   try:
#        filename = sys.argv[1]
#        S.filesimul= pyu.getshort(filename)
#        S.load()
        #S.filesimul='where2.simul'
        #S.filesimul='w2-Tx40.simul'
#    S.CIR(arange(Ntx),arange(Nrx),wav,store_level=16+8+4+2+1,alpha=alpha)
#        S.choose()
#        S.run()
#        S.save()
#    S.run()
#    param= {'band': 0.499,
#               'fc': 4.493,
#               'fe': 50,
#               'tresh': 3,
#               'Tw': 30}
#    wav = Waveform('impulse-generic',param)
##    S.save()
#    Ntx = S.tx.N
#    Nrx = len(S.filefield[0])
#    alpha = sqrt(1./30.0)
#    S.CIR(arange(Ntx),arange(Nrx),wav,store_level=16+8+4+2+1,alpha=alpha)
#   #S.indoor.interact()
#   #fitra = S.filetra[0][7]
#   #fitud = S.filetud[0][7]
#   #g0    = GrRay3D()
#   #g0.load(fitra,S.indoor)
#
#   #g1 = GrRay3D.GrRayTud()
#   #g1.load(S.filetud[0][7],S.sl,S.mat)
#   #S =Simulation(filestr='siradel-cut.str2',
#   #        fileslab='siradel.slab',
#   #        filemat='siradel.mat',
#   #        filepalch='def.palch',
#   #        filepatra='def.patra',
#   #        filefreq="def.freq",
#   #        filespaTx="Tx.spa",
#   #        filespaRx="Rx.spa",
#   #        fileantTx="dipole_XZ_NEC2.vsh3",
#   #        fileantRx="dipole_XZ_NEC2.vsh3")
#
#
##  # Propagation channel
##  VC=VectChannel(S,0,0,False)
##
##  # Transmission channel: SC0 (omnis antenna) & SCA (realistic antenna)
##  SCO = VC.vec2scal()
##  SCA = VC.vec2scalA(S.tx.A,S.rx.A)
##
##  # waveform to be transmitted
##  #param= {'band': 3,
##  #       'fc': 5,
##        #       'fe': 100,
##        #       'tresh': 10,
##        #       'Tw': 10}
##        #w = Waveform('impulse-generic',param)
##
##  # Wgamma
##  #f = w.sf.x
##  #ygamma = -1j*0.3/(4*pi*f)
##  #gamm  = bs.FUsignal(f,ygamma)
##  #w.sf=w.st.ftshift()
##  #Wgam  = w.sf*gamm
##  #wgam  = Wgam.ift()
##  #wgam  = Wgam.iftshift()
##
##  # simulated received signal
##  #ro    = SCO.applywavB(Wgam)  # omnis antenna
##  #    = SCA.applywavB(Wgam)  # monocone antenna
##
##  # Siradel measured received signal
##  #filemes = pyu.getlong('SIRADEL_08-07-31_P179.mat','measures')
##  M     = UWBMesure(1,1)
##  r3    = M.tdd.ch3
##
##  # construction of realistic transmitted waveform
##  rr    = M.RAW_DATA.tx
##  tv    = arange(1001)*0.005
##  wbis  = bs.TUsignal(tv,rr)
##  wbis.translate(-0.5)
##  Ewbis  = wbis.energy()
##  wbis.y = wbis.y/sqrt(Ewbis)
##  wbis.y = wbis.y-mean(wbis.y)
##  wbis.Yadd_zeros2l(1900)
##  wbis.Yadd_zeros2r(1100)
##  Wbis     = wbis.ft()
##  f        = Wbis.x
##  ygamma   = -1j*0.3/(4*pi*f)
##  gamm     = bs.FUsignal(f,ygamma)
##  Wgambis  = Wbis*gamm
##
##  robis   = SCO.applywavB(Wgambis)
##
##  rabis    = SCA.applywavB(Wgambis)  # monocone antenna
##  S.CIR(arange(3),arange(82),Wgambis,store_level=16+8+4+2+1,alpha=1.0)
##
##        for k in range(82):
##                I1 = hstack((II,S.CIRa.y))
##  #plot(ro.x,ro.y)
##  #plot(robis.x,robis.y,'r-')
##  #plot(ra.x,ra.y)
##  #plot(rabis.x,rabis.y,'b-')
##  #plot(t,r3,'g-')
##  #show()
##
##
############################
### Create Simulation object
### S = Simulation()
### S.run()
### VCg = VectChannel(S,0)
### alpha = pi*rand()
### beta  = pi*rand()
### gamma = pi*rand()
### Ta = MEulerAngle(alpha,beta,gamma)
### alpha = pi*rand()
### beta  = pi*rand()
### gamma = pi*rand()
### Tb = MEulerAngle(alpha,beta,gamma)
###
### VCl   = VCg2VCl(VCg,Ta,Tb)
###
### SCl   = VCl.vec2scal()
### SCAl  = VCl.vec2scalA(S.tx.A,S.rx.A)
###
### rAl = SCAl.applywavB(Wgam)
###
### figure()
### plot(rAl.x,rAl.y)
### ax=gca()
### ax.ticklabel_format(style='sci',axis='y')
### ax.yaxis.major.formatter.set_powerlimits((0,0))
### xlabel('retard (ns)',fontsize=28)
### ylabel('(V)',fontsize=28)
### xticks(fontsize=24)
### yticks(fontsize=24)
### title('r(t)',fontsize=30)
### show()
###
