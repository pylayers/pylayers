#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""

    This module run the electromagnetic simulation with an old
    version of the Ray Tracing tool

    Deprecated

Utility Functions
=================

.. autosummary::
    :toctree: generated/

     rename
     spafile

Pafreq class
=============

.. autosummary::
    :toctree: generated/

     Pafreq.__init__
     Pafreq.info
     Pafreq.load
     Pafreq.save
     Pafreq.gui




Simul class
============

.. autosummary::
    :toctree: generated/

    Simul.__init__
    Simul.gui
    Simul.updcfg
    Simul.clean
    Simul.clean_project
    Simul.save
    Simul.save_project
    Simul.load_project
    Simul.choose
    Simul.load
    Simul.layout
    Simul.show
    Simul.PL
    Simul.evalcir
    Simul.loadcir
    Simul.pltcir
    Simul.scatter
    Simul.info
    Simul.info2
    Simul.filtray
    Simul.showray
    Simul.show3l
    Simul._show3
    Simul.show3
    Simul.freq
    Simul.getlaunch
    Simul.gettra
    Simul.gettud
    Simul.launching
    Simul.tracing
    Simul.tratotud
    Simul.field
    Simul.run2
    Simul.run
    Simul.gt
    Simul.delay
    Simul.gr
    Simul.VC
    Simul.cir

"""
import doctest
import os
import re
import getopt
import sys
import shutil
import Tkinter, tkFileDialog
import time
import ConfigParser
import pdb
import cPickle
import numpy as np
import scipy as sp
import scipy.io as spio
import matplotlib.pylab as plt
import struct as stru
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
import pylayers.util.plotutil as plu
import pylayers.signal.waveform as wvf
import pylayers.signal.bsignal as bs
from pylayers.simul.radionode import RadioNode
from pylayers.util import easygui
from pylayers.antprop.slab import Slab, SlabDB, Mat, MatDB
# Handle Layout
from pylayers.gis.layout import Layout
# Handle Rays
#from   Channel import *
# Handle directory hierarchy
from pylayers.util.project import *
# Handle UWB measurements
from pylayers.measures import mesuwb as muwb
try:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi import mlab
except:
    print 'Layout:Mayavi is not installed'
import pdb


def rename(_filename,itx,irx,rep='./'):
    """ rename simulation file with proper itx/irx 

    Parameters
    ----------

    _filmename : short file name
    itx : transmitter index
    irx : transmitter index
    rep : directory

    See Also
    --------
    
    tratotud

    """
    filename = pyu.getlong(_filename,rep)
    f1  = re.sub('tx_[0-9]*','tx_'+str(itx),_filename)
    new = re.sub('rx_[0-9]*','rx_'+str(irx),f1)
    filenew = pyu.getlong(new,rep)
    os.rename(filename,filenew)
    return(new)

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



class Simul(PyLayers):
    """ Simulation Class

    Methods
    -------

    gui()
        graphical user interface
    choose()
        choose a simulation file in simuldir
    info()
        info about the simulation
    showray(itx,irx,iray)
        show a single ray
    help()
        help on using Simulation object
    freq()
        return the frequency base
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
    run(itx,irx)
        run simulation for links (itx,irx)
    cir(itx,irx,store_level=0,alpha=1.0)
        Channel Impulse Response calculation


    Attributes
    ----------

    fileconf

    filetauk
    filetang
    filerang

    tx
    rx

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

    """
    def __init__(self, _filesimul='default.ini'):
        self.filesimul = _filesimul
        self.config = ConfigParser.ConfigParser()
        self.config.add_section("files")
        self.config.add_section("frequency")
        self.config.add_section("waveform")
        self.config.add_section("output")

        self.dtang = {}
        self.drang = {}
        self.dtauk = {}
        self.dfield = {}
        self.dcir = {}
        self.output = {}

        #
        # Here was a nasty bug : Rule for the future
        #    "Always precise the key value of the passed argument"
        #
        # Mal nommer les choses, c'est ajouter au malheur du monde ( Albert Camus )
        #
        self.tx = RadioNode(name = '',
                            typ = 'tx',
                            _fileini = 'radiotx.ini',
                            _fileant = 'defant.vsh3',
                            )

        self.rx = RadioNode(name = '',
                            typ = 'rx',
                            _fileini = 'radiorx.ini',
                            _fileant = 'defant.vsh3',
                            )

        self.filefreq = "def.freq"

        self.progress = -1  # simulation not loaded

        self.filetang = []
        self.filerang = []
        self.filetauk = []
        self.filefield = []

        self.fileconf = "project.conf"
        self.cfield = []
        self.fGHz = np.linspace(2, 11, 181, endpoint=True)
        self.wav = wvf.Waveform()
        try:
            self.load(_filesimul)
        except:
            pass


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

    def updcfg(self):
        """ update simulation .ini config file 
        
        with values currently in use.
        """
        self.config.set("files", "struc", self.filestr)
        self.config.set("files", "conf", self.fileconf)
        self.config.set("files", "txant", self.tx.fileant)
        self.config.set("files", "rxant", self.rx.fileant)
        self.config.set("files", "tx", self.tx.fileini)
        self.config.set("files", "rx", self.rx.fileini)
        self.config.set("files", "mat", self.filematini)
        self.config.set("files", "slab", self.fileslabini)

        self.config.set("tud", "purc", str(self.patud.purc))
        self.config.set("tud", "nrmax", str(self.patud.nrmax))
        self.config.set("tud", "num", str(self.patud.num))

        #
        # frequency section
        #

        self.config.set("frequency", "fghzmin", self.fGHz[0])
        self.config.set("frequency", "fghzmax", self.fGHz[-1])
        self.config.set("frequency", "nf", str(len(self.fGHz)))
        #
        # waveform section
        #
        self.config.set("waveform", "tw", str(self.wav['twns']))
        self.config.set("waveform", "band", str(self.wav['bandGHz']))
        self.config.set("waveform", "fc", str(self.wav['fcGHz']))
        self.config.set("waveform", "thresh", str(self.wav['threshdB']))
        self.config.set("waveform", "type", str(self.wav['typ']))
        self.config.set("waveform", "fe", str(self.wav['feGHz']))
        #
        # output section
        #
        for k in self.output.keys():
            self.config.set("output",str(k),self.dout[k])

        # Initialize waveform
        self.wav = wvf.Waveform()
        # Update waveform 
        self.wav.read(self.config)
        self.save()

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
                filename = pyu.getlong(self.filelch[itx], pstruc['DIRLCH'])
                print filename
        if level > 1:
            for itx in range(self.tx.N):
                for irx in range(self.rx.N):
                    filename = pyu.getlong(self.filetra[itx][irx],
                                           pstruc(['DIRTRA']))
                    print filename
        if level > 2:
            for itx in range(self.tx.N):
                for irx in range(self.rx.N):
                    filename = pyu.getlong(self.filetud[itx][irx], pstruc['DIRTUD'])
                    print filename
                    filename = pyu.getlong(self.filetauk[itx][irx],pstruc['DIRTUD'])
                    print filename
        if level > 3:
            for itx in range(self.tx.N):
                for irx in range(self.rx.N):
                    filename = pyu.getlong(self.filetang[itx][irx], pstruc['DIRTUD'])
                    print filename
                    filename = pyu.getlong(self.filerang[itx][irx], pstruc['DIRTUD'])
                    print filename
                    filename = pyu.getlong(self.filefield[itx][irx], pstruc['DIRTUD'])
                    print filename


    def clean_project(self,verbose= True):
        """
        Clean Pyrpoject directory

        remove .lch, .tra, .field .tud, .tauk, .tang, .rang, <pyproject>/output

        remove [output] entries into .ini of self.filesimul


        Parameters
        ----------
 
        verbose : boolean
            Verbose mode on/off


        Returns
        -------
            Boolean:
                True if the project has been cleaned, False otherwise


        """

        if verbose:

            print "-----------------------------------"
            print "-----------------------------------"
            print "          WARNING                  "
            print "-----------------------------------"
            print "-----------------------------------"
            print "You are about to remove ALL previous computed raytracing files."
            print "If you decide to remove it, you will need to restart the entire \
    raytracing simulation to exploit simulation results"
            print "\n Do you want to remove these simulation files ? y/n"
            r=raw_input()

        else :
            r =='y'

        if r == 'y':
            inifile=self.filesimul
            try:
                path=os.getenv('BASENAME')
            except:
                print('Error : there is no project  directory in $BASENAME')



            dirlist=['output']
            extension=['.lch','.field','.tra','.tud','.tang','.rang','.tauk']
            rindic=False


            # remove file

            for d in dirlist:
                for ex in extension:
                    files = os.listdir(path +'/' +d )
                    for f in files:
                        if not os.path.isdir(path+'/'+d+'/'+f) and ex in f:
                            rindic=True
                            if verbose:
                                print f
                            os.remove(path+'/'+d+'/'+f)



                    if rindic:
                        if verbose:
                            print 'removed *' + ex +' from ' +d +'\n'
                        rindic=False


            # remove output into the self.filesimul ini file

            simcfg = ConfigParser.ConfigParser()
            simcfg.read(pyu.getlong(inifile,pstruc['DIRSIMUL']))
            simcfg.remove_section('output')
            f=open(pyu.getlong(inifile,pstruc['DIRSIMUL']),'wb')
            simcfg.write(f)
            f.close()
            self.dout = {}
            self.dlch = {}
            self.dtra = {}
            self.dtud = {}
            self.dtang = {}
            self.drang = {}
            self.dtauk = {}
            self.dfield = {}
            self.dcir = {}
            self.output = {}

            if verbose:
                print 'removed [output] entries into ' +inifile +'\n'
                print 'Project CLEANED'
            return True
        else :
            if verbose:
                print "clean project process ABORTED"
            return False


    def save(self):
        """ save simulation file

        Simulation files are .ini files which are saved in a dedicated
        directory basename/ini in the Project tree

        """
        filesimul = pyu.getlong(self.filesimul, "ini")
        fd = open(filesimul, "w")
        # getting current spa file if any
        try:
            self.config.set("files", "tx", self.tx.fileini)
        except:
            pass
        try:
            self.config.set("files", "rx", self.rx.fileini)
        except:
            pass
        self.config.write(fd)
        fd.close()
        # save tx
        self.tx.save()
        # save rx
        self.rx.save()
        # save slab and mat file
        # --
        # self.slab.save(self.fileslabini)
        # self.slab.mat.save(self.filematini)
        # --
        # fix bug #189
        #   slab is a member of S.L not of S anymore
        #   mat is a member of S.L.sl not of S.slab
        try:
            self.L.sl.save(self.fileslabini)
        except:
            pass
        try:
            self.L.sl.mat.save(self.filematini)
        except:
            pass

#    def saveold(self):
#        """ save simulation file

#        """
#        filesimul = pyu.getlong(self.filesimul, "ini")
#        fd = open(filesimul, "w")
#        config = ConfigParser.ConfigParser()

#        #
#        # files section
#        #
#        #config.add_section("files")
#        self.config.set("files", "conf", self.fileconf)
#        self.config.set("files", "struc", self.filestr)
#        self.config.set("files", "slab", self.fileslab)
#        self.config.set("files", "mat", self.filemat)

#        try:
#            self.config.set("files", "tx", self.tx.filespa)
#        except:
#            pass
#        try:
#            self.config.set("files", "rx", self.rx.filespa)
#        except:
#            pass
#        try:
#            self.config.set("files", "txant", self.tx.fileant)
#        except:
#            pass
#        try:
#            self.config.set("files", "rxant", self.rx.fileant)
#        except:
#            pass
#        self.config.set("files", "patra", self.filepatra)
#        self.config.set("files", "palch", self.filepalch)
#        self.palch.save()
#        self.patra.save()

#        #
#        # tud section
#        #

#        self.config.set("tud", "purc", self.patud.purc)
#        self.config.set("tud", "nrmax", self.patud.nrmax)
#        self.config.set("tud", "num", self.patud.num)

#        #
#        # frequency section
#        #
#        self.config.set("frequency", "fghzmin", self.freq[0])
#        self.config.set("frequency", "fghzmax", self.freq[-1])
#        self.config.set("frequency", "Nf", len(self.freq))

#        #
#        # output section
#        #
#        #filelch exists
#        if self.progress > 0:
#            #config.add_section("output")
#            for k in range(len(self.filelch)):
#                _fileout = "out" + "???"
#                filename = self.filelch[k]
#                self.config.set("launch", str(k + 1), filename)

#        # filetra exists
#        for k in range(len(self.filelch)):
#            if self.progress > 1:
#                #self.config.add_section("trace")
#                for l in arange(len(self.filetra[k])):
#                    filename = self.filetra[k][l]
#                    self.config.set("trace", "rx" + str(l + 1), filename)

#        # .tang exists
#        # .rang exists
#        # .tud exists
#            if self.progress > 2:
#                #config.add_section("tud")
#                #config.add_section("tang")
#                #config.add_section("rang")

#                for l in arange(len(self.filetud[k])):
#                    ftud = self.filetud[k][l]
#                    self.config.set("tud", "rx" + str(l + 1), ftud)

#                for l in arange(len(self.filetang[k])):
#                    ftang = self.filetang[k][l]
#                    self.config.set("tang", "rx" + str(l + 1), ftang)

#                for l in arange(len(self.filerang[k])):
#                    frang = self.filerang[k][l]
#                    self.config.set("rang", "rx" + str(l + 1), frang)

#        # .field exist
#        # .tauk exist
#            if self.progress > 3:
#                #config.add_section("tauk")
#                #config.add_section("field")
#                for l in arange(len(self.filetud[k])):
#                    ftauk = self.filetud[k][l]
#                    self.config.set("tauk", "rx" + str(l + 1), ftauk)

#                for l in arange(len(self.filefield[k])):
#                    ffield = self.filefield[k][l]
#                    self.config.set("field", "rx" + str(l + 1), ffield)

#        self.config.write(fd)
#        fd.close()

    def save_project(self):
        """ save Simulation files in a zipfile

        Simulation files are .ini files which are saved in a dedicated
        directory basename/ini in the Project tree

        """
        root = Tkinter.Tk()
        zipfileName = tkFileDialog.asksaveasfilename(parent=root,
                            filetypes = [("zipped file","zip")] ,
                            title="Save Project",
                            )
        pyu.zipd(basename,zipfileName)
        root.withdraw()
        print "Current project saved in", zipfileName

    def load_project(self):
        """ load Simulation files from a zipfile

        Simulation files are .ini files which are saved in a dedicated
        directory basename/ini in the Project tree

        """
        root = Tkinter.Tk()
        zipfileName= tkFileDialog.askopenfile(parent=root,
                                            mode='rb',
                                            title='Choose a project')
        dirname = tkFileDialog.askdirectory(parent=root,
                                    initialdir=basename,
                                    title='Please select a directory',
                                    mustexist=0)
        pyu.unzipd(dirname,zipfileName)
        root.withdraw()
        print "Current project loaded in", dirname

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
        try:
            _filetx = self.config.get("files", "tx")
        except:
            raise NameError('Error in section tx from '+ _filesimul)

        try:
            _filerx = self.config.get("files", "rx")
        except:
            raise NameError('Error in section rx from '+ _filesimul)


        try:
            _fileanttx = self.config.get("files", "txant")
        except:
            raise NameError('Error in section txant from '+ _filesimul)

        try:
           _fileantrx = self.config.get("files", "rxant")
        except:
            raise NameError('Error in section rxant from '+ _filesimul)

        try:
            self.tx = RadioNode(name = '',
                                typ = 'tx',
                                _fileini = _filetx,
                                _fileant = _fileanttx,
                                _filestr = self.filestr)

            self.rx = RadioNode(name = '',
                                typ = 'rx',
                                _fileini = _filerx,
                                _fileant = _fileantrx,
                                _filestr = self.filestr)
        except:
            raise NameError('Error during Radionode load')
#
# Load Layout
# 
        try:
            self.L = Layout(self.filestr)
        except:
            raise NameError('Layout load error')

#
# Frequency base
#
        if "frequency" in sections:
            try:
                self.fGHz = np.linspace(float(self.config.getfloat("frequency", "fghzmin")),
                                        float(self.config.getfloat("frequency", "fghzmax")),
                                        int(self.config.getint("frequency", "nf")),
                                        endpoint=True)
            except:
                raise NameError('Error in section frequency from '+ _filesimul)
            # update .freq file in tud directory

            filefreq = pyu.getlong(self.filefreq, pstruc['DIRTUD'])
            fd = open(filefreq, "w")
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
                _filename  =  self.config.get("output", itx)
                self.dout[int(itx)] = _filename
                filename = pyu.getlong(_filename, "output")
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

    def layout(self, _filestruc):
        """ load a layout in the simulation oject

        Parameters
        ----------

        _filestruc : string
            short file name of the Layout object

        Examples
        --------

        >>> from pylayers.simul.simulem import *
        >>> S = Simul()
        >>> S.layout('defstr.ini')

        """
        self.filestr = _filestruc

        self.L = Layout(_filestruc)
        # update config
        self.config.set("files", "struc", self.filestr)
        self.save()

    def show(self, itx=[-1], irx=[-1], furniture=True, s=8, c='b', traj=False, num=False,fig=[],ax=[]):
        """ show simulation

            Parameters
            -----------
            itx        : list of tx indices
            irx        : list of rx indices
            furniture  : boolean for METALIC furniture display
            s          : scale fir scatter plot  (default 8)
            c          : color for scatter plot  (default 'b')
            traj       : boolean  (def False)
            num        : boolean  (def False)
                display a number


            Examples
            --------
            >>> import matplotlib.pyplot as plt
            >>> from pylayers.simul.simulem import *
            >>> S = Simul()
            >>> S.load('w1.ini')
            >>> S.L.loadfur('FurW1.ini')
            >>> S.show()
            >>> plt.show()



        """
        if type(itx) == int:
            itx = [itx]
        if type(irx) == int:
            irx = [irx]


        if fig ==[]:
            fig = plt.gcf()
        if ax==[]:
            ax = fig.gca()

        #self.L.display['scaled']=False
        fig,ax=self.L.showG('s',fig=fig,ax=ax, aw=True)
        #
        if furniture:
            if 'lfur' in self.L.__dict__:
                for fur in self.L.lfur:
                    if fur.Matname == 'METAL':
                        fur.show(fig, ax)
            else:
                print "Warning : no furniture file loaded"

        if irx[0] == -1:
            ax.scatter(self.rx.position[0,:],
                       self.rx.position[1,:], c='b', s=s, alpha=0.5)
            #ax.scatter(self.rx.position[0,0],self.rx.position[0,1],c='k',s=s,linewidth=0)
            #ax.scatter(self.rx.position[1,0],self.rx.position[1,1],c='b',s=s,linewidth=0)
            #ax.scatter(self.rx.position[2,0],self.rx.position[2,1],c='g',s=s,linewidth=0)
            #ax.scatter(self.rx.position[3,0],self.rx.position[3,1],c='c',s=s,linewidth=0)
        else:
            for k in irx:
                ax.scatter(self.rx.position[0,k - 1],
                           self.rx.position[1,k - 1], c='b', s=s, alpha=0.5)
                if num:
                    ax.text(self.rx.position[0,k - 1],
                            self.rx.position[1,k - 1],
                            str(k), color='blue')

        if itx[0] == -1:
            ax.scatter(self.tx.position[0,:],
                       self.tx.position[1,:], c='r', s=s)
        else:
            if traj:
                cpt = 1
            for k in itx:
                ax.scatter(self.tx.position[0,k - 1],
                           self.tx.position[1,k - 1],
                           c=c, s=s, linewidth=0)
                if num:
                    if traj:
                        ax.text(self.tx.position[0,k - 1],
                                self.tx.position[1,k - 1],
                                str(cpt), color='black')
                        cpt = cpt + 1
                    else:
                        ax.text(self.tx.position[0,k - 1],
                                self.tx.position[1,k - 1],
                                str(k), color='black')

        return (fig,ax)
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

            td.append(d)
            tEa.append(Ea)
            tEo.append(Eo)

        plt.semilogx(td, 10 * np.log10(tEa), 'xr')
        plt.semilogx(td, 10 * np.log10(tEo), 'xb')
        plt.show()
        return td, tEa, tEo

    def evalcir(self,cutoff=4,algo='new'):
        """
        Parameters
        ----------

        S
        tx
        rx
        wav
        cutoff

        """

        crxp =-1
        ctxp =-1
        tcir = {}
        tx = self.tx.position
        Ntx = len(tx[0])
        rx = self.rx.position
        Nrx = len(rx[0])

        #for kt in range(1,Ntx-1):
        #print kt+1
        kt=0
        tcir[kt] = {}
        t = np.array([self.tx.position[0,kt],self.tx.position[1,kt],self.tx.position[2,kt]])
        for kr in range(Nrx):
            if (np.mod(kr,10)==0):
                print kr+1
            r = np.array([self.rx.position[0,kr],self.rx.position[1,kr],self.rx.position[2,kr]])
            ctx = self.L.pt2cy(t)
            crx = self.L.pt2cy(r)
            if (ctx<>ctxp)|(crx<>crxp):
                Si  = signature.Signatures(self.L,ctx,crx)
                ctxp = ctx
                crxp = crx
                Si.run4(cutoff=cutoff,algo=algo)
            r2d = Si.rays(t,r)
            #r2d.show(S.L)

            r3d = r2d.to3D(self.L)
            r3d.locbas(self.L)
            r3d.fillinter(self.L)
            Ct  = r3d.eval(self.fGHz)
            sca = Ct.prop2tran(self.tx.A,self.rx.A)
            cir = sca.applywavB(self.wav.sfg)
            tcir[kt][kr] = cir
        return(tcir)

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

        filecir = pyu.getlong(_filecir, pstruc['DIRCIR']+'/Tx' + ext)
        D = spio.loadmat(filecir)

        kxa = 'ta' + str(irx)
        kya = 'cira' + str(irx)

        kxo = 'to' + str(irx)
        kyo = 'ciro' + str(irx)

        cira = bs.TUsignal(D[kxa], D[kya][:, 0])
        ciro = bs.TUsignal(D[kxo], D[kyo][:, 0])

        return(cira, ciro)

    def pltcir(self, itx=1, irx=1, mode='linear', noise=False, color='b',format='a',fig=[],ax=[]):
        """ plot Channel Impulse Response

        Parameters
        ----------

        itx : Tx index
        irx : Rx index
        mode : str
            {'linear','dB'}
            noise : boolean
        color : string
            default 'b'

        >>> from pylayers.simul.simulem import *
        >>> S = Simul()
        >>> S.load('where2.ini')
        >>> S.run(1,1)
        >>> S.pltcir(1,1,mode='linear',noise=False,color='k')

        """

        if fig ==[]:
            fig = plt.gcf()
        #if ax==[]:
        #    ax = fig.gca()

        _filecir = self.dcir[itx][irx] + '.mat'
        filecir = pyu.getlong(_filecir, pstruc['DIRCIR']+'/Tx' + str('%0.3d' % itx))
        D = spio.loadmat(filecir)
        ax = fig.add_subplot('211')

        fig,ax=self.show(itx, irx,fig=fig,ax=ax)
        ax=fig.add_subplot('212')
        if 'a' in format :
            kxa = 't'
            kya = 'cir'
            ta = D[kxa]
            Tobs = ta[-1] - ta[0]
            te = ta[1] - ta[0]
            if noise:
                na = bs.Noise(Tobs + te, 1. / te)
                naf = na.gating(4.493, 0.5)
            cira = bs.TUsignal(ta, D[kya][:, 0])


        if 'o' in format:
            kxo = 'to' + str(irx)
            kyo = 'ciro' + str(irx)
            to = D[kxo]
            Tobs = to[-1] - to[0]
            te = to[1] - to[0]
            if noise:
                no = bs.Noise(Tobs + te, 1. / te)
                nof = no.gating(4.493, 0.5)
            ciro = bs.TUsignal(to, D[kyo][:, 0])

        if mode == 'linear':
            #plt.plot(ta,naf.y,color='k',label='Noise')
            plt.plot(ta, D[kya], label='Rx ' + str(irx), color=color)
            plt.xlabel('Time (ns)')

            '''if noise:
                naf.plot(col='k')
            cira.plot(col=color)'''
        else:
            '''if noise:
                naf.plotdB(col='k')
            cira.plotdB()'''
            plt.plot(ta, 20 * np.log10(abs(D[kya])), label='Rx ' + str(irx), color=color)
            plt.xlabel('Time (ns)')
#        plt.legend()
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
            print ">>> S.layout(filename.str,filematini,filematini) "
            print "default files exists for filematini and fileslabini "

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
            print "Distance (m) :", 0.3 / self.delay(itx, irx)
            print ""
            if itx in self.dlch.keys():
                print "-----"
                print "Launching "
                print "-----"
                print " ", self.dlch[itx]
            if irx in self.dtra[itx].keys():
                print "-----"
                print "Tracing "
                print "-----"
                print " ", self.dtra[itx][irx]
                gr = GrRay3D()
                gr.load(self.dtra[itx][irx], self.L)
                gr.info()
            if irx in self.dtud[itx].keys():
                print "-----"
                print "Tud parameters "
                print "-----"
                print " ", self.dtud[itx][irx]
                print " ", self.dtang[itx][irx]
                print " ", self.drang[itx][irx]
                gt = GrRay3D.GrRayTud()
                gt.load(self.dtud[itx][irx],
                        self.dtang[itx][irx],
                        self.drang[itx][irx], self.sl)
            if irx in self.dtauk[itx].keys():
                print self.dtauk[itx][irx]
                print self.dfield[itx][irx]
                VC = self.VC(itx, irx)
            if irx in self.dcir[itx].keys():
                print self.dcir[itx][irx]

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

    def _show3(self,rays=[],newfig = False,**kwargs):
        """ display of the simulation configuration
            using Mayavi

        Parameters
        ----------

        rays: Ray3d object :
            display the rays of the simulation
        newfig : boolean (default : False)
        kwargs of Rays.show3()


        see also
        --------

        pylayers.gis.layout
        pylayers.antprop.antenna
        pylayers.antprop.rays

        """
        Atx = self.tx.A
        Arx = self.rx.A
        Ttx = self.tx.orientation
        Trx = self.rx.orientation
        ptx = self.tx.position
        prx = self.rx.position

        self.L._show3(newfig=False,opacity=0.7)
        Atx._show3(T=Ttx.reshape(3,3),po=ptx,
            title=False,colorbar=False,newfig=False)
        Arx._show3(T=Trx.reshape(3,3),po=prx,
            title=False,colorbar=False,newfig=False)
        if rays != []:
            rays._show3(**kwargs)



    def show3(self,rays=[],**kwargs):
        """ geomview display of the simulation configuration

        Parameters
        ----------

        centered : boolean
            center the scene if True
        bdis  : boolean
            display local basis

        """
        try:
            self.tx.save()
        except:
            print('tx set is no defined')
        try:
            self.rx.save()
        except:
            print('rx set is no defined')
        _filename = self.filesimul.replace('.ini', '.off')
        filename = pyu.getlong(_filename, pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        try:
            sttx = "{<" + self.tx.filegeom + "}\n"
        except:
            sttx = "\n"
        try:
            strx = "{<" + self.rx.filegeom + "}\n"
        except:
            strx = "\n"
        try:
            stst = "{<" + self.L.filegeom + "}\n"
        except:
            stst = "\n"
        fo.write(sttx)
        fo.write(strx)
        fo.write(stst)
        fo.write("{</usr/share/geomview/geom/xyz.vect}\n")
        if rays !=[]:
            kwargs['bdis']=False
            kwargs['L']=self.L
            kwargs['centered']=False
            fo.write("{<" + rays.show3(**kwargs) + "}")

        fo.close()


        command = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
        os.system(command)

    def run(self, link, cirforce=True,verbose=False,cutoff=4):
        """ run the simulation for 1 tx and a set of rx

            Parameters
            ----------

            itx      : tx index
            srx      : list of rx index
            cirforce : boolean

            Warnings
            --------

            index point start with 1

            Example
            -------

            >>> from pylayers.simul.simulem import *
            >>> itx = 1
            >>> srx = [1,2,3]
            >>> S   = Simul()
            >>> S.load('where2.ini')
            >>> out = S.run(itx,srx)


        """
        
        # get file prefix

        link = DLink(force=True,L=self.L,fGHz=self.fGHz, verbose=False)
        prefix = self.filesimul.replace('.ini', '') 
        lsig = []
        for k,il in enumerate(link):
            tx = self.tx.points[il[0]]
            rx = self.rx.points[il[1]]
            ctx = S.L.pt2cy(tx)
            crx = S.L.pt2cy(rx)
            _filecir = prefix +'-cir-'+str(k)+'-'+str(link)+'-'+str((ctx,crx))
            D = {}
            D['Tx'] = tx
            D['Rx'] = rx

            
            link.a = tx
            link.b = rx
            ak,tauk = link.eval(verbose=False,diffrction=True)
            if (ctx,crx) not in lsig:
                Si  = signature.Signatures(S.L,ctx,crx)
                #
                # Change the run number depending on
                # the algorithm used for signature determination
                #
                Si.run4(cutoff=cutoff)
                # keep track and save signature
                _filesir = prefix + '-sig-'+str((ctx,crx))
                fd = open(filesig,'w')
                pickle.dump(Si,filesig)
                fd.close()
                lsig.appeng((ctx,crx))
                Si.dump(S.L,(ctx,crx))



            r2d = Si.rays(tx,rx)
            r2d.show(S.L)
    
            r3d = r2d.to3D()
            r3d.locbas(S.L)
            r3d.fillinter(S.L)

            Ct  = r3d.eval(S.freq)
            sco = Ct.prop2tran(a='theta',b='phi')
            sca = Ct.prop2tran(a=S.tx.A,b=S.rx.A)

            ciro = sco.applywavB(self.wav.sfg)
            cira = sca.applywavB(self.wav.sfg)

            D['to'] = ciro.x
            D['ciro'] = ciro.y
            D['t'] = cira.x
            D['cir'] = cira.y

            filename = pyu.getlong(_filename, cirdir)
            spio.savemat(filename, D)
            
    def delay(self, itx, irx):
        """ calculate LOS link delay

            Parameters
            ----------

            itx
            irx

            Returns
            -------

            delay : float 
                delay in ns

        """
        tx = self.tx.points[itx]
        rx = self.rx.points[irx]
        df = tx - rx
        dist = np.sqrt(np.dot(df, df))
        return(dist / 0.3)

if (__name__ == "__main__"):
    #plt.ion()
    doctest.testmod()
