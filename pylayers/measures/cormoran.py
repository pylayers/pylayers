# -*- coding:Utf-8 -*-
"""

This module handles CORMORAN measurement data

CorSer Class
============

.. autosummary::
    :toctree: generated/

    cor_log

CorSer init and load data
-------------------------

.. autosummary::
    :toctree: generated/

    CorSer.__init__
    CorSer.__repr__
    CorSer._loadBS
    CorSer._loadTCR
    CorSer._load_offset_dict
    CorSer._loadbody
    CorSer._loadcam
    CorSer._loadhkb
    CorSer._loadinfranodes
    CorSer.loadlog


CorSer exploit data
-------------------

.. autosummary::
    :toctree: generated/

    CorSer.getdevp
    CorSer.getlink
    CorSer.getlinkd
    CorSer.getlinkp
    CorSer.devmapper
    CorSer.export_csv


CorSer visualization tools
--------------------------

.. autosummary::
    :toctree: generated/


2D  vizualization
~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    CorSer.plot
    CorSer.pltmob
    CorSer.plttcr
    CorSer.pltvisi
    CorSer.animhkb
    CorSer.animhkbAP
    CorSer.imshow
    Corser.showlink
    CorSer.pltgt


3D vizualization
~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    CorSer._show3
    CorSer._show3i
    CorSer.imshowvisibility
    CorSer.imshowvisibility_i
    CorSer.anim



Video Vizualization
~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    CorSer.snapshot
    CorSer.snapshots
    CorSer.vlc


CorSer Internal processing
--------------------------

.. autosummary::
    :toctree: generated/


Time Alignement
~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/


    CorSer.align
    CorSer.offset_setter
    CorSer.offset_setter_video
    CorSer._align_on_devdf
    CorSer._align_devdf_on_hkb
    CorSer._apply_offset'
    CorSer._save_data_off_dict
    CorSer._save_offset_dict

Data frame transformation
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    CorSer.compute_visibility
    CorSer._computedevpdf
    CorSer._computedistdf
    CorSer._distancematrix
    CorSer._visiarray


Vizualizing processing
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    CorSer.visidev
    CorSer.visidev2
    CorSer.__refreshshow3i

Notes
-----

Useful members

distdf : distance between radio nodes (122 columns)
devdf  : device data frame

"""

import os
import pdb
import sys
import pandas as pd
import numpy as np
import numpy.ma as ma
import scipy.io as io
from pylayers.util.project import *
from pylayers.util.pyutil import *
#import  pylayers.util.mayautil as myu
from pylayers.mobility.ban.body import *
from pylayers.gis.layout import *
import pylayers.antprop.antenna as antenna

from matplotlib.widgets import Slider, CheckButtons, Button, Cursor
from pylayers.signal.DF import *

from moviepy.editor import *
from skimage import img_as_ubyte
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pickle

#Those lines handle incompatibility between mayavi and VTK
#and redirect noisy warning message into a log file
# import vtk
# output=vtk.vtkFileOutputWindow()
# output.SetFileName("mayaviwarninglog.tmp")
# vtk.vtkOutputWindow().SetInstance(output)


def cor_log(short=True):
    """ display cormoran measurement campaign logfile

    Parameters
    ----------

    short : boolean
        enable short version

    Examples
    --------

    >>> from pylayers.measures.cormoran import *
    >>> cor_log(short=True)

    """

    filelog = os.path.join(os.environ['CORMORAN'],'RAW','Doc','MeasurementLog.csv')
    log = pd.read_csv(filelog)
    if short :
        log['day'] =  [x.split('/')[0] for x in log['Date'].values]
        log['serie']=log['Meas Serie']
        return log[['serie','day','Subject','techno','Short Notes']]
    else:
        return log

def time2npa(lt):
    """ convert pd.datetime.time to numpy array

    Parameters
    ----------

    lt : pd.datetime.time

    Returns
    -------

    ta : numpy array
        time in seconds

    """
    ta = (lt.microsecond*1e-6+
    lt.second+
    lt.minute*60+
    lt.hour*3600)
    return(ta)

class CorSer(PyLayers):
    """ Handle CORMORAN measurement data

    Hikob data handling from CORMORAN measurement campaign

    11/06/2014
        single subject (Bernard and Nicolas)

    12/06/2014
        several subject (Jihad, Eric , Nicolas)
    """

    def __init__(self,serie=6,day=11,source='CITI'):
        """
        Parameters
        ----------

        serie : int
        day : int
        source : string

        Notes
        -----

        The environment variable CORMORAN is indicating the location of data directory

        """
        assert (day in [11,12]),"wrong day"

        try:
            self.rootdir = os.environ['CORMORAN']
        except:
            raise NameError('Please add a CORMORAN environement variable \
                            pointing to the data')

        # infos
        self.serie = serie
        self.day = day
        self.loadlog()

        if day == 11:
            if serie in [7,8]:
                print 'Serie '+str(serie) + ' has no hkb data and will not be loaded'
                exit
        if day ==12:
            if serie in [17,18,19,20]:
                raise AttributeError('Serie '+str(serie) + \
                                     ' has no hkb data and will not be loaded')
        #Measures
        if day==11:
            self.stcr = [1,2,3,4,10,11,12,32,33,34,35,9,17,18,19,20,25,26]
            self.shkb = [5,6,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
            self.sbs  = [5,6,7,8,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
            self.mocap = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]
            self.mocapinterf=[]
        if day==12:
            self.stcr = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
            self.shkb = [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
            self.sbs  = [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
            self.mocap =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
            self.mocapinterf = [5,6,7,8,13,14,15,16,21,22,23,24,]

        self.typ=''

        # HIKOB
        if serie in self.shkb:
            self._loadhkb(serie=serie,day=day,source=source)

        # IR-UWB  TCR
        if serie in self.stcr:
            self._loadTCR(serie=serie,day=day)

        # BeSpoon
        if serie in self.sbs:
            self._loadBS(serie=serie,day=day)

        # set filename
        if self.typ=='FULL':
            self._filename = 'Sc' + self.scenario + '_S' + str(self.serie) + '_R' + str(self.run) + '_' + self.typ.capitalize()
        else:
            self._filename = 'Sc' + self.scenario + '_S' + str(self.serie) + '_R' + str(self.run) + '_' + self.typ


        #Layout
        self.L= Layout('MOCAP-small2.ini')


        # Load Infrastructure Nodes
        self._loadinfranodes()

        # Load cameras
        self._loadcam()

        #BODY & interferers
        self.subject = str(self.log['Subject'].values[0]).split(' ')
        #filter typos in  self.subject
        self.subject = filter(lambda x : len(x)!=0,self.subject)
        if 'Jihad' in self.subject :
            uj = self.subject.index('Jihad')
            self.subject[uj]='Jihan'

        if serie in self.mocap :
            self._loadbody(serie=serie,day=day)
            self._distancematrix()
            self._computedevpdf()
            if isinstance(self.B,dict):
                for b in self.B:
                    self.B[b].traj.Lfilename=copy.copy(self.L.filename)
            else :
                self.B.traj.Lfilename=copy.copy(self.L.filename)

        # reference time is tmocap
        self.tmocap = self.B[self.subject[0]].time

        # load offset dict
        self.offset= self._load_offset_dict()

        ########################
        #realign Radio on mocap
        ########################
        # 1 - Resample radio time => mocap time
        # 2 - (if available) apply offset

        if ('BS' in self.typ) or ('FULL' in self.typ):
            print '\nBS data frame index: ',
            self._align_on_devdf(typ='BS')
            print 'Align on mocap OK...',
            try:
                self._apply_offset('BS')
                print 'time-offset applied OK'
            except:
                print 'WARNING time-offset NOT applied'
                print ('No BS offset not yet set => use self.offset_setter ')

        if ('TCR' in self.typ) or ('FULL' in self.typ):
            print '\nTCR data frame index:', 
            self._align_on_devdf(typ='TCR')
            print 'Align on mocap OK...',
            try:
                self._apply_offset('TCR')
                print 'time-offset applied OK'
            except:
                print 'WARNING time-offset NOT applied'
                print ('No TCR offset not yet set => use self.offset_setter')


        if ('HK' in self.typ) or ('FULL' in self.typ):
            print '\nHKB data frame index:',
            self._align_on_devdf(typ='HKB')
            print 'Align on mocap OK...',
            try:
                # self._apply_offset('HKB')
                print 'time-offset applied OK'
            except:
                print 'WARNING time-offset NOT applied'
                print ('No HKB offset not yet set => use self.offset_setter')

        print '\nCreate distance Dataframe...',
        self._computedistdf()
        print 'OK',

    def __repr__(self):
        st = ''
        st = st + 'filename : ' + self._filename + '\n'
        st = st + 'filewear : ' + self.filewear + '\n'
        st = st + 'filebody : ' + self.filebody + '\n'
        st = st + 'filemocap : ' + self.filemocap + '\n'
        st = st + 'Day : '+ str(self.day)+'/06/2014'+'\n'
        st = st + 'Serie : '+ str(self.serie)+'\n'
        st = st + 'Scenario : '+str(self.scenario)+'\n'
        st = st + 'Run : '+ str(self.run)+'\n'
        st = st + 'Type : '+ str(self.typ)+'\n'
        st = st + 'Original Video Id : '+ str(self.video)+'\n'
        st = st + 'Subject(s) : '

        for k in self.subject:
            st = st + k + ' '
        st = st + '\n\n'

        st = st+'Body available: ' + str('B' in dir(self)) + '\n\n'

        try :
            st = st+'BeSPoon : '+self._fileBS+'\n'
        except:
            pass
        try :
            st = st+'HIKOB : '+self._filehkb+'\n'
        except:
            pass
        try :
            st = st+'TCR : '+self._fileTCR+'\n'
        except:
            pass
        return(st)


    # @property
    # def dev(self):
    #     """ display device techno, id , id on body, body owner,...
    #     """
        
    #     title = '{0:21} | {1:7} | {2:8} | {3:10} '.format('Name in Dataframe', 'Real Id', 'Body Id', 'Subject')
    #     print title + '\n' + '-'*len(title) 
    #     if ('HK' in self.typ) or ('FULL' in self.typ):
    #         hkbkeys = self.idHKB.keys()
    #         hkbkeys.sort()
    #         for d in hkbkeys:
    #             dev = self.devmapper(self.idHKB[d],'HKB')
    #             print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
    #     if ('TCR' in self.typ) or ('FULL' in self.typ):
    #         tcrkeys = self.idTCR.keys()
    #         tcrkeys.sort()
    #         for d in tcrkeys:
    #             dev = self.devmapper(self.idTCR[d],'TCR')
    #             print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])

    @property
    def dev(self):
        """ display device techno, id , id on body, body owner,...
        """

        title = '{0:21} | {1:7} | {2:8} | {3:10} '.format('Name in Dataframe', 'Real Id', 'Body Id', 'Subject')
        print title + '\n' + '='*len(title) 
        # access points HKB
        for d in self.din:
            if ('HK' in d) :
                dev = self.devmapper(d,'HKB')
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
        if 'FULL' in self.typ:
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')
        for d in self.din:
            if ('BS' in d) :
                dev = self.devmapper(d,'BS')
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
        if 'FULL' in self.typ:
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')

        # access points TCR
        for d in self.din:
            if ('TCR' in d)  :
                dev = self.devmapper(d,'TCR')
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
        print '{0:66}'.format('-'*len(title) )
        #device per RAT per body
        for b in self.B:
            if b not in self.interf:
                #HKB per body
                for d in self.B[b].dev.keys():

                    if ('HK' in d):
                        dev = self.devmapper(d,'HKB')
                        print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
                #bespoon
                if ('FULL' in self.typ) or ('HKB' in self.typ):
                    print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')
                for d in self.B[b].dev.keys():
                    if ('BS' in d):
                        dev = self.devmapper(d,'BS')
                        print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
                # print '{0:66}'.format('-'*len(title) )
                #TCR per body
                if 'FULL' in self.typ:
                    print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')
                for d in self.B[b].dev.keys():
                    if ('TCR' in d):
                        dev = self.devmapper(d,'TCR')
                        print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
                print '{0:66}'.format('-'*len(title) )

    @property
    def ant(self):
        """ display device techno, id , id on body, body owner,...
        """

        title = '{0:21} | {1:7} | {2:8} | {3:10} '.format('Name in Dataframe', 'Real Id', 'Body Id', 'Subject')
        print title + '\n' + '='*len(title) 
        # access points HKB
        for d in self.din:
            if ('HK' in d) :
                dev = self.devmapper(d,'HKB')
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
        if 'FULL' in self.typ:
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')
        for d in self.din:
            if ('BS' in d) :
                dev = self.devmapper(d,'BS')
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
        if 'FULL' in self.typ:
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')

        # access points TCR
        for d in self.din:
            if ('TCR' in d)  :
                dev = self.devmapper(d,'TCR')
                print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
        print '{0:66}'.format('-'*len(title) )
        #device per RAT per body
        for b in self.B:
            if b not in self.interf:
                #HKB per body
                for d in self.B[b].dev.keys():

                    if ('HK' in d):
                        dev = self.devmapper(d,'HKB')
                        print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
                #bespoon
                if ('FULL' in self.typ) or ('HKB' in self.typ):
                    print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')
                for d in self.B[b].dev.keys():
                    if ('BS' in d):
                        dev = self.devmapper(d,'BS')
                        print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
                # print '{0:66}'.format('-'*len(title) )
                #TCR per body
                if 'FULL' in self.typ:
                    print '{0:21} | {1:7} | {2:8} | {3:10} '.format('','','','')
                for d in self.B[b].dev.keys():
                    if ('TCR' in d):
                        dev = self.devmapper(d,'TCR')
                        print '{0:21} | {1:7} | {2:8} | {3:10} '.format(dev[0],dev[1],dev[2],dev[3])
                print '{0:66}'.format('-'*len(title) )



    def _loadcam(self):
        """ load camera position 

        Returns
        -------

        update self.cam

        """

        self.cam = np.array([
            [-6502.16643961174,5440.97951452912,2296.44437108561],
            [-7782.34866625776,4998.47624994092,2417.5861326688],
            [8308.82897665828,3618.50516290547,2698.07710953287],
            [5606.68337709102,-6354.17891528277,2500.27779697402],
            [-8237.91886515041,-2332.98639475305,4765.31798299242],
            [5496.0942989988,6216.91946236788,2433.30012872688],
            [-8296.19706598514,2430.07325486109,4794.01607841197],
            [7718.37527064615,-4644.26760522485,2584.75330667172],
            [8471.27154730777,-3043.74550832061,2683.45089703377],
            [-8213.04824602894,-4034.57371591121,2368.54548665579],
            [-7184.66711497403,-4950.49444503781,2317.68563412347],
            [7531.66103727189,5279.02353243886,2479.36291603544],
            [-6303.08628709464,-7057.06193926342,2288.84938553817],
            [-5441.17834354692,6637.93014323586,2315.15657646861],
            [8287.79937470615,59.1614281340528,4809.14535447027]
            ])*1e-3


    def _loadinfranodes(self):
        """ load infrastructure nodes


nico

                        A4
                    mpts[6,7,8]
                        X

            A3                     A1
        mpts[9,10,11]        mpts[3,4,5]
            X                      X

                        A2
                    mpts[0,1,2]
                        X


        TCR = mpts[0,3,6,9]
        HKB = mpts[1,2,
                   4,5,
                   7,8,
                   10,11]

bernard


                        A3
                    mpts[3,4,5]
                        X

            A2                     A4
        mpts[6,7,8]        mpts[0,1,2]
            X                      X

                        A1
                    mpts[9,10,11]
                        X


        TCR = mpts[0,3,6,9]
        HKB = mpts[1,2,
                   4,5,
                   7,8,
                   10,11]


        """

        filename = os.path.join(self.rootdir,'RAW','11-06-2014','MOCAP','scene.c3d')

        print "\nload infrastructure node position:",
        a, self.infraname, pts, i = c3d.ReadC3d(filename)

        pts = pts/1000.
        mpts = np.mean(pts, axis=0)
        self.din={}
        if ('HK'  in self.typ) or ('FULL' in self.typ):
            uhkb = np.array([[1,2], [4,5], [7,8], [10,11]])
            mphkb = np.mean(mpts[uhkb], axis=1)

            self.din.update(
                {'HKB:1':{'p' : mphkb[3],
                          # 'T' : np.eye(3),
                          's3off' : 0.},

                 'HKB:2':{'p' : mphkb[2],
                          # 'T': np.array([[-0.44807362,  0.89399666,  0.],
                          #                [-0.89399666, -0.44807362,  0.],
                          #                [ 0.,0.,1.        ]]),
                          's3off':0.}      ,
                 'HKB:3':{'p':mphkb[1],
                          # 'T':array([[-0.59846007, -0.80115264,  0.],
                          #            [ 0.80115264, -0.59846007,  0.],
                          #            [ 0.,0.,  1.]]),
                          's3off':0.},
                 'HKB:4':{'p':mphkb[0],
                          # 'T':array([[-0.44807362, -0.89399666,  0.],
                          #            [ 0.89399666, -0.44807362,  0.],
                          #            [ 0.,0.,  1.]]),
                          's3off':0.}
                 })

        # TCR:31 is the coordinator which was not captured.
        # The position has been determined via optimization
        if ('TCR' in self.typ) or ('FULL' in self.typ):
            self.din.update({'TCR:32':{'p':mpts[9],
                                       'T':np.eye(3),
                                       's3off':0.1},
                 'TCR:24':{'p':mpts[6],
                           # 'T': np.array([[-0.44807362,  0.89399666,  0.],
                           #               [-0.89399666, -0.44807362,  0.],
                           #               [ 0.,0.,1.        ]]),
                           's3off':0.1},
                 'TCR:27':{'p':mpts[3],
                           # 'T':array([[-0.59846007, -0.80115264,  0.],
                           #           [ 0.80115264, -0.59846007,  0.],
                           #           [ 0.,0.,  1.]]),
                           's3off':0.1},
                 'TCR:28':{'p':mpts[0],
                           # 'T':array([[-0.44807362, -0.89399666,  0.],
                           #           [ 0.89399666, -0.44807362,  0.],
                           #           [ 0.,0.,  1.]]),
                           's3off':0.1},
                 'TCR:31':{'p':array([1.7719,-3.2655,1.74]),
                           # 'T':array([[-0.44807362, -0.89399666,  0.],
                           #           [ 0.89399666, -0.44807362,  0.],
                           #           [ 0.,0.,  1.]]),
                           's3off':0.0}
                 })

        if self.day == 12:
            #BS idem HKB:1 and HKB:2
            if ('BS'  in self.typ) or ('FULL' in self.typ):
                self.din.update(
                {'BS:74':{'p':mphkb[3],
                          # 'T':np.eye(3),
                          's3off':-0.2},
                 'BS:157':{'p':mphkb[2],
                          # 'T': np.array([[-0.44807362,  0.89399666,  0.],
                          #                [-0.89399666, -0.44807362,  0.],
                          #                [ 0.,0.,1.        ]]),
                          's3off':-0.2}      ,
                 })

        #load extra  information from inifile (antenna, rotation matrix,...)

        inifile = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','BodyandWear','AccesPoints.ini')
        config = ConfigParser.ConfigParser()
        config.read(inifile)

        for d in self.din:
            self.din[d]['antname']=config.get(d,'file')
            self.din[d]['ant']=antenna.Antenna(config.get(d,'file'))
            self.din[d]['T']=eval(config.get(d,'t'))
            self.din[d]['comment']=config.get(d,'comment')
        # self.pts= np.empty((12,3))
        # self.pts[:,0]= -mpts[:,1]
        # self.pts[:,1]= mpts[:,0]
        # self.pts[:,2]= mpts[:,2]
        # return mpts
        # self.dist = np.sqrt(np.sum((mpts[:,np.newaxis,:]-mpts[np.newaxis,:])**2,axis=2))


    def loadlog(self):
        """ load in self.log the log of the current serie
            from MeasurementLog.csv
        """

        filelog = os.path.join(self.rootdir,'RAW','Doc','MeasurementLog.csv')
        log = pd.read_csv(filelog)
        date = str(self.day)+'/06/14'
        self.log = log[(log['Meas Serie'] == self.serie) & (log['Date'] == date)]


    def _loadbody(self,day=11,serie=''):
        """ load motion capture file

        Parameters
        ----------

        day :
        serie :

        """
        assert day in [11,12],"wrong day in _loadbody"
        self.B={}
        color=['LightBlue','YellowGreen','PaleVioletRed','white','white','white','white','white','white','white']

        for us,subject in enumerate(self.subject):
            print "\nload ",subject, " body:",
            seriestr = str(self.serie).zfill(3)
            if day == 11:
                self.filemocap = os.path.join(self.rootdir,'RAW',str(self.day)+'-06-2014','MOCAP','serie_'+seriestr+'.c3d')
            elif day == 12:
                self.filemocap = os.path.join(self.rootdir,'RAW',str(self.day)+'-06-2014','MOCAP','Nav_serie_'+seriestr+'.c3d')

            # body and wear directory
            baw = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','BodyandWear')
            if subject =='Jihad':
                subject ='Jihan'


            #
            # Load body cylinder description : "Subject.ini"
            # Load wearable deviace description (contains antenna filename) :
            #
            self.filebody = os.path.join(baw, subject + '.ini')
            self.filewear = os.path.join(baw,subject + '_'  +str(self.day)+'-06-2014_' + self.typ + '.ini')

            if len(self.subject) >1 or self.mocapinterf:
                multi_subject=True
            else:
                multi_subject=False
            self.B.update({subject:Body(_filebody=self.filebody,
                             _filemocap=self.filemocap,unit = 'mm', loop=False,
                             _filewear=self.filewear,
                             centered=False,
                             multi_subject_mocap=multi_subject,
                             color=color[us])})

        if self.serie in self.mocapinterf:
            self.interf = ['Anis_Cylindre:',
                     'Benoit_Cylindre:',
                     'Bernard_Cylindre:',
                     'Claude_Cylindre:',
                     'Meriem_Cylindre:']
            intertmp=[]
            for ui,i in enumerate(self.interf):
                try:
                    print "load ",i, " interfering body:",

                    self.B.update({i:Cylinder(name=i,
                                              _filemocap=self.filemocap,
                                              unit = 'mm',
                                              color = color[ui])})
                    intertmp.append(i)
                except:
                    print "Warning ! load ",i, " FAIL !"
            self.interf=intertmp
        else :
            self.interf=[]
        # if len(self.subject) == 1:
        #     self.B = self.B[self.subject]


    def _loadTCR(self,day=11,serie='',scenario='20',run=1):
        """ load TCR data

        Parameters
        ----------

        day :
        serie :
        scenario :
        run :

        """

        #
        # TNET : (NodeId,MAC)
        #

        self.TNET={0:31,
        1:2,
        7:24,
        8:25,
        9:26,
        10:27,
        11:28,
        12:30,
        14:32,
        15:33,
        16:34,
        17:35,
        18:36,
        19:37,
        20:48,
        21:49}

        if day==11:
            self.dTCR ={'Unused':49,
                  'COORD':31,
                  'AP1':32,
                  'AP2':24,
                  'AP3':27,
                  'AP4':28,
                  'HeadRight':34,
                  'TorsoTopRight':25,
                  'TorsoTopLeft':30,
                  'BackCenter':35,
                  'HipRight':2,
                  'WristRight':26,
                  'WristLeft':48,
                  'KneeLeft':33,
                  'AnkleRight':36,
                  'AnkleLeft':37}
            dirname = os.path.join(self.rootdir,'POST-TREATED','11-06-2014','TCR')


        if day==12:
            dirname = os.path.join(self.rootdir,'POST-TREATED','12-06-2014','TCR')
            self.dTCR ={ 'COORD':31,
                        'AP1':32,
                        'AP2':24,
                        'AP3':27,
                        'AP4':28,
                   'Jihad:TorsoTopRight':35,
                   'Jihad:TorsoTopLeft':2,
                   'Jihad:BackCenter':33,
                   'Jihad:ShoulderLeft':37,
                   'Nicolas:TorsoTopRight':34,
                   'Nicolas:TorsoTopLeft':49,
                   'Nicolas:BackCenter':48,
                   'Nicolas:ShoulderLeft':36,
                   'Eric:TorsoCenter':30,
                   'Eric:BackCenter':25,
                   'Eric:ShoulderLeft':26}

        #
        # TCR  : (Name , MAC)
        # iTCR : (MAC , Name)
        # dTCR : (NodeId, Name)
        #

        self.idTCR={}
        for k in self.dTCR:
            self.idTCR[self.dTCR[k]]=k


        dTCRni={}
        for k in self.TNET.keys():
            dTCRni[k]=self.idTCR[self.TNET[k]]


        files = os.listdir(dirname)
        if serie != '':
            try:
                self._fileTCR = filter(lambda x : '_S'+str(serie)+'_' in x ,files)[0]
            except:
                self._fileTCR = filter(lambda x : '_s'+str(serie)+'_' in x ,files)[0]
            tt = self._fileTCR.split('_')
            self.scenario=tt[0].replace('Sc','')
            self.run = tt[2].replace('R','')
            self.typ = tt[3].replace('.csv','').upper()
            self.video = 'NA'
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            self._fileTCR = filter(lambda x : 'R'+str(run) in x ,filsc)[0]
            self.scenario= scenario
            self.run = str(run)

        filename = os.path.join(dirname,self._fileTCR)
        dtTCR = pd.read_csv(filename)
        tcr={}
        for k in dTCRni:
            for l in dTCRni:
                if k!=l:
                    d = dtTCR[((dtTCR['ida']==k) & (dtTCR['idb']==l))]
                    d.drop_duplicates('time',inplace=True)
                    del d['lqi']
                    del d['ida']
                    del d['idb']
                    d = d[d['time']!=-1]
                    d.index = d['time']
                    del d['time']
                    if len(d)!=0:
                        sr = pd.Series(d['dist']/1000,index=d.index)
                        tcr[dTCRni[k]+'-'+dTCRni[l]]= sr


        self.tcr = pd.DataFrame(tcr)
        self.tcr = self.tcr.fillna(0)
        ts = 75366400./1e9
        t = np.array(self.tcr.index)*ts
        t = t-t[0]
        self.tcr.index = t
        self.ttcr=self.tcr.index

    def _loadBS(self,day=11,serie='',scenario='20',run=1):
        """ load BeSpoon data

        Parameters
        ----------

        day : int
        serie : string
        scenario : string
        run : int

        """

        if day == 11:
            self.dBS = {'WristRight':157,'AnkleRight':74,'HandRight':0}
        elif day == 12:
            self.dBS = {'AP1':157,'AP2':74,'HandRight':0}
        self.idBS={}
        for k in self.dBS:
            self.idBS[self.dBS[k]]=k

        if day==11:
            dirname = os.path.join(self.rootdir,'POST-TREATED','11-06-2014','BeSpoon')
        if day==12:
            dirname = os.path.join(self.rootdir,'POST-TREATED','12-06-2014','BeSpoon')

        files = os.listdir(dirname)
        if serie != '':
            self._fileBS = filter(lambda x : 'S'+str(serie) in x ,files)[0]
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            self._fileBS = filter(lambda x : 'R'+str(run) in x ,filsc)[0]

        bespo = pd.read_csv(os.path.join(dirname,self._fileBS),index_col='ts')


        gb = bespo.groupby(['Sensor'])
        #get device id
        devid,idevid = np.unique(bespo['Sensor'],return_index=True)
        # get index of each group
        dgb={d:gb.get_group(d) for d in devid}
        lgb=[]
        for i in dgb:
            ind = dgb[i].index/1e3
            dti = pd.to_datetime(ind,unit='s')
            npai = time2npa(dti)
            npai = npai - npai[0]
            dgb[i].index=pd.Index(npai)
            lgb.append(pd.DataFrame(dgb[i]['d'].values,columns=[self.idBS[0]+'-'+self.idBS[i]],index=dgb[i].index))
        df = lgb[0].join(lgb[1])
        self.bespo = df


        #self.s157 = self.bespo[self.bespo['Sensor']==157]
        #self.s157.set_index(self.s157['tu'].values/1e9)
        #self.s74  = self.bespo[self.bespo['Sensor']==74]
        #self.s74.set_index(self.s74['tu'].values/1e9)
        #t157 = np.array(self.s157['tu']/(1e9))
        #self.t157 = t157-t157[0]
        #t74 = np.array(self.s74['tu']/(1e9))
        #self.t74 = t74 - t74[0]


    def _loadhkb(self,day=11,serie='',scenario='20',run=1,source='CITI'):
        """ load hkb measurement data

        Parameters
        ----------

        day : string
        serie : string
        scenario : string
        run : int
        source : 'string'

        Returns
        -------

        update self.hkb

        """

        if day == 11:
            if serie == 5:
                source = 'UR1'

        if day==11:
            self.dHKB ={'AP1':1,'AP2':2,'AP3':3,'AP4':4,
                       'HeadRight':5,'TorsoTopRight':6,'TorsoTopLeft':7,'BackCenter':8,'ElbowRight':9,'ElbowLeft':10,'HipRight':11,'WristRight':12,'WristLeft':13,'KneeLeft':14,'AnkleRight':16,'AnkleLeft':15}
            if source=='UR1' :
                dirname = os.path.join(self.rootdir,'POST-TREATED','11-06-2014','HIKOB')
            elif source=='CITI':
                dirname = os.path.join(self.rootdir,'POST-TREATED','11-06-2014','HIKOB','CITI')
        if day==12:
            self.dHKB= {'AP1':1,'AP2':2,'AP3':3,'AP4':4,'Jihad:TorsoTopRight':10,'Jihad:TorsoTopLeft':9,'Jihad:BackCenter':11,'JihadShoulderLeft':12,
             'Nicolas:TorsoTopRight':6,'Nicolas:TorsoTopLeft':5,'Nicolas:BackCenter':7,'Nicolas:ShoulderLeft':8,
             'Eric:TooTopRight':15,'Eric:TorsoTopLeft':13,'Eric:BackCenter':16,'Eric:ShoulderLeft':14}
            #if source=='UR1':
            dirname = os.path.join(self.rootdir,'POST-TREATED','12-06-2014','HIKOB')

        files = os.listdir(dirname)

        self.idHKB={}
        for k in self.dHKB:
            self.idHKB[self.dHKB[k]]=k

        if serie != '':
            self._filehkb = filter(lambda x : 'S'+str(serie) in x ,files)[0]
            tt = self._filehkb.split('_')
            if source == 'UR1':
                self.scenario=tt[0].replace('Sc','')
                self.run = tt[2].replace('R','')
                self.typ = tt[3]
                self.video = tt[4].replace('.mat','')
            elif source == 'CITI':
                self.scenario=tt[0].replace('Sc','')+tt[1]

                self.run = tt[3].replace('r','')
                self.typ = tt[4]
                if self.typ == 'HKB':
                    self.typ = 'HKBS'
                self.video = tt[5].replace('.mat','')

        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            if source=='UR1':
                self._filehkb = filter(lambda x : 'R'+str(run) in x ,filsc)[0]
            else:
                self._filehkb = filter(lambda x : 'r'+str(run) in x ,filsc)[0]


        data = io.loadmat(os.path.join(dirname,self._filehkb))
        if source=='UR1':
            self.rssi = data['rssi']
            self.thkb = data['t']
        else:
            self.rssi = data['val']
            self.thkb = np.arange(np.shape(self.rssi)[2])*25.832e-3

        def topandas():
            try:
                self.hkb = pd.DataFrame(index=self.thkb[0])
            except:
                self.hkb = pd.DataFrame(index=self.thkb)
            for k in self.idHKB:
                for l in self.idHKB:
                    if k!=l:
                        col  = self.idHKB[k]+'-'+self.idHKB[l]
                        rcol = self.idHKB[l]+'-'+self.idHKB[k]
                        if rcol not in self.hkb.columns:
                            rssi  = self.rssi[k-1,l-1,:]
                            self.hkb[col] = rssi

        topandas()
        self.hkb = self.hkb[self.hkb!=0]

    def compute_visibility(self,techno='HKB',square_mda=True,all_links=True):
        """ determine visibility of links for a given techno


            Parameters
            ----------

            techno  string
                select the given radio technology of the nodes to determine
                    the visibility matrix

            square_mda  boolean
                select ouput format
                    True : (device x device x timestamp)
                    False : (link x timestamp)

            all_links : bool
                compute all links or just those for which data is available

            Return
            ------

            if square_mda = True

            intersection : (ndevice x nbdevice x nb_timestamp)
                matrice of intersection (1 if link is cut 0 otherwise)
            links : (nbdevice)
                name of the links


            if square_mda = False

            intersection : (nblink x nb_timestamp)
                matrice of intersection (1 if link is cut 0 otherwise)
            links : (nblink x2)
                name of the links

            Example
            -------

            >>> from pylayers.measures.cormoran import *
            >>> import matplotlib.pyplot as plt
            >>> C=CorSer(serie=14,day=12)
            >>> inter,links=C.compute_visibility(techno='TCR',square_mda=True)
            >>> inter.shape
                (15, 15, 12473)
            >>>C.imshowvisibility_i(inter,links)


        """



        if techno == 'TCR':
            if not ((self.typ == 'TCR') or  (self.typ == 'FULL')):
                raise AttributeError('Serie has not data for techno: ',techno)
            hname = self.tcr.keys()
            dnode=copy.copy(self.dTCR)
            dnode.pop('COORD')
            prefix = 'TCR:'
        elif techno=='HKB':
            if not ((self.typ == 'HKBS') or  (self.typ == 'FULL')):
                raise AttributeError('Serie has not data for techno: '+techno)
            hname = self.hkb.keys()
            dnode=self.dHKB
            prefix = 'HKB:'
        # get link list
        if all_links:
            import itertools
            links =[l for l in itertools.combinations(dnode.keys(),2)]
        else:
            links=[n.split('-') for n in hname]
            links = [l for l in links if ('COORD' not in l[0]) and ('COORD' not in l[1])]
        #mapping between device name in self.hkb and on body/in self.devdf
        dev_bid = [self.devmapper(k,techno=techno)[2] for k in dnode.keys()]

        nb_totaldev=len(np.unique(self.devdf['id']))
        # extract all dev position on body
        # Mpdev : (3 x (nb devices and nb infra nodes) x nb_timestamp)
        Mpdev = np.empty((3,len(dev_bid),len(self.devdf.index)/nb_totaldev))

        # get all positions
        for ik,i in enumerate(dev_bid) :
            if i in self.din:
                Mpdev[:,ik,:] = self.din[i]['p'][:,np.newaxis]
            else:
                pts = self.devdf[self.devdf['id']==i][['x','y','z']].values.T
                if np.prod(pts.shape)!=0:
                    Mpdev[:,ik,:] = pts

        # create A and B from links
        nA = np.array([prefix+ str(dnode[l[0]]) for l in links])
        nB = np.array([prefix+ str(dnode[l[1]]) for l in links])

        dma = dict(zip(dev_bid,range(len(dev_bid))))
        mnA = [dma[n] for n in nA]
        mnB = [dma[n] for n in nB]

        A=Mpdev[:,mnA]
        B=Mpdev[:,mnB]


        # intersect2D matrix is
        # d_0: nb links
        #d_1: (cylinder number) * nb body + 1 * nb  cylinder_object
        # d_2 : nb frame
        intersect2D = np.zeros((len(links),
                                11*len(self.subject) + len(self.interf),
                                Mpdev.shape[-1]))
        # usub : index axes subject
        usub_start=0
        usub_stop=0
        # C-D correspond to bodies segments
        #C or D : 3 x 11 body segments x time
        # radius of cylinders are (nb_cylinder x time)
        for b in self.B:
            print 'processing shadowing from ',b
            # if b is a body not a cylinder
            if not 'Cylindre' in b:
                uta = self.B[b].sl[:,0].astype('int')
                uhe = self.B[b].sl[:,1].astype('int')
                rad = self.B[b].sl[:,2]

                C = self.B[b].d[:,uta,:]
                D = self.B[b].d[:,uhe,:]
                try:
                    radius = np.concatenate((radius,rad[:,np.newaxis]*np.ones((1,C.shape[2]))),axis=0)
                except:
                    radius = rad[:,np.newaxis]*np.ones((1,C.shape[2]))
                usub_start=usub_stop
                usub_stop=usub_stop+11
            else:

                cyl = self.B[b]
                # top of cylinder
                top = cyl.d[:,cyl.topnode,:]
                # bottom of cylinder =top with z =0
                bottom = copy.copy(cyl.d[:,cyl.topnode,:])
                bottom[2,:]=0.02
                #top 3 x 1 X time
                C=top[:,np.newaxis,:]
                D=bottom[:,np.newaxis,:]
                radius = np.concatenate((radius,cyl.radius[np.newaxis]))
                usub_start=usub_stop
                usub_stop=usub_stop+1

            f,g,X,Y,alpha,beta,dmin=seg.segdist(A,B,C,D,hard=True)

            intersect2D[:,usub_start:usub_stop,:]=g
            # import ipdb
            # ipdb.set_trace()
            #USEFUL Lines for debug
            #########################

            # def plt3d(ndev=53,ncyl=0,kl=11499):
            #     fig=plt.figure()
            #     ax=fig.add_subplot(111,projection='3d')
            #     if not isinstance(kl,list):
            #         kl=[kl]

            #     for ktime in kl:
            #         ax.plot([A[0,ndev,ktime],B[0,ndev,ktime]],[A[1,ndev,ktime],B[1,ndev,ktime]],[A[2,ndev,ktime],B[2,ndev,ktime]])
            #         [ax.plot([C[0,k,ktime],D[0,k,ktime]],[C[1,k,ktime],D[1,k,ktime]],[C[2,k,ktime],D[2,k,ktime]],'k') for k in range(11) ]
            #         ax.plot([X[0,ndev,ncyl,ktime],Y[0,ndev,ncyl,ktime]],[X[1,ndev,ncyl,ktime],Y[1,ndev,ncyl,ktime]],[X[2,ndev,ncyl,ktime],Y[2,ndev,ncyl,ktime]])
            #     ax.auto_scale_xyz([-5, 5], [-5, 5], [0, 2])
            #     plt.show()
            # import ipdb
            # ipdb.set_trace()




        uinter1 = np.where((intersect2D<=(radius-0.01)))
        uinter0 = np.where((intersect2D>(radius-0.01)))
        # intersect2D_=copy.copy(intersect2D)

        intersect2D[uinter1[0],uinter1[1],uinter1[2]]=1
        intersect2D[uinter0[0],uinter0[1],uinter0[2]]=0
        # #integrate the effect of all bodies by summing on axis 1
        intersect = np.sum(intersect2D,axis=1)>0

        if square_mda:
            dev= np.unique(links)
            ddev = dict(zip(dev,range(len(dev))))
            lmap = np.array(map(lambda x: (ddev[x[0]],ddev[x[1]]),links))
            M = np.nan*np.ones((len(dev),len(dev),intersect.shape[-1]))
            for i in range(len(intersect)):
                id1 = lmap[i][0]
                id2 = lmap[i][1]
                M[id1,id2,:]=intersect[i,:]
                M[id2,id1,:]=intersect[i,:]
            intersect=M
            links = dev

        self._visilinks = links
        self._visiintersect = intersect

        return intersect,links

    def imshowvisibility(self,techno='HKB',t=0,**kwargs):
        """  imshow visibility mda


        Parameters
        ----------

        techno : (HKB|TCR)
        t : float
            time in second

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> import matplotlib.pyplot as plt
        >>> C=CorSer(serie=6,day=12)
        >>> inter,links=C.compute_visibility(techno='TCR',square_mda=True)
        >>> i,l=C.imshowvisibility_i(inter,links)

        See Also
        --------

        pylayers.measures.CorSer.compute_visibility()

        """
        defaults = { 'grid':True,
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if 'fig' not in kwargs:
           fig = plt.figure()
        else:
           fig = kwargs.pop('fig')

        if 'ax' not in kwargs:
            ax = fig.add_subplot(111)
        else:
            ax = kwargs.pop('ax')

        if not '_visiintersect' in dir(self):
            print 'Visibility computed only once'
            self.compute_visibility(techno=techno)


        links = self._visilinks
        inter = self._visiintersect



        kt=np.where(self.tmocap <= t)[0][-1]
        plt.xticks(np.arange(0, len(links), 1.0))
        plt.yticks(np.arange(0, len(links), 1.0))
        ax.set_xlim([-0.5,len(links)-0.5])
        ax.set_ylim([len(links)-0.5,-0.5])
        ax.xaxis.set_ticks_position('top') 
        xtickNames = plt.setp(ax, xticklabels=links)
        ytickNames = plt.setp(ax, yticklabels=links)
        plt.setp(xtickNames, rotation=90, fontsize=8)
        plt.setp(ytickNames, rotation=0, fontsize=8)
        ims=[]
        ax.imshow(inter[:,:,kt],interpolation='nearest')
        if kwargs['grid']:
            ax.grid()
        return fig,ax

    def _show3i(self,t=0,**kwargs):
        """ show3 interactive
        """

        fig =plt.figure(num='Jog',figsize=(5,1.5))

        #set time to -10 is a trick to make appear interferers cylinder
        #because __refreshshow3i only update the data of the cylinder.
        # if cylinder is not present in the first _show3, they are not displayed
        # later.
        time=self.B[self.subject[0]].time

        fId = np.where(time<= t)[0][-1]

        kwargs['bodytime']=[self.tmocap[-10]]
        kwargs['returnfig']=True
        kwargs['tagtraj']=False
        mayafig = self._show3(**kwargs)
        self.__refreshshow3i(fId)
        # ax.grid()


        # matplotlib Widgets 

        slax=plt.axes([0.1, 0.5, 0.8, 0.3])
        slax.set_title('t='+str(time[fId]),loc='left')
        sliderx = Slider(slax, "time", 0, len(time),
                        valinit=fId, color='#AAAAAA')

        def update_x(val):
            value = int(sliderx.val)
            self.__refreshshow3i(val)
            slax.set_title('t='+str(time[val]),loc='left')
            fig.canvas.draw_idle()
        sliderx.on_changed(update_x)


        def plus(event):
            sliderx.set_val(sliderx.val +1)
            fig.canvas.draw_idle()

        def minus(event):
            sliderx.set_val(sliderx.val -1)
            fig.canvas.draw_idle()

        def pplus(event):
            sliderx.set_val(sliderx.val +10)
            fig.canvas.draw_idle()


        def mminus(event):
            sliderx.set_val(sliderx.val -10)
            fig.canvas.draw_idle()

        #QUIT by pressing 'q'
        def press(event):
            if event.key == 'q':
                mlab.close(mayafig)
                plt.close(fig)
        fig.canvas.mpl_connect('key_press_event', press)


        #-1 frame axes
        axm = plt.axes([0.2, 0.05, 0.1, 0.15])
        bm = Button(axm, '-1')
        bm.on_clicked(minus)
                #+1 frame axes
        axp = plt.axes([0.7, 0.05, 0.1, 0.15])
        bp = Button(axp, '+1')
        bp.on_clicked(plus)

        #-10 frames axes
        axmm = plt.axes([0.1, 0.05, 0.1, 0.15])
        bmm = Button(axmm, '-10')
        bmm.on_clicked(mminus)


        #+10 frames axes
        axpp = plt.axes([0.8, 0.05, 0.1, 0.15])
        bpp = Button(axpp, '+10')
        bpp.on_clicked(pplus)
        plt.show()


    def __refreshshow3i(self,kt):
        """ show3 update for interactive mode
            USED in imshowvisibility_i
        """
        t=self.tmocap[kt]

        for ib,b in enumerate(self.B):
            self.B[b].settopos(t=t,cs=True)

            try:
                # body
                X=np.hstack((self.B[b]._pta,self.B[b]._phe))
                self.B[b]._mayapts.mlab_source.set(x=X[0,:], y=X[1,:], z=X[2,:])
                # device
                udev = [self.B[b].dev[i]['uc3d'][0] for i in self.B[b].dev]
                Xd=self.B[b]._f[kt,udev,:].T
                self.B[b]._mayadev.mlab_source.set(x=Xd[0,:], y=Xd[1,:], z=Xd[2,:])
                # name
                uupper = np.where(X[2]==X[2].max())[0]
                self.B[b]._mayaname.actors.pop()
                self.B[b]._mayaname = mlab.text3d(X[0,uupper][0],X[1,uupper][0],X[2,uupper][0],self.B[b].name,scale=0.05,color=(1,0,0))
                # s = np.hstack((cylrad,cylrad))
            except:
                # cylinder
                X=np.vstack((self.B[b].top,self.B[b].bottom))
                self.B[b]._mayapts.mlab_source.set(x=X[:,0], y=X[:,1], z=X[:,2])
                # name
                self.B[b]._mayaname.actors.pop()
                self.B[b]._mayaname = mlab.text3d(self.B[b].top[0],self.B[b].top[1],self.B[b].top[2],self.B[b].name,scale=0.05,color=(1,0,0))
                #vdict
                V = self.B[b].traj[['vx','vy','vz']].iloc[self.B[b].toposFrameId].values

                self.B[b]._mayavdic.mlab_source.set(x= self.B[b].top[0],y=self.B[b].top[1],z=self.B[b].top[2],u=V[ 0],v=V[ 1],w=V[ 2])


    def imshowvisibility_i(self,techno='HKB',t=0,**kwargs):
        """  imshow visibility mda interactive

        Parameters
        ----------

        inter : (nb link x nb link x timestamps)
        links : (nblinks)
        time : intial time (s)


        Example
        -------

        >>> from pylayers.measures.cormoran import *
        >>> import matplotlib.pyplot as plt
        >>> C=CorSer(serie=6,day=12)
        >>> inter,links=C.visimda(techno='TCR',square_mda=True)
        >>> i,l=C.imshowvisibility_i(inter,links)

        """


        # if in_ipynb():
        #     notebook = False #program launch in ipyhon notebook
        #     from IPython.html import widgets # Widget definitions
        #     from IPython.display import display, clear_output# Used to display widgets in the notebook
        # else :
        #     notebook = False

        if not '_visiintersect' in dir(self):
            print 'Visibility is computed only once, Please wait\n'
            self.compute_visibility(techno=techno)


        links = self._visilinks
        inter = self._visiintersect


        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.3)


        time=self.tmocap

        fId = np.where(time<=t)[0][-1]

        vertc = [(0,-10),(0,-10),(0,10),(0,-10)]
        poly = plt.Polygon(vertc)
        pp = ax.add_patch(poly)

        plt.xticks(np.arange(0, len(links), 1.0))
        plt.yticks(np.arange(0, len(links), 1.0))
        ax.set_xlim([-0.5,len(links)-0.5])
        ax.set_ylim([len(links)-0.5,-0.5])
        ax.xaxis.set_ticks_position('top') 
        xtickNames = plt.setp(ax, xticklabels=links)
        ytickNames = plt.setp(ax, yticklabels=links)
        plt.setp(xtickNames, rotation=90, fontsize=8)
        plt.setp(ytickNames, rotation=0, fontsize=8)
        ims=[]
        l=ax.imshow(inter[:,:,fId],interpolation='nearest')
        #set time to -10 is a trick to make appear interferers cylinder
        #because __refreshshow3i only update the data of the cylinder.
        # if cylinder is not present in the first _show3, they are not displayed
        # later.
        kwargs['bodytime']=[self.tmocap[-10]]
        kwargs['returnfig']=True
        kwargs['tagtraj']=False
        mayafig = self._show3(**kwargs)
        self.__refreshshow3i(fId)
        # ax.grid()


        # matplotlib Widgets 
        slax=plt.axes([0.1, 0.15, 0.8, 0.05])
        slax.set_title('t='+str(time[fId]),loc='left')
        sliderx = Slider(slax, "time", 0, inter.shape[-1],
                        valinit=fId, color='#AAAAAA')

        # else :
        #     int_range = widgets.IntSliderWidget(min=0,max=inter.shape[-1],step=1,value=fId)
        #     display(int_range)


        def update_x(val):
            value = int(sliderx.val)
            sliderx.valtext.set_text('{}'.format(value))
            l.set_data(inter[:,:,value])
            self.__refreshshow3i(val)
            slax.set_title('t='+str(time[val]),loc='left')
            fig.canvas.draw_idle()
        sliderx.on_changed(update_x)
        # else:
        #     def update_x(name,value):

        #         clear_output(wait=True)
        #         display(plt.gcf())
        #         plt.imshow(inter[:,:,value],interpolation='nearest')
        #         # l.set_data(inter[:,:,value])
        #         kwargs['bodytime']=[self.tmocap[value]]
        #         self._show3(**kwargs)
        #         myu.inotshow('fig1',width=200,height=200,magnification=1)
        #         # slax.set_title('t='+str(time[val]),loc='left')
        #         # fig.canvas.draw_idle()
        #     int_range.on_trait_change(update_x, 'value')



        def plus(event):
            sliderx.set_val(sliderx.val +1)
            fig.canvas.draw_idle()
        # if not notebook:
        sliderx.on_changed(update_x)


        def minus(event):
            sliderx.set_val(sliderx.val -1)
            fig.canvas.draw_idle()
        # if not notebook:
        sliderx.on_changed(update_x)

        def pplus(event):
            sliderx.set_val(sliderx.val +10)
            fig.canvas.draw_idle()
        # if not notebook:
        sliderx.on_changed(update_x)


        def mminus(event):
            sliderx.set_val(sliderx.val -10)
            fig.canvas.draw_idle()
        # if not notebook:
        sliderx.on_changed(update_x)

        # #QUIT by pressing 'q'
        # def press(event):
        #     if event.key == 'q':
        #         mlab.close(mayafig)
        #         plt.close(fig)
        # fig.canvas.mpl_connect('key_press_event', press)

        # if not notebook:
        #-1 frame axes
        axm = plt.axes([0.3, 0.05, 0.1, 0.075])
        bm = Button(axm, '-1')
        bm.on_clicked(minus)

        #+1 frame axes
        axp = plt.axes([0.7, 0.05, 0.1, 0.075])
        bp = Button(axp, '+1')
        bp.on_clicked(plus)

        #-10 frames axes
        axmm = plt.axes([0.1, 0.05, 0.1, 0.075])
        bmm = Button(axmm, '-10')
        bmm.on_clicked(mminus)

        #+10 frames axes
        axpp = plt.axes([0.9, 0.05, 0.1, 0.075])
        bpp = Button(axpp, '+10')
        bpp.on_clicked(pplus)



        plt.show()




    def _distancematrix(self):
        """Compute the ditance matrix between the nodes

            self.dist : (nb frame x nb_node x nb_node)
            self.dist_nodesmap : list of used nodes (useful to make the association ;) )
        """

        if not isinstance(self.B,dict):
            B={self.subject[0]:self.B}
        else :
            B=self.B


        bn= []

        for b in B:
            if 'dev' in dir(B[b]):
                tdev=[]
                for k in B[b].dev:
                    bn.append(k)
                    tdev.append(B[b].dev[k]['uc3d'][0])
                tdev=np.array(tdev)
                try:
                    pnb = np.concatenate((pnb,B[b]._f[:,tdev,:]),axis=1)
                except:
                    pnb = B[b]._f[:,tdev,:]
        ln = []
        uin = []

        # infrastructure nodes
        if ('HK' in self.typ) or ('FULL' in self.typ):
            uin.extend(['HKB:1','HKB:2','HKB:3','HKB:4'])
        if ('TCR' in self.typ) or ('FULL' in self.typ):
            # TCR:31 is the coordinator (1.7719,-3.26)
            uin.extend(['TCR:32','TCR:24','TCR:27','TCR:28','TCR:31'])
        if self.day == 12:
            if ('BS' in self.typ) or ('FULL' in self.typ):
                uin.extend(['BS:74','BS:157'])
        ln = uin + bn
        pin = np.array([self.din[d]['p'] for d in uin])
        pin2 = np.empty((pnb.shape[0],pin.shape[0],pin.shape[1]))
        pin2[:,:,:]=pin
        p = np.concatenate((pin2,pnb),axis=1)
        self.dist = np.sqrt(np.sum((p[:,:,np.newaxis,:]-p[:,np.newaxis,:,:])**2,axis=3))
        self.dist_nodesmap = ln


    def _computedistdf(self):
        """Compute the distance dataframe from distance matrix
        """
        
        # HIKOB
        if ('HK' in self.typ) or ('FULL' in self.typ):
            devmap = {self.devmapper(k,'hkb')[0]:self.devmapper(k,'hkb')[2] for k in self.dHKB}
            udev = np.array([[self.dist_nodesmap.index(devmap[k.split('-')[0]]),self.dist_nodesmap.index(devmap[k.split('-')[1]])] for k in self.hkb.keys()])
            iudev =np.array([(self.dist_nodesmap[u[0]]+'-'+self.dist_nodesmap[u[1]]) for u in udev])
            df = pd.DataFrame(self.dist[:,udev[:,0],udev[:,1]],columns=iudev,index=self.tmocap)

        # BE Spoon
        if ('BS' in self.typ) or ('FULL' in self.typ):
            devmap = {self.devmapper(k,'BS')[0]:self.devmapper(k,'BS')[2] for k in self.dBS}
            udev = np.array([[self.dist_nodesmap.index(devmap[k.split('-')[0]]),self.dist_nodesmap.index(devmap[k.split('-')[1]])] for k in self.bespo.keys()])
            iudev =np.array([(self.dist_nodesmap[u[0]]+'-'+self.dist_nodesmap[u[1]]) for u in udev])
            dfb = pd.DataFrame(self.dist[:,udev[:,0],udev[:,1]],columns=iudev,index=self.tmocap)
            df = df.join(dfb)
            del dfb

        if ('TCR' in self.typ) or ('FULL' in self.typ):
            devmap = {self.devmapper(k,'tcr')[0]:self.devmapper(k,'tcr')[2] for k in self.dTCR}
            udev = np.array([[self.dist_nodesmap.index(devmap[k.split('-')[0]]),
                              self.dist_nodesmap.index(devmap[k.split('-')[1]])]
                             for k in self.tcr.keys() ])
            #                 for k in self.tcr.keys() if not 'COORD' in k])
            iudev =np.array([(self.dist_nodesmap[u[0]]+'-'+self.dist_nodesmap[u[1]]) for u in udev])
            dft = pd.DataFrame(self.dist[:,udev[:,0],udev[:,1]],columns=iudev,index=self.tmocap)
            if ('FULL' in self.typ):
                df = df.join(dft)
            else :
                df = dft
            del dft

        self.distdf=df

    # def accessdm(self,a,b,techno=''):
    #     """ access to the distance matrix

    #         give name|id of node a and b and a given techno. retrun Groung truth
    #         distance between the 2 nodes
    #     # """

    #     # a,ia,bia,subja=self.devmapper(a,techno)
    #     # b,ib,bib,subjb=self.devmapper(b,techno)

    #     if 'HKB' in techno :
    #         if isinstance(a,str):
    #             ia = self.dHKB[a]
    #         else:
    #             ia = a
    #             a = self.idHKB[a]

    #         if isinstance(b,str):
    #             ib = self.dHKB[b]
    #         else:
    #             ib = b
    #             b = self.idHKB[b]

    #     elif 'TCR' in techno :
    #         if isinstance(a,str):
    #             ia = self.dTCR[a]
    #         else:
    #             ia = a
    #             a = self.idTCR[a]

    #         if isinstance(b,str):
    #             ib = self.dTCR[b]
    #         else:
    #             ib = b
    #             b = self.idTCR[b]

    #     else :
    #         raise AttributeError('please give only 1 techno or radio node')

    #     ka = techno+':'+str(ia)
    #     kb = techno+':'+str(ib)

    #     ua = self.dist_nodesmap.index(ka)
    #     ub = self.dist_nodesmap.index(kb)

    #     return(ua,ub)





        # c3ds = self.B._f.shape
        # if 'Full' in self.typ:
        #     pdev= np.empty((c3ds[0],len(self.dHKB)+len(self.tcr)+len(bs),3))
        # elif 'HK' in self.typ:
        #     pdev= np.empty((c3ds[0],len(self.dHKB)+len(bs),3))
        # elif 'TCR' in self.typ:
        #     pdev= np.empty((c3ds[0],len(self.tcr),3))
        # else:
        #     raise AttributeError('invalid self.typ')

        # self.B.network()
        # DB = self.B.D2

        # ludev = np.array([[i,self.B.dev[i]['uc3d'][0]] for i in self.B.dev])
        # for i in ludev:
        #     pdev[:,eval(i[0])-1,:] = self.B._f[:,i[1],:]
        # # self.dist = np.sqrt(np.sum((mpts[:,np.newaxis,:]-mpts[np.newaxis,:])**2,axis=2))


    def vlc(self):
        """ play video of the associated serie
        """
        videofile = os.path.join(self.rootdir,'POST-TREATED', str(self.day)+'-06-2014','Videos')
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)

        try:
            uldir = luldir.index(True)
            _filename = ldir[uldir]
            filename = os.path.join(videofile,_filename)
            os.system('vlc '+filename +'&' )
        except:
            raise AttributeError('file '+ self._filename + ' not found')


    def snapshot(self,t0=0,offset=15.5,title=True,save=False,fig=[],ax=[],figsize=(10,10)):
        """ single snapshot plot

        Parameters
        ----------

        t0: float
        offset : float
        title : boolean
        save : boolean
        fig
        ax
        figsize : tuple

        Examples
        --------


        """

        if fig ==[]:
            fig=plt.figure(figsize=figsize)
        if ax == []:
            ax = fig.add_subplot(111)

        if self.offset[self._filename].has_key('video_sec'):
            offset = self.offset[self._filename]['video_sec']
        elif offset != '':
            offset = offset
        else:
            offset=0

        videofile = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','Videos')
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = os.path.join(videofile,_filename)
        vc = VideoFileClip(filename)
        F0 = vc.get_frame(t0+offset)
        I0 = img_as_ubyte(F0)
        ax.imshow(F0)
        if title:
            ax.set_title('t = '+str(t0)+'s')
        if save :
            plt.savefig(self._filename +'_'+str(t0) + '_snap.png',format='png')

        return fig,ax


    def snapshots(self,t0=0,t1=10,offset=15.5):
        """ take snapshots 

        Parameters
        ----------

        t0 : float
        t1 : float

        """

        if self.offset[self._filename].has_key('video_sec'):
            offset = self.offset[self._filename]['video_sec']
        elif offset != '':
            offset = offset
        else:
            offset=0


        videofile = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','Videos')
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = os.path.join(videofile,_filename)
        vc = VideoFileClip(filename)
        F0 = vc.get_frame(t0+offset)
        F1 = vc.get_frame(t1+offset)
        I0 = img_as_ubyte(F0)
        I1 = img_as_ubyte(F1)
        plt.subplot(121)
        plt.imshow(F0)
        plt.title('t = '+str(t0)+'s')
        plt.subplot(122)
        plt.imshow(F1)
        plt.title('t = '+str(t1)+'s')


    def _show3(self,**kwargs):
        """ mayavi 3d show of scenario

        Parameters
        ----------

        L : boolean
            display layout (True)

        body :boolean
            display bodytime(True)
        bodyname : boolean
            display body name
        bodytime: list
            list of time instant where body topos has to be shown


        devsize : float
            device on body size (100)
        devlist : list
            list of device name to show on body
        pattern : boolean
            display devices pattern

        trajectory : boolean
            display trajectory  (True)
        tagtraj : boolean
            tag on trajectory at the 'bodytime' instants (True)
        tagname : list
            name of the tagtrajs
        tagpoffset : ndarray
            offset of the tag positions (nb_of_tags x 3)
        fontsizetag : float
            size of the tag names


        inodes : boolean
            display infrastructure nodes
        inname : boolean
            display infra strucutre node name
        innamesize : float,
            size of name of infrastructure nodes (0.1)
        incolor: str
            color of infrastructure nodes ('r')
        insize
            size of infrastructure nodes (0.1)


        camera : boolean
            display Vicon camera position (True)
        cameracolor : str
            color of camera nodes ('b')
        camerasize  : float
            size of camera nodes (0.1)

        Examples
        --------

            >>> S  = Corser(6)
            >>> S._show3()



        """
        defaults = { 'L':True,
                     'body':True,
                     'bodyname':True,
                     'subject':[],
                     'interf':True,
                     'trajectory' :False,
                     'devsize':100,
                     'devlist':[],
                     'pattern':False,
                     'inodes' : True,
                     'inname' : True,
                     'innamesize' : 0.1,
                     'incolor' : 'r',
                     'insize' : 0.1,
                     'camera':True,
                     'cameracolor' :'k',
                     'camerasize' :0.1,
                     'bodytime':[],
                     'tagtraj':True,
                     'tagname':[],
                     'tagpoffset':[],
                     'fontsizetag':0.1,
                     'trajectory_color_range':True
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        cold = pyu.coldict()
        camhex = cold[kwargs['cameracolor']]
        cam_color = tuple(pyu.rgb(camhex)/255.)
        inhex = cold[kwargs['incolor']]
        in_color = tuple(pyu.rgb(inhex)/255.)

        if kwargs['subject'] == []:
            subject = self.subject
        else:
            subject = kwargs['subject']

        if kwargs['L']:
            self.L._show3(opacity=0.5)
        v = self.din.items()
        if kwargs['inodes']:
            X= np.array([v[i][1]['p'] for i in range(len(v))])
            mlab.points3d(X[:,0],X[:,1], X[:,2],scale_factor=kwargs['insize'],color=in_color)
            if kwargs['pattern']:
                for i in range(len(v)):
                    if not hasattr(self.din[v[i][0]]['ant'],'SqG'):
                        self.din[v[i][0]]['ant'].Fsynth()
                    self.din[v[i][0]]['ant']._show3(po=v[i][1]['p'],
                           T=self.din[v[i][0]]['T'],
                           ilog=False,
                           minr=0.01,
                           maxr=0.2,
                           newfig=False,
                           title=False,
                           colorbar=False,
                           )

        if kwargs['inname']:
            [mlab.text3d(v[i][1]['p'][0],
                        v[i][1]['p'][1],
                        v[i][1]['p'][2]+v[i][1]['s3off'],
                        v[i][0],
                        scale=kwargs['innamesize']) for i in range(len(v))]
        if kwargs['body']:

            if kwargs['bodytime']==[]:
                time =np.linspace(0,self.B[subject[0]].time[-1],5).astype(int)
                # time=range(10,100,20)
            else :
                time=kwargs['bodytime']

            for ki, i in enumerate(time):
                for ib,b in enumerate(subject):
                    self.B[b].settopos(t=i,cs=True)
                    self.B[b]._show3(dev=True,
                                    name = kwargs['bodyname'],
                                    devlist=kwargs['devlist'],
                                    devsize=kwargs['devsize'],
                                    tube_sides=12,
                                    pattern=kwargs['pattern'])
                    if kwargs['tagtraj']:
                        X=self.B[b].traj[['x','y','z']].values[self.B[b].toposFrameId]
                        if kwargs['tagpoffset']==[]:
                            X[2]=X[2]+0.2
                        else :
                            X=X+kwargs['tagpoffset'][ki]
                        if kwargs['tagname']==[]:
                            name = 't='+str(i)+'s'
                        else :
                            name = str(kwargs['tagname'][ki])
                        mlab.text3d(X[0],X[1],X[2],name,scale=kwargs['fontsizetag'])
                if kwargs['interf']:
                    for ib,b in enumerate(self.interf):
                        self.B[b].settopos(t=i,cs=True)
                        self.B[b]._show3(name=kwargs['bodyname'],tube_sides=12)

        if kwargs['trajectory']:
            for b in subject:
                self.B[b].traj._show3(kwargs['trajectory_color_range'])
        if kwargs['camera'] :
            mlab.points3d(self.cam[:,0],self.cam[:,1], self.cam[:,2],scale_factor=kwargs['camerasize'],color=cam_color)
        mlab.view(-111.44127634143871,
                    60.40674368088245,
                    24.492297713984197,
                    array([-0.07235499,  0.04868631, -0.00314969]))
        mlab.view(-128.66519195313163,
                   50.708933839573511,
                   24.492297713984247,
                   np.array([-0.07235499,  0.04868631, -0.00314969]))

    def anim(self):

        self._show3(body=False,inname=False,trajectory=False)
        [self.B[b].anim() for b in self.B]

        mlab.view(-43.413544538477254,
                    74.048193730704611,
                    11.425837641867618,
                    array([ 0.48298163,  0.67806043,  0.0987967 ]))


    def imshow(self,time=100,kind='time'):
        """ DEPRECATED

        Parameters
        ----------

        kind : string

            'mean','std'
        """
        fig = plt.figure(figsize=(10,10))
        self.D = self.rssi-self.rssi.swapaxes(0,1)

        try:
            timeindex = np.where(self.thkb[0]-time>0)[0][0]
        except:
            timeindex = np.where(self.thkb-time>0)[0][0]
        if kind=='time':
            dt1 = self.rssi[:,:,timeindex]
            dt2 = self.D[:,:,timeindex]

        if kind == 'mean':
            dt1 = ma.masked_invalid(self.rssi).mean(axis=2)
            dt2 = ma.masked_invalid(self.D).mean(axis=2)

        if kind == 'std':
            dt1 = ma.masked_invalid(self.rssi).std(axis=2)
            dt2 = ma.masked_invalid(self.D).std(axis=2)

        ax1 = fig.add_subplot(121)
        #img1 = ax1.imshow(self.rssi[:,:,timeindex],interpolation='nearest',origin='lower')
        img1 = ax1.imshow(dt1,interpolation='nearest')
        labels = map(lambda x : self.idHKB[x],range(1,17))
        plt.xticks(range(16),labels,rotation=80,fontsize=14)
        plt.yticks(range(16),labels,fontsize=14)
        if kind=='time':
            plt.title('t = '+str(time)+ ' s')
        if kind=='mean':
            plt.title(u'$mean(\mathbf{L})$')
        if kind=='std':
            plt.title(u'$std(\mathbf{L})$')
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        clb1 = fig.colorbar(img1,cax1)
        clb1.set_label('level dBm',fontsize=14)
        ax2 = fig.add_subplot(122)
        #img2 = ax2.imshow(self.D[:,:,timeindex],interpolation='nearest',origin='lower')
        img2 = ax2.imshow(dt2,interpolation='nearest')
        plt.title(u'$\mathbf{L}-\mathbf{L}^T$')
        divider = make_axes_locatable(ax2)
        plt.xticks(range(16),labels,rotation=80,fontsize=14)
        plt.yticks(range(16),labels,fontsize=14)
        cax2 = divider.append_axes("right", size="5%", pad=0.05)
        clb2 = fig.colorbar(img2,cax2)
        clb2.set_label('level dBm',fontsize=14)
        plt.tight_layout()
        plt.show()
        #for k in range(1,17):
        #    for l in range(1,17):
        #        self.dHKB[(k,l)]=iHKB[k]+' - '+iHKB[l]
        #        cpt = cpt + 1
        return fig,(ax1,ax2)

    def lk2nd(self,lk):
        """ transcode a lk from Id to real name

        Parameters
        ----------

        lk : string

        Examples
        --------

            >>> C=Corser(6)
            >>> lk = 'HKB:15-HKB:7'
            >>> C.lk2nd(lk)
        """
        u = lk.replace('HKB:','').split('-')
        v = map(lambda x : self.idHKB[int(x)],u)
        return(v)

    def _load_offset_dict(self):
        """ load offset_dictionnary.bin

        Returns
        -------

        d : dict

        {'Sc20_S5_R1_HKBS': {'hkb_index': -148, 'video_sec': 32.622087273809527},
        'Sc20_S6_R2_HKBS': {'bs_index': -124, 'hkb_index': -157},
        'Sc21a_S13_R1_HKBS': {'hkb_index': 537},
        'Sc21a_S14_R2_HKBS': {'hkb_index': 752},
        'Sc21a_S15_R3_HKBS': {'hkb_index': 438},
        'Sc21a_S16_R4_HKBS': {'hkb_index': 224},
        'Sc21b_S21_R1_HKBS': {'hkb_index': 368},
        'Sc21b_S22_R2_HKBS': {'hkb_index': -333},
        'Sc21b_S23_R3_HKBS': {'hkb_index': 136},
        'Sc22a_S9_R1_Full': {'hkb_index': 678}}

        Notes
        -----

        This is used for synchronization purpose

        """

        path = os.path.join(os.environ['CORMORAN'],'POST-TREATED')
        d = pickle.load( open( os.path.join(path,'offset_dictionnary.bin'), "rb" ) )

        return d

    def _save_offset_dict(self,d):

        path = os.path.join(os.environ['CORMORAN'],'POST-TREATED')
        d = pickle.dump( d, open( os.path.join(path,'offset_dictionnary.bin'), "wb" ) )

    def _save_data_off_dict(self,filename,typ,value):
        """ save
                - a given "value" of an for,
                - a serie/run "filename",
                - of a given typ (video|hkb|tcr|...)
        """

        d = self._load_offset_dict()
        try:
            d[filename].update({typ:value})
        except:
            d[filename]={}
            d[filename][typ]=value
        self._save_offset_dict(d)


    def offset_setter_video(self,a='AP1',b='WristRight',**kwargs):
        """ video offset setter
        """
        defaults = { 'inverse':True
                    }


        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        fig, axs = plt.subplots(nrows=2,ncols=1)
        fig.subplots_adjust(bottom=0.3)

        if isinstance(a,str):
            ia = self.dHKB[a]
        else:
            ia = a
            a = self.idHKB[a]

        if isinstance(b,str):
            ib = self.dHKB[b]
        else:
            ib = bq
            b = self.idHKB[b]


        time = self.thkb
        if len(time) == 1:
            time=time[0]


        sab = self.hkb[a+'-'+b].values
        sabt = self.hkb[a+'-'+b].index
        hkb = axs[1].plot(sabt,sab,label = a+'-'+b)
        axs[1].legend()


        try : 
            init = self.offset[self._filename]['video_sec']
        except:
            init=time[0]


        videofile = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','Videos')
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = os.path.join(videofile,_filename)
        vc = VideoFileClip(filename)
        F0 = vc.get_frame(init)
        I0 = img_as_ubyte(F0)
        axs[0].imshow(F0)



        ########
        # slider
        ########
        slide_xoffset_ax = plt.axes([0.1, 0.15, 0.8, 0.05])
        sliderx = Slider(slide_xoffset_ax, "video offset", 0, self.hkb.index[-1],
                        valinit=init, color='#AAAAAA')


        # vertc = [(0,-10),(0,-10),(0,10),(0,-10)]
        # poly = plt.Polygon(vertc)
        # pp = axs[1].add_patch(poly)


        def update_x(val):
            F0 = vc.get_frame(val)
            I0 = img_as_ubyte(F0)
            axs[0].imshow(F0)
            fig.canvas.draw_idle()
        sliderx.on_changed(update_x)

        # def cursor(val):
        #     try :
        #         pp.remove()
        #     except:
        #         pass
        #     vertc = [(sabt[0]+val,min(sab)-10),(sabt[0]+val,min(sab)-10),(sabt[0]+val,max(sab)+10),(sabt[0]+val,max(sab)-10)]
        #     poly = plt.Polygon(vertc)
        #     pp = axs[1].add_patch(poly)
        # sliderx.on_changed(cursor)

        def plus(event):
            sliderx.set_val(sliderx.val +0.2)
            fig.canvas.draw_idle()
        sliderx.on_changed(update_x)


        def minus(event):
            sliderx.set_val(sliderx.val -0.2)
            fig.canvas.draw_idle()
        sliderx.on_changed(update_x)


        def setter(event):
            self._save_data_off_dict(self._filename,'video_sec',sliderx.val)
            self.offset= self._load_offset_dict()
        axp = plt.axes([0.3, 0.05, 0.1, 0.075])
        axset = plt.axes([0.5, 0.05, 0.1, 0.075])
        axm = plt.axes([0.7, 0.05, 0.1, 0.075])

        bp = Button(axp, '<-')
        bp.on_clicked(minus)

        bset = Button(axset, 'SET offs.')
        bset.on_clicked(setter)

        bm = Button(axm, '->')
        bm.on_clicked(plus)

        plt.show()


    def offset_setter(self,a='HKB:1',b='HKB:12',techno='',**kwargs):
        """ offset setter
        """
        defaults = { 'inverse':True
                    }


        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        if plt.isinteractive():
            interactive = True
            plt.ioff()
        else :
            interactive = False

        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.2, left=0.3)

        a,ia,bia,subja,techno=self.devmapper(a,techno)
        b,ib,bib,subjb,techno=self.devmapper(b,techno)


        time = self.tmocap
        if len(time.shape) == 2:
            time = time[0,:]


        try : 
            init = time[0]#self.offset[self._filename]['hkb_index']
        except:
            init=time[0]


        var = self.getlinkd(ia,ib,techno).values
        if kwargs['inverse']:
            var = 10*np.log10(1./(var)**2)
        gt = ax.plot(time,var)


        ab = self.getlink(ia,ib,techno)
        sab = ab.values
        sabt = ab.index.values
        technoval = ax.plot(sabt,sab)


        ########
        # slider
        ########
        slide_xoffset_ax = plt.axes([0.1, 0.15, 0.8, 0.02])
        sliderx = Slider(slide_xoffset_ax, techno + " offset", -(len(sabt)/16), (len(sabt)/16),
                        valinit=init, color='#AAAAAA')

        slide_yoffset_ax = plt.axes([0.1, 0.10, 0.8, 0.02])
        slidery = Slider(slide_yoffset_ax, "gt_yoff", -100, 0,
                        valinit=0, color='#AAAAAA')

        slide_alpha_ax = plt.axes([0.1, 0.05, 0.8, 0.02])
        slideralpha = Slider(slide_alpha_ax, "gt_alpha", 0, 60,
                        valinit=30, color='#AAAAAA')

        def update_x(val):
            value = int(sliderx.val)
            rtechnoval = np.roll(sab,value)
            sliderx.valtext.set_text('{}'.format(value))
            technoval[0].set_xdata(sabt)
            technoval[0].set_ydata(rtechnoval)
            fig.canvas.draw_idle()
        sliderx.on_changed(update_x)
        sliderx.drawon = False


        def update_y(val):
            yoff = slidery.val
            alpha = slideralpha.val
            gt[0].set_ydata(alpha*var + yoff)
            fig.canvas.draw_idle()
        #initpurpose
        update_y(5)
        slidery.on_changed(update_y)
        slideralpha.on_changed(update_y)


        def setter(event):
            value = int(sliderx.val)
            try :
                nval = self.offset[self._filename][techno.lower()+'_index'] + value
            except : 
                nval = value
            self._save_data_off_dict(self._filename,techno.lower()+'_index',nval)
            self.offset= self._load_offset_dict()
            ax.set_title('WARNING : Please Reload serie to Valide offset change',color='r',weight='bold')

        axset = plt.axes([0.0, 0.5, 0.2, 0.05])
        bset = Button(axset, 'SET ' +techno+' offs.')
        bset.on_clicked(setter)

        plt.show()
        if interactive :
            plt.ion()


    # def offset_setter_hkb(self,a='AP1',b='WristRight',**kwargs):
    #     """ offset setter
    #     """
    #     defaults = { 'inverse':True
    #                 }


    #     for k in defaults:
    #         if k not in kwargs:
    #             kwargs[k] = defaults[k]

    #     if plt.isinteractive():
    #         interactive = True
    #         plt.ioff()
    #     else :
    #         interactive = False

    #     fig, ax = plt.subplots()
    #     fig.subplots_adjust(bottom=0.2, left=0.3)

    #     a,ia,bia,subja,techno=self.devmapper(a,'HKB')
    #     b,ib,bib,subjb,techno=self.devmapper(b,'HKB')


    #     time = self.thkb
    #     if len(time.shape) == 2:
    #         time = time[0,:]


    #     try : 
    #         init = time[0]#self.offset[self._filename]['hkb_index']
    #     except:
    #         init=time[0]


    #     var = self.getlinkd(ia,ib,'HKB').values
    #     if kwargs['inverse']:
    #         var = 10*np.log10(1./(var)**2)
    #     gt = ax.plot(self.B[self.B.keys()[0]].time,var)

    #     sab = self.hkb[a+'-'+b].values
    #     sabt = self.hkb[a+'-'+b].index
    #     hkb = ax.plot(sabt,sab)


    #     ########
    #     # slider
    #     ########
    #     slide_xoffset_ax = plt.axes([0.1, 0.15, 0.8, 0.02])
    #     sliderx = Slider(slide_xoffset_ax, "hkb offset", -(len(sabt)/16), (len(sabt)/16),
    #                     valinit=init, color='#AAAAAA')

    #     slide_yoffset_ax = plt.axes([0.1, 0.10, 0.8, 0.02])
    #     slidery = Slider(slide_yoffset_ax, "gt_yoff", -100, 0,
    #                     valinit=0, color='#AAAAAA')

    #     slide_alpha_ax = plt.axes([0.1, 0.05, 0.8, 0.02])
    #     slideralpha = Slider(slide_alpha_ax, "gt_alpha", 0, 10,
    #                     valinit=5, color='#AAAAAA')

    #     def update_x(val):
    #         value = int(sliderx.val)
    #         rhkb = np.roll(sab,value)
    #         sliderx.valtext.set_text('{}'.format(value))
    #         hkb[0].set_xdata(sabt)
    #         hkb[0].set_ydata(rhkb)
    #         fig.canvas.draw_idle()
    #     sliderx.on_changed(update_x)
    #     sliderx.drawon = False


    #     def update_y(val):
    #         yoff = slidery.val
    #         alpha = slideralpha.val
    #         gt[0].set_ydata(alpha*var + yoff)
    #         fig.canvas.draw_idle()
    #     #initpurpose
    #     update_y(5)
    #     slidery.on_changed(update_y)
    #     slideralpha.on_changed(update_y)


    #     def setter(event):
    #         value = int(sliderx.val)
    #         try :
    #             nval = self.offset[self._filename]['hkb_index'] + value
    #         except : 
    #             nval = value
    #         self._save_data_off_dict(self._filename,'hkb_index',nval)
    #         self.offset= self._load_offset_dict()
    #         ax.set_title('WARNING : Please Reload serie to Valide offset change',color='r',weight='bold')
    #     axset = plt.axes([0.0, 0.5, 0.2, 0.05])

    #     bset = Button(axset, 'SET offs.')
    #     bset.on_clicked(setter)


    #     plt.show()
    #     if interactive:
    #         plt.ion()


    def mtlbsave(self):
        """ Matlab format save

        """
        key = 'S'+str(self.day)+'_'+str(self.serie)
        filemat = key+'.mat'
        d = {}
        d[key]={}
        if ('HKB' in self.typ.upper()) or ('FULL' in self.typ.upper()):
            d[key]['HKB']={}
            links = list(self.hkb.columns)
            inter,lks = self.compute_visibility(techno='HKB')
            for l in links:
                ls   = l.split('-')
                nl = ls[0]+'_'+ls[1]
                d[key][nl] = {}
                ix0 = np.where(lks==ls[0])[0]
                ix1 = np.where(lks==ls[1])[0]
                Ssh = inter[ix0,ix1,:]
                Srssi= self.getlink(ls[0],ls[1],techno='HKB')
                # get distances between nodes
                Sdist = self.getlinkd(ls[0],ls[1],techno='HKB')
                dsh = rssi2dist(Sdist,Ssh,15)
                # rssi
                d[key]['HKB'][nl]['rssi'] = Srssi.values
                # dsh
                d[key]['HKB'][nl]['dsh'] = dsh
                #d['S6'][nl]['rssi_dec'] = np.roll(Srssi.values,-dec)
                d[key]['HKB'][nl]['sh'] = Ssh
                # time rssi
                d[key]['HKB'][nl]['tr'] = np.array(Srssi.index)
                # distance
                d[key]['HKB'][nl]['dist'] = Sdist.values
                # time mocap
                d[key]['HKB'][nl]['td'] = np.array(Sdist.index)

        if ('TCR' in self.typ.upper()) or ('FULL' in self.typ.upper()):
            d[key]['TCR']={}
            links = list(self.tcr.columns)
            inter,lks = self.compute_visibility(techno='TCR')
            for l in links:
                ls   = l.split('-')
                nl = ls[0]+'_'+ls[1]
                d[key]['TCR'][nl] = {}
                ix0 = np.where(lks==ls[0])[0]
                ix1 = np.where(lks==ls[1])[0]
                # intersection on the link
                Ssh = inter[ix0,ix1,:]
                Srange= self.getlink(ls[0],ls[1],techno='TCR')
                # get distances between nodes
                Sdist = self.getlinkd(ls[0],ls[1],techno='TCR')
                # rssi
                d[key]['TCR'][nl]['range'] = Srange.values
                # dsh
                #d['S6'][nl]['rssi_dec'] = np.roll(Srssi.values,-dec)
                d[key]['TCR'][nl]['sh'] = Ssh
                # time rssi
                d[key]['TCR'][nl]['tr'] = np.array(Srange.index)
                # distance
                d[key]['TCR'][nl]['dist'] = Sdist.values
                # time mocap
                d[key]['TCR'][nl]['td'] = np.array(Sdist.index)

        io.savemat(filemat,d)





    def pltvisi(self,a,b,techno='',**kwargs):
        """ plot visibility between link a and b


        Attributes
        ----------
        color:
            fill color
        hatch:
            hatch type
        label_pos: ('top'|'bottom'|'')
            postion of the label
        label_pos_off: float
            offset of postion of the label
        label_mob: str
            prefix of label in mobility
        label_stat: str
            prefix of label static


        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S = CorSer(6)
        >>> f,ax = S.plthkb('AP1','TorsoTopLeft',techno='HKB')
        >>> f,ax = S.pltvisi('AP1','TorsoTopLeft',techno='HKB',fig=f,ax=ax)
        >>> f,ax = S.pltmob(fig=f,ax=ax)
        >>> plt.title('hatch = visibility / gray= mobility')
        >>> plt.show()
        """


        defaults = { 'fig':[],
                     'figsize':(10,10),
                     'ax':[],
                     'color':'',
                     'hatch':'//',
                     'label_pos':'',
                     'label_pos_off':5,
                     'label_vis':'V',
                     'label_hide':'H'
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else :
            fig=kwargs['fig']

        if kwargs['ax'] ==[]:
            ax = fig.add_subplot(111)
        else :
            ax = kwargs['ax']



        aa= ax.axis()

        a,ia,nna,subjecta,technoa = self.devmapper(a,techno)
        b,ib,nnb,subjectb,technob = self.devmapper(b,techno)
        vv,tv,tseg,itseg = self._visiarray(nna,nnb)
        # vv.any : it exist NLOS regions
        if vv.any():
            if kwargs['color']=='':
                fig,ax=plu.rectplot(tv,tseg,ylim=aa[2:],
                                    fill=False,
                                    hatch=kwargs['hatch'],
                                    fig=fig,ax=ax)

            else :
                fig,ax=plu.rectplot(tv,tseg,ylim=aa[2:],
                                    color=kwargs['color'],
                                    hatch=kwargs['hatch'],
                                    fig=fig,ax=ax)

            if kwargs['label_pos']!='':
                if kwargs['label_pos'] == 'top':
                    yposV = aa[3]-kwargs['label_pos_off']+0.5
                    yposH = aa[3]-kwargs['label_pos_off']-0.5

                elif kwargs['label_pos'] == 'bottom':
                    yposV = aa[2]+kwargs['label_pos_off']+0.5
                    yposH = aa[2]+kwargs['label_pos_off']+0.5
                xposV= tv[tseg.mean(axis=1).astype(int)]
                xposH= tv[itseg.mean(axis=1).astype(int)]
                [ax.text(x,yposV,kwargs['label_vis']+str(ix+1)) for ix,x in enumerate(xposV)]
                [ax.text(x,yposH,kwargs['label_hide']+str(ix+1)) for ix,x in enumerate(xposH)]

        return fig,ax

    def pltmob(self,**kwargs):
        """ plot mobility

        Parameters
        ----------
        subject: str
            subject to display () if '', take the fist one from self.subject)
        showvel :  boolean
            display filtered velocity
        velth: float (0.7)
            velocity threshold
        fo : int (5)
            filter order
        fw: float (0.02)
            0 < fw < 1  (fN <=> 1)
        time_offset : int
            add time_offset to start later

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S = CorSer(6)
        >>> f,ax = S.plthkb('AP1','TorsoTopLeft',techno='HKB')
        >>> #f,ax = S.pltvisi('AP1','TorsoTopLeft',techno='HKB',fig=f,ax=ax)
        >>> f,ax = S.pltmob(fig=f,ax=ax)
        >>> plt.title('hatch = visibility / gray= mobility')
        >>> plt.show()
        """
        defaults = { 'subject':'',
                    'fig':[],
                    'figsize':(10,10),
                     'ax':[],
                     'showvel':False,
                     'velth':0.07,
                     'fo':5,
                     'fw':0.02,
                     'ylim':(),
                     'time_offset':0,
                     'color':'gray',
                     'hatch':'',
                     'label_pos':'top',
                     'label_pos_off':2,
                     'label_mob':'M',
                     'label_stat':'S'
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else :
            fig=kwargs['fig']

        if kwargs['ax'] ==[]:
            ax = fig.add_subplot(111)
        else :
            ax = kwargs['ax']

        if kwargs['subject']=='':
            subject=self.B.keys()[0]
        else:
            subject=kwargs['subject']

        V=self.B[subject].traj[['vx','vy']].values
        Vi=np.sqrt((V[:,0]**2+V[:,1]**2))
        f=DF()
        f.butter(kwargs['fo'],kwargs['fw'],'lowpass')
        Vif=f.filter(Vi)

        if kwargs['time_offset']>=0:
            zmo = np.zeros(kwargs['time_offset'])
            tmp = np.insert(Vif,zmo,0)
            Vif = tmp[:len(Vif)]
        else:
            zmo = np.zeros(-kwargs['time_offset'])
            tmp = np.concatenate((Vif,zmo))
            Vif = tmp[-kwargs['time_offset']:len(Vif)-kwargs['time_offset']]


        if kwargs['showvel']:
            fig2 = plt.figure()
            ax2=fig2.add_subplot(111)
            ax2.plot(self.B[subject].time[:-2],Vif)
            ax2.plot(Vif)
            cursor2 = Cursor(ax2, useblit=True, color='gray', linewidth=1)

        null = np.where(Vif<kwargs['velth'])[0]
        unu1 = np.where(np.diff(null)!=1)[0]
        unu2 = np.where(np.diff(null[::-1])!=-1)[0]
        unu2 = len(null)-unu2
        unu = np.concatenate((unu1,unu2))
        unu = np.sort(unu)
        sunu = unu.shape
        if sunu[0]%2:
            unu=np.insert(unu,-1,len(null)-1)
            sunu = unu.shape
        nullr=null[unu].reshape(sunu[0]/2,2)

        if kwargs['ylim'] != ():
            ylim = kwargs['ylim']
        else :
            axlim = ax.axis() 
            ylim = [axlim[2],axlim[3]]


        fig , ax =plu.rectplot(self.B[subject].time,nullr,ylim=ylim,
                                color=kwargs['color'],
                                hatch=kwargs['hatch'],
                                fig=fig,ax=ax)


        inullr = copy.copy(nullr)
        bb = np.insert(inullr[:,1],0,0)
        ee = np.hstack((inullr[:,0],null[-1]))
        inullr = np.array((bb,ee)).T
        # remove last
        inullr = inullr[:-1,:]

        if kwargs['label_pos']!='':
            if kwargs['label_pos'] == 'top':
                yposM = ylim[1]-kwargs['label_pos_off']+0.5
                yposS = ylim[1]-kwargs['label_pos_off']-0.5

            elif kwargs['label_pos'] == 'bottom':
                yposM = ylim[0]+kwargs['label_pos_off']+0.5
                yposS = ylim[0]+kwargs['label_pos_off']+0.5
            xposM= self.B[subject].time[nullr.mean(axis=1).astype(int)]
            xposS= self.B[subject].time[inullr.mean(axis=1).astype(int)]
            [ax.text(x,yposM,kwargs['label_mob']+str(ix+1),
                        horizontalalignment='center',
                        verticalalignment='center')
                        for ix,x in enumerate(xposM)]
            [ax.text(x,yposS,kwargs['label_stat']+str(ix+1),
                        horizontalalignment='center',
                        verticalalignment='center')
                        for ix,x in enumerate(xposS)]

        return fig,ax

    def animhkb(self,a,b,interval=10,save=False):
        """
        Parameters
        ----------

        a : node name |number
        b : node name | number
        save : bool
        """

        import matplotlib.animation as animation

        x = self.hkb.index
        link = a+'-'+b
        y = self.hkb[link].values

        fig, ax = plt.subplots()
        plt.xlim(0,x[-1])
        line = [ax.plot(x, y, animated=True)[0]]

        def animate(i):
            line[0].set_ydata(y[:i])
            line[0].set_xdata(x[:i])
            return line

        ani = animation.FuncAnimation(fig, animate, xrange(1, len(x)),
                                      interval=interval, blit=True)
        if save:
            ani.save(link+'.mp4')
        plt.title(link)
        plt.xlabel('time (s)')
        plt.ylabel('RSS (dBm)')
        plt.show()


    def animhkbAP(self,a,AP_list,interval=1,save=False,**kwargs):
        """
        Parameters
        ----------

        a : node name
        AP_nb=[]
        save : bool

        Example
        -------
            >>> from pylayers.measures.cormoran import *
            >>> S = CorSer(6)
            >>> S.animhkbAP('TorsoTopLeft',['AP1','AP2','AP3','AP4'],interval=100,xstart=58,figsize=(20,2))

        """
        import matplotlib.animation as animation

        defaults = {   'fig':[],
                       'figsize':(10,10),
                        'ax':[],
                        'label':'',
                        'xstart':0
        }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else :
            fig=kwargs['fig']

        if kwargs['ax'] ==[]:
            ax = fig.add_subplot(111)
        else :
            ax = kwargs['ax']


        ust = np.where(self.hkb.index>=kwargs['xstart'])[0][0]

        x = self.hkb.index[ust:]
        links = [l+'-'+a for l in AP_list]
        ly = [self.hkb[l].values[ust:] for l in links]

        color=['k','b','g','r']
        plt.xlim(kwargs['xstart'],x[-1]+3)
        line = [ax.plot(x, y, animated=True,
                        color=color[iy],
                        label=AP_list[iy]+'-'+kwargs['label'])[0] for iy,y in enumerate(ly)]

        def animate(i):
            for iy,y in enumerate(ly):
                line[iy].set_ydata(y[:i])
                line[iy].set_xdata(x[:i])

            return line
        plt.legend()
        plt.xlabel('time (s)')
        plt.ylabel('RSS (dBm)')
        ani = animation.FuncAnimation(fig, animate, xrange(0, len(x)),
                                      interval=interval, blit=True)
        if save:
            ani.save(a+'.mp4')
        #plt.title(links)
        plt.show()


    def plot(self,a,b,techno='',t='',**kwargs):
        """
        Parameters
        ----------

        a : str | int
            name |id
        b : str | int
            name |id
        techno : str (optional)
            radio techno
        t : float | list (optional)
            given time
            or [start,stop] time

        color : color
        distance : boolean (False)
            plot distance instead of value
        lin : boolean (False)
            display linear value instead of dB
        sqrtinv : boolean (False)
            apply : "sqrt (1/ dataset)"
        xoffset : float (0)
            add an offset on x axis
        yoffset :  float (1|1e3|1e6)
            add an offset on y axis
        
        title : boolean (True)
            display title
        shortlabel : boolean (True)
            enable short labelling
        fontsize : int (18)
            font size

        returnlines : boolean
            if True return the matplotlib ploted lines

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S = CorSer(6)
        >>> f,ax = S.plot('AP1','TorsoTopLeft',techno='HKB')
        >>> f,ax = S.pltvisi('AP1','TorsoTopLeft',techno='HKB',fig=f,ax=ax)
        >>> #f,ax = S.pltmob(fig=f,ax=ax)
        >>> #plt.title('hatch = visibility / gray= mobility')
        >>> plt.show()

        """

        defaults = { 'fig':[],
                     'ax':[],
                     'figsize':(6,4),
                     'color':'g',
                     'distance':False,
                     'lin':False,
                     'xoffset':0,
                     'yoffset': 1e6,
                     'sqrtinv':False,
                     'title':True,
                     'shortlabel':True,
                     'fontsize':18,
                     'returnlines':False
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        a,ia,bia,subja,techno=self.devmapper(a,techno)
        b,ib,bib,subjb,techno=self.devmapper(b,techno)

        ###create a short labeling
        if kwargs['shortlabel']:

            #find uppercase position
            uu =  np.nonzero([l.isupper() or l.isdigit() for l in a])[0]
            #cretae string from list
            labela = ''.join([a[i] for i in uu])

            uu =  np.nonzero([l.isupper() or l.isdigit() for l in b])[0]
            #cretae string from list
            labelb = ''.join([b[i] for i in uu])
            label = labela +'-'+labelb

            
        else:
            label = a+'-'+b

        if kwargs['distance']:
            label = 'dist ' + label

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else :
            fig=kwargs['fig']

        if kwargs['ax'] ==[]:
                ax = fig.add_subplot(111)
        else :
            ax = kwargs['ax']


        # get dataframe
        if not kwargs['distance']:
            df = self.getlink(a,b,techno,t)
            title = 'Received Power between ' + label
            ylabel = 'Received Power dBm'
        else :
            df = self.getlinkd(a,b,techno,t)
            title = 'Distance  between ' + label
            ylabel = 'distance (m)'


        #post processing on dataframe
        if kwargs['lin']:
            df = 10**(df/10) * kwargs['yoffset']
            
        if kwargs['sqrtinv']:
            df = np.sqrt(1./df)
            ylabel = u'$ (mW)^{-1/2} linear scale$'

        lines = df.plot(ax=ax,color=kwargs['color'],label=label)

        # Managing labelling
        if kwargs['title']:
            ax.set_title(label=title,fontsize=kwargs['fontsize'])
        if kwargs['lin']:
            if kwargs['yoffset']==1:
                ylabel = 'mW'
            if kwargs['yoffset']==1e3:
                ylabel = u'$\micro$W'
            if kwargs['yoffset']==1e6:
                ylabel = u'nW'

        ax.set_ylabel(ylabel)
        # if kwargs['data']==True:
        #     #ax.plot(self.thkb[0],self.rssi[ia,ib,:])
        #     #ax.plot(self.thkb[0],self.rssi[ib,ia,:])
        #     sab = self.hkb[a+'-'+b]

        #     if not(kwargs['dB']):
        #         sab = 10**(sab/10) * kwargs['yoffset']
        #         if kwargs['distance']:
        #             sab = np.sqrt(1/sab)
        #         if kwargs['reciprocal']:
        #             sba = 10**(sba/10 ) * kwargs['yoffset']
        #             sba = np.sqrt(1/sba)
        #     sab[t0:t1].plot(ax=ax,color=kwargs['colorab'],label=label,xlim=(t0,t1))
        #     if kwargs['reciprocal']:
        #         sba[t0:t1].plot(ax=ax,color=kwargs['colorba'],label=label)

        #     #title = 'Received Power   ' + self.title1
        #     if kwargs['dis_title']:
        #         #title = self.title1+kwargs['tit']
        #         title = kwargs['tit']
        #         ax.set_title(label=title,fontsize=kwargs['fontsize'])
        #     if not kwargs['distance']:
        #         if kwargs['dB']:
        #             ax.set_ylabel('Received Power dBm')
        #         else:
        #             if kwargs['yoffset']==1:
        #                 ax.set_ylabel('mW')
        #             if kwargs['yoffset']==1e3:
        #                 ax.set_ylabel(u'$\micro$W')
        #             if kwargs['yoffset']==1e6:
        #                 ax.set_ylabel(u'nW')

        #     else:
        #         ax.set_ylabel(u'$\prop (mW)^{-1/2} linear scale$')

        # if kwargs['reciprocal']==True:
        #     # if kwargs['data']==True:
        #     #     ax2=fig.add_subplot(212)
        #     r = self.hkb[a+'-'+b][self.hkb[a+'-'+b]!=0]- self.hkb[b+'-'+a][self.hkb[b+'-'+a]!=0]
        #     r[t0:t1].plot(ax=ax2)
        #     ax2.set_title('Reciprocity offset',fontsize=kwargs['fontsize'])
        if not kwargs['returnlines']:
            return fig,ax
        else:
            return fig,ax,lines


    def plthkb(self,a,b,techno='HKB',**kwargs):
        """
        DEPRECATED
        Parameters
        ----------

        a : node name |number
        b : node name | number
        t0 : start time
        t1 : stop time


        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S = CorSer(6)
        >>> f,ax = S.plthkb('AP1','TorsoTopLeft',techno='HKB')
        >>> f,ax = S.pltvisi('AP1','TorsoTopLeft',techno='HKB',fig=f,ax=ax)
        >>> f,ax = S.pltmob(fig=f,ax=ax)
        >>> plt.title('hatch = visibility / gray= mobility')
        >>> plt.show()

        """

        defaults = { 't0':0,
                     't1':-1,
                     'fig':[],
                     'ax':[],
                     'figsize':(8,8),
                     'xoffset':0,
                     'yoffset': 1e6,
                     'reciprocal':False,
                     'dB':True,
                     'data':True,
                     'colorab':'g',
                     'colorba':'b',
                     'distance':False,
                    'fontsize':18,
                    'shortlabel':True,
                    'dis_title':True,
                    'xlim':(),
                    'tit':''
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        t0 =kwargs['t0']
        t1 =kwargs['t1']
        if t1 ==-1:
            try:
                t1=self.thkb[0][-1]
            except:
                t1=self.thkb[-1]

        a,ia,bia,subja,technoa=self.devmapper(a,techno)
        b,ib,bib,subjb,technob=self.devmapper(b,techno)

        if kwargs['shortlabel']:

            #find uppercase position
            uu =  np.nonzero([l.isupper() or l.isdigit() for l in a])[0]
            #cretae string from list
            labela = ''.join([a[i] for i in uu])

            uu =  np.nonzero([l.isupper() or l.isdigit() for l in b])[0]
            #cretae string from list
            labelb = ''.join([b[i] for i in uu])

            label = labela +'-'+labelb
        else:
            label = a+'-'+b

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else :
            fig=kwargs['fig']

        if kwargs['ax'] ==[]:
            if kwargs['reciprocal']:
                ax = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
            else :
                ax = fig.add_subplot(111)
        else :
            ax = kwargs['ax']


        if kwargs['data']==True:
            #ax.plot(self.thkb[0],self.rssi[ia,ib,:])
            #ax.plot(self.thkb[0],self.rssi[ib,ia,:])
            sab = self.hkb[a+'-'+b]

            if not(kwargs['dB']):
                sab = 10**(sab/10) * kwargs['yoffset']
                if kwargs['distance']:
                    sab = np.sqrt(1/sab)
                if kwargs['reciprocal']:
                    sba = 10**(sba/10 ) * kwargs['yoffset']
                    sba = np.sqrt(1/sba)
            sab[t0:t1].plot(ax=ax,color=kwargs['colorab'],label=label,xlim=(t0,t1))
            if kwargs['reciprocal']:
                sba[t0:t1].plot(ax=ax,color=kwargs['colorba'],label=label)

            #title = 'Received Power   ' + self.title1
            if kwargs['dis_title']:
                #title = self.title1+kwargs['tit']
                title = kwargs['tit']
                ax.set_title(label=title,fontsize=kwargs['fontsize'])
            if not kwargs['distance']:
                if kwargs['dB']:
                    ax.set_ylabel('Received Power dBm')
                else:
                    if kwargs['yoffset']==1:
                        ax.set_ylabel('mW')
                    if kwargs['yoffset']==1e3:
                        ax.set_ylabel(u'$\micro$W')
                    if kwargs['yoffset']==1e6:
                        ax.set_ylabel(u'nW')

            else:
                ax.set_ylabel(u'$\prop (mW)^{-1/2} linear scale$')

        if kwargs['reciprocal']==True:
            # if kwargs['data']==True:
            #     ax2=fig.add_subplot(212)
            r = self.hkb[a+'-'+b][self.hkb[a+'-'+b]!=0]- self.hkb[b+'-'+a][self.hkb[b+'-'+a]!=0]
            r[t0:t1].plot(ax=ax2)
            ax2.set_title('Reciprocity offset',fontsize=kwargs['fontsize'])

        return fig,ax

    def plttcr(self,a,b,**kwargs):
        """
        Parameters
        ----------

        a : node name |number
        b : node name | number
        t0 : start time
        t1 : stop time

        """

        defaults = { 't0':0,
                     't1':-1,
                     'fig':[],
                     'ax':[],
                     'figsize':(8,8),
                     'data':True,
                     'colorab':'g',
                     'colorba':'b',
                     'linestyle':'default',
                     'inverse':False
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        t0 =kwargs['t0']
        t1 =kwargs['t1']
        if t1 ==-1:
            t1=self.ttcr[-1]

        if isinstance(a,str):
            ia = self.dTCR[a]
        else:
            ia = a
            a = self.idTCR[a]

        if isinstance(b,str):
            ib = self.dTCR[b]
        else:
            ib = b
            b = self.idTCR[b]


        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else:
            fig = kwargs['fig']

        if kwargs['ax'] ==[]:
            ax = fig.add_subplot(111)
        else :
            ax=kwargs['ax']

        if kwargs['data']==True:
            #ax.plot(self.thkb[0],self.rssi[ia,ib,:])
            #ax.plot(self.thkb[0],self.rssi[ib,ia,:])
            if kwargs['inverse']:
                sab = 1./(self.tcr[a+'-'+b])**2
                sba = 1./(self.tcr[b+'-'+a])**2
            else:
                sab = self.tcr[a+'-'+b]
                sba = self.tcr[b+'-'+a]
            sab[t0:t1].plot(ax=ax,color=kwargs['colorab'],marker='o',linestyle=kwargs['linestyle'])
            sba[t0:t1].plot(ax=ax,color=kwargs['colorba'],marker='o',linestyle=kwargs['linestyle'])
            ax.set_title(a+'-'+b)

        return fig,ax


    def pltgt(self,a,b,**kwargs):
        """ plt ground truth

        Parameters
        ----------

        t0
        t1
        fig
        ax
        figsize: tuple
        linestyle'
        inverse :False,
            display 1/distance  instead of distance
        log : boolean
            display log fo distance intead of distance
        gammma':1.,
            mulitplication factor for log : gamma*log(distance) 
            this can be used to fit RSS
        mode : string
            'HKB' | 'TCR' | 'FULL'
        visi : boolean,
            display visibility
        color: string color ('k'|'m'|'g'),
            color to display the visibility area
        hatch': strin hatch type ('//')
            hatch type to hatch visibility area
        fontsize: int
            title fontsize

        Example
        -------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(6)
        >>> S.pltgt('AP1','TorsoTopLeft')


        """

        defaults = { 'subject':'',
                     't0':0,
                     't1':-1,
                     'fig':[],
                     'ax':[],
                     'figsize':(8,8),
                     'linestyle':'default',
                     'inverse':False,
                     'log':True,
                     'gamma':-40,
                     'mode':'HKB',
                     'visi': True,
                     'fontsize': 14,
                     'color':'k',
                    'hatch':''
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]




        t0 =kwargs.pop('t0')
        t1 =kwargs.pop('t1')
        if t1 ==-1:
            #t1=self.thkb[-1]
            t1=self.ttcr[-1]


        label = a+'-'+b

        mode = kwargs.pop('mode')
        inverse = kwargs.pop('inverse')
        log = kwargs.pop('log')
        gamma = kwargs.pop('gamma')
        visibility = kwargs.pop('visi')
        fontsize = kwargs.pop('fontsize')
        hatch = kwargs.pop('hatch')
        subject = kwargs.pop('subject')


        if subject=='':
            subject=self.B.keys()[0]
        else:
            subject=subject

        if kwargs['fig']==[]:
            figsize = kwargs.pop('figsize')
            kwargs.pop('fig')
            fig = plt.figure(figsize=figsize)
        else:
            kwargs.pop('figsize')
            fig = kwargs.pop('fig')
        if kwargs['ax'] ==[]:
            kwargs.pop('ax')
            ax = fig.add_subplot(111)
        else :
            ax=kwargs.pop('ax')



        if mode == 'HKB' or mode == 'FULL':

            if isinstance(a,str):
                iahk = self.dHKB[a]
            else:
                iahk = a
                a = self.idHKB[a]

            if isinstance(b,str):
                ibhk = self.dHKB[b]
            else:
                ibhk = b
                b = self.idHKB[b]

            var = self.getlink(iahk,ibhk,'HKB').values
            if inverse:
                var = 1./(var)
                ax.set_ylabel(u'$m^{-2}$',fontsize=fontsize)
                if log :
                    #var = gamma*10*np.log10(var)
                    var = 20*np.log10(var)+gamma
                    ax.set_ylabel(u'$- 20 \log_{10}(d)'+str(gamma)+'$  (dB)',fontsize=fontsize)
                    plt.ylim(-65,-40)
            else:
                ax.set_ylabel(u'meters',fontsize=fontsize)
                if log :
                    var = gamma*10*np.log10(var)+gamma
                    ax.set_ylabel(u'$10log_{10}m^{-2}$',fontsize=fontsize)


            ax.plot(self.B[subject].time,var,label=label,**kwargs)
        #
        # TCR |Full
        #
        if mode == 'TCR' or mode == 'FULL':

            if isinstance(a,str):
                iatcr = self.dTCR[a]
            else:
                iatcr = a
                a = self.idTCR[a]

            if isinstance(b,str):
                ibtcr = self.dTCR[b]
            else:
                ibtcr = b
                b = self.idTCR[b]

            var = self.getlink(iatcr,ibtcr,'TCR').values

            #if inverse:
            #    var = 1./(var)**2
            #    if log :
            #        var = gamma*10*np.log10(var)
            #else:
            #    if log :
            #        var = gamma*10*np.log10(var)

            #pdb.set_trace()
            ax.plot(self.B[subject].time,var,**kwargs)


        if visibility:
            aa= ax.axis()
            vv,tv,tseg,itseg = self._visiarray(a,b)
            # vv.any : it exist NLOS regions
            if vv.any():
                fig,ax=plu.rectplot(tv,tseg,ylim=aa[2:],color=kwargs['color'],hatch=hatch,fig=fig,ax=ax)
                # for t in tseg:


        #axs[cptax].plot(visi.index.values,visi.values,'r')


        #if inverse:
        #    ax.set_title(u'Motion Capture Ground Truth :  inverse of squared distance',fontsize=fontsize+1)
        #else:
        #    ax.set_title('Motion Capture Ground Truth : evolution of distance (m)',fontsize=fontsize+1)

        ax.set_xlabel('Time (s)',fontsize=fontsize)
        plt.tight_layout()

        return fig, ax


    def pltlk(self,a,b,**kwargs):
        """ plt links

        Parameters
        ----------

        a : string
            node a name
        b : string
            node b name

        display: list
            techno to be displayed
        figsize
        t0: float
            time start
        t1 : float
            time stop
        colhk: plt.color

            color of hk curve
        colhk2:plt.color
            color of hk curve2 ( if recirpocal)
        linestylehk:
            linestyle hk

        coltcr:
            color tcr curve
        coltcr2:
            color of tcr curve2 ( if recirpocal)
        linestyletcr:
            linestyle tcr
        colgt:
            color ground truth
        inversegt:
            invert ground truth
        loggt: bool
            apply a log10 factor to ground truth
        gammagt:
            applly a gamma factor to ground truth (if loggt ! )
        fontsize:
            font size of legend
        visi:
            display visibility indicator
        axs :
            list of matplotlib axes

        Example
        -------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(6)
        >>> S.pltlk('AP1','TorsoTopLeft')
        """

        defaults = { 'display':[],
                     'figsize':(8,8),
                     't0':0,
                     't1':-1,
                     'colhk':'g',
                     'colhk2':'b',
                     'linestylehk':'default',
                     'coltcr':'g',
                     'coltcr2':'b',
                     'linestyletcr':'step',
                     'colgt': 'k',
                     'inversegt':True,
                     'loggt':True,
                     'gammagt':-40,
                     'fontsize':14,
                     'visi':True,
                     'axs' :[],
                     'gt':True,
                     'tit':''
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        display = kwargs.pop('display')

        if not isinstance(display,list):
            display=[display]


        if display == []:
            if ('tcr' in dir(self)) and ('hkb' in dir(self)):
                display.append('FULL')
            elif 'tcr' in dir(self):
                display.append('TCR')
            elif 'hkb' in dir(self):
                display.append('HKB')

        display = [t.upper() for t in display]

        if 'FULL' in display:
            ld = 2
        elif 'TCR' in display or 'HKB' in display:
            ld = 2

        #Axes management
        if kwargs['axs'] == []:
            kwargs.pop('axs')
            fig,axs = plt.subplots(nrows=ld,ncols=1,figsize=kwargs['figsize'],sharex=True)
        else :
            fig =plt.gcf()
            axs = kwargs.pop('axs')


        cptax= 0


        # HKB plot
        if 'HKB' in display or 'FULL' in display:
            if ('HKB' in self.typ.upper()) or ('FULL' in self.typ.upper()):
                if isinstance(a,str):
                    iahk = self.dHKB[a]
                else :
                    raise AttributeError('in self.pltlk, nodes id must be a string')
                if isinstance(b,str):
                    ibhk = self.dHKB[b]
                else :
                    raise AttributeError('in self.pltlk, nodes id must be a string')

            else :
                raise AttributeError('HK not available for the given scenario')






            kwargs['fig']=fig
            kwargs['ax']=axs[cptax]
            kwargs['colorab']=kwargs.pop('colhk')
            kwargs['colorba']=kwargs.pop('colhk2')
            kwargs['linestyle']=kwargs.pop('linestylehk')
            kwargs['tit']=kwargs.pop('tit')

            fig,axs[cptax]=self.plthkb(a,b,reciprocal=False,**kwargs)


            cptax+=1
        else :
            kwargs.pop('colhk')
            kwargs.pop('colhk2')
            kwargs.pop('linestylehk')


        #TCR plot
        if 'TCR' in display or 'FULL' in display:
            if ('TCR' in self.typ.upper()) or ('FULL' in self.typ.upper()):
                if isinstance(a,str):
                    iatcr = self.dTCR[a]
                else :
                    raise AttributeError('in self.pltlk, nodes id must be a string')
                if isinstance(b,str):
                    ibtcr = self.dTCR[b]
                else :
                    raise AttributeError('in self.pltlk, nodes id must be a string')
            else :
                raise AttributeError('TCR not available for the given scenario')

            kwargs['fig']=fig
            kwargs['ax']=axs[cptax]
            kwargs['colorab']=kwargs.pop('coltcr')
            kwargs['colorba']=kwargs.pop('coltcr2')
            kwargs['linestyle']=kwargs.pop('linestyletcr')
            tcrlink = a+'-'+b
            #plot only if link exist
            if tcrlink in self.tcr:
                fig,axs[cptax]=self.plttcr(a,b,**kwargs)
        else :
            kwargs.pop('coltcr')
            kwargs.pop('coltcr2')
            kwargs.pop('linestyletcr')
            #cptax+=1

        #
        # Ground Truth
        #
        #
        # HKB |Full
        #
        if kwargs.pop('gt'):
            kwargs['color'] = kwargs.pop('colgt')
            kwargs.pop('colorab')
            kwargs.pop('colorba')
            kwargs['ax']=axs[cptax]
            kwargs['inverse']=kwargs.pop('inversegt')
            kwargs['log']=kwargs.pop('loggt')
            kwargs['gamma']=kwargs.pop('gammagt')
            kwargs.pop('tit')

            if 'HKB' in display or 'FULL' in display:
                kwargs['mode']= 'HKB'
                fig,axs[cptax] = self.pltgt(a,b,**kwargs)
            elif 'TCR' in display or 'FULL' in display:
                kwargs['mode']= 'TCR'
                fig,axs[cptax] = self.pltgt(a,b,**kwargs)

        return fig,axs
        # aa = axs[cptax].axis()
        #
        # calculates visibility and display NLOS region
        # as a yellow patch over the shadowed region
        #


    def showpattern(self,a,techno='HKB',**kwargs):
        """ show pattern configuation for a given link and frame

        Parameters
        ----------

        a : int
            link index
        technoa : string
            'HKB'|'TCR'|'BS'
        technob
            default 'HKB'|'TCR'|'BS'
        phi : float
            antenna elevation in rad
        fig :
        ax  :
        t   : float
        phi : float
            pi/2
        ap  : boolean
            
        """

        defaults = { 'fig':[],
                     'ax':[],
                     't':0,
                     'phi':np.pi/2.,
                     'ap':False
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if kwargs['fig'] == []:
            fig=plt.figure()
        else :
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            ax=fig.add_subplot(111)
        else :
            ax = kwargs['ax']

        # display nodes
        #
        #
        #
        a,ia,ba,subjecta,techno = self.devmapper(a,techno)
        pa = self.getdevp(a,techno=techno,t=kwargs['t']).values


        if len(pa.shape) >1:
            pa=pa[0]
        ax.plot(pa[0],pa[1],'ob')
        ax.text(pa[0],pa[1],ba)

        if subjecta != '':
            self.B[subjecta].settopos(t=kwargs['t'])
            self.B[subjecta].dev[ba]['ant'].Fsynth()
            xa,ya,z,sa,v = self.B[subjecta].dev[ba]['ant']._computemesh(po=pa,T=self.B[subjecta].acs[ba],minr=0.01,maxr=0.1,ilog=False)
            p2 = np.where(self.B[subjecta].dev[ba]['ant'].phi<=kwargs['phi'])[0][-1]

            # ax.plot(xa[:,p2],ya[:,p2])
            ax.plot(xa[p2,:],ya[p2,:])

        else:
            self.din[ba]['ant'].Fsynth()
            xa,ya,z,sa,v = self.din[ba]['ant']._computemesh(po=self.din[ba]['p'],T=self.din[ba]['T'],minr=0.01,maxr=0.1,ilog=False)
            p2 = np.where(self.din[ba]['ant'].phi<=kwargs['phi'])[0][-1]
            ax.plot(xa[:,p2],ya[:,p2])

        return fig,ax

    def showlink(self,a='AP1',b='BackCenter',technoa='HKB',technob='HKB',**kwargs):
        """ show link configuation for a given frame

        Parameters
        ----------

        a : int
            link index
        b : int
            link index
        technoa : string 
            default 'HKB'|'TCR'|'BS'
        technob
            default 'HKB'|'TCR'|'BS'
        phi : float
            antenna elevation in rad

        """

        defaults = { 'fig':[],
                     'ax':[],
                     't':0,
                     'phi':np.pi/2.,
                     'ap':False
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if kwargs['fig'] == []:
            fig=plt.figure()
        else :
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            ax=fig.add_subplot(111)
        else :
            ax = kwargs['ax']

        # display nodes

        fig,ax=self.showpattern(a=a,techno=technoa,fig=fig,ax=ax)
        fig,ax=self.showpattern(a=b,techno=technob,fig=fig,ax=ax)

        plt.axis('equal')


        p1 = self.din['HKB:1']['p']
        p2 = self.din['HKB:2']['p']
        p3 = self.din['HKB:3']['p']
        p4 = self.din['HKB:4']['p']
        plt.plot(p1[0],p1[1],'og')
        plt.plot(p2[0],p2[1],'ob')
        plt.plot(p3[0],p3[1],'or')
        plt.plot(p4[0],p4[1],'ok')
        plt.axis('equal')
        
        # if A.ndim==2:
        #     plt.plot(A[iframe,0],A[iframe,1],'ob')
        #     plt.text(A[iframe,0],A[iframe,1],a)
        # else:
        #     plt.plot(A[0],A[1],'or')
        #     #plt.text(A[0],A[1],a)

        # if B.ndim==2:
        #     plt.plot(B[iframe,0],B[iframe,1],style)
        #     plt.text(B[iframe,0]+0.1,B[iframe,1]+0.1,b)
        # else:
        #     plt.plot(B[0],B[1],'ob')
        #     plt.text(B[0],B[1],b)
        # plt.xlim(-6,6)
        # plt.ylim(-5,5)


        # self.B[subjecta].settopos(t=t)
        # self.B[subjectb].settopos(t=t)

        # 

        
        # # display body

        # #pc = self.B.d[:,2,iframe] + self.B.pg[:,iframe].T
        # pc0 = self.B[subjecta].d[:,0,iframe] + self.B[subjecta].pg[:,iframe].T
        # pc1 = self.B[subjecta].d[:,1,iframe] + self.B[subjecta].pg[:,iframe].T
        # pc15 = self.B[subjecta].d[:,15,iframe] + self.B[subjecta].pg[:,iframe].T
        # #plt.plot(pc0[0],pc0[1],'og')
        # #plt.text(pc0[0]+0.1,pc0[1],str(iframe))
        # #plt.plot(pc1[0],pc1[1],'og')
        # #plt.plot(pc15[0],pc15[1],'og')
        # #ci00   = plt.Circle((pc0[0],pc0[1]),self.B[subjecta].sl[0,2],color='green',alpha=0.6)
        # #ci01   = plt.Circle((pc1[0],pc1[1]),self.B[subjecta].sl[0,2],color='green',alpha=0.1)
        # #ci100 = plt.Circle((pc0[0],pc0[1]),self.B[subjecta].sl[10,2],color='red',alpha=0.1)
        # ci1015 = plt.Circle((pc15[0],pc15[1]),self.B[subjecta].sl[10,2],color='green',alpha=0.5)
        # plt.axis('equal')
        # ax = plt.gca()
        # ax.add_patch(ci1015)
        # #ax.add_patch(ci01)
        # #ax.add_patch(ci100)
        # #ax.add_patch(ci1015)
        # #its = self.B[subjecta].intersectBody(A[iframe,:],B[iframe,:],topos=False,frameId=iframe)
        # #x.set_title('frameId :'+str(iframe)+' '+str(its.T))
        



    def visidev(self,a,b,technoa='HKB',technob='HKB',dsf=10):
        """ get link visibility status

        Returns
        -------

        visi : pandas Series
            0  : LOS
            1  : NLOS

        """

        A,B = self.getlinkp(a,b)
        A=A.values
        B=B.values
        aa,ia,ba,subjecta,technoa= self.devmapper(a,technoa)
        ab,ib,bb,subjectb,technob= self.devmapper(b,technob)


        if 'AP' not in aa:
            Nframe = A.shape[0]
        if 'AP' not in ab:
            Nframe = B.shape[0]
        else: 
            Nframe = len(self.B[self.B.keys()[0]].time)
        iframe = np.arange(0,Nframe-1,dsf)
        tvisi = []
        #
        # A : Nframe x 3
        # B : Nframe x 3
        # B.pg : 3 x Nframe
        #
        if subjecta != '':
            subject = subjecta
        elif subjectb != '':
            subject = subjectb
        else :
            raise AttributeError('Visibility can only be determine on a body for now')
        if self.B[subject].centered:
            A = A-self.B[subject].pg.T
            B = B-self.B[subject].pg.T


        for k in iframe:
            if len(np.shape(A))<2:
                A=A[np.newaxis,:]*np.ones((len(B),3))
            if len(np.shape(B))<2:
                B=B[np.newaxis,:]*np.ones((len(A),3))

            its = self.B[subject].intersectBody(A[k,:],B[k,:],topos=False,frameId=k)
            tvisi.append(its.any())
        visi = pd.Series(tvisi,index=iframe/100.)
        #return(visi,iframe)
        return(visi)

    def visidev2(self,a,b,technoa='HKB',technob='HKB',trange=[]):
        """ get link visibility status

        Returns
        -------
        trange : nd array
            time range
        visi : pandas Series
            0  : LOS
            1  : NLOS

        """

        A,B = self.getlinkp(a,b,technoa,technob)
        A=A.values
        B=B.values
        aa,ia,ba,subjecta,technoa= self.devmapper(a,technoa)
        ab,ib,bb,subjectb,technob= self.devmapper(b,technob)

        if 'AP' not in a:
            Nframe = A.shape[0]
        if 'AP' not in b:
            Nframe = B.shape[0]
        # iframe = np.arange(0,Nframe-1,dsf)
        tvisi = []
        #
        # A : Nframe x 3
        # B : Nframe x 3
        # B.pg : 3 x Nframe
        #
        if subjecta != '':
            subject = subjecta
        elif subjectb != '':
            subject = subjectb
        else :
            raise AttributeError('Visibility can only be determine on a body for now')

        if self.B[subject].centered:
            A = A-self.B[subject].pg.T
            B = B-self.B[subject].pg.T

        for t in trange:
            fid = self.B[subject].posvel(self.B[subjecta].traj,t)[0]
            its = self.B[subject].intersectBody(A[fid,:],B[fid,:],topos=False,frameId=fid)
            tvisi.append(its.any())
        visi = pd.Series(tvisi,index=trange)
        #return(visi,iframe)
        return(visi)



    def _visiarray(self,a,b,technoa='HKB',technob='HKB'):
        """ create entries for plu.rectplot
        """

        visi = self.visidev(a,b)
        tv = visi.index.values
        vv = visi.values.astype(int)
        if (not(vv.all()) and vv.any()):
            df = vv[1:]-vv[0:-1]

            um = np.where(df==1)[0]
            ud = np.where(df==-1)[0]
            lum = len(um)
            lud = len(ud)

            #
            # impose same size and starting
            # on leading edge um and endinf on
            # falling edge ud
            #
            if lum==lud:
                if ud[0]<um[0]:
                    um = np.hstack((np.array([0]),um))
                    ud = np.hstack((ud,np.array([len(vv)-1])))
            else:
                if ((lum<lud) & (vv[0]==1)):
                    um = np.hstack((np.array([0]),um))

                if ((lud<lum) & (vv[len(vv)-1]==1)):
                    ud = np.hstack((ud,np.array([len(vv)-1])))


            tseg = np.array(zip(um,ud))
            #else:
            #    tseg = np.array(zip(ud,um))
        else:
            if vv.all():
                tseg = np.array(zip(np.array([0]),np.array([len(vv)-1])))
            else :
                tseg = np.array([[0,0]])
        itseg = copy.copy(tseg)
        bb = np.insert(itseg[:,1],0,0)
        ee = np.hstack((itseg[:,0],len(vv)))
        itseg = np.array((bb,ee)).T
        # bb = np.hstack((bb,len(vv)))

        return vv,tv,tseg,itseg


    # def _computedevpdf(self):
    #     """ create a timestamped data frame
    #         with all positions of devices
    #     """
    #     t=self.B.traj.time()
    #     pos = np.empty((len(t),12,3))
    #     for ik,k in enumerate(t):
    #         self.B.settopos(t=k)
    #         pos[ik,:,:]=self.B.getlinkp()
    #     df=[]
    #     for d in range(pos.shape[1]):
    #         df_tmp=pd.DataFrame(pos[:,d,:],columns=['x','y','z'],index=t)
    #         df_tmp['id']=self.B.dev.keys()[d]
    #         try :
    #             df = pd.concat([df,df_tmp])
    #         except:
    #             df = df_tmp
    #     df = df.sort_index()
    #     cols=['id','x','y','z']
    #     self.devdf=df[cols]

    def _computedevpdf(self):
        """ create a timestamped data frame
            with positions of all devices


        """


        if not isinstance(self.B,dict):
            B={self.subject[0]:self.B}
        else :
            B=self.B

        for b in B:
            if 'dev' in dir(B[b]):
                dev = B[b].dev.keys()
                udev=[B[b].dev[d]['uc3d'] for d in dev]

                postmp = np.array([np.mean(B[b]._f[:,u,:],axis=1) for u in udev])
                pos = postmp.swapaxes(0,1)
                t = B[b].time
                for d in range(len(dev)):
                    df_tmp=pd.DataFrame(pos[:,d,:],columns=['x','y','z'],index=t)
                    df_tmp[['vx','vy','vz']]=df_tmp.diff()/(t[1]-t[0])
                    df_tmp['v']=np.sqrt(np.sum(df_tmp[['vx','vy','vz']]**2,axis=1))
                    df_tmp[['ax','ay','az']]=df_tmp[['vx','vy','vz']].diff()/(t[1]-t[0])
                    df_tmp['a']=np.sqrt(np.sum(df_tmp[['ax','ay','az']]**2,axis=1))
                    df_tmp['id']=B[b].dev.keys()[d]
                    df_tmp['subject']=B[b].name
                    try :
                        df = pd.concat([df,df_tmp])
                    except:
                        df = df_tmp

        for i in self.din:
            pos = self.din[i]['p']
            pos2 = pos[:,np.newaxis]*np.ones(len(t))
            df_tmp=pd.DataFrame(pos2.T,columns=['x','y','z'],index=t)
            df_tmp['v']=0.
            df_tmp['vx']=0.
            df_tmp['vy']=0.
            df_tmp['vz']=0.
            df_tmp['a']=0.
            df_tmp['ax']=0.
            df_tmp['ay']=0.
            df_tmp['az']=0.
            df_tmp['subject']=''
            df_tmp['id']=i
            df = pd.concat([df,df_tmp])

        df = df.sort_index()
        cols=['id','subject','x','y','z','v','vx','vy','vz','a','ax','ay','az']
        self.devdf=df[cols]

        # if ('HKB' in self.typ) or ('FULL' in selftyp):
        #     self.devdf=self._align_devdf_on_hkb(self.devdf,self.hkb)


    def export_csv(self,**kwargs):
        """ export to csv devices positions

        Parameters
        ----------
            unit : string ('mm'|'cm'|'m'),
                unit of positions in csv(default mm)
            tunit: string
                time unit in csv (default 'ns')
            'alias': dict
                dictionnary to replace name of the devices into the csv .
                example : if you want to replace a device id named 'TCR:34'
                to an id = 5, you have to add an entry in the alias dictionnary as :
                alias.update({'TCR34':5})
            offset : np.array
                apply an offset on positions

        Return
        ------

            a csv file into the folder <PylayersProject>/netsave

        """

        defaults={'unit' :'mm',
                  'tunit':'ns',
                  'offset':np.array(([0,0,0])),
                  'alias':{}}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        unit=kwargs.pop('unit')
        tunit=kwargs.pop('tunit')
        alias = kwargs.pop('alias')

        if alias == {}:

            alias={'TCR:49':4 #Nicolas TorsoLeft
            ,'TCR:34':5 #Nicolas TorsoRight
            ,'TCR:48':6 #Nicolas Back
            ,'TCR:36':7 #Nicolas Shoulder

            ,'TCR:2':8 # Jihad TorsoLeft
            ,'TCR:35':9 #Jihad TorsoRight
            ,'TCR:33':10 #Jihad Back
            ,'TCR:37':11 #Jihad Shoulder

            ,'TCR:30':12 #Eric Torso
            ,'TCR:25':13 #Eric Back
            ,'TCR:26':14 # Eric Shoulder
            }

        filename =pyu.getlong(self._filename,pstruc['DIRNETSAVE']) + '.csv'

        df = copy.deepcopy(self.devdf)

        ldf = df[['id','x','y','z']]

        #rename devices
        if alias != {}:
            for k in alias:
                u=ldf['id'] == k
                ldf.iloc[u.values,0]=str(alias[k])

        # fix position unit
        if unit == 'm':
            _unit = 1.
        if unit == 'cm':
            _unit = 1e2
        if unit == 'mm':
            _unit = 1e3

        ldf.loc[:,'x']=ldf.loc[:,'x']*_unit-kwargs['offset'][0]
        ldf.loc[:,'y']=ldf.loc[:,'y']*_unit-kwargs['offset'][1]
        ldf.loc[:,'z']=ldf.loc[:,'z']*_unit-kwargs['offset'][2]

        # fix time unit
        if tunit == 'ms':
            _tunit = 1e3
        if tunit == 'us':
            _tunit = 1e6
        if tunit == 'ns':
            _tunit = 1e9


        # add timestamp column

        ldf['Timestamp']=ldf.index*_tunit

        ldf.to_csv(filename, sep = ' ',index=False)


    def getlinkd(self,a,b,techno='',t=''):
        """    get the distance for a link between devices

        Parameters
        ----------

        a : str | int
            name |id
        b : str | int
            name |id

        optional

        techno : str
            radio techno

        t : float | list
            given time
            or [start,stop] time


        Returns
        -------

        dist : np.array()
            all distances for all timestamps for the given link

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S = CorSer(serie=6)
        >>> d = S.getlinkd('AP1','WristLeft',techno='HKB')

        """

        ra = a
        rb = b

        a,ia,nna,subjecta,techno = self.devmapper(a,techno)
        b,ib,nnb,subjectb,techno = self.devmapper(b,techno)


        df = self.distdf


        if (nna +'-' + nnb) in self.distdf.keys():
            link = nna +'-' + nnb
        elif (nnb +'-' + nna) in self.distdf.keys():
            link = nnb +'-' + nna
        else :
            raise AttributeError('Link between ' + str(ra) +' and ' + str(rb) + ' is not available in distdf dataframe')

        #determine time
        if isinstance(t,list):
            tstart = t[0]
            tstop = t[-1]
            val = df[(df.index >= tstart) & (df.index <= tstop)][link]
        elif t == '':
            val = df[link]
        else :
            hstep = (df.index[1]-df.index[0])/2.
            val = df[(df.index >= t-hstep) & (df.index <= t+hstep)][link]



        return val



    def getlinkp(self,a,b,technoa='',technob='',t='',fId=''):
        """    get a link devices positions

        Parameters
        ----------

        a : str | int
            name |id
        b : str | int
            name |id


        optional :

        technoa : str
            radio techno
        technob : str
            radio techno

        t : float | list
            given time | [time_start,time_stop]

        OR 

        fId : int
            frame id

        Returns
        -------

        pa,pb : np.array()

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(serie=34)
        >>> a,b=S.getlinkp('AP1','WristLeft')

        """



        pa = self.getdevp(a,technoa,t,fId)
        pb = self.getdevp(b,technob,t,fId)

        return pa,pb


    def getlink(self,a,b,techno='',t=''):
        """    get a link value

        Parameters
        ----------

        a : str | int
            name |id
        b : str | int
            name |id


        optional :


        techno : str
            radio techno

        t : float | list
            given time
            or [start,stop] time

        Returns
        -------

        Pandas Serie

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(serie=34)
        >>> S.getlink('AP1','WristLeft')
        """



        ra = a
        rb = b
        a,ia,nna,subjecta,techno = self.devmapper(a,techno)
        b,ib,nnb,subjectb,techno = self.devmapper(b,techno)

        if ('HK' in techno) :
            if (a +'-' + b) in self.hkb.keys():
                link = a +'-' + b
            elif (b +'-' + a) in self.hkb.keys():
                link = b +'-' + a
            else :
                raise AttributeError('Link between ' + str(ra) +' and ' + str(rb) + ' is not available in hkb dataframe')
            df = self.hkb

        elif ('BS' in techno):
            if (a +'-' + b) in self.bespo.keys():
                link = a +'-' + b
            elif (b +'-' + a) in self.bespo.keys():
                link = b +'-' + a
            else:
                raise AttributeError('Link between ' + str(ra) +' and ' + str(rb) + ' is not available in bespo dataframe')

            df = self.bespo

        elif ('TCR' in techno):
            if (a +'-' + b) in self.tcr.keys():
                link = a +'-' + b
            elif (b +'-' + a) in self.tcr.keys():
                link = b +'-' + a
            else:
                raise AttributeError('Link between ' + str(ra) +' and ' + str(rb) + ' is not available in tcr dataframe')

            df = self.tcr


        #determine time
        if isinstance(t,list):
            tstart = t[0]
            tstop = t[-1]
            val = df[(df.index >= tstart) & (df.index <= tstop)][link]
        elif t == '':
            val = df[link]
        else :
            hstep = (df.index[1]-df.index[0])/2.
            val = df[(df.index >= t-hstep) & (df.index <= t+hstep)][link]


        return val



    def getdevp(self,a,techno='',t='',fId=''):
        """ get a  device position

        Parameters
        ----------

        a : str | int
            name |id
        techno : str
            radio techno

        optional :

        t : float | list 
            given time |[time_start,time_stop]


        Returns
        -------

        pa : np.array()

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(serie=34)
        >>> a=S.getdevp('AP1','WristLeft')

        """

        a,ia,nna,subjecta,techno = self.devmapper(a,techno)

        device_select = self.devdf['id'] == nna


        if isinstance(t,list):
            tstart = t[0]
            tstop = t[-1]
            #findex = np.where((self.tmocap>=tstart) & (self.tmocap<=tstop))[0]


        elif t == '':
            tstart = 0.0
            tstop = 1000000
        else :
            hstep = (self.devdf[device_select].index[1]-self.devdf[device_select].index[0])/2.
            tstart = t-hstep
            tstop = t+hstep


        pa = self.devdf[(self.devdf.index >= tstart) &
                        (self.devdf.index <= tstop) & 
                        device_select][['x','y','z']]

        return pa




    def devmapper(self,a,techno=''):
        """  retrieve name of device if input is number
             or
             retrieve number of device if input is name

        Parameters
        ----------

        a : str | int
            name |id |bodyid
        techno : str
            radio techno


        Returns
        -------

        a : string
            dev name
        ia : int
            dev number
        ba : string
            dev refernce in body
        subject : string
            body owning the device


        """
        subject=''

        #if a is a bodyid (e.g. 'HKB:16') or a body part (e.g. AnkleRight)
        if isinstance(a,str):
            # case where body id is given as input
            if ('HKB' in a) or ('TCR' in a ) or ('BS' in a ) :

                ba = a
                techno, ia = a.split(':')
                ia=int(ia)
                try:
                    if techno.upper() == 'TCR':
                        a = self.idTCR[ia]
                    elif (techno.upper() == 'HKB'):
                        a = self.idHKB[ia]
                    elif (techno.upper() == 'BS'):
                        a = self.idBS[ia]
                except:
                    raise AttributeError('No device ' + a + ' for techno ' +techno)

                for b in self.B:
                    if not 'Cylindre' in b:
                        if ba in self.B[b].dev.keys():
                            subject = b
                            break

            #case where body part (e.g. AnkleRight) is given. Here techno is mandatory
            else :

                if techno == '':
                    if self.typ != 'FULL':
                        if self.typ == 'HKBS':
                            raise AttributeError('Please indicate the radio techno in argument : HKB or BS')
                        else :
                            techno = self.typ
                    else:
                        raise AttributeError('Please indicate the radio techno in argument : TCR, HKB, BS')


                try :
                    if techno.upper() == 'TCR':
                        ia = self.dTCR[a]
                        ba='TCR:'+str(ia)
                    elif (techno.upper() == 'HKB') :
                        ia = self.dHKB[a]
                        ba='HKB:'+str(ia)
                    elif techno.upper() == 'BS':
                        ia = self.dBS[a]
                        ba='BS:'+str(ia)
                except:
                    raise AttributeError('No device on body part: ' + a + ' for techno ' +techno)


                for b in self.B:
                    if not 'Cylindre' in b:
                        if ba in self.B[b].dev.keys():
                            subject = b
                            break

        # an id (number) is given
        else:
            # techno autodetection raise an error if conflict and invite to precise radio techno
            if techno == '':
                if hasattr(self,'idHKB'):
                    if a in self.idHKB.keys() :
                        if techno == '':
                            techno = 'HKB'
                        else :
                            raise AttributeError('Please indicate the radio techno in argument : TCR, HKB, BS')
                if hasattr(self,'idBS'):
                    if a in self.idBS.keys():
                        if techno == '':
                            techno = 'BS'
                        else :
                            raise AttributeError('Please indicate the radio techno in argument : TCR, HKB, BS')
                if hasattr(self,'idTCR'):
                    if a in self.idTCR.keys():
                        if techno == '':
                            techno = 'TCR'
                        else :
                            raise AttributeError('Please indicate the radio techno in argument : TCR, HKB, BS')

            try :
                if techno.upper() == 'TCR':
                    ia = a
                    a = self.idTCR[a]
                    ba='TCR:'+str(ia)
                elif (techno.upper() == 'HKB') :
                    ia = a
                    a = self.idHKB[a]
                    ba='HKB:'+str(ia)

                elif (techno.upper() == 'BS') :
                    ia = a
                    a = self.idBS[a]
                    ba='BS:'+str(ia)
            except:
                raise AttributeError('No device with ID: ' + str(a) + ' for techno ' +techno)


            for b in self.B:
                if not 'Cylindre' in b:
                    if ba in self.B[b].dev.keys():
                        subject = b
                        break

        return a,ia,ba,subject,techno


    def align(self,devdf,hkbdf):

        """ DEPRECATED align time of 2 data frames:

        the time delta of the second data frame is applyied on the first one
        (e.g. time for devdf donwsampled by hkb data frame time)


        Parameters
        ----------

        devdf : device dataframe
        hkbdf : hkbdataframe

        Return
        ------

        devdfc : 
            aligned copy device dataframe
        hkbdfc : 
            aligned copy hkbdataframe

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(6)
        >>> devdf = S.devdf[S.devdf['id']=='HKB:15']
        >>> hkbdf = S.hkb['AP1-AnkleLeft']
        >>> devdf2,hkbdf2 = S.align(devdf,hkbdf)


        """

        print ('warning DEPRECATED')

        devdfc=copy.deepcopy(devdf)
        hkbdfc=copy.deepcopy(hkbdf)
        idev = devdfc.index
        ihkb = hkbdfc.index

        devdfc.index = pd.to_datetime(idev,unit='s')
        hkbdfc.index = pd.to_datetime(ihkb,unit='s')
        import ipdb
        ipdb.set_trace()
        sf = (hkbdfc.index[2]-hkbdfc.index[1]).microseconds
        devdfc= devdfc.resample(str(sf)+'U')

        devdfc.index = pd.Series([val.time() for val in devdfc.index])
        hkbdfc.index = pd.Series([val.time() for val in hkbdfc.index])

        return devdfc,hkbdfc


    def _apply_offset(self,techno):
        """ apply offset from self.offset[self._filename][techno+'_index']

            if offset >0
                add np.nan at the begining
            if offset <0
                first values of self.hkb will be dropped
        """

        if techno == 'HKB':
            df = self.hkb
        elif techno == 'TCR':
            df = self.tcr
        elif techno == 'BS':
            df = self.bespo
        else :
            raise AttributeError('Unknown tecnology got applying offset')


        offset = self.offset[self._filename][techno.lower()+'_index']
        if offset <= 0 :
            index = df.index
            df = df.iloc[-offset:]
            df.index = index[0:offset]
        else :
            #extract time values
            npahkbi = df.index.values
            step = npahkbi[1]- npahkbi[0]
            nstart = npahkbi[0]+ (step * (offset))
            df.index = pd.Index(npahkbi + nstart)

            #add blank at begining
            df = pd.DataFrame({},columns=df.keys(),index=npahkbi[:offset])
            ndf=pd.concat([df,df])
            df=ndf
        if techno == 'HKB':
            self.thkb = df.index
        elif techno == 'TCR':
            self.ttcr = df.index
        elif techno == 'BS':
            self.tbs = df.index


    def _apply_hkb_offset(self):
        """ apply offset from self.offset[self._filename]['hkb_index']

            if offset >0
                add np.nan at the begining
            if offset <0
                first values of self.hkb will be dropped
        """


        # offset = self.offset[self._filename]['hkb_index']
        # if offset >=0:
        #     self.hkb = self.hkb.iloc[offset:]
        # else :
        #     # new length
        #     lhkb = len(self.hkb) + (-offset)
        #     #extract time values
        #     npahkbi = self.hkb.index.values
        #     #calculate new termianl time
        #     step = npahkbi[1]- npahkbi[0]
        #     nstop = npahkbi[-1]+ (step * (-offset))
        #     ni = np.linspace(0,nstop,lhkb)
        #     df = pd.DataFrame({},columns=self.hkb.keys(),index=ni[0:-offset])

        #     self.hkb.index=pd.Index(ni[-offset:])
        #     ndf=pd.concat([df,self.hkb])

        #     self.hkb=ndf

        offset = self.offset[self._filename]['hkb_index']
        if offset <= 0 :
            index = self.hkb.index
            self.hkb = self.hkb.iloc[-offset:]
            self.hkb.index = index[0:offset]
        else :
            #extract time values
            npahkbi = self.hkb.index.values
            step = npahkbi[1]- npahkbi[0]
            nstart = npahkbi[0]+ (step * (offset))
            self.hkb.index = pd.Index(npahkbi + nstart)

            #add blank at begining
            df = pd.DataFrame({},columns=self.hkb.keys(),index=npahkbi[:offset])
            ndf=pd.concat([df,self.hkb])
            self.hkb=ndf
            self.thkb = self.hkb.index

    def _align_on_devdf(self,typ=''):
        """ align hkb or bs time on device data frame (devdf) time index

        In place (a.k.a. replace old self.hkb by the resampled one)

        Parameters
        ----------

        typp : 'HKB' |'BS'

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S = CorSer(6)
        >>> devdf = S.devdf[S.devdf['id']=='HKB:15']
        >>> hkbdf = S.hkb['AP1-AnkleLeft']
        >>> devdf2 = S._align_on_hkb(devdf,hkbdf,typ ='HKB')


        """
        if typ == 'HKB':
            idf = self.hkb
        elif typ == 'BS':
            idf = self.bespo
        elif typ == 'TCR':
            idf = self.tcr

        # mocap time
        #
        # 0 0.010001 0.020002
        mocapindex = pd.to_datetime(self.tmocap,unit='s')
        # radio time
        # 0  0.023  0.0473
        idf.index = pd.to_datetime(idf.index,unit='s')


        sf = (mocapindex[2]-mocapindex[1]).microseconds
        df = idf.resample(str(sf)+'U',fill_method='ffill')

        nindex = time2npa(df.index)
        df.index = pd.Index(nindex)

        if typ == 'HKB':
            self.hkb = df
        elif typ == 'BS':
            self.bespo = df
        elif typ == 'TCR':
            self.tcr = df

    def _align_devdf_on_hkb(self,devdf,hkbdf):
        """ NOT USED Practically
            align time of 2 data frames:

        the time delta of the second data frame is applyied on the first one
        (e.g. time for devdf donwsampled by hkb data frame time)


        Parameters
        ----------

        devdf : device dataframe
        hkbdf : hkbdataframe

        Return
        ------

        devdfc : 
            aligned copy device dataframe
        hkbdfc : 
            aligned copy hkbdataframe

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(6)
        >>> devdf = S.devdf[S.devdf['id']=='HKB:15']
        >>> hkbdf = S.hkb['AP1-AnkleLeft']
        >>> devdf2 = S._align_devdf_on_hkb(devdf,hkbdf)


        """

        devdfc=copy.deepcopy(devdf)

        hkbdfc=copy.deepcopy(hkbdf)

        idev = devdfc.index
        ihkb = hkbdfc.index

        devdfc.index = pd.to_datetime(idev,unit='s')
        hkbdfc.index = pd.to_datetime(ihkb,unit='s')
        sf = (hkbdfc.index[2]-hkbdfc.index[1]).microseconds

        # cannot resapmple devdf directly because multiple similar index values
        #need to resampl each groupby separately
        gb = devdfc.groupby(['id'])
        #get device id 
        devid,idevid = np.unique(devdfc['id'],return_index=True)
        #save corresponding subject to each device
        subject = {devid[i]:devdfc['subject'].iloc[i] for i in idevid}
        # resample each group separatley
        dgb={d:gb.get_group(d).resample(str(sf)+'U') for d in devid}
        #re insert subject and device id information in each resampled group
        for d in dgb:
            dgb[d]['subject']=subject[d]
            dgb[d]['id']=d

        # create the realigned Dataframe

        lgb = [dgb[d] for d in dgb]
        df = pd.concat(lgb)
        df.sort_index(inplace=True)

        nindex = time2npa(df.index)
        df.index = pd.Index(nindex)
        cols=['id','subject','x','y','z','v','vx','vy','vz','a','ax','ay','az']
        df=df[cols]

        return df

def rssi2dist(dist,Ssh,offsetdB=15):
    """
    Parameters
    ----------

    dist
    Ssh 
    offsetdB : float

    """
    if type(dist)==pd.DataFrame:    
        z1 = 10*np.log10((1/dist**2)).values
    else:
        z1 = 10*np.log10((1/dist**2))

    u = np.where(Ssh[0]==1)[0]
    z1[u] = z1[u]-offsetdB
    z1 = z1-np.mean(z1)
    return(z1)

#    z2 = Srssi.values
#    z2m = np.mean(z2[~np.isnan(z2)])
#    z2[np.isnan(z2)]=z2m
#    z2 = z2-np.mean(z2)
#    z1n = z1/np.sqrt(np.sum(z1*z1))
#    z2n = z2/np.sqrt(np.sum(z2*z2))
#    cn,dec,ratio = resync(z1n,z2n)
#    tdec.append(dec)

#tratio.append(ratio)
#if ratio > maxratio:
#    maxratio =  ratio

    # def get_data(self,a,b):


    #     T=self.tcr[a+'-'+b]
    #     T.name=T.name+'-tcr'
    #     H=self.hkb[a+'-'+b]
    #     H.name=H.name+'-hkb'
    #     udhk = self.accessdm(a,b,'HKB')
    #     udtcr = self.accessdm(a,b,'HKB')
    #     dist_tcr=self.dist[:,udtcr[0],udtcr[1]]
    #     dist_hkb=self.dist[:,udhk[0],udhk[1]]
    #     tdist=np.linspace(0,self.dist.shape[0]/100.,self.dist.shape[0])
    #     D_tcr=pd.Series(dist_tcr,index=tdist)
    #     D_tcr.name = 'dist-tcr'
    #     D_hkb=pd.Series(dist_hkb,index=tdist)
    #     D_hkb.name = 'dist-hkb'

    #     return T,H,D_tcr,D_hkb


    # def get_dataframes(self,a,b):
    #     """ assemble all series in a DataFrame
    #     """

    #     T,H,DT,DH = self.get_data(a,b)
    #     NH=(np.sqrt(1/(10**(H/10)))/4e4)
    #     NHc=NH-NH.mean()
    #     DHc=DH-DH.mean()
    #     inh = NHc.index
    #     idh = DHc.index
    #     NHc.index = pd.to_datetime(inh,unit='m')
    #     DHc.index = pd.to_datetime(idh,unit='m')
    #     sD = (DHc.index[1]-DHc.index[0])
    #     sf= str(int(sD.microseconds*1e-3)) + 'ms'
    #     NHcr = NHc.resample(sf,fill_method='ffill')
    #     return NHcr,DHc
