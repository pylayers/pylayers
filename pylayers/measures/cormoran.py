# -*- coding:Utf-8 -*-

import os
import pdb
import sys
import pandas as pd
import numpy as np
import numpy.ma as ma
import scipy.io as io
from pylayers.util.project import *
from pylayers.util.pyutil import *
from pylayers.mobility.ban.body import *
from pylayers.gis.layout import *
from matplotlib.widgets import Slider, CheckButtons, Button, Cursor
from pylayers.signal.DF import *

from moviepy.editor import *
from skimage import img_as_ubyte
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def cor_ls():
    filelog = os.path.join(os.environ['CORMORAN'],'RAW','Doc','MeasurementLog.csv')
    log = pd.read_csv(filelog)
    return log


class CorSer(PyLayers):
    """ Hikob data handling from CORMORAN measurement campaign 11/06/2014

    """

    def __init__(self,serie=6,day=11,source='UR1'):

        try:
            self.rootdir = os.environ['CORMORAN']
        except:
            raise NameError('Please add a CORMORAN environement variable \
                            pointing to the data')

        # infos
        self.serie = serie
        self.day = day
        self.loadlog()

        # Measures

        if day==11:
            self.stcr = [1,2,3,4,10,11,12,32,33,34,35,9,17,18,19,20,25,26]
            self.shkb = [5,6,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
            self.sbs  = [5,6,7,8,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
            self.mocap = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]
            self.interferers =[]
            self.mocapinterf=[]
        if day==12:
            self.stcr = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
            self.shkb = []#[9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
            self.sbs  = []#[9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
            self.mocap =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
            self.mocapinterf = [5,6,7,8,13,14,15,16,21,22,23,24,]

        if serie in self.shkb:
            self.loadhkb(serie=serie,day=day,source=source)

        if serie in self.stcr:
            self.loadTCR(serie=serie,day=day)

        if serie in self.sbs:
            self.loadBS(serie=serie,day=day)


        if self.typ=='FULL':
            self._filename = 'Sc' + self.scenario + '_S' + str(self.serie) + '_R' + str(self.run) + '_' + self.typ.capitalize()
        else:
            self._filename = 'Sc' + self.scenario + '_S' + str(self.serie) + '_R' + str(self.run) + '_' + self.typ


        # Layout
        self.L= Layout('MOCAP-small2.ini')
        #self.L= Layout('MOCAP.ini')
        #self.L= Layout('MOCAPext.ini')

        # Infrastructure Nodes
        self.loadinfranodes()
        self.loadcam()

        # BODY
        self.subject = str(self.log['Subject'].values[0]).split(' ')
        # filter typos in  self.subject
        self.subject = filter(lambda x : len(x)!=0,self.subject)
        if 'Jihad' in self.subject:
            uj = self.subject.index('Jihad')
            self.subject[uj]='Jihan'
        
        if serie in self.mocap :
            self.loadbody(serie=serie,day=day)
            self._distancematrix()
            self.computedevpdf()
            if isinstance(self.B,dict):
                for b in self.B:
                    self.B[b].traj.Lfilename=copy.copy(self.L.filename)
            else :
                self.B.traj.Lfilename=copy.copy(self.L.filename)

        self.title1 = 'Scenario:'+str(self.scenario)+' Serie:'+str(self.serie)+' Run:'+str(self.run)
        self.title2 = 'Type:'+str(self.typ)+ ' Subject:'+str(self.subject[0])

        self.offset={'video':[],
                     'hkb':[],
                     'tcr':[],
                     'bs':[]
                     }

    def __repr__(self):
        st = ''
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


    def loadcam(self):

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


    def loadinfranodes(self):
        """ load infrastrucutre nodes



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
        a,self.infraname,pts,i = c3d.ReadC3d(filename)

        pts = pts/1000.
        mpts = np.mean(pts,axis=0)
        self.din={}
        if ('HK'  in self.typ) or ('FULL' in self.typ):
            uhkb = np.array([[1,2],[4,5],[7,8],[10,11]])
            mphkb = np.mean(mpts[uhkb],axis=1)

            self.din.update(
                {'HKB:1':{'p':mphkb[3],
                          'T':np.eye(3)},
                 'HKB:2':{'p':mphkb[2],
                          'T': np.array([[-0.44807362,  0.89399666,  0.],
                                         [-0.89399666, -0.44807362,  0.],
                                         [ 0.,0.,1.        ]])}      ,
                 'HKB:3':{'p':mphkb[1],
                          'T':array([[-0.59846007, -0.80115264,  0.],
                                     [ 0.80115264, -0.59846007,  0.],
                                     [ 0.,0.,  1.]])},
                 'HKB:4':{'p':mphkb[0],
                          'T':array([[-0.44807362, -0.89399666,  0.],
                                     [ 0.89399666, -0.44807362,  0.],
                                     [ 0.,0.,  1.]])}
                 })

        if ('TCR' in self.typ) or ('FULL' in self.typ):
            self.din.update({'TCR:32':{'p':mpts[9],
                                       'T':np.eye(3)},
                 'TCR:24':{'p':mpts[6],
                           'T': np.array([[-0.44807362,  0.89399666,  0.],
                                         [-0.89399666, -0.44807362,  0.],
                                         [ 0.,0.,1.        ]])},
                 'TCR:27':{'p':mpts[3],
                           'T':array([[-0.59846007, -0.80115264,  0.],
                                     [ 0.80115264, -0.59846007,  0.],
                                     [ 0.,0.,  1.]])},
                 'TCR:28':{'p':mpts[0],
                           'T':array([[-0.44807362, -0.89399666,  0.],
                                     [ 0.89399666, -0.44807362,  0.],
                                     [ 0.,0.,  1.]])}
                 })

        # self.pts= np.empty((12,3))
        # self.pts[:,0]= -mpts[:,1]
        # self.pts[:,1]= mpts[:,0]
        # self.pts[:,2]= mpts[:,2]
        # return mpts
        # self.dist = np.sqrt(np.sum((mpts[:,np.newaxis,:]-mpts[np.newaxis,:])**2,axis=2))


    def loadlog(self):
        """ load in self.log the log of current serie
            from MeasurementLog.csv
        """

        filelog = os.path.join(self.rootdir,'RAW','Doc','MeasurementLog.csv')
        log = pd.read_csv(filelog)
        date = str(self.day)+'/06/14'
        self.log = log[(log['Meas Serie'] == self.serie) & (log['Date'] == date)]


    def loadbody(self,day=11,serie=''):
        """ load log file
        """
        self.B={}
        color=['LightBlue','YellowGreen','PaleVioletRed','white','white','white','white','white','white','white']
        import ipdb
        ipdb.set_trace()
        for us,subject in enumerate(self.subject):
            print "load ",subject, " body\n"
            seriestr = str(self.serie).zfill(3)
            if day == 11:
                filemocap = os.path.join(self.rootdir,'RAW',str(self.day)+'-06-2014','MOCAP','serie_'+seriestr+'.c3d')
            elif day == 12:
                filemocap = os.path.join(self.rootdir,'RAW',str(self.day)+'-06-2014','MOCAP','Nav_serie_'+seriestr+'.c3d')
            baw = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','BodyandWear')
            if subject =='Jihad':
                subject ='Jihan'
            filebody = os.path.join(baw, subject + '.ini')
            filewear = os.path.join(baw,subject + '_'  +str(self.day)+'-06-2014_' + self.typ + '.ini')

            self.B.update({subject:Body(_filebody=filebody,
                             _filemocap=filemocap,unit = 'mm', loop=False,
                             _filewear=filewear,
                             centered=False,
                             multi_subject_mocap=True,
                             color=color[us])})

        if self.serie in self.mocapinterf:
            self.interf = ['Anis_Cylindre:',
                     'Benoit_Cylindre:',
                     'Bernard_Cylindre:',
                     'Claude_Cylindre:',
                     'Meriem_Cylindre:']

            for i in self.interf:
                print "load ",i, " interfering body"
                self.B.update({i:Cylinder(name=i,_filemocap=filemocap,unit = 'mm')})

        # if len(self.subject) == 1:
        #     self.B = self.B[self.subject]


    def loadTCR(self,day=11,serie='',scenario='20',run=1):
        """ load TCR data

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
                  'AnckleRight':36,
                  'AnckleLeft':37}
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

    def loadBS(self,day=11,serie='',scenario='20',run=1):
        """ load BeSpoon data

        """
        self.dBS = {'LeftWrist?':157,'RightAnckle?':74}
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

        self.bespo = pd.read_csv(os.path.join(dirname,self._fileBS),index_col='tu')
        self.s157 = self.bespo[self.bespo['Sensor']==157]
        #self.s157.set_index(self.s157['tu'].values/1e9)
        self.s74  = self.bespo[self.bespo['Sensor']==74]
        #self.s74.set_index(self.s74['tu'].values/1e9)
        #t157 = np.array(self.s157['tu']/(1e9))
        #self.t157 = t157-t157[0]
        #t74 = np.array(self.s74['tu']/(1e9))
        #self.t74 = t74 - t74[0]


    def loadhkb(self,day=11,serie='',scenario='20',run=1,source='UR1'):

        if day==11:
            self.dHKB ={'AP1':1,'AP2':2,'AP3':3,'AP4':4,
                       'HeadRight':5,'TorsoTopRight':6,'TorsoTopLeft':7,'BackCenter':8,'ElbowRight':9,'ElbowLeft':10,'HipRight':11,'WristRight':12,'WristLeft':13,'KneeLeft':14,'AnckleRight':16,'AnckleLeft':15}
            if source=='UR1':
                dirname = os.path.join(self.rootdir,'POST-TREATED','11-06-2014','HIKOB')
            elif source=='CITI':
                dirname = os.path.join(self.rootdir,'POST-TREATED','11-06-2014','HIKOB','CITI')
        if day==12:
            self.dHKB= {'AP1':1,'AP2':2,'AP3':3,'AP4':4,'Jihad:TorsoTopRight':10,'Jihad:TorsoTopLeft':9,'Jihad:BackCenter':11,'JihadShoulderLeft':12,
             'Nicolas:TorsoTopRight':6,'Nicolas:TorsoTopLeft':5,'Nicolas:BackCenter':7,'Nicolas:ShoulderLeft':8,
             'Eric:TooTopRight':15,'Eric:TorsoTopLeft':13,'Eric:BackCenter':16,'Eric:ShoulderLeft':14}
            if source=='UR1':
                dirname = os.path.join(self.rootdir,'POST-TREATED','12-06-2014','HIKOB')
            elif source=='CITI':
                dirname = os.path.join(self.rootdir,'POST-TREATED','12-06-2014','HIKOB','CITI')

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
                self.scenario=tt[0].replace('Sc','')
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

        self.topandas()
        self.hkb = self.hkb[self.hkb!=0]


    def _distancematrix(self):
        """Compute the ditance matrix between the nodes

            self.dist : (nb frame x nb_node x nb_node)
            self.dist_nodes : list of used nodes (useful to make the association ;) )
        """


        # tdev = []
        # for k in self.B.dev:
        #     tdev.append(self.B.dev[k]['uc3d'][0])
        # tdev = np.array(tdev)

        # # pnb : nframe x ndevices x 3
        # pnb = self.B._f[:,tdev,:]

        # ln = []
        # uin = []
        # if ('HK' in self.typ) or ('FULL' in self.typ):
        #     uin.extend(['HKB:1','HKB:2','HKB:3','HKB:4'])
        # if ('TCR' in self.typ) or ('FULL' in self.typ):
        #     uin.extend(['TCR:32','TCR:24','TCR:27','TCR:28'])
        # ln = uin + self.B.dev.keys()
        # pin = np.array([self.din[d]['p'] for d in uin])
        # pin2=np.empty((pnb.shape[0],pin.shape[0],pin.shape[1]))
        # pin2[:,:,:]=pin
        # p = np.concatenate((pin2,pnb),axis=1)
        # self.dist = np.sqrt(np.sum((p[:,:,np.newaxis,:]-p[:,np.newaxis,:,:])**2,axis=3))
        # self.dist_nodes = ln

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
        if ('HK' in self.typ) or ('FULL' in self.typ):
            uin.extend(['HKB:1','HKB:2','HKB:3','HKB:4'])
        if ('TCR' in self.typ) or ('FULL' in self.typ):
            uin.extend(['TCR:32','TCR:24','TCR:27','TCR:28'])
        ln = uin + bn
        pin = np.array([self.din[d]['p'] for d in uin])
        pin2=np.empty((pnb.shape[0],pin.shape[0],pin.shape[1]))
        pin2[:,:,:]=pin
        p = np.concatenate((pin2,pnb),axis=1)
        self.dist = np.sqrt(np.sum((p[:,:,np.newaxis,:]-p[:,np.newaxis,:,:])**2,axis=3))
        self.dist_nodes = ln

    def accessdm(self,a,b,techno):
        """ access to the distance matrix

            give name|id of node a and b and a given techno. retrun Groung truth
            distance between the 2 nodes
        """

        if 'HKB' in techno :
            if isinstance(a,str):
                ia = self.dHKB[a]
            else:
                ia = a
                a = self.idHKB[a]

            if isinstance(b,str):
                ib = self.dHKB[b]
            else:
                ib = b
                b = self.idHKB[b]

        elif 'TCR' in techno :
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

        else :
            raise AttributeError('please give only 1 techno or radio node')

        ka = techno+':'+str(ia)
        kb = techno+':'+str(ib)

        ua = self.dist_nodes.index(ka)
        ub = self.dist_nodes.index(kb)

        return(ua,ub)





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
            filename = videofile+_filename
            os.system('vlc '+filename +'&' )
        except:
            raise AttributeError('file '+ self._filename + ' not found')


    def snapshot(self,t0=0,offset=15.5,title=True,save=False,fig=[],ax=[],figsize=(10,10)):
        """ single snapshot plot
        """

        if fig ==[]:
            fig=plt.figure(figsize=figsize)
        if ax == []:
            ax = fig.add_subplot(111)

        videofile = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','Videos')
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = videofile+_filename
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
        """
        """
        videofile = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','Videos')
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = videofile+_filename
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
        devsize : float
            device on body size (100)
        bodytime: list
            list of time instant where body topos has to be shown

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
        incolor: str
            color of infrastructure nodes ('r')
        insize
            size of infrastrucutre nodes (0.1)


        camera : boolean
            display Vicon camera position (True)
        cameracolor : str
            color of camera nodes ('b')
        camerasize  : float
            size of camera nodes (0.1)




        """
        defaults = { 'L':True,
                     'body':True,
                     'subject':[],
                     'trajectory' :True,
                     'devsize':100,
                     'inodes' : True,
                     'inname' : True,
                     'incolor' : 'r',
                     'insize' : 0.1,
                     'camera':True,
                     'cameracolor' :'k',
                     'camerasize' :0.1,
                     'bodytime':[],
                     'tagtraj':True,
                     'tagname':[],
                     'tagpoffset':[],
                     'fontsizetag':0.5,
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
        if kwargs['inname']:
            [mlab.text3d(v[i][1]['p'][0],v[i][1]['p'][1],v[i][1]['p'][2],v[i][0],scale=0.5)
            for i in range(len(v))]
        if kwargs['body']:

            if kwargs['bodytime']==[]:
                time =np.linspace(0,self.B[subject[0]].time[-1],5).astype(int)
                # time=range(10,100,20)
            else :
                time=kwargs['bodytime']
            for ki, i in enumerate(time):
                for ib,b in enumerate(subject):
                    self.B[b].settopos(t=i,cs=True)
                    self.B[b]._show3(dev=True,devsize=kwargs['devsize'])
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

    def topandas(self):
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


    def imshow(self,time=100,kind='time'):
        """

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
            ib = b
            b = self.idHKB[b]


        time = self.thkb[0]
        if len(time) == 1:
            time=time[0]



        sab = self.hkb[a+'-'+b].values
        sabt = self.hkb[a+'-'+b].index
        hkb = axs[1].plot(sabt,sab)



        videofile = os.path.join(self.rootdir,'POST-TREATED',str(self.day)+'-06-2014','Videos')
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = videofile+_filename
        vc = VideoFileClip(filename)
        F0 = vc.get_frame(0)
        I0 = img_as_ubyte(F0)
        axs[0].imshow(F0)


        ########
        # slider
        ########
        slide_xoffset_ax = plt.axes([0.1, 0.15, 0.8, 0.05])
        sliderx = Slider(slide_xoffset_ax, "video offset", 0, self.hkb.index[-1],
                        valinit=time[0], color='#AAAAAA')


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
            self.offset['video']=sliderx.val

        axp = plt.axes([0.3, 0.05, 0.1, 0.075])
        axset = plt.axes([0.5, 0.05, 0.1, 0.075])
        axm = plt.axes([0.7, 0.05, 0.1, 0.075])

        bp = Button(axp, '<-')
        bp.on_clicked(minus)

        bs = Button(axset, 'SET offs.')
        bs.on_clicked(setter)

        bm = Button(axm, '->')
        bm.on_clicked(plus)

        plt.show()






    def offset_setter_hkb(self,a,b,**kwargs):
        """ offset setter
        """
        defaults = { 'inverse':True
                    }


        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.2, left=0.3)

        if isinstance(a,str):
            ia = self.dHKB[a]
        else:
            ia = a
            a = self.idHKB[a]

        if isinstance(b,str):
            ib = self.dHKB[b]
        else:
            ib = b
            b = self.idHKB[b]


        time = self.thkb[0]
        if len(time) == 1:
            time=time[0]


        dhk = self.accessdm(ia,ib,'HKB')
        if kwargs['inverse']:
            var = 10*np.log10(1./(self.dist[:,dhk[0],dhk[1]])**2)
        else :
            var = self.dist[:,dhk[0],dhk[1]]
        gt = ax.plot(self.B[self.B.keys()[0]].time,var)

        sab = self.hkb[a+'-'+b].values
        sabt = self.hkb[a+'-'+b].index
        hkb = ax.plot(sabt,sab)


        ########
        # slider
        ########
        slide_xoffset_ax = plt.axes([0.1, 0.15, 0.8, 0.02])
        sliderx = Slider(slide_xoffset_ax, "hkb offset", -(len(sabt)/16), (len(sabt)/16),
                        valinit=time[0], color='#AAAAAA')

        slide_yoffset_ax = plt.axes([0.1, 0.10, 0.8, 0.02])
        slidery = Slider(slide_yoffset_ax, "gt_yoff", -100, 0,
                        valinit=0, color='#AAAAAA')

        slide_alpha_ax = plt.axes([0.1, 0.05, 0.8, 0.02])
        slideralpha = Slider(slide_alpha_ax, "gt_alpha", 0, 10,
                        valinit=0, color='#AAAAAA')

        def update_x(val):
            value = int(sliderx.val)
            rhkb = np.roll(sab,value)
            sliderx.valtext.set_text('{}'.format(value))
            hkb[0].set_xdata(sabt)
            hkb[0].set_ydata(rhkb)
            fig.canvas.draw_idle()
        sliderx.on_changed(update_x)
        sliderx.drawon = False


        def update_y(val):
            yoff = slidery.val
            alpha = slideralpha.val
            gt[0].set_ydata(alpha*var + yoff)
            fig.canvas.draw_idle()
        slidery.on_changed(update_y)
        slideralpha.on_changed(update_y)

        plt.show()


    def pltvisi(self,a,b,**kwargs):
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
        >>> f,ax = S.plthkb('AP1','TorsoTopLeft')
        >>> f,ax = S.pltvisi('AP1','TorsoTopLeft',fig=f,ax=ax)
        >>> #f,ax = S.pltmob(showvel=False,ylim=([-100,-40]),fig=f,ax=ax)
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
        vv,tv,tseg,itseg = self.visiarray(a,b)
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
        >>> f,ax = S.plthkb('AP1','TorsoTopLeft')
        >>> #f,ax = S.pltvisi('AP1','TorsoTopLeft',fig=f,ax=ax)
        >>> f,ax = S.pltmob(showvel=False,ylim=([-100,-40]),fig=f,ax=ax)
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
                     'ylim':(-200,0),
                     'time_offset':0,
                     'color':'gray',
                     'hatch':'',
                     'label_pos':'top',
                     'label_pos_off':5,
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

        fig , ax =plu.rectplot(self.B[subject].time,nullr,ylim=kwargs['ylim'],
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
                yposM = kwargs['ylim'][1]-kwargs['label_pos_off']+0.5
                yposS = kwargs['ylim'][1]-kwargs['label_pos_off']-0.5

            elif kwargs['label_pos'] == 'bottom':
                yposM = kwargs['ylim'][0]+kwargs['label_pos_off']+0.5
                yposS = kwargs['ylim'][0]+kwargs['label_pos_off']+0.5
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

        a : node name | number
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


    def plthkb(self,a,b,**kwargs):
        """
        Parameters
        ----------

        a : node name | number
        b : node name | number
        t0 : start time
        t1 : stop time


        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S = CorSer(6)
        >>> f,ax = S.plthkb('AP1','TorsoTopLeft')
        >>> f,ax = S.pltvisi('AP1','TorsoTopLeft',fig=f,ax=ax)
        >>> f,ax = S.pltmob(showvel=False,ylim=([-100,-40]),fig=f,ax=ax)
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

        if isinstance(a,str):
            ia = self.dHKB[a]
        else:
            ia = a
            a = self.idHKB[a]

        if isinstance(b,str):
            ib = self.dHKB[b]
        else:
            ib = b
            b = self.idHKB[b]

        if kwargs['shortlabel']:

            #find uppercase position
            uu =  np.nonzero([l.isupper() or l.isdigit() for l in a])[0]
            # cretae string from list
            labela = ''.join([a[i] for i in uu])

            uu =  np.nonzero([l.isupper() or l.isdigit() for l in b])[0]
            # cretae string from list
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

        a : node name | number
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
            t1=self.thkb[0][-1]

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
        log : boolean
        gammma':1.,
        mode : string
            'HKB' | 'TCR' | 'FULL'
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

        if kwargs['subject']=='':
            subject=self.B.keys()[0]
        else:
            subject=kwargs['subject']


        t0 =kwargs.pop('t0')
        t1 =kwargs.pop('t1')
        if t1 ==-1:
            t1=self.thkb[0][-1]


        label = a+'-'+b

        mode = kwargs.pop('mode')
        inverse = kwargs.pop('inverse')
        log = kwargs.pop('log')
        gamma = kwargs.pop('gamma')
        visibility = kwargs.pop('visi')
        fontsize = kwargs.pop('fontsize')
        hatch = kwargs.pop('hatch')



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

            dhk = self.accessdm(iahk,ibhk,'HKB')

            if inverse:
                var = 1./(self.dist[:,dhk[0],dhk[1]])
                ax.set_ylabel(u'$m^{-2}$',fontsize=fontsize)
                if log :
                    #var = gamma*10*np.log10(var)
                    var = 20*np.log10(var)+gamma
                    ax.set_ylabel(u'$- 20 \log_{10}(d)'+str(gamma)+'$  (dB)',fontsize=fontsize)
                    plt.ylim(-65,-40)
            else:
                var = self.dist[:,dhk[0],dhk[1]]
                ax.set_ylabel(u'meters',fontsize=fontsize)
                if log :
                    var = gamma*10*np.log10(var)+gamma
                    ax.set_ylabel(u'$10log_{10}m^{-2}$',fontsize=fontsize)


            ax.plot(self.B[subject].time,var,label=label,**kwargs)
        #
        # TCR | Full
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

            dtcr = self.accessdm(iatcr,ibtcr,'TCR')
            if inverse:
                var = 1./(self.dist[:,dtcr[0],dtcr[1]])**2
                if log :
                    var = gamma*10*np.log10(var)
            else:
                var = self.dist[:,dtcr[0],dtcr[1]]
                if log :
                    var = gamma*10*np.log10(var)
            ax.plot(self.B[subject].time,var,**kwargs)


        if visibility:
            aa= ax.axis()
            vv,tv,tseg,itseg = self.visiarray(a,b)
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

        # Axes management
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


        # TCR plot
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
            # plot only if link exist
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
        # HKB | Full
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


    def showlink(self,a,b,technoa='HKB',technob='HKB',iframe=0,style='*b'):
        """ show link configuation for a given frame

        Parameters
        ----------

        a
        b
        technoa
        technob
        iframe
        style

        """
        # display nodes
        A,B = self.getdevp(a,b,technoa=technoa,technob=technob)
        a,ia,ba,subjecta = self.devmapper(a,technoa)
        b,ib,bb,subjectb = self.devmapper(b,technob)

        if A.ndim==2:
            plt.plot(A[iframe,0],A[iframe,1],'ob')
            plt.text(A[iframe,0],A[iframe,1],a)
        else:
            plt.plot(A[0],A[1],'or')
            #plt.text(A[0],A[1],a)

        if B.ndim==2:
            plt.plot(B[iframe,0],B[iframe,1],style)
            plt.text(B[iframe,0]+0.1,B[iframe,1]+0.1,b)
        else:
            plt.plot(B[0],B[1],'ob')
            plt.text(B[0],B[1],b)
        plt.xlim(-6,6)
        plt.ylim(-5,5)
        # display body

        #pc = self.B.d[:,2,iframe] + self.B.pg[:,iframe].T
        pc0 = self.B[subjecta].d[:,0,iframe] + self.B[subjecta].pg[:,iframe].T
        pc1 = self.B[subjecta].d[:,1,iframe] + self.B[subjecta].pg[:,iframe].T
        pc15 = self.B[subjecta].d[:,15,iframe] + self.B[subjecta].pg[:,iframe].T
        #plt.plot(pc0[0],pc0[1],'og')
        #plt.text(pc0[0]+0.1,pc0[1],str(iframe))
        #plt.plot(pc1[0],pc1[1],'og')
        #plt.plot(pc15[0],pc15[1],'og')
        #ci00   = plt.Circle((pc0[0],pc0[1]),self.B[subjecta].sl[0,2],color='green',alpha=0.6)
        #ci01   = plt.Circle((pc1[0],pc1[1]),self.B[subjecta].sl[0,2],color='green',alpha=0.1)
        #ci100 = plt.Circle((pc0[0],pc0[1]),self.B[subjecta].sl[10,2],color='red',alpha=0.1)
        ci1015 = plt.Circle((pc15[0],pc15[1]),self.B[subjecta].sl[10,2],color='green',alpha=0.5)
        plt.axis('equal')
        ax = plt.gca()
        ax.add_patch(ci1015)
        #ax.add_patch(ci01)
        #ax.add_patch(ci100)
        #ax.add_patch(ci1015)
        #its = self.B[subjecta].intersectBody(A[iframe,:],B[iframe,:],topos=False,frameId=iframe)
        #x.set_title('frameId :'+str(iframe)+' '+str(its.T))


    def visidev(self,a,b,technoa='HKB',technob='HKB',dsf=10):
        """ get link visibility status

        Returns
        -------

        visi : pandas Series
            0  : LOS
            1  : NLOS

        """

        A,B = self.getdevp(a,b,technoa,technob)
        aa,ia,ba,subjecta= self.devmapper(a,technoa)
        ab,ib,bb,subjectb= self.devmapper(b,technob)


        if 'AP' not in a:
            Nframe = A.shape[0]
        if 'AP' not in b:
            Nframe = B.shape[0]
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

        A,B = self.getdevp(a,b,technoa,technob)
        aa,ia,ba,subjecta= self.devmapper(a,technoa)
        ab,ib,bb,subjectb= self.devmapper(b,technob)

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



    def visiarray(self,a,b,technoa='HKB',technob='HKB'):
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

        itseg = copy.copy(tseg)
        bb = np.insert(itseg[:,1],0,0)
        ee = np.hstack((itseg[:,0],len(vv)))
        itseg = np.array((bb,ee)).T
        # bb = np.hstack((bb,len(vv)))

        return vv,tv,tseg,itseg


    # def computedevpdf(self):
    #     """ create a timestamped data frame
    #         with all positions of devices
    #     """
    #     t=self.B.traj.time()
    #     pos = np.empty((len(t),12,3))
    #     for ik,k in enumerate(t):
    #         self.B.settopos(t=k)
    #         pos[ik,:,:]=self.B.getdevp()
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

    def computedevpdf(self):
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
                udev=[B[b].dev[d]['uc3d'][0] for d in dev]
                pos = B[b]._f[:,udev,:]
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


        df = df.sort_index()
        cols=['id','subject','x','y','z','v','vx','vy','vz','a','ax','ay','az']
        self.devdf=df[cols]



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

        Return
        ------

            a csv file into the folder <PylayersProject>/netsave

        """

        defaults={'unit' :'mm',
                  'tunit':'ns',
                  'alias':{}}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        unit=kwargs.pop('unit')
        tunit=kwargs.pop('tunit')
        alias = kwargs.pop('alias')

        if alias == {}:

            alias={'TCR:49':4 # Nicolas TorsoLeft
            ,'TCR:34':5 # Nicolas TorsoRight
            ,'TCR:48':6 # Nicolas Back
            ,'TCR:36':7 # Nicolas Shoulder

            ,'TCR:2':8 # Jihad TorsoLeft
            ,'TCR:35':9 #Jihad TorsoRight
            ,'TCR:33':10 #Jihad Back
            ,'TCR:37':11 #Jihad Shoulder

            ,'TCR:30':12 # Eric Torso
            ,'TCR:25':13 # Eric Back
            ,'TCR:26':14 # Eric Shoulder
            }

        filename =pyu.getlong(self._filename,pstruc['DIRNETSAVE']) + '.csv'

        df = copy.deepcopy(self.devdf)

        ldf = df[['id','x','y','z']]

        # rename devices
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

        ldf.loc[:,'x']=ldf.loc[:,'x']*_unit
        ldf.loc[:,'y']=ldf.loc[:,'y']*_unit
        ldf.loc[:,'z']=ldf.loc[:,'z']*_unit

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


    def getdevp(self,a,b,technoa='HKB',technob='HKB'):
        """    get device position

        Parameters
        ----------

        a : str | int
            name | id
        b : str | int
            name | id
        technoa : str
            radio techno
        technob : str
            radio techno

        Returns
        -------

        pa,pb : np.array()

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(serie=34)
        >>> a,b=S.getdevp('AP1','WristLeft')

        """

        a,ia,nna,subjecta = self.devmapper(a,technoa)
        b,ib,nnb,subjectb = self.devmapper(b,technob)

        # node a
        # body node

        if subjecta != '':
            unna = self.B[subjecta].dev[nna]['uc3d'][0]
            pa = self.B[subjecta]._f[:,unna,:]
        # infra node
        else :
            pa = self.din[nna]['p']


        # node b
        # body node
        if subjectb != '':
            unnb = self.B[subjectb].dev[nnb]['uc3d'][0]
            pb = self.B[subjectb]._f[:,unnb,:]
        # infra node
        else :
            pb = self.din[nnb]['p']

        return pa,pb


    def devmapper(self,a,techno='HKB'):
        """  retrieve name of device if input is number
             or
             retrieve number of device if input is name

        Parameters
        ----------

        a : str | int
            name | id
        techno : str
            radio techno


        Return
        ------

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
        if isinstance(a,str):
            if techno == 'TCR':
                ia = self.dTCR[a]
                nna='TCR:'+str(ia)
            elif techno == 'HKB':
                ia = self.dHKB[a]
                ba='HKB:'+str(ia)


            for b in self.B:
                if ba in self.B[b].dev.keys():
                    subject = b
                    break

        else:
            if techno == 'TCR':
                ia = a
                a = self.idTCR[a]
                nna='TCR:'+str(ia)
            elif techno == 'HKB':
                ia = a
                a = self.idHKB[a]
                ba='HKB:'+str(ia)

            for b in self.B:
                if ba in self.B[b].dev.keys():
                    subject = b
                    break

        return a,ia,ba,subject


    def align(self,devdf,hkbdf):

        """ align time for device and hkb data frame

        devdf : device dataframe
        hkbdf : hkbdataframe

        Examples
        --------

        >>> from pylayers.measures.cormoran import *
        >>> S=CorSer(6)
        >>> devdf = S.devdf[S.devdf['id']=='HKB:15']
        >>> hkbdf = S.hkb['AP1-AnckleLeft']
        >>> devdf2,hkbdf2 = S.align(devdf,hkbdf)


        """


        idev = devdf.index
        ihkb = hkbdf.index
        devdf.index = pd.to_datetime(idev,unit='s')
        hkbdf.index = pd.to_datetime(ihkb,unit='s')
        sf = (hkbdf.index[2]-hkbdf.index[1]).microseconds
        devdf= devdf.resample(str(sf)+'U')
        return devdf,hkbdf

    def get_data(self,a,b):


        T=self.tcr[a+'-'+b]
        T.name=T.name+'-tcr'
        H=self.hkb[a+'-'+b]
        H.name=H.name+'-hkb'
        udhk = self.accessdm(a,b,'HKB')
        udtcr = self.accessdm(a,b,'HKB')
        dist_tcr=self.dist[:,udtcr[0],udtcr[1]]
        dist_hkb=self.dist[:,udhk[0],udhk[1]]
        tdist=np.linspace(0,self.dist.shape[0]/100.,self.dist.shape[0])
        D_tcr=pd.Series(dist_tcr,index=tdist)
        D_tcr.name = 'dist-tcr'
        D_hkb=pd.Series(dist_hkb,index=tdist)
        D_hkb.name = 'dist-hkb'

        return T,H,D_tcr,D_hkb


    def get_dataframes(self,a,b):
        """ assemble all series in a DataFrame
        """

        T,H,DT,DH = self.get_data(a,b)
        NH=(np.sqrt(1/(10**(H/10)))/4e4)
        NHc=NH-NH.mean()
        DHc=DH-DH.mean()
        inh = NHc.index
        idh = DHc.index
        NHc.index = pd.to_datetime(inh,unit='m')
        DHc.index = pd.to_datetime(idh,unit='m')
        sD = (DHc.index[1]-DHc.index[0])
        sf= str(int(sD.microseconds*1e-3)) + 'ms'
        NHcr = NHc.resample(sf,fill_method='ffill')
        return NHcr,DHc
