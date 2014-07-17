# -*- coding:Utf-8 -*-

import os
import pdb
import sys
import pandas as pd
import numpy as np
import numpy.ma as ma
import scipy.io as io
from pylayers.util.project import *
from pylayers.mobility.ban.body import *
from pylayers.gis.layout import *

from moviepy.editor import *
from skimage import img_as_ubyte
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class CorSer(PyLayers):
    """ Hikob data handling from CORMORAN measurement campaign 11/06/2014

    """

    def __init__(self,serie=6,day=11,source='UR1'):

        try:
            self.rootdir =os.environ['CORMORAN']
        except:
            raise NameError('Please add a CORMORAN environement variable \
                            pointing to the data')

        # infos
        self.serie = serie
        self.day = day
        self.loadlog()

        # Measures

        if day==11:
            stcr = [1,2,3,4,10,11,12,32,33,34,35,9,17,18,19,20,25,26]
            shkb = [5,6,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
            sbs  = [5,6,7,8,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
            mocap = [5,6,7,8,17,21,22,23,24,34]

        if day==12:
            stcr = []
            shkb = []
            sbs  = []

        if serie in shkb:
            self.loadhkb(serie=serie,day=day,source=source)

        if serie in stcr:
            self.loadTCR(serie=serie,day=day)

        if serie in sbs:
            self.loadBS(serie=serie,day=day)

        self._filename = 'Sc' + self.scenario + '_S' + str(self.serie) + '_R' + str(self.run) + '_' + self.typ.capitalize()


        # BODY
        self.subject = [str(self.log['Subject'].values[0])]
        if serie in mocap :
            self.loadbody(serie=serie,day=day)


        # Layout
        self.L= Layout('MOCAP.ini')

        self.loadinfranodes()
        self._distancematrix()


    def __repr__(self):
        st = ''
        st = st + 'Day : '+ str(self.day)+'/07/2014'+'\n'
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


    def loadinfranodes(self):
        """ load infrastrucutre nodes



                        A4 
                    mpts[6,7,8]
                        X 

            A3                     A3 
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


        """
        a,self.infraname,pts,i = c3d.ReadC3d('scene.c3d')

        pts=pts/1000.
        mpts = np.mean(pts,axis=0)
        self.din={}
        if ('HK'  in self.typ) or ('FULL' in self.typ):
            uhkb=np.array([[1,2],[4,5],[7,8],[10,11]])
            mphkb = np.mean(mpts[uhkb],axis=1)

            self.din.update({'HKB:1':mphkb[1],
                 'HKB:2':mphkb[0],
                 'HKB:3':mphkb[3],
                 'HKB:4':mphkb[2]})
        if ('TCR' in self.typ) or ('FULL' in self.typ):
            self.din.update({'TCR:32':mpts[3],
                 'TCR:24':mpts[0],
                 'TCR:27':mpts[9],
                 'TCR:28':mpts[6]})
        
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

        filelog = self.rootdir + '/RAW/Doc/MeasurementLog.csv'
        log = pd.read_csv(filelog)
        date = str(self.day)+'/06/14'
        self.log = log[(log['Meas Serie'] == self.serie) & (log['Date'] == date)]


    def loadbody(self,day=11,serie=''):
        """ load log file
        """
        self.B=[]
        for subject in self.subject:
            
            seriestr = str(self.serie).zfill(3)
            filemocap = self.rootdir + '/RAW/' + \
                        str(self.day)+'-06-2014/MOCAP/serie_' + seriestr + '.c3d'
            baw = self.rootdir + '/POST-TREATED/' + str(self.day)+'-06-2014/BodyandWear/'
            filebody = baw + subject + '.ini'
            filewear = baw + subject + '_'  +str(self.day)+'-06-2014_' + self.typ + '.ini'
            self.B.append(Body(_filebody=filebody,
                             _filemocap=filemocap,unit = 'mm',
                             _filewear=filewear))

        if len(self.subject) == 1:
            self.B = self.B[0]


        

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
            dirname = self.rootdir+'/POST-TREATED/11-06-2014/TCR'


        if day==12:
            dirname = self.rootdir+'/POST-TREATED/12-06-2014/TCR'
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

        filename = dirname + '/'+ self._fileTCR
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

    def loadBS(self,day=11,serie='',scenario='20',run=1):
        """ load BeSpoon data

        """
        self.dBS = {'LeftWrist?':157,'RightAnckle?':74}
        if day==11:
            dirname = self.rootdir+'/POST-TREATED/11-06-2014/BeSpoon'
        if day==12:
            dirname = self.rootdir+'/POST-TREATED/12-06-2014/BeSpoon'

        files = os.listdir(dirname)
        if serie != '':
            self._fileBS = filter(lambda x : 'S'+str(serie) in x ,files)[0]
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            self._fileBS = filter(lambda x : 'R'+str(run) in x ,filsc)[0]

        self.bespo = pd.read_csv(dirname+'/'+self._fileBS,index_col='tu')
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
                dirname = self.rootdir+'/POST-TREATED/11-06-2014/HIKOB'
            elif source=='CITI':
                dirname = self.rootdir+'/POST-TREATED/11-06-2014/HIKOB/CITI'
        if day==12:
            self.dHKB= {'AP1':1,'AP2':2,'AP3':3,'AP4':4,'Jihad:TorsoTopRight':10,'Jihad:TorsoTopLeft':9,'Jihad:BackCenter':11,'JihadShoulderLeft':12,
             'Nicolas:TorsoTopRight':6,'Nicolas:TorsoTopLeft':5,'Nicolas:BackCenter':7,'Nicolas:ShoulderLeft':8,
             'Eric:TorsoTopRight':15,'Eric:TorsoTopLeft':13,'Eric:BackCenter':16,'Eric:ShoulderLeft':14}
            if source=='UR1':
                dirname = self.rootdir+'/POST-TREATED/12-06-2014/HIKOB'
            elif source=='CITI':
                dirname = self.rootdir+'/POST-TREATED/12-06-2014/HIKOB/CITI'

        files = os.listdir(dirname)

        self.idHKB={}
        for k in self.dHKB:
            self.idHKB[self.dHKB[k]]=k

        if serie != '':
            self._filehkb = filter(lambda x : 'S'+str(serie) in x ,files)[0]
            tt = self._filehkb.split('_')
            self.scenario=tt[0].replace('Sc','')
            self.run = tt[2].replace('R','')
            self.typ = tt[3]
            self.video = tt[4].replace('.mat','')
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            if source=='UR1':
                self._filehkb = filter(lambda x : 'R'+str(run) in x ,filsc)[0]
            else:
                self._filehkb = filter(lambda x : 'r'+str(run) in x ,filsc)[0]


        data = io.loadmat(dirname+'/'+self._filehkb)
        if source=='UR1':
            self.rssi = data['rssi']
            self.t = data['t']
        else:
            self.rssi = data['val']
            self.t = np.arange(np.shape(self.rssi)[2])*25.832e-3

        self.topandas()
        self.hkb = self.hkb[self.hkb!=0]


    def _distancematrix(self):


        tdev = []
        for k in self.B.dev:
            tdev.append(self.B.dev[k]['uc3d'][0])
        tdev = np.array(tdev)

        # pnb : nframe x ndevices x 3
        pnb = self.B._f[:,tdev,:]
        ln = []
        uin = []
        if ('HK' in self.typ) or ('FULL' in self.typ):
            uin.extend(['HKB:1','HKB:2','HKB:3','HKB:4'])
        if ('TCR' in self.typ) or ('FULL' in self.typ):
            uin.extend(['TCR:32','TCR:24','TCR:27','TCR:28'])
        ln = uin + self.B.dev.keys()
        pin = np.array([self.din[d] for d in uin])
        pin2=np.empty((pnb.shape[0],pin.shape[0],pin.shape[1]))
        pin2[:,:,:]=pin
        p = np.concatenate((pin2,pnb),axis=1)
        self.dist = np.sqrt(np.sum((p[:,:,np.newaxis,:]-p[:,np.newaxis,:,:])**2,axis=3))
        self._lnd = ln


    def accessdm(self,a,b,techno):
        """ Access o the distance matrix

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

        ua = self._lnd.index(ka)
        ub = self._lnd.index(kb)

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
        videofile = self.rootdir+'/POST-TREATED/' +str(self.day) + '-06-2014/Videos/'
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = videofile+_filename
        os.system('vlc '+filename +'&' )



    def snapshot(self,t0=0,offset=15.5,save=False):
        """ single snapshot plot 
        """
        videofile = self.rootdir+'/POST-TREATED/' +str(self.day) + '-06-2014/Videos/'
        ldir = os.listdir(videofile)
        luldir = map(lambda x : self._filename in x,ldir)
        uldir = luldir.index(True)
        _filename = ldir[uldir]
        filename = videofile+_filename
        vc = VideoFileClip(filename)
        F0 = vc.get_frame(t0+offset)
        I0 = img_as_ubyte(F0)
        plt.subplot(111)
        plt.imshow(F0)
        plt.title('t = '+str(t0)+'s')
        if save :
            plt.savefig(self._filename +'_'+str(t0) + '_snap.png',format='png')



    def snapshots(self,t0=0,t1=10,offset=15.5):
        """
        """
        videofile = self.rootdir+'/POST-TREATED/' +str(self.day) + '-06-2014/Videos/'
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

    def show(self):

        fig,ax = self.L.showG('s')
        ax.plot(self.pts[:,0],self.pts[:,1],'o')

    def _show3(self):
        self.L._show3(opacity=0.5)
        v = self.din.items()
        X= np.array([v[i][1] for i in range(len(v))])

        mlab.points3d(X[:,0],X[:,1], X[:,2],scale_factor=0.1)
        [mlab.text3d(v[i][1][0],v[i][1][1],v[i][1][2],v[i][0],scale=0.5)
        for i in range(len(v))]
        for i in range(0,100,5):
            self.B.settopos(t=i,cs=True)
            self.B._show3(dev=True)
        mlab.view(54.989781407516112,
         64.187477298584483,
         20.433867676075128,
         np.array([-0.81123488, -1.65632874, -1.49091462]))



    def topandas(self):
        try:
            self.hkb = pd.DataFrame(index=self.t[0])
        except:
            self.hkb = pd.DataFrame(index=self.t)
        for k in self.idHKB:
            for l in self.idHKB:
                if k!=l:
                    column = self.idHKB[k]+'-'+self.idHKB[l]
                    rssi  = self.rssi[k-1,l-1,:]
                    self.hkb[column] = rssi


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
            timeindex = np.where(self.t[0]-time>0)[0][0]
        except:
            timeindex = np.where(self.t-time>0)[0][0]
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

    def plthkb(self,a,b,**kwargs):
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
                     'reciprocal':True,
                     'data':True
                    }
        
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]
        
        
        t0 =kwargs['t0']
        t1 =kwargs['t1']
        if t1 ==-1:
            t1=self.t[0][-1]

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


        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        if kwargs['ax'] ==[]:
            if kwargs['reciprocal']==True:
                ax = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
            else :
                ax = fig.add_subplot(111)

        if kwargs['data']==True:
            #ax.plot(self.t[0],self.rssi[ia,ib,:])
            #ax.plot(self.t[0],self.rssi[ib,ia,:])
            sab = self.hkb[a+'-'+b]
            sba = self.hkb[b+'-'+a]
            sab[t0:t1].plot(ax=ax)
            sba[t0:t1].plot(ax=ax)
            ax.set_title(a+'-'+b)
        if kwargs['reciprocal']==True:
            # if kwargs['data']==True:
            #     ax2=fig.add_subplot(212)
            r = self.hkb[a+'-'+b][self.hkb[a+'-'+b]!=0]- self.hkb[b+'-'+a][self.hkb[b+'-'+a]!=0]
            r[t0:t1].plot(ax=ax2)
            ax2.set_title('Reciprocity offset')

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
                     'data':True
                    }
        
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]
        
        
        t0 =kwargs['t0']
        t1 =kwargs['t1']
        if t1 ==-1:
            t1=self.t[0][-1]

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
        if kwargs['ax'] ==[]:
            ax = fig.add_subplot(111)

        if kwargs['data']==True:
            #ax.plot(self.t[0],self.rssi[ia,ib,:])
            #ax.plot(self.t[0],self.rssi[ib,ia,:])
            sab = self.tcr[a+'-'+b]
            sba = self.tcr[b+'-'+a]
            sab[t0:t1].plot(ax=ax)
            sba[t0:t1].plot(ax=ax)
            ax.set_title(a+'-'+b)

        return fig,ax

    def pltlk(self):

        display = [t.upper() for t in kwargs['display']]

        if 'HK' in display:
            if ('HK' in self.typ.upper()) or ('FULL' in self.typ.upper()):
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
            else :
                raise AttributeError('HK not available for the given scenario')

        if 'TCR' in display:
            if ('TCR' in self.typ.upper()) or ('FULL' in self.typ.upper()):
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
                raise AttributeError('TCR not available for the given scenario')


#s32 = Hikob(32)
#f,a = s32.imshow(40)
