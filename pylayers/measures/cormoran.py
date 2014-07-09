import os
import pdb
import sys
import pandas as pd
import numpy as np
import scipy.io as io
from pylayers.util.project import *
from moviepy.editor import *
from skimage import img_as_ubyte
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class CorSer(PyLayers):
    """ Hikob data handling from CORMORAN measurement campaign 11/06/2014

    """

    def __init__(self,serie=6,day=11,root='/home/uguen/svn2/measures/CORMORAN/',source='UR1'):

        self.root =root
        if day==11:
            stcr = [1,2,3,4,10,11,12,32,33,34,35,9,17,18,19,20,25,26]
            shkb = [5,6,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
            sbs  = [5,6,7,8,13,14,15,16,21,22,23,24,27,28,29,30,31,32,33,34,35]
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

    def __repr__(self):
        st = ''
        st = st + self._filename + '\n'
        st = st + self._fileBS + '\n'
        return(st)


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
            self.TCR ={'Unused':49,
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
            dirname = self.root+'/POST-TREATED/11-06-2014/TCR'


        if day==12:
            dirname = self.root+'/POST-TREATED/12-06-2014/TCR'
            self.TCR ={ 'COORD':31,
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
        self.iTCR={}
        for k in self.TCR:
            self.iTCR[self.TCR[k]]=k

        self.dTCR={}
        for k in self.TNET.keys():
            self.dTCR[k]=self.iTCR[self.TNET[k]]


        files = os.listdir(dirname)
        if serie != '':
            self._fileTCR = filter(lambda x : 'S'+str(serie) in x ,files)[0]
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            self._fileTCR = filter(lambda x : 'R'+str(run) in x ,filsc)[0]

        filename = dirname + '/'+ self._fileTCR
        dtTCR = pd.read_csv(filename)
        tcr={}
        for k in self.dTCR:
            for l in self.dTCR:
                if k!=l:
                    d = dtTCR[((dtTCR['ida']==k) & (dtTCR['idb']==l))]
                    del d['lqi']
                    del d['ida']
                    del d['idb']
                    d = d[d['time']!=-1]
                    d.index = d['time']
                    del d['time']
                    if len(d)!=0:
                        sr = pd.Series(d['dist']/1000,index=d.index)
                        sr = sr.drop_duplicates()
                        tcr[self.dTCR[k]+'_'+self.dTCR[l]]= sr

        self.tcr = pd.DataFrame(tcr)
        self.tcr = self.tcr.fillna(0)
        ts = 75366400./1e9
        t = np.array(self.tcr.index)*ts
        t = t-t[0]
        self.tcr.index = t

    def loadBS(self,day=11,serie='',scenario='20',run=1):
        """ load BeSpoon data

        """
        self.dBS = {157:'LeftWrist?',74:'RightAnckle?'}
        if day==11:
            dirname = self.root+'/POST-TREATED/11-06-2014/BeSpoon'
        if day==12:
            dirname = self.root+'/POST-TREATED/12-06-2014/BeSpoon'

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
            self.hkb ={'AP1':1,'AP2':2,'AP3':3,'AP4':4,
                       'HeadRight':5,'TorsoTopRight':6,'TorsoTopLeft':7,'BackCenter':8,'ElbowRight':9,'ElbowLeft':10,'HipRight':11,'WristRight':12,'WristLeft':13,'KneeLeft':14,'AnckleRight':16,'AnckleLeft':15}
            if source=='UR1':
                dirname = self.root+'/POST-TREATED/11-06-2014/HIKOB'
            elif source=='CITI':
                dirname = self.root+'/POST-TREATED/11-06-2014/HIKOB/CITI'
        if day==12:
            self.hkb= {'AP1':1,'AP2':2,'AP3':3,'AP4':4,'Jihad:TorsoTopRight':10,'Jihad:TorsoTopLeft':9,'Jihad:BackCenter':11,'JihadShoulderLeft':12,
             'Nicolas:TorsoTopRight':6,'Nicolas:TorsoTopLeft':5,'Nicolas:BackCenter':7,'Nicolas:ShoulderLeft':8,
             'Eric:TorsoTopRight':15,'Eric:TorsoTopLeft':13,'Eric:BackCenter':16,'Eric:ShoulderLeft':14}
            if source=='UR1':
                dirname = self.root+'/POST-TREATED/12-06-2014/HIKOB'
            elif source=='CITI':
                dirname = self.root+'/POST-TREATED/12-06-2014/HIKOB/CITI'

        files = os.listdir(dirname)

        self.ihkb={}
        for k in self.hkb:
            self.ihkb[self.hkb[k]]=k

        if serie != '':
            self._filename = filter(lambda x : 'S'+str(serie) in x ,files)[0]
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            if source=='UR1':
                self._filename = filter(lambda x : 'R'+str(run) in x ,filsc)[0]
            else:
                self._filename = filter(lambda x : 'r'+str(run) in x ,filsc)[0]


        data = io.loadmat(dirname+'/'+self._filename)
        if source=='UR1':
            self.rssi = data['rssi']
            self.t = data['t']
        else:
            self.rssi = data['val']
            self.t = np.arange(np.shape(self.rssi)[2])*25.832e-3

        self.topandas()
        self.df = self.df[self.df!=0]

    def snapshots(self,t0=0,t1=10,offset=15.5):
        """
        """
        videofile = self.root+'POST-TREATED/11-06-2014/Videos/'
        _filename = self._filename.replace('.mat','.mp4')
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



    def topandas(self):
        try:
            self.df = pd.DataFrame(index=self.t[0])
        except:
            self.df = pd.DataFrame(index=self.t)
        for k in self.ihkb:
            for l in self.ihkb:
                if k!=l:
                    column = self.ihkb[k]+'-'+self.ihkb[l]
                    rssi  = self.rssi[k-1,l-1,:]
                    self.df[column] = rssi


    def imshow(self,time):
        """
        Parameters
        ----------

        """
        fig = plt.figure(figsize=(10,10))
        self.D = self.rssi-self.rssi.swapaxes(0,1)
        try:
            timeindex = np.where(self.t[0]-time>0)[0][0]
        except:
            timeindex = np.where(self.t-time>0)[0][0]
        ax1 = fig.add_subplot(121)
        img1 = ax1.imshow(self.rssi[:,:,timeindex],interpolation='nearest',origin='lower')
        labels = map(lambda x : self.ihkb[x],range(1,17))
        plt.xticks(range(16),labels,rotation=80,fontsize=14)
        plt.yticks(range(16),labels,fontsize=14)
        plt.title('t = '+str(time)+ ' s')
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        clb1 = fig.colorbar(img1,cax1)
        clb1.set_label('level dBm')
        ax2 = fig.add_subplot(122)
        img2 = ax2.imshow(self.D[:,:,timeindex],interpolation='nearest',origin='lower')
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

    def pltlk(self,a,b,t0=0,t1=10,fig=[],ax=[],figsize=(8,8),reciprocal=True,data=True):
        """
        """
        ia = self.hkb[a]-1
        ib = self.hkb[b]-1

        if fig==[]:
            fig = plt.figure(figsize=figsize)
        if ax ==[]:
            ax = fig.add_subplot(111)

        if data==True:
            #ax.plot(self.t[0],self.rssi[ia,ib,:])
            #ax.plot(self.t[0],self.rssi[ib,ia,:])
            if reciprocal==True:
                plt.subplot(211)
            sab = self.df[a+'-'+b]
            sba = self.df[b+'-'+a]
            sab[t0:t1].plot()
            sba[t0:t1].plot()
            plt.title(a+'-'+b)
        if reciprocal==True:
            if data==True:
                plt.subplot(212)
            r = self.df[a+'-'+b][self.df[a+'-'+b]!=0]- self.df[b+'-'+a][self.df[b+'-'+a]!=0]
            r[t0:t1].plot()
            plt.title('Reciprocity offset')

        return fig,ax




#s32 = Hikob(32)
#f,a = s32.imshow(40)
