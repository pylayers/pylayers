import os
import pdb
import sys
import pandas as pd
import numpy as np
import scipy.io as io
from pylayers.util.project import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Hikob(PyLayers):
    """ Hikob data handling from CORMORAN measurement campaign 11/06/2014

    """

    def __init__(self,serie=6,root='/home/uguen/svn2/measures/CORMORAN/'):
        self.root =root
        self.load(serie=serie)

    def load(self,day=11,serie='',scenario='20',run=1):
        if day==11:
            self.hkb ={'AP1':1,'AP2':2,'AP3':3,'AP4':4,
                       'HeadRight':5,'TorsoTopRight':6,'TorsoTopLeft':7,'BackCenter':8,'ElbowRight':9,'ElbowLeft':10,'HipRight':11,'WristRight':12,'WristLeft':13,'KneeLeft':14,'AnckleRight':16,'AnckleLeft':15}
            dirname = self.root+'/POST-TREATED/11-06-2014/HIKOB'
        if day==12:
            self.hkb= {'AP1':1,'AP2':2,'AP3':3,'AP4':4,'Jihad:TorsoTopRight':10,'Jihad:TorsoTopLeft':9,'Jihad:BackCenter':11,'JihadShoulderLeft':12,
             'Nicolas:TorsoTopRight':6,'Nicolas:TorsoTopLeft':5,'Nicolas:BackCenter':7,'Nicolas:ShoulderLeft':8,
             'Eric:TorsoTopRight':15,'Eric:TorsoTopLeft':13,'Eric:BackCenter':16,'Eric:ShoulderLeft':14}
            dirname = self.root+'/POST-TREATED/11-06-2014/HIKOB'

        files = os.listdir(dirname)

        self.ihkb={}
        for k in self.hkb:
            self.ihkb[self.hkb[k]]=k

        if serie != '':
            _filename = filter(lambda x : 'S'+str(serie) in x ,files)[0]
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            _filename = filter(lambda x : 'R'+str(run) in x ,filsc)[0]

        data = io.loadmat(dirname+'/'+_filename)
        self.rssi = data['rssi']
        self.t = data['t']

    def topandas(self):
        self.df = pd.DataFrame(index=self.t[0])
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
        D = self.rssi-self.rssi.swapaxes(0,1)
        timeindex = np.where(self.t-time>0)[0][0]
        ax1 = fig.add_subplot(121)
        img1 = ax1.imshow(self.rssi[:,:,timeindex],interpolation='nearest',origin='lower')
        labels = map(lambda x : self.ihkb[x],range(1,17))
        plt.xticks(range(16),labels,rotation=80,fontsize=14)
        plt.yticks(range(16),labels,fontsize=14)
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        clb1 = fig.colorbar(img1,cax1)
        clb1.set_label('level dBm')
        ax2 = fig.add_subplot(122)
        img2 = ax2.imshow(D[:,:,timeindex],interpolation='nearest',origin='lower')
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

    def pltlk(self,a,b,fig=[],ax=[]):
        ia = self.hkb[a]-1
        ib = self.hkb[b]-1
        if fig==[]:
            fig = plt.figure()
        if ax ==[]:
            ax = fig.add_subplot(111)

        ax.plot(self.rssi[ia,ib,:])
        ax.plot(self.rssi[ib,ia,:])
        ax.plot(self.rssi[ib,ia,:]-self.rssi[ia,ib])

        return fig,ax




s32 = Hikob(32)
f,a = s32.imshow(40)
