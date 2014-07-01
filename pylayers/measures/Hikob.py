import os
import pdb
import sys
import pandas as pd
import numpy as np
import scipy.io as io
from pylayers.util.project import *
import matplotlib.pyplot as plt


class Hikob(PyLayers):
    """
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

        self.iHKB={}
        for k in self.hkb:
            self.iHKB[self.hkb[k]]=k

        if serie != '':
            _filename = filter(lambda x : 'S'+str(serie) in x ,files)[0]
        else:
            filesc = filter(lambda x : 'Sc'+scenario in x ,files)
            _filename = filter(lambda x : 'R'+str(run) in x ,filsc)[0]

        data = io.loadmat(dirname+'/'+_filename)
        self.rssi = data['rssi']
        self.t = data['t']
        #for k in range(1,17):
        #    for l in range(1,17):
        #        self.dHKB[(k,l)]=iHKB[k]+' - '+iHKB[l]
        #        cpt = cpt + 1






