#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""

    This module run simulation in complete trajectory


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
from pylayers.antprop.raysc import GrRay3D, GrRayTud
# Handle VectChannel and ScalChannel
from pylayers.antprop import antenna,channelc,signature
#from   Channel import *
# Handle directory hierarchy
from pylayers.util.project import *
# Handle UWB measurements
from pylayers.measures import mesuwb as muwb
import pdb
import pylayers.mobility.trajectory as tr
from pylayers.mobility.body.body1 import *
from pylayers.antprop.delaydispersion import *

class Simul(object):
	
    def __init__(self, _filesimul='simultraj.ini'):

        self.filesimul = _filesimul
        self.config = ConfigParser.ConfigParser()
        self.config.add_section("layout")
        self.config.add_section("person")
        self.config.add_section("trajectory")

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
        
        self.progress = -1  # simulation not loaded
        

        self.filetra = []
        self.filetud = []
        self.filetang = []
        self.filerang = []
        self.filetauk = []
        self.filefield = []

        self.claunching = []
        self.ctracing = []
        self.ctratotud = []
    
        self.cfield = []
        fmin = 3
        fmax = 6
        fstep = 0.05
        nf = (fmax-fmin)/fstep
        self.fGHz = np.linspace(fmin, fmax, nf, endpoint=True)
        self.wav = wvf.Waveform()
        
        self.load(_filesimul)
        
        
    def load (self,_filesimul):
        
        self.filesimul = _filesimul
        filesimul = pyu.getlong(self.filesimul, "ini")
     
       
        config = ConfigParser.ConfigParser()
        config.read(filesimul)
        sections = config.sections()
        di = {}
        for section in sections:
            
            di[section] = {}
            options = config.options(section)
            for option in options:
               di[section][option] = config.get(section,option)
        
        self.L = Layout(di['layout']['layout'])
        self.traj= tr.importsn(di['trajectory']['traj'])
        self.ap = di['acces_point']['ap']
        
        persons_files  = di['person'].values()
        self.dpersons = {}
        
        for person_file in persons_files:
            person_name  = person_file.replace('.ini','')
            self.dpersons.keys().append(person_name)
            person = Body(person_file)
            self.dpersons[person_name] =person
            
            person_name = Body(person_file)
            
    def links_generation(self, bB2B = False, bOB  = True, bB2I = False):

        lAP = self.ap
        lperson = self.dpersons.values()
        llink = []
        for person in lperson:
            
            if bOB:
                for dev1 in person.dev:
                    if person.dev[dev1]['typ']=='mobile':
                        for dev2 in person.dev:
                            if dev2<> dev1:
                                llink.append(((person.name,dev1),(person.name,dev2)))
                
            if bB2I: 
                for ap in lAP:
                    # may be criteria on distance between person and AP
                    for dev in person.dev:
                        #if person.dev[dev1]['infra']:
                        # in the future port would be an antenna port for MIMO AP
                        llink.append(((ap,'port'),(person.name,dev)))

            if bB2B:
                for dev1 in person.dev:
                    #if person.dev[dev1]['tobody']:
                    for alter in lperson:
                        if alter.name<>person.name:
                            for dev2 in alter.dev:
                                #if alter.dev[dev2]['tobody']:
                                llink.append(((person.name,dev1),(alter.name,dev2)))
        
        self.links  = llink
        
            
    
    def run(self,llink = []):
        
        if llink ==[]:
            llink = self.links
            
        time =  self.traj[0].time()[0:30]
        n_time =len(time)
        n_links = len(llink)
        Kmax = 100
        resultEnv  = np.zeros(shape=(n_time,n_links,Kmax,2))
        resultOb  = np.zeros(shape=(n_time,n_links,Kmax,2))
     
        
        for kt in range(0,n_time):
            
            if time[kt]%int(time[kt]) == 0:
                print 't = ', time[kt]           
            
            for kp, person in enumerate(self.dpersons.values()):
                person.settopos(self.traj[kp],t=time[kt],cs=True)
                
            
            for kl in range(0,n_links):
                link = llink[kl]
                A = link[0]
                B = link[1]
                pA = self.dpersons[A[0]].dcs[A[1]][:,0]
                TA = self.dpersons[A[0]].acs[A[1]]
                pB = self.dpersons[B[0]].dcs[B[1]][:,0]
                TB = self.dpersons[B[0]].acs[B[1]]
                cylA = self.dpersons[A[0]].dev[A[1]]['cyl']
                cylB = self.dpersons[B[0]].dev[B[1]]['cyl']
                
                interA = self.dpersons[A[0]].intersectBody3(pA,pB, topos = True)
                #interB = self.dpersons[B[0]].intersectBody2(pA,pB, topos = True)
                    
                   
                condition ='nlos'
                if interA==1:
                    condition = 'los'
                empA = A[1]
                empB = B[1]
                emp  = empA
                if empA == 'accelerometer':
                    emp = empB
                if emp == 'left_watch':
                    emp = 'right_watch'
                    #~ if condition == 'los':
                        #~ condition = 'nlos'
                    #~ else:
                        #~ condition = 'los'
                if emp == 'front_chest':
                    condition = 'los'
                
               
                alphakOb, taukOb = getchannel(emplacement = emp ,condition = condition, intersection = interA)
                alphak =  np.zeros(shape=(Kmax))
                tauk   =  np.zeros(shape=(Kmax))
                ntraj = len(taukOb)
                if ntraj < Kmax:               
                    
                    alphak[0:ntraj] =  alphakOb[0:ntraj]
                    tauk[0:ntraj] =  taukOb[0:ntraj]
                else:                    
                    alphak =  alphakOb[0:Kmax]
                    tauk =  taukOb[0:Kmax]                            
                    print ' warning ntraj > Kmax'
               
                tab = np.vstack((alphak,tauk)).T                   
                resultOb[kt,kl,:,:]= tab  
                  
                cycA =  self.L.pt2cy(pt = pA)                
                cycB =  self.L.pt2cy(pt = pB)
                sig = signature.Signatures(self.L,cycA,cycB)
                sig.run4(cutoff =1,algo='old')
                tx_2D = pA[0:2]
                rx_2D = pB[0:2]
                r2d = sig.rays(tx_2D,rx_2D)
                r2d.pTx = pA
                r2d.pRx = pB
                r3d = r2d.to3D(self.L)            
                r3d.locbas(self.L)
                r3d.fillinter(self.L)
                Cn = r3d.eval(fGHz=self.fGHz)
                Cn.locbas(Tt= TA,Tr =TB)
                AntA = antenna.Antenna(self.dpersons[A[0]].dev[A[1]]['file']) 
                AntB = antenna.Antenna(self.dpersons[B[0]].dev[B[1]]['file']) 
                H = Cn.prop2tran(a =AntA,b=AntB)
                H.applyFriis()
                ntraj = H.y.shape[0]
                alphak =  np.zeros(shape=(Kmax))
                tauk   =  np.zeros(shape=(Kmax))
                if ntraj < Kmax:
                    #pdb.set_trace()
                    
                    alphak[0:ntraj] =  np.real(np.sqrt(np.sum(H.y*np.conj(H.y),axis =1))/len(self.fGHz))
                    tauk[0:ntraj] =  H.tau0
                else:                    
                    alphak =  np.real(np.sqrt(np.sum(H.y*np.conj(H.y),axis =1))/len(self.fGHz))[0:Kmax]
                    tauk =  (H.tau0)[0:Kmax]                            
                    print ' warning ntraj > Kmax'
                 
                tab = np.vstack((alphak,tauk)).T
                   
                resultEnv[kt,kl,:,:]= tab
                
        
        return resultEnv, resultOb
        
        
               
            
       


        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


		
	
