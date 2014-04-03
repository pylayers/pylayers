#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""
   Run simulation with full human trajectory

.. currentmodule:: pylayers.simul.simultraj

Simul class
===========

.. autosummary::
    :toctree: generated/

     Simul.__init__
     Simul.load
     Simul.links_generation
     Simul.run

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
from pylayers.mobility.ban.body import *
from pylayers.antprop.statModel import *

class Simul(object):
    """
    Link oriented simulation

    A simulation requires :

        A Layout
        A Person
        A Trajectory

    """
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
        """  load

        Parameters
        ----------

        _filesimul :

        """
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

        # get the layout
        self.L = Layout(di['layout']['layout'])
        # get the trajectory
        traj  = tr.importsn(di['trajectory']['traj'])
        # resample trajectory
        for i in range(len(traj)):
            traj[i] = traj[i].resample(1)
        self.traj= traj
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
        """ links_generation

        Parameters
        ----------

        bB2B : boolean
        bOB  : boolean
        bB2I : boolean

        """

        lAP = self.ap
        lperson = self.dpersons.values()
        llink = []
        #pdb.set_trace()
        for person in lperson:
            if bOB:
                for dev1 in person.dev:
                    #if person.dev[dev1]['typ']=='dynamic':
                    for dev2 in person.dev:
                        if dev2<> dev1:
                            if ((person.name,dev2),(person.name,dev1)) not in llink :
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



    def run(self,llink = [], t =[]):
        """
        Parameters
        ----------

        llink : list

        """

        if llink ==[]:
            llink = self.links
        if t  == []:
            time =  self.traj[0].time()
        else:
            time =  self.traj[0].time()[t]


        n_time =len(time)
        n_links = len(llink)
        Kmax = 100
        resultEnv  = np.zeros(shape=(n_time,n_links,Kmax,2))
        resultOb  = np.zeros(shape=(n_time,n_links,Kmax,2))


        for kt in range(0,n_time):

            #~ if time[kt]%int(time[kt]) == 0:
                #~ print 't = ', time[kt]
            print 't = ', time[kt]
            for kp, person in enumerate(self.dpersons.values()):
                person.settopos(self.traj[kp],t=time[kt],cs=True)


            plt.figure()
            self.dpersons['Alex'].show(color='b',plane='xz',topos=True)
            plt.show()


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

                #interA = self.dpersons[A[0]].intersectBody3(pA,pB, topos = True)
                ##interB = self.dpersons[B[0]].intersectBody2(pA,pB, topos = True)


                #condition ='nlos'
                #if interA==1:
                    #condition = 'los'
                #empA = A[1]
                #empB = B[1]
                #emp  = empA
                #if empA == 'accelerometer':
                    #emp = empB
                #if emp == 'left_watch':
                    #emp = 'right_watch'
                    ##~ if condition == 'los':
                        ##~ condition = 'nlos'
                    ##~ else:
                        ##~ condition = 'los'
                #if emp == 'front_chest':
                    #condition = 'los'


                #alphakOb, taukOb = getchannel(emplacement = emp ,condition = condition, intersection = interA)
                #alphak =  np.zeros(shape=(Kmax))
                #tauk   =  np.zeros(shape=(Kmax))
                #ntraj = len(taukOb)
                #if ntraj < Kmax:
                    #alphak[0:ntraj] =  alphakOb[0:ntraj]
                    #tauk[0:ntraj] =  taukOb[0:ntraj]
                #else:
                    #alphak =  alphakOb[0:Kmax]
                    #tauk =  taukOb[0:Kmax]
                    #print ' warning ntraj > Kmax'

                #tab = np.vstack((alphak,tauk)).T
                #resultOb[kt,kl,:,:]= tab

                cycA =  self.L.pt2cy(pt = pA)
                cycB =  self.L.pt2cy(pt = pB)
                sig = signature.Signatures(self.L,cycA,cycB)
                #sig.run4(cutoff =1,algo='old')
                sig.run5(cutoff = 3)
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
                print 'link = ', link , '  pr  = ',10*np.log10(sum((abs(alphak)**2)))

        return resultEnv#, resultOb



    def runEnv(self,llink = [], t =[], show = False ):
        """
        Parameters
        ----------

        llink : list

        """

        if llink ==[]:
            llink = self.links
        if t  == []:
            time =  self.traj[0].time()
        else:
            time =  self.traj[0].time()[t]


        n_time =len(time)
        n_links = len(llink)
        Kmax = 100
        resultEnv  = np.zeros(shape=(n_time,n_links,Kmax,2))
        resultOb  = np.zeros(shape=(n_time,n_links,Kmax,2))
        taumin  = 0
        taumax  = 300
        taustep = 0.1
        x = np.arange(taumin,taumax,taustep)
        y = np.zeros(len(x))
        cira=bs.TUsignal(x,y)
        self.dist = np.zeros(shape=(n_time,n_links))
        self.vis  = np.zeros(shape=(n_time,n_links))   
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
                
                self.dist[kt,kl] = np.sqrt(sum((pA-pB)**2))
                interA = self.dpersons[A[0]].intersectBody3(pA,pB, topos = True)
               
                #interB = self.dpersons[B[0]].intersectBody3(pA,pB, topos = True)
                self.vis[kt,kl] = interA
                
                cylA = self.dpersons[A[0]].dev[A[1]]['cyl']
                cylB = self.dpersons[B[0]].dev[B[1]]['cyl']


                cycA =  self.L.pt2cy(pt = pA)
                cycB =  self.L.pt2cy(pt = pB)
                sig = signature.Signatures(self.L,cycA,cycB)
                sig.run5(cutoff=1,algo='old')
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
                if show:
                    print 'link  = ', link , '  pr env  = ',10*np.log10(sum((abs(alphak)**2)))
                #cira.aggcir(alphak,tauk)
                #self.cira  = cira
                self.chan = resultEnv
                
        #return resultEnv, cira

    def runOb(self,llink = [], t =[], show = False ):
        """
        Parameters
        ----------

        llink : list

        """

        if llink ==[]:
            llink = self.links
        if t  == []:
            time =  self.traj[0].time()
        else:
            time =  self.traj[0].time()[t]


        n_time =len(time)
        n_links = len(llink)
        Kmax =30
        resultEnv  = np.zeros(shape=(n_time,n_links,Kmax,2))
        resultOb  = np.zeros(shape=(n_time,n_links,Kmax,2))

        for kt in range(0,n_time):

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


                
                devIdA = A[1]
                devIdB = B[1]
                empA =self.dpersons[A[0]].dev[devIdA]['cyl']
                empB =self.dpersons[B[0]].dev[devIdB]['cyl']

                emp  = empA
                if empA == 'trunkb':
                    emp = empB
                if emp == 'forearml':
                    emp = 'forearmr'

                
                alphakOb, taukOb = getchannel(emplacement = emp, intersection = interA)


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
                if show:
                    print 'link  = ', link , '  pr ob  = ',10*np.log10(sum((abs(alphak)**2)))

        return resultOb
        
    def runObF(self,llink = [], t =[]):
        """
        Parameters
        ----------

        llink : list

        """

        if llink ==[]:
            llink = self.links
        if t  == []:
            time =  self.traj[0].time()
        else:
            time =  self.traj[0].time()[t]


        n_time =len(time)
        n_links = len(llink)
        Kmax = 30
        resultOb  = np.zeros(shape=(n_time,n_links,Kmax,2))
        slink  = np.zeros(shape=(n_time,n_links))
        tlink  = []
        for kl in range(0,n_links):
            link = llink[kl]
            A = link[0]
            B = link[1]
            typA = self.dpersons[A[0]].dev[A[1]]['typ']
            typB = self.dpersons[B[0]].dev[B[1]]['typ']
            if typA == 'dynamic':
                tlink.append(typA)
            else:
                if typB == 'dynamic':
                    tlink.append(typB)
                else:
                    tlink.append(typA)
           
        salpha = []
        stau   = []
        for kt in range(0,n_time):

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


                
                devIdA = A[1]
                devIdB = B[1]
                empA =self.dpersons[A[0]].dev[devIdA]['cyl']
                empB =self.dpersons[B[0]].dev[devIdB]['cyl']

                emp  = empA
                if empA == 'trunkb':
                    emp = empB
                if emp == 'forearml':
                    emp = 'forearmr'

                if kt ==0:
                    
                    a, d = getchannel(emplacement = emp, intersection = interA)
                    alpha_e0  = np.zeros(shape=(Kmax))
                    alpha_e0[0:len(a)]= a
                    tau_e0 = np.zeros(shape=(Kmax))
                    tau_e0[0:len(d)]= d
                    alpha =  np.zeros(shape=(Kmax))
                    tau   =  np.zeros(shape=(Kmax))
                    ntraj = len(a)
                    if ntraj < Kmax:
                        alpha[0:ntraj] =  alpha_e0[0:ntraj]
                        tau[0:ntraj] =  tau_e0[0:ntraj]
                    else:
                        alpha =  alpha_e0[0:Kmax]
                        tau =  tau_e0[0:Kmax]
                        print ' warning ntraj > Kmax'
                    
                    tab = np.vstack((alpha,tau)).T
                    resultOb[kt,kl,:,:]= tab
                    salpha.append(alpha_e0)
                    stau.append(tau_e0)
                else:
                    
                    if tlink[kl]== 'static':
                        bf = np.array([ 2*0.04343642])
                        af = np.array([ 1.,         -0.91312716])
                        gf =  0.2
                    else:
                        bf = np.array([ 2*0.40185134 ])
                        af = np.array([ 1. ,        -0.19629732])
                        gf = 0.63
                    
                    a, t = getchannel(emplacement = emp, intersection = interA)
                    alpha_e1  = np.zeros(shape=(Kmax))
                    alpha_e1[0:len(a)]= a
                    tau_e1 = np.zeros(shape=(Kmax))
                    tau_e1[0:len(t)]= t
                    
                    
                    alpha =  np.zeros(shape=(Kmax))
                    tau  =  np.zeros(shape=(Kmax))
                    ntraj = len(a)
                    if ntraj < Kmax:
                        alpha[0:ntraj] =  abs(bf[0]*(alpha_e1[0:ntraj]-(1-gf)*salpha[kl][0:ntraj]/kt)/gf-af[1]*resultOb[kt-1,kl,0:ntraj,0])
                        tau[0:ntraj] =  (bf[0]*(tau_e1[0:ntraj]-(1-gf)*stau[kl][0:ntraj]/kt)/gf-af[1]*resultOb[kt-1,kl,0:ntraj,1])
                        
                    else:
                        alpha =  (bf[0]*alpha_e1[0:Kmax]-af[1]*resultOb[kt-1,kl,0:ntraj,0])
                        tau =  (bf[0]*tau_e1[0:Kmax]-af[1]*resultOb[kt-1,kl,0:ntraj,1])
                        print ' warning ntraj > Kmax'
                    
                    tab = np.vstack((alpha,tau)).T
                    resultOb[kt,kl,:,:]= tab
                    salpha[kl] = salpha[kl]+alpha_e1
                    stau[kl] = stau[kl]+tau_e1
        
        
        return resultOb
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
