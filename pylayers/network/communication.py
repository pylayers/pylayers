# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of Network.

#Foobar is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Foobar is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Nicolas AMIOT        : nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN        : bernard.uguen@univ-rennes1.fr
#####################################################################

import numpy as np
#import scipy as sp 
#import networkx as nx
#import itertools
#import pickle as pk
#import pkgutil

#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#import Tkinter as tk

#import ConfigParser
#import copy
#import pdb
##from PyLayers.Network.Node import Node
#import pylayers.util.pyutil as pyu
#from pylayers.network.emsolver import EMSolver
#from pylayers.network.show import ShowNet,ShowTable
#from pylayers.util.pymysqldb import Database 
#import pylayers.util.pyutil as pyu

#import time
from pylayers.mobility.agent import Agent
from pylayers.network.emsolver import EMSolver
from pylayers.gis.layout import Layout
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import ConfigParser
import SimPy.Simulation
from SimPy.Simulation import Process,hold,SimEvent,Simulation,waitevent
from random import uniform,gauss
from pylayers.network.network import  Node,Network
import pdb




class dcond(dict):

    def __init__(self,ID=0):
        Cf = ConfigParser.ConfigParser()
        Cf.read(pyu.getlong('agent.ini','ini'))
        self.ag_opt = dict(Cf.items(ID))
        self.parse()

    def parse(self):
        string = 'condition'
        C=self.ag_opt['condition'].split('\n')
        for cit,c in enumerate(C):
            items=c.split(';')
            for k,it in enumerate(items):
                if k == 0:
                    self[str(cit)]={}
                    self[str(cit)][it.split(':')[0][1:]]=it.split(':')[1][1:-1]
#                elif it == 'message':
                elif it.split(':')[0]=='message':
                    self[str(cit)][it.split(':')[0]]=eval(it.split(':')[1])
                elif it=='and"' or it=='or"':
                    self[str(cit)]['and_or']=it[:-1]
                else:
                    try:
                        self[str(cit)][it.split(':')[0]]=it.split(':')[1][1:-1] 
                    except: 
                        pass




class TX(Process):
     def __init__(self,**args):
        defaults={'sim':None,
                  'net': Network(),
                  'ID': 0,
                  'dcond':{},
                  'levt': [],
                  'lcst': []
                  }
##       initialize attributes
        for key, value in defaults.items():
            if args.has_key(key):

                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key]=value  
        self.dcond = args['dcond']
        self.args=args
        self.net=args['net']
        self.PN=self.net.node[self.ID]['PN']

        Process.__init__(self,name='Tx'+str(self.ID),sim=self.sim)


     def run(self):
        while True:
            self.levt=[]
            for c in self.cond:
                if self.c_interpret(c):
                    self.levt.append(self.dcste[c])
            print 'Tx ', self.ID,' @',self.sim.now()
            yield hold, self, self.refresh


#    def c_interpret(self,c):
#        """
#        test si la constrainte est vrai a l'instant t
#         """
#        return True


#    def evt_create(self):




#class RX(Process):
#     def __init__(self,**args):
#        defaults={'sim':None,
#                  'ID':'A1'
#                  }
###       initialize attributes
#        for key, value in defaults.items():
#            if args.has_key(key):

#                setattr(self, key, args[key])
#            else:
#                setattr(self, key, value)
#                args[key]=value  
#        self.args=args
#        self.net=args['net']
#        self.sim=args['sim']
#        self.ID==args['ID']


#        Process.__init__(self,name='Rx-'+str(self.ID),sim=self.sim)

#        Cf = ConfigParser.ConfigParser()
#        Cf.read(pyu.getlong('agent.ini','ini'))
#        ag_opt = dict(Cf.items(self.ID))
#        self.cond = eval(ag_opt['condition'])
#        self.inode = eval(ag_opt['inode'])
#        self.irat = eval(ag_opt['irat'])
#        self.mess = eval(ag_opt['message'])


#     def run(self):
#        while True:
#            yield waitevent, self, self.event
#            print 'Tx ', self.ID,' @',self.sim.now()

if (__name__ == "__main__"):

    sim =Simulation()
    sim.initialize()


    L=Layout('TA-Office.str')
    L.build('str') # build 's'tructure, 't'opological and 'r'oom graphs
    N=Network()

    Ag=[]
    Cf = ConfigParser.ConfigParser()
    Cf.read(pyu.getlong('agent.ini','ini'))
    agents=eval(dict(Cf.items('used_agent'))['list'])
    for i, ag in enumerate(agents):
        ag_opt = dict(Cf.items(ag))
        print ag_opt['id']
        Ag.append(Agent(
                          ID=ag_opt['id'],
                          name=ag_opt['name'],
                          type=ag_opt['type'],
                          roomId=int(ag_opt['roomid']),
                          pos=np.array(eval(ag_opt['pos'])),
                          Layout=L,
                          net=N,
                          RAT=eval(ag_opt['rat']),
                          dcond=dcond(ag),
                          epwr=dict([(eval((ag_opt['rat']))[ep],eval((ag_opt['epwr']))[ep]) for ep in range(len(eval((ag_opt['rat']))))]),
                          sim=sim)
                  )


    N.create()
    N.EMS=EMSolver(L)

    for ldp in ['Pr','TOA']:
        N.compute_LDPs(N.nodes(),RAT='rat1',LDP=ldp)


    N.update_PN()



    tx1=TX(net=N,ID=Ag[0].ID,dcond=Ag[0].dcond)
    tx2=TX(net=N,ID=Ag[1].ID,dcond=Ag[1].dcond)
##    N.node[0]['PN'].node[0]['pe']=np.array((4,4))
##    N.node[0]['PN'].node[1]['pe']=np.array((8,8))
##    N.node[0]['PN'].node[2]['pe']=np.array((30,8))

##    nbtx=2
##    tx=[]
##    evt=[]
##    b=[]
##    for t in range(nbtx):
##        evt.append(SimEvent('pos'+str(t),sim=sim))
###        tx.append(TX(ID=t,sim=sim,evt=evt[t]))
##        b.append(Brain(ID=t,sim=sim,evt=evt[t]))
###        sim.activate (tx[t], tx[t].run(),0.0)
##        sim.activate (b[t], b[t].run(),0.0)
##    
##    T=TX(ID=t,sim=sim,evt=evt)
##    sim.activate (T, T.run(),0.0)
##    sim.simulate(until=10)

