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
#from pylayers.mobility.agent import Agent
from pylayers.network.emsolver import EMSolver
from pylayers.gis.layout import Layout
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import ConfigParser
import SimPy.Simulation
from SimPy.Simulation import Process,hold,SimEvent,Simulation,waitevent
from random import uniform,gauss
from pylayers.network.network import  Node,Network
import networkx as nx
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
                  'gcom': Gcom(),
                  'ID': 0,
                  'dcond':{},
                  'devt': {},
                  'lcst': []
                  }
##       initialize attributes
        for key, value in defaults.items():
            if args.has_key(key):

                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key]=value  
        self.args=args
        self.PN=self.net.node[self.ID]['PN']
        self.evt_create()
        self.c_init()
        Process.__init__(self,name='Tx'+str(self.ID),sim=self.sim)




    def run(self):
        while True:
            self.levt=[]
            for d in self.dcond.keys():
                if self.c_interpret(self.dcond[d]):
                    self.levt.append(self.dcste[c])
            print 'Tx ', self.ID,' @',self.sim.now()
            yield hold, self, self.refresh


    def evt_create(self):
#        for e in self.gcom.edges(self.ID,keys=True):
        for e in self.PN.edges(self.ID,keys=True):
            self.devt[e]=self.gcom.devt[e]


    def c_init(self):
        for dk in self.dcond.keys():
            d=self.dcond[dk]
            ### Rat
            if d['rat']=='all':
                lr=self.PN.SubNet.keys()
            else:
                try:
                    lr=d['rat'] # d['rat'] mustr contain a list of
                               # rat to be processed
                except: 
                    raise NameError('rat constraints must be a list of \
                    rat available in the personnal network of node' \
                    +str(self.ID) +'Please modify our agent.ini')


            ### Node
            ie =[]
            if d['node']=='all':
                if not isinstance(lr,list):
                    lr=[lr]
                for r in lr:
                    dg = nx.DiGraph(self.PN.SubNet[r])
                    ie.append(dg.edges_iter(self.ID))
            else:
                try:
                    for r in lr:
                        dg = nx.DiGraph(self.PN.SubNet[r])
                        iet.append(dg.edges(self.ID))
                        for i in iet:
                            for j in d['node']:
                                if j in i :
                                    ie.append(i)
                      # d['node'] must contain a list of node to be processed
                except:
                    raise NameError('node constraints must be a list of \
                    node available in the personnal network of node' \
                    +str(self.ID) +'Please modify our agent.ini')


            lmess = []
            di={}
            [di.update({it:[]}) for it in d['message'] ]
            for r in lr:
                if d['node']=='all':
                    lmess.append([{'message':di}] * len(self.PN.SubNet[r].edges()))
                else:
                    lmess.append([{'message':di}] * len(self.PN.SubNet[r].edges(d['node'])))

            self.gcom.fill_edge(ie,lr,lmess)




    def c_interpret(self,d):

        # first great criteria
        if 'distance' in d.keys():
            pass
        elif 'topology' in d.keys():
            pass
        else: 
            NameError('Type of mp constraint not yet taking into account')






class RX(Process):
    def __init__(self,**args):
        defaults={'sim':None,
                  'ID':'1',
                  'net':Network(),
                  'gcom': Gcom(),
                  'devt': {},
                  'refreshRSS':0.3,
                  'refreshTOA':0.3,
                  'mp':False
                  }
##       initialize attributes
        for key, value in defaults.items():
            if args.has_key(key):

                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key]=value  
        self.args=args
        self.PN=self.net.node[self.ID]['PN']
        self.create_evt()

        Process.__init__(self,name='Rx-'+str(self.ID),sim=self.sim)

        Cf = ConfigParser.ConfigParser()
        Cf.read(pyu.getlong('agent.ini','ini'))
        for s in Cf.sections():
            try:
                d= dict(Cf.items(s))
                if d['id']==self.ID:
                    self.refreshRSS=eval(d['refreshrss'])
                    self.refreshTOA=eval(d['refreshtoa'])
                break
            except:
                pass



    def swap_lt(self,lt):
        """
        swap edge tuple node number 
        """
        lto=[]
        for t in lt :
            if t[0] == self.ID:
                lto.append((t[1],t[0],t[2]))
            else :
                lto.append((t[0],t[1],t[2]))
        return lto



    def create_evt(self):
        rPNe = self.swap_lt(self.PN.edges(keys=True))
        for e in rPNe:
            self.devt[e]=self.gcom.devt[e]
        if self.mp :
            for e in rPNe:
                self.devt[e]=self.gcom.devt[e+str(mp)]


    def refresh_RSS(self):
        while 1:
            for rat in self.PN.SubNet.keys():
                [self.PN.edge[self.ID][n][rat].update(
                {'Pr':self.net.edge[self.ID][n][rat]['Pr'],'tPr':self.sim.now()})
                for n in self.PN.SubNet[rat].edge[self.ID].keys()]
#            print 'refresh RSS node', self.ID, ' @',self.sim.now()
            yield hold, self, self.refreshRSS


    def refresh_TOA(self):
        while 1:
            for rat in self.PN.SubNet.keys():
                [self.PN.edge[self.ID][n][rat].update(
                {'TOA':self.net.edge[self.ID][n][rat]['TOA'],'tTOA':self.sim.now()})
                for n in self.PN.SubNet[rat].edge[self.ID].keys()]
#            print 'refresh TOA node', self.ID, ' @',self.sim.now()
            yield hold, self, self.refreshTOA

#                [[self.SubNet[RAT].node[e]['PN'].edge[e][f][RAT].update(
#         {'TOA':self.SubNet[RAT].edge[e][f][RAT]['TOA'],'tTOA':self.sim.now()}) 
#         for f in self.SubNet[RAT].edge[e].keys()]
#         for e in self.SubNet[RAT].nodes()]


class Gcom(nx.MultiDiGraph):


    def __init__(self,net=Network(),sim=Simulation()):
        nx.MultiDiGraph.__init__(self)
        self.net=net
        self.sim=sim


    def create(self):
        self.create_graph()
        self.create_evt()

    def create_graph(self):
        for rat in self.net.SubNet:
            for n in self.net.SubNet[rat].nodes():
                G=nx.DiGraph(self.net.SubNet[rat])
                le = G.edges(n)
                ld = [{'message':[],'t':-1}] * len(le)
                try:
                    if le[0][0] == n :
                        self.add_edges_from(self.net.Gen_tuple(G.edges_iter(n),rat,ld))
                    else :
                        self.add_edges_from(self.net.Gen_tuple(nx.DiGraph(le).reverse().edges_iter(),rat,ld))
                except:
                    print 'WARNING : no edge on rat',rat


    def create_evt(self):
        self.devt={}
        for e in self.edges(keys=True):
            self.devt[e]=(SimEvent(e,sim=self.sim))


    def fill_edge(self,le,rat,mess):

        for i,r in enumerate(rat):
            self.add_edges_from(self.net.Gen_tuple(le[i],r,mess[i]))
#        Z=self.dcond['1']['message']*len(le[1])
#        self.gmp.add_edges_from(self.net.Gen_tuple(le,'rat1',Z))


#if (__name__ == "__main__"):

#    sim =Simulation()
#    sim.initialize()


#    L=Layout('TA-Office.str')
#    L.build('str') # build 's'tructure, 't'opological and 'r'oom graphs
#    N=Network()

#    Ag=[]
#    Cf = ConfigParser.ConfigParser()
#    Cf.read(pyu.getlong('agent.ini','ini'))
#    agents=eval(dict(Cf.items('used_agent'))['list'])
#    for i, ag in enumerate(agents):
#        ag_opt = dict(Cf.items(ag))
#        print ag_opt['id']
#        Ag.append(Agent(
#                          ID=ag_opt['id'],
#                          name=ag_opt['name'],
#                          type=ag_opt['type'],
#                          roomId=int(ag_opt['roomid']),
#                          pos=np.array(eval(ag_opt['pos'])),
#                          Layout=L,
#                          net=N,
#                          RAT=eval(ag_opt['rat']),
#                          dcond=dcond(ag),
#                          epwr=dict([(eval((ag_opt['rat']))[ep],eval((ag_opt['epwr']))[ep]) for ep in range(len(eval((ag_opt['rat']))))]),
#                          sim=sim)
#                  )


#    N.create()
#    N.EMS=EMSolver(L)

#    for ldp in ['Pr','TOA']:
#        N.compute_LDPs(N.nodes(),RAT='rat1',LDP=ldp)


#    N.update_PN()


#    gcom=Gcom(net=N,sim=sim)
#    gcom.create_graph()
#    gcom.create_evt()
#    tx=[]
#    rx=[]
#    for a in Ag:
##        tx.append(TX(net=N,ID=a.ID,dcond=a.dcond,gcom=gcom,sim=sim))
#        rx.append(RX(net=N,ID=a.ID,dcond=a.dcond,gcom=gcom,sim=sim))




###    N.node[0]['PN'].node[0]['pe']=np.array((4,4))
###    N.node[0]['PN'].node[1]['pe']=np.array((8,8))
###    N.node[0]['PN'].node[2]['pe']=np.array((30,8))

###    nbtx=2
###    tx=[]
###    evt=[]
###    b=[]
###    for t in range(nbtx):
###        evt.append(SimEvent('pos'+str(t),sim=sim))
####        tx.append(TX(ID=t,sim=sim,evt=evt[t]))
###        b.append(Brain(ID=t,sim=sim,evt=evt[t]))
####        sim.activate (tx[t], tx[t].run(),0.0)
###        sim.activate (b[t], b[t].run(),0.0)
###    
###    T=TX(ID=t,sim=sim,evt=evt)
###    sim.activate (T, T.run(),0.0)
###    sim.simulate(until=10)

