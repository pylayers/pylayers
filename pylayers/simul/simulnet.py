#   -*- coding:Utf-8 -*-
#####################################################################
#This file is part of LocaSimPy.

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

import pkgutil
import warnings
warnings.filterwarnings('ignore')

#import SimPy.Simulation
from SimPy.SimulationRT import SimulationRT, Process, hold
#import simpy # simpy 3
import numpy as np
import scipy as sp
import networkx as nx
import random
from random import seed

import time
import matplotlib.pyplot as plt
import ConfigParser

import pylayers.util.pyutil as pyu

from pylayers.network.network import Network, Node, PNetwork
from pylayers.network.communication import Gcom
from pylayers.network.show import ShowNet, ShowTable
from pylayers.mobility.agent import Agent
from pylayers.network.emsolver import EMSolver
from pylayers.gis.layout import Layout
from pylayers.antprop.slab import Slab
from pylayers.util.utilnet import str2bool
from pylayers.mobility.transit.World import world
#from pylayers.util.pymysqldb import Database as DB
from pylayers.util.project import *
from pylayers.util.save import *

import pdb
import os

class Simul(SimulationRT): # Sympy 2
#class Simul(sympy.RealtimeEnvironment):
    """

    Attributes
    ----------
    config  : config parser instance
    sim_opt : dictionary of configuration option for simulation
    ag_opt : dictionary of configuration option for agent
    lay_opt : dictionary of  configuration option for layout
    meca_opt : dictionary of  configuration option for mecanic
    net_opt : dictionary of  configuration option for network
    loc_opt : dictionary of  configuration option for localization
    save_opt : dictionary of  configuration option for save
    sql_opt : dictionary of  configuration option for sql


    Parameters
    ----------
    self.lAg  : list of Agent(Object)
        list of agents involved in simulation
    self.L : Layout
        Layout used in simulation   

    Notes
    ------
     All the prvious dictionnary are obtained from the chosen simulnet.ini file
    in the project directory

    """
    def __init__(self):
        SimulationRT.__init__(self) #Sympy 2
        #sympy.RealtimeEnvironment.__init__(self)  #simpy 3
        self.initialize()
        self.config = ConfigParser.ConfigParser()
        filename = pyu.getlong('simulnet.ini','ini')
        self.config.read(filename)
        self.sim_opt = dict(self.config.items('Simulation'))
        self.lay_opt = dict(self.config.items('Layout'))
        self.meca_opt = dict(self.config.items('Mechanics'))
        self.net_opt = dict(self.config.items('Network'))
        self.loc_opt = dict(self.config.items('Localization'))
        self.save_opt = dict(self.config.items('Save'))
        self.sql_opt = dict(self.config.items('Mysql'))

        self.verbose = str2bool(self.sim_opt['verbose'])
        if str2bool(self.net_opt['ipython_nb_show']):
            self.verbose = False
        self.roomlist=[]

        self.create()

    def __repr__(self):

        s = 'Simulation information' + '\n----------------------'
        s = s + '\nLayout: ' + self.lay_opt['filename'] 
        s = s + '\nSimulation duration: ' + self.sim_opt['duration'] 
        s = s + '\nRandom seed: ' + self.sim_opt['seed']
        s = s + '\nSave simulation: ' + self.save_opt['savep']

        s = s + '\n\nUpdate times' + '\n-------------' 
        s = s + '\nMechanical update: ' + self.meca_opt['mecanic_update_time'] 
        s = s + '\nNetwork update: ' + self.net_opt['network_update_time'] 
        s = s + '\nLocalization update: ' + self.net_opt['communication_mode'] 

        s = s + '\n\nAgents => self.lAg[i]' + '\n------' 
        s = s + '\nNumber of agents :' + str(len(self.lAg))
        s = s + '\nAgents IDs: ' + str([self.lAg[i].ID for i in range(len(self.lAg))])
        s = s + '\nDestination of agents choosed: ' + self.meca_opt['choose_destination']

        s = s + '\n\nNetwork' + '\n-------' 
        s = s + '\nNodes per RATs: ' + str(self.net.RAT)

        s = s + '\n\nLocalization'  + '------------' 
        s = s + '\nLocalization enable: ' + self.loc_opt['localization'] 
        s = s + '\nPostion estimation methods: ' + self.loc_opt['method'] 

        
        return s


    def create_layout(self):
        """
        Create Layout in Simpy the_world thanks to Tk backend

        """

        self.the_world = world(width = float(self.lay_opt['the_world_width']), 
                               height = float(self.lay_opt['the_world_height']),
                               scale=float(self.lay_opt['the_world_scale']))

        _filename = self.lay_opt['filename']
        self.L = Layout(_filename)



        try:
            self.L.dumpr()
            print 'Layout graphs are loaded from ',basename,'/struc/ini'
        except:
        #self.L.sl = sl
        #self.L.loadGr(G1)
            print 'This is the first time the layout file is used\
            Layout graphs are curently being built, it may take few minutes.'
            self.L.build()     
            self.L.dumpw()
        x_offset = 0  # float(self.lay_opt['x_offset'])
        y_offset = 0  # float(self.lay_opt['y_offset'])
        for ks in self.L.Gs.pos.keys():
            self.L.Gs.pos[ks] = (self.L.Gs.pos[ks][0] +
                                 x_offset, self.L.Gs.pos[ks][1] + y_offset)
        for ks in self.L.Gr.pos.keys():
            self.L.Gr.pos[ks] = (self.L.Gr.pos[ks][0] +
                                 x_offset, self.L.Gr.pos[ks][1] + y_offset)
        for ks in self.L.Gw.pos.keys():
            self.L.Gw.pos[ks] = (self.L.Gw.pos[ks][0] +
                                 x_offset, self.L.Gw.pos[ks][1] + y_offset)

        #
        # Create Layout
        #
        walls = self.L.thwall(0, 0)
        for wall in walls:
            points = []
            for ii in range(0, len(wall) - 1):
                self.the_world.add_wall(wall[ii], wall[ii + 1])

    def create_agent(self):
        """    
        create simulation's Agents

        ..todo:: change lAg list to a dictionnary ( modification in show.py too) 

        """


        self.lAg = []
        agents=[]
        Cf = ConfigParser.ConfigParser()
        Cf.read(pyu.getlong('agent.ini','ini'))
        agents=eval(dict(Cf.items('used_agent'))['list'])
        for i, ag in enumerate(agents):
            ag_opt = dict(Cf.items(ag))
            self.lAg.append(Agent(
                            ID=ag_opt['id'],
                            name=ag_opt['name'],
                            type=ag_opt['type'],
                            pos=np.array(eval(ag_opt['pos'])),
                            roomId=int(ag_opt['roomid']),
                            froom=eval(ag_opt['froom']),
                            meca_updt=float(self.meca_opt['mecanic_update_time']),
                            wait=float(ag_opt['wait']),
                            cdest=eval(self.meca_opt['choose_destination']),
                            loc=str2bool(self.loc_opt['localization']),
                            loc_updt=float(self.loc_opt['localization_update_time']),
                            loc_method=eval(self.loc_opt['method']),
                            L=self.L,
                            net=self.net,
                            epwr=dict([(eval((ag_opt['rat']))[ep],eval((ag_opt['epwr']))[ep]) for ep in range(len(eval((ag_opt['rat']))))]),
                            sens=dict([(eval((ag_opt['rat']))[ep],eval((ag_opt['sensitivity']))[ep]) for ep in range(len(eval((ag_opt['rat']))))]),
                            world=self.the_world,
                            RAT=eval(ag_opt['rat']),
                            save=eval(self.save_opt['save']),
                            gcom=self.gcom,
                            comm_mode=eval(self.net_opt['communication_mode']),
                            sim=self))
#                            
            if self.lAg[i].type == 'ag':
                self.activate(self.lAg[i].meca,
                              self.lAg[i].meca.move(), 0.0)

    def create_EMS(self):
        """
            Electromagnetic Solver object
        """
        self.EMS = EMSolver(L=self.L)

    def create_network(self):
        """
            create the whole network
        """

        self.net = Network(EMS=self.EMS)
        self.gcom=Gcom(net=self.net,sim=self)

        self.create_agent()
        # create network

        self.net.create()

        # create All Personnal networks
        for n in self.net.nodes():
            self.net.node[n]['PN'].get_RAT()
            self.net.node[n]['PN'].get_SubNet()
        self.gcom.create()


       # create Process Network
        self.Pnet = PNetwork(net=self.net,
                             net_updt_time=float(self.net_opt['network_update_time']),
                             L=self.L,
                             sim=self,
                             show_sg=str2bool(self.net_opt['show_sg']),
                             disp_inf=str2bool(self.net_opt['dispinfo']),
                             save=eval(self.save_opt['save']))
        self.activate(self.Pnet, self.Pnet.run(), 0.0)

    def create_visual(self):
        """ Create visual Tk process
        """

        self.visu = Updater(
            interval=float(self.sim_opt['show_interval']), sim=self)
        self.activate(self.visu, self.visu.execute(), 0.0)

    def create(self):
        """ Create the simulation 
            This method is called at the end of __init__

        """

        # this is just to redump the database at each simulation
        if 'mysql' in self.save_opt['save']:
            if str2bool(self.sql_opt['dumpdb']):
                os.popen('mysql -u ' + self.sql_opt['user'] + ' -p ' + self.sql_opt['dbname'] +\
                '< /private/staff/t/ot/niamiot/svn2/devel/simulator/pyray/SimWHERE2.sql' )

        ## TODO supprimer la ref en dur 
        if 'txt' in self.save_opt['save']:
            pyu.writeDetails(self)
            if os.path.isfile(basename+'/output/Nodes.txt'):
                print 'would you like to erase previous txt files ?'
                A=raw_input()
                if A == 'y':
                    for f in os.listdir(basename+'/output/'):
                        try:
                            fi,ext=f.split('.')
                            if ext == 'txt':
                                os.remove(basename+'/output/'+f)
                        except:
                            pass
                        
        self.create_layout()
        self.create_EMS()
        self.create_network()
        if str2bool(self.sim_opt['showtk']):
            self.create_visual()
        self.create_show()

        if str2bool(self.save_opt['savep']):
            self.save=Save(L=self.L,net=self.net,sim=self)
            self.activate(self.save,self.save.run(),0.0)

    def create_show(self):
        """ 
        """

        plt.ion()
        fig_net = 'network'
        fig_table = 'table'

        if str2bool(self.net_opt['show']):
            if str2bool(self.net_opt['ipython_nb_show']):
                notebook=True
            else:
                notebook =False
            self.sh = ShowNet(net=self.net, L=self.L,sim=self,fname=fig_net,notebook=notebook)
            self.activate(self.sh,self.sh.run(),1.0)

        if str2bool(self.net_opt['show_table']):
            self.sht = ShowTable(net=self.net,lAg=self.lAg,sim=self,fname=fig_table)
            self.activate(self.sht,self.sht.run(),1.0)


    def runsimul(self):
        """ Run simulation
        """

        seed(eval(self.sim_opt['seed']))
        self.simulate(until=float(self.sim_opt['duration']),
                      real_time=True,
                      rel_speed=float(self.sim_opt['speedratio']))
#        self.simulate(until=float(self.sim_opt['duration']))
        if self.save_opt['savep']:
            print 'Processing save results, please wait'
            self.save.mat_export()


if __name__ == '__main__':

    S = Simul()
    S.runsimul()

