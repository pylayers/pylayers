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

import SimPy.Simulation
from SimPy.Simulation import Simulation, Process, hold
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
from pylayers.network.show import ShowNet, ShowTable
from pylayers.mobility.agent import Agent
#import pylayers.Network.Localization
from pylayers.network.emsolver import EMSolver
from pylayers.gis.layout import Layout
from pylayers.antprop.slab import Slab
from pylayers.util.utilnet import str2bool
from pylayers.mobility.transit.World import world
from pylayers.util.pymysqldb import Database as DB
from pylayers.util.project import *

import pdb
import os

class Simul(Simulation):

    def __init__(self):
        Simulation.__init__(self)
        self.initialize()
        self.config = ConfigParser.ConfigParser()
        self.config.read(pyu.getlong('simulnet.ini','ini'))
        self.sim_opt = dict(self.config.items('Simulation'))
        self.ag_opt = dict(self.config.items('Agent'))
        self.lay_opt = dict(self.config.items('Layout'))
        self.meca_opt = dict(self.config.items('Mecanic'))
        self.net_opt = dict(self.config.items('Network'))
        self.loc_opt = dict(self.config.items('Localization'))
        self.save_opt = dict(self.config.items('Save'))
        self.sql_opt = dict(self.config.items('Mysql'))

    def create_layout(self):
        """
        Create Layout in Simpy the_world thantks to Tk backend

        TODO
        ----
        automatically update Simpy the_world in regards of the Pyray Layout

        """

        self.the_world = world(width=float(self.lay_opt['the_world_width']), height=float(self.lay_opt['the_world_height']), scale=float(self.lay_opt['the_world_scale']))
        tk = self.the_world.tk
        canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
        canvas.create_rectangle(x_(-1), y_(-1), x_(100), y_(100), fill='white')

        _filename = self.lay_opt['filename']
        #sl=Slab.SlabDB(self.lay_opt['slab'],self.lay_opt['slabmat'])
        #G1   = Graph.Graph(sl=sl,filename=_filename)
        self.L = Layout()
        if _filename.split('.')[1] == 'str':
            self.L.loadstr(_filename)
        else:
            self.L.loadstr2(_filename)

        try:
            self.L.dumpr()
        except:
        #self.L.sl = sl
        #self.L.loadGr(G1)
            self.L.buildGt()
            self.L.buildGr()
            self.L.buildGw()
            self.L.buildGv()

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
            for point in wall:
                points.append(x_(point[0]))
                points.append(y_(point[1]))
            canvas.create_polygon(points, fill='maroon', outline='black')
            for ii in range(0, len(wall) - 1):
                self.the_world.add_wall(wall[ii], wall[ii + 1])

    def create_agent(self, init='random'):
        """    create simulation's Agents

        """

        if init == 'random':
            self.lAg = []
            agents = ['A1', 'A2', 'A3','BS1','BS2']
#            agents=['A1','A2','A3','A4' ]
            Cf = ConfigParser.ConfigParser()
            Cf.read(pyu.getlong('agent.ini','ini'))


            for i, ag in enumerate(agents):
                ag_opt = dict(Cf.items(ag))
                self.lAg.append(Agent(
                                ID=ag_opt['id'],
                                name=ag_opt['name'],
                                type=ag_opt['type'],
                                pos=np.array(eval(ag_opt['pos'])),
                                roomId=int(ag_opt['roomid']),
                                meca_updt=float(self.meca_opt['mecanic_update_time']),
                                loc=str2bool(self.loc_opt['localization']),
                                loc_updt=float(self.loc_opt['localization_update_time']),
                                Layout=self.L,
                                net=self.net,
                                world=self.the_world,
                                RAT=eval(ag_opt['rat']),
                                save=eval(self.save_opt['save']),
                                sim=self))

                if self.lAg[i].type == 'ag':
                    self.activate(self.lAg[i].meca,
                                  self.lAg[i].meca.move(), 0.0)

#            R=[0,9,8,18,12,1,2,3,4,5,6,7,8,10,11,12,14,15]
#            #rat = [['rat1','rat6'],['rat1','rat2'],['rat2','rat3'],['rat3','rat4'],['rat4','rat5'],['rat5','rat6'],['rat6','rat1']]
#            rat = [['lte','bt'],['wifi','lte','bt'],['wifi','bt'],['wifi','bt'],['lte'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],['wifi','lte','bt'],]
#            type=['ag','ag','ap','ap','ap']
#            self.lAg=[]
#            for i in range(int(self.ag_opt['number_of_agent'])):
#                self.lAg.append(Agent.Agent(
#                                            ID=i,
#                                            type=type[i],
#                                            pos=[],
#                                            roomId=R[i],
#                                            meca_updt=float(self.meca_opt['mecanic_update_time']),
#                                            loc=Util.str2bool(self.loc_opt['localization']),
#                                            loc_updt=float(self.loc_opt['localization_update_time']),
#                                            Layout=self.L,
#                                            net=self.net,
#                                            RAT=rat[i],
#                                            sim=self))
#                if self.lAg[i].type == 'ag':
#                    self.activate(self.lAg[i].meca,self.lAg[i].meca.move(),0.0)

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

        self.create_agent()
        # create network
        self.net.create()
        # create All Personnal networks
        for n in self.net.nodes():
            self.net.node[n]['PN'].get_RAT()
            self.net.node[n]['PN'].get_SubNet()

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
        """ Create the simulation, to be ready to run

        """

        # this is just to redump the database at each simulation
        if 'mysql' in self.save_opt['save']:
            if str2bool(self.sql_opt['dumpdb']):
                os.popen('mysql -u ' + self.sql_opt['user'] + ' -p ' + self.sql_opt['dbname'] +\
                '< /private/staff/t/ot/niamiot/svn2/devel/simulator/pyray/SimWHERE2.sql' )

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
        

    def create_show(self):
        plt.ion()
        fig_net = 'network'
        fig_table = 'table'

        if str2bool(self.net_opt['show']):
            self.sh=ShowNet(net=self.net, L=self.L,sim=self,fname=fig_net)
            self.activate(self.sh,self.sh.run(),1.0)

        if str2bool(self.net_opt['show_table']):
            self.sht=ShowTable(net=self.net,lAg=self.lAg,sim=self,fname=fig_table)
            self.activate(self.sht,self.sht.run(),1.0)


    def runsimul(self):
        """ Run simulation
        """
        self.create()
        self.simulate(until=float(self.sim_opt['duration']))


if __name__ == '__main__':

    seed(0)
    S = Simul()
    S.runsimul()
