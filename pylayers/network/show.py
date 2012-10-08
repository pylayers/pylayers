# -*- coding:Utf-8 -*-
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

import SimPy.Simulation
from SimPy.Simulation import Process, hold
import numpy as np
import scipy as sp
import networkx as nx
import random
from random import seed

import time
import matplotlib.pyplot as plt
import ConfigParser
import pkgutil

from pylayers.gis.layout import Layout

import pylayers.util.pyutil as pyu
import pylayers.util
import pdb



class ShowNet(Process):
    '''
    Show network process
   
    Process which update the plot from the dynamic information from the network

    Method
    -------
        run : run the show process

    '''

    def __init__(self, **args):
        defaults = {'net': None,
                    'L': None,
                    'fname': 'network',
                    'fig': None,
                    'ax': None,
                    'sim': None}

##       initialize attributes
        for key, value in defaults.items():
            if key in args:
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key] = value
        self.args = args
        self.fname = self.args['fname']

        Process.__init__(self, name='shownet', sim=self.args['sim'])

        self.fig = plt.figure(self.fname, figsize=(20, 5), dpi=100)
        self.fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        self.fig, self.ax = self.L.showGs(fig=self.fig)

        self.legend = True
        self.ion = True
        self.info = False

        self.coll_plot = {}
        self.C = ConfigParser.ConfigParser()
        self.C.read(pyu.getlong('show.ini','ini'))
        self.RATcolor = dict(self.C.items('RATcolor'))
        self.RATes = dict(self.C.items('RATestyle'))
        self.update = dict(self.C.items('update'))
        self.option= dict(self.C.items('option'))
        self.cpt = self.sim.now()

    def run(self):
        """ 
        Display the network

        """

        while True:
            plt.figure(self.fname)
            self.ax.axis('scaled')

            try:
                self.coll_plot['node'][1] = []
                self.coll_plot['label'][1] = []
                self.coll_plot['edge'][1] = []
                Cl = []
            except:
                self.coll_plot['node'] = [[]]
                self.coll_plot['node'].append([])
                self.coll_plot['label'] = [[]]
                self.coll_plot['label'].append([])
                self.coll_plot['edge'] = [[]]
                self.coll_plot['edge'].append([])
                Cl = []
            rloop = self.net.RAT.keys()
            for ii, rl in enumerate(rloop):
                pos = self.net.get_pos(rl)
                self.coll_plot['node'][1].append(nx.draw_networkx_nodes(self.net, pos=pos, nodelist=self.net.SubNet[rl].nodes(), node_size=100., node_color='r'))
#                if self.option['estimate']:
#                    try:
#                        pose = self.net.get_pos_est(rl)
#                        pdb.set_trace()
#                        self.coll_plot['node'][1].append(nx.draw_networkx_nodes(self.net, pos=pose, nodelist=self.net.SubNet[rl].nodes(), node_size=100., node_color='b',alpha=0.5))
#                    except:
#                        pass
                Cl = nx.draw_networkx_labels(self.net.SubNet[rl], pos=pos, font_size=10)
                self.coll_plot['label'][1].extend(Cl.values())
                self.coll_plot['edge'][1].append((nx.draw_networkx_edges(self.net, pos=pos, edgelist=self.net.SubNet[rl].edges(), width=2., alpha=0.9, edge_color=self.RATcolor[rl], style=self.RATes[rl])))

            if self.legend:
                self.ax.legend((self.coll_plot['edge'][1]), (rloop), loc=3)
            if self.info:
                L = nx.get_edge_attributes(self.net, 'TOA')

            if self.ion:
                try:
                    [jj.remove() for jj in self.coll_plot['node'][0]]
                    [jj.remove() for jj in self.coll_plot[
                        'edge'][0] if jj is not None]
                    [jj.remove() for jj in self.coll_plot['label'][0]]
                except:
                    pass
                plt.draw()
                self.coll_plot['node'][0] = self.coll_plot['node'][1]
                self.coll_plot['edge'][0] = self.coll_plot['edge'][1]
                self.coll_plot['label'][0] = self.coll_plot['label'][1]

                yield hold, self, float(self.update['meca'])


class ShowTable(Process):
    '''
    Show table process
   
    Process which update the table from the dynamic information from the network

    Method
    -------
        run : run the show process

    '''
    def __init__(self, **args):

        defaults = {'net': None,
                    'fname': 'table',
                    'lAg': None,
                    'fig': None,
                    'ax': None,
                    'sim': None}

##       initialize attributes
        for key, value in defaults.items():
            if key in args:
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key] = value
        self.args = args
        self.fname = self.args['fname']

        self.fig = plt.figure(self.fname, figsize=(30, 15), dpi=50)
        self.fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        self.ax1 = self.fig.add_subplot(111, frame_on=False)
        self.ax1.xaxis.set_visible(False)
        self.ax1.yaxis.set_visible(False)

        self.fig2 = plt.figure(
            self.fname + str(2), figsize=(30, 15), dpi=50)
        self.ax2 = self.fig2.add_subplot(111, frame_on=False)
        self.ax2.xaxis.set_visible(False)
        self.ax2.yaxis.set_visible(False)

        self.coll_plot = {}
        self.colLabels1 = ('RAT', 'Node link', 'TOA ',
                           'TOA std', 'Pr', 'Pr std', 'distance')
        self.colLabels2 = ('name', 'pos x', 'pos y ',
                           'vel x', 'vel y', 'acc x', 'acc y')
        self.cellText = []
        self.C = ConfigParser.ConfigParser()
        self.C.read(pyu.getlong('show.ini','ini'))
        self.RATcolor = dict(self.C.items('RATcolor'))
        self.RATes = dict(self.C.items('RATestyle'))
        self.update = dict(self.C.items('update'))

        Process.__init__(self, name='PShowTable', sim=self.args['sim'])

    def run(self):
        while True:

            plt.figure(self.fname)
            self.ax1.axis('scaled')
            self.ax2.axis('scaled')

            try:
                self.coll_plot['table'][1] = []
            except:
                self.coll_plot['table'] = [[]]
                self.coll_plot['table'].append([])

            self.cellText1 = []
            for rat in self.net.RAT.keys():
                T = nx.get_edge_attributes(self.net.SubNet[rat], 'TOA')
                P = nx.get_edge_attributes(self.net.SubNet[rat], 'Pr')
                D = nx.get_edge_attributes(self.net.SubNet[rat], 'd')
                for i in self.net.SubNet[rat].edges():  # boucle sur toute les liaisons
                    r = [str(rat), str(i), "%.2f" % (T[i][0]), "%.2f" % (T[i][1]), "%.2f" % (P[i][0]), "%.2f" % (P[i][1]), "%.2f" % (D[i])]
                    self.cellText1.append(r)

            self.coll_plot['table'][1] = self.ax1.table(cellText=self.cellText1, colLabels=self.colLabels1, loc='center')
            self.coll_plot['table'][1].scale(2, 5)
            self.coll_plot['table'][1].auto_set_font_size(False)
            self.coll_plot['table'][1].set_fontsize(18)

#                self.coll_plot['table'][1].auto_set_font_size(False)
#                self.coll_plot['table'][1].set_fontsize(20)
            try:
                self.coll_plot['table'][0].remove()
            except:
                pass

            plt.draw()
            self.coll_plot['table'][0] = self.coll_plot['table'][1]

            try:
                self.coll_plot['table2'][1] = []
            except:
                self.coll_plot['table2'] = [[]]
                self.coll_plot['table2'].append([])

            plt.figure(self.fname + str(2))
            self.cellText2 = []

            for a in self.lAg:
                if a.type != 'ap':
                    r2 = [str(a.ID), "%.2f" % (a.meca.position[0]), "%.2f" % (a.meca.position[1]), "%.2f" % (a.meca.velocity[0]), "%.2f" % (a.meca.velocity[1]), "%.2f" % (a.meca.acceleration[0]), "%.2f" % (a.meca.acceleration[1])]
                    self.cellText2.append(r2)

            self.coll_plot['table2'][1] = self.ax2.table(cellText=self.cellText2, colLabels=self.colLabels2, loc='center')
            self.coll_plot['table2'][1].scale(2, 5)
            self.coll_plot['table2'][1].auto_set_font_size(False)
            self.coll_plot['table2'][1].set_fontsize(18)

            try:
                self.coll_plot['table2'][0].remove()
            except:
                pass

            plt.draw()
            self.coll_plot['table2'][0] = self.coll_plot['table2'][1]


            yield hold, self, float(self.update['table'])


