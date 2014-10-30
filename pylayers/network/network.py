# -*- coding:Utf-8 -*-
"""

.. currentmodule:: pylayers.network.network

Node Class
==========

.. autosummary::
    :toctree: generated/

     Node.__init__
     Node.randomMAC


Network Class
==============

.. autosummary::
    :toctree: generated/

    Network.__init__
    Network.__repr__


Network creation
----------------

.. autosummary::
    :toctree: generated/

    Network.add_devices
    Network.create
    Network._get_edges_typ
    Network._get_grp
    Network._get_llink
    Network._get_wstd
    Network._get_SubNet
    Network._connect
    Network._init_PN


Network attributes queries
--------------------------

.. autosummary::
    :toctree: generated/

    Network.get_pos
    Network.get_orient
    Network.get_pos_est
    Network.overview
    Network.haspe
    Network.pp


Network update
--------------

.. autosummary::
    :toctree: generated/

    Network.update_pos
    Network.update_orient
    Network.update_edges
    Network.update_PN
    Network.compute_LDPs
    Network.update_LDPs


Network Utilities
-----------------

.. autosummary::
    :toctree: generated/

    Network.perm
    Network.combi
    Network.Gen_tuple
    Network.dist_edge
    Network.show



Network save
------------

.. autosummary::
    :toctree: generated/

    Network.csv_save
    Network.init_save
    Network.mat_save
    Network.txt_save
    Network.loc_save
    Network.pyray_save
    Network.loc_save
    Network.ini_save


PNetwork Class
==============

    SimPy Process compliant version of the Network class

.. autosummary::
    :toctree: generated/

     PNetwork.__init__
     PNetwork.run

"""
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
import scipy as sp 
import networkx as nx
import itertools
import pickle as pk
import pkgutil
from copy import deepcopy


import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#import Tkinter as tk

import ConfigParser
import copy
import pdb
#from PyLayers.Network.Node import Node
import pylayers.util.pyutil as pyu
from pylayers.network.emsolver import EMSolver
from pylayers.network.show import ShowNet,ShowTable
#from pylayers.util.pymysqldb import Database
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.util.utilnet import str2bool
import time
import pylayers.util.pyutil as pyu


from SimPy.SimulationRT import Process,hold
import pprint
import select
import sys

try:
    from mayavi import mlab
    from tvtk.tools import visual

except:
    print 'mayavi not installed'


# How to take into account  1 specific key specifique for 1 MultiGraph
# MULTIGRAPH !!!! G.add_edge(10,11,key='wifi',attr_dict=dict(Pr=0,TOA=10))

class Node(PyLayers,nx.MultiGraph):
    """ Class Node

    inherit from networkx.MultiGraph()

    Attributes
    ----------

    Id    : float/hex/str/...
            node Id
    p    : np.array
            True position
    t    : time.time()
            Tag time
    wstd    : list
            available wstd of the node
    PN    : Network.Network
            Personal vision of the Network
    pos    : Dictionnary
            parser from Node.Node to networkx.node.pos

    Method
    ------

    RandomMac(): Generate a RAndom Mac adress

    """
    def __init__(self,**kwargs):
        nx.MultiGraph.__init__(self)

        defaults = { 'ID':0,
                     'name':'',
                     'p':np.array(()),
                     't':0.,
                     'pe':np.array(()),
                     'te':0.,
                     'wstd':[],
                     'epwr':{},
                     'sens':{},
                     'typ':'ag',
                     'grp':'',
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        # Personnal Network init
        self.ID = kwargs['ID']
        self.PN = Network(owner=self.ID,PN=True)
        self.PN.add_node(self.ID,dict(ID=kwargs['ID'],
                                      name=kwargs['name'],
                                      pe=kwargs['pe'],
                                      te=kwargs['te'],
                                      wstd=kwargs['wstd'],
                                      epwr=kwargs['epwr'],
                                      sens=kwargs['sens'],
                                      typ=kwargs['typ'],
                                      ))
        # Network init
        self.add_node(self.ID,dict(ID=kwargs['ID'],
                              name=kwargs['name'],
                              PN=self.PN,
                              p=kwargs['p'],
                              pe=self.PN.node[self.ID]['pe'],
                              t=kwargs['t'],
                              wstd=kwargs['wstd'],
                              epwr=kwargs['epwr'],
                              sens=kwargs['sens'],
                              typ=kwargs['typ'],
                              grp=kwargs['grp']))
        self.p   = self.node[self.ID]['p']
        self.pe  = self.PN.node[self.ID]['pe']
        self.t   = self.node[self.ID]['t']
        self.wstd = self.node[self.ID]['wstd']
        self.epwr = self.node[self.ID]['epwr']
        self.sens = self.node[self.ID]['sens']


    def randomMAC(self):
        """ Generate a random MAC address

        Returns
        -------

        macadress : string

        """
        mac = [ 0x00, 0x16, 0x3e,
        random.randint(0x00, 0x7f),
        random.randint(0x00, 0xff),
        random.randint(0x00, 0xff) ]
        return ':'.join(map(lambda x: "%02x" % x, mac))


class Network(PyLayers,nx.MultiDiGraph):
    """ Network class

    inherits from networkx.Graph()

    Attributes
    ----------

    wstd : dictionnary
        keys  = wstd
        value = list of nodes id
    wstde : dictionnary
        keys  = wstd
        value = list of edges id 
    SubNet : dictionnary
        keys  = wstd
        value = Subgraph of the given wstd
    pos : dictionnary
        keys  = node id
        value = node position

    Methods
    -------

    _get_wstd(self)  : Get wstd from nodes of the network
    _connect(self)  : Connect each node from a wireless standard
    create(self)   : compute get_wstd(),get_pos() and connect()
    update_LDP(self,n1,n2,wstd,LDP=None,value=[])    : update Location Dependent Parameter  
    compute_LDP(self,wstd) : compute the LDP value thanks to a ElectroMag Solver 
    update_pos(self,n,p=np.array)  : update node (or node list) position
    get_pos(self,wstd=None)         : get node positions
    pp(self)                    : pretty print on std out all edtges informations
    show(rat=None,legend=True)     : Display network for all rat or specified in Rat. 

    """



    def __init__(self,owner='sim',EMS=EMSolver(),PN=False):
        """ object constructor

        Parameters
        ----------

        owner : string
            'sim' |
        EMS : EMSolver
        PN : Boolean
            personal network activation

        """
        nx.MultiDiGraph.__init__(self)
        self.owner=owner
        self.wstd={}
        self.LDP = ['TOA','Pr']
        self.SubNet={}
        self.grp={}
        self.EMS=EMS
        self.coll_plot={}
        self.pos={}
        self.mat={}
        self.links={}
        self.relinks={}
        self.idx = 0
        self.lidx = 0
        self.isPN=PN

    def __repr__(self):

        if not self.isPN:
            s = 'Network information\n*******************\n'
            s = s + 'number of nodes: ' + str(len(self.nodes())) +'\n'
            title = '{0:7} | {1:15} |{2:7} | {3:4} | {4:17} | {5:10} '.format('ID', 'name', 'group', 'type', 'position (x,y,z)','wstd')
            s = s + title + '\n' + '-'*len(title) + '\n'
            for n in self.nodes():
                # for compliance with simulnet and simultraj
                # to be merged
                try:
                    wstd = self.node[n]['wstd'].keys()
                except:
                    wstd = self.node[n]['wstd']
                s = s + '{0:7} | {1:15} |{2:7} | {3:4} | {4:5.2f} {5:5.2f} {6:5.2f} | {7:10} '\
                .format(self.node[n]['ID'][:7], self.node[n]['name'][:15],
                self.node[n]['grp'][:7], self.node[n]['typ'][:4], self.node[n]['p'][0],
                self.node[n]['p'][1],self.node[n]['p'][2],wstd[:10]) + '\n'

    #             try:
    #                 s = s + 'node ID: ' + str(self.node[n]['ID']) + '\n'
    #             except:
    #                 s = s + 'node ID: ' + str(n) + '\n'
    #             try :
    #                 s = s + 'wstd: ' + str(self.node[n]['wstd'].keys()) + '\n'
    #             except:
    #                 s = s + 'wstd: ' + str(self.node[n]['wstd']) + '\n'
    #             try:
    #                 s = s + 'grp: ' + str(self.node[n]['grp']) + '\n'
    #             except:
    #                 s = s + 'type: ' + str(self.node[n]['typ']) + '\n'
    #             try:
    #                 s = s + 'pos: ' + str(self.node[n]['p']) + '\n'
    #             except:
    #                 pass
    #             s = s + '\n'
    #         # typ = nx.get_node_attributes(self,'typ').values() 

    #         # nodes = np.array(nx.get_node_attributes(self,'typ').items())

    #         # nb_ag = len(np.where(nodes=='ag')[0])
    #         # nb_ap = len(np.where(nodes=='ap')[0])

    #         # pag=np.where(nodes=='ag')
    #         # pap=np.where(nodes=='ap')

    #         # s = s +  '\n' + str(nb_ag) + ' Mobile Agents\n  -------------\n'
    #         # s = s + 'Agents IDs : ' + str([nodes[i,0] for i in pag[0]]) +'\n'


    #         # s = s +  '\n' + str(nb_ap) + ' Access points\n  -------------\n'
    #         # s = s + 'number of access point  : ' + '\n'
    #         # s = s + 'access points  IDs : ' + str([nodes[i,0] for i in pap[0]]) +'\n'

    #         # if len(self.SubNet.keys()) != 0 :
    #         #     s = s + '\n\nSubNetworks :' +str(self.SubNet.keys()) + '\n===========\n'
    #         #     for sub in self.SubNet.keys():
    #         #         s = s + '\t'+ sub + '\n' +  self.SubNet[sub].__repr__() + '\n'

        else:
            s = 'Personnal Network of node ' +str(self.owner)+ ' information\n***************************************\n'
            s = s + '{0:7} |{1:20} | {2:5} | {3:7}| {4:7}| {5:7}| {6:7}| {7:7}| {8:10}|'.format('peer','wstd', 'TOA','std TOA','tTOA', 'Pr', 'std Pr', 'tPr','visibility')
            for e1,e2 in self.edges():
                for r in self.edge[e1][e2].keys():
                    TOA = self.edge[e1][e2][r]['TOA'][0]
                    stdTOA = self.edge[e1][e2][r]['TOA'][1]
                    pr = self.edge[e1][e2][r]['Pr'][0]
                    stdpr = self.edge[e1][e2][r]['Pr'][1]
                    try :
                        tTOA = self.edge[e1][e2][r]['tTOA']
                    except:
                        tTOA = 'nan'
                    try :
                        tpr = self.edge[e1][e2][r]['tPr']
                    except:
                        tpr = 'nan'
                    vis = self.edge[e1][e2][r]['vis']
                    np.set_printoptions(precision=3)

                    s = s + '\n' + '{0:7} |{1:20} | {2:5.2f} | {3:7.2f}| {4:7}| {5:7.2f}| {6:7.2f}| {7:7}| {8:10}|'.format(e2 ,r ,TOA ,stdTOA ,tTOA ,pr , stdpr ,tpr, vis)
        return s



    def add_devices(self, dev, p=[], grp=''):
        """ add devices to the current network

        Parameters
        ----------

        dev : list
            list of Devices
        p : ndarray (Ndev x 3)
            np.array of devices' positions
        grp : string
            name of the group of device belong to.

        """

        if not isinstance(dev,list):
            dev=[dev]
        if p == []:
            p = np.nan*np.zeros((len(dev),3))
        elif len(p.shape) == 1:
            p = p.reshape(1,3)

        if (p.shape[0] != len(dev)):
            raise AttributeError('number of devices != nb pos')

        # check if unique ID (in dev and in network ) else raise error
        ids = [d.ID for d in dev]
        for d in dev:
            if d.ID in self:
                raise AttributeError('Devices must have a different ID')


        # add spectific node informations
        if 'ap' in grp:
            typ = 'ap'
        else :
            typ = 'ag'

        [d.__dict__.update({'p': p[ud, :],
              'T': np.eye(3),
              'grp':grp,
              'typ':typ,
              'dev':d,
                    }) for ud, d in enumerate(dev)]
# 
        # self.add_nodes_from([(d.ID, ldic[ud]) for ud,d in enumerate(dev)])
        self.add_nodes_from([(d.ID, d.__dict__) for d in dev])

        # create personnal network
        for ud, d in enumerate(dev):
            self.node[d.ID]['PN']= Network(owner=d.ID, PN=True)
            self.node[d.ID]['PN'].add_nodes_from([(d.ID,d.__dict__)])
        self._get_wstd()
        # for d in dev:
        #     for s in d.wstd.keys():
        #         try:
        #             self.wstd[s]
        #             if d.ID not in self.wstd[s]:
        #                 self.wstd[s].append(d.ID)
        #         except :
        #             self.wstd[s]=[d.ID]


    def perm(self,iterable,r,key,d=dict()):
        """  calculate permutation

        Notes
        -----

        combi = itertools.permutation(iterable,r) adapted

        This is an adapted version of itertools.permutations to
        comply with the networkx.add_edges_from method.

        itertools.permutations(range(4), 3) --> 012 013 021 023 031 302 102 103 ...

        self.perm([10,11],2,'wifi') -->     (10, 11, 'wifi', {'Pr': [], 'TOA': []}) 
                                (11, 10, 'wifi', {'Pr': [], 'TOA': []})

        Parameters
        ----------

        iterable : list
            list of node
        r     : int
            number of node gathered in the output tuple ( always set 2 ! )


        Returns
        --------

        out : tuple(node_list,r,wstd,d):
        node_list    : list of node1
        r        : gather r node in the tuple
        wstd        : the specified wstd
        d        : dictionnary of wstd attribute

        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>> l= [0,1,2]
        >>> key='toto'
        >>> d=dict(key1=1,key2=2)
        >>> perm=N.perm(l,2,key,d)
        >>> perm.next()
        (0, 1, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (0, 2, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (1, 0, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (1, 2, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (2, 0, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (2, 1, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        """

#        for l in self.LDP:
#            d[l]=[]
        pool = tuple(iterable)
        n = len(pool)
        r = n if r is None else r
        if r > n:
            return
        indices = range(n)
        cycles = range(n, n-r, -1)
        yield tuple((pool[indices[0]],pool[indices[1]],key,d))
        while n:
            for i in reversed(range(r)):
                cycles[i] -= 1
                if cycles[i] == 0:
                    indices[i:] = indices[i+1:] + indices[i:i+1]
                    cycles[i] = n - i
                else:
                    j = cycles[i]
                    indices[i], indices[-j] = indices[-j], indices[i]
                    yield tuple((pool[indices[0]],pool[indices[1]],key,d))
                    break
            else:
                return



    def combi(self,iterable,r,key,d=dict()):
        """ calculate combination

        Notes
        -----

        combi = itertools.combination(iterable,r) adapted

        This is an adapted version of itertools.combinations in order
        to comply with the networkx.add_edges_from method.
        itertools.combinations('ABCD', 2) --> AB AC AD BC BD CD
        itertools.combinations(range(4), 3) --> 012 013 023 123

        self.combi([10,11,12],2,'wifi') -->     (10, 11, 'wifi', {'Pr': [], 'TOA': []}) 
                                (10, 12, 'wifi', {'Pr': [], 'TOA': []})
                            (11, 12, 'wifi', {'Pr': [], 'TOA': []})

        Parameters
        ----------

        iterable : list
            list of node
        r     : int
            number of node gathered in the output tuple ( always set 2 ! )
        d  : dict


        Returns
        -------

        out : tuple(node_list,r,wstd,d):
        node_list    : list of node1
        r        : gather r node in the tuple
        wstd        : the specified wstd
        d        : dictionnary of wstd attribute

        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>> l= [0,1,2,3,4]
        >>> key='toto'
        >>> d=dict(key1=1,key2=2)
        >>> comb=N.combi(l,2,key,d)
        >>> comb.next()
        (0, 1, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> comb.next()
        (0, 2, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})

        """

#        for l in self.LDP:
#            d[l]=[]



        pool = iterable
        n = len(pool)
        if r > n:
            return
        indices = range(r)

        yield tuple((pool[indices[0]],pool[indices[1]],key,d))
        while True:

            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            for j in range(i+1, r):
                indices[j] = indices[j-1] + 1
            yield tuple((pool[indices[0]],pool[indices[1]],key,d))


    def Gen_tuple(self,gene,wstd,var):
        """ generate a specific tuple

        Parameters
        ----------

        gene : tuple(x,y) iterator
        wstd  : str
        var  : list
            len(var) = len(gene)

        Yield
        -----

        tuple : (gene[i][0],gene[i][1],wstd,var[i]) for iteration i

        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>> tup = zip(range(5),range(5))
        [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]
        >>> g = iter(tup)
        >>> wstd='string wstd'
        >>> var=[10,11,12,13,14]

        >>> T=N.Gen_tuple(g,wstd,var)
        >>> T.next()
        (0, 0, 'string wstd', 10)
        >>> T.next()
        (1, 1, 'string wstd', 11)

        """

        gvar=iter(var)
        while True:
            G=gene.next()
            Gvar=gvar.next()
            yield(tuple((G[0],G[1],wstd,Gvar)))


    def _get_wstd(self):
        """ get wireless standards from nodes of the network

        wstd argument specifies which wireless standard to append to the network.
        If None, all wireless standards are appended.

        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>>    N=Network.Network()
        >>> for i in range(3):
                no = Node(ID=i,wstd=['wifi','bt'])
                N.add_nodes_from(no.nodes(data=True))

        >>> N._get_wstd()
        {'bt': [0, 1, 2], 'wifi': [0, 1, 2]}

        """
#        if Rat !='None' :
        for no in self.nodes():
            for r in self.node[no]['wstd']:
                try:
                    self.wstd[r].append(no)
                except :
                    self.wstd[r]=[no]
#        else :
#            for no in self.nodes(): 
#                for r in self.node[no]['wstd']:
#                    try:
#                        self.wstd[r].extend(no)
#                    except :
#                        self.wstd[r]=[no]

        # uniquify results
        for ws in self.wstd.keys():
            self.wstd[ws]    = {}.fromkeys(self.wstd[ws]).keys()


    def update_edges(self, d , wstd, nodes=[]):
        """ update edges information for a given wstd

        Parameters
        ----------

        d: dict :
            dictionnary of information to be updated
        wstd : list | dict
            list of wstd where d has to be modified
        nodes : list
            list of nodes where information has to be applied 
            raise error if nodes in the list are not in wstd
        """
        if isinstance(wstd,dict):
            wstd = wstd.keys()
        elif not isinstance(wstd, list):
            wstd = [wstd]

        for w in wstd:
            if nodes == []:
                edges=self.perm(self.wstd[w], 2, w, d=d)
            else:
                nin = [n in self.wstd[w] for n in nodes]
                # raise error if some nodes are note in the wstd
                # no error raised if none nodes in wstd
                if sum(nin) != len(nodes) and (sum(nin) != 0):
                    unin = np.where(np.array(nin) == False)[0]
                    raise AttributeError(str(np.array(nodes)[unin]) +' are not in ' + w)
                else :
                    edges=self.perm(nodes, 2, w, d=d)
            try:
                self.SubNet[w].add_edges_from(edges)
            except:
                self.add_edges_from(edges)



    def _connect(self):
        """ connect nodes

        This method 
        1) Connect all nodes from the network sharing the same wstd 
        2) Create the associated SubNetworks
        3) Create lists of links : self.links and self.relinks

        """

        edge_dict={}
        for l in self.LDP:
            edge_dict[l]=np.array((np.nan, np.nan))
        edge_dict['vis'] = False


        for wstd in self.wstd.keys():
            self.update_edges(edge_dict,wstd)
            self._get_SubNet(wstd)

        # update  edges type informatiosn
        self._get_edges_typ()
        # create lists of links
        self._get_llinks()

    def _get_llinks(self):
        """ get list of links from the Network

        Notes
        -----

        Fill self.links and self.relinks
        """

        for wstd in self.wstd.keys():
            self.links[wstd]=[]
            self.relinks[wstd]=[]
            for i in itertools.combinations(self.wstd[wstd],2):
                self.links[wstd].append([i[0],i[1],self.edge[i[0]][i[1]][wstd]['typ']])
                # if self.node[i[0]]['grp'] == self.node[i[1]]['grp']\
                #     and (self.node[i[0]]['typ'] != 'ag'\
                #         or self.node[i[0]]['typ'] != 'ag'):
                #     self.links[wstd].append([i,'OB'])
                # else :
                #     nx.set_edge_attributes(self,i,{'typ':'OffB'})
                #     self.links[wstd].append([i,'OffB'])

            self.relinks[wstd]=[[i[1],i[0],i[2]] for i in self.links[wstd]]

    def _get_edges_typ(self):
        """ apply specific type on edges

        Notes
        -----

        types are :
            OB : On body
             when link' nodes of a link are:
                on the same agent
                and belong to the same group
            B2B : Body to Body
                when link' nodes of a link are:
                    between 2 agents
            B2I : Body to Infrastructure
                when link' nodes of a link are:
                    between an agent and an access point
            I2I : Infrastructure to Infrastructure
                when link' nodes of a link are:
                    between 2 access points
        """
        d = {}
        for n in self.SubNet:
            for e in self.SubNet[n].edges():
                e0 = self.node[e[0]]
                e1 = self.node[e[1]]
                if e0['typ'] == e1['typ'] == 'ag':
                    if e0['grp'] == e1['grp']:
                        self.update_edges({'typ': 'OB'}, n, e)
                    else :
                        self.update_edges({'typ': 'B2B'}, n, e)
                elif e0['typ'] == e1['typ'] == 'ap':
                    # if e0['grp'] == e1['grp']:
                    self.update_edges({'typ': 'I2I'}, n, e)
                        # print str(e0['ID']),str(e1['ID']),'I2I'
                else:    
                    self.update_edges({'typ': 'B2I'}, n, e)

    def _get_grp(self):
        """ 
            get group of the nodes of a network

        """

        for n in self.nodes():
            grp = self.node[n]['grp']
            if grp not in self.grp.keys():
                self.grp[grp] = []
            if n not in self.grp[grp]:
                self.grp[grp].extend([n])

    def _get_SubNet(self,wstd=[]):
        """
        get SubNetworks of a network

        Warnings
        --------

        ALWAYS use self._get_wstd() BEFORE !

        Parameters
        ----------

        wstd : specify which SubNet to create

        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>> for i in range(2):
                no = Node.Node(ID=i,wstd=['wifi','bt'])
                N.add_nodes_from(no.nodes(data=True))

        >>> no = Node.Node(ID=2,wstd=['wifi'])
        >>>    N.add_nodes_from(no.nodes(data=True))
        >>> N._get_wstd() # VERY IMPORTANT 

        >>> N._get_SubNet()
        >>> N.SubNet['bt'].nodes()
        [0, 1]
        >>> N.SubNet['wifi'].nodes()
        [0, 1, 2]


        """
        if wstd == []:
            for wstd in self.wstd:            
                # creating all SubNetworks 
                self.SubNet[wstd]= self.subgraph(self.wstd[wstd])
                # remove information from previous subnetwork (because subgraph copy the whole edge information)
                ek = self.SubNet[wstd].edges(keys=True)
                for e in ek :
                    if e[2] != wstd:
                        self.SubNet[wstd].remove_edge(e[0],e[1],e[2])
                for n in self.SubNet[wstd].nodes():
                    try:
                        self.SubNet[wstd].node[n]['epwr']=self.SubNet[wstd].node[n]['epwr']
                        self.SubNet[wstd].node[n]['sens']=self.SubNet[wstd].node[n]['sens']
                    except: 
                        pass


        elif wstd in self.wstd:
            # creating SubNetworks
            self.SubNet[wstd]= self.subgraph(self.wstd[wstd])
            # remove information from previous subnetwork (because subgraph copy the whole edge information)
            for k in self.wstd.keys():
                if k != wstd:
                    try:
                        self.SubNet[wstd].remove_edges_from(self.SubNet[k].edges(keys=True))
                    except :
                        pass
                for n in self.SubNet[wstd].nodes():
                    try:
                        self.SubNet[wstd].node[n]['epwr']=self.SubNet[wstd].node[n]['epwr']
                        self.SubNet[wstd].node[n]['sens']=self.SubNet[wstd].node[n]['sens']
                    except: 
                        pass


        else :
            raise AttributeError('invalid wstd name')


    def _init_PN(self):
        """
        Initializing personnal networks

        """


        for wstd, subnet in self.SubNet.iteritems():
            for n in subnet.nodes():
                for nn in subnet.nodes():
                    if nn != n:
                        try:
                            if wstd not in subnet.node[n]['PN'].node[nn]['wstd']: 
                                subnet.node[n]['PN'].node[nn]['wstd'].append(wstd)
                        except:
                            subnet.node[n]['PN'].add_node(nn,attr_dict=dict(wstd=[wstd],pe=np.array(()),te=0.),typ=subnet.node[nn]['typ'])


                Z= subnet.edges(n,keys=True,data=True)
                subnet.node[n]['PN'].add_edges_from(Z)


    def create(self):
        """ create the network

        This method computes :
            * _get_wstd()
            * _get_grp()
            * _connect()
            * _init_PN
            



        """

        self._get_wstd()
        self._get_grp()
        self._connect()
        self._init_PN()

    def update_PN(self):
        """ update personnal network
        """
        ####################################################################################
        # first iteration requested to correctely initiatilzing Personnal Networks's Subnets 
        for wstd in self.wstd.iterkeys():
            for ldp in self.LDP:
                self.compute_LDPs(self.nodes(),wstd)
        for n in self.nodes():
            self.node[n]['PN']._get_wstd()
            self.node[n]['PN']._get_SubNet()
            # Add access point position in each personal network (PN)
            [self.node[n]['PN'].node[n2].update({'pe':self.node[n2]['p']}) for n2 in self.node[n]['PN'].node.iterkeys() if self.node[n]['PN'].node[n2]['typ'] == 'ap']

        ####################################################################################


#    def visibility(func):
#        def wrapper(*args, **kwargs):
#            a = list(args)
#            pdb.set_trace()
#            print 'decorator',a
#            return func(*args, **kwargs)
#        return wrapper
    def dist_edge(self,e,dp):
        """ compute distance to edge

        Parameters
        ----------

        e :
        dp:

        """
        return(np.array([np.sqrt(np.sum((dp[i[0]]-dp[i[1]])**2)) for i in e]))


    def update_LDPs(self,ln,wstd,lD):   
        """Set a value between 2 nodes (n1 and n2) for a specific LDP from a wstd

        This method update :     * The network edges 
                    * The personal network (PN) of both n1 and n2

        Parameters
        ----------

        n1      : node ID
        n2      : node ID
        wstd     : string
            A specific wstd which exist in the network ( if not , raises an error)
        ln     : list 
            list of nodes
        lD     : list of dictionnary:
            [ {LDP1_1:[value , std],LDP2_1:[value , std] } , {LDPL_N:[value , std],LDPL_N:[value , std] }    ] for N nodes and L LDPS

        .. toto::

            Check if LDP value is compliant with the LDP


        """
        self.SubNet[wstd].add_edges_from(self.Gen_tuple(ln,wstd,lD))



    def compute_LDPs(self,wstd):
        """compute edge LDP

        Parameters
        ----------

        wstd     : string
            A specific wstd which exists in the network ( if not , raises an error)

        """
        # value    : list : [LDP value , LDP standard deviation] 
        # method    : ElectroMagnetic Solver method ( 'direct', 'Multiwall', 'PyRay'

        p=nx.get_node_attributes(self.SubNet[wstd],'p')
        epwr=nx.get_node_attributes(self.SubNet[wstd],'epwr')
        sens=nx.get_node_attributes(self.SubNet[wstd],'sens')
        e=self.links[wstd]#self.SubNet[wstd].edges()
        re=self.relinks[wstd] # reverse link aka other direction of link

        lp,lt, d, v= self.EMS.solve(p,e,'all',wstd,epwr,sens)
        lD=[{'Pr':lp[i],'TOA':lt[np.mod(i,len(e))] ,'d':d[np.mod(i,len(e))],'vis':v[i]} for i in range(len(d))]
        self.update_LDPs(iter(e+re),wstd,lD)


    def update_orient(self, n, T, now=0.):
        """
        Update Orientation(s) of a Device(s)/node(s)

        Parameters
        ----------

        n      : float/string (or a list of) 
            node ID (Nn x 3)
        T    : np.array  ( or a list of )
            node orientation (Nn x 3 x 3)

        Todo
        ----

        update the orientation of the antenna in the ACS (for now only DCS is updated)

        """

        if (isinstance(T,np.ndarray)) or (isinstance(n,list) and isinstance(T,list) ):
            # Tranfrom input as list
            if not(isinstance(n,list)):
                n=[n]
                T=[T]
            if len(n) == len(T):    
                d=dict(zip(n,T))    # transform data to be complient with nx.set_node_attributes            
                nowd=dict(zip(n,[now]*len(n)))
            else :
                raise TypeError('n and T must have the same length')
            # update position
            nx.set_node_attributes(self,'T',d)        
            # update time of ground truth position
            nx.set_node_attributes(self,'t',nowd)
        else :
            raise TypeError('n and p must be either: a key and a np.ndarray, or 2 lists')

    def update_pos(self, n, p, now=0., p_pe='p'):
        """
        Update Position(s) of Device(s)/node(s)

        Parameters
        ----------

        n      : float/string (or a list of)
            node ID
        p    : np.array  ( or a list of )
            node position 

        Todo
        ----

        update the position of the antenna in the ACS (for now only DCS is updated)

        """

        if (isinstance(p,np.ndarray)) or (isinstance(n,list) and isinstance(p,list) ):
            # Tranfrom input as list
            if not(isinstance(n,list)):
                n=[n]
                p=[p]
            if len(n) == len(p):    
                d=dict(zip(n,p))    # transform data to be complient with nx.set_node_attributes            
                nowd=dict(zip(n,[now]*len(n)))
            else :
                raise TypeError('n and p must have the same length')
            # update position
            nx.set_node_attributes(self,p_pe,d)        
            # update time of ground truth position
            if p_pe=='p':
                nx.set_node_attributes(self,'t',nowd)
        else :
            raise TypeError('n and p must be either: a key and a np.ndarray, or 2 lists')

    def update_dis(self):

        p = self.get_pos()
        e = self.edges()
        lp = np.array([np.array((p[e[i][0]],p[e[i][1]])) for i in range(len(e))])
        d = np.sqrt(np.sum((lp[:,0]-lp[:,1])**2,axis=1))
        [self.edge[ve[0]][ve[1]].update({'d':d[ie]}) for ie,ve in enumerate(self.edges())]



    def get_orient(self,wstd=None):
        """ get node orientations

        Parameters
        ----------

        wstd : specify a wstd to display node orientaion.
               If None, all wstd are displayed    

        Returns 
        -------

        dictionnary :     key     : node ID
        value     : np.array node position


        """
        if wstd == None:
            return nx.get_node_attributes(self,'T')
        else :
            try:
                return nx.get_node_attributes(self.SubNet[wstd],'T')
            except: 
                raise AttributeError('invalid wstd name')

    def get_pos(self,wstd=None):
        """ get node positions

        Parameters
        ----------

        wstd : specify a wstd to display node position. If None, all wstd are return

        Returns 
        -------

        dictionnary :     key     : node ID
        value     : np.array node position


        """
        if wstd == None:
            if self.node[self.nodes()[0]].has_key('p'):
                return nx.get_node_attributes(self,'p')
            else :
                return nx.get_node_attributes(self,'pe')
        else :
            try:
                if self.SubNet[wstd].node[self.SubNet[wstd].nodes()[0]].has_key('p'):
                    return nx.get_node_attributes(self.SubNet[wstd],'p')
                else :
                    return nx.get_node_attributes(self.SubNet[wstd],'pe')
            except: 
                raise AttributeError('invalid wstd name')

    def get_pos_est(self,wstd=None):
        """ get node estimated  positions ( only available in PN network)

        Parameters
        ----------

        wstd : specify a wstd to display node position. If None, all wstd are displayed    
        
        Returns 
        ------

        dictionnary :     key     : node ID
        value     : np.array node position
        

        """
        if wstd == None:
            return nx.get_node_attributes(self,'pe')
        else :
            try:
                return nx.get_node_attributes(self.SubNet[wstd],'pe')
            except: 
                raise AttributeError('invalid wstd name')



    def haspe(self,n):
        """
            Test if a node has an estimated point  pe key

        Parameters
        ----------

        n : int
            node numbner

        Returns
        -------

        Boolean : True if node n has a pe k

        """
        try:
            return  self.node[n]['pe'].any()
        except:
            return False


    def overview(self):
        """ overview of the network

        Returns
        -------

        O : dict


        """
        O={}
        for sn in self.SubNet.iteritems():
            for ldp in self.LDP:
                try:
                    O[sn[0]].update({ldp:nx.get_edge_attributes(sn[1],ldp)})
                except:
                    O[sn[0]]={ldp:nx.get_edge_attributes(sn[1],ldp)}

        return (O)


    def pp(self):
        """ pretty print information
            OBSOLETE

        Print information on edges connection and LDPs values and accuracy

        """
        for wstd in self.wstd.keys():
            print '-'*30
            print wstd
            print('{0:10} | {1:5} | {2:5} | {3:5} | {4:5} | {5:5} |'.format('Node link','TOA ','TOA std', 'Pr','Pr std', 'distance' ))
            print '-'*30
            T=nx.get_edge_attributes(self.SubNet[wstd],'TOA')
            P=nx.get_edge_attributes(self.SubNet[wstd],'Pr')
            D=nx.get_edge_attributes(self.SubNet[wstd],'d')
            for i in self.SubNet[wstd].edges(): # boucle sur toute les liaisons
                print('{0:10} | {1:1.4} | {2:7.4} | {3:1.4} | {4:7.4} | {5:7.4} |'.format(i,T[i][0],T[i][1],P[i][0],P[i][1],D[i]))




    def show(self, wstd=None, legend=False, ion=False, info=False, fig=plt.figure() ,ax=None, name=None):
        """ 
        Show the network

        Parameters 
        ----------

        wstd     : specify a wstd to display. If None, all wstd are displayed
        legend     : Bool. Toggle display edge legend
        ion     : interactive mode for matplotlib 
        info    : plot information on edges 
        fig     : plt.figure() to plot 
        ax      : plt.figure.ax to plot 
        name    : figure name


        """
        C = ConfigParser.ConfigParser()
        C.read(pyu.getlong('show.ini', 'ini'))
        color = ['r', 'g', 'b', 'm', 'y', 'c']*5
        style = ['-']*10

        wstdcolor = {k:color[uk] for uk, k in enumerate(self.SubNet.keys())}
        wstdes = {k:style[uk] for uk, k in enumerate(self.SubNet.keys())}

        # stdcolor = dict(C.items('wstdcolor'))
        # wstdes = dict(C.items('wstdestyle'))




        if wstd == None:
            rloop = self.wstd.keys()

        else :
            if isinstance(wstd,list):
                rloop = wstd    
            elif isinstance(wstd,str) :
                rloop=[wstd]    
            else :
                raise AttributeError('Arg must be a string or a string list')

        if fig==None:
            fig = plt.figure()
            ax=fig.add_subplot(111)
        elif ax== None:
            ax=fig.add_subplot(111)
        else:
            plt.figure(name)
            ax.axis('scaled')




        try:
            self.coll_plot['node'][1]=[]
            self.coll_plot['label'][1]=[]
            self.coll_plot['edge'][1]=[]
            Cl=[]
        except:
            self.coll_plot['node']=[[]]
            self.coll_plot['node'].append([])
            self.coll_plot['label']=[[]]
            self.coll_plot['label'].append([])
            self.coll_plot['edge']=[[]]
            self.coll_plot['edge'].append([])
            Cl=[]





        for ii,rl in enumerate(rloop):
            pos = self.get_pos(rl) 
            pos = {k:v[:2] for k,v in pos.items()}
            self.coll_plot['node'][1].append(nx.draw_networkx_nodes(
                                            self,
                                            pos=pos,
                                            nodelist=self.SubNet[rl].nodes(),
                                            node_size=100.,
                                            node_color='r',
                                            ax=ax))
            Cl=nx.draw_networkx_labels(self.SubNet[rl],
                                       pos=pos,
                                       font_size=10,
                                       ax=ax)
            self.coll_plot['label'][1].extend(Cl.values())
            self.coll_plot['edge'][1].append((nx.draw_networkx_edges(
                                              self,
                                              pos=pos,
                                              edgelist=self.SubNet[rl].edges(),
                                              arrows=False,
                                              width=2.,
                                              alpha=0.9,
                                              edge_color=wstdcolor[rl],
                                              style=wstdes[rl],
                                              ax=ax)))

        if legend:
            ax.legend((self.coll_plot['edge'][1]),(rloop),loc=3)
        if info :
            L=nx.get_edge_attributes(self,'TOA')


        if ion:
            try:
                [jj.remove() for jj in self.coll_plot['node'][0]]
                [jj.remove() for jj in self.coll_plot['edge'][0] if jj != None]
                [jj.remove() for jj in self.coll_plot['label'][0]]
            except:
                pass 
            plt.draw()
            self.coll_plot['node'][0]=self.coll_plot['node'][1]
            self.coll_plot['edge'][0]=self.coll_plot['edge'][1]
            self.coll_plot['label'][0]=self.coll_plot['label'][1]

        return fig, ax


    def _show3(self, wstd=None,newfig=False):
        """ Mayavi _show3

        Parameters
        ----------

        wstd : list
            list of wireless standards

        """

        color = ['r', 'g', 'b', 'm', 'y', 'c']*5
        wstdcolor = {k:color[uk] for uk, k in enumerate(self.SubNet.keys())}
        cold = pyu.coldict()

        if not newfig:
            f = mlab.gcf()

        if wstd == None:
            rloop = self.wstd.keys()

        else :
            if isinstance(wstd,list):
                rloop = wstd    
            elif isinstance(wstd,str) :
                rloop=[wstd]    
            else :
                raise AttributeError('Arg must be a string or a string list')

        for ii,rl in enumerate(rloop):
            
            pos = self.get_pos(rl)
            posv = pos.values()
            mp = dict(zip(pos.keys(),range(len(pos.keys()))))
            edg = self.SubNet[rl].edges()
            connect = [(mp[e[0]],mp[e[1]]) for e in edg]
            posv = np.array(posv)

            pts = mlab.points3d(posv[:,0], posv[:,1], posv[:,2],
                                   scale_factor=0.01, resolution=10)
            pts.mlab_source.dataset.lines = np.array(connect)
            tube = mlab.pipeline.tube(pts, tube_radius=0.01)
            colhex = cold[wstdcolor[rl]]
            col = tuple(pyu.rgb(colhex)/255.)
            mlab.pipeline.surface(tube, color=col)



    def csv_save(self,filename,S):
        """ save node positions into csv file

        Parameters 
        ----------

        filename : string 
                   name of the csv file
        S        : Simulation
                   Scipy.Simulation object

        """

        pos = np.array(nx.get_node_attributes(self,'p').values())
        pos = np.hstack((pos,np.zeros((len(self.nodes()),1))))  # passage en 3D
        pos = pos.reshape((1,len(self.nodes())*3))
        filecsv = pyu.getlong(filename,pstruc['DIRNETSAVE'])+'.csv'
        #file=open('../save_data/' +filename +'.csv','a')
        file = open(filecsv,'a')
        file.write(str(S.now()) +',')
        np.savetxt(file,pos,delimiter=',')
        file.close()


    def init_save(self,height=1.5):
        """
        Parameter
        ---------

        init_save

        """

        pos=nx.get_node_attributes(self,'p').items()

        AP=[]
        AG=[]
        api=1
        loc=False
        method = []
        # get methods for localization
        simcfg = ConfigParser.ConfigParser()
        simcfg.read(pyu.getlong('simulnet.ini','ini'))
        save =eval(simcfg.get('Save','save'))
        if 'loc' in save:
            loc = True
            method = eval(simcfg.get('Localization','method'))

        ## find Agent and Acces point
        for i in range(len(pos)):
            if self.node[pos[i][0]]['typ'] =='ap':
                AP.append(pos[i][0])
                if not os.path.isfile(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE'])):
                    file=open(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE']),'w')
                    config = ConfigParser.ConfigParser()
                    config.add_section('coordinates')
#                    config.set('coordinates',str(api), str(pos[i][1][0]) + ' ' + str(pos[i][1][1]) + ' '+str(height))
                    config.set('coordinates','1', str(pos[i][1][0]) + ' ' + str(pos[i][1][1]) + ' '+str(height))
                    api=api+1
                    config.write(file)
                    file.close()
            else:
                AG.append(pos[i][0])
                config = ConfigParser.ConfigParser()
                if not os.path.isfile(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE'])):
                    file=open(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE']),'w')
                    config.add_section('coordinates')
                    if loc :
                        if 'geo' in method:
                            config.add_section('geo_est')
                        if 'alg' in method:
                            config.add_section('alg_est')
                # if simulation has already been runed with localization, this
                # ensure that localization section will be created

                else :
                    file=open(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE']),'w')
                    config.read(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE']))
                    if 'coordinates' not in config.sections():
                        config.add_section('coordinates')
                    if 'geo_est' not in config.sections() and 'geo' in method:
                        config.add_section('geo_est')
                    if 'alg_est' not in config.sections() and 'alg' in method:
                        config.add_section('alg_est')

                config.write(file)
                file.close()



        if 'pyray' in save :

            file2=open(pyu.getlong('pyray.ini',pstruc['DIRNETSAVE']),'w')
            config = ConfigParser.ConfigParser()
            config.add_section('nodes')
            config.add_section('layout')
            config.add_section('simulation')
            config.set('nodes','AG',str(AG))
            config.set('nodes','AP',str(AP))
            config.set('simulation','updatetime',str(simcfg.get('Network','network_update_time')))
            config.set('layout','layoutname',str(simcfg.get('Layout','filename')))
            config.write(file2)
            file2.close()

        if 'loc' in save :

            file2=open(pyu.getlong('loc.ini',pstruc['DIRNETSAVE']),'w')
            config = ConfigParser.ConfigParser()
            config.add_section('nodes')
            config.add_section('simulation')
            config.set('nodes','AG',str(AG))
            config.set('nodes','AP',str(AP))
            config.set('simulation','loc_updatetime',str(simcfg.get('Localization','localization_update_time')))
            config.set('simulation','method',str(simcfg.get('Localization','method')))
            config.set('simulation','duration',str(simcfg.get('Simulation','duration')))
            config.write(file2)
            file2.close()

        return method





    def mat_save(self,S):
        """
        DEPRECATED 

        REPLACED BY pylayers.util.save

        DEPRECATED 

            save node positions into a matlab structure file


        Parameters
        ----------

        filename : string 
                   name of the mat file
        S        : Simulation
                   Scipy.Simulation object
        """

        pos=nx.get_node_attributes(self,'p').items()
        for i in range(len(pos)):
            if not 'BS' in pos[i][0]:
                try:
                    self.mat[pos[i][0]]['pos']=np.vstack((self.mat[pos[i][0]]['pos'],pos[i][1]))
                    self.mat[pos[i][0]]['time']=np.vstack((self.mat[pos[i][0]]['time'],S.now()))
                except:
                    self.mat[pos[i][0]]={}
                    self.mat[pos[i][0]]['pos']=pos[i][1]
                    self.mat[pos[i][0]]['time']=np.array(S.now())
            else :
                try:
                    self.mat[pos[i][0]]['pos']=pos[i][1]
                except:
                    self.mat[pos[i][0]]={}
                    self.mat[pos[i][0]]['pos']=pos[i][1]


        sp.io.savemat(pyu.getlong('mat.mat','save_data'),self.mat)

#    def sql_save(self,S):
#        """
#        save network state into mysqldatabase


#        Attributes:
#        ----------
#        
#        S        : Simulation
#                   Scipy.Simulation object

#        """
#        self.db.writenet(self,S.now())

    def txt_save(self,S):
        """
        DEPRECATED 

        REPLACED BY pylayers.util.save

        DEPRECATED 



        save network state into mysqldatabase


        Parameters 
        ----------

        S        : Simulation
                   Scipy.Simulation object

        """
        pyu.writenet(self,S)


    def loc_save(self,S):
        """
        DEPRECATED 

        REPLACED BY pylayers.util.save

        DEPRECATED 


        save txt 
        node ID , True pos x , True pos y , est pos x , est pos y , timestamp


        Parameters
        ----------

        S        : Simulation
                   Scipy.Simulation object

        """

        pos=nx.get_node_attributes(self,'p')
        pe=nx.get_node_attributes(self,'pe_alg')
        typ = nx.get_node_attributes(self,'typ')
        if self.idx == 0:
            entete = 'NodeID, True Position x, True Position y, Est Position x, Est Position y, Timestamp\n'
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/simulation.txt','write')
            file.write(entete)
            file.close()

        try:
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/simulation.txt','a')
            for n in self.nodes():
                if typ[n] != 'ap':
                    data = n + ',' + str(pos[n][0]) + ',' + str(pos[n][1]) + ',' + str(pe[n][0][0]) + ',' + str(pe[n][0][1]) + ',' +pyu.timestamp(S.now()) +',\n'
                    file.write(data)
            file.close()
            self.idx = self.idx +1
        except:
            pass


    # def dual_save(self,S):
    #     """
    #     DEPRECATED 

    #     REPLACED BY pylayers.util.save

    #     DEPRECATED 



    #     save txt 



    #     Parameters
    #     ----------

    #     S        : Simulation
    #                Scipy.Simulation object

    #     """

    #     pos=nx.get_node_attributes(self,'p')
    #     pclust = nx.get_node_attributes(self,'pe_clust')
    #     typ = nx.get_node_attributes(self,'typ')
    #     if self.idx == 0:
    #         entete = 'Timestamp, True Position x, True Position y, Est Position1 x, Est Position1 y,Est Position2 x, Est Position2 y\n'
    #         file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/pos.txt','write')
    #         file.write(entete)
    #         file.close()
    #         file2=open(basename+'/' + pstruc['DIRNETSAVE'] +'/rsslink.txt','write')
    #         entete2 = 'Timestamp, link, linkid, Pr, distance\n'
    #         file2.write(entete2)
    #         file2.close()
    #         file3=open(basename+'/' + pstruc['DIRNETSAVE'] +'/anchorposition.txt','write')
    #         data3 = 'node,pos x, pos y\n'
    #         file3.write(data3)
    #         for n in self.nodes():
    #             data3= n + ',' + str(self.node[n]['p'][0]) + ',' + str(self.node[n]['p'][1]) + '\n'
    #             file3.write(data3)
    #         file3.close()
    #         file4=open(basename+'/' + pstruc['DIRNETSAVE'] +'/toa.txt','w')
    #         entete4 = 'Timestamp, typ, toaid, toa,distance\n'
    #         file4.write(entete4)
    #         file4.close()

    #     try:

    #         file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/pos.txt','a')
    #         file2=open(basename+'/' + pstruc['DIRNETSAVE'] +'/rsslink.txt','a')
    #         file4=open(basename+'/' + pstruc['DIRNETSAVE'] +'/toa.txt','a')
    #         for n in self.nodes():
    #             if n == '1':
    #                 data =  pyu.timestamp(S.now()) +','+ str(pos[n][0]) + ',' + str(pos[n][1]) + ',' + str(pclust[n][0,0]) + ',' + str(pclust[n][0,1]) + ',' + str(pclust[n][1,0]) + ',' + str(pclust[n][1,1]) +'\n'
    #                 for e in self.edge[n].keys():
    #                     if e != '6' and e !='7':
    #                         try:
    #                             data2 = data2 +',link,' + str(e) + ',' + str(self.edge[n][e]['rat1']['Pr'][0]) +',' + str(np.sqrt(np.sum((pos[n]-pos[e])**2)))
    #                         except:
    #                             data2 = pyu.timestamp(S.now()) + ',link,' + str(e) + ',' + str(self.edge[n][e]['rat1']['Pr'][0]) +',' + str(np.sqrt(np.sum((pos[n]-pos[e])**2)))
    #                     else :
    #                         try:
    #                             data4 = data4 +',toa,' + str(e) + ',' + str(self.edge[n][e]['rat1']['TOA'][0]) +',' + str(np.sqrt(np.sum((pos[n]-pos[e])**2)))
    #                         except:

    #                             data4 = pyu.timestamp(S.now()) + ',toa,' + str(e) + ',' + str(self.edge[n][e]['rat1']['TOA'][0]) +',' +str(np.sqrt(np.sum((pos[n]-pos[e])**2)))

    #         data2=data2 + '\n'
    #         data4=data4 + '\n'
    #         file.write(data)
    #         file2.write(data2)
    #         file4.write(data4)


    #         file.close()
    #         file2.close()
    #         file4.close()
    #         self.idx = self.idx +1
    #     except:
    #         pass


    def pyray_save(self,S):
        """
            save node positions into ini file, compliant with pyray standard

        Parameters
        ----------

        filename : string 
                   name of the pyray file
        S        : Simulation
                   Scipy.Simulation object
        """


        assert len(self.SubNet.keys()) == 1 , NameError('when network.ini_save() \
        is used , only 1 wstd must be involved in the Network.\
        Please modify agent.ini')


        height= 1.5
        pos=nx.get_node_attributes(self,'p').items()

        ### create ini files
        if self.idx == 0:
            self.init_save(height=height)
        ### save agent positions
        for i in range(len(pos)):
            if self.node[pos[i][0]]['typ'] !='ap':
                config = ConfigParser.ConfigParser()
                config.read(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE']))
                config.set('coordinates',str(self.idx+1),value = str(pos[i][1][0]) + ' ' + str(pos[i][1][1]) + ' '+str(height))
                file=open(pyu.getlong(str(pos[i][0]) + '.ini',pstruc['DIRNETSAVE']),'w')
                config.write(file)
                file.close()


    def loc_save(self,S,node='all',p=False):
        """
            save node estimated positions into ini file,

        Attributes:
        ----------

        S        : Simulation
                   Scipy.Simulation object
        """

        if node  == 'all':
            node = self.nodes()
        elif not isinstance(node,list):
            node = [node]


        height=1.5
        ### create ini files
        if self.lidx == 0:
            self.init_save(height=height)

        pe_alg = nx.get_node_attributes(self,'pe_alg')
        pe_geo = nx.get_node_attributes(self,'pe_geo')
        p = nx.get_node_attributes(self,'p')
        ### save agent positions estimations
        for n in node:
            if self.node[n]['typ'] !='ap':
                config = ConfigParser.ConfigParser()
                config.read(pyu.getlong(str(n[0]) + '.ini',pstruc['DIRNETSAVE']))
                if pe_alg != {} :
                    config.set('alg_est',str(self.idx+1),value = str(pe_alg[n[0]][0]) + ' ' + str(pe_alg[n[0]][1]) + ' '+str(height))
                if pe_geo != {} :
                    config.set('geo_est',str(self.idx+1),value = str(pe_geo[n[0]][0]) + ' ' + str(pe_geo[n[0]][1]) + ' '+str(height))
                if p:
                    config.set('coordinates',str(self.idx+1),value = str(p[n[0]][0]) + ' ' + str(p[n[0]][1]) + ' '+str(height))
                file=open(pyu.getlong(str(n[0]) + '.ini',pstruc['DIRNETSAVE']),'w')
                config.write(file)
                file.close()
        self.lidx=self.lidx+1



    def ini_save(self,S,filename='simulnet_data.ini',height=1.5):
        """
        ----------
        DEPRECATED
        ----------
        Save an .ini file of node position . 
        Only links  which involve mobile nodes (typ 'ag') are kept.

        The produced init file is filled as follow:

            [timestamp]
            nodeID1_nodeID2 = x1,y1,z1,x2,y2,z2
            nodeID2_nodeID4 = x2,y2,z2,x4,y4,z4
            ....


        Attributes:
        ----------

        S        : Simulation
                   Scipy.Simulation object

        filename  : string
                   name of the saved ini file

        height    : float
                   height of the nodes





        """

        assert len(self.SubNet.keys()) == 1 , NameError('when network.ini_save() \
        is used , only 1 wstd must be involved in the Network.\
        Please modify agent.ini')


        if self.idx == 0:
            file=open(pyu.getlong(filename ,'output'),'w')
        else:
            file=open(pyu.getlong(filename ,'output'),'a')

        config = ConfigParser.ConfigParser()
        timestamp = pyu.timestamp(S.now())
        config.add_section(timestamp)
        for e in self.edges():
            if not ((self.node[e[0][0]]['typ'] == 'ap') and  (self.node[e[1][0]]['typ'] == 'ap')):
                key=str(e[0]) +'_' +str(e[1])
                value1 = str(self.node[e[0][0]]['p'][0])+ ',' +str(self.node[e[0][0]]['p'][1])+','+str(height)
                value2 = str(self.node[e[1][0]]['p'][0])+ ',' +str(self.node[e[1][0]]['p'][1])+','+str(height)
                config.set(timestamp, key, value1 + ' , ' + value2)

        config.write(file)
        file.close()

        self.idx=self.idx+1



#class PN(nx.MultiDiGraph):

#    def __init__(self,N):
#        nx.MultiDiGraph.__init__(self)
#        self.add_nodes_from(N)
#        pdb.set_trace()
#        self.add_edges_from( (u,v,key,deepcopy(datadict))
#                           for u,nbrs in self.adjacency_iter()
#                           for v,keydict in nbrs.items()
#                           for key,datadict in keydict.items() ) 
#        pdb.set_trace()
#        self.node=N.node












class PNetwork(Process):
    """
    Process version of the Network class
    """
    def __init__(self,**args):
        defaults={'net':Network(),
                  'L':[],
                  'net_updt_time':0.001,
                  'sim':None,
                  'show_sg':False,
                  'disp_inf':False,
                  'save':[]}

##       initialize attributes
        for key, value in defaults.items():
            if args.has_key(key):
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key]=value  
        self.args=args

        Process.__init__(self,name='PNetwork',sim=self.sim)
        self.cpt=self.sim.now()
        self.filename='pos'

        if 'mysql' in self.save:
           config = ConfigParser.ConfigParser()
           config.read(pyu.getlong('simulnet.ini','ini'))
           sql_opt = dict(config.items('Mysql'))
           self.net.db = Database(sql_opt['host'],sql_opt['user'],sql_opt['passwd'],sql_opt['dbname'])




    def run(self):


        ####################################################################################
        # first iteration requested to correctely initiatilzing Personnal Networks's Subnets 
        for wstd in self.net.wstd.iterkeys():
            self.net.compute_LDPs(wstd)
        for n in self.net.nodes():
            self.net.node[n]['PN']._get_wstd()
            self.net.node[n]['PN']._get_SubNet()
            # Add access point position in each personal network (PN)
            [self.net.node[n]['PN'].node[n2].update({'pe':self.net.node[n2]['p']}) for n2 in self.net.node[n]['PN'].node.iterkeys() if self.net.node[n]['PN'].node[n2]['typ'] == 'ap']
                
        ####################################################################################
        self.pos=self.net.get_pos()



        if 'csv' in self.save:
            nbnodes = len(self.net.nodes())
            entete = 'time'
            inode=self.net.nodes_iter()
            for i in inode:
                entete = entete +',x'+str(i) +',y'+str(i)+',z'+str(i)
            entete=entete +'\n'
            filecsv = pyu.getlong(self.filename,pstruc['DIRNETSAVE'])+'.csv'
            #file=open('../save_data/' +self.filename +'.csv','w')
            file = open(filecsv,'w')
            file.write(entete)
            file.close()



        while True:

            ############### compute LDP
            for wstd in self.net.wstd.iterkeys():
                self.net.compute_LDPs(wstd)
            
            if self.show_sg:
                ############### compute Signature (Sg)
                tx=self.net.node.keys()[0]
                rx=self.net.node.keys()[1]
                Sg=self.net.compute_Sg(tx,rx)

            ############## Show

            if self.show_sg:
                self.net.show_sig(Sg,tx,rx,ion=True,fig=fig,ax=ax)

            if self.disp_inf:
                self.net.pp()




#            ############# save network
#            REPLACED BY A SAVE PROCESS
            if 'csv' in self.save:
                self.net.csv_save(self.filename,self.sim)
#            if 'pyray' in self.save:
#                self.net.pyray_save(self.sim)
#            if 'matlab' in self.save:
#                self.net.mat_save(self.sim)
#            if 'msql' in self.save:
#                self.net.sql_save(self.sim)
#            if 'txt' in self.save:
#                self.net.txt_save(self.sim)

#            if 'ini' in self.save:
#                self.net.ini_save(self.sim)
#            if 'loc' in self.save:
#                self.net.loc_save(self.sim)
#            if 'dual' in self.save:
#                self.net.dual_save(self.sim)


            self.net.pos=self.net.get_pos()
            if self.sim.verbose:
                print 'network updated @',self.sim.now()
            self.net.idx=self.net.idx+1
            yield hold, self, self.net_updt_time



