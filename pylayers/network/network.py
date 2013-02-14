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
import scipy as sp 
import networkx as nx
import itertools
import pickle as pk
import pkgutil
from copy import deepcopy


import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import Tkinter as tk

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


from SimPy.SimulationRT import Process,hold
import pprint
import select
import sys


# How to take into account  1 specific key specifique for 1 MultiGraph
# MULTIGRAPH !!!! G.add_edge(10,11,key='wifi',attr_dict=dict(Pr=0,TOA=10))



        

class Node(nx.MultiGraph):
    """ Class Node
    inherit of networkx.Graph()
    
    Attributes
    ----------
        Id    : float/hex/str/...
                node Id
        p    : np.array
                True position
        t    : time.time()
                Tag time
        RAT    : list
                available RAT of the node
        PN    : Network.Network
                Personal vision of the Network 
        pos    : Dictionnary
                parser from Node.Node to networkx.node.pos

    Method
    ------
        RandomMac(): Generate a RAndom Mac adress    

    """
    def __init__(self,ID=0,p=np.array(()),t=time.time(),pe=np.array(()),te=time.time(),RAT=[],epwr={},sens={},type='ag'):
        nx.MultiGraph.__init__(self)

        # Personnal Network init
        self.ID=ID
        self.PN = Network(owner=self.ID,PN=True)
        self.PN.add_node(self.ID,dict(pe=pe,te=te,RAT=RAT,type=type))
        # Network init

        self.add_node(ID,dict(PN=self.PN,p=p,pe=self.PN.node[self.ID]['pe'],t=t,RAT=RAT,epwr=epwr,sens=sens,type=type))
        self.p    = self.node[self.ID]['p']
        self.pe    = self.PN.node[self.ID]['pe']
        self.t    = self.node[self.ID]['t']
        self.RAT = self.node[self.ID]['RAT']
        self.epwr = self.node[self.ID]['epwr']
        self.sens = self.node[self.ID]['sens']





    def randomMAC(self):
        mac = [ 0x00, 0x16, 0x3e,
        random.randint(0x00, 0x7f),
        random.randint(0x00, 0xff),
        random.randint(0x00, 0xff) ]
        return ':'.join(map(lambda x: "%02x" % x, mac))	




class Network(nx.MultiGraph):
    """Network class
    inherit of networkx.Graph()

    Attributes
    ----------
    RAT : dictionnary
        keys  = RAT 
        value = list of nodes id
    RATe : dictionnary
        keys  = RAT 
        value = list of edges id  
    SubNet : dictionnary
        keys  = RAT 
        value = Subgraph of the given RAT
    pos : dictionnary
        keys  = node id
        value = node position

    Methods
    -------
    get_RAT(self)                    : Get RAT from nodes of the network
    connect(self)                    : Connect each node from a rat together
    create(self)                    : compute get_RAT(),get_pos() and connect()
    update_LDP(self,n1,n2,RAT,LDP=None,value=[])    : update Location Dependent Parameter  
    compute_LDP(self,n1,n2,RAT,LDP,method='direct') : compute the LDP value thanks to a ElectroMag Solver 
    update_pos(self,n,p=np.array)            : update node (or node list) position
    get_pos(self,RAT=None)                : get node positions
    pp(self)                    : pretty print on std out all edtges informations
    show(Rat=None,legend=True)            : Display network for all rat or specified in Rat. 
    """



    def __init__(self,owner='sim',EMS=EMSolver(),PN=False):
        nx.MultiGraph.__init__(self)
        self.owner=owner
        self.RAT={}
        self.LDP = ['TOA','Pr'] 
        self.SubNet={}
        self.EMS=EMS
        self.coll_plot={}
        self.pos={}
        self.mat={}
        self.idx = 0
        self.lidx = 0

    def perm(self,iterable,r,key,d=dict()):
        """ combi = itertools.permutation(iterable,r) adapted 
    
        This is an adapted version of itertools.permutations in order to be complient with the networkx.add_edges_from method.
        itertools.permutations(range(4), 3) --> 012 013 021 023 031 302 102 103 ...
    
        self.perm([10,11],2,'wifi') -->     (10, 11, 'wifi', {'Pr': [], 'TOA': []}) 
                                (11, 10, 'wifi', {'Pr': [], 'TOA': []})

        Attributes
        ----------
        iterable : list
            list of node
        r     : int
            number of node gathered in the output tuple ( always set 2 ! )


        Returns
        ------     

        out : tuple(node_list,r,RAT,d):
        node_list    : list of node1        
        r        : gather r node in the tuple
        RAT        : the specified RAT
        d        : dictionnary of RAT attribute

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
        (0, 2, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})        for l in self.LDP:
        >>> perm.next()
        (1, 0, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (1, 2, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (2, 0, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        >>> perm.next()
        (2, 1, 'toto', {'Pr': [], 'TOA': [], 'key1': 1, 'key2': 2})
        """

        for l in self.LDP:
            d[l]=[]
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
        """ combi = itertools.combination(iterable,r) adapted 
    
        This is an adapted version of itertools.combinations in order to be complient with the networkx.add_edges_from method.
        itertools.combinations('ABCD', 2) --> AB AC AD BC BD CD
        itertools.combinations(range(4), 3) --> 012 013 023 123
    
        self.combi([10,11,12],2,'wifi') -->     (10, 11, 'wifi', {'Pr': [], 'TOA': []}) 
                                (10, 12, 'wifi', {'Pr': [], 'TOA': []})
                            (11, 12, 'wifi', {'Pr': [], 'TOA': []})

        Attributes
        ----------
        iterable : list
            list of node
        r     : int
            number of node gathered in the output tuple ( always set 2 ! )


        Returns
        ------     

        out : tuple(node_list,r,RAT,d):
        node_list    : list of node1        
        r        : gather r node in the tuple
        RAT        : the specified RAT
        d        : dictionnary of RAT attribute

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


    def Gen_tuple(self,gene,rat,var):
        """
        generate a specific tuple  

        Attributes
        ----------

        gene : tuple(x,y) iterator 
        rat  : str
        var  : list
            len(var) = len(gene)

        Yield
        -----
        tuple : (gene[i][0],gene[i][1],rat,var[i]) for iteration i

        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>> tup = zip(range(5),range(5))
        [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]
        >>> g = iter(tup)
        >>> rat='string rat'
        >>> var=[10,11,12,13,14]

        >>> T=N.Gen_tuple(g,rat,var)
        >>> T.next()
        (0, 0, 'string rat', 10)
        >>> T.next()
        (1, 1, 'string rat', 11)

            """

        gvar=iter(var)
        while True:
            G=gene.next()
            Gvar=gvar.next()
            yield(tuple((G[0],G[1],rat,Gvar)))


    def get_RAT(self):
        """ get rat from nodes of the network
        

        Attributes
        ----------
            Rat : specify which RAt you want to append to your network. If None, all rat are appended.

        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>>    N=Network.Network()
        >>> for i in range(3):
                no = Node.Node(ID=i,RAT=['wifi','bt'])
                N.add_nodes_from(no.nodes(data=True))

        >>> N.get_RAT()
        {'bt': [0, 1, 2], 'wifi': [0, 1, 2]}

        """
#        if Rat !='None' :

        for no in self.nodes():
            for r in self.node[no]['RAT']:
                try:
                    self.RAT[r].append(no)
                except :
                    self.RAT[r]=[no]
#        else :
#            for no in self.nodes(): 
#                for r in self.node[no]['RAT']:
#                    try:
#                        self.RAT[r].extend(no)
#                    except :
#                        self.RAT[r]=[no]

        # uniquify results
        for Rat in self.RAT.keys():
            self.RAT[Rat]    = {}.fromkeys(self.RAT[Rat]).keys()
        


    def connect(self):
        """ 
        Connect all nodes from the network sharing the same RAT and creating the associated Subnetwork

        """
        
        edge_dict={}
        for l in self.LDP:
            edge_dict[l]=[[]]
        edge_dict['vis']=False
        


        for ratnb,Rat in enumerate(self.RAT.keys()):
            edges=self.combi(self.RAT[Rat],2,Rat,d=edge_dict)
            self.add_edges_from(edges)    
            self.get_SubNet(Rat)

    def get_SubNet(self,Rat=None):
        """
        get SubNetworks of a network
        !!! ALWAYS use self.get_RAT() BEFORE !!!!!



        Attributes
        ----------
        Rat : specify which SubNet you want to create
        
        Examples
        --------

        >>> from pylayers.network.network import *
        >>> N=Network()
        >>> for i in range(2):
                no = Node.Node(ID=i,RAT=['wifi','bt'])
                N.add_nodes_from(no.nodes(data=True))

        >>> no = Node.Node(ID=2,RAT=['wifi'])
        >>>    N.add_nodes_from(no.nodes(data=True))
        >>> N.get_RAT() # VERY IMPORTANT 

        >>> N.get_SubNet()
        >>> N.SubNet['bt'].nodes()
        [0, 1]
        >>> N.SubNet['wifi'].nodes()
        [0, 1, 2]


        """

        if Rat == None:
        #    pdb.set_trace()
            for Rat in self.RAT:            
                # creating all SubNetworks 
                self.SubNet[Rat]= self.subgraph(self.RAT[Rat])
                # remove information from previous subnetwork (because subgraph copy the whole edge information)
                ek = self.SubNet[Rat].edges(keys=True)
                for e in ek :
                    if e[2] != Rat:
                        self.SubNet[Rat].remove_edge(e[0],e[1],e[2])
                for n in self.SubNet[Rat].nodes():
                    try:
                        self.SubNet[Rat].node[n]['epwr']=self.SubNet[Rat].node[n]['epwr']
                        self.SubNet[Rat].node[n]['sens']=self.SubNet[Rat].node[n]['sens']
                    except: 
                        pass


        elif Rat in self.RAT:
            # creating SubNetworks
            self.SubNet[Rat]= self.subgraph(self.RAT[Rat])
            # remove information from previous subnetwork (because subgraph copy the whole edge information)
            for k in self.RAT.keys():
                if k != Rat:
                    try:
                        self.SubNet[Rat].remove_edges_from(self.SubNet[k].edges(keys=True))
                    except :
                        pass
                for n in self.SubNet[Rat].nodes():
                    try:
                        self.SubNet[Rat].node[n]['epwr']=self.SubNet[Rat].node[n]['epwr']
                        self.SubNet[Rat].node[n]['sens']=self.SubNet[Rat].node[n]['sens']
                    except: 
                        pass


        else :
            raise NameError('invalid RAT name')


    def init_PN(self):
        """ 
        Initializing personnal networks

        """

#        for Sn in self.SubNet.iteritems():
#            print Sn
#            for n in Sn[1].nodes():
#                print Sn
#                pdb.set_trace()
#                try:
#                    [Sn[1].node[n]['PN'].node[nn]['RAT'].append(Sn[0]) for nn in Sn[1].nodes() if nn != n]
#                except:
#                    [Sn[1].node[n]['PN'].add_node(nn,attr_dict=dict(RAT=[Sn[0]],pe=np.array(()),te=time.time()),type=Sn[1].node[nn]['type']) for nn in Sn[1].nodes() if nn != n] 

        for Sn in self.SubNet.iteritems():
            for n in Sn[1].nodes():
                for nn in Sn[1].nodes():
                    if nn != n:
                        try:
                            Sn[1].node[n]['PN'].node[nn]['RAT'].append(Sn[0])
                        except:
                            Sn[1].node[n]['PN'].add_node(nn,attr_dict=dict(RAT=[Sn[0]],pe=np.array(()),te=time.time()),type=Sn[1].node[nn]['type'])
#                pdb.set_trace()

                ### init edge of PN

                Z= Sn[1].edges(n,keys=True,data=True)
                Sn[1].node[n]['PN'].add_edges_from(Z)
#                if Sn[1].node[n]['PN'].edge[n] == {}:
#                    Z1=self.swap_lt(Z,n)
#                    Sn[1].node[n]['PN'].add_edges_from(Z1)
#                print n,Sn[1].node[n]['PN'].edge[n]

    def create(self):
        """ create the network

        This method computes :     
            * get_RAT()
            * connect()
            * connect_PN()

        
        
        """
        self.get_RAT()
        self.connect()
        self.init_PN()

    def update_PN(self):
        """ update personnal network



        """
        ####################################################################################
        # first iteration requested to correctely initiatilzing Personnal Networks's Subnets 
        for rat in self.RAT.iterkeys():
            for ldp in self.LDP:
                self.compute_LDPs(self.nodes(),rat,ldp,method='direct')
        for n in self.nodes():
            self.node[n]['PN'].get_RAT()
            self.node[n]['PN'].get_SubNet()
            # Add access point position in each personal network (PN)
            [self.node[n]['PN'].node[n2].update({'pe':self.node[n2]['p']}) for n2 in self.node[n]['PN'].node.iterkeys() if self.node[n]['PN'].node[n2]['type'] == 'ap']

        ####################################################################################


#    def visibility(func):
#        def wrapper(*args, **kwargs):
#            a = list(args)
#            pdb.set_trace()
#            print 'decorator',a
#            return func(*args, **kwargs)
#        return wrapper
    def dist_edge(self,e,dp):
        return(np.array([np.sqrt(np.sum((dp[i[0]]-dp[i[1]])**2)) for i in e]))


    def update_LDPs(self,ln,RAT,lD):
        """Set a value between 2 nodes (n1 and n2) for a specific LDP from a RAT
            
        This method update :     * The network edges 
                    * The personal network (PN) of both n1 and n2

        Attributes
        ----------

        n1      : node ID
        n2      : node ID
        RAT     : string
            A specific RAT which exist in the network ( if not , raises an error)
        ln     : list 
            list of nodes
        lD     : list of dictionnary:
            [ {LDP1_1:[value , std],LDP2_1:[value , std] } , {LDPL_N:[value , std],LDPL_N:[value , std] }    ] for N nodes and L LDPS

        TODO 
        ----
        Check if LDP value is complient with the LDP
                 

        """
        
        self.SubNet[RAT].add_edges_from(self.Gen_tuple(self.SubNet[RAT].edges_iter(),RAT,lD))


    def compute_LDPs(self,ln,RAT):
        """compute edge LDP

        Attributes
        ----------

        n1      : float/string
            node ID
        n2      : float/string
            node ID
        RAT     : string
            A specific RAT which exist in the network ( if not , raises an error)
        value    : list : [LDP value , LDP standard deviation] 
        method    : ElectroMagnetic Solver method ( 'direct', 'Multiwall', 'PyRay'


        """

        p=nx.get_node_attributes(self.SubNet[RAT],'p')
        epwr=nx.get_node_attributes(self.SubNet[RAT],'epwr').values()

        e=self.SubNet[RAT].edges()
        lp,lt, d= self.EMS.solve(p,e,'all',RAT,epwr)
        lD=[{'Pr':lp[i],'TOA':lt[i] ,'d':d[i]} for i in range(len(d))]
        self.update_LDPs(ln,RAT,lD)



    def update_pos(self,n,p,p_pe='p'):
        """ 
        Update Position of a node

        Attributes
        ----------

        n      : float/string (or a list of)
            node ID
        p    : np.array  ( or a list of )
            node position 
                
        """

        if (isinstance(p,np.ndarray)) or (isinstance(n,list) and isinstance(p,list) ):
            # Tranfrom input as list
            if not(isinstance(n,list)):
                n=[n]
                p=[p]
            if len(n) == len(p):    
                d=dict(zip(n,p))    # transform data to be complient with nx.set_node_attributes            
            else :
                raise TypeError('n and p must have the same length')
            # update position
            nx.set_node_attributes(self,p_pe,d)        

        else :
            raise TypeError('n and p must be either: a key and a np.ndarray, or 2 lists')


    def get_pos(self,RAT=None):
        """ get node positions

        Attributes
        ----------
        RAT : specify a RAT to display node position. If None, all RAT are displayed    
        
        Returns 
        ------
        dictionnary :     key     : node ID
        value     : np.array node position
        

        """
        if RAT == None:
            if self.node[self.nodes()[0]].has_key('p'):
                return nx.get_node_attributes(self,'p')
            else :
                return nx.get_node_attributes(self,'pe')
        else :
            try:
                if self.SubNet[RAT].node[self.SubNet[RAT].nodes()[0]].has_key('p'):
                    return nx.get_node_attributes(self.SubNet[RAT],'p')
                else :
                    return nx.get_node_attributes(self.SubNet[RAT],'pe')
            except: 
                raise NameError('invalid RAT name')

    def get_pos_est(self,RAT=None):
        """ get node estimated  positions ( only available in PN network)

        Attributes
        ----------
        RAT : specify a RAT to display node position. If None, all RAT are displayed    
        
        Returns 
        ------
        dictionnary :     key     : node ID
        value     : np.array node position
        

        """
        if RAT == None:
            return nx.get_node_attributes(self,'pe')
        else :
            try:
                return nx.get_node_attributes(self.SubNet[RAT],'pe')
            except: 
                raise NameError('invalid RAT name')



    def haspe(self,n):
        """
            Test if a node has a pe key
        
        Returns
        -------
            Boolean : True if node n has a pe k
        """
        try:
            return  self.node[n]['pe'].any()
        except:
            return False


    def overview(self):
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
        for rat in self.RAT.keys():
            print '-'*30
            print rat
            print('{0:10} | {1:5} | {2:5} | {3:5} | {4:5} | {5:5} |'.format('Node link','TOA ','TOA std', 'Pr','Pr std', 'distance' ))
            print '-'*30
            T=nx.get_edge_attributes(self.SubNet[rat],'TOA')
            P=nx.get_edge_attributes(self.SubNet[rat],'Pr')
            D=nx.get_edge_attributes(self.SubNet[rat],'d')
            for i in self.SubNet[rat].edges(): # boucle sur toute les liaisons
                print('{0:10} | {1:1.4} | {2:7.4} | {3:1.4} | {4:7.4} | {5:7.4} |'.format(i,T[i][0],T[i][1],P[i][0],P[i][1],D[i]))




    def show(self,RAT=None,legend=False,ion=False,info=False,fig=plt.figure(),ax=None,name=None):
        """ 
        Show the network

        Attributes 
        ----------

        RAT     : specify a RAT to display. If None, all RAT are displayed
        legend     : Bool. Toggle display edge legend
        ion     : interactive mode for matplotlib 
        info    : plot information on edges 
        fig     : plt.figure() to plot 
        ax      : plt.figure.ax to plot 
        name    : figure name

            
        """
        C=ConfigParser.ConfigParser()
        C.read(pyu.getlong('show.ini','ini'))
        RATcolor=dict(C.items('RATcolor'))
        RATes    =dict(C.items('RATestyle'))

        
        

        if RAT == None:
            rloop = self.RAT.keys()

        else :
            if isinstance(RAT,list):
                rloop = RAT    
            elif isinstance(RAT,str) :
                rloop=[RAT]    
            else :
                raise NameError('Arg must be a string or a string list')

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
            self.coll_plot['node'][1].append(nx.draw_networkx_nodes(self,pos=pos,nodelist=self.SubNet[rl].nodes(),node_size=100.,node_color='r',ax=ax))
            Cl=nx.draw_networkx_labels(self.SubNet[rl],pos=pos,font_size=10,ax=ax)
            self.coll_plot['label'][1].extend(Cl.values())
            self.coll_plot['edge'][1].append((nx.draw_networkx_edges(self,pos=pos,edgelist=self.SubNet[rl].edges(),width=2.,alpha=0.9,edge_color=RATcolor[rl],style=RATes[rl],ax=ax)))
            
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

    def compute_Sg(self,tx,rx):
        if self.pos == {}:
            self.pos=self.get_pos()
        return (self.EMS.L.signature(self.pos[tx],self.pos[rx]))


    def show_sig(self,Sg,tx,rx,ion=False,fig=None,ax=None):
        if fig==None:
            fig = plt.figure()
            ax=fig.add_subplot(111)
        elif ax== None:
            ax=fig.add_subplot(111)

        try:
            self.coll_plot['Sg'][1]=[]
        except:
            self.coll_plot['Sg']=[[]]
            self.coll_plot['Sg'].append([])


        fig,ax,self.coll_plot['Sg'][1]=self.EMS.L.showSig(Sg,Tx=self.pos[tx],Rx=self.pos[rx],sr=True,fig=fig,ax=ax)


        if ion:
            try:
                [jj.remove() for jj in self.coll_plot['Sg'][0]]
            except:
                pass

            plt.draw()
            self.coll_plot['Sg'][0]=self.coll_plot['Sg'][1]



    def csv_save(self,filename,S):
        """
        save node positions into csv file

        Attributes:
        ----------
        
        filename : string 
                   name of the csv file
        S        : Simulation
                   Scipy.Simulation object

        """

        pos=np.array(nx.get_node_attributes(self,'p').values())
        pos=np.hstack((pos,np.zeros((len(self.nodes()),1))))  # passage en 3D
        pos=pos.reshape((1,len(self.nodes())*3))
        file=open('../save_data/' +filename +'.csv','a')
        file.write(str(S.now()) +',')
        np.savetxt(file,pos,delimiter=',')
        file.write('\n')
        file.close()


    def init_save(self,height=1.5):

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
            if self.node[pos[i][0]]['type'] =='ap':
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


        Attributes:
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


        Attributes:
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


        Attributes:
        ----------
        
        S        : Simulation
                   Scipy.Simulation object

        """

        pos=nx.get_node_attributes(self,'p')
        pe=nx.get_node_attributes(self,'pe_alg')
        typ = nx.get_node_attributes(self,'type')
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


    def dual_save(self,S):
        """
        DEPRECATED 

        REPLACED BY pylayers.util.save

        DEPRECATED 



        save txt 



        Attributes:
        ----------
        
        S        : Simulation
                   Scipy.Simulation object

        """

        pos=nx.get_node_attributes(self,'p')
        pclust = nx.get_node_attributes(self,'pe_clust')
        typ = nx.get_node_attributes(self,'type')
        if self.idx == 0:
            entete = 'Timestamp, True Position x, True Position y, Est Position1 x, Est Position1 y,Est Position2 x, Est Position2 y\n'
            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/pos.txt','write')
            file.write(entete)
            file.close()
            file2=open(basename+'/' + pstruc['DIRNETSAVE'] +'/rsslink.txt','write')
            entete2 = 'Timestamp, link, linkid, Pr, distance\n'
            file2.write(entete2)
            file2.close()
            file3=open(basename+'/' + pstruc['DIRNETSAVE'] +'/anchorposition.txt','write')
            data3 = 'node,pos x, pos y\n'
            file3.write(data3)
            for n in self.nodes():
                data3= n + ',' + str(self.node[n]['p'][0]) + ',' + str(self.node[n]['p'][1]) + '\n'
                file3.write(data3)
            file3.close()
            file4=open(basename+'/' + pstruc['DIRNETSAVE'] +'/toa.txt','w')
            entete4 = 'Timestamp, type, toaid, toa,distance\n'
            file4.write(entete4)
            file4.close()

        try:

            file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/pos.txt','a')
            file2=open(basename+'/' + pstruc['DIRNETSAVE'] +'/rsslink.txt','a')
            file4=open(basename+'/' + pstruc['DIRNETSAVE'] +'/toa.txt','a')
            for n in self.nodes():
                if n == '1':
                    data =  pyu.timestamp(S.now()) +','+ str(pos[n][0]) + ',' + str(pos[n][1]) + ',' + str(pclust[n][0,0]) + ',' + str(pclust[n][0,1]) + ',' + str(pclust[n][1,0]) + ',' + str(pclust[n][1,1]) +'\n'
                    for e in self.edge[n].keys():
                        if e != '6' and e !='7':
                            try:
                                data2 = data2 +',link,' + str(e) + ',' + str(self.edge[n][e]['rat1']['Pr'][0]) +',' + str(np.sqrt(np.sum((pos[n]-pos[e])**2)))
                            except:
                                data2 = pyu.timestamp(S.now()) + ',link,' + str(e) + ',' + str(self.edge[n][e]['rat1']['Pr'][0]) +',' + str(np.sqrt(np.sum((pos[n]-pos[e])**2)))
                        else :
                            try:
                                data4 = data4 +',toa,' + str(e) + ',' + str(self.edge[n][e]['rat1']['TOA'][0]) +',' + str(np.sqrt(np.sum((pos[n]-pos[e])**2)))
                            except:

                                data4 = pyu.timestamp(S.now()) + ',toa,' + str(e) + ',' + str(self.edge[n][e]['rat1']['TOA'][0]) +',' +str(np.sqrt(np.sum((pos[n]-pos[e])**2)))

            data2=data2 + '\n'
            data4=data4 + '\n'
            file.write(data)
            file2.write(data2)
            file4.write(data4)


            file.close()
            file2.close()
            file4.close()
            self.idx = self.idx +1
        except:
            pass


    def pyray_save(self,S):
        """
            save node positions into ini file, compliant with pyray standard

        Attributes:
        ----------
        
        filename : string 
                   name of the pyray file
        S        : Simulation
                   Scipy.Simulation object
        """


        assert len(self.SubNet.keys()) == 1 , NameError('when network.ini_save() \
        is used , only 1 rat must be involved in the Network.\
        Please modify agent.ini')


        height= 1.5
        pos=nx.get_node_attributes(self,'p').items()

        ### create ini files
        if self.idx == 0:
            self.init_save(height=height)
        ### save agent positions
        for i in range(len(pos)):
            if self.node[pos[i][0]]['type'] !='ap':
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
            if self.node[n]['type'] !='ap':
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
        Only links  which involve mobile nodes (type 'ag') are kept.

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
        is used , only 1 rat must be involved in the Network.\
        Please modify agent.ini')


        if self.idx == 0:
            file=open(pyu.getlong(filename ,'output'),'w')
        else:
            file=open(pyu.getlong(filename ,'output'),'a')

        config = ConfigParser.ConfigParser()
        timestamp = pyu.timestamp(S.now())
        config.add_section(timestamp)
        for e in self.edges():
            if not ((self.node[e[0][0]]['type'] == 'ap') and  (self.node[e[1][0]]['type'] == 'ap')):
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
        for rat in self.net.RAT.iterkeys():
            self.net.compute_LDPs(self.net.nodes(),rat)
        for n in self.net.nodes():
            self.net.node[n]['PN'].get_RAT()
            self.net.node[n]['PN'].get_SubNet()
            # Add access point position in each personal network (PN)
            [self.net.node[n]['PN'].node[n2].update({'pe':self.net.node[n2]['p']}) for n2 in self.net.node[n]['PN'].node.iterkeys() if self.net.node[n]['PN'].node[n2]['type'] == 'ap']
                
        ####################################################################################
        self.pos=self.net.get_pos()



        if 'csv' in self.save:
            nbnodes=len(self.net.nodes())
            entete = 'time'
            inode=self.net.nodes_iter()
            for i in inode:
                entete=entete + ',' + i +',,'
            entete=entete +'\n'
            file=open('../save_data/' +self.filename +'.csv','w')
            file.write(entete)
            file.close()



        while True:
            ############### compute LDP
            for rat in self.net.RAT.iterkeys():
                self.net.compute_LDPs(self.net.nodes(),rat)
            
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
#            if 'csv' in self.save:
#                self.net.csv_save(self.filename,self.sim)
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

            print 'network updated @',self.sim.now()

            self.net.idx=self.net.idx+1
            yield hold, self, self.net_updt_time



