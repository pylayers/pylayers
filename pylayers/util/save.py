from SimPy.SimulationRT import Simulation, Process, hold
import numpy as np
import scipy as sp
import scipy.io as spio
import networkx as nx
import matplotlib.pyplot as plt
import ConfigParser

from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.network.network import Network, Node, PNetwork
from pylayers.gis.layout import Layout
import copy


import pickle

import pdb
import os




class Save(Process):
    """
    Save all variable of a simulnet simulation.
    Save process can be setup with the save.ini file from /<project>/ini
        
    Attributes
    ----------
    net : pylayers.network.network()
    sim : SimPy.SimulationRT()

    savemat : dictionnary with all the saved results from a simulation
               ( obtained after self.export() )

    Methods
    -------

    run ():
        save the current simulation evey k step (setup into save.ini)
    load():
        Load saved results of a simulation. file extension  .pck
        
    export(etype) :
        export the results into the etype format.
        available format : 
        - 'python'
        - 'matlab'


    """
    def __init__(self, **args):
        defaults = {'L': None,
                    'net': None,
                    'sim': None}


##       initialize attributes
        for key, value in defaults.items():
            if key in args:
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key] = value
        self.args = args

        Process.__init__(self, name='save', sim=self.args['sim'])


        self.C = ConfigParser.ConfigParser()
        self.C.read(pyu.getlong('save.ini','ini'))
        self.opt = dict(self.C.items('config'))
        self.pos = dict(self.C.items('position'))
        self.ldp = dict(self.C.items('ldp'))
        self.rat = dict(self.C.items('rat'))
        self.lpos=eval(self.pos['position'])
        self.lldp=eval(self.ldp['ldp'])
        self.lrat=eval(self.rat['rat'])


        self.sim=args['sim']
        self.net=args['net']


        self.save={}
        self.filename = eval(self.opt['filename'])
        self.file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/' +self.filename,'write')
        self.save['saveopt']={}
        self.save['saveopt']['lpos']=self.lpos
        self.save['saveopt']['lldp']=self.lldp
        self.save['saveopt']['lrat']=self.lrat
        pickle.dump(self.save, self.file)
        self.file.close()
        self.idx=0

    def load(self,filename=[]):
        """
        Load a saved trace simulation

        Examples
        --------

        >>> from pylayers.util.save import *
        >>> S=Save()
        >>> S.load()
        """

        if filename == []:
            filename = self.filename

        out=[0]
        infile = open(basename+'/' + pstruc['DIRNETSAVE'] +'/'+filename, 'r')
        while 1:
            try:
                out.append(pickle.load(infile))
            except (EOFError, pickle.UnpicklingError):
                break
            out.pop(0)
        infile.close()
        dout= dict(out[-1])
        return dout


    def mat_export(self):
        """
            export save simulation to a matlab file

        Examples
        --------
        >>> from pylayers.util.save import *
        >>> S=Save()
        >>> S.mat_export()

        """
        self.save=self.load()
        self.savemat=copy.deepcopy(self.save)
        nodes=self.save['saveopt']['type'].keys()
        for inn,n in enumerate(nodes):
            self.savemat['node_'+n]=self.save[n]
            for n2 in nodes:
                if n2 != n:
                    self.savemat['node_'+n]['node_'+n2]=self.save[n][n2]
                    del self.savemat[n][n2]
            del self.savemat[n]

            for o in self.save['saveopt']:
                if o =='subnet' and inn == 0:
                    for r in self.save['saveopt']['lrat']:
                        li=self.save['saveopt'][o][r]
                        self.savemat['saveopt'][o][r]=['node_'+l for l in li]
                
                else :
                    try:
                        self.savemat['saveopt'][o]['node_'+n]=self.save['saveopt'][o][n]
                        del self.savemat['saveopt'][o][n]
                    except:
                        pass

        spio.savemat(basename+'/' + pstruc['DIRNETSAVE'] +'/' +self.filename,self.savemat)


    def run(self):
        """
            Run the save Result process
        """


        ### init save dictionnary
        self.save['saveopt']['Layout'] = self.L.filename
        self.save['saveopt']['type']=nx.get_node_attributes(self.net,'type')
        self.save['saveopt']['epwr']=nx.get_node_attributes(self.net,'epwr')
        self.save['saveopt']['sens']=nx.get_node_attributes(self.net,'sens')


        self.save['saveopt']['subnet']={}
        for rat in self.lrat:
            self.save['saveopt']['subnet'][rat]=self.net.SubNet[rat].nodes()

        [self.save.update({n:{}}) for n in self.net.nodes()]

        # find the size of save array regarding the simulation duration and the save sample time
        nb_sample=np.ceil(eval(self.sim.sim_opt['duration'])/eval(self.opt['save_update_time']))+1


        # create void array to be fill with simulation data
        for n in self.net.nodes():
            for position in self.lpos:
                self.save[n][position]=np.zeros((nb_sample,2))*np.nan


        for e in self.net.edges():
            self.save[e[0]][e[1]]={}
            self.save[e[1]][e[0]]={}
            for rat in self.lrat:
                self.save[e[0]][e[1]][rat]={}
                self.save[e[1]][e[0]][rat]={}
                for ldp in self.lldp:
                    self.save[e[0]][e[1]][rat][ldp]=np.zeros((nb_sample,2))*np.nan
                    self.save[e[1]][e[0]][rat][ldp]=np.zeros((nb_sample,2))*np.nan


        while True:
            rl={}
            for rat in self.lrat:
                for ldp in self.lldp:
                    rl[rat+ldp]=nx.get_edge_attributes(self.net.SubNet[rat],ldp)

            for n in self.net.nodes():
                for position in self.lpos:
                    try:
                        p = nx.get_node_attributes(self.net,position)
                        self.save[n][position][self.idx]=p[n]
                    except:
                        pass

            for e in self.net.edges():
                for rat in self.lrat:
                    for ldp in self.lldp:
                        self.save[e[0]][e[1]][rat][ldp][self.idx]=rl[rat+ldp][e]
                        self.save[e[1]][e[0]][rat][ldp][self.idx]=rl[rat+ldp][e]

            self.file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/' +self.filename,'a')
            pickle.dump(self.save, self.file)
            self.file.close()
            self.idx=self.idx+1
            yield hold, self, eval(self.opt['save_update_time'])


