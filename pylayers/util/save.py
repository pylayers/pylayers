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
        defaults = {'net': None,
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
        saveopt=dout['saveopt']
        dout.pop('saveopt')
        return saveopt,dout

    def export(self,etype='python'):
        """
        Export simulation results into a given extension

        Parameters
        ----------
        etype:
            'python' or ''matlab

        Examples
        --------
        >>> from pylayers.util.save import *
        >>> S=Save()
        >>> S.export('matlab')
        """


        savecfg,d=self.load(self.filename+'.tmp')
        # need to sort the time stamp 
        TS = np.array((d.keys()))
        sTSi=np.argsort(TS)
        sTS=np.sort(TS)

        self.savemat={}
        self.savemat['timestamp']=sTS


        ##### init
        for n in savecfg['type'].keys():
            if etype == 'matlab':
                dkey='node_'+n
            else :
                dkey=n
            self.savemat[dkey]={}
            self.savemat[dkey]['type']=savecfg['type'][n]
            self.savemat[dkey]['sens']=savecfg['sens'][n]
            self.savemat[dkey]['epwr']=savecfg['epwr'][n]
            for p in savecfg['lpos']:
                self.savemat[dkey][p]=[]
            for r in savecfg['lrat']:
                # test if node has the rat and if it is not the only one on that rat
                if (n in savecfg['subnet'][r]) and  (len(savecfg['subnet'][r]) > 1):
                    self.savemat[dkey][r]={}
                    for l in savecfg['lldp']:
                        self.savemat[dkey][r][l]={}

    
        ### fill in dict
        for t in sTS:
            for p in savecfg['lpos']:
                for n in d[t][p].keys():
                    # handle matlab struct name cannot be a number
                    if etype == 'matlab':
                        dkey='node_'+n
                    else :
                        dkey=n


                    ########### Position
                    try :
                        if p != 'pe_clust':
                            self.savemat[dkey][p]=np.vstack((self.savemat[dkey][p],d[t][p][n]))
                        else :
                            self.savemat[dkey][p]=np.dstack((self.savemat[dkey][p],d[t][p][n]))
                    except:
                        self.savemat[dkey][p]=d[t][p][n]
#                        except:
#                            self.savemat[dkey]={}
#                            for r in savecfg['lrat']:
#                                self.savemat[dkey][r]={}
#                                for l in savecfg['lldp']:
#                                    self.savemat[dkey][r][l]={}
#                            self.savemat[dkey][p]=d[t][p][n]

                   ############### links
        for t in sTS:
            for r in savecfg['lrat']:
                for l in savecfg['lldp']:
                    for ii in d[t][r][l].keys():
                        for n in savecfg['type'].keys():
                            # edge value (node1,node2) if n == node1
                            if etype == 'matlab':
                                dkey='node_'+n
                            else :
                                dkey=n
                            if n in ii[0]:
                                if etype == 'matlab':
                                    dkey2='node_'+ii[1]
                                else :
                                    dkey2=ii[1]

                                try:
                                    self.savemat[dkey][r][l][dkey2]=np.vstack((self.savemat[dkey][r][l][dkey2],d[t][r][l][ii]))
                                except:
                                    self.savemat[dkey][r][l][dkey2]=d[t][r][l][ii]

                            # if n == node2
                            elif n in ii[1]:
                                if etype == 'matlab':
                                    dkey2='node_'+ii[0]
                                else :
                                    dkey2=ii[0]

                                try:
                                    self.savemat[dkey][r][l][dkey2]=np.vstack((self.savemat[dkey][r][l][dkey2],d[t][r][l][ii]))
                                except:
                                    self.savemat[dkey][r][l][dkey2]=d[t][r][l][ii]
                    

        if  etype == 'matlab':
            spio.savemat(basename+'/' + pstruc['DIRNETSAVE'] +'/' +self.filename,self.savemat)
        if  etype == 'python':

            self.savemat['saveopt']=savecfg
            self.file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/' +self.filename+'.pck','w')
            pickle.dump(self.savemat, self.file)
            self.file.close()



    def run(self):
        """
            Run the save Result process
        """
        self.save['saveopt']['type']=nx.get_node_attributes(self.net,'type')
        self.save['saveopt']['epwr']=nx.get_node_attributes(self.net,'epwr')
        self.save['saveopt']['sens']=nx.get_node_attributes(self.net,'sens')
        self.save['saveopt']['subnet']={}
        for rat in self.lrat:
            self.save['saveopt']['subnet'][rat]=self.net.SubNet[rat].nodes()

        while True:
            self.save[self.sim.now()]={}
            for position in self.lpos:
                self.save[self.sim.now()][position]=nx.get_node_attributes(self.net,position)
            for rat in self.lrat:
                self.save[self.sim.now()][rat]={}
                for ldp in self.lldp:
                    self.save[self.sim.now()][rat][ldp]=nx.get_edge_attributes(self.net.SubNet[rat],ldp)
            self.file=open(basename+'/' + pstruc['DIRNETSAVE'] +'/' +self.filename+'.tmp','a')
            pickle.dump(self.save, self.file)
            self.file.close()
            yield hold, self, eval(self.opt['save_update_time'])


