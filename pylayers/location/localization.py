from pylayers.util.project import *
import sys

import SimPy.Simulation
from SimPy.Simulation import Process,hold

from pylayers.util import utilnet
from pylayers.network.network import Network, Node
from pylayers.location.locarule import Take_all,  merge_rules

import pdb

#mport CSP10
from pylayers.location.geometric.constraints.cla import *
from pylayers.location.geometric.constraints.rss import *
from pylayers.location.geometric.constraints.toa import *
from pylayers.location.geometric.constraints.tdoa import *
from pylayers.location.geometric.constraints.exclude import *


class Localization(object):

    def __init__(self,**args):
        
        defaults={'PN':Network(),'method':'RGPA','rule':[Take_all()],'dc':[],'ID':'0'}

        for key, value in defaults.items():
            if args.has_key(key):
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key]=value  
        self.args=args
        self.cla = CLA()

#    def get_const(self,RAT=None,LDP=None):
#        """ get constraints
#            
#            get the constraint of the network following rule given in self.rule list.
#            These rules are defined in Loca_Rule
#        
#        """
#        self.dc= merge_rules(self,RAT=RAT,LDP=LDP)



    def fill_cla(self):
        """
            Fill the constraint layer array
        """

        for e in self.PN.edge[self.ID].keys():
            for rat in self.PN.edge[self.ID][e].keys():
                try:
                    self.cla.append(
                        RSS(id = rat+'-'+self.ID+'-'+e,
                            value = self.PN.edge[self.ID][e][rat]['Pr'][0],
                            std = self.PN.edge[self.ID][e][rat]['Pr'][1],
                            vcw = 1.0,
                            model={},
                            p = self.PN.node[e]['pe'],
                            origin={'id':self.ID,'link':[e],'rat':rat,'ldp':'Pr'}
                            )
                                    )
                except:
                    pass


#                elif ldp == 'TOA':
#                    pass


#                elif ldp == 'TDOA':
#                    pass




    def update(self,RAT=None,LDP=None):
        for c in self.cla.c:
            rat,ldp,e,own=c.origin.values()
            c.p = self.PN.node[e[0]]['pe']
            c.value = self.PN.edge[e[0]][self.ID][rat][ldp][0]
            c.std = self.PN.edge[e[0]][self.ID][rat][ldp][1]
        self.cla.update()


class PLocalization(Process):
    def __init__(self,loc=Localization(),loc_updt_time=.5,sim=None):
        Process.__init__(self,name='Location',sim=sim)
        self.loc = loc  
        self.loc_updt_time = loc_updt_time    
        
    def run(self):
#        self.loc.get_const()
        self.loc.fill_cla()
        while True:
            self.loc.update()
            print 'localization update @',self.sim.now()
            yield hold, self, self.loc_updt_time
