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
        
        defaults={'PN':Network(),'method':'RGPA','rule':[Take_all()],'dc':[]}

        for key, value in defaults.items():
            if args.has_key(key):
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key]=value  
        self.args=args

        self.cla = CLA()

    def get_const(self,RAT=None,LDP=None):
        """ get constraints
            
            get the constraint of the networl followinf rule given in self.rule list.
            These rules are defined in Loca_Rule
        
        """
        self.dc= merge_rules(self,RAT=RAT,LDP=LDP)



    def fill_cla(self):
        for ldp in self.dc.keys():
            for rat in self.dc[ldp].keys():
                if ldp == 'Pr':
                    [self.cla.append(
                          RSS(value=self.dc[ldp][rat][i][1][0],
                          std=self.dc[ldp][rat][i][1][1],
                          vcw=1.0,
                          model={},
                          p=pan ) 
                                    ) for i in range(len(self.dc[ldp][rat]))]



                elif ldp == 'TOA':
                    pass


                elif ldp == 'TDOA':
                    pass



class PLocalization(Process):
    def __init__(self,loc=Localization(),loc_updt_time=.5,sim=None):
        Process.__init__(self,name='Location',sim=sim)
        self.loc = loc  
        self.loc_updt_time = loc_updt_time    
        
    def run(self):
        while True:

            self.loc.get_const()
            self.loc.fill_cla()

            print 'localization update @',self.sim.now()
            yield hold, self, self.loc_updt_time
