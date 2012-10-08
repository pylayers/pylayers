from pylayers.util.project import *
import sys

import SimPy.Simulation
from SimPy.Simulation import Process,hold

from pylayers.util import utilnet
from pylayers.network.network import Network, Node
from pylayers.location.locarule import Take_all,  merge_rules

import pdb

from pylayers.location.geometric.constraints.cla import *
from pylayers.location.geometric.constraints.rss import *
from pylayers.location.geometric.constraints.toa import *
from pylayers.location.geometric.constraints.tdoa import *
from pylayers.location.geometric.constraints.exclude import *


class Localization(object):

    def __init__(self,**args):
        
        defaults={'PN':Network(),'net':Network(),'method':'RGPA','rule':[Take_all()],'dc':[],'ID':'0'}

        for key, value in defaults.items():
            if key in args:
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key] = value
        self.args = args
        self.cla = CLA()

#    def get_const(self, RAT=None, LDP=None):
#        """ get constraints
#
#   get the constraint of the networl followinf rule given in self.rule list.
#            These rules are defined in Loca_Rule
#
#                args[key]=value  
#
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

        for e in self.net.node[self.ID]['PN'].edge[self.ID].keys():
            for rat in self.net.node[self.ID]['PN'].edge[self.ID][e].keys():
                try:
                    self.cla.append(
                        RSS(id = rat+'-Pr-'+self.ID+'-'+e,
                            value = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['Pr'][0],
                            std = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['Pr'][1],
                            model={},
                            p = self.net.node[self.ID]['PN'].node[e]['pe'],
                            origin={'id':self.ID,'link':[e],'rat':rat,'ldp':'Pr'}
                            )
                                    )
                except:
                    pass

                try:
                    self.cla.append(
                        TOA(id = rat+'-TOA-'+self.ID+'-'+e,
                            value = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['TOA'][0]*0.3,
                            std = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['TOA'][1]*0.3,
                            p= self.net.node[self.ID]['PN'].node[e]['pe'],
                            origin={'id':self.ID,'link':[e],'rat':rat,'ldp':'TOA'}
                            )
                                    )
                except:
                    pass

#                elif ldp == 'TDOA':
#                    pass


    def update(self,rat='all',ldp='all'):
        if rat == 'all':
            rat=self.net.node[self.ID]['PN'].SubNet.keys()
        if ldp == 'all':
            ldp=['Pr','TOA','TDOA']
        for c in self.cla.c:
            crat,cldp,e,own=c.origin.values()
            if (crat in rat) and (cldp in ldp) :
                c.p = self.net.node[self.ID]['PN'].node[e[0]]['pe']
                c.value = self.net.node[self.ID]['PN'].edge[e[0]][self.ID][crat][cldp][0]
                c.std = self.net.node[self.ID]['PN'].edge[e[0]][self.ID][crat][cldp][1]
        self.cla.update()


    def savep(self):
        self.net.node[self.ID]['PN'].update_pos(self.ID,self.cla.pe,p_pe='pe')
        self.net.update_pos(self.ID,self.cla.pe,p_pe='pe')

class PLocalization(Process):
    def __init__(self, loc=Localization(), loc_updt_time=.5, sim=None):
        Process.__init__(self, name='Location', sim=sim)
        self.loc = loc
        self.loc_updt_time = loc_updt_time

    def run(self):
#        self.loc.get_const()
        self.loc.fill_cla()
        while True:

            self.loc.update()
            self.loc.cla.merge2()
            self.loc.cla.refine(self.loc.cla.Nc)
            self.loc.cla.estpos2()
            self.loc.savep()
            self.loc.cla.Nc=len(self.loc.cla.c)

            print '-------------LOC--------------------------'
            for i in self.loc.cla.c:
                print i.id,i.value
            print '--------------------------------------------'
            print 'localization update @',self.sim.now()
            yield hold, self, self.loc_updt_time
