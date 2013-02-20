from pylayers.util.project import *
import sys


from SimPy.SimulationRT import Process,hold

from pylayers.util import utilnet
from pylayers.network.network import Network, Node
from pylayers.location.locarule import Take_all,  merge_rules

import pdb

from pylayers.location.geometric.constraints.cla import *
from pylayers.location.geometric.constraints.rss import *
from pylayers.location.geometric.constraints.toa import *
from pylayers.location.geometric.constraints.tdoa import *
from pylayers.location.geometric.constraints.exclude import *
from pylayers.location.algebraic.algebraic import *

from   pylayers.network.model import Model
import networkx as nx

class Localization(object):

    def __init__(self,**args):
        
        defaults={'PN':Network(),'net':Network(),'method':['geo','alg'],'model':{},'rule':[Take_all()],'dc':[],'ID':'0','save':[]}

        for key, value in defaults.items():
            if key in args:
                setattr(self, key, args[key])
            else:
                setattr(self, key, value)
                args[key] = value
        self.args = args
        self.config = ConfigParser.ConfigParser()
        self.config.read(pyu.getlong('EMSolver.ini', 'ini'))

        self.cla = CLA()
        self.algloc=algloc()
        self.idx = 0
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
        ## loop on edges
        for e in self.net.node[self.ID]['PN'].edge[self.ID].keys():
            ## loop on rat
            for rat in self.net.node[self.ID]['PN'].edge[self.ID][e].keys():
                try:
                    param = dict(self.config.items(rat+'_PLM'))
                    self.cla.append(
                        RSS(id = rat+'-Pr-'+self.ID+'-'+e,
                            value = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['Pr'][0],
                            std = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['Pr'][1],
                            model=Model(f=eval(param['f']), rssnp=eval(param['rssnp']), d0=eval(param['d0']), method=param['method']),
                            p = self.net.node[self.ID]['PN'].node[e]['pe'],
                            origin={'id':self.ID,'link':[e],'rat':rat,'ldp':'Pr'}
                            )
                                    )
                except:
                    param = dict(self.config.items(rat+'_PLM'))
                    self.cla.append(
                        RSS(id = rat+'-Pr-'+self.ID+'-'+e,
                            p = self.net.node[self.ID]['PN'].node[e]['pe'],
                            origin={'id':self.ID,'link':[e],'rat':rat,'ldp':'Pr'}
                            )
                                    )

                try:
                    self.cla.append(
                        TOA(id = rat+'-TOA-'+self.ID+'-'+e,
                            value = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['TOA'][0],
                            std = self.net.node[self.ID]['PN'].edge[self.ID][e][rat]['TOA'][1],
                            p= self.net.node[self.ID]['PN'].node[e]['pe'],
                            origin={'id':self.ID,'link':[e],'rat':rat,'ldp':'TOA'}
                            )
                                    )
                except:
                    self.cla.append(
                        TOA(id = rat+'-TOA-'+self.ID+'-'+e,
                            p= self.net.node[self.ID]['PN'].node[e]['pe'],
                            origin={'id':self.ID,'link':[e],'rat':rat,'ldp':'TOA'}
                            )
                                    )


#                 elif ldp == 'TDOA':
#                    pass




    def update(self,rat='all',ldp='all'):
        """
            update constraints information (anchor position, value and std) 
            from the network graph
        """
        if rat == 'all':
            rat=self.net.node[self.ID]['PN'].SubNet.keys()
        if ldp == 'all':
            ldp=['Pr','TOA','TDOA']
        else:
            if not isinstance(ldp,list):
                ldp =[ldp]

        self.algloc.nodes={}
        self.algloc.ldp={}
        for c in self.cla.c:
            crat,cldp,e,own=c.origin.values()
            if (crat in rat) and (cldp in ldp) :
                pos = self.net.node[self.ID]['PN'].node[e[0]]['pe']
                value = self.net.node[self.ID]['PN'].edge[self.ID][e[0]][crat][cldp][0]
                std = self.net.node[self.ID]['PN'].edge[self.ID][e[0]][crat][cldp][1]
                if 'geo' in self.method :
                    c.p = pos
                    c.value = value
                    c.std = std
                if 'alg' in self.method :
                    if c.runable:
                        if c.type == 'TOA':
                            try :
                                self.algloc.nodes['RN_TOA'] = np.vstack((self.algloc.nodes['RN_TOA'],pos))
                                self.algloc.ldp['TOA'] = np.hstack((self.algloc.ldp['TOA'],value))
                                self.algloc.ldp['TOA_std'] = np.hstack((self.algloc.ldp['TOA_std'],std))

                            except:
                                self.algloc.nodes['RN_TOA'] = pos
                                self.algloc.ldp['TOA'] = value
                                self.algloc.ldp['TOA_std'] = std
                        if c.type == 'Pr':
                            self.algloc.nodes['RN_RSS'] = pos.T
                            self.algloc.ldp['RSS'] = value
                            self.algloc.ldp['RSS_std'] = std.T
                            self.algloc.ldp['d0'] = c.param['d0']
                            ####### -pl0 from alg loc ############
                            self.algloc.ldp['PL0'] = -c.param['PL0']
            else:
                c.runable = False

        try:
            self.algloc.nodes['RN_TOA']=self.algloc.nodes['RN_TOA'].T
        except:
            pass
        self.cla.update()

    def savep(self,value,name='pe'):
        """
            write an estimated position into self.net
        """
        self.net.node[self.ID]['PN'].update_pos(self.ID,value,p_pe=name)
        self.net.update_pos(self.ID,value,p_pe=name)



    def compute_geo(self,rat='all',ldp='all',pe=True):
        """
            Compute postion with the geometric algorithm
        """


        if sum(self.cla.runable) >= 2:
            cpe = self.cla.compute(pe=pe)
            if cpe:
                self.savep(self.cla.pe,name='pe')
                self.savep(self.cla.pe,name='pe_geo')

        # in case of lack of observables
        elif sum(self.cla.runable) >= 1:
            cpe = self.cla.compute_amb(pe=pe)
            if cpe:
                self.savep(np.array(self.cla.pecluster),name='pe_clust')

    def compute_alg(self,rat='all',ldp='all',pe=True):
        """
            Compute postion with the algebraic algorithm
        """

        if len(self.cla.c) !=0:
            if ldp == 'all':
                ldp=['Pr','TOA','TDOA']
            elif not isinstance(ldp,list):
                ldp=[ldp]

            pe_alg = self.algloc.wls_locate('Pr' in ldp, 'TOA' in ldp, 'TDOA' in ldp, 'mode')
            self.savep(pe_alg.T,name='pe_alg')




    def compute_crb(self,rat='all',ldp='all',pe=True):
        """
            Compute CramerRao bound
        """
        
        if ldp == 'all':
            ldp=['Pr','TOA','TDOA']
        crb = self.algloc.crb(np.zeros((2,1)),'Pr' in ldp, 'TOA' in ldp, 'TDOA' in ldp)
        self.savep(np.array(crb),name='crb')




class PLocalization(Process):
    def __init__(self, loc=Localization(), loc_updt_time=.5, sim=None):
        Process.__init__(self, name='Location', sim=sim)
        self.loc = loc
        self.loc_updt_time = loc_updt_time
        self.method = self.loc.method
        self.sim = sim
    def run(self):
#        self.loc.get_const()
        self.loc.fill_cla()
        while True:
            self.loc.update(ldp='TOA')
            if 'geo'in self.method :
                self.loc.compute_geo(ldp='TOA')
            if 'alg'in self.method :
                self.loc.compute_alg(ldp='TOA')
            if self.sim.verbose:
                print 'localization node',self.loc.ID, ' update @',self.sim.now()
            yield hold, self, self.loc_updt_time
