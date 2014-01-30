from pylayers.util.project import *
import sys


from SimPy.SimulationRT import Process,hold

from pylayers.util import utilnet
from pylayers.network.network import Network, Node
#from pylayers.location.locarule import Take_all,  merge_rules

import pdb

from pylayers.location.geometric.constraints.cla import *
from pylayers.location.geometric.constraints.rss import *
from pylayers.location.geometric.constraints.toa import *
from pylayers.location.geometric.constraints.tdoa import *
from pylayers.location.geometric.constraints.exclude import *
from pylayers.location.algebraic.algebraic import *

from pylayers.network.communication import Gcom, TX, RX

from   pylayers.network.model import PLSmodel
import networkx as nx

class Localization(object):
    """

    """

    def __init__(self,**args):

        defaults={'PN':Network(),'net':Network(),'method':['geo','alg'],'model':{},'ID':'0','save':[]}

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

    def __repr__(self):
        s = 'Localization information\n*************************\n'
        s = s + '\nNode ID: ' + str(self.ID)
        s = s + '\nLocalization methods: ' + str(self.method)
        s = s + '\n\n CLA:\n'
        s = s + self.cla.__repr__() + '\n'

        return s

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
                            model = PLSmodel(f = eval(param['f']), 
                                             rssnp = eval(param['rssnp']),
                                             d0 = eval(param['d0']), 
                                             method = param['method']),
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
                if self.net.node[self.ID]['PN'].edge[self.ID][e[0]][crat]['vis']:
                    c.visible = True
                    pos = self.net.node[self.ID]['PN'].node[e[0]]['pe']
                    value = self.net.node[self.ID]['PN'].edge[self.ID][e[0]][crat][cldp][0]
                    std = self.net.node[self.ID]['PN'].edge[self.ID][e[0]][crat][cldp][1]
                    if 'geo' in self.method :
                        # methods from constraint.py
                        c.updc('p',pos)
                        c.updc('value',value)
                        c.updc('std',std)

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
                else :
                    c.visible = False
            else :
                c.usable = False

        try:
            self.algloc.nodes['RN_TOA']=self.algloc.nodes['RN_TOA'].T
        except:
            pass
        self.cla.update()

        # ############"" FOR DEBUG
        # prss=np.where(np.array(self.cla.type)=='RSS')[0]
        # for x in range(len(self.cla.c)):
        #     if x in prss:
        #         self.cla.usable[x]=False


    def savep(self,value,now=0.,name='pe'):
        """
            write an estimated position into self.net
        """
        self.net.node[self.ID]['PN'].update_pos(self.ID,value,p_pe=name)
        self.net.update_pos(self.ID,value,now=now,p_pe=name)



    def compute_geo(self,rat='all',ldp='all',now=0.,pe=True):
        """
            Compute position with the geometric algorithm

        Returns
        -------
            True if estimated position has been computed
        """

        if sum(self.cla.usable) >= 2:
            cpe = self.cla.compute(pe=pe)
            if cpe:
                self.savep(self.cla.pe,name='pe')
                self.savep(self.cla.pe,name='pe_geo')
                self.net.node[self.ID]['PN'].node[self.ID]['te']=now # estimated time
                return True
            return False
        # in case of lack of observables
        elif sum(self.cla.usable) >= 1:
            cpe = self.cla.compute(pe=pe)
            if cpe:
                self.savep(np.array(self.cla.pecluster),now=now,name='pe_clust')
                self.net.node[self.ID]['PN'].node[self.ID]['te']=now # estimated time
                return True
            return False
        else :
            return False

    def compute_alg(self,rat='all',ldp='all',now=0.,pe=True):
        """
            Compute position with the algebraic algorithm
        Returns
        -------
            True if estimatited position has been computed
        """

        if sum(self.cla.usable) >= 2 :
            if ldp == 'all':
                ldp=['Pr','TOA','TDOA']
            elif not isinstance(ldp,list):
                ldp=[ldp]

            pe_alg = self.algloc.wls_locate('Pr' in ldp, 'TOA' in ldp, 'TDOA' in ldp, 'mode')
            self.savep(pe_alg.T,now=now,name='pe_alg')
            self.net.node[self.ID]['PN'].node[self.ID]['te']=now
            return True
        else:
            return False



#    def compute_crb(self,rat='all',ldp='all',now=0.,pe=True):
#        """
#            Compute CramerRao bound
#        """
#
#        if ldp == 'all':
#            ldp=['Pr','TOA','TDOA']
#        crb = self.algloc.crb(np.zeros((2,1)),'Pr' in ldp, 'TOA' in ldp, 'TDOA' in ldp)
#        self.savep(np.array(crb),now=now,name='crb')




class PLocalization(Process):
    def __init__(self, loc=Localization(), tx=TX() ,loc_updt_time=.5, sim=None):
        Process.__init__(self, name='Location', sim=sim)
        self.loc = loc
        self.tx=tx
        self.loc_updt_time = loc_updt_time
        self.method = self.loc.method
        self.sim = sim




    def run(self):
#        self.loc.get_const()
        self.loc.fill_cla()
        self.loc.update(ldp='TOA')
        while True:


            # if no previous position have been computed or if position is obsolete
            bep=False
            if self.loc.net.node[self.loc.ID]['pe'].size == 0 :
                self.tx.cmdrq.signal()
                self.loc.update(ldp='TOA')
            # try to obtain an estimated position
            if 'geo' in self.method :
                bep = self.loc.compute_geo(ldp='TOA',now=self.sim.now())
#                    if 'alg' in self.method and bep:
#                        bep = self.loc.compute_alg(ldp='TOA',now=self.sim.now())

#                     if no position has been computed

#            if not bep or (self.sim.now() - self.loc.net.node[self.loc.ID]['PN'].node[self.loc.ID]['te']>self.loc_updt_time):
            if self.sim.verbose:
                print 'localization request communication from node',self.loc.ID, '@',self.sim.now()
            self.tx.cmdrq.signal()
            self.loc.update(ldp='TOA')
            

            if bep and self.sim.verbose :
                print 'LOCALIZATION node',self.loc.ID, ' update @',self.sim.now()
            yield hold, self, self.loc_updt_time
