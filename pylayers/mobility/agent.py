from SimPy.SimulationRT import Simulation
from pylayers.mobility.transit.Person import Person
from pylayers.mobility.transit.World import world
from pylayers.mobility.transit.SteeringBehavior import Seek, Separation, Containment, InterpenetrationConstraint, queue_steering_mind, default_steering_mind


import numpy as np
import networkx as nx
import time
import ConfigParser
import pandas as pd
import pylayers.util.pyutil as pyu
from pylayers.network.network import Node, Network
from pylayers.network.communication import Gcom, TX, RX
from pylayers.location.localization import Localization, PLocalization
from pylayers.gis.layout import Layout
from pylayers.util.utilnet import *
#from pylayers.util.pymysqldb import Database
import pdb


class Agent(object):

    """   Agent 

    Members
    -------

    args
    ID
    name
    typ
    net 
    epwr 
    gcom
    sim 
    wstd 
    sens
    dcond 
    meca : transit.Person
    net : pylayers.network.Network 
    sim : 
    PN :    
    rxt
    rxr 



    """

    def __init__(self, **args):
        """ Mobile Agent Init

           Parameters
           ----------

           'ID': string
                agent ID
           'name': string
                Agent name
           'typ': string
                agent typ . 'ag' for moving agent, 'ap' for static acces point
           'pos' : np.array([])
                numpy array containing the initial position of the agent
           'roomId': int
                Room number where the agent is initialized (Layout.Gr)
           'meca_updt': float
                update time interval for the mechanical process
           'loc': bool
                enable/disable localization process of the agent
           'loc_updt': float
                update time interval for localization process
           'L': pylayers.gis.Layout()
           'net':pylayers.network.Network(),
           'wstd': list of string
                list of used radio access techology of the agent
           'world': transit.world()
                Soon deprecated
           'save': list of string
                list of save method ( soon deprecated)
           'sim':Simpy.SimulationRT.Simulation(),
           'epwr': dictionnary
                dictionnary of emmited power of transsmitter{'wstd#':epwr value}
           'sens': dictionnary
                dictionnary of sensitivity of reveicer {'wstd#':sens value}
           'dcond': dictionnary
                Not used yet
           'gcom':pylayers.communication.Gcom()
                Communication graph
           'comm_mod': string
                Communication between nodes mode:
                'autonomous': all TOAs are refreshed regulary
                'synchro' : only visilbe TOAs are refreshed
        """
        defaults = {'ID': '0',
                    'name': 'johndoe',
                    'typ': 'ag',
                    'color': 'k',
                    'pdshow': False,
                    'pos': np.array([]),
                    'roomId': -1,
                    'froom': [],
                    'wait': [],
                    'seed': 0,
                    'cdest': 'random',
                    'meca_updt': 0.1,
                    'loc': False,
                    'loc_updt': 0.5,
                    'loc_method': ['geo'],
                    'L': Layout(),
                    'network': True,
                    'net': Network(),
                    'wstd': ['rat1'],
                    'world': world(),
                    'save': [],
                    'sim': Simulation(),
                    'epwr': {},
                    'sens': {},
                    'dcond': {},
                    'gcom': Gcom(),
                    'comm_mode': 'autonomous'}

        for key, value in defaults.items():
            if key not in args:
                args[key] = value

        self.args = args
        self.ID = args['ID']
        self.name = args['name']
        self.typ = args['typ']
        # Create Network
        self.net = args['net']
        self.epwr = args['epwr']
        self.gcom = args['gcom']
        self.sim = args['sim']
        self.wstd = args['wstd']
        if args['epwr'] == {}:
            self.epwr = {x: 0 for x in self.wstd}
        else:
            self.epwr = args['epwr']

        if args['sens'] == {}:
            self.sens = {x: -180 for x in self.wstd}
        else:
            self.sens = args['sens']

        try:
            self.dcond = args['dcond']
        except:
            pass

        # check if node id already given
        if self.ID in self.net.nodes():
            raise NameError(
                'another agent has the ID: ' + self.ID + ' .Please use an other ID')

        if self.typ == 'ag':
            # mechanical init
            self.meca = Person(ID=self.ID,
                               color=args['color'],
                               pdshow=args['pdshow'],
                               roomId=args['roomId'],
                               L=args['L'],
                               net=self.net,
                               interval=args['meca_updt'],
                               wld=args['world'],
                               sim=args['sim'],
                               seed=args['seed'],
                               moving=True,
                               froom=args['froom'],
                               wait=args['wait'],
                               cdest=args['cdest'],
                               save=args['save']
                               )
            self.meca.behaviors = [Seek(), Containment(),
                                   Separation(), InterpenetrationConstraint()]
            self.meca.steering_mind = queue_steering_mind
            # Network init
            self.node = Node(ID=self.ID,name=self.name, p=conv_vecarr(self.meca.position),
                             t=self.sim.now(), wstd=args['wstd'],
                             epwr=self.epwr, sens=self.sens, typ=self.typ)
            self.net.add_nodes_from(self.node.nodes(data=True))

            self.sim.activate(self.meca, self.meca.move(), 0.0)
            self.PN = self.net.node[self.ID]['PN']

            # Communication init

            if args['comm_mode'] == 'synchro' and args['network']:
                # The TOA requests are made every refreshTOA time ( can be modified in agent.ini)
                # This Mode will be deprecated in future version

                self.rxr = RX(net=self.net,
                              ID=self.ID,
                              dcond=self.dcond,
                              gcom=self.gcom,
                              sim=self.sim)

                self.rxt = RX(net=self.net,
                              ID=self.ID,
                              dcond=self.dcond,
                              gcom=self.gcom,
                              sim=self.sim)

                self.sim.activate(self.rxr, self.rxr.refresh_RSS(), 0.0)
                self.sim.activate(self.rxt, self.rxt.refresh_TOA(), 0.0)

            elif args['comm_mode'] == 'autonomous' and args['network']:
                # The requests are made by node only when they are in
                # visibility of pairs.

                # self.rxr only manage a refresh RSS process
                self.rxr = RX(net=self.net, ID=self.ID,
                              gcom=self.gcom, sim=self.sim)
                # self.tx manage all requests to other nodes
                self.tx = TX(net=self.net, ID=self.ID,
                             gcom=self.gcom, sim=self.sim)
                # self.tx replies to  requests from self.tx
                self.rx = RX(net=self.net, ID=self.ID,
                             gcom=self.gcom, sim=self.sim)

                self.sim.activate(self.rxr, self.rxr.refresh_RSS(), 0.0)
                self.sim.activate(self.tx, self.tx.request(), 0.0)
                self.sim.activate(self.rx, self.rx.wait_request(), 0.0)

        elif self.typ == 'ap':
            if args['roomId'] == -1:
                self.node = Node(ID=self.ID, p=self.args['pos'],
                                 t=self.sim.now(), wstd=args['wstd'],
                                 epwr=self.epwr, sens=self.sens, typ=self.typ)
            else:
                pp = np.array(args['L'].Gr.pos[self.args['roomId']])
                self.node = Node(
                    ID=self.ID, p=pp, t=self.sim.now(), wstd=args['wstd'],
                    epwr=self.epwr, sens=self.sens, typ=self.typ)
            self.net.add_nodes_from(self.node.nodes(data=True))
            self.sim = args['sim']

            self.PN = self.net.node[self.ID]['PN']
            self.PN.node[self.ID]['pe'] = self.net.node[self.ID]['p']
            if args['comm_mode'] == 'autonomous' and args['network']:
                self.rx = RX(net=self.net, ID=self.ID,
                             gcom=self.gcom, sim=self.sim)
                self.sim.activate(self.rx, self.rx.wait_request(), 0.0)

            p = self.args['pos']
            self.posdf = pd.DataFrame(
                {'t': pd.Timestamp(0), 'x': p[0], 'y': p[1], 'z': p[2],
                                       'vx': np.array([0.0]), 'vy': np.array([0.0]),
                                       'ax': np.array([0.0]), 'ay': np.array([0.0]),
                 }, columns=['t', 'x', 'y', 'z', 'vx', 'vy', 'ax', 'ay'], index=np.array([0]))

        else:
            raise NameError(
                'wrong agent typ, it must be either agent (ag) or acces point (ap) ')

        if self.typ == 'ap':
            self.MoA = 1
        else:
            self.MoA = 0

        if 'mysql' in args['save']:
            config = ConfigParser.ConfigParser()
            config.read(pyu.getlong('simulnet.ini', 'ini'))
            sql_opt = dict(config.items('Mysql'))
            db = Database(sql_opt['host'], sql_opt['user'],
                          sql_opt['passwd'], sql_opt['dbname'])
            db.writenode(self.ID, self.name, self.MoA)

        if 'txt' in args['save']:
            pyu.writenode(self)
        if self.typ != 'ap' and args['loc']:

            self.loc = Localization(net=self.net, ID=self.ID,
                                    method=args['loc_method'])
            self.Ploc = PLocalization(loc=self.loc,
                                      loc_updt_time=args['loc_updt'],
                                      tx=self.tx,
                                      sim=args['sim'])

            self.sim.activate(self.Ploc, self.Ploc.run(), 1.5)

    def __repr__(self):
        s = 'General Agent info \n********************\n'
        s = s + 'name : ' + self.name + '\n'
        s = s + 'ID: ' + self.ID + '\n'
        s = s + 'typ: ' + self.typ

        s = s + '\n\n More Agent information about:'
        s = s + '\n+ Mecanichal => self.meca'

        s = s + '\n+ Network => self.net'
        s = s + '\n+ Personnal Network => self.PN'
        s = s + '\n+ Localization => self.loc\n\n'

        try:
            s = s + self.PN.__repr__() + '\n\n'
        except:
            s = s + 'No network simulated'

        if self.typ != 'ap':
            s = s + self.meca.__repr__() + '\n\n'
            try:
                s = s + self.loc.__repr__() + '\n\n'
            except:
                s = s + 'no localization simulated'

        return s
