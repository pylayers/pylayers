from SimPy.SimulationRT import Process, Simulation
from pylayers.mobility.transit.Person import Person
from pylayers.mobility.transit.vec3 import vec3
from pylayers.mobility.transit.World import world
from pylayers.mobility.transit.SteeringBehavior import Seek, Wander, Queuing, FollowWaypoints, Separation, Containment, InterpenetrationConstraint, queue_steering_mind, default_steering_mind


from random import normalvariate, uniform, gauss
import numpy as np
import networkx as nx
import time
import ConfigParser
import pylayers.util.pyutil as pyu
from pylayers.network.network import Node, Network
from pylayers.network.communication import Gcom, TX, RX
from pylayers.location.localization import Localization, PLocalization
from pylayers.gis.layout import Layout
from pylayers.util.utilnet import *
#from pylayers.util.pymysqldb import Database
import pdb


class Agent(object):
    def __init__(self, **args):
        """ Mobile Agent Init

           Parameters
           ----------
           'ID': string
                agent ID
           'name': string
                Agent name
           'type': string
                agent type . 'ag' for moving agent, 'ap' for static acces point
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
           'Layout': pylayers.gis.Layout()
           'net':pylayers.network.Network(),
           'RAT': list of string
                list of used radio access techology of the agent
           'world': transit.world()
                Soon deprecated
           'save': list of string
                list of save method ( soon deprecated)
           'sim':Simpy.SimulationRT.Simulation(),
           'epwr': dictionnary
                dictionnary of emmited power of transsmitter{'rat#':epwr value}
           'sens': dictionnary
                dictionnary of sensitivity of reveicer {'rat#':sens value}
           'dcond': dictionnary
                Not used yet
           'gcom':pylayers.communication.Gcom()
                Communication graph
           'comm_mod': string
                Communication between nodes mode: 
                'autonomous': all TOAs are refreshed regulary
                'synchro' : only visilbe TOAs are refreshed 
        """
        defaults = {'ID': 0,
                    'name': 'johndoe',
                    'type': 'ag',
                    'pos': np.array([]),
                    'roomId': 0,
                    'froom':[],
                    'wait':[],
                    'cdest':'random',
                    'meca_updt': 0.1,
                    'loc': False,
                    'loc_updt': 0.5,
                    'loc_method': ['geo'],
                    'Layout': Layout(),
                    'net': Network(),
                    'RAT': ['wifi'],
                    'world': world(),
                    'save': [],
                    'sim': Simulation(),
                    'epwr':{},
                    'sens': {},
                    'dcond': {},
                    'gcom': Gcom(),
                    'comm_mode':'autonomous'}

        for key, value in defaults.items():
            if key not in args:
                args[key] = value

        self.args = args
        self.ID = args['ID']
        self.name = args['name']
        self.type = args['type']
        # Create Network
        self.net = args['net']
        self.epwr = args['epwr']
        self.gcom = args['gcom']
        try:
            self.dcond = args['dcond']
        except:
            pass

        if self.type == 'ag':
            # mechanical init
            self.meca = Person(ID=self.ID,
                                roomId=args['roomId'],
                                L=args['Layout'],
                                net=self.net,
                                interval=args['meca_updt'],
                                wld=args['world'],
                                sim=args['sim'],
                                moving=True,
                                froom=args['froom'],
                                wait=args['wait'],
                                cdest=args['cdest'],
                                save=args['save'])
            self.meca.behaviors = [Queuing(),Seek(), Containment(),\
                                   Separation(), InterpenetrationConstraint()]
            self.meca.steering_mind = queue_steering_mind
#            self.meca.steering_mind = queue_steering_mind
        # filll in network

            ## Network init
            self.node = Node(ID=self.ID, p=conv_vecarr(self.meca.position),
                             t=time.time(), RAT=args['RAT'],
                             epwr=args['epwr'], sens=args['sens'], type=self.type)
            self.net.add_nodes_from(self.node.nodes(data=True))
            self.sim = args['sim']
            self.sim.activate(self.meca, self.meca.move(), 0.0)
            self.PN = self.net.node[self.ID]['PN']

            ## Communication init

            if args['comm_mode'] == 'synchro':
                ## The TOA requests are made every refreshTOA time ( can be modified in agent.ini)

                self.rxr = RX(net=self.net, ID=self.ID,
                              dcond=self.dcond, gcom=self.gcom, sim=self.sim)
                self.rxt = RX(net=self.net, ID=self.ID,
                              dcond=self.dcond, gcom=self.gcom, sim=self.sim)

                self.sim.activate(self.rxr, self.rxr.refresh_RSS(), 0.0)
                self.sim.activate(self.rxt, self.rxt.refresh_TOA(), 0.0)


            elif args['comm_mode'] == 'autonomous':
                ## The TOA requests are made by node only when they are in visibility of pairs.

                self.rxr = RX(net=self.net, ID=self.ID,
                              dcond=self.dcond, gcom=self.gcom, sim=self.sim)
                self.rxt = RX(net=self.net, ID=self.ID,
                              dcond=self.dcond, gcom=self.gcom, sim=self.sim)
                self.txt = TX(net=self.net, ID=self.ID,
                              dcond=self.dcond, gcom=self.gcom, sim=self.sim)

                self.sim.activate(self.rxr, self.rxr.refresh_RSS(), 0.0)
                self.sim.activate(self.rxt, self.rxt.wait_TOArq(), 0.0)
                self.sim.activate(self.txt, self.txt.request_TOA(), 0.0)



        elif self.type == 'ap':
            if args['roomId'] == -1:
                self.node = Node(ID=self.ID, p=self.args['pos'],
                                 t=time.time(), RAT=args['RAT'],
                                 epwr=args['epwr'], sens=args['sens'], type=self.type)
            else:
                pp = np.array(args['Layout'].Gr.pos[self.args['roomId']])
                self.node = Node(ID=self.ID, p=pp, t=time.time(), RAT=args['RAT'],
                                 epwr=args['epwr'], sens=args['sens'], type=self.type)
            self.net.add_nodes_from(self.node.nodes(data=True))
            self.sim = args['sim']
#            self.sim.activate(self.meca, self.meca.move(),0.0)
            self.PN = self.net.node[self.ID]['PN']
            self.PN.node[self.ID]['pe'] = self.net.node[self.ID]['p']

        else:
            raise NameError('wrong agent type, it must be either agent (ag) or acces point (ap) ')

        if self.type == 'ap':
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
        if args['loc'] and self.type != 'ap':

            self.loc = Localization(net=self.net, ID=self.ID,
                                    method=args['loc_method'])
            self.Ploc = PLocalization(loc=self.loc,
                                      loc_updt_time=args['loc_updt'], sim=args['sim'])
            self.sim.activate(self.Ploc, self.Ploc.run(), 1.5)

