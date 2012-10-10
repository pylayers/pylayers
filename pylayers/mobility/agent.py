from SimPy.Simulation import Process,Simulation
from pylayers.mobility.transit.Person3 import Person3
from pylayers.mobility.transit.vec3 import vec3
from pylayers.mobility.transit.World import world
from pylayers.mobility.transit.SteeringBehavior2 import Seek, Wander, Queuing, FollowWaypoints, Separation, Containment, InterpenetrationConstraint, queue_steering_mind, default_steering_mind


from random import normalvariate,uniform,gauss
import numpy as np
import networkx as nx
import time
import ConfigParser
import pylayers.util.pyutil as pyu
from pylayers.network.network import  Node,Network
from pylayers.location.localization import Localization,PLocalization
from pylayers.gis.layout import Layout
from pylayers.util.utilnet import *
from pylayers.util.pymysqldb import Database 
import pdb


class Agent(object):
    def __init__(self,**args):
        defaults = {'ID': 0,'name': 'johndoe','type':'ag','pos':np.array([]),'roomId':0, 'meca_updt':0.1,'loc':False,'loc_updt':0.5,'Layout':Layout(),'net':Network(),'RAT':['wifi'],'world':world(),'save':False, 'sim':Simulation}

        for key, value in defaults.items():
            if not args.has_key(key):
                args[key]=value  

        self.args = args
        self.ID = args['ID']
        self.name = args['name']
        self.type=args['type']
        # Create Network
        self.net=args['net']
        # mecanique
        if self.type == 'ag':
            self.meca=Person3( ID=self.ID,
                               roomId=args['roomId'],
                               L=args['Layout'],
                               net=self.net,
                               interval=args['meca_updt'],
                               wld=args['world'],
                               sim=args['sim'],
                               moving=True,
                               save=args['save'])
            self.meca.behaviors  = [Seek(),Containment(),Separation(),InterpenetrationConstraint()]
            self.meca.steering_mind = queue_steering_mind
#            self.meca.steering_mind = queue_steering_mind
        # filll in network
            self.node = Node(ID=self.ID,p=conv_vecarr(self.meca.position),t=time.time(),RAT=args['RAT'],type=self.type)
            self.net.add_nodes_from(self.node.nodes(data=True))
            self.sim=args['sim']
            self.sim.activate(self.meca, self.meca.move(),0.0)
            self.PN=self.net.node[self.ID]['PN']

        elif self.type== 'ap':
#            self.meca=Person3(ID=self.ID,roomId=args['roomId'],L=args['Layout'],net=self.net,interval=args['meca_updt'],sim=args['sim'],moving=False)
#            self.meca.behaviors  = []
            if args['roomId'] == -1:
                self.node = Node(ID=self.ID,p=self.args['pos'],t=time.time(),RAT=args['RAT'],type=self.type)
            else:
                pp = np.array(args['Layout'].Gr.pos[self.args['roomId']])
                self.node = Node(ID=self.ID,p=pp,t=time.time(),RAT=args['RAT'],type=self.type)
            self.net.add_nodes_from(self.node.nodes(data=True))
            self.sim=args['sim']
#            self.sim.activate(self.meca, self.meca.move(),0.0)
            self.PN=self.net.node[self.ID]['PN']
            self.PN[self.ID]['pe']=self.net.node[self.ID]['p']

        else :
            raise NameError('wrong agent type, it must be either agent (ag) or acces point (ap) ')


        if self.type == 'ap':
            self.MoA=1
        else :
            self.MoA=0

        if 'mysql' in args['save']:
            config = ConfigParser.ConfigParser()
            config.read(pyu.getlong('simulnet.ini','ini'))
            sql_opt = dict(config.items('Mysql'))
            db = Database(sql_opt['host'],sql_opt['user'],sql_opt['passwd'],sql_opt['dbname'])
            db.writenode(self.ID,self.name,self.MoA)

        if 'txt' in args['save']:
            pyu.writenode(self)
        if args['loc'] and self.type != 'ap':
            self.loc=Localization(net=self.net,ID=self.ID)
            self.Ploc = PLocalization(loc=self.loc,loc_updt_time=args['loc_updt'],sim=args['sim'])
            self.sim.activate(self.Ploc,self.Ploc.run(),1.0)


