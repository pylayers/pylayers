from SimPy.Simulation import *
from math import *
from random import normalvariate,uniform
from pylayers.mobility.transit.vec3 import vec3
from pylayers.mobility.transit.World import world
from pylayers.transit.SteeringBehavior2 import default_steering_mind
from random import uniform,gauss
import numpy as np

# "near collision avoidance" inspired from
#   http://people.revoledu.com/kardi/publication/Kouchi2001.pdf

def truncate(self, max):
    if self.length() > max:
        return self.normalize() * max
    else:
        return vec3(self)
vec3.truncate = truncate

def scale(self, size):
    return self.normalize() * size
vec3.scale = scale

def copy(self):
    return vec3(self)
vec3.copy = copy

class Person2(Process):
    """ Person2 
    """
    max_acceleration = 2.0 # m/s/s
    max_speed        = 1.2 # m/s
    #radius = 0.2106  # if one person takes 1.5 feet^2 of space, per traffic stds
    # 2r = 0.5 to 0.7 for "sports fans", per the Helbing, Farkas, Vicsek paper
    radius           = 2.85   # per the Teknomo, et al, paper 
    mass             = 80 # kg
    average_radius   = 0.6
    npers            = 0
#GeomNet          = np.array((0,0,[[1,2,3]],[[1,0,0]],[[0,0,1]]),dtype=GeomNetType)
    def __init__(self, interval=0.5,roomId=0, L=[],sim=None):
        """
        boid is initialized in a room of Layout L 
        """                
        #GeomNetType = np.dtype([('Id',int), 
        #                ('time',int), 
        #		('p',float,(1,3)),
        #		('v',float,(1,3)),
        #		('a',float,(1,3))])
        Person2.npers +=1
        Process.__init__(self,sim=sim)
        self.L = L 
        self.world = world()
        self.interval = interval
        self.manager = None
        self.manager_args = []
        self.waypoints = []
        self.roomId    = roomId
        self.nextroomId   = int(np.floor(uniform(0,self.L.Gr.size())))
        while self.nextroomId == self.roomId : # test destination different de l'arrive
            self.nextroomId   = int(np.floor(uniform(0,self.L.Gr.size())))
            print "nextroom : ", self.nextroomId
            self.wp           =  self.L.waypoint(roomId,self.nextroomId)
            for tup in self.wp[1:]:
                self.waypoints.append(vec3(tup)) 
            try:
                self.position = vec3(L.Gr.pos[roomId])
            except:         
                self.position = vec3()
#		self.old_pos = vec3()
		    self.stuck = 0               
            self.destination = self.waypoints[0]
            self.velocity = vec3()
            self.localx = vec3(1, 0)
            self.localy = vec3(0, 1)
            self.world.add_boid(self)
               
                #        GeomNet = np.vstack((GeomNet,np.array((Nnodes,0,[list(self.position)],[list(self.velocity)],[list(self.velocity)]),dtype=GeomNetType)))

                # from Helbing, et al "Self-organizing pedestrian movement"
                self.max_speed = normalvariate(1.36, 0.26)
                self.desired_speed = self.max_speed
                self.radius = normalvariate(self.average_radius, 0.025) / 2
                self.intersection = vec3()
                self.arrived = False 
                self.endpoint = False 
                self.behaviors = []
                self.steering_mind = default_steering_mind
                self.cancelled = 0


        def move(self):
                """
                Gr : Graph of rooms  
                """
                while True:


                        while self.cancelled:
                                yield passivate, self
                                print "Person.move: activated after being cancelled" 
                        checked = []
                        for zone in self.world.zones(self):
                                if zone not in checked:
                                        checked.append(zone)
                                        zone(self)
                        acceleration = self.steering_mind(self) 
                        acceleration = acceleration.truncate(self.max_acceleration)
                        self.acceleration = acceleration 
                        velocity = self.velocity + acceleration * self.interval
                        self.velocity = velocity.truncate(self.max_speed) 
                        if velocity.length() > 0.2:
                                # record direction only when we've really had some
                                self.localy = velocity.normalize()
                                self.localx = vec3(self.localy.y, -self.localy.x)

                        self.position = self.position + self.velocity * self.interval
                        self.update()
                        self.world.update_boid(self)


                        if self.arrived:

                                self.arrived = False
                                if self.endpoint:
                                        self.endpoint=False
                                        self.roomId = self.nextroomId
                                #
                                # Si porte on continue
                                #
                                #
                                # ig destination --> next room
                                #
                                #adjroom  = self.L.Gr.neighbors(self.roomId)
                                #Nadjroom = len(adjroom)
				        self.nextroomId   = int(np.floor(uniform(0,self.L.Gr.size())))
					while (self.nextroomId == self.roomId) or (self.nextroomId in self.forbid) : # test destination different de l'arrive
	                                        self.nextroomId   = int(np.floor(uniform(0,self.L.Gr.size())))
				        wp                =  self.L.waypoint(self.roomId,self.nextroomId)
                                        for tup in wp[1:]:
                                                self.waypoints.append(vec3(tup)) 
                                #nextroom = adjroom[k]
                                        print "room : ",self.roomId
                                        print "nextroom : ",self.nextroomId
                                #p_nextroom = self.L.Gr.pos[self.nextroomId]
                                #setdoors1  = self.L.Gr.node[self.roomId]['doors']
                                #setdoors2  = self.L.Gr.node[nextroom]['doors']
                                #doorId     = np.intersect1d(setdoors1,setdoors2)[0]
                                #
                                # coord door
                                #
                                #unode = self.L.Gs.neighbors(doorId)	
                                #p1    = self.L.Gs.pos[unode[0]]
                                #p2    = self.L.Gs.pos[unode[1]]
                                #print p1
                                #print p2
                                #pdoor = (np.array(p1)+np.array(p2))/2
                                        self.destination = self.waypoints[0]
                                #waittime = uniform(0,10)

                                #if self.manager:
                                #    if self.manager(self, *self.manager_args):
                                #        yield hold , self , waittime
                                #else:
                                #    yield hold, self , waittime 
			                yield hold, self, 50*abs(gauss(5,2.5))

                                else:	
                                        del self.waypoints[0]
                                #print "wp : ", self.waypoints
                                        if len(self.waypoints)==1:
                                                self.endpoint=True
                                        self.destination = self.waypoints[0]
                                #print "dest : ", self.destination        
                        else:
                                yield hold, self, self.interval

        def delete(self):
                tk = self.world.tk
                tk.canvas.delete(self.graphic)
                self.world.remove_boid(self)
                self.cancelled = 1

        def update(self):
                tk = self.world.tk
                if tk is None: return
                if not hasattr(self, 'graphic'):
                    self.graphic = tk.canvas.create_polygon(0, 0, 0, 1, 1, 1, 1, 0, smooth=True, fill='red', outline='black', tag='person')
                    if tk.collision_vectors:
                        self.graphic_front_collision = tk.canvas.create_line(0, 0, 1, 0, fill='white')
                        self.graphic_left_collision = tk.canvas.create_line(0, 0, 1, 0, fill='white')
                        self.graphic_right_collision = tk.canvas.create_line(0, 0, 1, 0, fill='white')
                        self.graphic_intersection = tk.canvas.create_oval(0, 0, 1, 1, fill='red')
                        self.graphic_intersection_normal = tk.canvas.create_line(0, 0, 1, 1, fill='red')
                    if tk.vectors:
                        self.graphic_velocity = tk.canvas.create_line(0, 0, 1, 0, fill='green')
                        self.graphic_acceleration = tk.canvas.create_line(0, 0, 1, 0, fill='blue')
                x_, y_ = tk.x_, tk.y_
                position = self.position
                width, depth = self.localx * (self.radius), self.localy * (self.radius * 0.65)
                ul = position - depth + width
                ur = position + depth + width
                lr = position + depth - width
                ll = position - depth - width
                tk.canvas.coords(self.graphic,
                                 x_(ul.x), y_(ul.y), x_(ur.x), y_(ur.y),
                                 x_(lr.x), y_(lr.y), x_(ll.x), y_(ll.y))
                tk.canvas.create_line(x_(ur.x),y_(ur.y),x_(ul.x),y_(ul.y))
        #        if tk.vectors:
        #            xx, yy, unused_z = position
        #            velocity = self.velocity
        #            tk.canvas.coords(self.graphic_velocity,
        #                             x_(xx), y_(yy), x_(xx + velocity.x), y_(yy + velocity.y))
        #            acceleration = self.acceleration
        #            tk.canvas.coords(self.graphic_acceleration,
        #                             x_(xx), y_(yy), x_(xx + acceleration.x), y_(yy + acceleration.y))
        #        if tk.collision_vectors:
        #            xx, yy, unused_z = position
        #            speed = self.velocity.length()
        #            front_check = 0.5 + speed * 1.5
        #            side_check = 0.5 + speed * 0.5
        #            point = self.localy.scale(front_check)
        #            tk.canvas.coords(self.graphic_front_collision,
        #                             x_(xx), y_(yy), x_(xx + point.x), y_(yy + point.y))
        #            point = (self.localy + self.localx).scale(side_check)
        #            tk.canvas.coords(self.graphic_left_collision,
        #                             x_(xx), y_(yy), x_(xx + point.x), y_(yy + point.y))
        #            point = (self.localy - self.localx).scale(side_check)
        #            tk.canvas.coords(self.graphic_right_collision,
        #                             x_(xx), y_(yy), x_(xx + point.x), y_(yy + point.y))
        #            if self.intersection:
        #                xx, yy, unused_z = self.intersection
        #                tk.canvas.coords(self.graphic_intersection,
        #                                 x_(xx-0.2), y_(yy-0.2), x_(xx + 0.2), y_(yy + 0.2))
        #                point = self.intersection_normal
        #                tk.canvas.coords(self.graphic_intersection_normal,
        #                                 x_(xx), y_(yy), x_(xx + point.x), y_(yy + point.y))
        #            else:
        #                tk.canvas.coords(self.graphic_intersection,
        #                                 0, 0, 0, 0)
        #                tk.canvas.coords(self.graphic_intersection_normal,
        #                                 0, 0, 0, 0)
                tk.canvas.lower(self.graphic, 'vehicle')
