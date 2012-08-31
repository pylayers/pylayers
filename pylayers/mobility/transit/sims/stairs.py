import SimPy.Simulation
from SimPy.Simulation import *
from Transit.Person import Person
from Transit.StateProcess import Updater
from Transit.World import world
from Transit.vec3 import vec3
from Transit.SteeringBehavior import Seek, Queuing, FollowWaypoints, Separation, Containment, InterpenetrationConstraint, queue_steering_mind
from Transit.Zones import Stairs
from math import *
from random import uniform

class PersonGenerator(Process):
    def generate(self, interval, waypoints):
        people = []
        rate = 0.5
        #while True:
	for k in range(10):
            start = uniform(0, 5)

            p = Person(interval=interval)
	    #p.destination= vec3(12.5,uniform(2,20))
            p.behaviors = [Containment(), Separation(), Queuing(), FollowWaypoints(), InterpenetrationConstraint()]
 	    #p.behaviors = [FollowWaypoints(),InterpenetrationConstraint()]
 	    #p.behaviors = [Seek(),InterpenetrationConstraint()]
            p.steering_mind = queue_steering_mind
            p.position = vec3(start, 0.0)
            p.waypoints = waypoints[0:]
            activate(p, p.move(), 0.0)
            people.append(p)

            count = 0
            for p in people:
                if p.position[1] > 25:
                    self.cancel(p)
                    p.world.remove_boid(p)
                    people.remove(p)
                if p.position[1] < 7:
                    count += 1

            if count > 30:
                rate = 1.5
            elif count < 20:
                rate = 0.5

            yield hold, self, rate

def main(mode='path'):
    interval = 0.1
    the_world = world(width=25, height=25, scale=25)
    initialize()

    xx = 12.5
    door_width = 3.0
    corner_size = 1

    waypoints = [vec3(xx, 4.0), vec3(xx, 5.0), vec3(xx, 6.0),
                 vec3(xx, 40.0)]

    tk = the_world.tk
    canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
    canvas.create_rectangle(x_(-1), y_(-1), x_(100), y_(100), fill='blue')

    # stairs are added here so the graphics appear below the walls
    if mode == 'stairs':
        steps = 12
        step_length = 0.28
        zone = Stairs(lower_left=vec3(xx - door_width / 2, 7),
                      upper_right=vec3(xx + door_width / 2, 7 + steps * step_length))
        the_world.add_zone(zone)
        zone.draw()
    wall1 = ( (-100, 5.0), (12.5 - door_width / 2 - corner_size, 5.0),
                    (12.5 - door_width / 2, 5.0 + corner_size), (12.5 - door_width / 2, 100) )
    wall2 = ( (12.5 + door_width / 2, 100), (12.5 + door_width / 2, 5.0 + corner_size),
                    (12.5 + door_width / 2 + corner_size, 5.0), (100, 5.0) )
    wall3 = ( (12,10),(13,10),(13,11),(12,11) )
    for wall in ( wall1 , wall2 ,wall3 ):
        points = []
        for point in wall:
            points.append(x_(point[0]))
            points.append(y_(point[1]))
        canvas.create_polygon(points, fill='red', outline='black')
        for ii in range(0, len(wall)-1):
            the_world.add_wall(wall[ii], wall[ii+1])

    pg = PersonGenerator()
    activate(pg, pg.generate(interval, waypoints), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    SimPy.Simulation.simulate(until=100)

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3:
        main(sys.argv[1], int(sys.argv[2]))
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main('stairs')
