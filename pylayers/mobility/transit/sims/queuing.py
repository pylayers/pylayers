from SimPy.Simulation import *
from Transit.Person import Person
from Transit.StateProcess import Updater
from Transit.World import world
from Transit.vec3 import vec3
from Transit.SteeringBehavior import Queuing, FollowWaypoints, Separation, Containment, queue_steering_mind
from math import *
from random import uniform

def random_people(cx, cy, waypoints, cluster_size, interval, people):
    def random_position():
        range = cluster_size * (pi * Person.average_radius**2 + 0.2)
        offset = range / 2
        return uniform(0.0, range) - offset
    for ii in range(0, cluster_size):
        p = Person(interval)
        p.behaviors = [Containment(), Separation(), Queuing(), FollowWaypoints()]
        p.steering_mind = queue_steering_mind
        try_again = True
        while try_again:
            px, py = cx + random_position(),  cy + random_position()
            if py > (21.0 - Person.radius):
                continue
            try_again = False
            for person in people:
                jx, jy, unused_z = person.position
                if sqrt((jx-px)**2 + (jy-py)**2) < Person.average_radius*2.5:
                    try_again = True
        p.position = vec3(px, py)
        p.waypoints = waypoints[0:]
        activate(p, p.move(), 0.0)
        people.append(p)
    return people

class Restarter(Process):
    def __init__(self, people, waypoints):
        Process.__init__(self)
        self.people = people
        self.waypoints = waypoints

    def execute(self):
        waypoints = self.waypoints
        tk = world().tk
        while True:
            for person in self.people:
                if person.position[1] > 25:
                    start_angle = uniform(0, pi)
                    xx, yy = cos(start_angle) * 25 + 12.5, 21.0 - person.radius * 1.2 - sin(start_angle) * 25
                    person.position = vec3(xx, yy)
                    person.waypoints = waypoints[0:]
            yield hold, self, 0.1

def main(cluster_size=30):
    interval = 0.1
    the_world = world(width=25, height=25, scale=25)
    initialize()

    xx = 12.5
    door_width = 4
    corner_size = 1
    people = []

    waypoints = [vec3(xx, 21.0), vec3(xx, 22.0), vec3(xx, 23.0),
                 vec3(xx, 40.0)]

    people = random_people(xx, 4.0, waypoints, cluster_size, interval, people)

    tk = the_world.tk
    canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
    canvas.create_rectangle(x_(-1), y_(-1), x_(100), y_(100), fill='moccasin')

    for wall in ( ( (-100, 21.0), (12.5 - door_width / 2 - corner_size, 21.0),
                    (12.5 - door_width / 2, 21.0 + corner_size), (12.5 - door_width / 2, 100) ),
                  ( (12.5 + door_width / 2, 100), (12.5 + door_width / 2, 21.0 + corner_size),
                    (12.5 + door_width / 2 + corner_size, 21.0), (100, 21.0) ), ):
        points = []
        for point in wall:
            points.append(x_(point[0]))
            points.append(y_(point[1]))
        canvas.create_polygon(points, fill='gray', outline='black')
        for ii in range(0, len(wall)-1):
            the_world.add_wall(wall[ii], wall[ii+1])

    restarter = Restarter(people, waypoints)
    activate(restarter, restarter.execute(), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    simulate(until=10000)

main()
