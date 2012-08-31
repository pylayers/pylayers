from SimPy.Simulation import *
from Transit.Person import Person
from Transit.StateProcess import Updater
from Transit.World import world
from Transit.vec3 import vec3
from math import *
from random import uniform
from Transit.SteeringBehavior import Wander, Separation, Containment

class Restarter(Process):
    def __init__(self, people):
        Process.__init__(self)
        self.people = people

    def execute(self):
        tk = world().tk
        while True:
            for person in self.people:
                for ii in 0, 1:
                    if person.position[ii] > 25:
                        person.position[ii] = 0
                    elif person.position[ii] < 0:
                        person.position[ii] = 25
            yield hold, self, 0.1

def main(cluster_size=20):
    interval = 0.1
    the_world = world(width=25, height=25, scale=25)
    initialize()

    xx = 12.5

    p = Person(interval)
    p.behaviors = [Wander(), Separation(), Containment()]
    p.position = vec3(12.5, 12.5)
    p.destination = vec3(-100, 100)
    activate(p, p.move(), 0.0)
    people = [p]

    tk = the_world.tk
    canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
    canvas.create_rectangle(x_(-1), y_(-1), x_(100), y_(100), fill='moccasin')

    for wall in ( ( (-100, 14.5), (10.5, 14.5), (10.5, 100) ),
                  ( (10.5, -100), (10.5, 8.5), (8.5, 10.5), (-100, 10.5) ),
                  ( (14.5, 100), (14.5, 14.5), ( 100, 14.5) ),
                  ( ( 100, 10.5), (16.5, 10.5), (14.5, 8.5), (14.5, -100) ),
                  ):
        points = []
        for point in wall:
            points.append(x_(point[0]))
            points.append(y_(point[1]))
        canvas.create_polygon(points, fill='gray', outline='black')
        for ii in range(0, len(wall)-1):
            the_world.add_wall(wall[ii], wall[ii+1])

    restarter = Restarter(people)
    activate(restarter, restarter.execute(), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    simulate(until=10000)

main()

