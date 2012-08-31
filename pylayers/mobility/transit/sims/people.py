"""
There are four example modes here: solo, immovable object, clump, and two clumps
The mode is the first command line argument, ie:

    people.py 'immovable object'

The clump modes take an additional argument which is the number of
people in each clump, ie:

    people.py 'two clumps' 20

"""

from SimPy.Simulation import *
from Transit.Person import Person
from Transit.StateProcess import Updater
from Transit.World import world
from Transit.vec3 import vec3
from Transit.SteeringBehavior import Arrive, Separation
from math import *
from random import random

def random_people(cx, cy, dstx, dsty, cluster_size, interval, people):
    def random_position():
        range = cluster_size * (pi * Person.radius**2 + 0.2)
        offset = range / 2
        return random() * range - offset
    for ii in range(0, cluster_size):
        p = Person(interval)
        p.behaviors = [Arrive(), Separation()]
        try_again = True
        while try_again:
            px, py = cx + random_position(),  cy + random_position()
            try_again = False
            for person in people:
                jx, jy, unused_z = person.position
                if sqrt((jx-px)**2 + (jy-py)**2) < Person.radius*2.5:
                    try_again = True
        p.position = vec3(px, py)
        p.destination = vec3(dstx, dsty)
        activate(p, p.move(), 0.0)
        people.append(p)

def main(testing='two_clumps', cluster_size=20):
    interval = 0.1
    world(width=25, height=25, scale=25)
    initialize()

    xx = 12.5
    people = []

    p = Person(interval)
    p.behaviors = [Arrive(), Separation()]
    p.position = vec3(xx, 4.0)
    p.destination = vec3(xx, 21.0)
    people.append(p)
    activate(p, p.move(), 0.0)

    if testing == 'solo':
        pass

    if testing == 'immovable object':
        p = Person(interval)
        p.behaviors = [Arrive(), Separation()]
        p.position = vec3(xx - 0.08, 14.0)
        p.destination = vec3(xx - 0.08, 14.0)
        people.append(p)
        p.update()
        world().update_boid(p)
        # DON'T activate(p, p.move(), 0.0)

    if testing in ('clump', 'two clumps'):
        random_people(xx, 21.0, xx, 4.0, cluster_size, interval, people)

    if testing == 'two clumps':
        random_people(xx, 4.0, xx, 21.0, cluster_size, interval, people)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    simulate(until=240)

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3:
        main(sys.argv[1], int(sys.argv[2]))
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main()

