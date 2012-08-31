from SimPy.Simulation import *
from Transit.StateProcess import StateProcess, Updater, new_state
from Transit.PRTStation import PRTStationManager
from Transit.Path import *
from Transit.PRT import PRT, Diverge, station, prt_draw
from Transit.vec3 import vec3
from Transit.Person import Person
from Transit.SteeringBehavior import Seek, Separation
from CommonUtils import PRTGenerator
from random import random, uniform

interval = 0.1

class PersonGenerator(Process):
    def __init__(self, station1):
        Process.__init__(self)
        self.station1 = station1

    def generate(self):
        while True:
            start = uniform(5, 40)

            p = Person(interval=interval)
            p.behaviors = [Seek(), Separation()]
            p.position = vec3(100.0, start)
            platform_position = start * 8 / 199
            p.destination = self.destination1.copy()
            p.manager = self.station1.person_at_station
            p.manager_args = []
            activate(p, p.move(), 0.0)

            yield hold, self, 3.5

if __name__ == '__main__':
    initialize()

    berths = 10
    psm1 = PRTStationManager(berths)
    activate(psm1, psm1.check(), 0.0)
    pg = PersonGenerator(psm1)
    activate(pg, pg.generate(), 0.0)

    s1 = Straight((70., -30.), (70., 5.0))
    s1.next, last_left, last_right, length = station((70.0, 5.0), PRT.max_speed, 'right', berths, psm1)
    last_left.next = last_right.next = s2 = Straight((70.0, 5.0 + length), (70., 10000.0))
    psm1.mainline = last_left
    prt_draw(s1, 70.0, -31.0)
    pg.destination1 = vec3(s1.next.destinations[0].x + 8, (s1.next.platform[0].y + s1.next.platform[1].y) / 2)

    tg = PRTGenerator()
    activate(tg, tg.generate(s1, interval, percent_empty=0.60), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)

    simulate(until=1000)
    #simulate(until=400)
