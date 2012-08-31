from SimPy.Simulation import *
from Transit.Railcar import Railcar, rail_draw
from Transit.Person import Person
from Transit.PRT import PRT, Diverge, station, prt_draw
from Transit.Path import *
from Transit.StateProcess import StateProcess, Updater, new_state
from Transit.vec3 import vec3
from Transit.PRTStation import PRTStationManager
from Transit.RailStation import RailStationManager
from Transit.SteeringBehavior import Seek, Separation
from CommonUtils import PRTGenerator
from Tkinter import *
from random import random, uniform

interval = 0.1
prts_only = False

class RailcarGenerator(Process):
    def __init__(self, rail_station):
        Process.__init__(self)
        self.rail_station = rail_station

    def generate(self):
        while True:
            r1 = r = Railcar('railcar', interval)
            r.x = 30
            r.position = -120
            r.destination = 114
            r.station_manager = self.rail_station
            activate(r, r.execute(), 0.0)
            r2 = r = Railcar('railcar', interval)
            r.x = 30
            r.position = r1.position - r.length - 1
            r.destination = r1.destination - r.length - 1
            r.station_manager = self.rail_station
            activate(r, r.execute(), 0.0)
            yield hold, self, 100
            self.cancel(r1)
            self.cancel(r2)
            yield hold, self, 200

class PersonGenerator(Process):
    def __init__(self, station1, station2, rail_station):
        Process.__init__(self)
        self.station1 = station1
        self.station2 = station2
        self.rail_station = rail_station

    def generate(self):
        while True:
            start = uniform(1, 199)

            if not prts_only:
                p = Person(interval=interval)
                p.behaviors = [Seek(), Separation()]
                p.position = vec3(0.0, start)
                platform_position = start * 54 / 199 - 27
                p.destination = vec3(28, 100 + platform_position)
                p.manager = self.rail_station.person_at_platform
                activate(p, p.move(), 0.0)

            p = Person(interval=interval)
            p.behaviors = [Seek(), Separation()]
            p.position = vec3(100.0, start)
            platform_position = start * 8 / 199
            if start < 100:
                p.destination = self.destination1.copy()
                p.manager = self.station1.person_at_station
                p.manager_args = []
            else:
                p.destination = self.destination2.copy()
                p.manager = self.station2.person_at_station
                p.manager_args = []
            activate(p, p.move(), 0.0)

            yield hold, self, 7

if __name__ == '__main__':
    initialize()

    rail_draw(30., 0.)

    # PRT segments
    psm1 = PRTStationManager(3)
    activate(psm1, psm1.check(), 0.0)
    psm2 = PRTStationManager(3)
    activate(psm2, psm2.check(), 0.0)

    rsm = RailStationManager()
    rsm.destination_x = 30.0
    activate(rsm, rsm.check(), 0.0)


    pg = PersonGenerator(psm1, psm2, rsm)
    activate(pg, pg.generate(), 0.0)

    s1 = Straight((70., -30.), (70., 5.0))
    s1.next, last_left, last_right, length \
        = station((70.0, 5.0), PRT.max_speed, 'right', 3, psm1)
    last_left.next = last_right.next = s2 \
        = Straight((70.0, 5.0 + length), (70., 125.0))
    psm1.mainline = last_left
    s2.next, last_left, last_right, length \
        = station((70.0, 125.0), PRT.max_speed, 'right', 3, psm2)
    last_left.next = last_right.next = s7 = Straight((70.0, 125.0 + length), (70., 10000.))
    psm2.mainline = last_left
    prt_draw(s1, 70.0, -31.0)
    pg.destination1 = vec3(s1.next.destinations[0].x + 2.2, s1.next.destinations[0].y)
    pg.destination2 = vec3(s2.next.destinations[0].x + 2.2, s2.next.destinations[0].y)

    rg = RailcarGenerator(rsm)
    activate(rg, rg.generate(), 0.0)

    tg = PRTGenerator()
    activate(tg, tg.generate(s1, interval), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)

    simulate(until=1000)
    #simulate(until=400)
