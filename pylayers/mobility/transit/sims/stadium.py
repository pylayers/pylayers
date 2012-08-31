from SimPy.Simulation import *
from Transit.StateProcess import StateProcess, Updater, new_state
from Transit.PRTStation import PRTStationManager
from Transit.RailStation import RailStationManager
from Transit.Path import *
from Transit.PRT import PRT, Diverge, station, prt_draw
from Transit.Railcar import Railcar, rail_draw
from Transit.World import world
from Transit.Person import Person
from Transit.vec3 import vec3
from Transit.SteeringBehavior import FollowWaypoints, Separation, Containment, InterpenetrationConstraint, queue_steering_mind
from Transit.SteeringBehavior import Seek, Separation
from CommonUtils import PRTGenerator, PRTChecker
from random import random, uniform
import os

interval = 0.1

class RailcarGenerator(Process):
    def __init__(self, rx, rail_station):
        self.rx = rx
        self.rail_station = rail_station
        Process.__init__(self)

    def generate(self):
        while True:
            r1 = r = Railcar('railcar', interval)
            r.x = self.rx
            r.position = -120
            r.destination = 114
            activate(r, r.execute(), 0.0)
            r.station_manager = self.rail_station
            r2 = r = Railcar('railcar', interval)
            r.x = self.rx
            r.position = r1.position - r.length - 1
            r.destination = r1.destination - r.length - 1
            activate(r, r.execute(), 0.0)
            r.station_manager = self.rail_station
            while r1.position < 250:
                yield hold, self, 1
            self.cancel(r1)
            self.cancel(r2)

class PersonGenerator(Process):
    def __init__(self,
                 rail_station, r_x, r_y1, r_y2,
                 prt_station1, p1_x, p1_y1, p1_y2,
                 prt_station2, p2_x, p2_y1, p2_y2):
        Process.__init__(self)
        self.rail_station = rail_station
        self.r_x = r_x
        self.r_y1 = r_y1
        self.r_y2 = r_y2
        self.prt_station1 = prt_station1
        self.p1_x = p1_x
        self.p1_y1 = p1_y1
        self.p1_y2 = p1_y2
        self.prt_station2 = prt_station2
        self.p2_x = p2_x
        self.p2_y1 = p2_y1
        self.p2_y2 = p2_y2

    def generate(self):
        """
        creates a Person at the same time and same relative place on each
        of the PRT and LRT lines.
        """
        prt_platform_length = self.p1_y2 - self.p1_y1
        rail_platform_length = self.r_y2 - self.r_y1
        people = 0
        # http://findarticles.com/p/articles/mi_m1215/is_9_204/ai_108788436
        event_max = 5000
        while people < event_max:
            people += 1
            if (people % 100) == 0:
                print "people: ", people

            start = uniform(0, 100)

            if start < 33.33:
                p = Person(interval=interval)
                p.behaviors = [Seek(), Separation()]
                p.manager_args = []
                if start < 16.66:
                    p_start = self.p1_y1 + start * 6 * prt_platform_length / 100
                    p.position = vec3(self.p1_x, p_start)
                    platform_position = (p_start - self.p1_y1) * 0.7
                    p.destination = vec3(28., self.p1_y1 + 4 + platform_position)
                    p.manager = self.prt_station1.person_at_station
                else:
                    p_start = self.p2_y1 + (start - 16.66) * 6 * prt_platform_length / 100
                    p.position = vec3(self.p2_x, p_start)
                    platform_position = (p_start - self.p2_y1) * 0.7
                    p.destination = vec3(28., self.p2_y1 + 4 + platform_position)
                    p.manager = self.prt_station2.person_at_station
                activate(p, p.move(), 0.0)

            p = Person(interval=interval)
            p.behaviors = [Seek(), Separation()]
            p.position = vec3(self.r_x, self.r_y1 + 12 + start * 0.6 * rail_platform_length / 100)
            platform_position = start * (self.r_y2 - self.r_y1) / 100
            p.destination = vec3(78., self.r_y1 + platform_position)
            p.manager = self.rail_station.person_at_platform
            p.manager_args = []
            activate(p, p.move(), 0.0)

            yield hold, self, 0.25

        yield passivate, self

def main(mode='straight'):
    the_world = world()
    initialize()

    tk = the_world.tk
    canvas, x_, y_, s_ = tk.canvas, tk.x_, tk.y_, tk.s_

    ###
    ### PRT
    ###
    berths = 6
    psm1 = PRTStationManager(berths)
    psm2 = PRTStationManager(berths)
    activate(psm1, psm1.check(), 0.0)
    activate(psm2, psm2.check(), 0.0)

    s1 = Straight((20., -30.), (20., -5.0))
    s1.next, last_left, last_right, length \
        = station((20.0, -5.0), PRT.max_speed, 'right', berths, psm1)
    last_left.next = last_right.next = s2 \
        = Straight((20., -5. + length), (20., -5. + length + 5.))
    psm1.mainline = last_left
    s2.next, last_left, last_right, length \
        = station((20.0, -5.0 + length + 5), PRT.max_speed, 'right', berths, psm2)
    last_left.next = last_right.next \
        = Straight((20.0, 500.), (20., 10000.))
    psm2.mainline = last_left
    prt_draw(s1, 20.0, -31.0)

    pc = PRTChecker()
    activate(pc, pc.check(s2.members), 0.0)

    p1_ul, p1_lr = prt_platform = s1.next.platform
    prt_platform_center = (p1_ul.y + p1_lr.y) / 2
    psm1.platform = prt_platform
    p2_ul, p2_lr = s2.next.platform

    canvas.create_polygon( ((x_(p1_lr.x), y_(p1_ul.y)), (x_(p1_lr.x + 15), y_(p1_ul.y + 4)),
                            (x_(p1_lr.x + 15), y_(p1_lr.y - 4)), (x_(p1_lr.x), y_(p1_lr.y))),
                           fill='#ccc', outline='#ccc')

    canvas.create_polygon( ((x_(p2_lr.x), y_(p2_ul.y)), (x_(p2_lr.x + 15), y_(p2_ul.y + 4)),
                            (x_(p2_lr.x + 15), y_(p2_lr.y - 4)), (x_(p2_lr.x), y_(p2_lr.y))),
                           fill='#ccc', outline='#ccc')

    tg = PRTGenerator()
    # 75% empty to simulate sending extras to the stadium
    activate(tg, tg.generate(s1, interval, percent_empty=0.75), 0.0)

    ###
    ### Rail
    ###
    rx, ry = 80., 0.

    rsm = RailStationManager()
    activate(rsm, rsm.check(), 0.0)

    rail_stop_center = 100

    rail_draw(rx, ry, rail_stop_center)

    rg = RailcarGenerator(rx, rsm)
    rsm.mode = "stadium"
    activate(rg, rg.generate(), 20.0 + 30)

    rg.waypoint_high = vec3(rx - 4, p1_ul.y + 3)
    rg.waypoint_low = vec3(rx - 4, p1_lr.y - 3)

    canvas.create_polygon( ((x_(rx - 20), y_(rail_stop_center + 20)),
                            (x_(rx - 5), y_(rail_stop_center + Railcar.length)),
                            (x_(rx - 5), y_(rail_stop_center - Railcar.length)),
                            (x_(rx - 20), y_(rail_stop_center - 20))),
                           fill='#ccc', outline='#ccc')
        
    ###
    ### People
    ###
    pg = PersonGenerator(rsm, rx - 20, rail_stop_center - Railcar.length + 2,
                         rail_stop_center + Railcar.length - 2,
                         psm1, p1_lr.x + 15, p1_lr.y - 4, p1_ul.y + 4,
                         psm2, p2_lr.x + 15, p2_lr.y - 4, p2_ul.y + 4)
    activate(pg, pg.generate(), 45.0)

    ###
    ###
    ###
    u = Updater(interval)
    activate(u, u.execute(), 0.0)

    simulate(until=2000)

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3:
        main(sys.argv[1], int(sys.argv[2]))
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main()
