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
from CommonUtils import PRTGenerator, PRTChecker, Clock
from random import random, uniform
import os

interval = 0.1

class RailcarGenerator(Process):
    def __init__(self, mode, platform_center_y, doorway_width, prt_station, rail_x, headway):
        Process.__init__(self)
        self.world = world()
        self.mode = mode
        self.platform_center_y = platform_center_y
        self.doorway_width = doorway_width
        self.prt_station = prt_station
        self.rail_x = rail_x
        self.headway = headway
        self.rail_waypoints = []
        self.rail_position = []
        self.rail_destination = []
        self.rows = 3
        self.columns = 2
        self.compactness = 1.2
        if self.rows > 6 or self.columns > 4:
            self.compactness = 0.8
        rx = self.rail_x
        railside_x = rx - Railcar.width / 2
        outside_x = rx + Railcar.width / 2
        top_rail_stop = platform_center_y + Railcar.length / 2 + 0.5
        for kk in (0, 1):
            position = -140 - (Railcar.length + 1) * kk
            destination = top_rail_stop - (Railcar.length + 1) * kk
            self.rail_position.append(position)
            self.rail_destination.append(destination)
            top = destination + Railcar.length / 2
            bottom = destination - Railcar.length / 2
            self.world.add_wall( (railside_x, top), (outside_x, top),
                                 (outside_x, bottom), (railside_x, bottom) )
            rail_waypoints = []
            for ii in range(0, len(Railcar.doors)):
                rail_waypoints.append( [vec3(railside_x, top - Railcar.doors[ii]),
                                        vec3(rx - 3.5, top - Railcar.doors[ii])] )
            self.rail_waypoints.append(rail_waypoints)
        print "passegers transferring:", self.rows * self.columns * len(Railcar.doors) * 2

    def generate(self):
        sm = RailStationManager()
        activate(sm, sm.check(), 0.0)
        rows = self.rows
        compactness = self.compactness
        while True:
            cars = []
            for kk in (0, 1):
                r = Railcar('railcar', interval)
                cars.append(r)
                r.x = self.rail_x
                r.position = self.rail_position[kk]
                r.destination = self.rail_destination[kk]
                r.station_manager = sm
                for ii in range(0, len(Railcar.doors)):
                    for jj in range(0, rows):
                        door = self.rail_destination[kk] + Railcar.length/2 - Railcar.doors[ii]
                        yy = door - (jj - (rows - 1) / 2) * (Person.average_radius * 2.0 * compactness)
                        self.person(r, yy, self.rail_waypoints[kk][ii])
                activate(r, r.execute(), 0.0)

            yield hold, self, 100
            self.cancel(cars[0])
            self.cancel(cars[1])
            yield hold, self, self.headway - 100

    def person(self, railcar, yy, railcar_waypoints):
        columns = self.columns
        compactness = self.compactness
        for ii in range(0, columns):
            xx = self.rail_x + (ii - (float(columns) - 1.0) / 2.0) * (Person.average_radius * 2.0 * compactness)
            p = Person(interval=interval)
            p.behaviors = [Containment(), Separation(), FollowWaypoints(), InterpenetrationConstraint()]
            p.steering_mind = queue_steering_mind
            p.position = vec3(xx, yy)
            p.waypoints = railcar_waypoints[0:]
            if self.mode == 'straight':
                if p.waypoints[0].y - self.waypoint_high.y > 3:
                    p.waypoints.append(self.waypoint_high)
                if self.waypoint_low.y - p.waypoints[0].y > 3:
                    p.waypoints.append(self.waypoint_low)
            elif self.mode == 'doorway':
                if p.waypoints[0].y > self.platform_center_y:
                    p.waypoints.append(self.waypoint_high)
                if p.waypoints[0].y <= self.platform_center_y:
                    p.waypoints.append(self.waypoint_low)
                p.waypoints.append(vec3(p.waypoints[-1]))
                p.waypoints[-1].x -= 7.0
            elif self.mode == 'stair':
                p.waypoints.append(self.waypoint_low)
                p.waypoints.append(vec3(self.waypoint_low))
                p.waypoints[-1].x -= self.doorway_width * 0.75
                p.waypoints.append(self.waypoint_prt_door)
                p.waypoints.append(vec3(self.waypoint_prt_door))
                p.waypoints[-1].x -= 1.0
            p.manager = self.prt_station.person_at_station
            railcar.departing.append(p)

def main(mode='straight'):
    the_world = world()
    initialize()

    rx, ry = 35., 0.

    berths = 9
    psm1 = PRTStationManager(berths)
    activate(psm1, psm1.check(), 0.0)

    door_width = 4.0
    stair_width = door_width

    s1 = Straight((20., -30.), (20., 35.0))
    s1.next, last_left, last_right, length = station((20.0, 35.0), PRT.max_speed, 'right', berths, psm1)
    last_left.next = last_right.next = c1 = Straight((20., 35. + length), (20., 10000.))
    #last_left.next = last_right.next = c1 = Curve(0, pi/2, (-10., 35.0 + length), 30.)
    #c1.next = Straight((-10., 35.0 + length + 30.), (-10000., 35.0 + length + 30.))
    psm1.mainline = last_left
    prt_draw(s1, 20.0, -31.0)

    pc = PRTChecker()
    activate(pc, pc.check(c1.members), 0.0)

    pp_ul, pp_lr = prt_platform = s1.next.platform
    prt_platform_center = (pp_ul.y + pp_lr.y) / 2
    psm1.platform = prt_platform

    if mode in ('straight', 'doorway'):
        rail_stop_center = prt_platform_center
    elif mode == 'stair':
        platform_height = 6
        step_height = 0.19
        step_length = 0.28
        steps = platform_height / step_height
        stair_length = steps * step_length
        rail_stop_center = prt_platform_center - stair_length - door_width

    #rail_draw(rx, ry, prt_platform_center_y)
    rail_draw(rx, ry, rail_stop_center)

    rail_headway = 450
    rg = RailcarGenerator(mode, rail_stop_center, door_width, s1.next, rx, rail_headway)
    activate(rg, rg.generate(), 20.0 + 30)

    # prt platform walls, clockwise from top left
    the_world.add_wall( (pp_lr.x, pp_ul.y), (pp_ul.x, pp_ul.y),       # top right and left
                        (pp_ul.x, pp_lr.y), (pp_lr.x, pp_lr.y) )       # bottom left and right

    railside_x = rx - Railcar.width / 2
    r1_top = rg.rail_destination[0] + Railcar.length / 2
    r2_top = rg.rail_destination[1] + Railcar.length / 2
    r1_doors = []
    r2_doors = []
    for door in Railcar.doors:
        r1_doors.append( ((railside_x, r1_top - door + 1), (railside_x, r1_top - door - 1)) )
        r2_doors.append( ((railside_x, r2_top - door + 1), (railside_x, r2_top - door - 1)) )

    railside_x = rx - Railcar.width / 2

    # rail platform walls
    the_world.add_wall( (rx - 5, ry + 155), (railside_x, ry + 155), # top left and right
                        r1_doors[0][0] )                            # first car, top door
    the_world.add_wall( r1_doors[0][1], r1_doors[1][0] )
    the_world.add_wall( r1_doors[1][1], r1_doors[2][0] )
    the_world.add_wall( r1_doors[2][1], r1_doors[3][0] )
    the_world.add_wall( r1_doors[3][1], r2_doors[0][0] )
    the_world.add_wall( r2_doors[0][1], r2_doors[1][0] )
    the_world.add_wall( r2_doors[1][1], r2_doors[2][0] )
    the_world.add_wall( r2_doors[2][1], r2_doors[3][0] )
    the_world.add_wall( r2_doors[3][1],                             # second car, bottom door
                        (railside_x, ry + 45), (rx - 5, ry + 45))   # bottom right and left

    tk = the_world.tk
    canvas, x_, y_, s_ = tk.canvas, tk.x_, tk.y_, tk.s_
    if mode == 'straight':
        the_world.add_wall( (rx - 5, pp_ul.y + 4),                 # prt/rail door top
                            (rx - 5, ry + 155) )                    # top right
        the_world.add_wall( (rx - 5, ry + 45),                      # bottom left
                            (rx - 5, pp_lr.y - 4) )                # prt/rail door bottom

        rg.waypoint_high = vec3(rx - 4, pp_ul.y + 3)
        rg.waypoint_low = vec3(rx - 4, pp_lr.y - 3)

        # between rail and prt platforms
        the_world.add_wall((pp_lr.x, pp_ul.y), (rx - 5, pp_ul.y + 4))
        the_world.add_wall((rx - 5, pp_lr.y - 4), (pp_lr.x, pp_lr.y))

        canvas.create_polygon( ((x_(pp_lr.x), y_(pp_ul.y)), (x_(rx - 5), y_(pp_ul.y + 4)),
                                (x_(rx - 5), y_(pp_lr.y - 4)), (x_(pp_lr.x), y_(pp_lr.y))),
                               fill='#ccc', outline='#ccc')

    elif mode == 'doorway':
        the_world.add_wall( (pp_lr.x, pp_ul.y), (pp_lr.x, prt_platform_center + door_width / 2),
                            (rx - 5, prt_platform_center + door_width / 2),
                            (rx - 5, ry + 155) )
        the_world.add_wall( (rx - 5, ry + 45), (rx - 5, prt_platform_center - door_width / 2),
                            (pp_lr.x, prt_platform_center - door_width / 2),
                            (pp_lr.x, pp_lr.y) )
        rg.waypoint_high = vec3(rx - 4, prt_platform_center + door_width / 3)
        rg.waypoint_low = vec3(rx - 4, prt_platform_center - door_width / 3)
        canvas.create_polygon( ((x_(pp_lr.x), y_(prt_platform_center + door_width / 2)),
                                (x_(rx - 5), y_(prt_platform_center + door_width / 2)),
                                (x_(rx - 5), y_(prt_platform_center - door_width / 2)),
                                (x_(pp_lr.x), y_(prt_platform_center - door_width / 2))),
                               fill='#ccc', outline='#ccc')
    elif mode == 'stair':
        epd = extra_platform_depth = 3
        wall1 = ( (pp_lr.x, pp_ul.y), (pp_lr.x + epd, prt_platform_center + door_width / 2),
                  (pp_lr.x + epd + stair_width, prt_platform_center + door_width / 2),
                  (pp_lr.x + epd + stair_width, prt_platform_center - door_width / 2),
                  (rx - 5, prt_platform_center - door_width / 2 - stair_length), 
                  (rx - 5, ry + 155) )
        wall2 = ( (rx - 5, ry + 45), (rx - 5, prt_platform_center - door_width / 2 - stair_length - door_width),
                  (rx - 5 - door_width, prt_platform_center - door_width / 2 - stair_length - door_width),
                  (rx - 5 - door_width, prt_platform_center - door_width / 2 - stair_length),
                  (pp_lr.x + epd, prt_platform_center - door_width / 2),
                  (pp_lr.x, pp_lr.y) )
        the_world.add_wall( *wall1 )
        the_world.add_wall( *wall2 )
        rg.waypoint_low = vec3(rx - 4, prt_platform_center - stair_length - door_width)
        rg.waypoint_prt_door = vec3(pp_lr.x + door_width / 2, prt_platform_center)
        points = []
        for point in wall1 + wall2:
            points.append(x_(point[0]))
            points.append(y_(point[1]))
        canvas.create_polygon(points, fill='gray')
        
    tg = PRTGenerator()
    # 50% empty to simulate sending extras to the transit station
    activate(tg, tg.generate(s1, interval, percent_empty=0.50), 0.0)

    clock = Clock()
    # "0" time taken from when Railcar called station_manager.doors_open
    activate(clock, clock.execute(), 43.684 + 30)

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
