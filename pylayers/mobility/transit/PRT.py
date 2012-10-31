from SimPy.Simulation import *
from pylayers.mobility.transit.StateProcess import StateProcess, Updater, new_state, next_state
from pylayers.mobility.transit.Path import *
from pylayers.mobility.transit.World import world
from pylayers.mobility.transit.vec3 import vec3
from math import *
from random import lognormvariate, uniform
import os

PARKING = 'parking'

class PRT(StateProcess):
    # station thruput, 300 for two, 1300 for twelve; probably not be linear as assumed
    # berths:  2    3    4    5    6    7  ...  12
    #          300  437  566  687  799  903     1300
    max_jerk = 2.5 # m/s^3
    max_acceleration = 2.5 # m/s^2, rounded just above 0.25G
    max_speed = 11.2    # m/s, 25mph
    #station_speed = 1.25 # m/s, just under 3mph
    station_speed = 2.5 # m/s, just under 3mph
    comfortable_jerk = 1.25 # m/s^3
    comfortable_acceleration = 1.25 # m/s^2
    headway = 0.5 # s
    dwell = 20
    width = 1.27        # 50 inches
    length = 2.54       # 100 inches; alternatively 8.5 feet, or 102 inches
    berth_length = 2.67 # 105 inches; alternatively 9 feet, or 108 inches
    door_time = 2
    
    def __init__(self, interval, segment, distance):
        Process.__init__(self)
        self.world = world()
        self.interval = interval
        self.segment = segment
        self.distance = distance # along segment
        self.speed = PRT.max_speed
        self.desired_speed = PRT.max_speed
        self.acceleration = 0.0
        self.switch = 'right'
        self.occupied = False
        self.ready_to_depart = False
        self.waiting_to_depart = False
        self.boarding = False
        self.tangent = 0.0
        self.door_opening = False
        self.door_closing = False
        self.prt_ahead = None
        self.point_segment = None
        self.point_distance = 0.0
        self.manager = None
        self.manager_args = []

    def execute(self):
        interval = self.interval
        while True:
            if self.door_opening:
                yield hold, self, self.door_time
                self.manager(self, *self.manager_args)
                self.door_opening = False
                continue

            if self.door_closing:
                # This is intended to simulate the person's boarding
                # and situating themselves and pressing "go", beyond
                # the walk time from the platform (which is handled by
                # the person movement and arrival).  However, I don't
                # know where this factor of "2 to 20 seconds" comes
                # from (cause or observed).  lognormvariate is used to
                # give a peak around 3 seconds (plus walk time,
                # recall) with more on the upper side.
                yield hold, self, self.door_time + lognormvariate(1.3, 0.5) + 2.0
                #yield hold, self, self.door_time + uniform(2, 20)
                self.manager(self, *self.manager_args)
                self.door_closing = False
                continue

            if self.boarding or self.waiting_to_depart:
                yield hold, self, 0.01
                continue

            interval = self.interval
            prt_ahead = self.prt_ahead
            if prt_ahead is PARKING:
                interval = 0.0020
                separation = self.point_distance - self.distance
                distance_to_point = separation - self.berth_length
                self.acceleration = 0.0
                if self.speed > self.station_speed:
                    self.acceleration = -self.max_acceleration
                else:
                    time_to_stop = self.speed / self.comfortable_acceleration
                    distance_to_stop = self.speed / 2 * time_to_stop
                    if distance_to_point <= distance_to_stop:
                        self.acceleration = -self.comfortable_acceleration
                    elif self.speed < self.station_speed:
                        self.acceleration = self.comfortable_acceleration
                self.speed += self.acceleration * interval
                if self.speed < 0:
                    self.speed = 0.0
            elif prt_ahead is None:
                speed_difference = self.speed - self.desired_speed
                interval_acceleration = min(abs(speed_difference), self.max_acceleration * interval)
                if speed_difference < 0:
                    self.speed += interval_acceleration
                else:
                    self.speed -= interval_acceleration
#                 if self.speed > self.desired_speed:
#                     self.acceleration = -self.max_acceleration
#                 elif self.speed < self.desired_speed:
#                     self.acceleration = self.max_acceleration
#                 self.speed += self.acceleration * interval
            else:
                if self.segment is prt_ahead.segment:
                    separation = prt_ahead.distance - self.distance
                else:
                    separation = prt_ahead.distance + self.segment.length - self.distance
                ahead_speed = prt_ahead.speed
                ahead_acceleration = prt_ahead.acceleration
                following_distance = max(ahead_speed * 0.5, self.berth_length)
                distance_to_point = separation - following_distance
                self.acceleration = 0.0
                if self.speed > self.station_speed:
                    self.acceleration = -self.max_acceleration
                else:
                    time_to_stop = self.speed / self.comfortable_acceleration
                    distance_to_stop = self.speed / 2 * time_to_stop
                    if distance_to_point <= distance_to_stop:
                        self.acceleration = -self.comfortable_acceleration
                self.speed += self.acceleration * interval
                if self.speed < ahead_speed:
                    self.speed = ahead_speed
                elif self.speed < 0.0:
                    self.speed = 0.0
            self.move(self.speed * interval)
            self.x, self.y, self.tangent = self.segment.position(self.distance)
            self.update()
            if abs(self.speed) == 0.0:
                self.speed = 0.0
                if self.manager:
                    self.manager(self, *self.manager_args)
            elif prt_ahead is not PARKING:
                if self.segment is self.point_segment and self.distance > self.point_distance:
                    if self.manager:
                        self.manager(self, *self.manager_args)

            yield hold, self, interval

    def move(self, delta_x):
        self.distance += delta_x
        if self.distance > self.segment.length:
            self.segment.members.remove(self)
            self.distance -= self.segment.length
            self.segment = self.segment.next
            while self.segment.type in ('diverge', 'decision'):
                if self.segment.type == 'diverge':
                    self.segment = getattr(self.segment, self.switch)
                elif self.segment.type == 'decision':
                    self.switch = self.segment.prt_switch(self)
                    self.segment = self.segment.next
            self.segment.members.append(self)

    def delete(self):
        self.world.tk.canvas.delete(self.graphic)

    # GUI
    def update(self):
        tk = self.world.tk
        if tk is None: return
        if not hasattr(self, 'graphic'):
            self.graphic = tk.canvas.create_polygon(0, 0,  1, 0, 1, 1, 0, 1, fill='black', outline='#333')
            tk.canvas.addtag('vehicle', 'withtag', self.graphic)
        x_, y_ = tk.x_, tk.y_
        tangent = self.tangent
        x, y = self.x, self.y
        x1, y1 = cos(tangent)*-self.width/2 - sin(tangent)*self.length/2, sin(tangent)*-self.width/2 + cos(tangent)*self.length/2
        x2, y2 = cos(tangent)*self.width/2 - sin(tangent)*self.length/2, sin(tangent)*self.width/2 + cos(tangent)*self.length/2
        tk.canvas.coords(self.graphic,
                         x_(x - x1), y_(y + y1), x_(x - x2), y_(y + y2),
                         x_(x + x1), y_(y - y1), x_(x + x2), y_(y - y2))
        if hasattr(self, 'color'):
            tk.canvas.itemconfigure(self.graphic, fill=self.color)
        elif self.boarding:
            tk.canvas.itemconfigure(self.graphic, fill='orange')
        elif self.occupied:
            tk.canvas.itemconfigure(self.graphic, fill='red')
        else:
            tk.canvas.itemconfigure(self.graphic, fill='green')

class Diverge:
    type = 'diverge'

# 1.8 is the guideway separation, but is wrong for the length (too
# sharp of turn).  Of course, splines are the wrong curve anyway.
# Should be using easements.
spline_height = 7.5
spline_length = Spline((0.0, 0.0), (0.9, 0.9), (1.8, spline_height-0.9), (1.8, spline_height)).length

# This is the simplest possible station on a straight guideway
def station(start, mainline_speed, side, berths, decision, queue_berths=None):
    if side == 'left':
        side = -1
    else:
        side = 1
    if queue_berths is None:
        queue_berths = berths
    seconds = mainline_speed / PRT.max_acceleration
    average_ramp_speed = (mainline_speed + 0) / 2
    ramp_length = average_ramp_speed * seconds
    offline_guideway_straight_length = ramp_length * 2 - spline_length * 2 + (berths + queue_berths) * PRT.berth_length
    online_guideway_straight_length = offline_guideway_straight_length + spline_height * 2
    xx, yy = start
    d1 = decision
    d1.next = d2 = Diverge()
    d2.left = s1 = Straight((xx, yy), (xx, yy+online_guideway_straight_length))

    d2.right = sp1 = Spline((xx, yy), (xx, yy+spline_height / 2),
                            (xx+(1.8*side), yy+spline_height / 2), (xx+(1.8*side), yy+spline_height))
    yy_top = yy+spline_height+offline_guideway_straight_length
    sp1.next = s2 = Straight((xx+(1.8*side), yy+spline_height), (xx+(1.8*side), yy_top))
    s2.destinations = []
    for ii in range(0, berths + queue_berths):
        s2.destinations.append(vec3(xx+(1.8*side), yy_top - (ramp_length - spline_length) - ii * PRT.berth_length))
    s2.platform = (vec3(xx+(1.8+PRT.width/2)*side, yy_top - (ramp_length - spline_length) + PRT.berth_length / 2),
                   vec3(xx+(1.8+PRT.width/2+3)*side, yy_top - (ramp_length - spline_length) - berths * PRT.berth_length
                        + PRT.berth_length / 2))
    d1.destinations = s2.destinations
    d1.platform = s2.platform
    s2.next = sp2 = Spline((xx+(1.8*side), yy_top), (xx+(1.8*side), yy_top+spline_height/2),
                           (xx, yy_top+spline_height/2), (xx, yy_top+spline_height))
    return d1, s1, sp2, online_guideway_straight_length

def prt_draw(segment, last_x, last_y):
    tk = world().tk
    canvas, x_, y_, s_ = tk.canvas, tk.x_, tk.y_, tk.s_
    if segment.type == 'decision':
        prt_draw(segment.next, last_x, last_y)
    elif segment.type == 'diverge':
        prt_draw(segment.left, last_x, last_y)
        prt_draw(segment.right, last_x, last_y)
    else:
        for ii in range(0, 10):
            xx, yy, tangent = segment.position(segment.length * float(ii) / 10)
            canvas.create_line(x_(last_x), y_(last_y), x_(xx), y_(yy), fill='#888', width=s_(0.914))
            last_x, last_y = xx, yy
        if hasattr(segment, 'platform'):
            platform = segment.platform
            canvas.create_rectangle(x_(platform[0].x), y_(platform[0].y),
                                    x_(platform[1].x), y_(platform[1].y), fill='#ccc', outline="")
            if os.environ.has_key('DESTINATIONS'):
                for destination in segment.destinations:
                    canvas.create_line(x_(destination.x-1), y_(destination.y),
                                       x_(destination.x+1), y_(destination.y))
        if hasattr(segment, 'next'):
            prt_draw(segment.next, last_x, last_y)

if __name__ == '__main__':
    interval = 1./30
    initialize()
    from PRTStation import PRTStationManager

    side = 'left'
    berths = 10
    s1 = Straight((30., 10.), (30., 50.))
    s1.next, last_left, last_right, length = station((30.0, 50.0), PRT.max_speed, side, berths, PRTStationManager(berths))
    last_left.next = last_right.next = Straight((30., 50. + length), (30., 10000.))
    prt_draw(s1, 30.0, 0.0)

    for ii in range(0, 20):
        p = PRT(interval, s1, 0.)
	s1.members.append(p)
        activate(p, p.execute(), ii)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    simulate(until=240)
