from SimPy.Simulation import *
from Tkinter import *
from math import *
from pylayers.mobility.transit.StateProcess import StateProcess, Updater, new_state
from pylayers.mobility.transit.World import world


# Easements
#   http://prtnews.com/wiki/EaseMents

class Straight:
    type = 'straight'

    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.length = sqrt((end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2)
        self.members = []

    def position(self, distance):
        along = distance / self.length
        xx = (self.end[0] - self.start[0]) * along + self.start[0]
        yy = (self.end[1] - self.start[1]) * along + self.start[1]
        return xx, yy, 0.0


class Curve:
    type = 'curve'

    def __init__(self, start, end, center, radius):
        self.start = start
        self.end = end
        self.center = center
        self.radius = radius
        self.length = self.radius * abs(end - start)
        self.members = []

    def position(self, distance):
        if self.start < self.end:
            along = self.start + distance / self.length \
                            * abs(self.end - self.start)
        else:
            along = self.start - distance / self.length \
                            * abs(self.end - self.start)
        xx = self.center[0] + cos(along) * self.radius
        yy = self.center[1] + sin(along) * self.radius
        return xx, yy, pi - along


def B0(tt):
    return tt ** 3


def B1(tt):
    return 3 * tt ** 2 * (1 - tt)


def B2(tt):
    return 3 * tt * (1 - tt) ** 2


def B3(tt):
    return (1 - tt) ** 3


class Spline:
    type = 'spline'

    def __init__(self, p0, p1, p2, p3):
        self.members = []
        self.x0, self.y0 = p0
        self.x1, self.y1 = p1
        self.x2, self.y2 = p2
        self.x3, self.y3 = p3
        self.segments = [(0., self.x0, self.y0)]
        last_x, last_y = p0
        length = 0
        for ii in range(1, 21):
            xx, yy = self.spline(1. - float(ii) / 20)
            length += sqrt((xx - last_x) ** 2 + (yy - last_y) ** 2)
            self.segments.append((length, xx, yy))
            last_x, last_y = xx, yy
        self.length = length

    def spline(self, tt):
        xx = self.x0 * B0(tt) + self.x1 * B1(tt) \
            + self.x2 * B2(tt) + self.x3 * B3(tt)
        yy = self.y0 * B0(tt) + self.y1 * B1(tt) \
            + self.y2 * B2(tt) + self.y3 * B3(tt)
        return xx, yy

    def position(self, distance):
        segments = self.segments
        this_distance, this_x, this_y = segments[0]
        for ii in range(1, len(segments)):
            next_distance, next_x, next_y = segments[ii]
            if distance - next_distance < 1e-10:
                break
            this_distance, this_x, this_y = segments[ii]
        segment_length = next_distance - this_distance
        distance_along_segment = distance - this_distance
        ratio = distance_along_segment / segment_length
        dx = next_x - this_x
        dy = next_y - this_y
        return this_x + dx * ratio, this_y + dy * ratio, 0.0  # atan2(dx, dy), looks odd with splines instead of easements


# used only in the __main__ test below
class _Mover(StateProcess):
    def __init__(self, interval):
        StateProcess.__init__(self)
        self.interval = interval

    def move(self):
        speed = 6.
        for segment in self.segments:
            seconds_to_distance = segment.length / speed
            segment_distance = 0
            for ii, seconds in self.pace(segment.length / speed):
                segment_distance += speed * seconds
                self.x, self.y, tangent = segment.position(segment_distance)
                self.update()
                yield hold, self, seconds
            assert abs(segment_distance - segment.length) < 0.00001

    def update(self):
        tk = world().tk
        if not hasattr(self, 'graphic'):
            self.graphic = tk.canvas.create_oval(0, 0, 1, 1,)
        x_, y_ = tk.x_, tk.y_
        tk.canvas.coords(self.graphic,
                         x_(self.x - 0.21), y_(self.y - 0.21),
                         x_(self.x + 0.21), y_(self.y + 0.21))
        if hasattr(self, 'last_x'):
            tk.canvas.create_line(x_(self.last_x), y_(
                self.last_y), x_(self.x), y_(self.y))
        self.last_x, self.last_y = self.x, self.y

if __name__ == '__main__':
    interval = 0.1
    initialize()

    m = _Mover(interval)
    m.segments = (Straight((30., 10.), (30., 40.)),
                  Curve(pi, pi / 2, (60., 40.), 30),
                  Curve(-pi / 2, 0, (60., 100.), 30),
                  Straight((90., 100.), (90., 150.)),
                  Spline((90., 150.), (90., 200.), (30., 200.),
                         (30., 150.))
                  )
    m.x, m.y = (75.0, 600.0)
    activate(m, m.move(), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    simulate(until=50)
    import sys
    sys.stdin.readline()
