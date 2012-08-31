from SimPy.Simulation import Process, activate, hold
from Transit.PRT import PRT
from Transit.World import world
from random import random, uniform
from math import *

class PRTGenerator(Process):
    def __init__(self):
        Process.__init__(self)
        self.prts = []

    def generate(self, segment, interval, percent_empty=0.30):
        while True:
            p = PRT(interval, segment, 0.)
            segment.members.append(p)
            self.prts.append(p)
            p.occupied = random() > percent_empty
            activate(p, p.execute(), 0.0)
            yield hold, self, uniform(0.5, 1.5)
            for prt in self.prts:
                if prt.distance > 200:
                    self.cancel(prt)
                    prt.delete()

class PRTChecker(Process):
    def __init__(self):
        Process.__init__(self)
        self.prts = []

    def check(self, prts):
        while True:
            last_prt = None
            for prt in prts:
                #if last_prt is not None and abs(prt.distance - last_prt.distance) < (prt.max_speed * 0.45):
                #    print prt.distance, last_prt.distance, prt.segment, last_prt.segment, prt.speed, last_prt.speed
                if (prt.distance < 150
                    and last_prt is not None
                    and prt.segment is last_prt.segment
                    and abs(prt.distance - last_prt.distance) < (prt.max_speed * 0.45)):
                    if not hasattr(prt, 'color'):
                        print abs(prt.distance - last_prt.distance), prt.speed
                    prt.color = 'purple'
                elif prt.speed < prt.desired_speed:
                    print 'too slow!', prt.speed, prt.desired_speed
                last_prt = prt
            yield hold, self, 0.01
                    
class Clock(Process):
    def __init__(self, stopwatch_seconds=60):
        Process.__init__(self)
        self.stopwatch_seconds = stopwatch_seconds
        tk = world().tk
        canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
        cx, cy = 27, 136
        self.cx, self.cy = cx, cy
        canvas.create_oval(x_(cx - 5), y_(cy - 5), x_(cx + 5), y_(cy + 5), fill='white')
        steps = 60
        if stopwatch_seconds == 60:
            steps = 5
        for ii in range(0, stopwatch_seconds, steps):
            inner_radius = 4.2
            if steps == 5 and ii % 15 in (5, 10):
                inner_radius = 4.6
            radians = self.seconds_to_radians(ii)
            canvas.create_line(x_(cx + cos(radians) * inner_radius), y_(cy + sin(radians) * inner_radius),
                               x_(cx + cos(radians) * 5), y_(cy + sin(radians) * 5))
        self.hand = canvas.create_line(0, 0, 1, 1)

    def seconds_to_radians(self, seconds):
        # rotate counterclockwise so 0 points up
        return -seconds * 2 * pi / self.stopwatch_seconds + pi / 2

    def execute(self):
        stopwatch = 0.0
        tk = world().tk
        canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
        hand, cx, cy = self.hand, self.cx, self.cy
        while True:
            radians = self.seconds_to_radians(stopwatch)
            canvas.coords(hand, x_(cx), y_(cy),
                          x_(cx + cos(radians) * 4), y_(cy + sin(radians) * 4))
            stopwatch += 1.0
            yield hold, self, 1.0
