from SimPy.Simulation import *
from Transit.StateProcess import StateProcess, Updater, new_state
from World import world
import os

class Railcar(StateProcess):
    distance_to_stop = 225.61739 # m
    acceleration = 1.34112 # m/s/s
    max_speed = 24.6 # m/s
    dwell = 60 # s
    width = 2.75 # m
    length = 28.65 # m
    doors = [5.0, 11.0, 17.65, 23.65]
    sections = [1.0, 7.5, 14.33, 21.15, 27.65]
    platform_length = 67 # 220 feet


    def __init__(self, name, interval):
        Process.__init__(self)
        self.world = world()
        self.interval = interval
        self.name = name
        self.x = 35
        self.position = 0
        self.speed = 24.6 # m/s
        self.destination = 240.0
        self.departing = []
        self.station_manager = None
        self.state_generator = self.cruising()

    def accelerate(self, seconds, acceleration):
        new_speed = self.speed + acceleration * seconds
        average_speed = (self.speed + new_speed) / 2
        self.position += average_speed * seconds
        self.speed = new_speed
        self.update()

    # states

    def cruising(self):
        distance_to_braking = abs(self.destination - self.position) - self.distance_to_stop
        for ii, seconds, distance in self.pace(distance_to_braking / self.speed, distance_to_braking):
            yield hold, self, seconds
            self.position += distance
            self.update()
        yield new_state, self.braking()

    def braking(self):
        for ii, seconds in self.pace((self.speed - 0) / self.acceleration):
            yield hold, self, seconds
            self.accelerate(seconds, -self.acceleration)
        assert abs(self.speed - 0) < 0.00001
        yield new_state, self.in_station()

    def in_station(self):
        # @@ station manager should be based on guideway, not "fixed" on 'self'
        if self.station_manager: self.station_manager.arrived(self)
        yield hold, self, 5
        if self.station_manager: self.station_manager.doors_open(self)
        if not self.station_manager or self.station_manager.mode is "normal":
            yield hold, self, self.dwell - 10
        elif self.station_manager and self.station_manager.mode is "stadium":
            while self.station_manager.boarded < 400:
                yield hold, self, 1
        if self.station_manager: self.station_manager.closing_doors(self)
        yield hold, self, 5
        if self.station_manager: self.station_manager.leaving(self)
        yield new_state, self.accelerating()

    def accelerating(self):
        for ii, seconds in self.pace((self.max_speed - 0) / self.acceleration):
            yield hold, self, seconds
            self.accelerate(seconds, self.acceleration)
        assert abs(self.speed - self.max_speed) < 0.00001
        yield new_state, self.cruising()

    # views

    def update(self):
        tk = self.world.tk
        if tk is None: return
        if not hasattr(self, 'graphic'):
            # these colors are light steel blue and steel blue in Unix X11
            self.graphic = tk.canvas.create_rectangle(0, 0,  1, 1, fill='#B0C4DE', outline='#4682B4')
            tk.canvas.addtag('vehicle', 'withtag', self.graphic)
        x_, y_ = tk.x_, tk.y_
        tk.canvas.coords(self.graphic,
                         x_(self.x - self.width/2), y_(self.position - self.length/2),
                         x_(self.x + self.width/2), y_(self.position + self.length/2))

def rail_draw(last_x, last_y, platform_center_y=100):
    tk = world().tk
    canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
    # destinations (car 1 and 2)
    if os.environ.has_key('DESTINATIONS'):
        tk.canvas.create_line(x_(last_x-2), y_(last_y+86), x_(last_x+2), y_(last_y+86), fill='#000')
        tk.canvas.create_line(x_(last_x-2), y_(last_y+114), x_(last_x+2), y_(last_y+114), fill='#000')
    # rails
    tk.canvas.create_line(x_(last_x-0.717), y_(last_y),
                                x_(last_x-0.717), y_(last_y+200), fill='#888')
    tk.canvas.create_line(x_(last_x+0.717), y_(last_y),
                                x_(last_x+0.717), y_(last_y+200), fill='#888')
    # platform
    warning_strip_width = 0.610 # 24 inches
    warning_strip_color = '#ccc' # '#ff6' for a darker yellow, '#ccc' to hide
    plat_top = platform_center_y + Railcar.platform_length / 2
    plat_bot = platform_center_y - Railcar.platform_length / 2
    tk.canvas.create_rectangle(x_(last_x - 5), y_(plat_top),
                                     x_(last_x - Railcar.width / 2 - warning_strip_width), y_(plat_bot),
                               fill='#ccc', outline="")
    tk.canvas.create_rectangle(x_(last_x - Railcar.width / 2 - warning_strip_width), y_(plat_bot),
                                     x_(last_x - Railcar.width / 2), y_(plat_top), fill=warning_strip_color, outline="")


if __name__ == '__main__':
    interval = 0.5
    initialize()

    r = Railcar('railcar', interval)
    r.position = -120
    r.destination = 114
    activate(r, r.execute(), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    simulate(until=240)
