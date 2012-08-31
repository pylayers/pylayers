from SimPy.Simulation import *
from Transit.StateProcess import PASSIVATE, DONT_PASSIVATE
from Transit.World import world

# this should be handled by notifications, not polling
class RailStationManager(Process):
    def __init__(self):
        Process.__init__(self)
        self.world = world()
        self.state = "waiting for train"
        self.boarders = []
        self.boarded = 0
        self.mode = "normal" # or "stadium"
        self.destination_x = 81

    def check(self):
        tk = self.world.tk
        while True:
            if self.state == "move'em in":
                for person in self.boarders:
                    person.destination.x = self.destination_x
                    person.manager = self.person_in_vehicle
                    activate(person, person.move(), 0.0)
            yield hold, self, 0.5

    ###
    ### Interactions with Person
    ###
    def person_at_platform(self, person):
        person.manager = None
        person.world.tk.canvas.itemconfig(person.graphic, outline='red')
        self.boarders.append(person)
        return PASSIVATE

    def person_in_vehicle(self, person):
        self.boarded += 1
        person.delete()
        return PASSIVATE

    ###
    ### Interactions with rail car
    ###
    def arrived(self, railcar):
        self.state = "waiting for doors to open"

    def doors_open(self, railcar):
        self.state = "move'em in"
        for person in railcar.departing:
            activate(person, person.move(), 0.0)
        reactivate(self)

    def closing_doors(self, railcar):
        self.state = "hold'em off"

    def leaving(self, railcar):
        print "boarded: ", self.boarded
        self.boarded = 0
        reactivate(self)
        self.state = "waiting for train"
