from SimPy.Simulation import *
from Transit.StateProcess import Updater
from Transit.Person import Person
from Transit.World import world
from Transit.SteeringBehavior import Seek, Containment, Wander
from Transit.vec3 import vec3

this is broken!

This simulation hasn't been updated to support waypoints.

walls = [
    [ (( 2.25, 20.25), (13.50,  9.00)), ((15.75,  6.75), (18.00,  9.00), ( 4.50, 22.50)), ],
    [ ((15.75,  6.75), ( 9.00,  0.00), ( 6.75,  2.25)), (( 9.00,  4.50), (13.50,  9.00)), ],
    [ (( 6.75,  2.25), ( 0.00,  9.00), ( 2.25, 11.25)), (( 4.50,  9.00), ( 9.00,  4.50)), ],
    [ (( 2.25, 11.25), (13.50, 22.50)), ((15.75, 24.75), (18.00, 22.50), ( 4.50,  9.00)), ],
    [ ((15.75, 24.75), ( 9.00, 31.50), ( 6.75, 29.25)), (( 9.00, 27.00), (13.50, 22.50)), ],
    [ (( 6.75, 29.25), ( 0.00, 22.50), ( 2.25, 20.25)), (( 4.50, 22.50), ( 9.00, 27.00)), ],
    ]

def draw(walls):
    the_world = world()
    tk = the_world.tk
    canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
    for wall in walls:
        the_world.add_wall(wall)
        for lines in wall:
            scaled_lines = []
            for xx, yy in lines:
                scaled_lines.append((x_(xx), y_(yy)))
            canvas.create_line(scaled_lines, fill='#ccc')

def more(person):
    pass

def main():
    interval = 0.1
    initialize()
    draw(walls)

    xx = 4.0
    people = []

    p = Person(interval)
    p.behaviors = [Wander()]
    p.position = vec3(xx, 21)
    p.destination = vec3(xx, 21)
    people.append(p)
    p.people = people
    activate(p, p.move(), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    simulate(until=240)

if __name__ == '__main__':
    main()
