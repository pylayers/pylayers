"""

Simulates a "packed" station to get throughput statistics.

The number of berths is parameter, and defaults to 9.

    thruput 5

With the '--analyze' option and the name of a file that contains the
output of the simulation, statistics are generated.

    thruput --analyze thruput.log

Each run in the log file must begin with a line of the form:

    ---- 5 ----

where the number is the number of berths for the following run.

Multiple runs can be generated using these commads to the shell:

    for ii in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
        echo ---- $ii ----;
        SCALE=5 PYTHONPATH=`pwd` python sims/thruput.py $ii;
    done 2>&1 | tee /tmp/thruput.log

"""

from Transit.Path import *
from Transit.Person import Person
from Transit.SteeringBehavior import Arrive
from Transit.StateProcess import Updater
from Transit.PRT import PRT, station, prt_draw
from Transit.vec3 import vec3
from SimPy.Simulation import *
from CommonUtils import PRTGenerator, PRTChecker
from random import uniform

class PersonGenerator(Process):
    def __init__(self, station):
        Process.__init__(self)
        self.station = station

    def execute(self):
        station = self.station
        platform = station.platform
        while True:
            for berth in range(0, len(station.berth_queues)):
                queue = station.berth_queues[berth]
                if len(queue) is 0:
                    p = Person(interval=interval)
                    p.behaviors = [Arrive()]
                    p.position = vec3(platform[0].x + 2,
                                      platform[0].y - PRT.berth_length / 2 - berth * PRT.berth_length)
                    p.destination = vec3(platform[0].x + 0.5,
                                      platform[0].y - PRT.berth_length / 2 - berth * PRT.berth_length)
                    p.velocity = vec3(-0.1, 0.0)
                    activate(p, p.move(), 0.0)
                    queue.append(p)
            yield hold, self, 0.01

sample_time = 900

def analyze(log_file):
    log = open(log_file, 'r')
    first_time = None
    count = 0
    def summary():
        per_hour = count * (3600 / sample_time)
        per_minute = float(per_hour) / 60.
        per_berth = per_minute / berths
        print per_hour, "%.2f" % per_minute, "%.2f" % per_berth
    for line in log:
        if line[0:4] == '----':
            if count != 0:
                summary()
            print line,
            count = 0
            first_time = None
            dash, berths, dash = line.split()
            berths = int(berths)
            continue
        values = line.split()
        departure_time = float(values[0])
        if first_time is None:
            first_time = departure_time
        if departure_time < first_time + sample_time:
            count += 1
    summary()
        

if __name__ == '__main__':
    interval = 1./30
    initialize()
    from Transit.PRTStation import PRTStationManager
    import sys

    if len(sys.argv) == 3 and sys.argv[1] == '--analyze':
        analyze(sys.argv[2])
        sys.exit(0)

    side = 'right'
    if len(sys.argv) == 2:
        berths = int(sys.argv[1])
    else:
        berths = 9
    queue_berths = berths
    psm1 = PRTStationManager(berths, queue_berths)
    psm1.launch_statistics = True
    activate(psm1, psm1.check(), 0.0)
    s1 = Straight((30., 10.), (30., 50.))
    s1.next, last_left, last_right, length = station((30.0, 50.0), PRT.max_speed, side, berths, psm1, queue_berths)
    last_left.next = last_right.next = Straight((30., 50. + length), (30., 10000.))
    psm1.mainline = last_left
    prt_draw(s1, 30.0, 0.0)

    tg = PRTGenerator()
    activate(tg, tg.generate(s1, interval, percent_empty=1.0), 0.0)

    pg = PersonGenerator(psm1)
    activate(pg, pg.execute(), 0.0)

    pc = PRTChecker()
    activate(pc, pc.check(last_left.next.members), 0.0)

    u = Updater(interval)
    activate(u, u.execute(), 0.0)
    # extra minute to allow for start when first PRT leaves berth
    simulate(until=sample_time + 60)
