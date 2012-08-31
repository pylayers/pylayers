"""
Most forward berth is 0.

"""

# Due to the"looking backwards" calculations in PRTLaunch, the minimum
# station size must be 3 berths so that the"look back" distance is
# less than the length of the associated mainline.


from SimPy.Simulation import *
from Transit.PRT import PRT, spline_length, PARKING
from Transit.StateProcess import PASSIVATE, DONT_PASSIVATE
from Transit.SteeringBehavior import Arrive, FollowWaypoints, Separation, Queuing, InterpenetrationConstraint
from Transit.vec3 import vec3

class PRTLaunch(Process):
    """Sub-process to launch a PRT."""
    def execute(self, station):
        prts = station.prts
        headway = PRT.max_speed * PRT.headway
        while True:
            if not (len(prts) and prts[0].ready_to_depart):
                yield hold, self, 0.5
                continue

            prt = prts.pop(0)
            prt.prt_ahead = None

            front_berth_distance = 17.315 + PRT.berth_length * (station.total_berths)
            if prt.distance < front_berth_distance:
                if len(prts):
                    prts[0].prt_ahead = PARKING
                    prts[0].point_distance = 17.315 + PRT.berth_length * (station.total_berths + 1)
                reactivate(prt)
                while prt.distance < front_berth_distance:
                    yield hold, self, 0.01
            elif prt.speed < 0.01:
                prt.waiting_to_depart = True
            else:
                if len(prts):
                    prts[0].prt_ahead = PARKING
                    prts[0].point_distance = 17.315 + PRT.berth_length * (station.total_berths + 1)

            blocking = True
            while blocking:
                blocking = False
                speed_difference = PRT.max_speed - prt.speed
                time_to_reach_speed = speed_difference / PRT.max_acceleration
                average_speed = (PRT.max_speed + prt.speed) / 2
                distance_to_reach_speed = average_speed * time_to_reach_speed
                distance_to_merge = prt.segment.length - prt.distance + prt.segment.next.length
                distance_at_speed = distance_to_merge - distance_to_reach_speed
                distance_to_check = station.mainline.length - (distance_at_speed + time_to_reach_speed * PRT.max_speed)
                for check_prt in station.mainline.members:
                    if (distance_to_check - headway) < check_prt.distance < (distance_to_check + headway):
                        blocking = True
                        break
                if blocking:
                    yield hold, self, 0.01
            if hasattr(station, 'launch_statistics'):
                print now(), "launch"
            prt.ready_to_depart = False
            if prt.waiting_to_depart:
                prt.waiting_to_depart = False
                if len(prts):
                    prts[0].prt_ahead = PARKING
                    prts[0].point_distance = 17.315 + PRT.berth_length * (station.total_berths + 1)
            prt.desired_speed = PRT.max_speed
            reactivate(prt)

class PRTStationManager(Process):
    """
    view_queues are people waiting to arrive at a berth, as a type of"view" of which berths may be busy.
    berth_queues are people at a berth waiting to load.
    """
    # PRTStationManager also acts as a PRT line segment
    type = 'decision'

    def __init__(self, berths, queue_berths=None):
        Process.__init__(self)
        self.berths = berths
        if queue_berths is None:
            self.queue_berths = berths
        else:
            self.queue_berths = queue_berths
        self.total_berths = self.berths + self.queue_berths
        self.view_queues = []
        self.berth_queues = []
        for ii in range(0, berths):
            self.view_queues.append([])
            self.berth_queues.append([])
        self.prts = []
        self.prt_launch = PRTLaunch()
        activate(self.prt_launch, self.prt_launch.execute(self), 0.0)

    def check(self):
        while True:
            # push out empty vehicles in front of a PRT ready to depart
            for check_prt in self.prts:
                if not (check_prt.boarding or check_prt.ready_to_depart):
                    continue
                if check_prt.ready_to_depart:
                    for push_prt in self.prts:
                        if push_prt is check_prt:
                            break
                        push_prt.ready_to_depart = True
                break

            yield hold, self, 0.5

    ###
    ### Interactions with Person
    ###
    def person_at_station(self, person):
        def by_weighted_distance(left, right):
            return cmp(left[1], right[1])
        pp = self.platform
        queue_depths = []
        line_weighting = 10.0
        forward_berth_weighting = 4.0
        for berth in range(0, self.berths):
            queue_depths.append((berth, (float(len(self.view_queues[berth]) + len(self.berth_queues[berth])) * line_weighting
                                         + (person.position - vec3(pp[0].x, pp[0].y - PRT.berth_length / 2 - berth * PRT.berth_length)).length()
                                         - float(self.berths - berth) / forward_berth_weighting)))
        queue_depths.sort(by_weighted_distance)
        berth = queue_depths[0][0]
        self.view_queues[berth].append(person)
        person.manager = self.person_at_berth
        person.manager_args = [berth]
        person.waypoints = [vec3(pp[1].x - 0.5,
                                 pp[0].y - PRT.berth_length / 2 - berth * PRT.berth_length),
                            vec3(pp[0].x + 0.5,
                                 pp[0].y - PRT.berth_length / 2 - berth * PRT.berth_length)]
        person.behaviors = [Separation(), FollowWaypoints(), Queuing(), InterpenetrationConstraint()]
        return DONT_PASSIVATE

    def person_at_berth(self, person, berth):
        self.view_queues[berth].remove(person)
        self.berth_queues[berth].append(person)
        person.behaviors = [Arrive()]
        person.manager = None
        pp = self.platform
        person.destination = vec3(pp[0].x + 0.5 + (len(self.berth_queues[berth]) - 1) * person.radius  * 2 * 1.2,
                                  pp[0].y - PRT.berth_length / 2 - berth * PRT.berth_length)
        person.world.tk.canvas.itemconfig(person.graphic, outline='red')
        return DONT_PASSIVATE

    def person_in_vehicle(self, person, prt):
        person.delete()
        prt.entered += 1
        if prt.entered == prt.occupied:
            prt.door_closing = True
            prt.manager = self.prt_door_closed
            prt.manager_args = []
        return PASSIVATE

    ###
    ### Interactions with PRT
    ###
    def prt_switch(self, prt):
        next_right = self.next.right
        # 6 meters into the straight segment was determind empirically
        if (prt.occupied
            or len(next_right.members)
            or (len(next_right.next.members)
                and (next_right.next.members[-1].distance < 6))):
            return 'left'
        count = 0
        for check_prt in self.prts:
            if check_prt.distance < 17.315 + PRT.berth_length * (self.queue_berths + 0.5):
                count += 1
        if count == self.queue_berths:
            return 'left'
        if len(self.prts):
            prt.prt_ahead = self.prts[-1]
        else:
            prt.prt_ahead = PARKING
            prt.point_distance = 17.315 + PRT.berth_length * (self.total_berths + 1)
        prt.manager = self.prt_stopped_at_berth
        prt.manager_args = []
        prt.desired_speed = prt.station_speed
        self.prts.append(prt)
        return 'right'

    def prt_door_open(self, prt, berth):
        prt.occupied = 0
        def load_prt(prt, depth):
            person = self.berth_queues[berth].pop(0)
            person.destination.x -= 1.2
            person.manager = self.person_in_vehicle
            person.manager_args = [prt]
            reactivate(person)
            prt.occupied += 1
        load_prt(prt, 0)
        if len(self.berth_queues[berth]): load_prt(prt, 1)
        if len(self.berth_queues[berth]): load_prt(prt, 2)
        if len(self.berth_queues[berth]): load_prt(prt, 3)
        prt.entered = 0
        prt.update()
        pp = self.platform
        destination_x = pp[0].x + 0.5
        for person in self.berth_queues[berth]:
            person.destination.x = destination_x
            destination_x += person.radius * 2 * 1.2
            reactivate(person)

    def prt_door_closed(self, prt):
        prt.ready_to_depart = True
        prt.boarding = False

    def prt_stopped_at_berth(self, prt):
        # @@ "standard" ramp @ 11 m/s is 25.088m, - spline_length is 17.315
        berth = self.total_berths - int(round((prt.distance - 17.315) / PRT.berth_length))
        if berth < self.berths and len(self.berth_queues[berth]):
            prt.door_opening = True
            prt.boarding = True
            prt.manager = self.prt_door_open
            prt.manager_args = [berth]
