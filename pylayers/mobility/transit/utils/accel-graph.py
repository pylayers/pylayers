import biggles
import Numeric, math
from Tkinter import *

prt_berth_advance_trajectory = [
    [0.0, 0.0, 0.0, 0.0, 1.25],
    [1.0, 0.208333, 0.625000, 1.250000, -1.25],
    [2.0, 1.250000, 1.250000, 0.000000, 0.0],
    [2.136000, 1.420000, 1.250000, 0.000000, -1.25],
    [3.136000, 2.461667, 0.625000, -1.250000, 1.25],
    [4.136000, 2.670000, 0.000000, 0.000000, 0.0],
    ]

unit_trajectory = [
    [0.0, 0.0, 0.0, 0.0, 1.0],
    [1.0, 0.1666667, 0.5, 1.0, -1.0],
    [2.0, 1.0, 1.0, 0.0, -1.0],
    [3.0, 1.8333333, 0.5, -1.0, 1.0],
    [4.0, 2.0, 0.0, 0.0, 0.0],
    ]

prt_to_max_jerk = [
    [0.0, 0.0, 0.0, 0.0, 1.25],
    [1.0, 0.208333333333, 0.625, 1.25, -1.25],
    [2.0, 1.25, 1.25, 0.0, -1.25],
    [3.0, 2.29166666667, 0.625, -1.25, 1.25],
    [4.0, 2.5, 0.0, 0.0, 0.0],
    ]


def accel_plot(trajectory, t_step=0.01):
    t_total = trajectory[-1][0] + t_step
    t = Numeric.arange(0, t_total, t_step)
    d = Numeric.arange(0, t_total, t_step)
    v = Numeric.arange(0, t_total, t_step)
    a = Numeric.arange(0, t_total, t_step)
    j = Numeric.arange(0, t_total, t_step)
    traj = list(trajectory)
    print "  time   distance   veloc    accel    jerk"
    for ii in range(0, len(t)):
        if t[ii] > traj[1][0]:
            del traj[0]
            # tt is relative to the trajectory segments t_0
            print "popping trajectory"
        t_0, d_0, v_0, a_0, j_0 = traj[0]
        tt = t[ii] - t_0
        d[ii] = d_0 + v_0*tt + a_0*tt**2/2 + j_0*tt**3/6
        v[ii] = v_0 + a_0*tt + j_0*tt**2/2
        a[ii] = a_0 + j_0*tt
        j[ii] = j_0
        #print "%8.4f %8.4f %8.4f %8.4f %8.4f" % (t[ii], d[ii], v[ii], a[ii], j[ii])

    return t, d, v, a, j

def main(argv):
    (t, d, v, a, j) = accel_plot(prt_berth_advance_trajectory)
    p = biggles.FramedPlot()
    p.title = 'Acceleration Profile'
    p.xlabel = 'Time'
    p.add( biggles.Curve(t, d, color='red') )
    p.add( biggles.Curve(t, v, color='yellow') )
    p.add( biggles.Curve(t, a, color='blue') )
    p.add( biggles.Curve(t, j, color='green') )
    p.show()

if __name__ == '__main__':
    main(sys.argv)
