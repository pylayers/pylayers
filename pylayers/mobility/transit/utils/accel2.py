"""
This script calculates a stop-to-stop trajectory by binary mid-point
search.

still broken.

"""

distance = 2.67 / 2
comfortable_jerk = 1.25 # m/s^3
max_acceleration = 2.5 # m/s^2, rounded just above .25g
max_velocity = 1.34 # 3 miles/hour in station
time_to_max_accel = max_acceleration / comfortable_jerk
print "max_time = %8.5f" % time_to_max_accel

upper_time = time_to_max_accel
lower_time = 0.0

d_0 = v_0 = a_0 = 0.0
j_0 = comfortable_jerk

def distance_at_time(t):
    t_m = t / 2
    t_f = t - t_m # same as midway, since it's half
    d_m = d_0 + v_0*t_m + a_0*t_m**2/2 + j_0*t_m**3/6
    v_m = v_0 + a_0*t_m + j_0*t_m**2/2
    a_m = a_0 + j_0*t_m
    print t_m, d_m, v_m, a_m, -j_0
    d_f = d_m + v_m*t_f + a_m*t_f**2/2 + -j_0*t_f**3/6
    v_f = v_m + a_m*t_f + -j_0*t_f**2/2
    a_f = a_m + -j_0*t_f
    print t, d_f, v_f, a_f, -j_0
    return d_f

print "dinstance at upper_time = %8.5f" % distance_at_time(time_to_max_accel)

while abs(upper_time - lower_time) > 0.00001:
    mid_time = lower_time + (upper_time - lower_time) / 2
    check_distance = distance_at_time(mid_time)
    if check_distance > distance:
        lower_time = mid_time
    else:
        upper_time = mid_time

print upper_time, lower_time
