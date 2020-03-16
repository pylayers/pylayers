def accel(t, d_0, v_0, a_0, j_0):
    d_t = d_0 + v_0 * t + a_0 * t ** 2.0 / 2.0 + j_0 * t ** 3.0 / 6.0
    v_t = v_0 + a_0 * t + j_0 * t ** 2.0 / 2.0
    a_t = a_0 + j_0 * t
    return d_t, v_t, a_t

def prt_max_jerk():
    print("prt_max_jerk")
    print("[0.0, 0.0, 0.0, 0.0, 1.25],")
    d_max_a, v_max_a, a_max_a = accel(1.0, 0.0, 0.0, 0.0, 1.25)
    print("[1.0, " + str(d_max_a) + ", "  + str(v_vmax_a) + ", " + str(a_max_a) + ", -1.25]")
    d_mid, v_mid, a_mid = accel(1.0, d_max_a, v_max_a, a_max_a, -1.25)
    print("[2.0, " + str(d_mid) + ", " + str(v_mid) + ", " + str(a_mid) + "-1.25],")
    d_max_a2, v_max_a2, a_max_a2 = accel(1.0, d_mid, v_mid, a_mid, -1.25)
    print("[3.0, " + str(d_max_a2) + ", " + str(v_max_a2) + ", " + str(a_max_a2) + ", 1.25],")
    d_fin, v_fin, a_fin = accel(1.0, d_max_a2, v_max_a2, a_max_a2, 1.25)
    print("[4.0, " + str(d_fin) + ", " + str(v_fin) + ", " + str(a_fin) + ", 0.0],")

def prt_to_one_berth_higher_velocity():
    # ie. does not reach max_station_velocity
    print("prt_to_one_berth_higher_velocity")
    print("[0.0, 0.0, 0.0, 0.0, 1.25],")
    d_max_a, v_max_a, a_max_a = accel(1.0, 0.0, 0.0, 0.0, 1.25)
    print("[1.0, " + str(d_max_a) + ", " + str(v_max_a) + ", " + str(a_max_a) + ", 0.0],")
    t_a = 0.0446683
    d_aft, v_aft, a_aft = accel(t_a, d_max_a, v_max_a, a_max_a, 0.0)
    print("[" + str(t_a+1.0= + ", " + str(d_aft) + ", " + str(v_aft) + ", " + str(a_aft) + ", -1.25],")
    d_mid, v_mid, a_mid = accel(1.0, d_aft, v_aft, a_aft, -1.25)
    print("[" + str(t_a+2.0) + ", " + str(d_mid) + ", " + str(v_mid) + ", " + str(a_mid) + ", -1.25],")
    d_aft2, v_aft2, a_aft2 = accel(1.0, d_mid, v_mid, a_mid, -1.25)
    print("[" + str(t_a+3.0) + ", " + str(d_aft2) + ", " + str(v_aft2) + ", " + str(a_aft2) + ", 0.0],")
    d_max_a2, v_max_a2, a_max_a2 = accel(t_a, d_aft2, v_aft2, a_aft2, 0.0)
    print("[" + str(t_a*2 + 3.0) + ", " + str(d_max_a2) + ", " + str(v_max_a2) + ", " + str(a_max_a2) + ", 1.25],")
    d_fin, v_fin, a_fin = accel(1.0, d_max_a2, v_max_a2, a_max_a2, 1.25)
    print("[" + str(t_a*2 + 4.0) + ", " + str(d_fin) + ", " + str(v_fin) + ", " + str(a_fin) + ", 0.0],")

def prt_to_one_berth_max_velocity():
    # in this one, max_v = max_a = max_j = 1.25
    # max_v of 1.25 is just short of 3mph (2.8mph)
    print("prt_to_one_berth_max_velocity")
    print("[0.0, 0.0, 0.0, 0.0, 1.25],")
    d_max_a, v_max_a, a_max_a = accel(1.0, 0.0, 0.0, 0.0, 1.25)
    print("[1.0, " + str(d_max_a) + ", " + str(v_max_a) + ", " + str(a_max_a) + ", -1.25],")
    d_max_v, v_max_v, a_max_v = accel(1.0, d_max_a, v_max_a, a_max_a, -1.25)
    print("[2.0, " + str(d_max_v) + ", " + str(v_max_v) + ", " + str(a_max_v) + ", 0.0],")
    t_v = 0.136
    d_aft, v_aft, a_aft = accel(t_v, d_max_v, v_max_v, a_max_v, 0.0)
    print("[" + str(t_v+2.0) + ", " + str(d_aft) + ", " + str(v_aft) + ", " + str(a_aft) + ", -1.25],")
    d_max_a2, v_max_a2, a_max_a2 = accel(1.0, d_aft, v_aft, a_aft, -1.25)
    print("[" + str(t_v+3.0) + ", " + str(d_max_a2) + ", " + str(v_max_a2) + ", " + str(a_max_a2) + ", 1.25],")
    d_fin, v_fin, a_fin = accel(1.0, d_max_a2, v_max_a2, a_max_a2, 1.25)
    print("[" + str(t_v+4.0) + ", " + str(d_fin) + ", " + str(v_fin) + ", " + str(a_fin) + ", 0.0],")

prt_max_jerk()
prt_to_one_berth_higher_velocity()
prt_to_one_berth_max_velocity()
