import numpy as np
import pylayers.mobility.trajectory as tr
import pylayers.mobility.body.body as body
import matplotlib.pyplot as plt
time = np.arange(0,10,0.1)
v = 4000/3600.
x = v*time
y = np.zeros(len(time))
traj = tr.Trajectory()
John = body.Body()
John.settopos(traj,2.3)
John.show()
