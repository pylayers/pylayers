import numpy as np
import pylayers.mobility.trajectory as tr
import pylayers.mobility.ban.body as body
import matplotlib.pyplot as plt
time = np.arange(0,10,0.1)
v = 4000/3600.
x = v*time
y = np.zeros(len(time))
traj = tr.Trajectory()
traj.generate()
bc = body.Body()
bc.settopos(traj,2.3,2)
bc.setccs(topos=True)
bc.setdcs()
bc.show(plane='yz',color='b',widthfactor=80)
plt.show()
