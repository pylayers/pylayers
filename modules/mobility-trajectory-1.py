from pylayers.mobility.trajectory import *
import matplotlib.pyplot as plt 
import numpy as np 
t = np.arange(0,10,0.01)
x = 2*t*np.cos(t)
y = 3*t*np.sin(t) 
z = 0*t
pt =np.vstack((x,y,z)).T
traj = Trajectory(t,pt)
f,a = traj.plot()
plt.show()
