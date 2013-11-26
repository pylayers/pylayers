from pylayers.mobility.body.body import *
from pylayers.mobility.trajectory import Trajectory
from IPython.display import Image

plt.ion()
bc = Body()
np.shape(bc.d)
v = bc.vmocap
time = np.arange(0,10,0.01)
x = v*time
y = np.zeros(len(time))
#
traj = Trajectory(time,np.vstack((x,y)).T)
#
bc.settopos(traj=traj,tk=7)
bc.setccs(topos=True)
bc.setaccs()
bc.show3(topos=True,ccs=True,wire=True,pattern=True)
#bc.show3(topos=True,ccs=True,wire=True,pattern=True,velocity=True)
#
