from pylayers.mobility.trajectory import *
from pylayers.mobility.ban.body import *
b = Body()
traj = Trajectory()
traj.generate()
b.settopos(traj,t=3,cs=True)
# #>>> b._show3(topos=True,pattern=True)
