from pylayers.util.geomutil import *
from pylayers.util.plotutil import *
import numpy as np
import scipy as sp

vl = np.array([0,0,1])  # beam in z direction 
pl = np.array([1,0,0])  # polar along x

phi = np.pi/2  # beam in y direction 
tilt = 0       # no tilt
polar = 'H'    # polar H 

M = MAzTiltPol(vl,pl,phi,tilt,polar)
assert False
print(M)




