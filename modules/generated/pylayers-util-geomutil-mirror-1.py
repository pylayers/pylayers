from pylayers.util.geomutil import *
from pylayers.util.plotutil import *
import matplotlib.pyplot as plt
import numpy as np
p = np.random.randn(2,1000)
pa  = np.array([0,0])
pb  = np.array([0,1])
M = mirror(p,pa,pb)
plt.plot(p[0,:],p[1,:],'or',alpha=0.2)
plt.plot(M[0,:],M[1,:],'ob',alpha=0.2)
