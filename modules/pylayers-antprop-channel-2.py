from pylayers.signal.bsignal import *
import numpy as np
N = 20
fGHz = np.arange(1,3,1)
taud = np.sort(np.random.rand(N))
alpha = np.random.rand(N,len(fGHz))
#s = Tchannel(x=fGHz,y=alpha,tau=taud)
#s.plot3d()
