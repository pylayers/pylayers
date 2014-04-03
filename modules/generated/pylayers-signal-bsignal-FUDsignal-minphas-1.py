from pylayers.signal.bsignal import *
import numpy as np
fGHz = np.arange(2,11,0.1)
tau1 = np.array([1,2,3])[:,np.newaxis]
y = np.exp(-2*1j*np.pi*fGHz[np.newaxis,:]*tau1)/fGHz[np.newaxis,:]
H = FUDsignal(x=fGHz,y=y,taud=np.array([15,17,18]))
f,a = H.plot(typ=['ru'])
H.minphas()
H.taue
# array([ 1.,  2.,  3.])
f,a = H.plot(typ=['ru'])
