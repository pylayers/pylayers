from pylayers.signal.bsignal import *
import numpy as np
x = np.arange(0,1,0.01)
y = np.sin(2*np.pi*x)
s = Bsignal(x,y)
su = s.extract(np.arange(4,20))
f,a = s.plot()
f,a = su.plot()
