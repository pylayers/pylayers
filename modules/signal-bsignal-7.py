from pylayers.signal import *
from matplotlib.pylab import *
ip = EnImpulse()
fig,ax = ip.plot(types=['v'])
ip.zlr(-10,10)
fig,ax = ip.plot(types=['v'],fig=fig,ax=ax)
show()
