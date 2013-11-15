from pylayers.signal.bsignal import *
from matplotlib.pylab import *
ip = EnImpulse()
ip.translate(-10)
fig,ax=ip.plot(types=['v'])
show()
