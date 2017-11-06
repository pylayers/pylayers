from pylayers.signal.bsignal import *
from matplotlib.pylab import *
ip = TUsignal()
ip.EnImpulse()
f,a = ip.plot(typ=['v'])
ip.zlr(-10,10)
f,a = ip.plot(typ=['v'])
