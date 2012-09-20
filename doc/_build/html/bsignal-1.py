from pylayers.signal.bsignal import *
from pylayers.simul.simulem import *
from matplotlib.pylab import *
fc     = 4
band   = 2
thresh = 10
fe     = 100
ip     = EnImpulse([],fc,band,thresh,fe)
ip.plot()
show()