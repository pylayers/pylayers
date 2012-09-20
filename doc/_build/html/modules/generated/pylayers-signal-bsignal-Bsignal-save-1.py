from pylayers.signal.bsignal import *
import matplotlib.pyplot as plt
e = EnImpulse()
e.plot()
e.save('impulse.mat')
del e
h = TUsignal()
h.load('impulse.mat')
h.plot()
