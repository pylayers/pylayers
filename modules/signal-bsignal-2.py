from pylayers.signal.bsignal import *
import matplotlib.pyplot as plt
e = EnImpulse()
fig,ax = e.plot(types=['v'])
e.save('impulse.mat')
del e
h = TUsignal()
h.load('impulse.mat')
fig,ax = h.plot(types=['v'])
