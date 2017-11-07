from pylayers.signal.bsignal import *
import matplotlib.pyplot as plt
e = TUsignal()
e.EnImpulse(feGHz=100)
fig,ax = e.plot(typ=['v'])
tit1 = plt.title('original waveform')
e.save('impulse.mat')
del e
h = TUsignal()
h.load('impulse.mat')
fig,ax = h.plot(typ=['v'])
tit2 = plt.title('retrieved waveform')
