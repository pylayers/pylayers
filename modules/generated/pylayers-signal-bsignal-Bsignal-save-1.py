from pylayers.signal.bsignal import *
import matplotlib.pyplot as plt
e = EnImpulse(fe=100)
fig,ax = e.plot(typ=['v'])
plt.title('original waveform')
e.save('impulse.mat')
del e
h = TUsignal()
h.load('impulse.mat')
fig,ax = h.plot(typ=['v'])
plt.title('retrieved waveform')
