from pylayers.measures.mesuwb import *
import matplotlib.pyplot as plt
M = UWBMeasure(1)
ch4 = M.tdd.ch4
ch4.plot(color='k')
plt.title("WHERE1 M1 UWB Channel impulse response")
ch4.plot(color='k')
plt.title("WHERE1 M1 UWB Channel impulse response")
ax1=plt.axis([10,160,-90,-50])
ch4.plot(color='k')
plt.title("WHERE1 M1 UWB Channel impulse response")
ax2=plt.axis([20,120,-80,-50])
plt.show()
tau_moy = ch4.tau_moy()
print "tau_moy: %2.2f ns" % tau_moy
tau_rms = ch4.tau_rms()
print "tau_rms: %2.2f" % tau_rms
