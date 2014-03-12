from pylayers.measures.mesuwb import *
M1 = UWBmeasure(1)
ch4 = M.tdd.ch4
ch4.plot(color='k')
ch4.plot(color='k')
plt.axis([10,160,-90,-50])
ch4.plot((color='k')
plt.axis([20,120,-80,-50])
tau_moy = ch4.tau_moy()
print "tau_moy: %2.2f" % tau_moy
# 38.09
tau_rms = ch4.tau_rms()
print "tau_rms: %2.2f" % tau_rms
# 13.79
