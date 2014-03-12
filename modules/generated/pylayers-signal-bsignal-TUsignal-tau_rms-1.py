from pylayers.measures.mesuwb import *
import matplotlib.pyplot as plt
M = UWBMeasure(1)
ch4 = M.tdd.ch4
f1,a1=ch4.plot(color='k')
plt.title("WHERE1 M1 UWB Channel impulse response")
f2,a2=ch4.plot(color='k')
plt.title("WHERE1 M1 UWB Channel impulse response")
ax1=plt.axis([10,160,-90,-50])
f3,a3=ch4.plot(color='k')
plt.title("WHERE1 M1 UWB Channel impulse response")
ax2=plt.axis([20,120,-80,-50])
