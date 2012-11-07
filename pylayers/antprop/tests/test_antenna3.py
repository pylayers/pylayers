from pylayers.antprop.antenna import *
import matplotlib.pylab as plt 
from numpy import *
import pdb
#
# Visalization of coefficient in different shapes
#
filename = 'S1R1.mat'

A = Antenna('mat',filename,'ant/UWBAN/Matfile')

#plot(freq,angle(A.Ftheta[:,maxPowerInd[1],maxPowerInd[2]]*exp(2j*pi*freq.reshape(len(freq))*electricalDelay)))
freq = A.fa.reshape(104,1,1)/1e9
delayCandidates = arange(-10,10,0.001)
electricalDelay = A.getdelay(freq,delayCandidates)
disp('Electrical Delay = ' + str(electricalDelay)+' ns') 


A.Ftheta = A.Ftheta*exp(2*1j*pi*freq*electricalDelay)
A.Fphi   = A.Fphi*exp(2*1j*pi*freq*electricalDelay)

dsf = 2
#
# Calculate Vector Spherical Harmonics
#
A.vshd(dsf)
A.C.info()
A.C.s1tos2(15)
A.C.info()
A.C.s2tos3(1e-2)
A.C.info()
A.C.show(typ='s2')
plt.show()
A.C.show(typ='s3')
plt.show()
