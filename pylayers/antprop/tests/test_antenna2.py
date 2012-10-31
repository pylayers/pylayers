from pylayers.antprop.antenna import *
import matplotlib.pylab as plt 
from numpy import *
import pdb

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

v  = np.abs(A.C.Br.s1)
u  = np.nonzero(v==v.max())

plt.figure(figsize=(15,15))
for l in range(16):
    plt.subplot(4,4,l+1)
    plt.plot(np.real(A.C.Br.s1[:,l,0]),np.imag(A.C.Br.s1[:,l,0]),'k')
    plt.plot(np.real(A.C.Br.s1[:,l,1]),np.imag(A.C.Br.s1[:,l,1]),'b')
    plt.plot(np.real(A.C.Br.s1[:,l,2]),np.imag(A.C.Br.s1[:,l,2]),'r')
    plt.plot(np.real(A.C.Br.s1[:,l,2]),np.imag(A.C.Br.s1[:,l,2]),'g')
    plt.axis([-0.6,0.6,-0.6,0.6])
    plt.title('l='+str(l))

plt.show()
