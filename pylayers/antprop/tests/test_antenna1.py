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


# An example of broadcasting

A.Ftheta = A.Ftheta*exp(2*1j*pi*freq*electricalDelay)

A.Fphi   = A.Fphi*exp(2*1j*pi*freq*electricalDelay)


dsf = 2

A.vshd(dsf)

tn  = []
tet = []
tep = []
te  = []
tmse = []
tl  = range(1,15)
for l in tl:
    print 'l : ',l
    A.C.s1tos2(l)
    errelTh,errelPh,errel = A.errel(l,20,dsf,typ='s2')
    tet.append(errelTh)
    tep.append(errelPh)
    te.append(errel)

line1 = plt.plot(array(tl),10*log10(array(tep)),'b')
line2 = plt.plot(array(tl),10*log10(array(tet)),'r')
line3 = plt.plot(array(tl),10*log10(array(te)),'g')
plt.xlabel('order l')
plt.ylabel(u'$\epsilon_{rel}$  (dB)',fontsize=18)
plt.title('Evolution of reconstruction relative error wrt order')
plt.legend((u'$\epsilon_{rel}^{\phi}$',u'$\epsilon_{rel}^{\\theta}$',u'$\epsilon_{rel}^{total}$'))
#plt.legend((line1,line2,line3),('a','b','c'))
plt.show()
#legend(('errel_phi','errel_theta','errel'))
