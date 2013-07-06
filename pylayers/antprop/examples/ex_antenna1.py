from pylayers.antprop.antenna import *
from pylayers.antprop.antvsh import *
import matplotlib.pylab as plt
from numpy import *
import pdb
"""
This test : 
    
    1 : loads a measured antenna
    2 : applies an electrical delay obtained from data with getdelay method
    3 : evaluate the antenna vsh coefficient with a downsampling factor of 2
    4 : evaluates the relative error of reconstruction (vsh3) for various values of order l
    5 : display the results

"""


filename = 'S1R6.mat'

A = Antenna(filename,'ant/UWBAN/Matfile')

#plot(freq,angle(A.Ftheta[:,maxPowerInd[1],maxPowerInd[2]]*exp(2j*pi*freq.reshape(len(freq))*electricalDelay)))
freq = A.fa.reshape(104,1,1)
delayCandidates = arange(-10,10,0.001)
electricalDelay = A.getdelay(freq,delayCandidates)
disp('Electrical Delay = ' + str(electricalDelay)+' ns') 


A.Ftheta = A.Ftheta*exp(2*1j*pi*freq*electricalDelay)
A.Fphi   = A.Fphi*exp(2*1j*pi*freq*electricalDelay)


dsf = 2
A = vsh(A,dsf)

tn  = []
tet = []
tep = []
te  = []
tmse = []
Lmax = 40
A.C.s1tos2(40)
A.C.s2tos3(1e-5)
errelTh,errelPh,errel = A.errel(20,dsf,typ='s3')

#tl = arange(1,6,0.01)
#for l in tl:
#    thresh = 10**(-l)
#    A.C.s2tos3(thresh)
#    print errel
#    tet.append(errelTh)
#    tep.append(errelPh)
#    te.append(errel)
##
##line1 = plt.plot(array(tl),10*log10(array(tep)),'b')
##line2 = plt.plot(array(tl),10*log10(array(tet)),'r')
##line3 = plt.plot(array(tl),10*log10(array(te)),'g')
##plt.xlabel('threshold')
##plt.ylabel(u'$\epsilon_{rel}$  (dB)',fontsize=18)
##plt.title('Evolution of reconstruction relative error wrt order')
##plt.legend((u'$\epsilon_{rel}^{\phi}$',u'$\epsilon_{rel}^{\\theta}$',u'$\epsilon_{rel}^{total}$'))
###plt.legend((line1,line2,line3),('a','b','c'))
##plt.show()
###legend(('errel_phi','errel_theta','errel'))
