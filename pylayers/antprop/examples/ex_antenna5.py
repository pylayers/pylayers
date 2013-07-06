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


filename = 'S1R1.mat'

A = Antenna(filename,'ant/UWBAN/Matfile')
B = Antenna(filename,'ant/UWBAN/Matfile')

#plot(freq,angle(A.Ftheta[:,maxPowerInd[1],maxPowerInd[2]]*exp(2j*pi*freq.reshape(len(freq))*electricalDelay)))
freq = A.fa.reshape(104,1,1)
delayCandidates = arange(-10,10,0.001)
electricalDelay = A.getdelay(freq,delayCandidates)
disp('Electrical Delay = ' + str(electricalDelay)+' ns') 


A.Ftheta = A.Ftheta*exp(2*1j*pi*freq*electricalDelay)
B.Ftheta = B.Ftheta*exp(2*1j*pi*freq*electricalDelay)
A.Fphi   = A.Fphi*exp(2*1j*pi*freq*electricalDelay)
B.Fphi   = B.Fphi*exp(2*1j*pi*freq*electricalDelay)


dsf = 2
A = vsh(A,dsf)
B = vsh(B,dsf)

tn  = []
tet = []
tep = []
te  = []
tmse = []
l  = 20
A.C.s1tos2(l)
B.C.s1tos2(l)
u = np.shape(A.C.Br.s2)
Nf = u[0]
Nk = u[1]
tr = np.arange(2,Nk)
A.C.s2tos3_new(Nk)
B.C.s2tos3(1e-6)

UA = np.sum(A.C.Cr.s3*np.conj(A.C.Cr.s3),axis=0)
UB = np.sum(B.C.Cr.s3*np.conj(B.C.Cr.s3),axis=0)
ua = A.C.Cr.ind3
ub = B.C.Cr.ind3

da ={}
db ={}

for k in range(Nk):
    da[str(ua[k])]=UA[k]
    db[str(ub[k])]=UB[k]

tu = []

for t in sort(da.keys()):
    tu.append(da[t] - db[t])

errelTha,errelPha,errela = A.errel(l,20,dsf,typ='s3')
errelThb,errelPhb,errelb = B.errel(l,20,dsf,typ='s3')

print "a: nok",errela,errelPha,errelTha
print "b: ok ",errelb,errelPhb,errelThb

for r in tr:
    E = A.C.s2tos3_new(r)
    errelTh,errelPh,errel = A.errel(l,20,dsf,typ='s3')
    print 'r : ',r,errel,E
    tet.append(errelTh)
    tep.append(errelPh)
    te.append(errel)
#
line1 = plt.plot(array(tr),10*log10(array(tep)),'b')
line2 = plt.plot(array(tr),10*log10(array(tet)),'r')
line3 = plt.plot(array(tr),10*log10(array(te)),'g')
#
plt.xlabel('order l')
plt.ylabel(u'$\epsilon_{rel}$  (dB)',fontsize=18)
plt.title('Evolution of reconstruction relative error wrt order')
plt.legend((u'$\epsilon_{rel}^{\phi}$',u'$\epsilon_{rel}^{\\theta}$',u'$\epsilon_{rel}^{total}$'))
plt.legend((line1,line2,line3),('a','b','c'))
plt.show()
plt.legend(('errel_phi','errel_theta','errel'))
