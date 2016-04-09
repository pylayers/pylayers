from pylayers.measures.mesuwb import *
from pylayers.simul.link import *
from pylayers.signal.waveform import *
import pdb

L = Layout('WHERE1.ini')
tx_id = 59
rx_id = [1,2,3,4]
M = UWBMeasure(tx_id)
Tx = M.tx
Rx1 = M.rx[1]
Rx2 = M.rx[2]
Rx3 = M.rx[3]
Rx4 = M.rx[4]
mex1 = M.tdd.ch1.x
mey1 = M.tdd.ch1.y
mex2 = M.tdd.ch2.x
mey2 = M.tdd.ch2.y
mex3 = M.tdd.ch3.x
mey3 = M.tdd.ch3.y
mex4 = M.tdd.ch4.x
mey4 = M.tdd.ch4.y
Lk1 = DLink(L=L,a=Tx,b=Rx1,cutoff=4,verbose=False)
Lk2 = DLink(L=L,a=Tx,b=Rx2,cutoff=4,verbose=False)
Lk3 = DLink(L=L,a=Tx,b=Rx3,cutoff=4,verbose=False)
Lk4 = DLink(L=L,a=Tx,b=Rx4,cutoff=4,verbose=False)
Lk1.Aa = Antenna('defant.vsh3')
Lk1.Ab = Antenna('defant.vsh3')
Lk2.Aa = Antenna('defant.vsh3')
Lk2.Ab = Antenna('defant.vsh3')
Lk3.Aa = Antenna('defant.vsh3')
Lk3.Ab = Antenna('defant.vsh3')
Lk4.Aa = Antenna('defant.vsh3')
Lk4.Ab = Antenna('defant.vsh3')
Lk1.eval(alg=5)
Lk2.eval(alg=5)
Lk3.eval(alg=5)
Lk4.eval(alg=5)
wav = Waveform(typ='W1compensate')
ir1 = Lk1.H.applywavB(wav.sf)
ir2 = Lk2.H.applywavB(wav.sf)
ir3 = Lk3.H.applywavB(wav.sf)
ir4 = Lk4.H.applywavB(wav.sf)
plt.subplot(411)
plt.semilogy(mex1,mey1[0,:],'b',linewidth=2,alpha=0.3)
plt.semilogy(ir1.x,ir1.y[0,:],'r')
plt.xlim(0,100)
plt.ylim(1e-6,1e-2)
plt.subplot(412)
plt.semilogy(mex2,mey2[0,:],'b',linewidth=2,alpha=0.3)
plt.semilogy(ir2.x,ir2.y[0,:],'r')
plt.xlim(0,100)
plt.ylim(1e-6,1e-2)
plt.subplot(413)
plt.semilogy(mex3,mey3[0,:],'b',linewidth=2,alpha=0.3)
plt.semilogy(ir3.x,ir3.y[0,:],'r')
plt.xlim(0,100)
plt.ylim(1e-6,1e-2)
plt.subplot(414)
plt.semilogy(mex4,mey4[0,:],'b',linewidth=2,alpha=0.3)
plt.semilogy(ir4.x,ir4.y[0,:],'r')
plt.xlim(0,100)
plt.ylim(1e-6,1e-2)
plt.show()
