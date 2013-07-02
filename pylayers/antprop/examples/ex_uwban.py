from pylayers.antprop.antenna import *
from pylayers.antprop.spharm import *
from pylayers.antprop.antvsh import *
from pylayers.util.pyutil import *
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.pyplot as plt
import os

_filename = 'S1R2.mat'
A = Antenna(_filename,'ant/UWBAN/Matfile')
filename=getlong(_filename,'ant/UWBAN/Matfile')
Norig = os.path.getsize(filename)
freq = A.fa.reshape(104,1,1)
ed = A.getdelay(freq)
A.Ftheta = A.Ftheta*exp(2*1j*pi*freq*ed)
A.Fphi   = A.Fphi*exp(2*1j*pi*freq*ed)


A = vsh(A,dsf=2)
A.C.s1tos2(20)
A.C.s2tos3(1e-5)
A.savevsh3()
filevsh3 = getlong(_filename.replace('.mat','.vsh3'),'ant')
Nvsh3 =  os.path.getsize(filevsh3)
ratio = Norig/(1.*Nvsh3)
et1,et2,et3 =A.errel(dsf=1,typ='s3')
et3l = 10*log10(et3)
Nc = len(A.C.Br.ind3)
Nf = A.Nf
csize = 4*Nc*Nf
ch1 = _filename.replace('.mat','')
ch2 = ', Nf ='+str(Nf)
ch3 = ', ['+str(A.fa[0])+','+str(A.fa[-1])+' ] GHz'
ch4 = ', Nc = '+ str(Nc)
ch5 = ', size = '+ str(csize)
ch6 = ', compress = '+ str(ratio)[0:5]
ch7 = ', relative error ='+str(et3l)[0:5]+' dB'
fig1 = plt.figure(figsize=(10,10))
A.C.plot(subp=False,titre=ch1+ch2+ch3+ch4+ch5+ch6)
fig1.savefig(_filename.replace('.mat','')+'-pl.png')
fig2 = plt.figure(figsize=(10,10))
A.C.show(typ='s3',titre=ch1+ch2+ch3+ch4+ch5+ch6)
fig2.savefig(_filename.replace('.mat','')+'-im.png')
