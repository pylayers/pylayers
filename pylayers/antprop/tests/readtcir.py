import pickle
import numpy as np
import  matplotlib.pyplot as plt 
from pylayers.signal.bsignal import *
#
file=open("tcir5.pickle","r")
tcir=pickle.load(file)
file.close()
#
for i in tcir[1].keys():
    cir = tcir[1][i]
    cir.zlr(0,150)
    try:
        ttcir=np.vstack((ttcir,cir.y))
    except:
        ttcir=cir.y
#    ttcir.append(cir)
#attcir=np.array(ttcir)
dmax=150*0.3
plt.imshow(20*np.log10(ttcir),vmax=-40,vmin=-150,origin='lower',extent=[0,dmax,1,69])
plt.xlabel(r'delay $\times$ c (meters)',fontsize=20)
#plt.ylabel(r'distance along trajectory (meters)',fontsize=20)
plt.ylabel(r'trajectory index number',fontsize=20)
clb=plt.colorbar()
clb.set_label('level (dB)',fontsize=20)

plt.axis('tight')
#
