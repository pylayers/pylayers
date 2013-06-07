import pickle
import numpy as np
import  matplotlib.pyplot as plt 
from pylayers.signal.bsignal import *
#
file=open("tcir2.pickle","r")
tcir=pickle.load(file)
file.close()
#
for i in tcir[1].keys():
    cir = tcir[1][i]
    cir.zlr(0,200)
    try:
        ttcir=np.vstack((ttcir,cir.y))
    except:
        ttcir=cir.y
#    ttcir.append(cir)
#attcir=np.array(ttcir)
plt.imshow(20*np.log10(ttcir),vmax=-40,vmin=-150,origin='lower')
plt.colorbar()
plt.axis('tight')
#
