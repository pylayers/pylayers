import matplotlib.pyplot as plt
import numpy as np
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
L = Layout('defstr.ini')
s = Signature(seq)
tx = np.array([760,1113])
rx = np.array([762,1114])
s.ev(L)
M = s.image(tx)
isvalid,Y = s.backtrace(tx,rx,M)

fig,ax = L.showG('s',labels=1,aw=1,axes=1)
l1 = ax.plot(tx[0],tx[1],'or')
l2 = ax.plot(rx[0],rx[1],'og')
l3 = ax.plot(M[0,:],M[1,:],'ob')
l4 = ax.plot(Y[0,:],Y[1,:],'xk')
ray = np.hstack((np.hstack((tx.reshape(2,1),Y)),rx.reshape(2,1)))
l5 = ax.plot(ray[0,:],ray[1,:],color='#999999',alpha=0.6,linewidth=0.6)

plt.show()
