import matplotlib.pyplot as plt
import numpy as np
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
L = Layout()
L.buildGt()
L.buildGr()
seq = [1,5,1]
s = Signature(seq)
tx = np.array([4,-1])
rx = np.array([1,1])
s.ev(L)
M = s.image(tx)
Y = s.backtrace(tx,rx,M)
fig = plt.figure()
ax = fig.add_subplot(111)
l1 = ax.plot(tx[0],tx[1],'or')
l2 = ax.plot(rx[0],rx[1],'og')
l3 = ax.plot(M[0,:],M[1,:],'ob')
l4 = ax.plot(Y[0,:],Y[1,:],'xk')
ray = np.hstack((np.hstack((tx.reshape(2,1),Y)),rx.reshape(2,1)))
l5 = ax.plot(ray[0,:],ray[1,:],color='#999999',alpha=0.6,linewidth=0.6)
fig,ax = L.showGs(fig,ax)
plt.show()
