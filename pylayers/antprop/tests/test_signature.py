import matplotlib.pyplot as plt
import numpy as np
from pylayers.gis.layout import *
from pylayers.antprop.signature import *

# load the layout graphs

L = Layout('defstr.str')
#L.loadstr('defstr.str')
#L.buildGt()
#L.buildGr()
L.build()
s1,s2 = L.signature(0,0)
# define a sequence of interactions

sig = s1[:,3:4]

s = Signature(sig)

#L.showSig(s1)

s.ev(L)
tx = np.array([4,-1])
rx = np.array([1,1])
M, Y = s.sig2ray(L, tx, rx)

#M = s.image(tx)
#Y = s.backtrace(tx,rx,M)
fig = plt.figure()
ax = fig.add_subplot(111)
l1 = ax.plot(tx[0],tx[1],'or')
l2 = ax.plot(rx[0],rx[1],'og')
l3 = ax.plot(M[0,:],M[1,:],'ob')
if Y<> None:
    l4 = ax.plot(Y[0,:],Y[1,:],'xk')
    ray = np.hstack((np.hstack((tx.reshape(2,1),Y)),rx.reshape(2,1)))
    l5 = ax.plot(ray[0,:],ray[1,:],color='#999999',alpha=0.6,linewidth=0.6)
fig,ax = L.showGs(fig,ax)
plt.show()

"""

pa = self.pa
pb = self.pb
typ = self.typ

N = np.shape(pa)[1]
I2 = np.eye(2)
z0 = np.zeros((2, 1))

pkm1 = rx.reshape(2, 1)
Y = pkm1
k = 0
beta = .5
while (((beta <= 1) & (beta >= 0)) & (k < N)):
l0 = np.hstack((I2, pkm1 - M[:, -(k + 1)].reshape(2, 1), z0))
l1 = np.hstack((I2, z0, pa[:, -(k + 1)].reshape(
2, 1) - pb[:, -(k + 1)].reshape(2, 1)))
T = np.vstack((l0, l1))
yk = np.hstack((pkm1[:, 0].T, pa[:, -(k + 1)].T))
xk = la.solve(T, yk)
pkm1 = xk[0:2].reshape(2, 1)
gk = xk[2::]
beta = gk[1]

"""
