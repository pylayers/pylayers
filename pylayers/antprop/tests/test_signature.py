import matplotlib.pyplot as plt
import numpy as np
from pylayers.gis.layout import *
from pylayers.antprop.signature import *

# load the layout graphs

L = Layout('defstr.str')
L.build()
tx = np.array([2, -1])
rx = np.array([1, 1])

#L = Layout('TA-Office.str')
#L.build()
#tx = np.array([20, 8])
#rx = np.array([35, 6])

#s1,s2 = L.signature(0,0)
S = Signatures(L, tx, rx)
s1 = S.get_sigarr()

# define a sequence of interactions
#sig = s1[:, 69:72]
#sig =np.array([[ 4. , 7.],  [ 1.,  2.]])

#s = Signature(sig)
#s.ev(L)

#L.showSig(s1)

#s.ev(L)
#M, Y = s.sig2ray(L, tx, rx)

#M = s.image(tx)
#Y = s.backtrace(tx,rx,M)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#l1 = ax.plot(tx[0], tx[1], 'or')
#l2 = ax.plot(rx[0], rx[1], 'og')
#l3 = ax.plot(M[0, :], M[1, :], 'ob')
#if Y is not None:
    #l4 = ax.plot(Y[0, :], Y[1, :], 'xk')
    #ray = np.hstack((np.hstack((tx.reshape(2, 1), Y)), rx.reshape(2, 1)))
    #l5 = ax.plot(ray[0, :], ray[1, :], color='#999999', alpha=0.6,
        #linewidth=0.6)
#fig, ax = L.showGs(fig, ax)
#plt.show()
#
Mi, Yi = S.sigs2rays(L, tx, rx, s1)
