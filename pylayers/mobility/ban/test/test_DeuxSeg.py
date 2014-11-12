import numpy as np
from pylayers.mobility.ban.DeuxSeg import *



nseg=10

A = np.random.rand(3,nseg)
B = np.random.rand(3,nseg)
C = np.random.rand(3,nseg)
D = np.random.rand(3,nseg)

aa=[]
bb=[]
dd=[]
ff=[]
gg=[]
for i in range(nseg):
    a,b,d=dmin3d_old(A[:,i],B[:,i],C[:,i],D[:,i])
    f,g = dist_old(A[:,i],B[:,i],C[:,i],D[:,i],a,b)
    aa.append(a)
    bb.append(b)
    dd.append(d)
    ff.append(f)
    gg.append(g)

aa=np.array(aa)
bb=np.array(bb)
dd=np.array(dd)
ff=np.array(ff)
gg=np.array(gg)


a,b,d=dmin3d(A,B,C,D)
f,g = dist(A,B,C,D,a,b)

assert (aa - a).any()==0
assert (bb - b).any()==0
assert (dd - d).any()==0
assert (ff - f).any()==0
assert (gg - g).any()==0


###############
# DEBUG 1
###############

# np.random.seed(0)

# A = np.random.rand(3)
# B = np.random.rand(3)
# C = np.random.rand(3)
# D = np.random.rand(3)

# X=np.dot(A,B)
# Y=np.dot(C,D)

# A2=np.vstack((A,C)).T
# B2=np.vstack((B,D)).T

# XY=np.einsum('ij,ij->j',A2,B2)

# print 'X == XY[0]',np.allclose(X,XY[0])
# print 'Y == XY[1]',np.allclose(Y,XY[1])

# assert np.allclose(X,XY[0])
# assert np.allclose(Y,XY[1])


###############
# DEBUG 2
###############

# a,b,d=dmin3d(A,B,C,D)
# f,g = dist(A,B,C,D,a,b)


