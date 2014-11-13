import numpy as np
from pylayers.mobility.ban.DeuxSeg import *



N=5
M=10

A = np.random.rand(3,N)
B = np.random.rand(3,N)

C = np.random.rand(3,M)
D = np.random.rand(3,M)


alpha,beta,dmin = dmin3d(A,B,C,D)
f,g = dist(A,B,C,D,alpha,beta)



## verif:
aa=np.empty((N,M))
bb=np.empty((N,M))
dd=np.empty((N,M))
ff=np.empty((N,M))
gg=np.empty((N,M))
for i in range(N):
    for j in range(M):
        aa[i,j],bb[i,j],dd[i,j]=dmin3d_nonvectorized(A[:,i],B[:,i],C[:,j],D[:,j])
        ff[i,j],gg[i,j] = dist_nonvectorized(A[:,i],B[:,i],C[:,j],D[:,j],aa[i,j],bb[i,j])

assert (aa-alpha).any()==0
assert (bb-beta).any()==0
assert (dd-dmin).any()==0
assert (ff-f).any()==0
assert (gg-g).any()==0

######################



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



# nseg1=5
# nseg2=10

# A = np.random.rand(3,nseg1)
# B = np.random.rand(3,nseg1)
# C = np.random.rand(3,nseg2)
# D = np.random.rand(3,nseg2)

# aa=[]
# bb=[]
# dd=[]
# ff=[]
# gg=[]
# for i in range(nseg):
#     a,b,d=dmin3d_old(A[:,i],B[:,i],C[:,i],D[:,i])
#     f,g = dist_old(A[:,i],B[:,i],C[:,i],D[:,i],a,b)
#     aa.append(a)
#     bb.append(b)
#     dd.append(d)
#     ff.append(f)
#     gg.append(g)

# aa=np.array(aa)
# bb=np.array(bb)
# dd=np.array(dd)
# ff=np.array(ff)
# gg=np.array(gg)

# import ipdb
# ipdb.set_trace()
# a,b,d=dmin3d(A,B,C,D)
# f,g = dist(A,B,C,D,a,b)