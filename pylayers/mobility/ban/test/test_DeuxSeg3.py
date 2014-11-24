# -*- coding:Utf-8 -*-

import numpy as np
from pylayers.mobility.ban.DeuxSeg import *




# link
N=4
# body cylinders
M=10
# time
K=5

A = 50*np.random.rand(3,N,K)
B = 5*np.random.rand(3,N,K)
C = 20*np.random.rand(3,M,K)
D = 120*np.random.rand(3,M,K)

# # 3 x N x M x K
# AC = C[:,np.newaxis,:]-A[:,:,np.newaxis]
# # 3 x M 
# CD = D-C
# # 3 x N 
# BA = A-B

# # u0 : N x M
# u0 = np.einsum('ijk...,ijk...->jk...',AC,AC)#np.dot(AC,AC)
# # u4 : N 
# u4 = np.einsum('ij...,ij...->j...',BA,BA)[:,np.newaxis]#np.dot(BA,BA)
# # u5 : M 
# u5 = np.einsum('ij...,ij...->j...',CD,CD)[np.newaxis,:]#np.dot(CD,CD)
# # u1 : N x M
# u1 = np.einsum('ij...,ijk...->jk...',BA,AC)#np.dot(BA,AC)
# # u2 : N x M
# u2 = np.einsum('ik...,ijk...->jk...',CD,AC)#np.dot(CD,AC)
# # u3 : N x M
# u3 = np.einsum('ik...,ij...->jk...',CD,BA)#np.dot(CD,BA) 

# # den : N x M
# den   = u4*u5-u3*u3
# # alpha = N x M
# alpha = (u2*u3-u1*u5)/(1.*den)
# # beta = N x M
# beta  = (u1*u3-u2*u4)/(1.*den)
# # dmin : N x M
# dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5) 

# f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
# X  = A[:,:,np.newaxis]-alpha[np.newaxis,:,:]*BA[:,:,np.newaxis] # A - alpha*BA
# Y  = C[:,np.newaxis,:] + beta[np.newaxis,:,:]*CD[:,np.newaxis,:]# C + beta*CD

# g =np.einsum('ijk...,ijk...->jk...',X-Y,X-Y)


# alpha,beta,dmin = dmin3d(A,B,C,D)
# f,g = dist(A,B,C,D,alpha,beta)

f,g,alpha,beta,dmin=segdist(A,B,C,D,hard=False)

## verif:
aa=np.empty((N,M,K))
bb=np.empty((N,M,K))
dd=np.empty((N,M,K))
ff=np.empty((N,M,K))
gg=np.empty((N,M,K))

for i in range(N):
    for j in range(M):
        for k in range(K):
            aa[i,j,k],bb[i,j,k],dd[i,j,k]=dmin3d_nonvectorized(A[:,i,k],B[:,i,k],C[:,j,k],D[:,j,k])
            ff[i,j,k],gg[i,j,k] = dist_nonvectorized(A[:,i,k],B[:,i,k],C[:,j,k],D[:,j,k],aa[i,j,k],bb[i,j,k])

assert (aa-alpha).any()==0
assert (bb-beta).any()==0
assert (dd-dmin).any()==0
assert (ff-f).any()==0
assert (gg-g).any()==0

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
