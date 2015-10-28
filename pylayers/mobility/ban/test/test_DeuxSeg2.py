# -*- coding:Utf-8 -*-

import numpy as np
from pylayers.mobility.ban.DeuxSeg import *





N=5
M=10

A = 50*np.random.rand(3,N)
B = 5*np.random.rand(3,N)
C = 20*np.random.rand(3,M)
D = 120*np.random.rand(3,M)

# #3 x N x M
# AC = C[:,np.newaxis,:]-A[:,:,np.newaxis]
# # 3 x M 
# CD = D-C
# # 3 x N 
# BA = A-B

# #u0 : N x M
# u0 = np.einsum('ijk,ijk->jk',AC,AC)#np.dot(AC,AC)
# #u4 : N 
# u4 = np.einsum('ij,ij->j',BA,BA)[:,np.newaxis]#np.dot(BA,BA)
# # u5 : M 
# u5 = np.einsum('ij,ij->j',CD,CD)[np.newaxis,:]#np.dot(CD,CD)
# # u1 : N x M
# u1 = np.einsum('ij,ijk->jk',BA,AC)#np.dot(BA,AC)
# # u2 : N x M
# u2 = np.einsum('ik,ijk->jk',CD,AC)#np.dot(CD,AC)
# # u3 : N x M
# u3 = np.einsum('ik,ij->jk',CD,BA)#np.dot(CD,BA) 

# # den : N x M
# den   = u4*u5-u3*u3
# #alpha = N x M
# alpha = (u2*u3-u1*u5)/(1.*den)
# # beta = N x M
# beta  = (u1*u3-u2*u4)/(1.*den)
# # dmin : N x M
# dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5) 

# f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
# X  = A[:,:,np.newaxis]-alpha[np.newaxis,:,:]*BA[:,:,np.newaxis] # A - alpha*BA
# Y  = C[:,np.newaxis,:] + beta[np.newaxis,:,:]*CD[:,np.newaxis,:]# C + beta*CD

# g =np.einsum('ijk,ijk->jk',X-Y,X-Y)


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
