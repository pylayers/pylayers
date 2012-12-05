from numpy import *
import pdb
import timeit
import time

def mul(A,k,axis=1):
    if k>0:
        Ak = A[:,k-1,:,:]
        return(einsum('ikl,iln->ikn',Ak,mul(A,k-1,axis=1)))
    else:
        return(eye(2).reshape(1,2,2))

A = array([[1,2],[3,4]]).reshape(1,1,2,2)
B = array([[1,-1],[2,3]]).reshape(1,1,2,2)
A = concatenate((A,B),axis=1)
N = 30
Nf = 3000
A = random.rand(Nf*N*2*2).reshape(Nf,N,2,2)
timeit.timeit()
start_time = time.clock()
C = mul(A,N)
print time.clock()-start_time
B = eye(2).reshape(1,2,2)
start_time = time.clock()
for k in range(N):
    B = einsum('ikl,ilm->ikm',A[:,k,:,:],B)
