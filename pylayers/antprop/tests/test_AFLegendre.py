from pylayers.antprop.antenna import AFLegendre, AFLegendre2 ,VW,VW2,index_vsh
import numpy as np 
import time


L = 2
M = 2 
ind = index_vsh(L, M)
l = ind[:, 0]
m = ind[:, 1]
x = np.arange(0,1,0.5)
phi = np.linspace(0,2*np.pi,len(x))

tic = time.clock()
Pmm1l, Pmp1l = AFLegendre2(L,M,x)
toc = time.clock()
print toc -tic 

tic = time.clock()
Pmm1n, Pmp1n = AFLegendre(L,M,x)
toc = time.clock()
print toc-tic 

d1 = Pmm1l - Pmm1n
d2 = Pmp1l - Pmp1n


tic = time.clock()
V,W = VW(l,m,phi,x,Pmm1l,Pmp1l)
toc = time.clock()
print "Old : ", toc-tic 

tic = time.clock()
V2,W2 = VW2(l,m,phi,x,Pmm1l,Pmp1l)
toc = time.clock()
print "New : ", toc-tic 

u = np.ravel(V-V2)
v = np.ravel(W-W2)
