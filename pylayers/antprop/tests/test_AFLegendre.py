from pylayers.antprop.antenna import AFLegendre, AFLegendre2, AFLegendre3, VW,VW2,index_vsh
import numpy as np 
import time


L = 2
M = 2 
ind = index_vsh(L, M)
l = ind[:, 0]
m = ind[:, 1]
x = np.arange(0,1,0.0001)
phi = np.linspace(0,2*np.pi,len(x))

#tic = time.clock()
Pmm1l, Pmp1l = AFLegendre2(L,M,x)
#toc = time.clock()
#print toc -tic 

#tic = time.clock()
Pmm1n, Pmp1n = AFLegendre(L,M,x)
#toc = time.clock()
#print toc-tic 
P, Q = AFLegendre3(L,M,x)

d1 = Pmm1l - Pmm1n
d2 = Pmp1l - Pmp1n
d3 = Pmm1l - P
d4 = Pmp1l - Q
print d1,d2,d3,d4


