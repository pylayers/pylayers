from pylayers.antprop.antenna import *
from scipy.special import lpmv
from numpy import *

"""
    theta and phi should have the same size

"""

Nt = 90
Np = 180
th = linspace(0,pi,Nt)
ph = linspace(0,2*pi,Np)

theta = kron(th,ones(Np))
phi   = kron(ones(Nt), ph)

L = 5
M = 5
t = indexvsh(5)
l = t[:,0]
m = t[:,1]
x = -cos(theta)

Pmm1l, Pmp1l = AFLegendre(L, M, x)
LG = lpmv(m.reshape(21,1,1),l.reshape(1,21,1),x.reshape(1,1,16200))
#V1,W1 = VW(l,m,x,phi,Pmm1l,Pmp1l)
V2,W2 = VW2(l,m,x,phi,Pmm1l,Pmp1l)
V3,W3 = VW3(l,m,theta,phi)

