from pylayers.antprop.spharm import *
from scipy.special import lpmv
from numpy import *

"""
    theta and phi should have the same size

    This test  evaluates Legendre Polynomials for L=5 and M = 5

    lpmv is a scipy function implementing Legendre polyomials

    VW VW2 VW3 are 3 implementation of evaluation of vector spherical harmonics
"""

Nt = 90
Np = 180
th = linspace(0,pi,Nt)     # 1 x Nt 
ph = linspace(0,2*pi,Np)   # 1 x Np 

theta = kron(th,ones(Np))  # 1 x Nt*Np = 1 x Ndir
phi   = kron(ones(Nt), ph) # 1 x Nt*Np = 1 x Ndir

L = 5
M = 5
t = indexvsh(5)    # K(L,M) x 2 
l = t[:,0]         # K(L,M) x 1
m = t[:,1]         # K(L,M) x 1
x = -cos(theta)    # 1 x Ndir

Pmm1l, Pmp1l = AFLegendre(L, M, x)
LG = lpmv(m.reshape(21,1,1),l.reshape(1,21,1),x.reshape(1,1,16200))

V1,W1 = VW(l,m,x,phi,Pmm1l,Pmp1l)
V2,W2 = VW2(l,m,x,phi,Pmm1l,Pmp1l)
V3,W3 = VW3(l,m,theta,phi)

