
# coding: utf-8

# This notebook is an implementation from the paper "Fast UTD Diffraction Coefficient Using Only One Suare Root" from Jean Francçois Legendre and Thierry Marsault.

# $N$ est compris entre 0 et 1. 0 est la lame de coupeau et 1 est le demi-plan 

# In[13]:

import numpy as np 
import matplotlib.pyplot as plt
from pylayers.antprop.slab import *


# In[4]:

Nf=1
Nr=1
fGHz = np.linspace(1,10,Nf)
N = 320/180.*np.ones(Nr)#320/180.
phi0 = np.ones(Nr)#np.linspace(0.01,2*np.pi-0.01,Nr)#40*np.pi/180.
phi = np.linspace(0.01,2*np.pi-0.01,Nr)
#phi = np.linspace(0,3*np.pi/2,10)
dm = MatDB()
mat0 = dm['METAL']
matN = dm['METAL']
si = 10000.*np.ones(Nr)
sd = 1.*np.ones(Nr)
beta=np.array([np.pi/2.])



# def quickD(fGHz,phi0,phi,si,sd,N,mat0,matN,beta=np.pi/2,debug=False):


def qsin(x):
    return x *(1.0- (0.1666666666*x**2) *(1.0- (0.05*x**2) *(1.* - (0.1428571429*x**2) )))
def qcos(x):
    return 1.0- (0.5*x**2) *(1.0- (0.0833333333*x**2)* (1. - (0.0333333333*x**2) *(1.0- (0.0178571428*x**2) )))



fGHz  = fGHz[:,None]
phi0  = phi0[None,:]
phi   = phi[None,:]
si    = si[None,:]
sd    = sd[None,:]
N     = N[None,:]
beta  = beta[None,:]

L     = si*sd/(si+sd)




# tho  = np.empty((fGHz.shape[0],phi.shape[1]))
# thn  = np.empty((fGHz.shape[0],phi.shape[1]))
# # PHI0 = phi0 * np.ones(phi.shape)
# # PHI  = np.ones(phi0.shape)*phi
# # BN   = np.ones(phi0.shape)*N



# c1 = phi>phi0
# c2 = ~c1
# tho[:,c1[0,:]] = phi0[:,c1[0,:]]
# thn[:,c1[0,:]] = N[:,c1[0,:]]*np.pi-phi[:,c1[0,:]]
# tho[:,c2[0,:]] = phi[:,c2[0,:]]
# thn[:,c2[0,:]] = N[:,c2[0,:]]*np.pi-phi0[:,c2[0,:]]



# er0  = np.real(mat0['epr'])
# err0 = np.imag(mat0['epr'])
# ur0  = np.real(mat0['mur'])
# urr0 = np.imag(mat0['mur'])
# sigma0 = mat0['sigma']
# deltah0 = mat0['roughness']

# erN  = np.real(matN['epr'])
# errN = np.imag(matN['epr'])
# urN  = np.real(mat0['mur'])
# urrN = np.imag(mat0['mur'])
# sigmaN = matN['sigma']
# deltahN = matN['roughness']





# In[7]:
k = 0.3/fGHz
alpha = N*np.sqrt(L*2*k)


def gfunc(xsi):
    sg = np.ones(xsi.shape)
    c1m = xsi<0
    c1p = xsi>(np.pi/4)
    xsi[c1m]=-xsi[c1m]
    xsi[c1p]=xsi[c1p]-np.pi/4
    sg[c1m]=-1.
    sg[c1p]=-1.
    y = alpha * np.sin(xsi)
    ny = y<0
    y[ny]=-y[ny]
    sg[ny]=-sg[ny]
    uy = (y**2)<0.005
    g = np.zeros(y.shape)
    g[uy] = -0.5 + 0.39384228*y[uy]*(1.+1j)
    g[~uy] = -0.1994711402/y[~uy] * ((0.5/(y[~uy]**2+1.))+1j*(0.5/(y[~uy]**2-1.)))
    g = g * np.cos(xsi) * sg * (1.+1j)

    return g
# In[10]:

xsi11 = (np.pi/4+phi-phi0)/(2*N)
xsi12 = (np.pi/4-phi+phi0)/(2*N)
xsi13 = (np.pi/4-phi-phi0)/(2*N)
xsi14 = (phi+phi0-(2*N-1)*np.pi/4)/(2*N)

Df1 = gfunc(xsi11)*np.sqrt(L)
Df2 = gfunc(xsi12)*np.sqrt(L)
Df3 = gfunc(xsi13)*np.sqrt(L)
Df4 = gfunc(xsi14)*np.sqrt(L)

# sg11 = np.ones(xsi11.shape)
# c11m = np.where(xsi11<0)
# c11p = np.where(xsi11>(np.pi/4))

# xsi11[c11m]=xsi11[c11m]*(-1)
# xsi11[c11p]=xsi11[c11p]-np.pi/4
# sg11[c11m]=-1.
# sg11[c11p]=-1.

# # In[124]:

# sg12 = np.ones(xsi12.shape)
# c12m = np.where(xsi12<0)
# c12p = np.where(xsi12>(np.pi/4))
# xsi12[c12m]=xsi12[c12m]*(-1)
# xsi12[c12p]=xsi12[c12p]-np.pi/4
# sg12[c12m]=-1.
# sg12[c12p]=-1.

# # In[125]:
# sg13 = np.ones(xsi13.shape)
# c13m = np.where(xsi13<0)
# c13p = np.where(xsi13>(np.pi/4))
# xsi13[c13m]=xsi13[c13m]*(-1)
# xsi13[c13p]=xsi13[c13p]-np.pi/4
# sg13[c13m]=-1.
# sg13[c13p]=-1.

# # In[126]:
# sg14 = np.ones(xsi14.shape)
# c14m = np.where(xsi14<0)
# xsi14[c14m]=xsi14[c14m]*(-1)
# c14p = np.where(xsi14>(np.pi/4))
# xsi14[c14p]=xsi14[c14p]-np.pi/4
# sg14[c14m]=-1.
# sg14[c14p]=-1.


# y11 = alpha * qsin(xsi11)
# y12 = alpha * qsin(xsi12)
# y13 = alpha * qsin(xsi13)
# y14 = alpha * qsin(xsi14)

# ny11 = np.where(y11<0)[0]
# ny12 = np.where(y12<0)[0]
# ny13 = np.where(y13<0)[0]
# ny14 = np.where(y14<0)[0]

# y11[ny11]=-y11[ny11]
# y12[ny11]=-y12[ny11]
# y13[ny11]=-y13[ny11]
# y14[ny11]=-y14[ny11]

# sg11[ny11]=-sg11[ny11]
# sg12[ny11]=-sg12[ny11]
# sg13[ny11]=-sg13[ny11]
# sg14[ny11]=-sg14[ny11]


# uy11 = y11**2<0.005
# uy12 = y12**2<0.005
# uy13 = y13**2<0.005
# uy14 = y14**2<0.005

# g11 = np.zeros(y11.shape)
# g12 = np.zeros(y12.shape)
# g13 = np.zeros(y13.shape)
# g14 = np.zeros(y14.shape)

# g11[uy11] = -0.5 + 0.39384228*y11[uy11]*(1.+1j)
# g11[~uy11] = -1994711402/y11[~uy11] * ((0.5/(y11[~uy11]**2+1.))+1j*(0.5/(y11[~uy11]**2+1.)))
# g12[uy12] = -0.5 + 0.39384228*y12[uy12]*(1.+1j)
# g12[~uy12] = -1994711402/y12[~uy12] * ((0.5/(y12[~uy12]**2+1.))+1j*(0.5/(y12[~uy12]**2+1.)))
# g13[uy13] = -0.5 + 0.39384228*y13[uy13]*(1.+1j)
# g13[~uy13] = -1994711402/y13[~uy13] * ((0.5/(y13[~uy13]**2+1.))+1j*(0.5/(y13[~uy13]**2+1.)))
# g14[uy14] = -0.5 + 0.39384228*y14[uy14]*(1.+1j)
# g14[~uy14] = -1994711402/y14[~uy14] * ((0.5/(y14[~uy14]**2+1.))+1j*(0.5/(y14[~uy14]**2+1.)))


# g11 = g11 * qcos(xsi11) * sg11 * (1.+1j)
# g12 = g12 * qcos(xsi12) * sg12 * (1.+1j)
# g13 = g13 * qcos(xsi13) * sg13 * (1.+1j)
# g14 = g14 * qcos(xsi14) * sg14 * (1.+1j)



# D1f = np.sqrt(L) * g11
# D2f = np.sqrt(L) * g12
# D2f = np.sqrt(L) * g13
# D4f = np.sqrt(L) * g14



