import numpy as np
from pylayers.antprop.slab import *
Nf=3
Nr=10
Nb=5
fGHz = np.linspace(0,10,Nf)
N = np.linspace(1,10,Nb)#320/180.
phi0 = np.linspace(0.01,2*np.pi-0.01,Nr)#40*np.pi/180.
phi = np.linspace(0.01,2*np.pi-0.01,Nr)
dm = MatDB()
mat0 = dm['METAL']
matN = dm['METAL']
si = 10000.*np.ones(Nr)
sd = 1.*np.ones(Nr)
plt.ion()
Ds,Dh,D1,D2,D3,D4 = diff(fGHz,phi0,phi,si,sd,N,mat0,matN)
