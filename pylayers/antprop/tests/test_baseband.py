import numpy as np
from pylayers.simul.link import *
import matplotlib.pyplot as plt
import time
import pdb 
print("=======================")
print(" start test_baseband.py ")
print("=======================")


# Only 1 frequency point for the propagation channel
force = True
WMHz = 100
Nf   = 1024
fcGHz = 5
fGHz1 = np.array([fcGHz])
fGHz2 = np.linspace(fcGHz-WMHz*1e-3/2.,fcGHz+WMHz*1e-3/2,Nf)
DL1=DLink(L='defstr.lay',fGHz=fGHz1)
DL1.a=np.array([2,2,1.0])
DL1.b=np.array([8,4,1.5])
DL2=DLink(L='defstr.lay',fGHz=fGHz2)
DL2.a=np.array([2,3.5,1.0])
DL2.b=np.array([8,3.5,1.5])
DL1.eval(diffraction=True,force=force,cutoff=6,ra_vectorized=True,applywav=False)
DL2.eval(diffraction=True,force=force,cutoff=6,ra_vectorized=True,applywav=False)
H1 = DL1.H.baseband(fcGHz=fcGHz,WMHz=WMHz,Nf=Nf)
H2 = DL2.H.baseband(fcGHz=fcGHz,WMHz=WMHz,Nf=Nf)

plt.plot(np.real(H1.y[0,0,:]),np.imag(H1.y[0,0,:]),'ob')
plt.plot(np.real(H2.y[0,0,:]),np.imag(H2.y[0,0,:]),'or')
