import matplotlib.pyplot as plt
from pylayers.antprop.antenna import *
ifreq = 0
A     = Antenna('defant.vsh3')
A.Nt  = 30
A.Np  = 60
A.Nf  = len(A.fa)
theta = np.linspace(0,np.pi,A.Nt)
phi   = np.linspace(0,2*np.pi,A.Np)
Fth3,Fph3 = A.Fsynth3(theta,phi)
FTh3 = Fth3.reshape(A.Nf,A.Nt,A.Np)
FPh3 = Fph3.reshape(A.Nf,A.Nt,A.Np)
show3D(FTh3,theta,phi,ifreq)
txt = plt.title('show3D example')
