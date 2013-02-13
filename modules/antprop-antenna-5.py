import matplotlib.pyplot as plt
from pylayers.antprop.antenna import *
ifreq = 0
A     = Antenna('vsh3','defant.vsh3')
A.Nt  = 30
A.Np  = 60
A.Nf  = len(A.fa)
theta = np.linspace(0,np.pi,A.Nt)
phi   = np.linspace(0,2*np.pi,A.Np)
th    = np.kron(theta,np.ones(A.Np))
ph    = np.kron(np.ones(A.Nt),phi)
Fth3,Fph3 = A.Fsynth3(th,ph)
FTh3 = Fth3.reshape(A.Nf,A.Nt,A.Np)
FPh3 = Fph3.reshape(A.Nf,A.Nt,A.Np)
show3D(FTh3,theta,phi,ifreq)
txt = plt.title('show3D example')
plt.show()
