from pylayers.antprop.antenna import *
import numpy as np
import matplotlib.pylab as plt
A = Antenna('vsh3','defant.vsh3')
theta = np.linspace(0,np.pi,70)
phi = np.linspace(0,2*np.pi,180)
th = np.kron(theta,np.ones(len(phi)))
ph = np.kron(np.ones(len(theta)),phi)
Fth,Fph = A.Fsynth3(th,ph)
