from pylayers.antprop.antenna import *
import numpy as np
import matplotlib.pylab as plt
A = Antenna('defant.vsh3')
theta = np.linspace(0,np.pi,70)
phi = np.linspace(0,2*np.pi,180)
F = A.Fsynth3(theta,phi,pattern=True)
