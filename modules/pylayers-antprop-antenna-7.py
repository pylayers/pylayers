from pylayers.antprop.antenna import *
import  matplotlib.pyplot as plt
import numpy as np
n=5
m=3
theta = np.linspace(0,np.pi,30)
phi   = np.linspace(0,2*np.pi,60)
plotVW(n,m,theta,phi)
