from pylayers.util.plotutil import *
import numpy as np
theta = np.linspace(0,np.pi,90)
phi = np.linspace(0,2*np.pi,180)
rho = np.ones((len(theta),len(phi)))
fig=plt.figure()
pol3D(fig,rho,theta,phi)
plt.show()
