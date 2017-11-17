import numpy as np 
import matplotlib.pyplot as plt
from pylayers.antprop.antssh import *
nth = 90
nph = 180
L  = 20
theta = np.linspace(0,np.pi,nth)
phi = np.linspace(0,2*np.pi,nph)
Y,idx = SSHFunc(L,theta,phi)
plt.subplot(121)
plt.imshow(np.abs(Y.T),cmap='jet')
plt.axis('auto')
plt.subplot(122)
plt.imshow(np.angle(Y.T),cmap='jet')
plt.axis('auto')

