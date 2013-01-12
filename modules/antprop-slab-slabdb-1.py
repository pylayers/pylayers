from pylayers.antprop.slab import *
import numpy as np
import matplotlib.pyplot as plt
theta = np.arange(0,np.pi/2,0.01)
fGHz  = np.arange(0.1,10,0.2)
sl  = SlabDB('matDB.ini','slabDB.ini')
mat = sl.mat
air = mat['AIR']
brick  = mat['BRICK']
II  = MatInterface([air,brick],0,fGHz,theta)
II.RT()
fig = plt.figure()
II.plotwrta(0)
plt.show()
