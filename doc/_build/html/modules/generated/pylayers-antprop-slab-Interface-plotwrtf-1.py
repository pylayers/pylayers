from pylayers.antprop.slab import *
import matplotlib.pylab as plt
import numpy as np
theta = np.arange(0,np.pi/2,0.01)
fGHz  = np.arange(0.1,10,0.2)
sl    = SlabDB('def.mat','def.slab')
mat   = sl.mat
lmat  = [mat['AIR'],mat['WOOD']]
II    = MatInterface(lmat,0,fGHz,theta)
II.RT()
fig = plt.figure()
II.plotwrtf(10)
plt.show()
