from pylayers.antprop.slab import *
import matplotlib.pylab as plt
import numpy as np
theta = np.arange(0,np.pi/2,0.01)
fGHz = np.arange(0.1,10,0.2)
sl = SlabDB('matDB.ini','slabDB.ini')
mat   = sl.mat
lmat  = [mat['AIR'],mat['WOOD']]
II    = MatInterface(lmat,0,fGHz,theta)
II.RT()
fig,ax = II.plotwrt(var='a',kv=10,types=['m'])
air = mat['AIR']
brick  = mat['BRICK']
II  = MatInterface([air,brick],0,fGHz,theta)
II.RT()
fig,ax = II.plotwrt(var='f',color='k',types=['m'])
plt.show()
