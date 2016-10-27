from pylayers.antprop.slab import *
import matplotlib.pylab as plt
import numpy as np
theta = np.arange(0,np.pi/2,0.01)
fGHz = np.arange(1,100,0.2)
sl = SlabDB('matDB.ini','slabDB.ini')
mat = sl.mat
lmatname = [mat['WOOD_ITUR'],mat['AIR_ITUR']]
II = MatInterface(lmatname,0,fGHz,theta)
II.RT()
fig,ax = II.plotwrt(var='a',kv=10,typ=['m'])
air = mat['AIR_ITUR']
wood = mat['WOOD_ITUR']
II = MatInterface([air,wood],0,fGHz,theta)
II.RT()
fig,ax = II.plotwrt(var='f',color='k',typ=['m'])
plt.ion()
plt.show()