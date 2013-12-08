import numpy as np
import matplotlib.pyplot as plt
from pylayers.antprop.slab import *
sl = SlabDB('matDB.ini','slabDB.ini')
lmatname = ['PLATRE-57GHz','AIR','PLATRE-57GHz']
lthick = [0.018,0.03,0.018]
sl.add('placo',lmatname,lthick)
theta = np.arange(0,np.pi/2,0.01)
fGHz = np.array([57.5])
sl['placo'].ev(fGHz,theta)
fig,ax=sl['placo'].plotwrt(var='a',types=['m'])
plt.show()
