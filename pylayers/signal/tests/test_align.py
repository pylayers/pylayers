from pylayers.util import project
from pylayers.signal.bsignal import *
import matplotlib.pylab as plt
import numpy as np
import pdb
import glob


plt.ion()
i1 = EnImpulse()
i2 = EnImpulse()
i2.translate(-10.4)
i3  = i1.align(i2)
#i3 = TUsignal()
#i3.x = L[0].x
#i3.y = np.vstack((L[0].y,L[1].y))
#i3.imshow()
#i3.flatteny(reversible=True)
#i3.plot(iy=0,col='r')
#i3.plot(iy=1,col='b')
#plt.show()

