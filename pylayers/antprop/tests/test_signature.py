import matplotlib.pyplot as plt
import numpy as np
from pylayers.gis.layout import *
from pylayers.antprop.signature import *

# load the layout graphs

L = Layout('defstr.str')
L.build()
tx = np.array([8, -1])
rx = np.array([1, 1])

#L = Layout('TA-Office.str')
#L.build()
#tx = np.array([20, 8])
#rx = np.array([35, 6])


S = Signatures(L, tx, rx)
s1 = S.get_sigarr()
rays = S.sigs2rays(L, tx, rx, s1)
S.show_rays2D(L, rays, tx, rx)
plt.show()
