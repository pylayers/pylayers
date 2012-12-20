import matplotlib.pyplot as plt
import numpy as np
from pylayers.gis.layout import *
from pylayers.antprop.signature import *

# load the layout graphs

L = Layout('defstr.str')
L.build()
tx = np.array([8., -1., 1.])
rx = np.array([1., 1., 2.])

#L = Layout('TA-Office.str')
#L.build()
#tx = np.array([20, 8, 1])
#rx = np.array([35, 6, 2])


S = Signatures(L, tx, rx)

s1 = S.get_sigslist(tx, rx)

rays2d = S.sigs2rays(s1)

rays3d = S.ray2D3D(rays2d)

#S.show3(rays=rays3d)



