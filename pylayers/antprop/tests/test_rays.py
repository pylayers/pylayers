from pylayers.antprop.rays import *
from pylayers.simul.simulem import *
from numpy import *

# load default simuluation
S     = Simul()
# get the tud file for link 1 -> 1
grTud = GrRayTud()
grTud.load(S.dtud[1][1],S.dtang[1][1],S.drang[1][1],S.L.sl)
# get the first ray of the tay cluster 
g0 = grTud.rayTud[0]
g0.info()
g0.inter
g0.inter[0]
# get the first 3 interactions 
i0 = g0.inter[0]
i1 = g0.inter[1]
i2 = g0.inter[2]

I = Interactions()
I.C = np.array([[1,-1],[2,3]]).reshape(1,1,2,2)
M = np.array([[1,2],[3,4]])
I.addB(M)
U  = np.eye(2).reshape(1,1,2,2)
