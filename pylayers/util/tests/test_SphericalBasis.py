#
# Test SphericalBasis
# BTB_tx 
#
from pylayers.util.geomutil import *


a_g = np.array([[np.pi/2,np.pi/2],[np.pi/2,np.pi/2],[np.pi/2,np.pi/2],[np.pi/2,np.pi/2]])
G = SphericalBasis(a_g)
assert(len(G.shape)==3)
assert(G.shape[0]==3)
assert(G.shape[1]==3)

# Identity test
T = np.eye(3)
a_l,R = BTB(a_g,T)  
assert np.allclose(a_l,a_g)

