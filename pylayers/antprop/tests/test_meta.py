from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import networkx as nx
import numpy as np

L=Layout('DLR.ini')
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()
#L.build()
#L.dumpw()
nc1 = 1
nc2 = 15

poly1 = L.Gt.node[nc1]['polyg']
cp1 = poly1.centroid.xy

poly2 = L.Gt.node[nc2]['polyg']
cp2 = poly2.centroid.xy
ptx = np.array([cp1[0][0],cp1[1][0],1.5])
prx = np.array([cp2[0][0],cp2[1][0],1.5])
print ptx
print prx
d = np.sqrt(np.dot((ptx-prx),(ptx-prx)))
tau = d/0.3
print d,tau




S   = Signatures(L,nc1,nc2)
metasig = [S.meta()]
print "S.run"
S.run(metasig,cutoff=5)


#S.run(L,metasig,cutoff=3)
print "r = S.rays "
r = S.rays(ptx,prx)
print "r3 = r.to3D "
r3 = r.to3D()
print "r3.locbas "
r3.locbas(L)
print "r3.fillinter "
r3.fillinter(L)
r3.show(L)
plt.show()
