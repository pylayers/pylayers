from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import networkx as nx
from IPython.display import Image,HTML,Latex

L=Layout('DLR.osm')
#L.build()
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()

plt.ion()
L.showG('t',labels=True,figsize=(10,10))

nc1 = 1
nc2 = 14
S=Signatures(L,nc1,nc2)
S.run1(cutoff=4)


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


r=S.rays(ptx,prx)
figsize=(10,10)
r.show(L,i=[3],r=[0,1])
r3 = r.to3D()
r3.locbas(L)
r3.fillinter(L)
r3.show(L,i=[3],r=[2])
pg=np.array([[760.59],[1124.911],[0]])
r3.show3(strucname='DLR',pg=pg)

