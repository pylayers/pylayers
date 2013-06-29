from pylayers.antprop.rays import *
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import pylayers.signal.bsignal as bs
from pylayers.simul.simulem import *

S = Simul()
filestr = 'defstr3'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
S.L.build()
S.tx.clear()
S.rx.clear()
tx=array([760,1114,1.0])
rx=array([766,1114,1.5])
#tx=array([763,1120,1.0])
#rx=array([750,1130,1.5])
S.tx.point(tx)
S.rx.point(rx)
Ctx = S.L.pt2cy(S.tx.position[:,0])
Crx = S.L.pt2cy(S.rx.position[:,0])
Si = Signatures(S.L,Ctx,Crx)
Si.run1(cutoff=5)
#Si.run2(cutoff=3,dcut=2)
r2d = Si.rays(tx,rx)
r3d = r2d.to3D(lsss=S.L.lsss)
r3d.locbas(S.L)
r3d.fillinter(S.L)
pg = np.sum(S.L.pt,axis=1)/np.shape(S.L.pt)[1]
pg = array([pg[0],pg[1],0]).reshape(3,1)
r3d.show3(strucname='defstr3',pg=pg)
