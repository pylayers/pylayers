from pylayers.antprop.rays import *
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import pylayers.signal.bsignal as bs
import pylayers.signal.waveform as wvf 
from pylayers.simul.simulem import *
import matplotlib.pyplot as plt 
import time
print "======================="
print " start test_rays.py (Ray Tracing numpy) "
print "======================="
S = Simul()
filestr = 'defstr'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
S.L.build()
tx=array([759,1114,1.0])
rx=array([767,1114,1.5])
S.tx.clear()
S.rx.clear()
S.tx.point(tx)
S.rx.point(rx)
Ctx = S.L.pt2cy(S.tx.position[:,0])
Crx = S.L.pt2cy(S.rx.position[:,0])
fGHz=np.arange(2,11,0.1)
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
#
# Dans un sens
#
Si1 = Signatures(S.L,Ctx,Crx)
Si1.run4(cutoff=5,algo='old')
r2d = Si1.rays(tx,rx)
r3d1 = r2d.to3D(S.L)
r3d1.locbas(S.L)
r3d1.fillinter(S.L)
C1 = r3d1.eval(fGHz)
###C1.sort()
sc1 = C1.prop2tran(a='theta',b='theta')
cir1 = sc1.applywavB(wav.sfg)
#####
###### puis dans l'autre 
######
#####print "second rayon"
r2d2 = r2d.reciprocal()
###### get new reciprocal r3d
r3d2 = r2d2.to3D(S.L)
r3d2.locbas(S.L)
r3d2.fillinter(S.L)
C2=r3d2.eval(fGHz)
#####C2.sort()
sc2=C2.prop2tran(a='theta',b='theta')
cir2 = sc2.applywavB(wav.sfg)
######
######print r3d1[2]['sig'][:,:,0]
######print r3d2[2]['sig'][:,:,1]
######
######
r3d1.check_reciprocity(r3d2)
C1.check_reciprocity(C2)
#plt.plot(cir1.x,cir1.y,'b',cir2.x,cir2.y,'r')
