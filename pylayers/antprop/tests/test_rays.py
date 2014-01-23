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
filestr = 'defstr3'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
S.L.Gs.node[1]['ss_name']=['WOOD','AIR','METAL']
S.L.build()
tx=array([759,1114,1.0])
rx=array([767,1114,1.5])
S.tx.clear()
S.rx.clear()
S.tx.point(tx)
S.rx.point(rx)
Ctx = S.L.pt2cy(S.tx.position[:,0])
Crx = S.L.pt2cy(S.rx.position[:,0])
Si = Signatures(S.L,Ctx,Crx)
print "Signature : run4" 
tic = time.time()
Si.run4(cutoff=5,algo='old')
toc = time.time()
print "Elapsed Time : ",toc-tic
#Si.run2(cutoff=3,dcut=2)
print "Signature : rays" 
r2d = Si.rays(tx,rx)
print "Rays : to3D" 
r3d = r2d.to3D(S.L)
print "Rays : locbas" 
r3d.locbas(S.L)
print "Rays : fillinter" 
r3d.fillinter(S.L)

fGHz=np.arange(2,11,0.1)
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
print "Rays : eval" 
Cwood=r3d.eval(fGHz)
scwood=Cwood.prop2tran(a='theta',b='theta')
cirwood = scwood.applywavB(wav.sfg)

S.L.Gs.node[1]['ss_name']=['METAL','AIR','WOOD']
# graph to numpy 
print "Layout : g2npy" 
S.L.g2npy()
r3d.fillinter(S.L,append=True)
Cmetal=r3d.eval(fGHz)
scmetal=Cmetal.prop2tran(a='theta',b='theta')
cirmetal = scmetal.applywavB(wav.sfg)



S.L.Gs.node[1]['ss_name']=['AIR','AIR','WOOD']
# graph to numpy 
S.L.g2npy()
r3d.fillinter(S.L,append=True)
Cair=r3d.eval(fGHz)
scair=Cair.prop2tran(a='theta',b='theta')
cirair = scair.applywavB(wav.sfg)

print "======================="
print " stop test_rays.py (Ray Tracing numpy) "
print "======================="
