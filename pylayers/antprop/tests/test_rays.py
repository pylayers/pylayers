from pylayers.antprop.rays import *
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import pylayers.signal.bsignal as bs
import pylayers.signal.waveform as wvf 
from pylayers.simul.simulem import *
import matplotlib.pyplot as plt 

S = Simul()
filestr = 'defstr3'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
S.L.Gs.node[1]['ss_name']=['WOOD','AIR','METAL']
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
r3d = r2d.to3D(S.L)
r3d.locbas(S.L)
r3d.fillinter(S.L)
pg = np.sum(S.L.pt,axis=1)/np.shape(S.L.pt)[1]
pg = array([pg[0],pg[1],0]).reshape(3,1)
#r3d.show3(strucname='defstr3',pg=pg)
fGHz=np.arange(2,11,0.1)
#Cn=r3d.eval(fGHz,ib=[1])

Cwood=r3d.eval(fGHz)
scwood=Cwood.prop2tran(a='theta',b='theta')
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
cirwood = scwood.applywavB(wav.sfg)
print(S.L.Gs.node[1])
S.L.Gs.node[1]['ss_name']=['METAL','AIR','METAL']
print(S.L.Gs.node[1])
r3d.fillinter(S.L)
Cmetal=r3d.eval(fGHz)
scmetal=Cmetal.prop2tran(a='theta',b='theta')
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
cirmetal = scmetal.applywavB(wav.sfg)
plt.subplot(211)
cirwood.plot()
plt.subplot(212)
cirmetal.plot()
plt.show()
