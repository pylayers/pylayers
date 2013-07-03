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
S.L.Gs.node[1]['ss_name']=['AIR','AIR','WOOD']
S.L.g2npy()
S.L.build()
S.tx.clear()
S.rx.clear()
#tx=array([1,0.3,1.0])
#rx=array([7,0,1.0])
tx=array([759,1114,1.0])
rx=array([767,1114,1.5])
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
r3d.show3(strucname='defstr3',pg=pg)
fGHz=np.arange(2,11,0.5)
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
r3d.fillinter(S.L,append=True)
Cair=r3d.eval(fGHz)
scair=Cair.prop2tran(a='theta',b='theta')
cirair = scair.applywavB(wav.sfg)

S.L.Gs.node[1]['ss_name']=['WOOD','AIR','METAL']
S.L.g2npy()
# graph to numpy 
Cwood=r3d.eval(fGHz)
scwood=Cwood.prop2tran(a='theta',b='theta')
cirwood = scwood.applywavB(wav.sfg)

S.L.Gs.node[1]['ss_name']=['METAL','AIR','WOOD']
# graph to numpy 
S.L.g2npy()
r3d.fillinter(S.L,append=True)
Cmetal=r3d.eval(fGHz)
scmetal=Cmetal.prop2tran(a='theta',b='theta')
cirmetal = scmetal.applywavB(wav.sfg)


fig=plt.figure()
ax1 = fig.add_subplot(311)
cirair.plot()
plt.title("['AIR','AIR','WOOD']")
ax2 = fig.add_subplot(312,sharey=ax1)
cirwood.plot()
plt.title("['WOOD','AIR','WOOD']")
ax3 = fig.add_subplot(313,sharey=ax1)
cirmetal.plot()
plt.title("['METAL','AIR','WOOD']")
plt.show()
