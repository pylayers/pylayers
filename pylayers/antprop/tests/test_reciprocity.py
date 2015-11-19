from pylayers.antprop.rays import *
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import pylayers.signal.bsignal as bs
import pylayers.signal.waveform as wvf 
from pylayers.simul.simulem import *
import matplotlib.pyplot as plt 
import time
print "======================="
print " start test_reciprocity.py "
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
# get cycle from coordinates
Ctx = S.L.pt2cy(S.tx.position[:,0])
Crx = S.L.pt2cy(S.rx.position[:,0])
fGHz=np.arange(2,11,0.1)
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
#
# Dans un sens
#
tic = time.time()
Si1 = Signatures(S.L,Ctx,Crx)
#Si1.run4(cutoff=5,algo='old')
Si1.run5()
toc = time.time()
print "signature ",toc-tic
tic = time.time()
r2d = Si1.rays(tx,rx)
print "2D rays ",tic-toc
r3d1 = r2d.to3D(S.L)
toc = time.time()
print "3D rays ",toc-tic
r3d1.locbas(S.L)
tic = time.time()
print "3D rays locbas",tic-toc
r3d1.fillinter(S.L)
toc = time.time()
print "3D fill interaction ",toc-tic
C1 = r3d1.eval(fGHz)
tic = time.time()
print "eval field ",tic-toc
###C1.sort()
<<<<<<< HEAD
sc1 = C1.prop2tran()
chw1 = sc1.apply(wav.sfg)
rir1 = chw1.rir(Nz=500,ffts=1)
=======
# sc1 = C1.prop2tran()
# cir1 = sc1.applywavB(wav.sfg)
>>>>>>> 0ecb3dc98d4ba79bcea9a4a4d9753a9ec1fec2a7
#####
###### puis dans l'autre
######
#####print "second rayon"
print '##############'
print '# reciprocal #'
print '##############'

r2d2 = r2d.reciprocal()
###### get new reciprocal r3d
r3d2 = r2d2.to3D(S.L)
r3d2.locbas(S.L)

r3d2.fillinter(S.L)
C2=r3d2.eval(fGHz)
#####C2.sort()
sc2=C2.prop2tran()
sc2.sort()
chw2 = sc2.apply(wav.sfg)
rir2 = chw2.rir(Nz=500,ffts=1)
plt.imshow(rir2,interpolation='nearest',cmap=plt.cm.jet)
plt.axis('auto')
plt.figure()
plt.imshow(np.log10(abs(rir2)),interpolation='nearest',cmap=plt.cm.jet)
plt.axis('auto')
# sc2=C2.prop2tran()
# chw = sc2.apply(wav.sfg)
# cir = chw.rir(Nz=500,ffts=1)
# plt.imshow(cir.y[:,0,0,:],interpolation='nearest')
# plt.axis('auto')
# cir2 = sc2.applywavB(wav.sfg)
######
######print r3d1[2]['sig'][:,:,0]
######print r3d2[2]['sig'][:,:,1]
######
######
r3d1.check_reciprocity(r3d2)
C1.check_reciprocity(C2)
# plt.figure()
# plt.plot(cir1.x,cir1.y[0,0,:],'b',cir2.x,cir2.y[0,0,:],'r')
# plt.axis('auto')
# plt.figure()
# plt.plot(cir1.x,cir1.y[0,0,:],'b',cir2.x,cir2.y[0,0,:],'r')
# plt.axis('auto')
