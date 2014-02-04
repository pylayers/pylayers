from pylayers.antprop.rays import *
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import pylayers.signal.bsignal as bs
import scipy.linalg as la
import pylayers.signal.waveform as wvf
from pylayers.simul.simulem import *
import matplotlib.pyplot as plt
import time
import copy
print "======================="
print " start test_rays.py (Ray Tracing numpy) "
print "======================="
S = Simul()
filestr = 'defstr'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
S.L.build()
tx=array([759,1114,1.0])
Ctx = S.L.pt2cy(S.tx.position[:,0])
Crx = S.L.pt2cy(S.rx.position[:,0])
rx=array([760,1114,1.0])
S.tx.clear()
S.rx.clear()
S.tx.point(tx)
S.rx.point(rx)
fGHz=np.arange(2,11,0.1)
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
Si1 = Signatures(S.L,Ctx,Crx)
Si1.run4(cutoff=5,algo='old')
#
# To test Friis Formula we loop on the x coordinate of the receiver.
# We start at 1 m from the transmitter and move away from it along x axis
# on 9 m
# loop on rx
#
tEtt =[]
tEpp =[]
tEtp =[]
tEpt =[]
dist = np.arange(1,9,0.2)
vx    = dist+759
for x in vx:
    rx=array([x,1114,1.0])
    print "Rx :",rx
    S.rx.clear()
    S.rx.point(rx)

    Crx2 = S.L.pt2cy(S.rx.position[:,0])
    if Crx2<>Crx:
        print "new signature"
        Crx = Crx2
        Si1 = Signatures(S.L,Ctx,Crx)
        Si1.run4(cutoff=5,algo='old')
    r2d = Si1.rays(tx,rx)
    r3d1 = r2d.to3D(S.L)
    r3d1.locbas(S.L)
    r3d1.fillinter(S.L)
    C1 = r3d1.eval(fGHz)
    Ett,Epp,Etp,Ept = C1.energy(sumray=True,Friis=True)
    tEtt.append(Ett)
    tEpp.append(Epp)
    tEtp.append(Etp)
    tEpt.append(Ept)
np.save('tEtt',inp.array(tEtt))
np.save('tEtp',inp.array(tEtp))
np.save('tEpt',inp.array(tEpt))
np.save('tEpp',inp.array(tEpp))
