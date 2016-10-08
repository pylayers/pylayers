from pylayers.antprop.rays import *
from pylayers.gis.layout import *
from pylayers.antprop.signature import Signatures
import pylayers.signal.bsignal as bs
import pylayers.signal.waveform as wvf 
from pylayers.simul.link import *
import matplotlib.pyplot as plt
import time
print "======================="
print " start test_reciprocity.py "
print "======================="
S1 = DLink(L='defstr.ini')
S2 = DLink(L='defstr.ini')
S1.a=array([759,1114,1.0])
S1.b=array([767,1114,1.5])
S2.a=array([767,1114,1.5])
S2.b=array([759,1114,1.0])
fGHz =np.arange(2,11,0.1)
wav = wvf.Waveform(fcGHz=5,bandGHz=3)
#
# Dans un sens
#
S1.eval(force=1)
S2.eval(force=1)
######
####### puis dans l'autre
#######
######print "second rayon"
#print '##############'
#print '# reciprocal #'
#print '##############'
#
#r2d2 = r2d.reciprocal()
####### get new reciprocal r3d
#r3d2 = r2d2.to3D(S.L)
#r3d2.locbas(S.L)
#
#r3d2.fillinter(S.L)
#C2=r3d2.eval(fGHz)
######C2.sort()
#sc2=C2.prop2tran()
#sc2.sort()
#chw2 = sc2.apply(wav.sfg)
#rir2 = chw2.rir(Nz=500,ffts=1)
#plt.imshow(rir2,interpolation='nearest',cmap=plt.cm.jet)
#plt.axis('auto')
#plt.figure()
#plt.imshow(np.log10(abs(rir2)),interpolation='nearest',cmap=plt.cm.jet)
#plt.axis('auto')
## sc2=C2.prop2tran()
## chw = sc2.apply(wav.sfg)
## cir = chw.rir(Nz=500,ffts=1)
## plt.imshow(cir.y[:,0,0,:],interpolation='nearest')
## plt.axis('auto')
## cir2 = sc2.applywavB(wav.sfg)
#######
#######print r3d1[2]['sig'][:,:,0]
#######print r3d2[2]['sig'][:,:,1]
#######
#######
#r3d1.check_reciprocity(r3d2)
#C1.check_reciprocity(C2)
## plt.figure()
## plt.plot(cir1.x,cir1.y[0,0,:],'b',cir2.x,cir2.y[0,0,:],'r')
## plt.axis('auto')
## plt.figure()
## plt.plot(cir1.x,cir1.y[0,0,:],'b',cir2.x,cir2.y[0,0,:],'r')
## plt.axis('auto')
