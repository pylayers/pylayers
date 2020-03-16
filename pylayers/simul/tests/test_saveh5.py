from pylayers.antprop.rays import *
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import pylayers.signal.bsignal as bs
import pylayers.signal.waveform as wvf 
import matplotlib.pyplot as plt 
from pylayers.antprop.antenna import *
from pylayers.antprop.channel import *

import copy
import time

from numpy import array
from tvtk.api import tvtk
from mayavi.scripts import mayavi2


S = Simul()
filestr = 'TA-Office'
S.layout(filestr+'.lay','matDB.ini','slabDB.ini')
try:
    S.L.dumpr()
except:
    S.L.build()
    S.L.dumpw()

S.tx.clear()
S.rx.clear()

# set tx and rx position from cycles
ctx = 0
crx = 1
tx=S.L.cy2pt(ctx)
rx=S.L.cy2pt(crx)+0.5
S.tx.point(tx)
S.rx.point(rx)

# set antenna orientation
Ttx=np.eye(3)
Trx=np.eye(3)
S.tx.orientation = Ttx
S.rx.orientation = Trx

# set frequency range
fGHz=np.arange(2,11,0.1)


###### RUN TEST
Si1=Signatures(S.L,ctx,crx,cutoff=2)
Si1.run5()
Si1.saveh5()
print('Sig : saveh5 passed ')
Si1=Signatures(S.L,ctx,crx,cutoff=2)
Si1.loadh5()
print('Sig : loadh5 passed ')

r2d = Si1.rays(tx,rx)
r3d1 = r2d.to3D(S.L)
r3d1.saveh5()
print('Ray3d : saveh5 passed ')
r3d1=Rays(tx,rx)
r3d1.loadh5('TA-Office_70_0.h5')
print('Ray3d : loadh5 passed ')
r3d1.locbas(S.L)
r3d1.saveh5()
print('Ray3d locbas: saveh5 passed ')
r3d1=Rays(tx,rx)
r3d1.loadh5('TA-Office_70_0.h5')
print('Ray3d locbas: loadh5 passed ')

r3d1.fillinter(S.L)
r3d1.eval(fGHz)
r3d1.saveh5()
print('Ray3d eval: saveh5 passed ')
r3d1=Rays(tx,rx)
C=r3d1.loadh5('TA-Office_70_0.h5')
print('Ray3d eval: loadh5 passed ')
C.saveh5(S.L.filename,0,tx,rx)
print('Ctilde: saveh5 passed ')
C=Ctilde()
C.loadh5(S.L.filename,0)
print('Ctilde: loadh5 passed ')
sca=C.prop2tran()
sca.saveh5(S.L.filename,0,tx,rx,Ttx,Trx)
print('Tchannel: saveh5 passed ')
sca = Tchannel()
sca.loadh5(S.L.filename,0)
print('Tchannel: loadh5 passed ')
