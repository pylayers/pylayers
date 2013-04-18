# -*- coding: utf-8 -*-

from pylayers.simul.simulem import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
from pylayers.antprop.signature import *
import pylayers.util.pyutil as pyu
from pylayers.gis.layout import *
import pylayers.signal.bsignal as bs
from datetime import datetime
import time
import pickle



S = Simul()
filestr = 'DLR'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
try:
    S.L.dumpr()
except:
    S.L.build()
    S.L.dumpw()


tx = S.tx.position[:,4]
Rtx = S.L.pt2ro(tx)
print "transmitter :",tx," is in room ",Rtx

rx = array([15,3,2.5])
#rx = S.rx.position[:,28]
Rrx = S.L.pt2ro(rx)
print "mobile node :",rx," is in room ",Rrx


Si = Signatures(S.L,tx,rx)
Si.run(tx,rx,2)
Si.showi()
