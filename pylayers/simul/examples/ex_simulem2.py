from pylayers.simul.simulem import *

# create a Simul object

S = Simul()

# loading a layout 
filestr = 'defstr'
S.layout(filestr+'.str','matDB.ini','slabDB.ini')
# setting transmitter

S.tx = RadioNode(typ='tx')
S.tx.point([1.2,1,1.4])

# setting receiver

S.rx = RadioNode(typ='rx')
S.rx.line(100,ptt=[2,1.2,1.5],pth=[8.2,1.1,1.6])

S.save()

# print launching parameters
S.palch.info()

# ang Tx : angular step from Tx
S.palch.angTx  = 1

# ISB ang Incident Shadow Boundary angle (degree) 
S.palch.ISBang = 90  

# ray elimination Threshold 
S.palch.ethreshold = 0.001

# maximum depth
S.palch.maxdeep  = 10

# typealgo = 0 (include diffraction) 1 (no diffraction)
S.palch.typalgo = 1
title = str(S.palch.angTx) + '-' +\
        str(S.palch.ISBang) + '-' +\
        str(S.palch.ethreshold) + '-' + \
        str(S.palch.maxdeep) + '-' + \
        str(S.palch.typalgo)

S.palch.save()
S.pafreq.fghzmin=2
S.pafreq.fghzmax=11
S.pafreq.nf=181
S.pafreq.save()

S.run(1,verbose=True)
