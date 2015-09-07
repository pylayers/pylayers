# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pylayers.gis.layout import *
from pylayers.util.geomutil import *
from pylayers.antprop.signature import *
from pylayers.antprop.channel import *
import pylayers.signal.waveform as wvf
import networkx as nx
import numpy as np
import time
import logging
from IPython.display import Image,HTML,Latex


# <codecell>


L = Layout('WHERE1.ini')
#L = Layout('defstr2.ini')
#try:
#    L.dumpr()
#except:
L.build()
L.dumpw()


# <markdowncell>

# A signature is calculated between two cycles. This presents an interest for evaluating coverage maps. All links which have their termination in the same cycle can be obtained with the same signature. 

# <codecell>

nc1 = 7
nc2 = 2


# <markdowncell>

# The topological graph $\mathcal{G}_t$ contains all the cycles of the layout. A member data of the node of $\mathcal{G}_t$ is the shapely polygon associated with the cycle. 

# <codecell>

poly1 = L.Gt.node[nc1]['polyg']
cp1 = poly1.centroid.xy

poly2 = L.Gt.node[nc2]['polyg']
cp2 = poly2.centroid.xy

# <markdowncell>

# Here the centrod of each cycle is chosen in order to become respectively the transmitter and the receiver of the link. 

# <codecell>

ptx = np.array([cp1[0][0],cp1[1][0],1.5])
prx = np.array([cp2[0][0],cp2[1][0],1.5])
print 'Tx : ',ptx
print 'Rx : ',prx
v   = prx-ptx
mv  = np.sqrt(np.sum(v*v,axis=0))
vn  = v/mv
print vn


# <codecell>

d = np.sqrt(np.dot((ptx-prx),(ptx-prx)))
tau = d/0.3
print "Distance (m): ,",d
print "Delay (ns): ,",tau


S=Signatures(L,nc1,nc2)

t1 = time.time()
S.run3(cutoff=2,dcut=3)
t2 = time.time()
print "elapsed time :",t2-t1


t1 = time.time()
r=S.rays(ptx,prx)
r3=r.to3D()
r3.locbas(L)
r3.fillinter(L)
t2 = time.time()
print "elapsed time :",t2-t1

r3.show3(strucname='WHERE1')
