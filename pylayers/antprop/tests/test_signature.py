#!/usr/bin/python
#-*- coding:Utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from pylayers.gis.layout import *
from pylayers.antprop.signature import *

# load the layout graphs

L = Layout('TA-Office.str')
L.boundary()
print L.ax
try:
    print " Début Construction"
    L.dumpr()
    print " Fin Construction"
except:
    print " Début lecture"
    L.build()
    L.dumpw()
    print " Fin lecture"
tx = np.array([8., 8., 1.])
rx = np.array([30., 11., 2.])

#L = Layout('TA-Office.str')
#L.build()
#tx = np.array([20, 8, 1])
#rx = np.array([35, 6, 2])


S = Signatures(L, tx, rx)

print "Calcul signatures"
s1 = S.get_sigslist(tx, rx)
print "Fin calcul signatures"

print "signatures --> rayons "
r2d = S.sigs2rays(s1)
print "fin signatures --> rayons "

r22 = r2d['2']
pt2 = r22['pt']
sig2 = r22['sig']
pt2 = np.swapaxes(pt2,0,2)
pt2 = np.swapaxes(pt2,1,2)
tx2 = np.kron(np.ones(2),tx).reshape(2,3,1)
rx2 = np.kron(np.ones(2),rx).reshape(2,3,1)
tx2[:,2,:]=0
rx2[:,2,:]=0
#pt  = np.concatenate((tx2,pt2,rx2),axis=2)
#vsi = pt[:, :, 1:] - pt[:,:,:-1]
#si  = np.sqrt(np.sum(vsi*vsi, axis=1))
#alpha = np.cumsum(si,axis=1)
#c  = alpha[:,-1].reshape(2,1)
#alpha = alpha/c
#pt[:,2,1:]= tx[2]+alpha*(rx[2]-tx[2])


print "rayons 2D --> rayons3D "
rays3d = S.ray2D3D(r2d)
print "fin rayons 2D --> rayons3D "

S.show3(rays=rays3d)



