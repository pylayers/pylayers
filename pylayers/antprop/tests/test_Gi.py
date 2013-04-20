#-*- coding:Utf-8 -*-
from pylayers.gis.layout import *
from pylayers.antprop.signature import *
import networkx as nx
import numpy as np
#
# supprimer les diffractions
# enrichir Gi pour gerer le second ordre
# utiliser key dans add_edge
#

L=Layout('DLR.ini')
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()
#L.build()
#L.dumpw()
g0 = L.dGi[0]
#
#
#
for e in g0.edges():
    # extract  both interactions
    i0 = eval(e[0])
    i1 = eval(e[1])
    try:
        nstr0 = i0[0]
    except:
        nstr0 = i0

    try:
        nstr1 = i1[0]
    except:
        nstr1 = i1
    p0 = np.array(L.Gs.pos[nstr0])
    p1 = np.array(L.Gs.pos[nstr1])
    v01  = p1-p0
    v01m = np.sqrt(np.dot(v01,v01))
    v01n = v01/v01m
    v10n = -v01n
    # next interaction
    for i2 in g0.nodes():
        i2 = eval(i2)
        # in the sense i0->i1->i2
        output01 = [ ]
        if type(i2)==int: # diffraction
            print "Dif :" ,i2
            nstr2 = i2
        else:
            if len(i2)>2: #
                print "T :" ,i2
            else:
                print "R :" ,i2
            nstr2 = i2[0]
            norm = L.Gs.node[i2[0]]['norm']

        p2 = np.array(L.Gs.pos[nstr2])
        v12 = p2-p1
        v12m = np.sqrt(np.dot(v12,v12))
        v12n = v12/v12m
        # Heuristic rules
        # if 0 = R 1 = R
        #    p12.p01 > 0
        # if 0 = T 1 = R
        # if 0 = T 1 = T
        # if 0 = R 1 = R
        # if 0 = R 1 = R

        print norm


        # in the sense 1->0->2
        v02 = p2-p0
        v02m = np.sqrt(np.dot(v02,v02))
        v02n = v02/v02m
        output10 = [ ]

    #print e,vn
