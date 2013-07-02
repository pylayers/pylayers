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

L=Layout('defstr.ini')
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()
#L.build()
#L.dumpw()
g0 = L.Gi
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

    output = []
    if nstr1>0:
        # segment unitary vector
        l1 = L.seguv(np.array([nstr1]))
        p0 = np.array(L.Gs.pos[nstr0])
        p1 = np.array(L.Gs.pos[nstr1])
        v01  = p1-p0
        v01m = np.sqrt(np.dot(v01,v01))
        v01n = v01/v01m
        v10n = -v01n
        # next interaction
        for i2 in nx.neighbors(g0,str(i1)):
            i2 = eval(i2)
            if type(i2)==int:
                nstr2 = i2
            else:
                nstr2 = i2[0]
            p2 = np.array(L.Gs.pos[nstr2])
            v12 = p2-p1
            v12m = np.sqrt(np.dot(v12,v12))
            v12n = v12/v12m
            d1 = np.dot(v01n,l1)
            d2 = np.dot(l1,v12n)
            if d1*d2>=0:
                output.append(i2)
            else:
                pass
                #print i0,i1,i2,v01n,v12n,l2,d1,d2

            # in the sense 1->0->2
            # v02 = p2-p0
            # v02m = np.sqrt(np.dot(v02,v02))
            # v02n = v02/v02m
            # output10 = [ ]
    L.Gi.add_edge(i0,i1,output=output)
    #print i0,i1,output
