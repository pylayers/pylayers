from pylayers.gis.cycles import *
from pylayers.gis.layout import *
import networkx as nx

#L = Layout('TA-Office2.ini')
#L = Layout('DLR.ini')
#L = Layout('defstr.ini')
L = Layout('exemple.str')
L.build()
L.dumpw()
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()


