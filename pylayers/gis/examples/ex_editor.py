from pylayers.gis.layout import *
import numpy as np
import networkx as nx 
from pylayers.util.geomutil import *

#L = Layout('DLR.ini','matDB.ini','slabDB.ini')
L = Layout()
L.display['fileoverlay']="TA-Office.png"
#L.saveini('TA-Office.ini')
# L.build()
#C = nx.algorithms.cycles.cycle_basis(L.Gs)
#for n in C[4]:
#    try:
#        p = np.vstack((p,array(L.Gs.pos[n]).T))
#    except:
#        p = array(L.Gs.pos[n]).T
#P = Polygon(p) 
#P.plot()
#plt.show()
##L.editor()
#
