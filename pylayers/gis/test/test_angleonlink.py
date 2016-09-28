from pylayers.gis.layout import *
L = Layout('TA-Office.ini')
# single point 
p1 = np.array([3,3])
# list of points 
p2 = 10*np.random.rand(2,3)
sgl = L.angleonlink(p1,p2)
L.sla[sgl['s']]
