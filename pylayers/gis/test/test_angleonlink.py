from pylayers.gis.layout import *
L = Layout('TA-Office.ini')
p1 = np.array([3,3])
p2 = 10*np.random.rand(2,1000)
sgl = L.angleonlink(p1,p2)
L.sla[sgl['s']]
