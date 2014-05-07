from pylayers.gis.layout  import *
L = Layout('scattering.ini')
#L.build()
L.build(verbose=True)
#L.dumpr()
poly =  L.Gc.node[0]['polyg']
cycle = L.Gc.node[0]['cycle']
poly.vnodes=cycle.cycle
Gv = poly.buildGv(show=True,eded=False)
