from pylayers.gis.layout  import *
from pylayers.simul.link  import *
L = Layout('scattering.ini')
L.build()
Lk = DLink(L=L,cutoff=3,force=True)
Lk.a[0]=Lk.a[0]-4
Lk.b[1]=Lk.b[1]+60
Lk.eval(force=True)

# list of points
#lpnt = np.array(filter(lambda x : x <0,L.Gs.node))
# degree of points
#ldeg = np.array(map(lambda x : nx.degree(L.Gs,x),lpnt))
# index of points of degree > 2
#u = np.where(ldeg>2)[0]
# list of point of degree > 2
#ldegsup2 = lpnt[u]
# list of point of degree 2 including airwalls
#ldeg2 = L.degree[2]

#ldeg2f = filter(lambda x : x in ldeg2 ,ldegsup2)

#aseg = map(lambda x : filter(lambda y : y not in L.name['AIR'],nx.neighbors(L.Gs,x)),ldeg2f)

