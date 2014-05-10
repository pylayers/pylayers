from pylayers.gis.layout  import *
from pylayers.simul.link  import *
from pylayers.antprop.signature import *
L = Layout('scattering.ini')
L.build()
#Si = Signatures(L,0,1)
#Si.run5(cutoff=4)
#Si.show(L,ctx=0,crx=1)
Lk = DLink(L=L,cutoff=4,force=True)
pta = Lk.a
ptb = Lk.b
pta[0] = pta[0]-4
ptb[1] = ptb[1]+60
ptb[0] = ptb[0]-12
Lk.a = pta
Lk.b = ptb
Lk.eval(si_algo='old',force=True)
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

