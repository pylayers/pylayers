from pylayers.gis.layout  import *
from pylayers.simul.link  import *
L = Layout('scattering.ini')
L.build()
Lk = DLink(L=L,cutoff=3)
Lk.a[0]=Lk.a[0]-4
Lk.b[1]=Lk.b[1]+60
Lk.eval()
