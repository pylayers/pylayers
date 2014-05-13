from pylayers.gis.layout  import *
from pylayers.simul.link  import *
from pylayers.antprop.signature import *
L = Layout('scattering.ini',force=True)
L.build()
L.dumpw()
#Si = Signatures(L,0,1)
#Si.run5(cutoff=4)
#Si.show(L,ctx=0,crx=1)
Lk = DLink(L=L,cutoff=4,force=True)
pta = Lk.a
ptb = Lk.b
pta[0] = pta[0]-4
ptb[1] = ptb[1]+20
ptb[0] = ptb[0]-18
Lk.a = pta
Lk.b = ptb
Lk.eval(si_algo='old',force=True)
#fig = plt.figure()
#for i in Lk.R.keys():
#    ax = fig.add_subplot(3,3,i)
#    fig,ax = Lk.R.show(L=L,i=i,fig=fig,ax=ax)


