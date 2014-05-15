from pylayers.gis.layout  import *
from pylayers.simul.link  import *
from pylayers.antprop.signature import *
L = Layout('scattering.ini',force=True)
L.build()
L.dumpw()
#Si = Signatures(L,0,1)
#Si.run5(cutoff=4)
#Si.show(L,ctx=0,crx=1)
Lk = DLink(L=L,cutoff=3,force=True)
pta = np.array([-20,-20,1.2])
ptb = np.array([-30,0,1.2])
Lk.a = pta
Lk.b = ptb
Lk.eval(si_algo='old',force=True)
#fig = plt.figure()
#for i in Lk.R.keys():
#    ax = fig.add_subplot(3,3,i)
#    fig,ax = Lk.R.show(L=L,i=i,fig=fig,ax=ax)


