from pylayers.gis.layout import *
from pylayers.simul.link import *

L = Layout('Luebbers.ini',force=True,check=False)
L.build()
L.dumpw()
plt.ion()
#L.showG('st',aw=True,labels=True,nodelist=L.ldiffout)
f,lax= plt.subplots(2,2)
L.showG('s',aw=True,labels=True,fig=f,ax=lax[0][0])
lax[0][0].set_title('Gs',fontsize=18)
L.showG('st',aw=True,labels=True,fig=f,ax=lax[0][1])
lax[0][1].set_title('Gt',fontsize=18)
L.showG('v',aw=True,labels=True,fig=f,ax=lax[1][0])
lax[1][0].set_title('Gv',fontsize=18)
L.showG('i',aw=True,labels=True,fig=f,ax=lax[1][1])
lax[1][1].set_title('Gi',fontsize=18)

#DL = DLink(L=L)
#DL.a = np.array([37.5,6.2,1.5])
#DL.b = np.array([12.5,30,1.5])
#DL.eval()
#DL.R.show(L=L)
