import mayavi.mlab as mlab
from pylayers.gis.layout import *
from pylayers.simul.link import *

L = Layout('Luebbers.ini')
#L.showG('st',aw=True,labels=True,nodelist=L.ldiffout)
#f,lax= plt.subplots(2,2)
#L.showG('s',aw=True,labels=True,fig=f,ax=lax[0][0])
#lax[0][0].set_title('Gs',fontsize=18)
#L.showG('st',aw=True,labels=True,fig=f,ax=lax[0][1])
#lax[0][1].set_title('Gt',fontsize=18)
#L.showG('v',aw=True,labels=True,fig=f,ax=lax[1][0])
#lax[1][0].set_title('Gv',fontsize=18)
#L.showG('i',aw=True,labels=True,fig=f,ax=lax[1][1])
#lax[1][1].set_title('Gi',fontsize=18)
#
fGHz=np.arange(0.5,1,0.01)
DL = DLink(L=L,fGHz=fGHz)
DL.Aa = Antenna('Omni',fGHz=fGHz)
DL.Ab = Antenna('Omni',fGHz=fGHz)
DL.a = np.array([37.5,6.2,1.5])
DL.b = np.array([12.5,30,1.5])
DL.eval(force=True,cutoff=4,verbose=False,ra_ceil_H=0)
#DL.R.show(L=L)
