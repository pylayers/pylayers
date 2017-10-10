from pylayers.simul.link import *
import pdb
#fGHz=np.arange(59,61,0.01)
fGHz=np.array([32])
L = Layout('TC2_METIS.lay',bbuild=True)
a = np.array([200,300,15])
b = np.array([130,80,1.2])
ca = L.pt2cy(a)
cb = L.pt2cy(b)
DL=DLink(L=L,fGHz=fGHz)
DL.a = a
DL.b = b
DL.Aa = Antenna(typ='Omni',param={'pol':'t','GmaxdB':2})
DL.Ab = Antenna(typ='hplanesectoralhorn',fGHz=32)
DL.eval(verbose=0,force=True,ra_ceil_H=0,applywav=False,cutoff=4)
DL.show(bsig=False,rays=True)
plt.show()
