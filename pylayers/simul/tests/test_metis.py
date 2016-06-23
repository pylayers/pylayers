from pylayers.simul.link import *
import pdb
fGHz=np.arange(59,61,0.01)
L = Layout('TC2_METIS_new.ini')
a = np.array([200,300,15])
b = np.array([20,80,1.2])
ca = L.pt2cy(a)
cb = L.pt2cy(b)
DL=DLink(L=L,fGHz=fGHz)
DL.a = a
DL.b = b
#print ca
#print cb
DL.eval(verbose=0,force=True,ra_ceil_H=0,applywav=False,cutoff=6)
#Si = Signatures(L,ca,cb,cutoff=6)
#Si.run()
