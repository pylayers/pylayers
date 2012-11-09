from pylayers.util import project
from pylayers.signal.bsignal import *
import matplotlib.pylab as plt
import numpy as np
import pdb
import glob

tcir = []
for a,b,c in glob.os.walk(basename+'/output/Tx001'):
    c.sort()
    for s in c:
        cir = TUsignal()
        cir.readcir(s,'Tx001')
        tcir.append(cir)
        print s
plt.ion()
#
tx = tcir[0].x
ty = tcir[0].y
for k,cir in enumerate(tcir[1:]):
    print "k:",k
    x = cir.x
    y = cir.y

    # expand ty to length (tx)
    Ny   = len(y)
    Mty,Nty = np.shape(ty)
    pdb.set_trace()
    ndiff = len(ty)-len(y)
    print Ny,Nty,ndiff
    # x plus court que tx
    if ndiff ==0:
        ty = np.vstack((ty,y))
    if ndiff > 0:
        # indice de tx qui contient x[0]
        u  = np.nonzero(tx==x[0])[0][0]
        # indice de ty qui contient x[-1]
        v  = np.nonzero(tx==x[-1])[0][0]
        ty = np.vstack((ty,np.zeros(Ntx)))
        print ndiff,np.shape(ty)
        if ((u>=0) & (v<=Ntx)):
            print u,v
            ty[-1,u:v+1]=y
        elif u >=0:
            print u,v
            ty[-1,u:len(y)+1]=y
        elif v <Ntx:
            print u,v
            ty[-1,v-len(y):v+1]=y
    # x plus long que tx
    # il  faut ajouter des zeros
    if ndiff < 0:
        # indice de tx qui contient x[0]
        u  = np.nonzero(tx==x[0])[0][0]
        # indice de ty qui contient x[-1]
        v  = np.nonzero(tx==x[-1])[0][0]
        if abs(tx[0]-x[0])<1e-7:
            zr = zeros(Ntx,ndiff)
            ty = np.hstack((ty,zr))

        pdb.set_trace()
        ty = np.hstack((zl,ty))

        tx = np.union1d(tx,x)

#plt.pcolor(tx,np.arange(5),ty)
#plt.axis('auto')
