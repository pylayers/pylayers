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
tc = tcir[1]
for k,cir in enumerate(tcir[2:]):
    tc = tc.align(cir)
tc.imshow(dB=False)
plt.xlabel('ns')
plt.ylabel('trajectory index')
plt.title('CIR along a trajectory - crossing a wall')
#    x = cir.x
#    y = cir.y
#
#    # expand ty to length (tx)
#    Ny = len(y)
#    shty = np.shape(ty)
#    #pdb.set_trace()
#    if len(shty)==2:
#        Nty = shty[1]
#        Nlignes = shty[0]
#    else:
#        Nty = shty[0]
#
#    ndiff = Nty - Ny 
#    print 'Nlignes :',Nlignes
#    print 'longueur cir:',Ny
#    print 'longueur tableau :',Nty
#    print 'ecart :',ndiff
#    # x shorter than tx
#    if ndiff ==0:
#        ty = np.vstack((ty,y))
#    if ndiff > 0:
#        # indice de tx qui contient x[0]
#        u  = np.nonzero(tx==x[0])[0][0]
#        # indice de ty qui contient x[-1]
#        v  = np.nonzero(tx==x[-1])[0][0]
#        ty = np.vstack((ty,np.zeros(Ntx)))
#        print ndiff,np.shape(ty)
#        if ((u>=0) & (v<=Ntx)):
#            print u,v
#            ty[-1,u:v+1]=y
#        elif u >=0:
#            print u,v
#            ty[-1,u:len(y)+1]=y
#        elif v <Ntx:
#            print u,v
#            ty[-1,v-len(y):v+1]=y
#    # x plus long que tx
#    # il  faut ajouter des zeros
#    if ndiff < 0:
#        # indice de tx qui contient x[0]
#        u  = np.nonzero(tx==x[0])
#        # indice de ty qui contient x[-1]
#        v  = np.nonzero(tx==x[-1])
#        print u
#        print v
#        if abs(tx[0]-x[0])<1e-7:
#            zr = np.zeros(Nlignes,ndiff)
#            ty = np.hstack((ty,zr))
#
#        pdb.set_trace()
#        ty = np.hstack((zl,ty))
#
#        tx = np.union1d(tx,x)
#
