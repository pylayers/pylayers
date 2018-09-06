from pylayers.simul.link import *
import pdb

DL0 = DLink(L='testair0.lay')
DL1 = DLink(L='testair1.lay')
DL0.a = np.array([1,3,1])
DL0.b = np.array([8,1,2.5])
DL1.a = np.array([1,3,1])
DL1.b = np.array([8,1,2.5])
DL0.eval(force=1,cutoff=1,threshold=0.1)
DL1.eval(force=1,cutoff=1,threshold=0.1)
DL0.plt_cir()
DL1.plt_cir()
#pdb.set_trace()
#B = np.array([[0,1,0],[0,0,1],[1,0,0]])
#
## Sans mur d'air 
## Bo0 3 x 3 
#print("without air wall")
Bo0_0 = DL0.R[1]['Bo0'][:,:,0]
Bi_0 = DL0.R[1]['Bi'][:,:,:,0]
Bo_0 = DL0.R[1]['Bo'][:,:,:,0]
BiN_0 = DL0.R[1]['BiN'][:,:,0]
#print "Bo0"
#print Bo0
#print DL0.R[1]['B'][:,:,0,0]
##"print np.dot(Bo0[:,1:].T,Bi[:,1:,0])
##print np.dot(Bi[:,1:,0].T,Bo0[:,1:])
## get the indices of interactions
#linter = DL0.R[1]['rays'][0]
#for k in range(Bi.shape[2]):
#    print "Bi"+str(k)
#    print "  ",Bi[:,:,k]
#    print "Interaction"+str(k)
#
#    print "  ",DL0.R.I.I[0,linter[k],:,:]
#    print "Bo"+str(k)
#    print "  ",Bo[:,:,k] 
#print DL0.R[1]['B'][:,:,1,0]
##print np.dot(Bo[:,1:,0].T,BiN[:,1:])
##print np.dot(BiN[:,1:].T,Bo[:,1:,0])
#print "BiN"
#print BiN
#print ("\n")
#print("with air wall")
## Avec mur d'air 
Bo0_1 = DL1.R[1]['Bo0'][:,:,0]
Bi_1 = DL1.R[1]['Bi'][:,:,:,0]
Bo_1 = DL1.R[1]['Bo'][:,:,:,0]
BiN_1 = DL1.R[1]['BiN'][:,:,0]
#print "Bo0"
#print Bo0
#print DL1.R[2]['B'][:,:,0,0]
#print np.dot(Bo0[:,1:].T,Bi[:,1:,0])
#print np.dot(Bi[:,1:,0].T,Bo0[:,1:])
## get the indices of interactions
#linter = DL1.R[2]['rays']
#for k in range(Bi.shape[2]):
#    if k>0:
#        print DL1.R[2]['B'][:,:,k,0]
#    print "Bi"+str(k)
#    print "  ",Bi[:,:,k]
#    print "Interaction"+str(k)
#
#    print "  ",DL1.R.I.I[0,linter[0,k],:,:]
#    print "Bo"+str(k)
#    print "  ",Bo[:,:,k] 
#print k
#print DL1.R[2]['B'][:,:,2,0]
##print np.dot(Bo[:,1:,0].T,BiN[:,1:])
##print np.dot(BiN[:,1:].T,Bo[:,1:,0])
#print "BiN"
#print BiN
##-------------------
## The problem 
N0 = Bo0_0
N1 = Bi_0[:,:,0]
A = np.dot(N1.T,N0)
print( A)
M0 = Bo0_1
I  = np.eye(3)
M1 = Bi_1[:,:,0]
M2 = Bo_1[:,:,0] 
M3 = Bi_1[:,:,1] 
T1 = np.dot(M1.T,M0)
T2 = np.dot(I,T1)
T3 = np.dot(M2,T2)
T4 = np.dot(M3.T,T3)
print( T4 )
#print np.dot(M2.T,M1)
A1 = DL0.R[1]['B'][:,:,0,0]
A2 = DL0.R[1]['B'][:,:,1,0]
R0  = DL0.R.I.I[0,1,:,:]
res0 = np.dot(A2,np.dot(R0,A1))
assert np.allclose(N0,M0)
assert np.allclose(N1,M3)
B1 = DL1.R[2]['B'][:,:,0,0]
B2 = DL1.R[2]['B'][:,:,1,0]
B3 = DL1.R[2]['B'][:,:,2,0]
R1  = DL1.R.I.I[0,2,:,:]
res1 = np.dot(B3,np.dot(R1,np.dot(B2,B1)))
print( res0-res1)
#U = Bi_1[:,:,0]
#V = Bo_1[:,:,0]
#B = np.dot(Bi_1[:,:,0].T,Bo0_1)
#C = np.dot(Bo_1[:,:,0].T,B.T)
#print np.dot(U,V.T)
#print C
#print "Bi",A
#print "Bo",B
#print "BiBo",np.dot(A,B)
#
##for k in DL1.R:
##    print "groupe d'interactions ",k
##    Bo0 = DL1.R[k]['Bo0']
##    Bi = DL1.R[k]['Bi']
##    Bo = DL1.R[k]['Bo']
##    BiN = DL1.R[k]['BiN']
##    nray = Bi.shape[3]
##    for ir in range(nray):
##        for il in range(k):
##            mBi = Bi[:,:,il,ir]
##            mBo = Bo[:,:,il,ir]
##            print np.dot(mBo[:,1:].T,mBi[:,1:])
##            #dmBi = np.linalg.det(mBi)
##            #dmBo = np.linalg.det(mBo)
##            #print dmBi,dmBo
##            #print np.linalg.det(np.dot(Bi[:,:,il,ir],B.T))
##            #print np.linalg.det(np.dot(Bo[:,:,il,ir],B.T))
##
###tind = []
###for tau in DL0.H.taud:
###    if tau in DL1.H.taud:
###        u = np.where(DL1.H.taud==tau)[0]
###        v = np.where(DL0.H.taud==tau)[0]
###        try:
###            diff = np.abs(DL1.H.y[u]-DL0.H.y[v]).squeeze()
###            if diff > 1e-10: 
###                tind.append(zip(u,v))
###        except:
###            pass
