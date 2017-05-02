#-*- coding:Utf-8 -*-
import numpy as np 
import scipy.stats as st 
from pylayers.util.pyutil import coldict
from itertools import permutations
import  matplotlib.pyplot as plt 
import pdb


def delay_association(s):
    """
    """


if __name__=="__main__":

    Nsnap = 100
    u = np.arange(Nsnap)
    np.random.seed(12)
    #E1=st.expon(4)
    #E2=st.expon(0.2)
    #A1=st.norm(1,0.2)

#
#
#  Construct a delay/time test signal with N path and Nsnap snapshots
#
# ta : array of amplitude 
# tt : array of delays 
#
    N = 6 
    ta = np.zeros((1,Nsnap))
    tt = np.zeros((1,Nsnap))
    for k in range(N):
        tk = (4+10*np.random.rand(1))*np.cos(2*np.pi*4*np.random.rand(1)*u/100.+2*np.pi*np.random.rand(1))+0.0*np.random.rand(Nsnap)
        ak = 10*(np.sin(2*np.pi*6*np.random.rand(1)*u/100.+2*np.pi*np.random.rand(1))**2+0.4)
        ta= np.vstack((ta,ak))
        tt= np.vstack((tt,tk))

W = 4
bpermut = np.arange(N)
s = np.arange(N)
for w in range(2,W):
    for k1 in range(N-(w-1)):
        d = s[k1:k1+w]
        for k2 in permutations(d):
            p = np.hstack((s[0:k1],k2,s[k1+w:]))
            diff = p[None,:]-bpermut
            sdiff = np.sum(np.abs(diff),axis=1)
            if sdiff.all()!=0:
                bpermut = np.vstack((bpermut,p))
#bpermut = np.arange(N)
#for k in permutations(np.arange(N)):
#     bpermut=np.vstack((bpermut,np.array(k)))

#for k in range(Nsnap):
#    xk = E1.rvs(N)
#    dk = E2.rvs(N)
#    ak = 10**(A1.rvs(N))
#    tk = np.cumsum(xk)
#    akn = 0.1*ta[:,-1][:,None]+0.9*ak[:,None]
#    tkn = 0.1*tt[:,-1][:,None]+0.9*tk[:,None]
#    ta= np.hstack((ta,akn))
#    tt= np.hstack((tt,tkn))

# 5x100x2
#
# Ground truth signal 
#
tp_gt = np.dstack((ta[1:,:],tt[1:,:]))

ta_ = np.zeros((N,1))
tt_ = np.zeros((N,1))
tu_ = np.zeros((N,1))
for k in range(Nsnap):
    uk  = np.argsort(tp_gt[:,k,1])
    tu_ = np.hstack((tu_,uk[:,None]))
    ta_ = np.hstack((ta_,tp_gt[uk,k,0][:,None]))
    tt_ = np.hstack((tt_,tp_gt[uk,k,1][:,None]))
#
# Sorted signal : this is the signal the algorithm
# start with in order to recover at best the ground truth
#
tp_s = np.dstack((ta_[:,1:],tt_[:,1:],tu_[:,1:]))
def distance(v1,v2,p,debug=False):
    """
    Calulate distance between v1 and v2[p] 
    p is a permutation 

    Parameters
    ----------
    v1 : array N,
    v2 : array N, 
    p  : a N permutation 

    """
    #
    # Apply permutation on vector v2 
    #
    v2p = v2[np.argsort(p),...]
    #dist =  np.sum(tp2[0]-tp1[0])+np.sum(tp2[1]-tp2[1])
    #dist =  np.sqrt(np.sum(np.abs(v2p[:,0]-p1[:,0])**2)+np.sum(np.abs(v2p[:,1]-v1[:,1])**2)+ha)
    dist =  np.sqrt(np.sum(np.abs(v2p[:,0]-v1[:,0])**2)*np.sum(np.abs(v2p[:,1]-v1[:,1])**2))
    if debug: 
        #print v2p[:,0]
        #print v1[:,0]
        #print "----"
        #if ((p==np.array([4,2,0,3,1])).all() | (p==np.array([0,1,2,3,4])).all()):
        print "current[]",v2p[:,1]
        print "expected", v1[:,1]

    #dist =  np.prod(np.abs(v2p[:,0]-v1[:,0])**2)*np.prod(np.abs(v2p[:,1]-v1[:,1])**2)
    #dist =  np.sqrt(np.sum(np.abs(v2p[:,1]-v1[:,1])**2))
    return(dist) 

#plt.figure()
#for k in range(N):
    #plt.scatter(np.arange(Nsnap),np.arange(tp[k,:,1]),s=np.arange(tp[k,:,0]),c=color[k])
#    plt.scatter(np.arange(Nsnap),tp_gt[k,:,1],s=tp_gt[k,:,0]*3,c=np.random.rand(1,3))
#    plt.title('Ground Truth')
#plt.figure()
#for k in range(N):
    #plt.scatter(np.arange(Nsnap),np.arange(tp[k,:,1]),s=np.arange(tp[k,:,0]),c=color[k])
#    plt.scatter(np.arange(Nsnap),tp_s[k,:,1],s=tp_s[k,:,0]*3,c=np.random.rand(1,3))
#    plt.title('Sorted paths')
#plt.figure()
#ddist = {} 
#dhamm = {} 
#p1=tp[:,0,:]
#p2=tp[:,1,:]
#for p in lpermut:
#    dhamm[tuple(p)] = np.max(np.abs(p-np.arange(N)))
#    ddist[tuple(p)]=distance(p1,p2,p)

def p_opt(v1,v2,bpermut,debug=False):
    """ find the optimum permutation 
    w.r.t the distance function 

    Parameters
    ----------
    v1 : array 1 (N, 
    v2 : array 2 (N,
    bpermut : array of permutations 

    """
    d_min = 1e10
    N = np.shape(v1)[0] 
    popt = np.arange(N)
    Npermut = bpermut.shape[0]
    if debug:
        #print bpermut
        pass
    for k in range(Npermut):
        p = bpermut[k,:]
        d = distance(v1,v2,p,debug)
        if d<d_min:
            d_min = d
            popt = p
        if debug:
            #if ((p==np.array([4,2,0,3,1])).all() | (p==np.array([0,1,2,3,4])).all()):
            print p,d
    return popt,d_min

#imin = np.where(ddist.values()==np.min(ddist.values()))
#print("imin : ",imin)
#plt.subplot(211)
#plt.plot(np.arange(len(lpermut)),ddist.values())
#for k in imin[0]:
#    plt.scatter(k,ddist.values()[k],c='r',s=20)
#plt.subplot(212)
#plt.plot(np.arange(len(lpermut)),dhamm.values(),'k')

print "start reconstruction"
tp_r = np.zeros((N,Nsnap,5))
popt = np.arange(N).astype(int)
#pprev = np.arange(N).astype(int)
#pprev = tp_s[:,0,2].astype(int)
#
# The inverse permutation is obtained with function argsort 
#      p        p^-1 = argsort(p) 
#  gt -----> s  ----------------> gt
#
pprev = tp_s[:,0,2].astype(int)
for k in range(Nsnap):
    if k==87:
        debug=True
    else:
        debug=False
    if k>4: 
        tp_expected = tp_r[:,k-1,:]+(tp_r[:,k-1,:]-tp_r[:,k-5,:])/4.
        #pprev = tp_r[:,k-1,2].astype(int)
        if k==87:
            print "oo :",tp_s[pprev,k,1]
        #popt,dmin = p_opt(tp_expected,tp_s[pprev,k,:],bpermut,debug)
        popt,dmin = p_opt(tp_expected[pprev,:],tp_s[:,k,:],bpermut,debug)
        #pprev = pprev[popt]
    elif k>3:
        tp_expected = tp_r[:,k-1,:]+(tp_r[:,k-1,:]-tp_r[:,k-4,:])/3.
        #pprev = tp_r[:,k-1,2].astype(int)
        #popt,dmin = p_opt(tp_expected,tp_s[pprev,k,:],bpermut,debug)
        popt,dmin = p_opt(tp_expected[pprev,:],tp_s[:,k,:],bpermut,debug)
        #pprev= pprev[popt]
    elif k>2: 
        tp_expected = tp_r[:,k-1,:]+(tp_r[:,k-1,:]-tp_r[:,k-3,:])/2.
        #pprev = tp_r[:,k-1,2].astype(int)
        #popt,dmin = p_opt(tp_expected,tp_s[pprev,k,:],bpermut,debug)
        popt,dmin = p_opt(tp_expected[pprev,:],tp_s[:,k,:],bpermut,debug)
        #pprev = pprev[popt] 
    elif k>1:
        tp_expected = tp_r[:,k-1,:]+(tp_r[:,k-1,:]-tp_r[:,k-2,:])
        #pprev = tp_r[:,k-1,2].astype(int)
        #popt,dmin = p_opt(tp_expected,tp_s[pprev,k,:],bpermut,debug)
        popt,dmin = p_opt(tp_expected[pprev,:],tp_s[:,k,:],bpermut,debug)
        #print k,popt
        #pprev = pprev[popt] 
    elif k>0:
        popt,dmin = p_opt(tp_s[:,k-1,:],tp_s[:,k,:],bpermut,debug)
        #print k,popt
        #pprev = pprev[popt]
        #tp_r[:,k-1,2].astype(int)
    else:
        pass
        #u_expected_sorted = np.argsort(tp_expected[:,1])
        #tp_expected_sorted = tp_expected[u_expected_sorted,:]
        #popt,dmin = p_opt(tp_expected_sorted,tp_s[:,k,:],bpermut,N)
          #print(popt)
    tp_r[:,k,0] = tp_s[np.argsort(pprev[popt]),k,0]
    tp_r[:,k,1] = tp_s[np.argsort(pprev[popt]),k,1]
    #tp_r[:,k,2] = pprev[popt]
    tp_r[:,k,4] = pprev
    tp_r[:,k,2] = popt
    if k>1:
        tp_r[:,k,3] = tp_expected[:,1]
    print "~~~~~~~~~~~"
    print "k : ",k
    print "gt delay :",tp_gt[:,k,1]
    print "recovered delay :",tp_r[:,k,1]
    try:
        print "tp_expected :",tp_expected[:,1] 
    except:
        pass
    print "gt pk :",tp_s[:,k,2]
    print "pprev :",pprev
    print "popt :",popt
    try:
        print "tp_expected sorted pprev :",tp_expected[pprev,1]
    except:
        print "problem",pprev
    print "s[k] :",tp_s[:,k,1]
    pbest_ = np.argsort(pprev)[tp_s[:,k,2].astype(int)] 
    print "pbest :",pbest_
    #print "popt  :",tp_r[:,k,2]
    #pprev_ = tp_r[:,k,4].astype(int)
    #popt_ = tp_r[:,k,2].astype(int)
    #print "recovered delay from sorted (1):",tp_s[pprev_[popt_],k,1]
    #print "recovered delay from sorted (2):",tp_s[tp_s[:,k,2].astype(int),k,1]
    #print "recovered delay from sorted (3):",tp_s[np.argsort(tp_s[:,k,2].astype(int)),k,1]
    #print "pbest pk -1",np.argsort(tp_s[:,k,2].astype(int))
    #print "Good permutation ",tp_s[:,k,2].astype(int)
    #print "pprev : ",pprev_
    #print "popt : ",popt_
    #print "pprev[popt]",pprev_[popt_]
    #print "pprev",pprev
    #print "pprev[pbest]",pprev_[pbest_]
    #print "delay (pprev[pbest])",tp_s[np.argsort(tp_s[:,k,2]),k,1]
    #print "delay (pprev[pbest-1])",tp_s[pprev_[np.argsort(pbest_)],k,1]
    #print "delay (pprev[popt])",tp_s[pprev_[popt_],k,1]
    #print "pprev[popt]^-1",np.argsort(pprev_[popt_])
    pprev = pprev[popt]
    print "pprev (output) :",pprev
    if (pprev!=tp_s[:,k,2]).any():
        print "Erreur en ",k
    print "=================="
#
#plt.figure()
#for k in range(N):
#    #plt.scatter(np.arange(Nsnap),np.arange(tp[k,:,1]),s=np.arange(tp[k,:,0]),c=color[k])
#    plt.scatter(np.arange(Nsnap),tp_r[k,:,1],s=tp_r[k,:,0]*3,c=np.random.rand(1,3),label=str(k))
#    plt.title('Reconstructed paths')
#plt.legend(loc='best')
#for k in range(N): 
#    plt.subplot(211)
#    plt.plot(tp_r[k,:,0],label=str(k))
#    plt.plot(tp_gt[k,:,0],'o')
#    plt.subplot(212)
#    plt.plot(tp_r[k,:,1],label=str(k))
#    plt.plot(tp_gt[k,:,1],'o')
#E = np.sqrt((tp_gt[:,:,1]-tp_r[:,:,1])**2)
#print np.sum(E)
#plt.legend()
plt.figure()
imin=0
#imax = 10
imax = Nsnap
for k in range(N):
    plt.subplot(N,1,k+1)
    plt.plot(tp_gt[k,imin:imax,1],'o',label='gt : '+str(k))
    plt.plot(tp_r[k,imin:imax,1],'o',label='r : '+str(k))
    plt.plot(tp_r[k,imin:imax,3],'or',label='e : '+str(k))
#    
plt.legend()
plt.figure()
y = np.argsort(tp_r[:,:,2],axis=0)
for k in range(N):
    plt.subplot(N,1,k+1)
    plt.plot(tp_s[k,imin:imax,2],'o',label=str(k))
    plt.plot(y[k,imin:imax],'o',label=str(k))
plt.show()
