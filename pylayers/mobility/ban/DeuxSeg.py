# -*- coding:Utf-8 -*-
import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import scipy.linalg as la

def dist_nonvectorized(A,B,C,D,alpha,beta):
    """

    Parameters
    ----------

    A
    B
    C
    D
    alpha
    beta

    """
    AC=C-A
    CD=D-C
    BA=A-B

    u0 = np.dot(AC,AC)
    u4 = np.dot(BA,BA)
    u5 = np.dot(CD,CD)
    u1 = np.dot(BA,AC)
    u2 = np.dot(CD,AC)
    u3 = np.dot(CD,BA)

    f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
    M  = A - alpha*BA
    N  = C + beta*CD
    g  = np.dot(M-N,M-N)
    return(f,g)

def dmin3d_nonvectorized(A,B,C,D):
    """
     dmin3d evaluate the minimal distance between 2 set of segments 

     this should be vectorized 

      A : (3xN) initial point segment 1
      B   (3xN) end point segment 1
      C   (3xN) starting point segment 2
      D   (3xN) end point segment 2  
    """

    AC=C-A
    CD=D-C
    BA=A-B

    u0 = np.dot(AC,AC)
    u4 = np.dot(BA,BA)
    u5 = np.dot(CD,CD)
    u1 = np.dot(BA,AC)
    u2 = np.dot(CD,AC)
    u3 = np.dot(CD,BA) 


    den   = u4*u5-u3*u3
    alpha = (u2*u3-u1*u5)/(1.*den)
    beta  = (u1*u3-u2*u4)/(1.*den)
    #~ print ' dmin**2 ', u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
    dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5) 
    #~ print 'dmin', dmin
    return(alpha,beta,dmin)


# def dist_old(A,B,C,D,alpha,beta):
#     """

#     Parameters
#     ----------

#     A
#     B
#     C
#     D
#     alpha
#     beta

#     """
#     if len(A.shape) ==1 :
#         A=A.reshape(3,1) 
#     if len(B.shape) ==1 :
#         B=B.reshape(3,1) 
#     if len(C.shape) ==1 :
#         C=C.reshape(3,1) 
#     if len(D.shape) ==1 :
#         D=D.reshape(3,1) 


#     AC=C-A
#     CD=D-C
#     BA=A-B

#     u0 = np.einsum('ij,ij->j',AC,AC)#np.dot(AC,AC)
#     u4 = np.einsum('ij,ij->j',BA,BA)#np.dot(BA,BA)
#     u5 = np.einsum('ij,ij->j',CD,CD)#np.dot(CD,CD)
#     u1 = np.einsum('ij,ij->j',BA,AC)#np.dot(BA,AC)
#     u2 = np.einsum('ij,ij->j',CD,AC)#np.dot(CD,AC)
#     u3 = np.einsum('ij,ij->j',CD,BA)#np.dot(CD,BA)

#     f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
#     M  = A - alpha*BA
#     N  = C + beta*CD
#     g  = np.einsum('ij,ij->j',M-N,M-N)#np.dot(M-N,M-N)
#     return(f,g)


# def dmin3d_old(A,B,C,D):
#     """
#      dmin3d evaluate the minimal distance between 2 set of segments 

#      this should be vectorized 

#       A : (3xN) initial point segment 1
#       B   (3xN) end point segment 1
#       C   (3xN) starting point segment 2
#       D   (3xN) end point segment 2  
#     """

#     if len(A.shape) ==1 :
#         A=A.reshape(3,1) 
#     if len(B.shape) ==1 :
#         B=B.reshape(3,1) 
#     if len(C.shape) ==1 :
#         C=C.reshape(3,1) 
#     if len(D.shape) ==1 :
#         D=D.reshape(3,1) 

#     AC=C-A
#     CD=D-C
#     BA=A-B

#     u0 = np.einsum('ij,ij->j',AC,AC)#np.dot(AC,AC)
#     u4 = np.einsum('ij,ij->j',BA,BA)#np.dot(BA,BA)
#     u5 = np.einsum('ij,ij->j',CD,CD)#np.dot(CD,CD)
#     u1 = np.einsum('ij,ij->j',BA,AC)#np.dot(BA,AC)
#     u2 = np.einsum('ij,ij->j',CD,AC)#np.dot(CD,AC)
#     u3 = np.einsum('ij,ij->j',CD,BA)#np.dot(CD,BA) 


#     den   = u4*u5-u3*u3
#     alpha = (u2*u3-u1*u5)/(1.*den)
#     beta  = (u1*u3-u2*u4)/(1.*den)
#     #~ print ' dmin**2 ', u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
#     dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5) 
#     #~ print 'dmin', dmin
#     return(alpha,beta,dmin)


def dist(A,B,C,D,alpha,beta):
    """
     distance between AB-CD

     AB (3xN) or (3xNxK) to add a time axis
     CD (3xM) or (3xMxK) to add a time axis


     Parameters
     ----------

        A : (3xN) initial point segment 1
            [or (3xNxK) initial point segment 1 for K realizations]
        B   (3xN) end point segment 1
            [or (3xNxK) end point segment 1 for K realizations]
        C   (3xM) starting point segment 2
            [or (3xMxK) initial point segment 2 for K realizations]
        D   (3xM) end point segment 2  
            [or (3xMxK) end point segment 2 for K realizations]
        alpha : (N x M) parametrization
            [or (NxMxK)] 
        beta : (N x M) parametrization 
            [or (NxMxK)] 
    Returns
    -------
        f : N x M
            [or (NxMxK)] 
        g : N x M
            [or (NxMxK)] 
    """
    if len(A.shape) ==1 :
        A=A.reshape(3,1) 
    if len(B.shape) ==1 :
        B=B.reshape(3,1) 
    if len(C.shape) ==1 :
        C=C.reshape(3,1) 
    if len(D.shape) ==1 :
        D=D.reshape(3,1) 
    
    assert alpha.shape[0] == A.shape[1]
    assert alpha.shape[1] == C.shape[1]
    assert beta.shape[0] == A.shape[1]
    assert beta.shape[1] == C.shape[1]



    # 3 x N x M
    AC = C[:,np.newaxis,:]-A[:,:,np.newaxis]
    # 3 x M 
    CD = D-C
    # 3 x N 
    BA = A-B

    # u0 : N x M
    u0 = np.einsum('ijk...,ijk...->jk...',AC,AC)#np.dot(AC,AC)
    # u4 : N 
    u4 = np.einsum('ij...,ij...->j...',BA,BA)[:,np.newaxis]#np.dot(BA,BA)
    # u5 : M 
    u5 = np.einsum('ij...,ij...->j...',CD,CD)[np.newaxis,:]#np.dot(CD,CD)
    # u1 : N x M
    u1 = np.einsum('ij...,ijk...->jk...',BA,AC)#np.dot(BA,AC)
    # u2 : N x M
    u2 = np.einsum('ik...,ijk...->jk...',CD,AC)#np.dot(CD,AC)
    # u3 : N x M
    u3 = np.einsum('ik...,ij...->jk...',CD,BA)#np.dot(CD,BA) 


    # f : N x M
    f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
    # X : 3 x N x M
    X  = A[:,:,np.newaxis]-alpha[np.newaxis,:,:]*BA[:,:,np.newaxis] # A - alpha*BA
    # Y : 3 x N x M
    Y  = C[:,np.newaxis,:] + beta[np.newaxis,:,:]*CD[:,np.newaxis,:]# C + beta*CD
    # g : N x M
    g =np.einsum('ijk...,ijk...->jk...',X-Y,X-Y)

    return(f,g)

def dmin3d(A,B,C,D):
    """
     dmin3d evaluate the minimal distance between 2 set of segments 
     Note that the number of segment of AB is NOT NECESSARILY the same than CD

     AB (3xN) or (3xNxK) to add a time axis
     CD (3xM) or (3xMxK) to add a time axis

     Parameters
     ----------

        A : (3xN) initial point segment 1
            [or (3xNxK) initial point segment 1 for K realizations]
        B   (3xN) end point segment 1
            [or (3xNxK) end point segment 1 for K realizations]
        C   (3xM) starting point segment 2
            [or (3xMxK) initial point segment 2 for K realizations]
        D   (3xM) end point segment 2  
            [or (3xMxK) end point segment 2 for K realizations]

    Returns
    -------
        alpha 
            parametrization N x M 
            [or (NxMxK) ]
        # beta 
            parametrization N x M 
            [or (NxMxK)]
        dmin 
            minimal distance N x M 
            [or (NxMxK)]
    """

    if len(A.shape) ==1 :
        A=A.reshape(3,1) 
    if len(B.shape) ==1 :
        B=B.reshape(3,1) 
    if len(C.shape) ==1 :
        C=C.reshape(3,1) 
    if len(D.shape) ==1 :
        D=D.reshape(3,1) 


    # 3 x N x M
    AC = C[:,np.newaxis,:]-A[:,:,np.newaxis]
    # 3 x M 
    CD = D-C
    # 3 x N 
    BA = A-B

    # u0 : N x M
    u0 = np.einsum('ijk...,ijk...->jk...',AC,AC)#np.dot(AC,AC)
    # u4 : N 
    u4 = np.einsum('ij...,ij...->j...',BA,BA)[:,np.newaxis]#np.dot(BA,BA)
    # u5 : M 
    u5 = np.einsum('ij...,ij...->j...',CD,CD)[np.newaxis,:]#np.dot(CD,CD)
    # u1 : N x M
    u1 = np.einsum('ij...,ijk...->jk...',BA,AC)#np.dot(BA,AC)
    # u2 : N x M
    u2 = np.einsum('ik...,ijk...->jk...',CD,AC)#np.dot(CD,AC)
    # u3 : N x M
    u3 = np.einsum('ik...,ij...->jk...',CD,BA)#np.dot(CD,BA) 

    # den : N x M
    den   = u4*u5-u3*u3
    # alpha = N x M
    alpha = (u2*u3-u1*u5)/(1.*den)
    # beta = N x M
    beta  = (u1*u3-u2*u4)/(1.*den)
    # dmin : N x M
    dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5)

    return(alpha,beta,dmin)


def segdist(A,B,C,D,hard=True):
    """
     distance between AB-CD
     Note that the number of segment of AB is NOT NECESSARILY the same than CD
     (This function gathers dmin3d and dist from DeuxSeg)

     AB (3xN) or (3xNxK) to add a time axis
     CD (3xM) or (3xMxK) to add a time axis


     Parameters
     ----------

        A : (3xN) initial point segment 1
            [or (3xNxK) initial point segment 1 for K realizations]
        B   (3xN) end point segment 1
            [or (3xNxK) end point segment 1 for K realizations]
        C   (3xM) starting point segment 2
            [or (3xMxK) initial point segment 2 for K realizations]
        D   (3xM) end point segment 2  
            [or (3xMxK) end point segment 2 for K realizations]
        hard : boolean:
            if True when alpha<0 and beta<0 then alpha=0 and beta =0
                    and when alpha>1 and beta>1 then alpha=1 and beta=1
            elif False: use True alpha beta values in computation
    Returns
    -------

        f : N x M
            [or (NxMxK)] 
        g : N x M
            [or (NxMxK)] 
        alpha : (N x M) parametrization
            [or (NxMxK)] 
        beta : (N x M) parametrization 
            [or (NxMxK)] 
        dmin 
            minimal distance N x M 
            [or (NxMxK)]
    """

    if len(A.shape) ==1 :
        A=A.reshape(3,1) 
    if len(B.shape) ==1 :
        B=B.reshape(3,1) 
    if len(C.shape) ==1 :
        C=C.reshape(3,1) 
    if len(D.shape) ==1 :
        D=D.reshape(3,1) 


    # 3 x N x M
    AC = C[:,np.newaxis,:]-A[:,:,np.newaxis]
    # 3 x M 
    CD = D-C
    # 3 x N 
    BA = A-B

    # u0 : N x M
    u0 = np.einsum('ijk...,ijk...->jk...',AC,AC)#np.dot(AC,AC)
    # u4 : N 
    u4 = np.einsum('ij...,ij...->j...',BA,BA)[:,np.newaxis]#np.dot(BA,BA)
    # u5 : M 
    u5 = np.einsum('ij...,ij...->j...',CD,CD)[np.newaxis,:]#np.dot(CD,CD)
    # u1 : N x M
    u1 = np.einsum('ij...,ijk...->jk...',BA,AC)#np.dot(BA,AC)
    # u2 : N x M
    u2 = np.einsum('ik...,ijk...->jk...',CD,AC)#np.dot(CD,AC)
    # u3 : N x M
    u3 = np.einsum('ik...,ij...->jk...',CD,BA)#np.dot(CD,BA) 

    # den : N x M
    den   = u4*u5-u3*u3
    # alpha = N x M
    alpha = (u2*u3-u1*u5)/(1.*den)
    # beta = N x M
    beta  = (u1*u3-u2*u4)/(1.*den)
    # dmin : N x M
    dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5)

    if hard :
        ap = np.where(alpha>1)[0]
        am = np.where(alpha<0)[0]
        bp = np.where(beta>1)[0]
        bm = np.where(beta<0)[0]
        alpha[ap]=np.ones(len(ap))
        alpha[am]=np.zeros(len(am))
        beta[bp]=np.ones(len(bp))
        beta[bm]=np.zeros(len(bm))



    # f : N x M
    f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
    # X : 3 x N x M
    X  = A[:,:,np.newaxis]-alpha[np.newaxis,:,:]*BA[:,:,np.newaxis] # A - alpha*BA
    # Y : 3 x N x M
    Y  = C[:,np.newaxis,:] + beta[np.newaxis,:,:]*CD[:,np.newaxis,:]# C + beta*CD
    # g : N x M
    g =np.einsum('ijk...,ijk...->jk...',X-Y,X-Y)

    return(f,g,alpha,beta,dmin)


# <codecell>
if (__name__=="__main__"):

    A = np.random.rand(3)
    B = np.random.rand(3)
    C = np.random.rand(3)
    D = np.random.rand(3)

    a,b,d=dmin3d_nonvectorized(A,B,C,D)

    if a < 0:
        a = 0

    if a > 1:
        a = 1
    if b < 0:
        b = 0
    if b > 1:
        b = 1

    f, g   = dist_nonvectorized(A,B,C,D,a,b)

    print a,b,d

    print f,g

    print 'sqrt ' ,  np.sqrt(f), np.sqrt(g)



