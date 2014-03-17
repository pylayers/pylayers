import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import scipy.linalg as la

def dist(A,B,C,D,alpha,beta):
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

def dmin3d(A,B,C,D):
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

# <codecell>
if (__name__=="__main__"):

    A = np.random.rand(3)
    B = np.random.rand(3)
    C = np.random.rand(3)
    D = np.random.rand(3)

    a,b,d=dmin3d(A,B,C,D)

    if a < 0:
        a = 0

    if a > 1:
        a = 1
    if b < 0:
        b = 0
    if b > 1:
        b = 1

    f, g   = dist(A,B,C,D,a,b)

    print a,b,d

    print f,g

    print 'sqrt ' ,  np.sqrt(f), np.sqrt(g)



