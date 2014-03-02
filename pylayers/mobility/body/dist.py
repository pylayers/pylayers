import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
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

    Return
    ------
    f
    g

    """
    AC=C-A
    CD=D-C
    BA=A-B

    u0 = np.dot(AC,AC)
    u4 = np.dot(BA,BA)
    u5 = np.dot(CD,CD)
    u1 = np.dot(CD,AC)
    u2 = np.dot(BA,AC)
    u3 = np.dot(CD,BA)

    f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
    M  = A - alpha*BA
    N  = C + beta*CD
    g  = np.dot(M-N,M-N)
    return(f,g)


#A = rand(3)
#B = rand(3)
#C = rand(3)
#D = rand(3)

A =  np.array([-199.92987677,   19.85458989,   78.87541506])
B =  np.array([-138.412703  ,   74.27783155,   96.88739959])
C =  np.array([-170.65697276,   65.39777929,   86.91978279])
D =  np.array([-194.55622849,   65.83291724,   42.14616457])

AC = C-A
CD = D-C
BA = A-B
u0 = np.dot(AC,AC)
u4 = np.dot(BA,BA)
u5 = np.dot(CD,CD)
u1 = np.dot(CD,AC)
u2 = np.dot(BA,AC)
u3 = np.dot(CD,BA)

a = np.linspace(-2,2,100)
b = np.linspace(-2.5,2.5,100)
alpha,beta = np.meshgrid(a,b)
f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5

Z = np.array([[u4,u3],[u3,u5]])
y = np.array([-u1,-u2])
print Z
print y
#print det(Z)





d2,d3 = dist(A,B,C,D,0.5,0.5)
print d2,d3

# <codecell>


# <codecell>

x = la.solve(Z,y)
yr = np.dot(Z,x)


print yr
print "alpha : ",x[0]
print "beta : ",x[1]

al = x[0]
be = x[1]
if al<0:
    al = 0
if al>1:
    al = 1

if be<0:
    be = 0
if be>1:
    be = 1

M = A + al*(B-A)
N = C + be*(D-C)

MO = A - 2*(B-A)
MI = A + 2*(B-A)
NO = C - 2*(D-C)
NI = C + 2*(D-C)

print "alpha : ",al
print "beta : ",be

# <codecell>



# <codecell>


# <markdowncell>

# \\( \pmatrix{\alpha \\\\ \beta }= \frac{1}{u_4 u_5 - u_3^2} \pmatrix{u_2 u_3 - u_1 u_5 \\\\ u_1 u_3-u_2 u_4 }\\)

# <codecell>

def dmin3d(A,B,C,D):
    """
     dmin3d evaluate the minimal distance between 2 set of segments

     this should be vectorized

      A : (3xN) initial point segment 1
      B   (3xN) end point segment 1
      C   (3xN) starting point segment 2
      D   (3xN) starting point segment 2
    """

    AC=C-A
    CD=D-C
    BA=A-B

    u0 = np.dot(AC,AC)
    u4 = np.dot(BA,BA)
    u5 = np.dot(CD,CD)
    u1 = np.dot(CD,AC)
    u2 = np.dot(BA,AC)
    u3 = np.dot(CD,BA)

    den   = u4*u5-u3*u3
    alpha = (u2*u3-u1*u5)/(1.*den)
    beta  = (u1*u3-u2*u4)/(1.*den)

    dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5) 

    return(alpha,beta,dmin)


