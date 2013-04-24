from numpy import *
from scipy import *
from matplotlib.pylab import *
from scipy.linalg import *



def dist(A,B,C,D,alpha,beta):
    AC=C-A
    CD=D-C
    BA=A-B
    
    u0 = dot(AC,AC)
    u4 = dot(BA,BA)
    u5 = dot(CD,CD)
    u1 = dot(CD,AC)
    u2 = dot(BA,AC)
    u3 = dot(CD,BA)
    
    f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5
    M  = A - alpha*BA
    N  = C + beta*CD
    g  = dot(M-N,M-N)
    return(f,g)


#A = rand(3)
#B = rand(3)
#C = rand(3)
#D = rand(3)

A =  array([-199.92987677,   19.85458989,   78.87541506])
B =  array([-138.412703  ,   74.27783155,   96.88739959])
C =  array([-170.65697276,   65.39777929,   86.91978279])
D = array([-194.55622849,   65.83291724,   42.14616457])

AC = C-A
CD = D-C
BA = A-B
u0 = dot(AC,AC)
u4 = dot(BA,BA)
u5 = dot(CD,CD)
u1 = dot(CD,AC)
u2 = dot(BA,AC)
u3 = dot(CD,BA)

a = linspace(-2,2,100)
b = linspace(-2.5,2.5,100)
alpha,beta = meshgrid(a,b)
f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5

Z = array([[u4,u3],[u3,u5]])
y = array([-u1,-u2])
print Z 
print y
#print det(Z)





d2,d3 = dist(A,B,C,D,0.5,0.5)
print d2,d3

# <codecell>


# <codecell>

x = solve(Z,y)
yr = dot(Z,x)


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
    
    u0 = dot(AC,AC)
    u4 = dot(BA,BA)
    u5 = dot(CD,CD)
    u1 = dot(CD,AC)
    u2 = dot(BA,AC)
    u3 = dot(CD,BA)
    
    den   = u4*u5-u3*u3
    alpha = (u2*u3-u1*u5)/(1.*den)
    beta  = (u1*u3-u2*u4)/(1.*den)
 #   if (alpha > 0):
 #       alpha = 0
 #   if (alpha >1):
 #       alpha = 1
 #   if (beta > 0):
 #       beta = 0
 #   if (beta >1):
 #       beta = 1  
        
    dmin = sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5) 
   
    return(alpha,beta,dmin)

# <codecell>

a,b,d=dmin3d(A,B,C,D)
print a,b,d

# <codecell>

a

# <codecell>

b

# <codecell>

