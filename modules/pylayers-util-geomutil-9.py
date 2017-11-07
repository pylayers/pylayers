import scipy as sp
import numpy as np
from pylayers.util.geomutil import *
from pylayers.util.plotutil import *
import matplotlib.pylab as plt
N = 10
A = sp.rand(2,N)
B = sp.rand(2,N)
C = sp.rand(2,N)
D = sp.rand(2,N)
b1 = intersect(A,B,C,D)
pt1 = A[:,b1]
ph1 = B[:,b1]
pt2 = C[:,b1]
ph2 = D[:,b1]
f1,a1 = displot(pt1,ph1,'r')
f2,a2 = displot(pt2,ph2,'b')
ti = plt.title('test intersect')
A = np.array([[0],[0]])
B = np.array([[1],[1]])
C = np.array([[1],[0]])
D = np.array([[0],[1]])
intersect(A,B,C,D)
# array([ True], dtype=bool)
intersect(A,B,C,D)[0]
# True
