# >>>from pylayers.util.plotutil import *
import scipy as sp
import numpy as np
from pylayers.util.geomutil import *
from pylayers.util.plotutil import *
import matplotlib.pylab as plot
N = 20
A = sp.rand(2,N)
B = sp.rand(2,N)
C = np.array(([0.5,0.5])).reshape(2,1)
left=isleft(A,B,C)
il = np.where(left)[0]
inl = np.where(~left)[0]
plt.scatter(C[0],C[1],color='b',s=10)
displot(A[:,il],B[:,il],arrow=True,color='g')
displot(A[:,inl],B[:,inl],arrow=True,color='r')
