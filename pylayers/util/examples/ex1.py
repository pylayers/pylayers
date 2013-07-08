import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import scipy.linalg as la
from pylayers.util.geomutil import *
import numpy as np 

# cylindre original 
X = np.array([[0,0,0],[0,0,1],[0,0,-1],[1,0,0]]).T

# cylindre final 
Y = np.array([[3,3,3],[3,4,5],[1,3,1],[4,2,3]]).T

B = Y[:,0]
Yc = Y-B[:,np.newaxis]
pX = la.pinv(X)
A = np.dot(Yc,pX)
Xn = np.dot(A,X)+B[:,np.newaxis]

cyl = Geomoff('cylinder')
#cyl.show3()
pt= cyl.loadpt()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xs = pt[:,0]
ys = pt[:,1]
zs = pt[:,2]
#ax.scatter(xs=xs,ys=ys,zs=zs,zdir='z')
#plt.axis('scaled')
#plt.show()


#
#
#

pA = np.array([[3],[3],[3]])
pB = np.array([[10],[10],[10]])
pM = (pA+pB)/2.
T = onbfromaxe(pA,pB)
R = 0.3
Y = np.hstack((pM,pA,pB,pM+R*T[:,0].reshape(3,1),pM+R*T[:,1].reshape(3,1)))
A,B = cylmap(Y)
# 
ptn = np.dot(A,pt.T)+B

xs = ptn[0,:]
ys = ptn[1,:]
zs = ptn[2,:]
ax.scatter(xs=xs,ys=ys,zs=zs,zdir='z')
plt.axis('scaled')
plt.show()
