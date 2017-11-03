from pylayers.util.geomutil import *
import numpy as np 
import matplotlib.pyplot as plt 
import pdb
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
Ndim = 2
np.random.seed(4)
Npoints = 100000

if Ndim>2:
    ax = fig.add_subplot(111, projection='3d')
    p = 10*np.random.rand(Npoints,3)-5
    #p = np.array([[1,2,0.1],[1,3,10],[1,.2,0.4]])
    v = np.array([[1,0,0],[1,1,0],[0,0,1]]).T
    apex = np.array([0,0,0]) 
else:
    apex = np.array([0,0]) 
    v = np.array([[1,1],[-1,1]]).T
    p = 10*np.random.rand(Npoints,2)-5
    ax = fig.add_subplot(111)
    ax.scatter(p[:,0],p[:,1],s=100,c='blue',alpha=0.05,linewidth=0)

bb = are_points_inside_cone(p,apex,v,radius=np.inf)

if Ndim>2:
    ax.scatter(p[~bb,0],p[~bb,1],p[~bb,2],s=100,c='blue',alpha=0.4,linewidth=0)
    ax.scatter(p[bb,0],p[bb,1],p[bb,2],s=100,c='red',alpha=1,linewidth=0)
    #ax.plot([0,10*v[0,0]],[0,10*v[0,1]],[0,10*v[0,2]],color='blue')
    #ax.plot([0,10*v[1,0]],[0,10*v[1,1]],[0,10*v[1,2]],color='blue')
    #ax.plot([0,10*v[2,0]],[0,10*v[2,1]],[0,10*v[2,2]],color='blue')
else:
    ax.scatter(p[bb,0],p[bb,1],s=100,c='red',alpha=1,linewidth=0)
#plt.figure()
#plt.scatter(p[bb,0],p[bb,1],s=100,c='red',alpha=1,linewidth=0)

#plt.axis('equal')
plt.show()
