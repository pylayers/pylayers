import numpy as np
import scipy.stats as sp
from pylayers.mobility.body import c3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
from visual import *
import  pdb as pdb
import pylayers.util.plotutil as pltu


def dist (A,B):
	"""
	evaluate the distance between two points 1 and B
	"""
	
	d  = sqrt((A[0] -B[0])**2 +(A[1] -B[1])**2+ (A[2] -B[2])**2)
	return d 
	

class Body(object):
    """
       Class to manage c3d files
    """

    def __init__(self):

        self.g = nx.Graph()
        nodes_ID = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
        self.marker_set = ['STRN', 'CLAV', 'RFHD', 'RSHO', 'LSHO', 'RELB', 'LELB', 'RWRB', 'LWRB', 'RFWT', 'LFWT', 'RKNE', 'LKNE', 'RANK', 'LANK']
        self.g.add_nodes_from(nodes_ID)
        self.g.add_edge(0, 1)
        self.g.add_edge(0, 9)
        self.g.add_edge(0, 10)
        self.g.add_edge(1, 2)
        self.g.add_edge(1, 3)
        self.g.add_edge(1, 4)
        self.g.add_edge(3, 5)
        self.g.add_edge(4, 6)
        self.g.add_edge(5, 7)
        self.g.add_edge(6, 8)
        self.g.add_edge(9, 11)
        self.g.add_edge(10, 12)
        self.g.add_edge(11, 13)
        self.g.add_edge(12, 14)

    def LoadMotion(self, filename='07_01.c3d'):

        s,p,f = c3d.read_c3d(filename)

        # self.d 3 x np x nf
        self.d = np.ndarray(shape=(3, 15, np.shape(f)[0]))
        ind = []
        for i in range(len(self.marker_set)):
            ind.append(p.index(s[0] + self.marker_set[i]))

        # f.T : 3 x np x nf 
        self.d = f[0:nframes, ind, :].T
        self.g.pos={}
        for i in range(15):
            self.g.pos[i]=(self.d[1,i,0],self.d[2,i,0])

    def CylinderModel(self, nc=10, frameID  = 0 ):
		
		"""
		
		nc  =  cylinder number 
		c : array(shape  =  (nc,4)), Cylinder Id , A coordiante, B Coordinate , cylinder radius 
		
		"""
		
		self.c  = np.ndarray(shape = (nc,8))
		i = 0
		self.c[i,0] = i
		self.c[i,1:4] = list((self.d[:,9,frameID] +self.d[:,10,frameID])/2 )
		self.c[i,4:7] = (self.d[:,3,frameID] +self.d[:,4,frameID])/2
		self.c[i,7] = dist(self.d[:,3,frameID],self.d[:,4,frameID])/2
		
		i = 1
		self.c[i,0] = i
		self.c[i,1:4] = (self.d[:,3,frameID] +self.d[:,4,frameID])/2
		self.c[i,4:6] = self.c[0,4:6]
		self.c[i,6] = self.d[2,2,frameID]
		self.c[i,7] = 5
		
		i = 2 
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,6,frameID]
		self.c[i,4:7] = self.d[:,4,frameID]
		self.c[i,7] = 5
		
		i = 3
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,5,frameID]
		self.c[i,4:7] = self.d[:,3,frameID]
		self.c[i,7] = 5
		
		i = 4
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,8,frameID]
		self.c[i,4:7] = self.d[:,6,frameID]
		self.c[i,7] = 5
		
		i = 5 
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,7,frameID]
		self.c[i,4:7] = self.d[:,5,frameID]
		self.c[i,7] = 5
			 
			 
		i = 6
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,12,frameID]
		self.c[i,4:7] = self.d[:,10,frameID]
		self.c[i,7] = 5
		
		i = 7 
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,11,frameID]
		self.c[i,4:7] = self.d[:,9,frameID]
		self.c[i,7] = 5
		
		i = 8
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,14,frameID]
		self.c[i,4:7] = self.d[:,12,frameID]
		self.c[i,7] = 5
		
		i = 9
		self.c[i,0] = i
		self.c[i,1:4] = self.d[:,13,frameID]
		self.c[i,4:7] = self.d[:,11,frameID]
		self.c[i,7] = 5

       


if __name__ == '__main__':

    nframes = 126
    B = Body()
    B.LoadMotion()
    c10_15 = B.d
    #nx.draw(B.g)
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    frameID = 75
    ax.scatter(c10_15[0, :, frameID], c10_15[1, :, frameID], c10_15[2, :, frameID])
    ax.axis('scaled')
    B.CylinderModel(frameID)
    fig  = plt.figure()
    for i in range(2,B.c.shape[0]):
		pltu.cylinder(fig,B.c[i,1:4],B.c[i,4:7],B.c[i,7])	

    plt.axis('scaled')
    plt.show()
