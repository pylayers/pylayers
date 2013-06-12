import numpy as np
import scipy.stats as sp
from pylayers.mobility.body import c3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
#from visual import *
import  pdb as pdb


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

    def CylinderModel(self):

        return 0


if __name__ == '__main__':

    nframes = 126
    B = Body()
    B.LoadMotion()
    c10_15 = B.d
    #nx.draw(B.g)
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    frameID = 30
    ax.scatter(
        c10_15[0, :, frameID], c10_15[1, :, frameID], c10_15[2, :, frameID])
    ax.axis('scaled')
    #ax.scatter(c10_15[frameID,:,0],c10_15[frameID,:,1],c10_15[frameID,:,2] )
    pointID = 13
    plt.figure()
    plt.plot(c10_15[0, pointID].T, '*--', label='x')
    plt.plot(c10_15[1, pointID].T, '*--', label='y')
    plt.plot(c10_15[2, pointID].T, '*--', label='z')
    plt.legend()

    plt.show()
