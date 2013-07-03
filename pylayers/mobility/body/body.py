import numpy as np
g = plt.figure()                                                                                     
                                                                                import scipy.stats as sp
from pylayers.mobility.body import c3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import  pdb as pdb
import pylayers.util.plotutil as pltu
import doctest


def Basis_generation(A, B):
    v = (B - A) / np.linalg.norm(B - A)
    test = True
    while test:
        random_vector = np.random.random(3)
        u = random_vector - np.dot(np.dot(random_vector, v), v)
        if(u.any()):
            test = False
    u = u / np.linalg.norm(u)
    w = np.cross(u, v)
    return u, v, w


def ChangeBasis(u0, v0, w0, v1):

    # Rotation autour de l'axe w

    v2 = v1 - np.dot(np.dot(v1, w0), w0)  # projection de v1 sur le plan (u,v)
    v2 = v2 / np.linalg.norm(v2)
    c = np.dot(v2, u0)
    s = np.dot(v2, v0)
    u1 = np.dot(s, u0) - np.dot(c, v0)
    u1 = u1 / np.linalg.norm(u1)
    w1 = np.cross(u1, v1)
    return u1, v1, w1


def dist(A, B):
    """
    evaluate the distance between two points 1 and B
    """

    d = np.sqrt((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
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

    def LoadMotion(self, filename='07_01.c3d', nframes=126):

        s, p, f = c3d.read_c3d(filename)

        # self.d 3 x np x nf
        self.d = np.ndarray(shape=(3, 15, np.shape(f)[0]))
        ind = []
        for i in range(len(self.marker_set)):
            ind.append(p.index(s[0] + self.marker_set[i]))

        # f.T : 3 x np x nf
        self.d = f[0:nframes, ind, :].T
        self.g.pos = {}
        for i in range(15):
            self.g.pos[i] = (self.d[1, i, 0], self.d[2, i, 0])

    def CylinderModel(self, nc=10):

        """

        nc  =  cylinder number
        c : array(shape  =  (nc,8)), Cylinder Id , A coordiante, B Coordinate , cylinder radius

        """

        self.c = np.ndarray(shape=(nc, 8, 126))
        i = 0
        self.c[i, 0] = i
        self.c[i, 1:4] = (self.d[:, 9] + self.d[:, 10]) / 2
        self.c[i, 4:6] = self.c[i, 1:3]
        self.c[i, 6] = (self.d[2, 3] + self.d[2, 4]) / 2
        self.c[i, 7] = dist(self.d[:, 9], self.d[:, 10]) / 2

        i = 1
        self.c[i, 0] = i
        self.c[i, 1:4] = self.c[0, 4:7]
        self.c[i, 4:6] = self.c[0, 4:6]
        self.c[i, 6] = self.d[2, 2, ]
        self.c[i, 7] = 2

        i = 2
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 6]
        self.c[i, 4:7] = self.d[:, 4]
        self.c[i, 7] = 5

        i = 3
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 5]
        self.c[i, 4:7] = self.d[:, 3]
        self.c[i, 7] = 5

        i = 4
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 8]
        self.c[i, 4:7] = self.d[:, 6]
        self.c[i, 7] = 5

        i = 5
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 7]
        self.c[i, 4:7] = self.d[:, 5]
        self.c[i, 7] = 5

        i = 6
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 12, ]
        self.c[i, 4:7] = self.d[:, 10]
        self.c[i, 7] = 5

        i = 7
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 11]
        self.c[i, 4:7] = self.d[:, 9]
        self.c[i, 7] = 5

        i = 8
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 14]
        self.c[i, 4:7] = self.d[:, 12]
        self.c[i, 7] = 5

        i = 9
        self.c[i, 0] = i
        self.c[i, 1:4] = self.d[:, 13]
        self.c[i, 4:7] = self.d[:, 11]
        self.c[i, 7] = 5

    def cylinder_basis0(self, frameId=0):
        nc = self.c.shape[0]
        self.basis0 = np.ndarray(shape=(nc, 9))
        for i in range(nc):
            Ai = self.c[i, 1:4, frameId]
            Bi = self.c[i, 4:7, frameId]
            u, v, w = Basis_generation(Ai, Bi)
            #pdb.set_trace()
            self.basis0[i, 0:3] = u
            self.basis0[i, 3:6] = v
            self.basis0[i, 6:] = w

    def cylinder_basis_k(self, frameId):

        nc = self.c.shape[0]
        self.basisk = np.ndarray(shape=(nc, 9))
        for i in range(nc):
            u0 = self.basis0[i, 0:3]
            v0 = self.basis0[i, 3:6]
            w0 = self.basis0[i, 6:]
            v1 = self.c[i, 4:7, frameId] - self.c[i, 1:4, frameId]
            v1 = v1 / np.linalg.norm(v1)
            uk, vk, wk = ChangeBasis(u0, v0, w0, v1)
            self.basisk[i, 0:3] = uk
            self.basisk[i, 3:6] = vk
            self.basisk[i, 6:] = wk

    def cyl_antenna(self, cylinderId, l, alpha, frameId=0):
        r = self.c[cylinderId, 7, frameId]

        x = r * np.cos(alpha)
        y = r * np.sin(alpha)
        z = l
        if frameId == 0:
            u0 = self.basis0[cylinderId, 0:3]
            v0 = self.basis0[cylinderId, 3:6]
            w0 = self.basis0[cylinderId, 6:]

        else:
            self.cylinder_basis_k(frameId)
            u0 = self.basisk[cylinderId, 0:3]
            v0 = self.basisk[cylinderId, 3:6]
            w0 = self.basisk[cylinderId, 6:]
        #~ #pdb.set_trace()
        self.ant = x.reshape((len(x)), 1) * u0 + y.reshape(
            (len(y)), 1) * w0 + z.reshape((len(z)), 1) * v0


def translate(cycle, new_origin):
    """ rotate a cycle of frames by an angle alpha

    Parameters
    ----------

    cycle :  np.array
            3 x np x nf
    alpha : float
        angle in radians


    Returns
    -------

    cycle modified : np.array
        3 x np x nf
    """


    cycle_tr = np.ndarray(shape=cycle.shape)
    old_origin = cycle[:, 0, 0]
    cycle_tr = cycle + new_origin.reshape((3, 1, 1)) - \
        old_origin.reshape((3, 1, 1))

    return cycle_tr


def rotation(cycle, alpha=np.pi/2):
    """ rotate a cycle of frames by an angle alpha

    Parameters
    ----------

    cycle :  np.array
            3 x np x nf
    alpha : float
        angle in radians


    Returns
    -------

    cycle modified : np.array
        3 x np x nf
    """


    cycle_rot = np.ndarray(shape=cycle.shape)
    cycle_rot[0, :, :] = (cycle[0, :, :]) * np.cos(alpha) + (cycle[1, :, :]) * np.sin(alpha)
    cycle_rot[1, :, :] = -(cycle[0, :, :]) * np.sin(alpha) + (cycle[1, :, :]) * np.cos(alpha)
    cycle_rot[2, :, :] = cycle[2, :, :]
    cycle_rot = translate(cycle_rot, cycle[:, 0, 0])

    return cycle_rot


def Global_Trajectory(cycle, traj):
    """

    Parameters
    ----------

    cycle :  walking step cycle (2 step), shape = (3,npoints  = 15, nframes = 126)
    traj  : trajectory described by the gravity center, shape =(3,nposition)

    We assume that the body moves straight between two successive positions

    Returns
    -------

    data : list
        list
    """

    data = []

    fr_start_index = 0
    ref_fr = cycle[:, :, 0]
    vect_ortho = ref_fr[:, 3] - ref_fr[:, 4]
    vect_ortho = vect_ortho / np.linalg.norm(vect_ortho)
    v = np.random.random(3)
    v = v - np.dot(np.dot(v, vect_ortho), vect_ortho)
    v[2] = 0
    vect_ant = v / np.linalg.norm(v)

    for i in range(1, traj.shape[1]):

        print 'i = ', i

        vect_depl = traj.T[i] - traj.T[i - 1]
        vect_depl = vect_depl / np.linalg.norm(vect_depl)

        alpha = np.arccos(np.dot(vect_ant, vect_depl))

        cycle_i = rotation(cycle, alpha)

        dist_inter = dist(traj.T[i], traj.T[i - 1])
        Nfr = int(dist_inter * 126.0 / 140.0)
        cycle_i = translate(cycle_i, traj.T[i - 1])
        if Nfr < 126:
            if fr_start_index + Nfr < 126:
                data.append(cycle_i[:, :, fr_start_index:fr_start_index + Nfr])
                fr_start_index = fr_start_index + Nfr
            else:
                data.append(cycle_i[:, :, fr_start_index:])
                cycle_i = translate(cycle_i, cycle_i[:, 0, -1] + vect_depl)
                data.append(cycle_i[:, :, 0:fr_start_index + Nfr - 126])
                fr_start_index = fr_start_index + Nfr - 126

    return data


if __name__ == '__main__':
    plt.ion()
    doctest.testmod()

