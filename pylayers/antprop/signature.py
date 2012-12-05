#-*- coding:Utf-8 -*-
import doctest
import numpy as np
import scipy as sp
import scipy.linalg as la
import pdb
#import GrRay3D
#import Graph
#import pylayers.antprop.slab
#import pylayers.gis.layout
#import pylayers.util.geomutil as geu
import matplotlib.pyplot as plt


class Signature(object):
    """ class Signature

    A signature contains two lists

    seq : list of interaction numbers
    typ : list of interaction type
    """
    def __init__(self, seq):
        """
        pa  : tail point of intercation segment
        pb  : head point of intrcation segement
        pc  : middle point  of interaction segment
        typ : type of interaction 1-R 2-T 3-D
        seq : sequence of interaction point (edges (>0)  or vertices (<0)
        """
        self.seq = seq
        self.typ = []

    def info(self):
        """
        """
        for k in self.__dict__.keys():
            print k, ':', self.__dict__[k]

    def split(self):
        """
        split signature
        """
        pass

    def ev(self, L):
        """  evaluation of Signature

        Parameters
        ----------
            L : Layout

        Notes
        -----
        Le type des interactions d extremite reste indetermine a ce stade
        """
        N = len(self.seq)
        self.pa = np.zeros((2, N))  # tail
        self.pb = np.zeros((2, N))  # head
        self.pc = np.zeros((2, N))  # center
        self.typ = np.zeros(N)
        self.norm = np.zeros((2, N))

        for n in range(N):
            k = self.seq[n]
            if k > 0:  # segment
                ta, he = L.Gs.neighbors(k)
                norm1 = L.Gs.node[k]['norm']
                norm = np.array([norm1[0], norm1[1]])
                self.pa[:, n] = np.array(L.Gs.pos[ta])
                self.pb[:, n] = np.array(L.Gs.pos[he])
                self.pc[:, n] = np.array(L.Gs.pos[k])
                self.norm[:, n] = norm
                self.typ[n] = 1
            else:      # node
                pa = np.array(L.Gs.pos[k])
                norm = array([0, 0])
                self.pa[:, n] = pa
                self.pb[:, n] = pa
                self.pc[:, n] = pa
                self.norm[:, n] = norm
                self.typ[n] = 3
        #
        #  vecteurs entre deux points adjascents de la signature
        #
        #self.v   = s.pc[:,1:]-s.pc[:,:-1]
        #self.vn  = self.v /np.sqrt(np.sum(self.v*self.v,axis=0))
        #u1       = np.sum(self.norm*self.vn[:,0:-1],axis=0)
        #u2       = np.sum(self.norm*self.vn[:,1:],axis=0)
        #self.typ = sign(u1*u2)
        #return(vn)
        #return(typ)

    def evtx(self, L, tx):
        """ evtx

        Parameters
        ----------
            L  : Layout
            tx : np.array (2xN)
            rx : np.array (2xM)

        """
        self.pa = tx.reshape(2, 1)
        self.pb = tx.reshape(2, 1)
        self.pc = tx.reshape(2, 1)
        self.typ = np.array([0])
        for k in self.seq:
            if k > 0:
                ta, he = L.Gs.neighbors(k)
                norm1 = L.Gs.node[k]['norm']
                norm = np.array([norm1[0], norm1[1]]).reshape(2, 1)
                pa = np.array(L.Gs.pos[ta]).reshape(2, 1)
                pb = np.array(L.Gs.pos[he]).reshape(2, 1)
                pc = np.array(L.Gs.pos[k]).reshape(2, 1)
                self.pa = np.hstack((self.pa, pa))
                self.pb = np.hstack((self.pb, pb))
                self.pc = np.hstack((self.pc, pc))
                try:
                    self.norm = np.hstack((self.norm, norm))
                except:
                    self.norm = norm
                self.typ = np.hstack((self.typ, np.array([1])))
            else:
                pa = np.array(L.Gs.pos[k]).reshape(2, 1)
                norm = np.array([0, 0]).reshape(2, 1)
                self.pa = np.hstack((self.pa, pa))
                self.pb = np.hstack((self.pb, pa))
                self.pc = np.hstack((self.pc, pa))
                try:
                    self.norm = np.hstack((self.norm, norm))
                except:
                    self.norm = norm
                self.typ = np.hstack((self.typ, np.array([3])))
        self.pa = np.hstack((self.pa, rx.reshape(2, 1)))
        self.pb = np.hstack((self.pb, rx.reshape(2, 1)))
        self.pc = np.hstack((self.pc, rx.reshape(2, 1)))
        self.typ = np.hstack((self.typ, np.array([0])))
        #
        #  vecteur entre deux points adjascents de la signature
        #
        self.v = s.pc[:, 1:] - s.pc[:, :-1]
        self.vn = self.v / np.sqrt(sum(self.v * self.v, axis=0))
        u1 = sum(self.norm * self.vn[:, 0:-1], axis=0)
        u2 = sum(self.norm * self.vn[:, 1:], axis=0)
        self.typ = sign(u1 * u2)
        #return(vn)
        #return(typ)

    def image(self, tx):
        """

        Parameters
        ----------
            tx : np.array (2xN)

        Returns
        -------
            M : Matrix for ray calculation  (to be stored)
        Notes
        ------

        pa  : first point of segment or diffraction point
        pb  : last point of segment or diffraction point
        typ :
        """
        pa = self.pa
        pb = self.pb
        typ = self.typ
        # number of interactions
        N = np.shape(pa)[1]
        # detect diffraction
        usig = np.nonzero(typ == 3)[0]
        # detect transmission
        vsig = np.nonzero(typ == 2)[0]
        dsig = np.ones(N)
        rsig = np.ones(N)
        # dsig = 1 R|T
        # dsig = 0 D
        rsig[usig] = 0
        rsig[vsig] = 0
        dsig[usig] = 1

        Ntx = np.shape(tx)
        pab = pb - pa
        pabo = np.array([pa[1, :] - pb[1, :], pb[0, :] - pa[0, :]])
        alpha = np.sum(pab * pab, axis=0)

        a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
        b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
        c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
                            pa[1, :] * (pa[1, :] - pb[1, :]) * (pb[0, :] - pa[0, :]))
        d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
                            pa[0, :] * (pa[1, :] - pb[1, :]) * (pb[0, :] - pa[0, :]))

        I2 = np.eye(2)
        Z2 = np.zeros((2, 2))
        S = np.zeros((N, 2, 2))
        D = np.zeros((N, 2, 2))
        D[:, 0, 0] = 1
        D[:, 1, 1] = 1
        S[:, 0, 0] = -a
        S[:, 0, 1] = b
        S[:, 1, 0] = b
        S[:, 1, 1] = a
        #
        # La matrice contient les termes à partir de l'interaction 1
        #
        a1 = a[1::] * rsig[1::]
        a2 = -np.ones(N - 1) * (1 - rsig[1::])
        b1 = b[1::] * rsig[1::]
        c1 = c[1::] * rsig[1::]
        d1 = d[1::] * rsig[1::]
        #
        # 1 sur la diagonale sauf pour les diffractions
        #
        # N : Nombre d'interactions
        # A est une matrice 2*N x 2*N
        #
        A = np.diag(np.kron(dsig, [1, 1]))
        # dm1  : Sous diagonale 1 : 2*N - 1
        dm1 = np.zeros(2 * N - 1)
        #  0 b 0 b 0 b
        kdm1 = range(1, 2 * N - 1, 2)
        dm1[kdm1] = b1
        # dm2 : Sous diagonale 2 : 2*(N-1)
        # -a a -a a -a a  ...
        pdb.set_trace()
        dm2 = -np.kron(a1, np.array([1, 0])) + np.kron(a1, np.array([0, 1]))
        dm2 = dm2 + np.kron(a2, [1, 1])
        # annulation de la sous diganale 2 pour la diffraction
        # if diffraction
        if len(usig) > 0:
            v = np.kron(2 * usig, np.array(
                [1, 0])) + np.kron(2 * usig + 1, np.array([0, 1]))
            dm2[v] = 0
        # dm3 : sous diagonale 3 : 2*N-3
        #
        dm3 = np.zeros(2 * N - 3)
        kdm3 = range(0, 2 * N - 3, 2)
        dm3[kdm3] = b1
        #A     = np.eye(2*N)
        Am1 = np.diagflat(dm1, -1)
        Am2 = np.diagflat(dm2, -2)
        Am3 = np.diagflat(dm3, -3)
        A = A + Am1 + Am2 + Am3
        #IA    = la.inv(A)
        #
        # Evaluate y
        #
        if typ[0] == 1:
            vc0 = np.array([c[0], d[0]])
            v0 = np.dot(-S[0, :, :], tx) + vc0   # TODO vérifier signe de S
        if typ[0] == 2:
            v0 = tx
        if typ[0] == 3:
            v0 = pa[:, 0]
        #c0    = c[0]+tx[0]*a[0]-tx[1]*b[0]
        #d0    = d[0]-tx[1]*a[0]-tx[0]*b[0]

        #y    = np.hstack((np.array([c0,d0]),np.kron(c1,np.array([1,0]))+np.kron(d1,np.array([0,1]))))
        y = np.hstack((v0, np.kron(
            c1, np.array([1, 0])) + np.kron(d1, np.array([0, 1]))))
        #v1   = np.hstack((0,0,y[0:-2]))
        #v2   = np.hstack((0,0,np.kron(d[0:-1],np.array([1,0]))+np.kron(c[0:-1],np.array([0,1]))))
        #alpha  = np.hstack((0,0,np.kron(a1,np.array([1,-1]))))
        #beta   = np.hstack((0,0,np.kron(b1,np.array([1,1]))))
        x = la.solve(A, y)
        #xv   = x * np.kron(2./alpha,np.array([1,1]))
        M = np.vstack((x[0::2], x[1::2]))
        #return(M,A,y,dm1,dm2,dm3,a1,b1,c1,d1)
        return(M)

    def backtrace(self, tx, rx, M):
        """ backtracing step

        Parameters
        ----------
        tx :  transmitter
        rx :  receiver
        M  :  M

        Returns

        Examples
        --------

        .. plot::
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> from pylayers.gis.layout import *
            >>> from pylayers.antprop.signature import *
            >>> L = Layout()
            >>> L.loadstr('exemple.str')
            >>> L.buildGt()
            >>> L.buildGr()
            >>> seq = [1,5,1]
            >>> s = Signature(seq)
            >>> tx = np.array([4,-1])
            >>> rx = np.array([1,1])
            >>> s.ev(L)
            >>> M = s.image(tx)
            >>> Y = s.backtrace(tx,rx,M)
            >>> fig = plt.figure()
            >>> ax = fig.add_subplot(111)
            >>> l1 = ax.plot(tx[0],tx[1],'or')
            >>> l2 = ax.plot(rx[0],rx[1],'og')
            >>> l3 = ax.plot(M[0,:],M[1,:],'ob')
            >>> l4 = ax.plot(Y[0,:],Y[1,:],'xk')
            >>> ray = np.hstack((np.hstack((tx.reshape(2,1),Y)),rx.reshape(2,1)))
            >>> l5 = ax.plot(ray[0,:],ray[1,:],color='#999999',alpha=0.6,linewidth=0.6)
            >>> fig,ax = L.showGs(fig,ax)
            >>> plt.show()

        """

        pa = self.pa
        pb = self.pb
        typ = self.typ

        N = np.shape(pa)[1]
        I2 = np.eye(2)
        z0 = np.zeros((2, 1))

        pkm1 = rx.reshape(2, 1)
        Y = pkm1
        k = 0
        beta = .5
        while (((beta <= 1) & (beta >= 0)) & (k < N)):
            l0 = np.hstack((I2, pkm1 - M[:, -(k + 1)].reshape(2, 1), z0))
            l1 = np.hstack((I2, z0, pa[:, -(k + 1)].reshape(
                2, 1) - pb[:, -(k + 1)].reshape(2, 1)))
            T = np.vstack((l0, l1))
            yk = np.hstack((pkm1[:, 0].T, pa[:, -(k + 1)].T))
            xk = la.solve(T, yk)
            pkm1 = xk[0:2].reshape(2, 1)
            gk = xk[2::]
            beta = gk[1]
            Y = np.hstack((Y, pkm1))
            k = k + 1

        if ((k == N) & ((beta > 0) & (beta < 1))):
            Y = np.hstack((Y, tx.reshape(2, 1)))

        return(Y)

#pa    = np.array([[2,6],[8,7]])
#pb    = np.array([[5,8],[8,2]])
if __name__ == "__main__":
    doctest.testmod()
