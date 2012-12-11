#-*- coding:Utf-8 -*-
import doctest
import numpy as np
#import scipy as sp
import scipy.linalg as la
import pdb
import networkx as nx
#import GrRay3D
#import Graph
#import pylayers.antprop.slab
#import pylayers.gis.layout
#import pylayers.util.geomutil as geu
#import pylayers.util.graphutil as gph
import pylayers.util.pyutil as pyu
import matplotlib.pyplot as plt

class Signatures(object):
    """
    gather all signatures from a layout given tx and rx
    Attributes
    ----------
        L : gis.Layout
        pTx : numpy.ndarray
            position of Tx
        pRx : numpy.ndarray
            position of Rx
    """

    def __init__(self, L, pTx, pRx):
        """
        """
        self.L = L
        self.pTx = pTx
        self.pRx = pRx

    def info(self):
        """
        """
        print "Signatures for scenario defined by :"
        print "Layout"
        print "======"
        self.L.info()
        print "================================"
        print "Transmitter position: ", self.pTx
        print "Receiver position: ", self.pRx

    def get_sigarr(self):
        """
        get signatures (in one array) between iTx and iRx
        signatures are separated by zeros
        Parameters
        ----------
            L : Layout
            iTx : integer
            iRx : integer
        Returns
        -------
            sigarr = numpy.ndarray

        Warnings
        --------
        This a temporary function
            There is some algorithmic work to find the best way to determine
            signature
            T4 : limit the ndt to only edges and nodes in visibility from Tx

        """
        # Here we take all the vnodes >0  from the room
        #
        # Practically those list of nodes should depend on pTx , pRx
        #
        try:
            self.L.Gi
        except:
            self.L.build()

        NroomTx = self.L.pt2ro(self.pTx)
        NroomRx = self.L.pt2ro(self.pRx)

        if not self.L.Gr.has_node(NroomTx) or not self.L.Gr.has_node(NroomRx):
            raise AttributeError('Tx or Rx is not in Gr')

        #
        # .. todo:: modifier inter afin de ne pas retenir les points non
        # diffractants
        #
        ndt = self.L.Gt.node[self.L.Gr.node[NroomTx]['cycle']]['inter']
        ndr = self.L.Gt.node[self.L.Gr.node[NroomRx]['cycle']]['inter']
        ntr = np.intersect1d(ndt,ndr)
        sigarr = np.array([]).reshape(2, 0)
        for nt in ndt:
            for nr in ndr:
                addpath = False
                if (type(nt) != type(nr)):
                    try:
                        path = nx.dijkstra_path(self.L.Gi, nt, nr)
                        addpath = True
                    except:
                        pass
                        #print 'no path between ',nt,nr
                elif (nt != nr):
                    try:
                        path = nx.dijkstra_path(self.L.Gi, nt, nr)
                        addpath = True
                    except:
                        pass
                        #print 'no path between ',nt,nr
                else:
                    addpath = True
                    path = [nt]
                if addpath:
                    sigarr = np.hstack((sigarr, np.array([[0], [0]])))
                    for interaction in path:
                        it = eval(interaction)
                        if type(it) == tuple:
                            sigarr = np.hstack((sigarr,
                                    np.array([[it[0]], [1]])))
                        elif it < 0:
                            sigarr = np.hstack((sigarr,
                                    np.array([[it], [3]])))
                        else:
                            sigarr = np.hstack((sigarr,
                                    np.array([[it], [2]])))
        sigarr = np.hstack((sigarr, np.array([[0], [0]])))
        return sigarr

    def sigs2rays(self, L, pTx, pRx, sig_arr):
        """
        from signature arrays to rays
        Parameters
        ----------
            L : Layout
            pTx : ndarray
                2D transmitter position
            pRx : ndarray
                2D receiver position
            sig_arr : ndarray
        Returns
        -------
            M : numpy.ndarray
            Y : numpy.ndarray
        """
        zsig = np.array(np.where(sig_arr[0,:]==0.))
        fig=plt.figure()
        ax = fig.add_subplot(111)
        L.showGs(fig, ax)
        for i in range(np.shape(zsig)[1]-1):
            print i
            i1 = zsig[:,i][0]+1
            i2 = zsig[:,i+1][0]
            sig = sig_arr[:, i1:i2]
            s = Signature(sig)
            print 'sig'
            print sig
            Mi, Yi = s.sig2ray(L, pTx, pRx)
            print np.shape(Mi)
            print np.shape(Yi)
            if Yi is not None:
                rayi = np.hstack((np.hstack((pTx.reshape(2, 1), Yi)), pRx.reshape(2, 1)))
                ax.plot(rayi[0, :], rayi[1, :], alpha=0.6, linewidth=1.)
            


class Signature(object):
    """ class Signature

    A signature contains two lists

    seq : list of interaction numbers
    typ : list of interaction type
    """
    def __init__(self, sig):
        """
        pa  : tail point of intercation segment
        pb  : head point of intrcation segement
        pc  : middle point  of interaction segment
        typ : type of interaction 1-R 2-T 3-D
        seq : sequence of interaction point (edges (>0)  or vertices (<0)
        """
        self.seq = sig[0, :]
        self.typ = sig[1, :]

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
        #self.typ = np.zeros(N)
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
                #self.typ[n] = 1
            else:      # node
                pa = np.array(L.Gs.pos[k])
                norm = np.array([0, 0])
                self.pa[:, n] = pa
                self.pb[:, n] = pa
                self.pc[:, n] = pa
                self.norm[:, n] = norm
                #self.typ[n] = 3
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

    def evtx(self, L, tx, rx):
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
        self.typ = np.sign(u1 * u2)
        #return(vn)
        #return(typ)

    def image2(self, tx):
        """
        """
        pa = self.pa
        pb = self.pb
        pab = pb - pa
        alpha = np.sum(pab * pab, axis=0)
        zalpha = np.where(alpha==0.)
        alpha[zalpha] = 1.

        a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
        b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
        c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
                            pa[1, :] * (pa[1, :] - pb[1, :]) *
                            (pb[0, :] - pa[0, :]))
        d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
                            pa[0, :] * (pa[1, :] - pb[1, :]) *
                            (pb[0, :] - pa[0, :]))
        
        typ = self.typ
        # number of interactions
        N = np.shape(pa)[1]

        S = np.zeros((N, 2, 2))
        S[:, 0, 0] = -a
        S[:, 0, 1] = b
        S[:, 1, 0] = b
        S[:, 1, 1] = a
        blocks = np.zeros((N-1,2,2))
        A = np.eye(N*2)
        
        # detect diffraction
        usig = np.nonzero(typ[1:] == 3)[0]
        if len(usig)>0:
            blocks[usig,:,:] = np.zeros((2,2))
        # detect transmission
        tsig = np.nonzero(typ[1:] == 2)[0]
        if len(tsig)>0:
            blocks[tsig,:,:] = np.zeros((2,2))
        # detect reflexion
        rsig = np.nonzero(typ[1:] == 1)[0]
        if len(rsig)>0:
            blocks[rsig,:,:] = S[rsig+1,:,:]
        A = pyu.fill_block_diag(A, blocks, 2, -1)

        y = np.zeros(2*N)
        if typ[0] == 1:
            vc0 = np.array([c[0], d[0]])
            v0 = np.dot(-S[0, :, :], tx) + vc0
        if typ[0] == 2:
            v0 = tx
        if typ[0] == 3:
            v0 = pa[:, 0]
        y[0:2] = v0
        for i in range(len(typ[1:])):
            if typ[i+1] == 1:
                y[2*(i+1):2*(i+1)+2] = np.array([c[i+1],d[i+1]])
            if typ[i+1] == 2:
                y[2*(i+1):2*(i+1)+2] = np.zeros(2)
            if typ[i+1] == 3:
                y[2*(i+1):2*(i+1)+2] = pa[:,i+1]


        x = la.solve(A, y)
        M = np.vstack((x[0::2], x[1::2]))
        return M
        

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

        #Ntx = np.shape(tx)
        pab = pb - pa
        #pabo = np.array([pa[1, :] - pb[1, :], pb[0, :] - pa[0, :]])
        alpha = np.sum(pab * pab, axis=0)
        zalpha = np.where(alpha==0.)
        alpha[zalpha] = 1.

        a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
        b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
        c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
                            pa[1, :] * (pa[1, :] - pb[1, :]) *
                            (pb[0, :] - pa[0, :]))
        d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
                            pa[0, :] * (pa[1, :] - pb[1, :]) *
                            (pb[0, :] - pa[0, :]))
        #
        # introduce the diffraction point coords in vk (using c and d)
        #
        c[zalpha] = pa[0,zalpha]
        d[zalpha] = pa[1,zalpha]
        #I2 = np.eye(2)
        #Z2 = np.zeros((2, 2))
        S = np.zeros((N, 2, 2))
        D = np.zeros((N, 2, 2))
        D[:, 0, 0] = 1
        D[:, 1, 1] = 1
        S[:, 0, 0] = -a
        S[:, 0, 1] = b
        S[:, 1, 0] = b
        S[:, 1, 1] = a
        # N : Nombre d'interactions
        # A est une matrice 2*N x 2*N
        #
        A = np.diag(np.kron(dsig, [1, 1]))
        if N > 1:
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
            # dm1  : Sous diagonale 1 : 2*N - 1
            dm1 = np.zeros(2 * N - 1)
            #  0 b 0 b 0 b
            kdm1 = range(1, 2 * N - 1, 2)
            dm1[kdm1] = b1
            # dm2 : Sous diagonale 2 : 2*(N-1)
            # -a a -a a -a a  ...
            dm2 = -np.kron(a1, np.array([1, 0])) \
                  + np.kron(a1, np.array([0, 1]))
            dm2 = dm2 + np.kron(a2, [1, 1])
            # annulation de la sous diganale 2 pour la diffraction
            # if diffraction
            if len(usig) > 0:
                v = np.kron(2 * (usig-1), np.array([1, 0])) +\
                    np.kron(2 * usig -1, np.array([0, 1]))
                dm2[v] = 0
            dm3 = np.zeros(2 * N - 3)
            kdm3 = range(0, 2 * N - 3, 2)
            dm3[kdm3] = b1
            Am1 = np.diagflat(dm1, -1)
            Am2 = np.diagflat(dm2, -2)
            Am3 = np.diagflat(dm3, -3)
            A = A + Am1 + Am2 + Am3
        #
        # Evaluate y
        #
        # depending on the first interaction the 2 first terms of y are
        #   1 - (R) the image of tx wrt first interaction segment  (far)
        #   2 - (T) the transmitter tx itself                    (close)
        #   3 - (D) the diffraction point itself                (middle)
        if typ[0] == 1:
            vc0 = np.array([c[0], d[0]])
            v0 = np.dot(-S[0, :, :], tx) + vc0
        if typ[0] == 2:
            v0 = tx
        if typ[0] == 3:
            v0 = pa[:, 0]
        # for following intercations
        #  1 - (R)    vk
        #  2 - (T)    z2
        #  3 - (D)    the diffarction point itself
        
        if N > 1:
            y = np.hstack((v0, np.kron(c1, np.array([1, 0])) +
                               np.kron(d1, np.array([0, 1]))))
        else:
            y = v0
        x = la.solve(A, y)
        M = np.vstack((x[0::2], x[1::2]))
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
        cpt = 0
        while (((beta <= 1) & (beta >= 0)) & (k < N)):
            if int(typ[k]) != 3:
                # Formula (30) of paper Eucap 2012
                l0 = np.hstack((I2, pkm1 - M[:, -(k + 1)].reshape(2, 1), z0
                              ))
                l1 = np.hstack((I2, z0,
                     pa[:, -(k + 1)].reshape(2, 1) -
                     pb[:, -(k + 1)].reshape(2, 1)
                     ))
                
                T = np.vstack((l0, l1))
                yk = np.hstack((pkm1[:, 0].T, pa[:, -(k + 1)].T))
                deT = np.linalg.det(T)
                if abs(deT)<1e-15:
                    return(None)
                xk = la.solve(T, yk)
                pkm1 = xk[0:2].reshape(2, 1)
                gk = xk[2::]
                beta = gk[1]
                Y = np.hstack((Y, pkm1))
            else:
                Y = np.hstack((Y, pa[:,k].reshape((2,1))))
                pkm1 = pa[:,k].reshape((2,1))
            k = k + 1

        if ((k == N) & ((beta > 0) & (beta < 1))):
            Y = np.hstack((Y, tx.reshape(2, 1)))
            return(Y)
        else:
            return(None)

    def sig2ray(self, L, pTx, pRx):
        """
        from signature to rays
        Parameters
        ----------
            L : Layout
            pTx : ndarray
                2D transmitter position
            pRx : ndarray
                2D receiver position
        Returns
        -------
            M : numpy.ndarray
            Y : numpy.ndarray
        """
        try:
            L.Gr
        except:
            L.build()

        self.ev(L)
        M = self.image2(pTx)
        Y = self.backtrace(pTx, pRx, M)
        return M, Y
        


if __name__ == "__main__":
    doctest.testmod()
