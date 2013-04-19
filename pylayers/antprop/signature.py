#-*- coding:Utf-8 -*-
import doctest
import numpy as np
#import scipy as sp
import scipy.linalg as la
import pdb
import networkx as nx
#import GrRay3D
#import Graph
#import pylayers.gis.layout
import pylayers.util.geomutil as geu
#import pylayers.util.graphutil as gph
import pylayers.util.pyutil as pyu
import matplotlib.pyplot as plt
from pylayers.util.project import *
from mpl_toolkits.mplot3d import Axes3D
from pylayers.antprop.rays import Rays
#from numba import autojit

def showsig(L,s,tx,rx):
    """
    """
    L.display['thin']=True
    fig,ax = L.showGs()
    L.display['thin']=False
    L.display['edlabel']=True
    L.showGs(fig=fig,ax=ax,edlist=s,width=4)
    plt.plot(tx[0],tx[1],'x')
    plt.plot(rx[0],rx[1],'+')
    plt.title(str(s))
    plt.show()
    L.display['edlabel']=False

class Signatures(dict):
    """
    gathers all signatures from a layout given tx and rx

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
        Parameters
        ----------
        L : Layout
        pTx :
        pRx : 
        """
        self.L = L
        self.pTx = pTx
        self.pRx = pRx

#    def __repr__(self):
#        s = 
#        return(s)
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


    def all_simple_paths(self,G, source, target, cutoff=None):
        if cutoff < 1:
            return

        visited = [source]
        stack = [iter(G[source])]
        while stack:
            # stack is a list of iterators

            children = stack[-1]
            child = next(children, None)
            if child is None:
                stack.pop()
                visited.pop()
            elif len(visited) < cutoff:
                if child == target:
                    yield visited + [target]
                elif child not in visited:
                    visited.append(child)
                    # explore all child connexion
                    # TODO : limit the explorable childs
                    stack.append(iter(G[child]))
            else: #len(visited) == cutoff:
                if child == target or target in children:
                    yield visited + [target]
                stack.pop()
                visited.pop()


    def run(self, tx, rx,cutoff=1):
        """ get signatures (in one list of arrays) between tx and rx

        Parameters
        ----------

            tx : numpy.ndarray
            rx : numpy.ndarray
            cutoff :

        Returns
        -------

            sigslist = numpy.ndarray

        """
        try:
            self.L.Gi
        except:
            self.L.build()
        # all the vnodes >0  from the room
        #
        NroomTx = self.L.pt2ro(tx)
        NroomRx = self.L.pt2ro(rx)
        #print NroomTx,NroomRx

        if not self.L.Gr.has_node(NroomTx) or not self.L.Gr.has_node(NroomRx):
            raise AttributeError('Tx or Rx is not in Gr')

        # list of interaction in roomTx 
        # list of interaction in roomRx
        ndt = self.L.Gt.node[self.L.Gr.node[NroomTx]['cycle']]['inter']
        ndr = self.L.Gt.node[self.L.Gr.node[NroomRx]['cycle']]['inter']

        # transmitter
        ndt1 = filter(lambda l: len(eval(l))>2,ndt) # Transmission
        ndt2 = filter(lambda l: len(eval(l))<3,ndt) # Reflexion

        # receiver
        ndr1 = filter(lambda l: len(eval(l))>2,ndr) # Transmission
        ndr2 = filter(lambda l: len(eval(l))<3,ndr) # Reflexion

        # tx,rx : attaching rule
        #
        # tx attachs to out transmisision point 
        # rx attachs to in transmission point

        #
        # WARNING : room number <> cycle number
        #

        ncytx = self.L.Gr.node[NroomTx]['cycle']
        ncyrx = self.L.Gr.node[NroomRx]['cycle']

        ndt1 = filter(lambda l: eval(l)[2]<>ncytx,ndt1)
        ndr1 = filter(lambda l: eval(l)[1]<>ncyrx,ndr1)


        ndt = ndt1 + ndt2
        ndr = ndr1 + ndr2
        #print ndt
        #print ndr
        ntr = np.intersect1d(ndt, ndr)

        for nt in ndt:
            for nr in ndr:

                if (nt != nr):
                    paths = list(self.all_simple_paths(self.L.Gi,source=nt,target=nr,cutoff=cutoff))
                else:
                    paths = [[nt]]
                ### supress the followinfg loops .
                for path in paths:
                    sigarr = np.array([],dtype=int).reshape(2, 0)
                    for interaction in path:

                        it = eval(interaction)
                        if type(it) == tuple:
                            if len(it)==2: #reflexion
                                sigarr = np.hstack((sigarr,
                                                np.array([[it[0]],[1]],dtype=int)))
                            if len(it)==3: #transmission
                                sigarr = np.hstack((sigarr,
                                                np.array([[it[0]],[2]],dtype=int)))
                        elif it < 0: #diffraction
                            sigarr = np.hstack((sigarr,
                                                np.array([[it],[3]],dtype=int)))
                    #print sigarr
                    try:
                        self[len(path)] = np.vstack((self[len(path)],sigarr))
                    except:
                        self[len(path)] = sigarr


    def showi(self):
        """ interactive show
        
        press n to visit rays

        Attributes
        ----------
            ni : number of interaction
            us : signature index
        """
        plt.ion()
        fig=plt.figure()
        ax=fig.add_subplot(111)
#        fig,ax=self.L.showG(fig=fig,ax=ax,graph='s')
#        plt.draw()

        nit = self.keys()
        uni = 0
        ni = nit[uni]
        
        ust = len(self[ni])/2
        us = 0

        st='a'
        while st != 'q':
            inter=[]
            ax=fig.add_subplot(111)
            fig,ax=self.L.showG(fig=fig,ax=ax,graph='s')
            title = '# interaction :', ni, 'signature #',us,'/',ust
             
            ax.set_title(title)

            line = self.pTx[:2]
            ax.plot(self.pRx[0],self.pRx[1],'xb')
            ax.plot(self.pTx[0],self.pTx[1],'xr')

            if ni not in self.keys():
                print "incorrect number of interactions"
                
            pos={}

            try:
                for u in self[ni][us*2]:
                    pos.update({u:self.L.Gs.pos[u]})
                    line = np.vstack((line,np.array((self.L.Gs.pos[u]))))
                nx.draw_networkx_nodes(self.L.Gs,pos=pos,nodelist=pos.keys(),node_color='r',ax=ax)

                for ii in self[ni][(us*2)+1]:
                    if ii == 1:
                        inter.append('R')
                    if ii == 2:
                        inter.append('T')




            except:
                print "signature index out of bounds of signature"
            line = np.vstack((line,self.pRx[:2]))
            ax.plot(line[:,0],line[:,1])
            
            plt.draw()
            print inter
            st = raw_input()
            ax.cla()
            if st == 'n':
                if us+2 <= ust:
                    us=us+2

                else:
                    uni = uni+1
                    try:
                        ni = nit[uni]
                        ust = len(self[ni])/2
                        us=0
                    except:
                        uni=0
                        ni=nit[uni]
                        us = 0
                


            else:
                print 'press n for next signature'    
    


    def rays(self):
        """ from signatures dict to 2D rays

        Parameters
        ----------

            dsig : dict 

        Returns
        -------

            rays : dict

        """
        rays = Rays(self.pTx,self.pRx)
        for k in self:
            tsig = self[k]
            shsig = np.shape(tsig)
            for l in range(shsig[0]/2):
                sig = tsig[2*l:2*l+2,:]
                s   = Signature(sig)
                Yi  = s.sig2ray(self.L, self.pTx[:2], self.pRx[:2])
                if Yi is not None:
                    Yi = np.fliplr(Yi)
                    nint = len(sig[0, :])
                    if nint in rays.keys():
                        Yi3d = np.vstack((Yi[:, 1:-1], np.zeros((1, nint))))
                        Yi3d = Yi3d.reshape(3, nint, 1)
                        rays[nint]['pt'] = np.dstack(( rays[nint]['pt'], Yi3d))
                        rays[nint]['sig'] = np.dstack(( rays[nint]['sig'], sig.reshape(2, nint, 1)))
                    else:
                        rays[nint] = {'pt': np.zeros((3, nint, 1)),
                                      'sig': np.zeros((2, nint,
                                                            1),dtype=int)}
                        rays[nint]['pt'][0:2, :, 0] = Yi[:, 1:-1]
                        rays[nint]['sig'][:, :, 0] = sig
        return rays

class Signature(object):
    """ class Signature

    A signature contains two lists

    seq : list of interaction numbers
    typ : list of interaction type
    """
    def __init__(self, sig):
        """
        pa  : tail point of interaction segment
        pb  : head point of interaction segment
        pc  : center point of interaction segment
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


    def image(self, tx):
        """
        Compute the images of tx with respect to the signature segments
        Parameters
        ----------
            tx : numpy.ndarray
        Returns
        -------
            M : numpy.ndarray
        """

        pa = self.pa
        pb = self.pb
        pab = pb - pa
        alpha = np.sum(pab * pab, axis=0)
        zalpha = np.where(alpha == 0.)
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
        blocks = np.zeros((N - 1, 2, 2))
        A = np.eye(N * 2)

        # detect diffraction
        usig = np.nonzero(typ[1:] == 3)[0]
        if len(usig) > 0:
            blocks[usig, :, :] = np.zeros((2, 2))
        # detect transmission
        tsig = np.nonzero(typ[1:] == 2)[0]
        if len(tsig) > 0:
            #blocks[tsig, :, :] = np.zeros((2, 2))
            blocks[tsig, :, :] = -np.eye(2)
        # detect reflexion
        rsig = np.nonzero(typ[1:] == 1)[0]
        if len(rsig) > 0:
            blocks[rsig, :, :] = S[rsig + 1, :, :]
        A = pyu.fill_block_diag(A, blocks, 2, -1)

        y = np.zeros(2 * N)
        if typ[0] == 1:
            vc0 = np.array([c[0], d[0]])
            v0 = np.dot(-S[0, :, :], tx) + vc0
        if typ[0] == 2:
            v0 = tx
        if typ[0] == 3:
            v0 = pa[:, 0]
        y[0:2] = v0
        for i in range(len(typ[1:])):
            if typ[i + 1] == 1:
                y[2 * (i + 1):2 * (i + 1) + 2] = np.array([c[i + 1], d[i + 1]])
            if typ[i + 1] == 2:
                #y[2 * (i + 1):2 * (i + 1) + 2] = y[2*i:2*i+2]
                y[2 * (i + 1):2 * (i + 1) + 2] = np.array([0,0]) 
            if typ[i + 1] == 3:
                y[2 * (i + 1):2 * (i + 1) + 2] = pa[:, i + 1]

        x = la.solve(A, y)
        M = np.vstack((x[0::2], x[1::2]))
        return M
    def backtrace(self, tx, rx, M):
        """
        backtracing step: given the image, tx, and rx, this function
        traces the 2D ray.

        Parameters
        ----------
            tx :  numpy.ndarray
                  transmitter
            rx :  numpy.ndarray
                  receiver
            M  :  numpy.ndarray
                  images obtained using image()

        Returns
        -------
            Y : numpy.ndarray
                2D ray

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
                l0 = np.hstack((I2, pkm1 - M[:, N - (k + 1)].reshape(2, 1), z0
                                ))
                l1 = np.hstack((I2, z0,
                                pa[:, N - (k + 1)].reshape(2, 1) -
                                pb[:, N - (k + 1)].reshape(2, 1)
                                ))

                T = np.vstack((l0, l1))
                yk = np.hstack((pkm1[:, 0].T, pa[:, N - (k + 1)].T))
                deT = np.linalg.det(T)
                if abs(deT) < 1e-15:
                    return(None)
                xk = la.solve(T, yk)
                pkm1 = xk[0:2].reshape(2, 1)
                gk = xk[2::]
                alpha = gk[0]
                beta = gk[1]
                Y = np.hstack((Y, pkm1))
            else:
                Y = np.hstack((Y, pa[:, k].reshape((2, 1))))
                pkm1 = pa[:, k].reshape((2, 1))
            k = k + 1
        if ((k == N) & ((beta > 0) & (beta < 1)) & ((alpha > 0) & (alpha < 1))):
            Y = np.hstack((Y, tx.reshape(2, 1)))
            return(Y)
        else:
            return(None)

    def sig2ray(self, L, pTx, pRx):
        """
        convert a signature to a 2D ray
        Parameters
        ----------
            L : Layout
            pTx : ndarray
                2D transmitter position
            pRx : ndarray
                2D receiver position
        Returns
        -------
            Y : numpy.ndarray
        """
        try:
            L.Gr
        except:
            L.build()

        self.ev(L)
        M = self.image(pTx)
        Y = self.backtrace(pTx, pRx, M)
        return Y
 
# def get_sigslist(self, tx, rx):
#        """
#        get signatures (in one list of arrays) between tx and rx
#        Parameters
#        ----------
#            tx : numpy.ndarray
#            rx : numpy.ndarray
#        Returns
#        -------
#            sigslist = numpy.ndarray
#        """
#        try:
#            self.L.Gi
#        except:
#            self.L.build()
#        # all the vnodes >0  from the room
#        #
#        NroomTx = self.L.pt2ro(tx)
#        NroomRx = self.L.pt2ro(rx)
#        print NroomTx,NroomRx
#
#        if not self.L.Gr.has_node(NroomTx) or not self.L.Gr.has_node(NroomRx):
#            raise AttributeError('Tx or Rx is not in Gr')
#
#        #list of interaction 
#        ndt = self.L.Gt.node[self.L.Gr.node[NroomTx]['cycle']]['inter']
#        ndr = self.L.Gt.node[self.L.Gr.node[NroomRx]['cycle']]['inter']
#
#        ndt1 = filter(lambda l: len(eval(l))>2,ndt)
#        ndt2 = filter(lambda l: len(eval(l))<3,ndt)
#        ndr1 = filter(lambda l: len(eval(l))>2,ndr)
#        ndr2 = filter(lambda l: len(eval(l))<3,ndr)
#
#        print ndt1
#        print ndr1
#        ndt1 = filter(lambda l: eval(l)[2]<>NroomTx,ndt1)
#        ndr1 = filter(lambda l: eval(l)[1]<>NroomRx,ndr1)
#
#        ndt = ndt1 + ndt2
#        ndr = ndr1 + ndr2
#
#        ntr = np.intersect1d(ndt, ndr)
#        sigslist = []
#
#        for nt in ndt:
#            print nt
#            for nr in ndr:
#                addpath = False
#                print nr
#                if (nt != nr):
#                    try:
#                        path = nx.dijkstra_path(self.L.Gi, nt, nr)
#                        #paths = nx.all_simple_paths(self.L.Gi,source=nt,target=nr)
#                        addpath = True
#                        showsig(self.L,path,tx,rx)
#                    except:
#                        pass
#                if addpath:
#                    sigarr = np.array([]).reshape(2, 0)
#                    for interaction in path:
#                        it = eval(interaction)
#                        if type(it) == tuple:
#                            if len(it)==2: #reflexion
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[1]])))
#                            if len(it)==3: #transmission
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]], [2]])))
#                        elif it < 0: #diffraction
#                            sigarr = np.hstack((sigarr,
#                                                np.array([[it], [3]])))
#                    sigslist.append(sigarr)
#
#        return sigslist
#
#    def update_sigslist(self):
#        """
#        get signatures taking into account reverberations
#
#        Returns
#        -------
#            sigslist: numpy.ndarry
#
#        Notes
#        -----
#        This is a preliminary function need more investigations
#
#        """
#        pTx = self.pTx
#        pRx = self.pRx
#        NroomTx = self.L.pt2ro(pTx)
#        NroomRx = self.L.pt2ro(pRx)
#        if NroomTx == NroomRx:
#            sigslist = self.get_sigslist(pTx, pRx)
#        else:
#            sigslist = []
#            sigtx = self.get_sigslist(pTx, pTx)
#            sigrx = self.get_sigslist(pRx, pRx)
#            sigtxrx = self.get_sigslist(pTx, pRx)
#            sigslist = sigslist + sigtxrx
#            for sigtr in sigtxrx:
#                for sigt in sigtx:
#                    if (sigt[:, -1] == sigtr[:, 0]).all():
#                        if np.shape(sigtr)[1] == 1 or np.shape(sigt)[1] == 1:
#                            pass
#                        else:
#                            sigslist.append(np.hstack((sigt, sigtr[:, 1:])))
#                for sigr in sigrx:
#                    if (sigr[:, 0] == sigtr[:, -1]).all():
#                        if np.shape(sigtr)[1] == 1 or np.shape(sigr)[1] == 1:
#                            pass
#                        else:
#                            sigslist.append(np.hstack((sigtr, sigr[:, 1:])))
#
#        return sigslist
#
#    def image_ceilfloor(self, tx, pa, pb):
#        """
#        Compute the images of tx with respect to ceil or floor
#        Parameters
#        ----------
#            tx : numpy.ndarray
#            pa : numpy.ndarray
#            pb : numpy.ndarray
#        Returns
#        -------
#            M : numpy.ndarray
#        """
#
#        pab = pb - pa
#        alpha = np.sum(pab * pab, axis=0)
#        zalpha = np.where(alpha == 0.)
#        alpha[zalpha] = 1.
#
#        a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
#        b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
#        c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
#                            pa[1, :] * (pa[1, :] - pb[1, :]) *
#                            (pb[0, :] - pa[0, :]))
#        d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
#                            pa[0, :] * (pa[1, :] - pb[1, :]) *
#                            (pb[0, :] - pa[0, :]))
#
#        S = np.zeros((1, 2, 2))
#        S[:, 0, 0] = -a
#        S[:, 0, 1] = b
#        S[:, 1, 0] = b
#        S[:, 1, 1] = a
#        A = np.eye(2)
#
#        vc0 = np.array([c[0], d[0]])
#        y = np.dot(-S[0, :, :], tx) + vc0
#
#        x = la.solve(A, y)
#        M = np.vstack((x[0::2], x[1::2]))
#        return M
#
#    def backtrace_ceilfloor(self, tx, rx, pa, pb, M):
#        """
#        backtracing step: given the image, tx, and rx, this function
#        traces the 2D ray.
#
#        Parameters
#        ----------
#            tx :  numpy.ndarray
#                  transmitter
#            rx :  numpy.ndarray
#                  receiver
#            M  :  numpy.ndarray
#                  images obtained using image()
#
#        Returns
#        -------
#            Y : numpy.ndarray
#                2D ray
#
#
#        """
#        N = np.shape(pa)[1]
#        I2 = np.eye(2)
#        z0 = np.zeros((2, 1))
#
#        pkm1 = rx.reshape(2, 1)
#        Y = pkm1
#        k = 0
#        beta = .5
#        cpt = 0
#        while (((beta <= 1) & (beta >= 0)) & (k < N)):
#            l0 = np.hstack((I2, pkm1 - M[:, N - (k + 1)].reshape(2, 1), z0
#                            ))
#            l1 = np.hstack((I2, z0,
#                            pa[:, N - (k + 1)].reshape(2, 1) -
#                            pb[:, N - (k + 1)].reshape(2, 1)
#                            ))
#
#            T = np.vstack((l0, l1))
#            yk = np.hstack((pkm1[:, 0].T, pa[:, N - (k + 1)].T))
#            deT = np.linalg.det(T)
#            if abs(deT) < 1e-15:
#                return(None)
#            xk = la.solve(T, yk)
#            pkm1 = xk[0:2].reshape(2, 1)
#            gk = xk[2::]
#            alpha = gk[0]
#            beta = gk[1]
#            Y = np.hstack((Y, pkm1))
#            k += 1
#        if ((k == N) & ((beta > 0) & (beta < 1))):  # & ((alpha > 0) & (alpha < 1))):
#            Y = np.hstack((Y, tx.reshape(2, 1)))
#            return(Y)
#        else:
#            return(None)
#   def sigs2rays(self, sigslist):
#        """ from signatures list to 2D rays
#
#        Parameters
#        ----------
#
#            sigslist : list
#
#        Returns
#        -------
#
#            rays : dict
#
#        """
#        rays = {}
#        for sig in sigslist:
#            s = Signature(sig)
#            Yi = s.sig2ray(self.L, self.pTx[:2], self.pRx[:2])
#            if Yi is not None:
#                #pdb.set_trace()
#                Yi = np.fliplr(Yi)
#                nint = len(sig[0, :])
#                if str(nint) in rays.keys():
#                    Yi3d = np.vstack((Yi[:, 1:-1], np.zeros((1, nint))))
#                    Yi3d = Yi3d.reshape(3, nint, 1)
#                    rays[str(nint)]['pt'] = np.dstack((
#                                                      rays[str(nint)]['pt'], Yi3d))
#                    rays[str(nint)]['sig'] = np.dstack((
#                                                       rays[str(nint)]['sig'],
#                                                       sig.reshape(2, nint, 1)))
#                else:
#                    rays[str(nint)] = {'pt': np.zeros((3, nint, 1)),
#                                       'sig': np.zeros((2, nint, 1))}
#                    rays[str(nint)]['pt'][0:2, :, 0] = Yi[:, 1:-1]
#                    rays[str(nint)]['sig'][:, :, 0] = sig
#        return rays


if __name__ == "__main__":
    doctest.testmod()
