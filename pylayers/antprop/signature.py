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
import pylayers.util.geomutil as geu
#import pylayers.util.graphutil as gph
import pylayers.util.pyutil as pyu
import matplotlib.pyplot as plt
from pylayers.util.project import *
from mpl_toolkits.mplot3d import Axes3D

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

    def get_sigslist(self, tx, rx):
        """
        get signatures (in one list of arrays) between tx and rx
        Parameters
        ----------
            tx : numpy.ndarray
            rx : numpy.ndarray
        Returns
        -------
            sigslist = numpy.ndarray
        """
        try:
            self.L.Gi
        except:
            self.L.build()
        # Here we take all the vnodes >0  from the room
        #
        # Practically those list of nodes should depend on tx , rx
        #
        NroomTx = self.L.pt2ro(tx)
        NroomRx = self.L.pt2ro(rx)

        if not self.L.Gr.has_node(NroomTx) or not self.L.Gr.has_node(NroomRx):
            raise AttributeError('Tx or Rx is not in Gr')
        ndt = self.L.Gt.node[self.L.Gr.node[NroomTx]['cycle']]['inter']
        ndr = self.L.Gt.node[self.L.Gr.node[NroomRx]['cycle']]['inter']
        ntr = np.intersect1d(ndt,ndr)
        sigslist=[]
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
                    sigarr = np.array([]).reshape(2, 0)
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
                    sigslist.append(sigarr)
        
        return sigslist

    def update_sigslist(self):
        """
        get sigantures taking into account reverberations
        Returns
        -------
            sigslist: numpy.ndarry

        Notes
        -----
        This is a prelminary function need more investigations
        """
        pTx = self.pTx
        pRx = self.pRx
        NroomTx = self.L.pt2ro(pTx)
        NroomRx = self.L.pt2ro(pRx)
        if NroomTx == NroomRx:
            sigslist = self.get_sigslist(pTx, pRx)
        else:
            sigslist=[]
            sigtx = self.get_sigslist(pTx, pTx)
            sigrx = self.get_sigslist(pRx, pRx)
            sigtxrx = self.get_sigslist(pTx, pRx)
            sigslist=sigslist+sigtxrx
            for sigtr in sigtxrx:
                for sigt in sigtx:
                    if (sigt[:,-1] == sigtr[:,0]).all():
                        if np.shape(sigtr)[1]==1 or np.shape(sigt)[1]==1:
                            pass
                        else:
                            sigslist.append(np.hstack((sigt,sigtr[:,1:])))
                for sigr in sigrx:
                    if (sigr[:,0] == sigtr[:,-1]).all():
                        if np.shape(sigtr)[1]==1 or np.shape(sigr)[1]==1:
                            pass
                        else:
                            sigslist.append(np.hstack((sigtr,sigr[:,1:]))) 
                 
        return sigslist
            
    def image_ceilfloor(self, tx, pa, pb):
        """
        Compute the images of tx with respect to ceil or floor
        Parameters
        ----------
            tx : numpy.ndarray
            pa : numpy.ndarray
            pb : numpy.ndarray
        Returns
        -------
            M : numpy.ndarray
        """
        
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

        S = np.zeros((1, 2, 2))
        S[:, 0, 0] = -a
        S[:, 0, 1] = b
        S[:, 1, 0] = b
        S[:, 1, 1] = a
        A = np.eye(2)

        vc0 = np.array([c[0], d[0]])
        y = np.dot(-S[0, :, :], tx) + vc0


        x = la.solve(A, y)
        M = np.vstack((x[0::2], x[1::2]))
        return M

    def backtrace_ceilfloor(self, tx, rx, pa, pb, M):
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


        """
        N = np.shape(pa)[1]
        I2 = np.eye(2)
        z0 = np.zeros((2, 1))

        pkm1 = rx.reshape(2, 1)
        Y = pkm1
        k = 0
        beta = .5
        cpt = 0
        while (((beta <= 1) & (beta >= 0)) & (k < N)):
            l0 = np.hstack((I2, pkm1 - M[:, N-(k + 1)].reshape(2, 1), z0
                          ))
            l1 = np.hstack((I2, z0,
                 pa[:, N-(k + 1)].reshape(2, 1) -
                 pb[:, N-(k + 1)].reshape(2, 1)
                 ))
            
            T = np.vstack((l0, l1))
            yk = np.hstack((pkm1[:, 0].T, pa[:, N-(k + 1)].T))
            deT = np.linalg.det(T)
            if abs(deT)<1e-15:
                return(None)
            xk = la.solve(T, yk)
            pkm1 = xk[0:2].reshape(2, 1)
            gk = xk[2::]
            alpha = gk[0]
            beta = gk[1]
            Y = np.hstack((Y, pkm1))
            k+=1
        if ((k == N) & ((beta > 0) & (beta < 1))):# & ((alpha > 0) & (alpha < 1))):
            Y = np.hstack((Y, tx.reshape(2, 1)))
            return(Y)
        else:
            return(None)
        
    def sigs2rays(self, sigslist):
        """
        from signatures list to 2D rays
        Parameters
        ----------
            sigslist : list
        Returns
        -------
            rays : dict
        """
        rays = {}
        for sig in sigslist:
            s = Signature(sig)
            Yi = s.sig2ray(self.L, self.pTx[:2], self.pRx[:2])
            if Yi is not None:
                Yi = np.fliplr(Yi)
                nint = len(sig[0,:])
                if str(nint) in rays.keys():
                    Yi3d = np.vstack((Yi[:,1:-1],np.zeros((1,nint))))
                    Yi3d = Yi3d.reshape(3,nint,1)
                    rays[str(nint)]['pt'] = np.dstack((
                            rays[str(nint)]['pt'],Yi3d))
                    rays[str(nint)]['sig'] = np.dstack((
                            rays[str(nint)]['sig'],
                            sig.reshape(2,nint,1)))
                else:
                    rays[str(nint)] = {'pt': np.zeros((3,nint,1)),
                                       'sig': np.zeros((2,nint,1))}
                    rays[str(nint)]['pt'][0:2,:,0] = Yi[:,1:-1]
                    rays[str(nint)]['sig'][:,:,0] = sig
        return rays

    def show_rays2D(self, rays):
        """
        plot 2D rays within the simulated environment
        Parameters
        ----------
            rays: dict
        """

        fig=plt.figure()
        ax = fig.add_subplot(111)
        self.L.showGs(fig, ax)
        ax.plot(self.pTx[0], self.pTx[1], 'or')
        ax.plot(self.pRx[0], self.pRx[1], 'og')
        for i in rays.keys():
            for j in range(len(rays[i]['pt'][0,0,:])):
                ray = np.hstack((self.pTx[0:2].reshape((2,1)),
                                 np.hstack((rays[i]['pt'][0:2,:,j],
                                          self.pRx[0:2].reshape((2,1))))
                               ))
                ax.plot(ray[0, :], ray[1, :], alpha=0.6, linewidth=1.)


    def ray2D3D(self, rays):
        """
        transform 2D ray to 3D ray (no ceil no floor here)
        Parameters
        ----------
            rays : dict
        
        Returns
        -------
            rays : dict
        """
        pTx = self.pTx
        pRx = self.pRx
        for i in rays.keys():
            pts = rays[i]['pt'][0:2,:,:]
            sig = rays[i]['sig']
            t = self.pTx[0:2].reshape((2,1,1))*\
                np.ones((1,1,len(pts[0,0,:])))
            r = self.pRx[0:2].reshape((2,1,1))*\
                np.ones((1,1,len(pts[0,0,:])))
            pts1 = np.hstack((t,np.hstack((pts,r))))
            si1 = pts1[:,1:,:]-pts1[:,:-1,:]
            si = np.sqrt(np.sum(si1*si1,axis=0))
            alpha = np.zeros(np.shape(si[:-1,:]))
            for j in range(len(alpha[:,0])):
                alpha[j,:] = np.sum(si[0:j+1,:],axis=0)/\
                             np.sum(si,axis=0)
                rays[i]['pt'][2,j,:]=  pTx[2]+alpha[j,:]*(pRx[2]-pTx[2])
                
        rays = self.ray_ceilfloor(rays=rays, nr=1)
        return rays


    def ray_ceilfloor(self, rays, nr=1):
        """
        compute 3D rays reflected nr times on ceil and floor
        Parameters
        ----------
            rays : dict
            nr : int
        
        Returns
        -------
            rays : dict
        """
        #
        # Compute for floor
        #
        pax = np.array([[0.],[3.]])
        pbx = np.array([[10.],[3.]])
        pay = np.array([[-2.],[3.]])
        pby = np.array([[2.],[3.]])

        txx = np.array([self.pTx[0],self.pTx[2]])
        txy = self.pTx[1:]

        rxx = np.array([self.pRx[0],self.pRx[2]])
        rxy = self.pRx[1:]
        
        Mx = self.image_ceilfloor(txx, pax, pbx)
        My = self.image_ceilfloor(txy, pay, pby)

        Yx = self.backtrace_ceilfloor(txx, rxx, pax, pbx, Mx)
        Yy = self.backtrace_ceilfloor(txy, rxy, pax, pbx, My)

        Yxy = np.vstack((Yx[0:1,:],Yy[0:1,:]))
        pts = np.array([]).reshape(3,0)
        sig = np.array([]).reshape(2,0)
        Ii = []
        for i in range(len(Yxy[0,:])-1):
            p1 = Yxy[:,i]
            p2 = Yxy[:,i+1]
            I = self.L.seginline(p1,p2)
            #I = np.hstack((I,it))
            if np.shape(I)[1] != 0:
                print np.shape(I)
                Iz = np.nan*np.ones((1,np.shape(I)[1]))
                I = np.vstack((I[1:,:],Iz))
                print I 
                pts = np.hstack((pts, I))
                pts = np.hstack((pts, np.vstack((Yx[0:1,i],Yy[:,i]))))
        print pts
        rayf = np.vstack((Yx[0:1,1:-1],Yy[:,1:-1]))
        #print rayf
        rays[str(nr)]['pt'] = np.dstack((rays[str(nr)]['pt'],rayf))

        return rays

        
        

        

    def show_ray3d(self, _filestr='defstr', ray = np.array([]), bdis=True
            , bbas=False, bstruc=True, col=np.array([1, 0, 1]), id=0,
            linewidth=1):
        """
        plot a 3D ray
        Parameters
        ----------

        bdis :
            display boolean - if False return .vect filename
        bbas :
            display local basis
        bstruc :
            display structure
        col  :
            color of the ray
        id   :
            id of the ray
        linewidth :
        """
        
        filerac = pyu.getlong("ray" + str(id), pstruc['DIRGEOM'])
        _filerac = pyu.getshort(filerac)
        filename_list = filerac + '.list'
        filename_vect = filerac + '.vect'
        try:
            fo = open(filename_vect, "w")
        except:
            raise NameError(filename)

        fo.write("appearance { linewidth %d }\n" % linewidth)

        fo.write("VECT\n")

        fo.write("1 %d 1\n\n" % len(ray[0,:]))
        fo.write("%d\n" % len(ray[0,:]))
        fo.write("1\n")
        for i in range(len(ray[0,:])):
            fo.write("%g %g %g\n" % (ray[0, i], ray[1,i],
                                    ray[2, i]))
        #fo.write("%d %d %d 0\n" % (col[0],col[1],col[2]))
        fo.write("%g %g %g 0\n" % (col[0], col[1], col[2]))
        fo.close()

        #
        # Ajout des bases locales
        #

        fo = open(filename_list, "w")
        fo.write("LIST\n")
        fo.write("{<" + filename_vect + "}\n")
        if (bstruc):
            #fo.write("{<strucTxRx.off}\n")
            fo.write("{<" +_filestr +".off}\n")

        filename = filename_list
        fo.close()

        if (bdis):
        #
        # Geomview Visualisation
        #
            chaine = "geomview -nopanel -b 1 1 1 " + filename + \
                " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)

        
    def show3(self, rays={}, bdis=True, bstruc=True, id=0):
        """
        plot 3D rays within the simulated environment
        Parameters
        ----------
            raysarr: numpy.ndarray
        """
        pTx = self.pTx.reshape((3,1))
        pRx = self.pRx.reshape((3,1))
        filename = pyu.getlong("grRay" + str(id) + ".list", pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        if bstruc:
            fo.write("{<defstr.off}\n")
            #fo.write("{<strucTxRx.off}\n")
            k = 0
            for i in rays.keys():
                for j in range(np.shape(rays[i]['pt'])[2]):
                    ray = np.hstack((pTx,
                            np.hstack((rays[i]['pt'][:,:,j],pRx))))
                    #ray = rays[i]['pt'][:,:,j]
                    col = np.array([2, 0, 1])
                    fileray = self.show_ray3d(ray=ray, bdis=False, bstruc=False, col=col, id=k)
                    k+=1
                    fo.write("{< " + fileray + " }\n")
        fo.close()
        if (bdis):
            chaine = "geomview " + filename + " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)
        
        

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
                l0 = np.hstack((I2, pkm1 - M[:, N-(k + 1)].reshape(2, 1), z0
                              ))
                l1 = np.hstack((I2, z0,
                     pa[:, N-(k + 1)].reshape(2, 1) -
                     pb[:, N-(k + 1)].reshape(2, 1)
                     ))
                
                T = np.vstack((l0, l1))
                yk = np.hstack((pkm1[:, 0].T, pa[:, N-(k + 1)].T))
                deT = np.linalg.det(T)
                if abs(deT)<1e-15:
                    return(None)
                xk = la.solve(T, yk)
                pkm1 = xk[0:2].reshape(2, 1)
                gk = xk[2::]
                alpha = gk[0]
                beta = gk[1]
                Y = np.hstack((Y, pkm1))
            else:
                Y = np.hstack((Y, pa[:,k].reshape((2,1))))
                pkm1 = pa[:,k].reshape((2,1))
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

if __name__ == "__main__":
    doctest.testmod()
