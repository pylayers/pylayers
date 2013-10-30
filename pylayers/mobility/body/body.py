import numpy as np
import scipy.stats as sp

from pylayers.mobility.body import c3d
from pylayers.mobility import trajectory 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import pdb as pdb
import pylayers.util.pyutil as pyu
import pylayers.util.plotutil as pltu
import pylayers.util.geomutil as geu
import doctest

def ChangeBasis(u0, v0, w0, v1):
    """

    Parameters
    ----------

    u0
    v0
    w0
    v1

    """

    # Rotate with respect to axe w

    v2 = v1 - np.dot(np.dot(v1, w0), w0)  # projection of v1 on plan (u,v)
    v2 = v2 / np.linalg.norm(v2)
    c = np.dot(v2, u0)
    s = np.dot(v2, v0)
    u1 = np.dot(s, u0) - np.dot(c, v0)
    u1 = u1 / np.linalg.norm(u1)
    w1 = np.cross(u1, v1)
    return u1, v1, w1


def dist(A, B):
    """
    evaluate the distance between two points A and B
    """

    d = np.sqrt((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    return d


class BodyCylinder(object):
    """ Class to manage the Body model 
    """

    def __init__(self):
        self.g = nx.Graph()
        self.npoints = 15 
        self.nodes_Id = {
            0:'STRN', 
            1:'CLAV',
            2:'RFHD', 
            3:'RSHO', 
            4:'LSHO',
            5:'RELB', 
            6:'LELB', 
            7:'RWRB', 
            8:'LWRB',
            9:'RFWT',
            10:'LFWT', 
            11:'RKNE',
            12:'LKNE', 
            13:'RANK', 
            14:'LANK'}
#            15:'LFIN'}
#            16:'LELB',
#            17:'RUPA',
#            18:'LUPA',
#            19:'RFRM',
#            20:'LFRM',
#            21:'LSHN',
#            22:'RSHN',
#            23:'LTHI',
#            24:'RTHI',
#            25:'LHEE',
#            26:'RHEE',
#            27:'LTOE',
#            28:'RTOE',
#            29:'RMT5'}
        self.g.add_nodes_from(self.nodes_Id.keys())
        self.g.add_edge(0, 1,radius =1)
        self.g[0][1]['radius']=0.18
        #self.g.add_edge(0, 9)
        #self.g.add_edge(0, 10)
        self.g.add_edge(1, 2)
        self.g[1][2]['radius']=0.12
        #self.g.add_edge(1, 3)
        #self.g.add_edge(1, 4)
        self.g.add_edge(3, 5)
        self.g[3][5]['radius']=0.05
        self.g.add_edge(4, 6)
        self.g[4][6]['radius']=0.05
        self.g.add_edge(5, 7)
        self.g[5][7]['radius']=0.05
        self.g.add_edge(6, 8)
        self.g[6][8]['radius']=0.05
        self.g.add_edge(9, 11)
        self.g[9][11]['radius']=0.05
        self.g.add_edge(10, 12)
        self.g[10][12]['radius']=0.05
        self.g.add_edge(11, 13)
        self.g[11][13]['radius']=0.05
        self.g.add_edge(12, 14)
        self.g[12][14]['radius']=0.05

    def __repr__(self):
        st = ''
        
        if 'filename' in dir(self):
            st = st +'filename : '+ self.filename +'\n'
        if 'nframes' in dir(self):    
            st = st +'nframes : ' + str(self.nframes) +'\n'
        if 'pg' in dir(self):
            st = st + 'Centered : True'+'\n'
        if 'topos' in dir(self):
            st = st + 'topos : True'+'\n'
        else:    
            st = st + 'topos : False'+'\n'
        return(st)    

    def center(self):
        """ centering the body 

        Returns
        -------

        self.pg : center of gravity 
        self.d  : set of centered frames


        Notes
        -----

        Here only the projection of the body centroid in the
        plane 0xy is calculated

        """
        # self.d  : 3 x 16 x Nf
        # self.pg : 3 x Nf
        self.pg = np.sum(self.d,axis=1)/self.npoints
        self.pg[2,:] = 0
        self.d = self.d - self.pg[:,np.newaxis,:]


    def posvel(self,traj,tk,Tstep):
        """ position and velocity

        Parameters
        ----------

        traj : ndarray
            nx3
        tk : float 
            time for evaluation of topos
        Tstep : flloat 
        
        Returns
        -------
        
        kf 
        kt 
        vsn
        wsn 
        vtn 
        wtn 

        """
        # tk should be in the trajectory time range
        assert ((tk<traj[-1,0]) & (tk>traj[0,0])),'posvel: tk not in trajectory time range'

        tf = Tstep/(1.0*self.nframes) # frame sampling period  
        kt = int(np.floor(tk/tf))
        kf = int(np.floor(np.mod(tk,Tstep)/tf))
        # self.pg : 3 x Nframes 
        # traj : Nptraj x 3 (t,x,y)

        # vs  : speed vector along motion capture frame
        # vsn : unitary speed vector along motion capture frame
        vs = self.pg[0:-1,kf] - self.pg[0:-1,kf-1]
        vsn = vs/np.sqrt(np.dot(vs,vs))
        wsn = np.array([vsn[1],-vsn[0]])
        #
        # vt : speed vector along trajectory 
        #
        vt = traj[kt+1,1:] - traj[kt,1:]
        vtn = vt/np.sqrt(np.dot(vt,vt))
        wtn = np.array([vtn[1],-vtn[0]])
        
        return(kf,kt,vsn,wsn,vtn,wtn)

    def settopos(self,traj,tk=0,Tstep=3):
        """ translate the body on a time stamped trajectory

        Parameters
        ----------

        traj : ndarray (3,N)
            t,x,y
        tk : float 
            time for evaluation of topos (seconds) this value should be in the
            range of the trajectory timestamp
        Tstep : float 
           duration of the periodic motion sequence (seconds)

        Examples
        --------

        >>> import numpy as np 
        >>> time = np.arange(0,10,0.1)
        >>> v = 4000/3600.
        >>> x = v*time
        >>> y = np.zeros(len(time))
        >>> traj = np.vstack((time.T,x.T,y.T))
        >>> bc = BodyCylinder()
        >>> bc.loadC3D(filename='07_01.c3d')
        >>> bc.center()
        >>> bc.settopos(traj,2.3,2)

        Notes
        -----

        topos is the current spatial global position of a body configuration.


        See Also 
        --------

        pylayers.util.geomutil.affine

        """

        #
        #
        # psa : origin source
        # psb = psa+vsn : a point in the direction of pedestrian motion 
        #
        # pta : target translation 
        # ptb = pta+vtn : a point in the direction of trajectory 
        
        kf,kt,vsn,wsn,vtn,wtn = self.posvel(traj,tk,Tstep)

        psa = np.array([0,0])
        psb = psa + vsn
        psc = psa + wsn

        pta = traj[kt,1:]
        ptb = pta + vtn
        ptc = pta + wtn

        X   = np.array([[0,0],[psb[0],psb[1]],[psc[0],psc[1]]]).T
        Y   = np.array([[pta[0],pta[1]],[ptb[0],ptb[1]],[ptc[0],ptc[1]]]).T

        a,b = geu.affine(X,Y)
        A = np.eye(3)                
        B = np.zeros((3,1))                
        A[0:-1,0:-1] = a                
        B[0:-1,:] = b

        self.topos = (np.dot(A,self.d[:,:,kf])+B)
        
    
    def loadC3D(self, filename='07_01.c3d', nframes=126 ,unit='cm'):
        """ load nframes of motion capture C3D file 

        Parameters
        ----------

        filename : string
            file name 
        nframes  : int 
            number of frames 

        """

        if 'pg' in dir(self):
            del self.pg

        self.nframes = nframes
        self.filename = filename

        s, p, f = c3d.read_c3d(filename)

        CM_TO_M = 0.01
        #
        # self.d 3 x np x nf
        # 
        self.d = np.ndarray(shape=(3, self.npoints, np.shape(f)[0]))
        #if self.d[2,:,:].max()>50:
        ind = []
        for i in self.nodes_Id:
            ind.append(p.index(s[0] + self.nodes_Id[i]))

        # f.T : 3 x np x nf
        self.d = f[0:nframes, ind, :].T
        if unit=='cm':
            self.d = self.d*CM_TO_M
        self.g.pos = {}
        for i in range(self.npoints):
            self.g.pos[i] = (self.d[1, i, 0], self.d[2, i, 0])
        #
        # Extension of cylinder
        #

        self.g.add_node(15)
        self.npoints = 16

        pm  = (self.d[:, 9, 0] + self.d[:, 10, 0])/2.
        pmf = (self.d[:, 9, :] + self.d[:, 10, :])/2.
        pmf = pmf[:,np.newaxis,:]

        self.d = np.concatenate((self.d,pmf),axis=1)

        self.g.add_edge(0, 15)
        self.g.pos[15] = (pm[1],pm[2])
        self.g[0][15]['radius']=0.1

    def movie(self,topos=False,tk=[],traj=[]):
        """ Create a geomview movie

        Parameters
        ----------

        topos : Boolean
        tk    : np.array time index in s 
        traj  : np.array Npt x 3 (t,x,y)    

        See Also
        --------

        BodyCylinder.geomfile

        """

        if not topos:
            for k in range(self.nframes):
                self.geomfile(iframe=k,verbose=True)
        else:
            for k,ttk in enumerate(tk):
                stk = str(k).zfill(6) # for string alignement 
                self.settopos(traj=traj,tk=ttk,Tstep=1)
                self.geomfile(topos=True,verbose=False,tag=stk)

    def plot3d(self,iframe=0,topos=False,fig=[],ax=[],col='b'):
        """ scatter 3d plot 
        
        Parameters
        ----------

        iframe : int 
        topos : boolean    
        fig : 
        ax  :

        Returns
        -------
        
        fig,ax 

        """
        if fig == []:
            fig = plt.figure()
        if ax == []:
            ax = fig.add_subplot(111, projection='3d')
        if not topos:    
            ax.scatter(self.d[0, :, iframe], self.d[1, :, iframe], self.d[2, :, iframe],color=col)
        else:
            ax.scatter(self.topos[0, :], self.topos[1, :], self.topos[2, :],color=col)

        for k,e in enumerate(self.g.edges()):
            e0 = e[0]
            e1 = e[1]
            if not topos:
                pA = self.d[:,e0,iframe].reshape(3,1)
                pB = self.d[:,e1,iframe].reshape(3,1)
            else:    
                pA = self.topos[:,e0].reshape(3,1)
                pB = self.topos[:,e1].reshape(3,1)
            ax.plot(np.array([pA[0][0],pB[0][0]]),np.array([pA[1][0],pB[1][0]]),
                    np.array([pA[2][0],pB[2][0]]),zdir='z',c=col)
        #ax.auto_scale_xyz([-2,2], [-2, 1], [-2, 2])
        return(fig,ax)
                    
    def show3(self,iframe=0,topos=True,tag=''): 
        """ create geomfile for frame iframe 

        Parameters
        ----------

        iframe : int 
            frame number (useless if topos == True)
        topos : boolean 
            if True show the current body topos
        tag : aditional string for naming file .off (useless if topos==False)

        """
        self.geomfile(iframe=0,topos=topos,tag=tag)
        bdy = geu.Geomlist('body'+str(iframe))
        bdy.show3()

    def geomfile(self,iframe=0,verbose=False,topos=False,tag=''):
        """ create a geomview file from a body configuration 

        Parameters
        ----------

        iframe : int 
        verbose : boolean
        topos : boolean 
        tag : string 

        Notes
        -----



        """
        cyl = geu.Geomoff('cylinder')
        pt = cyl.loadpt()
        if not topos:
            _filebody = str(iframe).zfill(4)+'body.list'
        else:    
            _filebody = tag+'-body.list'
        filebody = pyu.getlong(_filebody,"geom")
        #
        # To be change : The defauly layout is DLR.off
        # 
        filestruc = pyu.getlong('DLR.off',"geom")
        fo = open(filebody,"w")
        fo.write("LIST\n")
        fo.write('{<'+filestruc+'}\n')
        if verbose:
            print ("LIST\n")
        for k,e in enumerate(self.g.edges()):
            e0 = e[0]
            e1 = e[1]
            if not topos:
                pA = self.d[:,e0,iframe].reshape(3,1)
                pB = self.d[:,e1,iframe].reshape(3,1)
            else:    
                pA = self.topos[:,e0].reshape(3,1)
                pB = self.topos[:,e1].reshape(3,1)
            pM = (pA+pB)/2.
            T = geu.onbfromaxe(pA,pB)
            R = self.g[e0][e1]['radius']
            Y = np.hstack((pM,pA,pB,pM+R*T[0,:,0].reshape(3,1),pM+R*T[0,:,1].reshape(3,1),pB+R*T[0,:,0].reshape(3,1)))
            A,B = geu.cylmap(Y)
            ptn = np.dot(A,pt.T)+B
            if not topos:
                _filename = 'edge'+str(k)+'-'+str(iframe)+'.off'
            else:
                _filename = tag+'-edge'+str(k)+'.off'
            filename = pyu.getlong(_filename,"geom")
            cyl.savept(ptn.T,_filename)
            fo.write('{<'+filename+'}\n')
            if verbose:
                print('{<'+filename+'}\n')
        fo.close()



    def updbasis0(self,frameId=0,topos=True):
        """ update basis0

        Parameters
        ----------

        frameId : int 
            default 0 
        topos : boolean     
            default True

        Returns
        -------

        self.basis0 : ndarray (nc,3,3)

        Notes
        -----

        There are as many basis as cylinders (body graph edges) 

        """

        nc = len(self.g.edges())
        #
        # basis0 : nc x 9 
        #
        self.basis0 = np.ndarray(shape=(nc,3,3))

        for k,e in enumerate(self.g.edges()):
            e0 = e[0]
            e1 = e[1]
            if not topos:
                pA = self.d[:,e0,iframe].reshape(3,1)
                pB = self.d[:,e1,iframe].reshape(3,1)
            else:    
                pA = self.topos[:,e0].reshape(3,1)
                pB = self.topos[:,e1].reshape(3,1)
            pM = (pA+pB)/2.
            T = geu.onbfromaxe(pA,pB)
            self.basis0[k,:,:] = T 

    def cylinder_basis_k(self, frameId):
        """
        """
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
        """
        """
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
    """  rotate a cycle of frames by an angle alpha

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
    """  rotate a cycle of frames by an angle alpha

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

    cycle :  walking step cycle (2 step), shape = (3,npoints  = 16, nframes = 126)
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

#    nframes = 126
#    Bc = BodyCylinder()
#    Bc.loadC3D()
#    c10_1t = Bc.d
#    #nx.draw(B.g)
#    fig = plt.figure()
#
#    ax = fig.add_subplot(111, projection='3d')
#    frameID = 30
#    ax.scatter(c10_15[0, :, frameID], c10_15[1, :, frameID], c10_15[2, :, frameID])
#    ax.axis('scaled')
#
#    #ax.scatter(c10_15[frameID,:,0],c10_15[frameID,:,1],c10_15[frameID,:,2] )
#
#    pointID = 13
#
#    plt.figure()
#    plt.plot(c10_15[0, pointID].T, '*--', label='x')
#    plt.plot(c10_15[1, pointID].T, '*--', label='y')
#    plt.plot(c10_15[2, pointID].T, '*--', label='z')
#    plt.legend()
#
#    plt.show()

