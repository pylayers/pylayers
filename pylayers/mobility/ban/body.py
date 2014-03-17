"""

Body Class
===========

This class implements the body model

.. autosummary::
    :toctree: /generated

    Body.__init__
    Body.__repr__
    Body.load
    Body.center
    Body.posvel
    Body.settopos
    Body.setccs
    Body.setdcs
    Body.setacs
    Body.loadC3D
    Body.plot3d
    Body._show3
    Body.show
    Body.show3
    Body.geomfile
    Body.movie
    Body.intersectBody
    Body.body_link
    Body.cylinder_basis_k
    Body.cyl_antenna

Miscelianous Functions
======================

.. autosummary::
    :toctree: /generated

    ChangeBasis
    translate
    rotation
    dist
    Global_Trajectory

"""
import numpy as np
import scipy.stats as sp
import ConfigParser
from pylayers.mobility.ban import c3d
import pylayers.mobility.trajectory as tr
import matplotlib.pyplot as plt
import pylayers.antprop.antenna as ant
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import pdb as pdb
import pylayers.util.pyutil as pyu
import pylayers.util.plotutil as plu
import pylayers.util.geomutil as geu
import pylayers.mobility.ban.DeuxSeg as seg
import doctest
import itertools as itt
try:
    from mayavi import mlab
    from tvtk.tools import visual

except:
    print 'mayavi not installed'


class Body(object):
    """ Class to manage a Body model

    Members
    -------

    ncyl : number of cylinder
    sl   :
    d    :
    topos :
    vtopos :

    Methods
    -------

    load
    center
    posvel
    loadC3D
    settopos
    setccs
    setdcs
    geomfile
    plot3d
    movie
    cylinder_basis_k
    cyl_antenna


    """

    def __init__(self,_filebody='John.ini',_filemocap='07_01.c3d'):
        """
        Parameters
        ----------

        _filebody : string
        _filemocap : string

        """

        self.name = _filebody.replace('.ini','')
        di = self.load(_filebody)
        self.loadC3D(filename=_filemocap,centered=True)

    def __repr__(self):
        st = ''

        st = "My name is : " + self.name + '\n\n'

        for k in self.dev.keys():
            st = st + 'I have a '+self.dev[k]['name']+' device '
            side = str(self.dev[k]['cyl'])[-1]
            if side=='l':
                st = st+'on the left '
            if side=='r':
                st = st+'on the right '
            if side=='u':
                st = st+'on the upper part of '
            if side=='b':
                st = st+'on the lower part of '
            st = st + str(self.dev[k]['cyl'])[0:-1]+'\n'



        if 'topos' not in dir(self):
            st = st+ 'I am nowhere yet\n\n'
        else :
            st = st + 'My centroid position is \n'+ str(self.centroid)+"\n"
        if 'filename' in dir(self):
            st = st +'filename : '+ self.filename +'\n'
        if 'nframes' in dir(self):
            st = st +'nframes : ' + str(self.nframes) +'\n'
        if 'pg' in dir(self):
            st = st + 'Centered : True'+'\n'
        #if 'mocapinfo' in dir(self):
        #    st = st + str(self.mocapinfo)+'\n'
        if 'tmocap' in dir(self):
            st = st + 'Mocap Duration : ' + str(self.Tmocap)+'\n'
        if 'vmocap' in dir(self):
            st = st + 'Mocap Speed : ' + str(self.vmocap)+'\n'

        return(st)


    def load(self,_filebody='John.ini'):
        """ load a body ini file

        Parameters
        ----------

        _filebody : body short filename

        Notes
        -----

        A body .ini file contains 4 sections

        + section [nodes]
        Node number = Node name
        + section [cylinder]
        CylinderId = {'t':tail node number, 'h':head node number , 'r': cylinder' radius}
        + section [antennas]
        Antenna number = {'cylId' cylinder Id, 'l': ,'h': parameter,'a':angle,'filename': filename}

        """
        filebody = pyu.getlong(_filebody,'ini')
        config = ConfigParser.ConfigParser()
        config.read(filebody)
        sections = config.sections()
        di = {}
        for section in sections:
            di[section] = {}
            options = config.options(section)
            for option in options:
                if section=='nodes':
                    di[section][option] = config.get(section,option)
                else:
                    di[section][option] = eval(config.get(section,option))

        keys = map(lambda x : eval(x),di['nodes'].keys())
        self.nodes_Id = {k:v for (k,v) in zip(keys,di['nodes'].values())}


        self.sl = np.ndarray(shape=(len(di['cylinder'].keys()),3))
        self.dcyl = {}
        for cyl in di['cylinder'].keys():
            t = di['cylinder'][cyl]['t']
            h = di['cylinder'][cyl]['h']
            r = di['cylinder'][cyl]['r']
            i = di['cylinder'][cyl]['i']
            self.dcyl[cyl]=i
            #pdb.set_trace()
            #
            # sl : segment list of the body
            # line index of sl corresponds to cylinder id from .ini file
            #
            self.sl[i,:] = np.array([t,h,r])

        self.ncyl = len(di['cylinder'].values())

        # update devices dict
        self.dev={}
        for dev in di['device'].keys():
            self.dev[dev]=di['device'][dev]

        return(di)

    def center(self):
        """ centering the body

        Returns
        -------

        self.pg : center of gravity
        self.vg : velocity
        self.d : set of centered frames
        self.smocap : integrated distance
        self.vmocap : averaged velocity

        Notes
        -----

        The center method creates a centered version of the motion capture data stored in
        self.d
        It also calculates :
        self.smocap : total distance along trajectory
        self.vmocap : averaged speed along trajectory
        Here only the projection of the body centroid in the plan 0xy is calculated

        """
        # self.d : 3 x 16 x Nf
        # self.pg : 3 x Nf
        if not self.centered:
            self.pg = np.sum(self.d,axis=1)/self.npoints
            self.pg[2,:] = 0
            self.d = self.d - self.pg[:,np.newaxis,:]
            # speed vector of the gravity centernp.
            self.vg = self.pg[:,1:]-self.pg[:,0:-1]
            # duplicate last spped vector for size homogeneity
            self.vg = np.hstack((self.vg,self.vg[:,-1][:,np.newaxis]))
            # length of trajectory
            d = self.pg[0:-1,1:]-self.pg[0:-1,0:-1]
            self.smocap = np.cumsum(np.sqrt(np.sum(d*d,axis=0)))
            self.vmocap = self.smocap[-1]/self.Tmocap
            self.centered = True

    def posvel(self,traj,t):
        """
        traj : Tajectory DataFrame
            nx3
        t : float
            trajectory time for evaluation of topos

        Returns
        -------

        kf  : frame integer index
        kt  : trajectory integer index
        vsn : normalized speed vector along motion capture trajectory (source)
        wsn : planar vector orthogonal to vsn
        vtn : normalized speed vector along motion trajectory (target)
        wtn : planar vector orthogonal to wtn

        """
        # t should be in the trajectory time range
        assert ((t>=traj.tmin) & (t<=traj.tmax)),'posvel: t not in trajectory time range'

        sk = traj.distance(t) # covered distance along trajectory at time t
        smax = self.smocap[-1]
        ks = int(np.floor(sk/smax)) # number of sequences
        df = sk - ks*smax # covered distance into the sequence
        kf = np.where(self.smocap>=df)[0][0]

        #tf = self.Tmocap/(1.0*self.nframes) # frame body time sampling period
        #timetraj = traj.time()
        #tt = timetraj[1]-timetraj[0]        # trajectory time sampling period

        kt = int(np.floor(t/traj.ts))        # trajectory time integer index
        # self.pg : 3 x Nframes
        # traj : Nptraj x 3 (t,x,y)


        #
        #  BODY SIDE
        #
        # vs  : speed vector along motion capture frame
        # vsn : unitary speed vector along motion capture frame
        #

        vs = self.pg[0:-1,kf] - self.pg[0:-1,kf-1]
        vsn = vs/np.sqrt(np.dot(vs,vs))
        wsn = np.array([vsn[1],-vsn[0]])

        #
        #
        #  TRAJECTORY SIDE  (Topos)
        #
        #
        # vt : speed vector along trajectory
        #

        vt = np.array([traj['vx'][kt],traj['vy'][kt]])
        vtn = vt/np.sqrt(np.dot(vt,vt))
        wtn = np.array([vtn[1],-vtn[0]])

        # vt = traj[kt+1,1:] - traj[kt,1:]
        # vt = traj[kt+1,1:] - traj[kt,1:]

        return(kf,kt,vsn,wsn,vtn,wtn)

    def settopos(self,traj,t=0,cs=False):
        """ translate the body on a time stamped trajectory

        Parameters
        ----------

        traj : ndarray (3,N)
        t,x,y
        t : float
        time for evaluation of topos (seconds) this value should be in the
        range of the trajectory timestamp

        Returns
        -------

        self.topos
        self.vtopos

        Examples
        --------

        .. plot::
            :include-source:

            >>> import numpy as np
            >>> import pylayers.mobility.trajectory as tr
            >>> import pylayers.mobility.ban.body as body
            >>> import matplotlib.pyplot as plt
            >>> time = np.arange(0,10,0.1)
            >>> v = 4000/3600.
            >>> x = v*time
            >>> y = np.zeros(len(time))
            >>> traj = tr.Trajectory()
            >>> traj.generate()
            >>> John = body.Body()
            >>> John.settopos(traj,2.3)
            >>> fig,ax = John.show(plane='xz',color='b')
            >>> plt.title('xz')
            >>> plt.show()

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
        #
        # kt : trajectory integer index  
        # kf : frame integer index  

        kf,kt,vsn,wsn,vtn,wtn = self.posvel(traj,t)

        psa = np.array([0,0])
        psb = psa + vsn
        psc = psa + wsn

        pta = np.hstack((traj['x'].values[kt],traj['y'].values[kt]))
        ptb = pta + vtn
        ptc = pta + wtn

        self.centroid = pta

        X = np.array([[0,0],[psb[0],psb[1]],[psc[0],psc[1]]]).T
        Y = np.array([[pta[0],pta[1]],[ptb[0],ptb[1]],[ptc[0],ptc[1]]]).T

        a,b = geu.affine(X,Y)

        A = np.eye(3)
        B = np.zeros((3,1))
        A[0:-1,0:-1] = a
        B[0:-1,:] = b

        self.toposFrameId = kf
        #
        # TOPOS = A d + B     d == BODY at kf frame
        #
        self.topos = (np.dot(A,self.d[:,:,kf])+B)

        self.vtopos = np.hstack((vtn,np.array([0])))[:,np.newaxis]

        # if asked for calculation of coordinates systems
        if cs:
            # calculate cylinder coordinate system 
            self.setccs(topos=True)
            # calculate device coordinate system 
            self.setdcs(topos=True)
            # calculate antenna coordinate system 
            self.setacs()




    def setdcs(self, topos = True, frameId =0):
        """ set device coordinate system (dcs) from a topos

        This method evaluates the set of all dcs.
        It provides the information necessary for device placement on
        the body.

        If N is the number of antenna an dcs is an MDA of size 3x4xN

        Parameters
        ----------

        topos : boolean
                default : True
        frameId : int
                default 0

        Returns
        -------

        self.dcs : dictionnary

        Examples
        --------

        .. plot::
            :include-source:

            >>> import numpy as np
            >>> import pylayers.mobility.trajectory as tr
            >>> import pylayers.mobility.ban.body as body
            >>> import matplotlib.pyplot as plt
            >>> time = np.arange(0,10,0.1)
            >>> v = 4000/3600.
            >>> x = v*time
            >>> y = np.zeros(len(time))
            >>> traj = tr.Trajectory()
            >>> traj.generate()
            >>> bc = body.Body()
            >>> bc.settopos(traj,2.3,2)
            >>> bc.setccs(topos=True)
            >>> bc.setdcs()
            >>> bc.show(plane='yz',color='b',widthfactor=80)
            >>> plt.show()

        """
        self.dcs = {}
        for dev in self.dev.keys():

            # retrieving antenna placement information from dictionnary ant
            cylname = self.dev[dev]['cyl']
            Id = self.dcyl[cylname]
            alpha = self.dev[dev]['a']*np.pi/180.
            l = self.dev[dev]['l']
            h = self.dev[dev]['h']

            # getting cylinder information

            #pdb.set_trace()
            #~ a_edge = np.array(l_edge)
            #~ a_data = a_edge[:,2]
            #~ dtk = filter(lambda x: x['id']==str(Id),a_data)
            #~ k = np.where(a_data ==dtk)[0][0]

            #~ kta = ed[0]
            #~ khe = ed[1]
            kta = int(self.sl[int(Id),0])
            khe = int(self.sl[int(Id),1])
            Rcyl = self.sl[int(Id),2]

            if topos == True :
                pta = np.array(self.topos[:,kta])
                phe = np.array(self.topos[:,khe])
            else:
                pta = np.array(self.d[:,kta, frameId])
                phe = np.array(self.d[:,khe, frameId])


            vl = phe - pta
            lmax = np.sqrt(np.dot(vl,vl))

            #CCS = self.ccs[k,:,:]
            CCS = self.ccs[Id,:,:]
            #self.nodes_Id[kta],self.nodes_Id[khe]

            # applying rotation and translation

            Rot = np.array([[np.cos(alpha),-np.sin(alpha),0],[np.sin(alpha),np.cos(alpha),0],[0,0,1]])
            CCSr = np.dot(CCS,Rot)
            neworigin = pta + CCSr[:,2]*(l*lmax) + CCSr[:,0]*(Rcyl+h)
            self.dcs[dev] = np.hstack((neworigin[:,np.newaxis],CCSr))

    def setacs(self):
        """ set antenna coordinate system (acs) from a topos or a set of frames

        """

        self.acs = {}
        for dev in self.dev.keys():
            Rab = self.dev[dev]['T']
            U = self.dcs[dev]
            # extract only orthonormal basis
            Rbg = U[:,1:]
            self.acs[dev]  = np.dot(Rbg,Rab)

    def loadC3D(self, filename='07_01.c3d', nframes=300 ,unit='cm',centered = False):
        """ load nframes of motion capture C3D file

        Parameters
        ----------

        filename : string
        file name
        nframes : int
        number of frames

        Notes
        -----

        The body is centered at the

        """


        #if 'pg' in dir(self):
        # del self.pg

        s, p, f, info = c3d.read_c3d(filename)

        self.mocapinfo = info

        self.filename = filename
        if nframes<>-1:
            self.nframes = nframes
        else:
            self.nframes = np.shape(f)[0]
        #
        # s : prefix
        # p : list of points name
        # f : nframe x npoints x 3
        #

        CM_TO_M = 0.01

        # duration of the motion capture snapshot

        self.Tmocap = self.nframes / info['VideoFrameRate']

        #
        # motion capture data
        #
        # self.d : 3 x npoints x nframes
        #

        self.npoints = 15

        self.d = np.ndarray(shape=(3, self.npoints, self.nframes))

        #if self.d[2,:,:].max()>50:
        # extract only known nodes in nodes_Id
        ind = []
        for i in self.nodes_Id:
            if self.nodes_Id[i]<>'BOTT':
                ind.append(p.index(s[0] + self.nodes_Id[i]))

        # f.T : 3 x npoints x nframe
        #
        # cm to meter conversion if required
        #
        self.d = f[0:nframes, ind, :].T
        if unit=='cm':
            self.d = self.d*CM_TO_M

        #
        # Extension of cylinder
        #

        self.npoints = 16

        pm = (self.d[:, 9, 0] + self.d[:, 10, 0])/2.
        pmf = (self.d[:, 9, :] + self.d[:, 10, :])/2.
        pmf = pmf[:,np.newaxis,:]

        self.d = np.concatenate((self.d,pmf),axis=1)

        #self.nodes_Id[15]='bottom'
        if centered:
            self.centered = False
            self.center()



    def movie(self,**kwargs):
        """ creates a geomview movie

        Parameters
        ----------

        lframe : []
        verbose : False
        topos : True
        wire : True
        ccs : False
        dcs : False
        struc : True
        traj : []
        filestruc:'DLR.off'

        See Also
        --------

        Body.geomfile

        """

        defaults = {'lframe':[],
                    'verbose':False,
                    'topos': True,
                    'wire': True,
                    'ccs': False,
                    'dcs': False,
                    'struc':True,
                    'pattern':False,
                    'traj':[],
                    'filestruc':'DLR.off'
                   }


        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if not kwargs['topos']:
            for k in range(self.nframes):
                self.geomfile(iframe=k,verbose=True)
        else:
            t = kwargs['traj'].time()
            for k,tt in enumerate(t):
                stk = str(k).zfill(6) # for string alignement
                self.settopos(traj=kwargs['traj'],t=tt)
                if kwargs['ccs']:
                    self.setccs(topos=True)
                if kwargs['dcs']:
                    self.setdcs()
                kwargs['tag']=stk
                self.geomfile(**kwargs)

    def plot3d(self,iframe=0,topos=False,fig=[],ax=[],col='b'):
        """ scatter 3d plot
        Parameters
        ----------

        iframe : int
        topos : boolean
        fig :
        ax :

        Returns
        -------
        fig,ax

        """
        if fig == []:
            fig = plt.figure()
        if ax == []:
            ax = fig.add_subplot(111, projection='3d')
        if not topos:
            ax.scatter(self.d[0, :, iframe],
                       self.d[1, :, iframe],
                       self.d[2, :, iframe],color=col)
        else:
            ax.scatter(self.topos[0, :], self.topos[1, :], self.topos[2, :],color=col)

        for k in range(self.sl.shape[0]):
            # e0 : tail node of cylinder segment
            e0 = self.sl[k,0]
            # e1 : head node of cylinder segment
            e1 = self.sl[k,1]
            if not topos:
                pA = self.d[:,e0,iframe].reshape(3,1)
                pB = self.d[:,e1,iframe].reshape(3,1)
            else:
                pA = self.topos[:,e0].reshape(3,1)
                pB = self.topos[:,e1].reshape(3,1)

            ax.plot(np.array([pA[0][0],pB[0][0]]),
                    np.array([pA[1][0],pB[1][0]]),
                    np.array([pA[2][0],pB[2][0]]),zdir='z',c=col)

        #ax.auto_scale_xyz([-2,2], [-2, 1], [-2, 2])
        ax.autoscale(enable=True)
        return(fig,ax)

    def _show3(self,**kwargs):
        """ mayavi visualization

        Parameters
        ----------

        iframe : int
            frame index (default 0 )
        widthfactor : float
            cylinder scaling factor (default 1.0)
        topos : boolean
            show topos if True
        pattern : boolean
            show pattern if True
        ccs : boolean
            show ccs if True
        k : frequency index
            select frequency index for displaying antenna pattern

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.mobility.trajectory import *
            >>> from pylayers.mobility.ban.body import *
            >>> b = Body()
            >>> traj = Trajectory()
            >>> traj.generate()
            >>> b.settopos(traj,t=3,cs=True)
            #>>> b._show3(topos=True,pattern=True)

        """
        defaults = {'iframe' : 0,
                    'widthfactor' : 1.,
                    'topos':False,
                    'pattern':False,
                    'ccs':False,
                    'k':0}

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        args = {}
        for k in kwargs:
            if k not in defaults:
                args[k] = kwargs[k]

        visual.set_viewer(mlab.gcf())

        fId = kwargs['iframe']


        cold = pyu.coldict()
        colhex = cold[kwargs['color']]
        body_color = tuple(pyu.rgb(colhex)/255.)


        for k in range(self.ncyl):

            kta = int(self.sl[k,0])
            khe = int(self.sl[k,1])
            cylrad = self.sl[k,2]
            if kwargs['topos']:
                pta =  np.array([self.topos[0, kta], self.topos[1, kta], self.topos[2, kta]])
                phe =  np.array([self.topos[0, khe], self.topos[1, khe], self.topos[2, khe]])
            else:
                pta =  np.array([self.d[0, kta, fId], self.d[1, kta, fId], self.d[2, kta, fId]])
                phe =  np.array([self.d[0, khe, fId], self.d[1, khe, fId], self.d[2, khe, fId]])

            ax = phe-pta
            cc = (pta+phe)/2.
            l = np.sqrt(np.sum(ax**2))
            cyl = visual.Cylinder(pos=(pta[0],pta[1],pta[2]),axis=(ax[0],ax[1],ax[2]), radius=cylrad*kwargs['widthfactor'],length=l)

            mlab.pipeline.surface(cyl.polydata,color=body_color)
            f.children[-1].name=self.name +' ' +self.idcyl[k]            


            if kwargs['ccs']:

                pt = pta+cylrad*kwargs['widthfactor']*self.ccs[k,:,0]
                pte = np.repeat(pt[:,np.newaxis],3,axis=1)
                mlab.quiver3d(pte[0],pte[1],pte[2],self.ccs[k,0],self.ccs[k,1],self.ccs[k,2],scale_factor=0.2)




        if kwargs['pattern']:
            self.setacs()
            for key in self.dcs.keys():
                Ant =  ant.Antenna(self.dev[key]['file'])

                if not hasattr(Ant,'SqG'):
                    Ant.Fsynth()

                U = self.dcs[key]
                V = Ant.SqG[kwargs['k'],:,:]
                T = self.acs[key]

                Ant._show3(po=U[:,0],
                           T=T,
                           ilog=False,
                           minr=0.01,
                           maxr=0.2,
                           newfig=False,
                           title=False,
                           colorbar=False)

        fig = mlab.gcf()
        mlab.savefig('Body.png',figure=fig)


    def show(self,**kwargs):
        """ show a 2D plane projection of the body

        Parameters
        ----------

        frameiId : int
        plane : string
            'yz' | 'xz' | 'xy'
        widthfactor : int
        topos : boolean
            default False



        """

        defaults = {'frameId' : 0,
                    'plane': 'yz',
                    'widthfactor' : 40,
                    'topos':False}

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        args = {}
        for k in kwargs:
            if k not in defaults:
                args[k] = kwargs[k]

        if kwargs['plane'] == 'yz':
            ax1 = 1
            ax2 = 2
        if kwargs['plane'] == 'xz':
            ax1 = 0
            ax2 = 2
        if kwargs['plane'] == 'xy':
            ax1 = 0
            ax2 = 1

        fId = kwargs['frameId']

        for k in range(self.ncyl):

            kta = self.sl[k,0]
            khe = self.sl[k,1]
            cylrad = self.sl[k,2]
            if kwargs['topos']:
                pta =  np.array([self.topos[ax1, kta], self.topos[ax2, kta]])[:,np.newaxis]
                phe =  np.array([self.topos[ax1, khe], self.topos[ax2, khe]])[:,np.newaxis]
            else:
                pta =  np.array([self.d[ax1, kta, fId], self.d[ax2, kta, fId]])[:,np.newaxis]
                phe =  np.array([self.d[ax1, khe, fId], self.d[ax2, khe, fId]])[:,np.newaxis]

            fig,ax = plu.displot(pta,phe,linewidth = cylrad*kwargs['widthfactor'],**args)
            #try:
            #
            #    pta = np.vstack((pta,np.array([self.d[ax1, kta, fId], self.d[ax2, kta, fId]])[np.newaxis,:]))
            #    phe = np.vstack((phe,np.array([self.d[ax1, khe, fId], self.d[ax2, khe, fId]])[np.newaxis,:]))
            #except:
            #    pta =  np.array([self.d[ax1, kta, fId], self.d[ax2, kta, fId]])[np.newaxis,:]
            #    phe =  np.array([self.d[ax1, khe, fId], self.d[ax2, khe, fId]])[np.newaxis,:]
            #print pta
            #print phe
            #fig,ax = plu.displot(pta,phe,linewidth = cylrad*kwargs['widthfactor'],**args)
        #pdb.set_trace()
        #fig,ax = plu.displot(pta.T,phe.T,**args)

        plt.axis('scaled')
        return(fig,ax)


    def show3(self,**kwargs):
        """ create geomfile for frame iframe

        Parameters
        ----------

        iframe : int
        frame number (useless if topos == True)
        topos : boolean
        if True shows the current body topos
        tag : aditional string for naming file .off (useless if topos==False)

        """

        defaults = { 'iframe': 0,
                    'verbose':False,
                    'topos':False,
                    'tag':'',
                    'wire':False,
                    'ccs':False,
                    'lccs':[],
                    'dcs':False,
                    'struc':False,
                    'pattern':False,
                    'filestruc':'DLR.off'
                  }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        bdy = self.geomfile(**kwargs)
        bdy.show3()

    def geomfile(self,**kwargs):
        """ create a geomview file from a body configuration

        Parameters
        ----------

        iframe : int
        frame id (useless if topos==True)
        verbose : boolean
        topos : boolean
        frame id or topos
        wire : boolean
        body as a wire or cylinder
        ccs : boolean
        display cylinder coordinate system
        cacs : boolean
        display cylinder antenna coordinate system
        acs : boolean
        display antenna coordinate system
        struc : boolean
        displat structure layout
        tag : string
        filestruc : string
        name of the Layout

        Notes
        -----

        This function creates either a 3d representation of the frame iframe
        or if topos==True a representation of the current topos.



        """
        defaults = { 'iframe': 0,
                    'verbose':False,
                    'topos':False,
                    'tag':'',
                    'wire': False,
                    'ccs': False,
                    'lccs': [],
                    'dcs': False,
                    'ldcs': [],
                    'struc':False,
                    'pattern':False,
                    'velocity':False,
                    'filestruc':'DLR.off',
                    'fileant':'defant.vsh3',
                    'k':0 }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        iframe = kwargs['iframe']

        Ant = ant.Antenna(kwargs['fileant'])
        Ant.Fsynth3()

        if kwargs['lccs']==[]:
            lccs = np.arange(11)
        else:
            lccs = kwargs['lccs']

        if not kwargs['wire']:
            # load reference cylinder
            #cyl = geu.Geomoff('cylinder')
            cyl = geu.Geomoff('cylinder')
            ptc = cyl.loadpt()

        if not kwargs['topos']:
            _filebody = str(iframe).zfill(4)+'body'
        else:
            if kwargs['tag']<>'':
                _filebody = kwargs['tag']+'-body'
            else:
                _filebody = 'body'

        bodylist = geu.Geomlist(_filebody,clear=True)
        filestruc = pyu.getlong(kwargs['filestruc'],"geom")

        bodylist.append("LIST\n")

        # add layout
        if kwargs['struc']:
            bodylist.append('{<'+filestruc+'}\n')
        if kwargs['verbose']:
            print ("LIST\n")

        dbody = {}
        for k in range(self.sl.shape[0]):
            # e0 : tail node of cylinder segment
            e0 = self.sl[k,0]
            # e1 : head node of cylinder segment
            e1 = self.sl[k,1]
            Rcyl = self.sl[k,2]

            if not kwargs['topos']:
                pA = self.d[:,e0,iframe].reshape(3,1)
                pB = self.d[:,e1,iframe].reshape(3,1)
                vg = self.vg[:,iframe][:,np.newaxis]
            else:
                pA = self.topos[:,e0].reshape(3,1)
                pB = self.topos[:,e1].reshape(3,1)
                vg = self.vtopos
            pM = (pA+pB)/2.

            if kwargs['wire']:
                dbody[k]=(pA,pB)
            else:
                # affine transformation of cylinder
                T = geu.onb(pA,pB,vg)
                Y = np.hstack((pM,pA,pB,pM+Rcyl*T[0,:,0].reshape(3,1),
                                        pM+Rcyl*T[0,:,1].reshape(3,1),
                                        pB+Rcyl*T[0,:,0].reshape(3,1)))
                # idem geu.affine for a specific cylinder
                #A,B = geu.cylmap(Y,r=2,l=6)
                A,B = geu.cylmap(Y)
                ptn = np.dot(A,ptc.T)+B

                if not kwargs['topos']:
                    _filename = 'edge'+str(k)+'-'+str(kwargs['iframe'])+'.off'
                else:
                    _filename = kwargs['tag']+'-edge'+str(k)+'.off'

                filename = pyu.getlong(_filename,"geom")
                cyl.savept(ptn.T,_filename)
                bodylist.append('{<'+filename+'}\n')
                if kwargs['verbose']:
                    print('{<'+filename+'}\n')

            # display selected cylinder coordinate system
            if kwargs['ccs']:
                if k in lccs:
                    fileccs = kwargs['tag']+'ccs'+str(k)
                    geov = geu.GeomVect(fileccs)
                    pt = pA[:,0]+Rcyl*self.ccs[k,:,0]

                    geov.geomBase(self.ccs[k,:,:],pt=pt,scale=0.1)
                    bodylist.append('{<'+fileccs+'.vect'+"}\n")

        # display antenna cylinder coordinate system
        if kwargs['dcs']:
            for key in self.dcs.keys():
                filedcs = kwargs['tag']+'dcs-'+key
                U = self.dcs[key]
                geoa = geu.GeomVect(filedcs)

                geoa.geomBase(U[:,1:],pt=U[:,0],scale=0.1)
                bodylist.append('{<'+filedcs+'.vect'+"}\n")

        # display antenna pattern

        if kwargs['pattern']:
            self.setacs()
            for key in self.dcs.keys():
                Ant =  ant.Antenna(self.dev[key]['file'])
                if not hasattr(Ant,'SqG'):
                    Ant.Fsynth3()
                U = self.dcs[key]
                _filepatt = kwargs['tag']+'patt-'+key
                geo = geu.Geomoff(_filepatt)
                V = Ant.SqG[kwargs['k'],:,:]
                #T = U[:,1:]
                #Rab = self.dev[key]['T']
                #T = np.vstack((U[:,1+DT[0]],U[:,1+DT[1]],U[:,1+DT[2]]))
                #Rbg = U[:,1:]
                # combine rotation antenna -> body -> global
                #T = np.dot(Rbg,Rab)
                #T = np.eye(3)
                T  = self.acs[key]
                geo.pattern(Ant.theta,Ant.phi,V,po=U[:,0],T=T,ilog=False,minr=0.01,maxr=0.2)
                bodylist.append('{<'+_filepatt+'.off'+"}\n")

        # wireframe body
        if kwargs['wire']:
            if not kwargs['topos']:
                _filebody = 'body'+str(kwargs['iframe'])
            else:
                _filebody = kwargs['tag']+'bwire'
            bodygv = geu.GeomVect(_filebody,clear=True)
            bodygv.segments(dbody,i2d=False,linewidth=5)
            bodylist.append('{<'+_filebody+'.vect}\n')


        if kwargs['velocity']:
            _filevelo = kwargs['tag']+'velo'
            velo = geu.GeomVect(_filevelo,clear=True)
            pta = self.pg
            phe = self.pg+self.vtopos
            ds = {0:(pta,phe)}
            velo.segments(ds,i2d=True,linewidth=3)
            bodylist.append('{<'+_filevelo+'.vect}\n')

        return(bodylist)


    def setccs(self,frameId=0,topos=False):
        """ set cylinder coordinate system

        Parameters
        ----------

        frameId : int
        frame id in the mocap dataframe (default 0)
        topos : boolean
        default False

        Returns
        -------

        self.ccs : ndarray (nc,3,3)
        Notes
        -----

        There are as many frames as cylinders (body graph edges)

        ccs is a MDA (nc x 3 x 3 ) where nc denotes the number of cylinders
        For each cylinder there is an attached coordinate systems

        1st vector
        2nd
        """

        nc = self.ncyl
        #
        # ccs : nc x 3 x 3
        #
        self.ccs = np.empty((nc,3,3))
        for k in range(self.sl.shape[0]):
            # e0 : tail node of cylinder segment
            e0 = int(self.sl[k,0])
            # e1 : head node of cylinder segment
            e1 = int(self.sl[k,1])

            if not topos:
                # pA : tail point
                pA = self.d[:,e0,frameId].reshape(3,1)
                # pB : head point
                pB = self.d[:,e1,frameId].reshape(3,1)
                vg = self.vg[:,frameId][:,np.newaxis]
            else:
                # pA : tail point
                pA = self.topos[:,e0].reshape(3,1)
                # pB : head point
                pB = self.topos[:,e1].reshape(3,1)
                # vtopos : mean topos velociy
                vg = self.vtopos
            pM = (pA+pB)/2.
            # create an orthonormal basis
            # 1 st vector : vg normalized (blue)
            # 2 nd vector : 3 x 1
            # 3 rd vector : PB-PA normalized
            T = geu.onb(pA,pB,vg)
            self.ccs[k,:,:] = T
            #~ # e0 : tail node of cylinder segment
            #~ e0 = e[0]
            #~ # e1 : head node of cylinder segment
            #~ e1 = e[1]
            #~ if not topos:
                #~ # pA : tail point
                #~ pA = self.d[:,e0,frameId].reshape(3,1)
                #~ # pB : head point
                #~ pB = self.d[:,e1,frameId].reshape(3,1)
                #~ vg = self.vg[:,frameId][:,np.newaxis]
            #~ else:
                #~ # pA : tail point
                #~ pA = self.topos[:,e0].reshape(3,1)
                #~ # pB : head point
                #~ pB = self.topos[:,e1].reshape(3,1)
                #~ # vtopos : mean topos velociy
                #~ vg = self.vtopos
            #~ pM = (pA+pB)/2.
            #~ # create an orthonormal basis
            #~ # 1 st vector : vg normalized (blue)
            #~ # 2 nd vector : 3 x 1
            #~ # 3 rd vector : PA-PB normalized
            #~ T = geu.onb(pA,pB,vg)
            #~ self.ccs[k,:,:] = T


    def intersectBody(self,A,B, topos = True, frameId = 0, cyl =[]):
        """

        Parameters
        ----------

        A
        B
        topos
        frameId
        cyl

        Returns
        -------

        intersect : np.array (,ncyl)
            O : AB not intersected by cylinder
            1 : AB intersected by cylinder

        """

        intersect = np.zeros((self.ncyl,1))
        mu = np.zeros((self.ncyl,1))
        lmd = 0.075
        for k in range (self.ncyl):
            if k not in cyl:

                if topos  == True:
                    kta  = self.sl[k,0]
                    khe  = self.sl[k,1]
                    C = self.topos[:,kta]
                    D = self.topos[:,khe]
                else:
                    kta  = self.sl[k,0]
                    khe  = self.sl[k,1]
                    C = self.d[:,kta,frameId]
                    D = self.d[:,khe,frameId]

                alpha, beta,dmin = seg.dmin3d(A,B,C,D)
                if alpha < 0:
                    alpha = 0
                if alpha > 1 :
                    alpha  = 1
                if beta < 0:
                    beta = 0
                if beta > 1:
                    beta = 1
                dmin = np.sqrt(seg.dist (A,B,C,D,alpha,beta)[1])

                if dmin  < self.sl[k,2]:
                    intersect[k]=1


                if 0 < alpha < 1 and 0 < beta < 1 :
                    #print 'dmin = ', dmin  
                    #print 'r = ', self.sl[k,2]  
                    dAB = np.sqrt(sum((A-B)**2))
                    if alpha <> 0:
                        mu[k] =(dmin-self.sl[k,2])*np.sqrt(2/(lmd*dAB*abs(alpha)*abs(1-alpha)))

        return intersect


    def body_link(self, topos = True,frameId = 0):
        """

        Parameters
        ----------
        topos :  boolean
            default True
        frameId : int
            used in case topos == False. Indicates the frame Id.

        Returns
        -------

        link_vis : np.array (,nlinks)
            number of intersected cylinder on the link

        links is a list of couple of strings inticating the different links
        between devices.

        """

        self.links = list(itt.combinations(self.dev.keys(),2))
        n_link = len(self.links)
        link_vis = np.ndarray(shape = (n_link))

        for k,link in enumerate(self.links):
            A = self.dcs[link[0]][:,0]
            B = self.dcs[link[1]][:,0]

            inter  = self.intersectBody(A,B, topos=topos,frameId = frameId, cyl =[])[0]


            link_vis[k] =  sum(inter)
        return link_vis


    def cylinder_basis_k(self, frameId):
        """

        Parameters
        ----------

        frameId : int

        """
        nc = self.c.shape[0]
        self.basisk = np.ndarray(shape=(nc, 9))
        for i in range(nc):
            u0 = self.ccs[i, 0:3]
            v0 = self.ccs[i, 3:6]
            w0 = self.ccs[i, 6:]
            v1 = self.c[i, 4:7, frameId] - self.c[i, 1:4, frameId]
            v1 = v1 / np.linalg.norm(v1)
            uk, vk, wk = ChangeBasis(u0, v0, w0, v1)
            self.basisk[i, 0:3] = uk
            self.basisk[i, 3:6] = vk
            self.basisk[i, 6:] = wk

    def cyl_antenna(self, cylinderId, l, alpha, frameId=0):
        """

        Parameters
        ----------

        cylinderId : int
        index of cylinder
        l : distance from origin of cylider
        alpha : angle from reference direction
        frameId : frameId

        """
        r = self.c[cylinderId, 7, frameId]

        x = r * np.cos(alpha)
        y = r * np.sin(alpha)
        z = l

        if frameId == 0:
            u0 = self.ccs[cylinderId, 0:3]
            v0 = self.ccs[cylinderId, 3:6]
            w0 = self.ccs[cylinderId, 6:]

        else:
            self.cylinder_basis_k(frameId)
            u0 = self.basisk[cylinderId, 0:3]
            v0 = self.basisk[cylinderId, 3:6]
            w0 = self.basisk[cylinderId, 6:]

        self.ant = x.reshape((len(x)), 1) * u0 + \
                   y.reshape((len(y)), 1) * w0 + \
                   z.reshape((len(z)), 1) * v0


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
if __name__ == '__main__':
    # plt.ion()
    # doctest.testmod()
    bd = Body(_filebody='John.ini')
    lt = tr.importsn()
    traj = tr.Trajectory()
    traj.generate()
    bd.settopos(traj,0.3,cs=True)
    #bd.setccs(topos=True)
    #bd.setdcs()
    #bd.show3(k=46,wire=True,dcs=True,topos=True,pattern=True)
    bd._show3(k=46,ccs=True,topos=True,pattern=True,newfig=False)
    #bd.show3(wire=True,dcs=True,topos=True)
    #bd.show3(wire=False,dcs=True,topos=True)
    #bd.movie(traj=lt[0],wire=True,dcs=False,pattern=False,filestruc='TA-Office.off')
