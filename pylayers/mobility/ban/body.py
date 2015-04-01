# -*- coding:Utf-8 -*-
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
import os
import copy
from pylayers.mobility.ban import c3d
import pylayers.mobility.trajectory as tr
import matplotlib.pyplot as plt
import pylayers.antprop.antenna as ant
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import networkx as nx
import pdb as pdb
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.plotutil as plu
import pylayers.util.geomutil as geu
import pylayers.mobility.ban.DeuxSeg as seg
import doctest
import itertools as itt
from pylayers.util.project import *
try:
    from mayavi import mlab
    from tvtk.tools import visual

except:
    print 'mayavi not installed'


class Body(PyLayers):
    """ Class  to manage a Body model

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

    def __init__(self,
                 _filebody='John.ini',
                 _filemocap=[],
                 _filewear = [],
                  traj=[],
                  unit='cm',
                  loop=True,
                 centered=True,
                 multi_subject_mocap=False,
                 color='white'):
        """ object constructor

        Parameters
        ----------

        _filebody : string
        _filemocap : string
        unit : str
            unit of the mocap file 'm'|'cm'|'mm'
        _filewear : string
        traj : tr.Trajectory
        loop : bool
            True : indicate if the mocap file is used as a sequence to be looped on a trajectory (default)
            False: the mocap is self sufficient and describes the complete body movement

        See Also
        --------

        pylayers.mobility.trajectory

        """
        self._multi_subject_mocap=multi_subject_mocap

        # extract name from _filebody

        self.name = _filebody.replace('.ini','')
        di = self.load(_filebody,_filemocap,unit,_filewear)
        # if _filemocap != []:
        #     if unit==[]:
        #         raise AttributeError('Please set the unit of the mocap file mm|cm|m')
        #     self.loadC3D(filename=_filemocap,unit=unit)
        #
        #  When the motion capture is on the correct glabal coordinate system
        #  as for example for data coming from the CORMORAN measurement
        #  campaign it is not required to center the body.
        #  centering makes sense only when using the topos projection
        #
        #
        self.cylfromc3d(centered=centered)

        try:
            self.ccsfromc3d(di)
            self.mocapccs=True
        except:
            self.mocapccs=False

        if isinstance(traj,tr.Trajectory):
            self.traj=traj
        self.centered=centered
        self.mocap_loop=loop

        self.color=color
        # otherwise self.traj use values from c3d file
        # obtain in self.loadC3D

    def __repr__(self):
        st = ''

        st = "My name is : " + self.name + '\n\n'

        for k in self.dev.keys():
            if self.dev[k]['status']=='simulated':
                st = st + 'I have a '+self.dev[k]['name']+' device with id #'+k+' '
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
            else :
                st = st + 'I have a '+self.dev[k]['name']+' device with id #'+k+' on '+\
                            self.dev[k]['radiomarkname']+'\n'




        if 'topos' not in dir(self):
            st = st+ '\nI am nowhere yet\n\n'
        else :
            st = st + '\n@ t=' +str(self.time[self.toposFrameId]) +' (frameID='+ str(self.toposFrameId) +'),\n'+'My centroid position is ' +str(self.pg[:2,self.toposFrameId])+"\n\n"
        if 'filewear' in dir(self):
            st = st +'filewear : '+ self.filewear +'\n'
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
            st = st + 'Mocap Speed : ' + "%2.3f" % self.vmocap+' m/s \n'

        st = st + '\n'

        return(st)


    def load(self,_filebody='John.ini',_filemocap=[],unit=[],_filewear=[]):
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
        + section [wearable]
        + section [mocap]

        """

         # check if local or global path
        if ('/'  in _filebody) or ('\\' in _filebody):
            filebody = _filebody
            ne = os.path.basename(_filebody)
            self.name = os.path.splitext(ne)[0]
        else :
            filebody = pyu.getlong(_filebody,pstruc['DIRBODY'])
        if not os.path.isfile(filebody):
            raise NameError(_filebody + ' cannot be found in'
                             + filebody)

        config = ConfigParser.ConfigParser()
        config.read(filebody)
        sections = config.sections()
        di = {}
        for section in sections:
            di[section] = {}
            options = config.options(section)
            for option in options:
                if section=='cylinder' or option =='nframes':
                    di[section][option] = eval(config.get(section,option))
                else:
                    di[section][option] = config.get(section,option)

        keys = map(lambda x : eval(x),di['nodes'].keys())
        nodes_Id = {k:v for (k,v) in zip(keys,di['nodes'].values())}
        # identifier are always 4 character. otherwise its a list
        fnid = filter(lambda x: len(x[1])>4 , nodes_Id.items())
        for k,v in fnid:
            # clean bracket and coma
            vc = v.split('[')[1].split(']')[0].split(',')
            nodes_Id.update({k:vc})

        self.nodes_Id=nodes_Id

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

        self.idcyl={}
        [self.idcyl.update({v:k}) for k,v in self.dcyl.items()]



        # if a mocap file is given in the config file

        if _filemocap == []:
            unit = di['mocap']['unit']
            nframes = di['mocap']['nframes']
            self.loadC3D(di['mocap']['file'],nframes = nframes, unit = unit)
        else:
            if unit  == []:
                raise AttributeError('Please indicate the unit of the motion capture file')
            self.loadC3D(_filemocap, unit = unit)

        #
        # update devices dict from wearable file
        #
        try :
            del self.dev
        except:
            pass
        self.dev={}


        # read default in ini file
        if _filewear == []:
            devfilename = pyu.getlong(di['wearable']['file'],pstruc['DIRWEAR'])
            self.filewear = di['wearable']['file']
            if not os.path.exists(devfilename):
                raise AttributeError('the wareable file '+di['wearable']['file']+
                                 ' cannot be found in $BASENAME/'+pstruc['DIRWEAR'])
        else : 
            # check if local or global path
            if ('/' or '\\') in _filewear:
                devfilename = _filewear
            else :
                devfilename = pyu.getlong(_filewear,pstruc['DIRWEAR'])
            self.filewear = devfilename

            if not os.path.exists(devfilename):
                raise AttributeError('the wareable file '+ devfilename +
                                 ' cannot be found')

        devconf = ConfigParser.ConfigParser()
        devconf.read(devfilename)
        sections = devconf.sections()
        self.dev = {}
        for section in sections:
            self.dev[section] = {}
            options = devconf.options(section)
            for option in options:
                # non case sensitive in .ini file
                if option=='t':
                    option=option.upper()
                #manage non string data
                try:
                    self.dev[section][option]=eval(devconf.get(section,option))
                except:
                    self.dev[section][option]=devconf.get(section,option)
                if option == 'file':
                    self.dev[section]['ant']=ant.Antenna(self.dev[section]['file'])


        try:
            # mocapprefix : retrieve where the prefix is the body name
            ump = [self.name.lower() in p.lower() for p in self._s]
            if sum(ump) >1:
                # Handle case CorSer (serie=3,day=11)
                self._mocap_prefix='Bernard:'
            else:
                self._mocap_prefix = self._s[ump.index(True)]
        except:
            self._mocap_prefix = self._p[-1].split(':')[0]+':'


        #
        # filter real device and get devices
        #
        rd = dict(filter(lambda x: x[1]['status']== 'real',self.dev.items()))


        # for d in rd :
        #     if self.dev[d]['name'] == 'hikob':
        #         bd = [self.dev[d]['radiomarkname'] in n for n in self._p if not 'TCR' in n]
        #         self.dev[d]['uc3d'] = np.where(bd)[0]
        #     else :
        #         bd = [self.dev[d]['radiomarkname'] in n for n in self._p]
        #         self.dev[d]['uc3d'] = np.where(bd)[0]
        #     if len(self.dev[d]['uc3d']) == 0:
        #         print 'Warning : device ',d, 'not present in mocap'
        #         import ipdb
        #         ipdb.set_trace()



        prefix = ['Bernard:','Bernard_','NicolasCormoran:',
                  'Nicolas_FullBody_ClusterOnly:',
                  'Eric_FullBody_ClusterOnly:',
                  'Jihan_FullBody_ClusterOnly2:',
                  'Jihan_FullBody_ClusterOnly:', # j12s8
                  'Jihan_FullBody_sansClusters:', # j12s21
                  'Nicolas_FullBody_sansClusters:', # j12s21
                  'Eric_FullBody_sansClusters:', # j12s21
                  'Nicolas_FullBody:',
                  'Jihan_FullBody:',
                  'Eric_FullBody:']

        self._mocanodes = self._p


        for p in prefix:
            tmpnode=[]
            for n in self._mocanodes:
                tmpnode.append(n.replace(p,''))
                self._mocanodes = tmpnode


        # 2 remove multiple entries due to orientation marker
        self._mocanodes = [n.split(':')[0] for n in self._mocanodes]

        for d in rd :
            if self.dev[d]['name'] == 'hikob':
                bd = [self.dev[d]['radiomarkname'] in n for n in self._mocanodes if not 'TCR' in n]
                self.dev[d]['uc3d'] = np.where(bd)[0]
            else :
                bd = [self.dev[d]['radiomarkname'] in n for n in self._mocanodes]
                self.dev[d]['uc3d'] = np.where(bd)[0]



        return(di)


    def loadC3D(self, filename='07_01.c3d', nframes=-1 ,unit='cm'):
        """ load nframes of motion capture C3D file

        Parameters
        ----------

        filename : string
            file name
        nframes : int
            number of frames
        unit : str (mm|cm|mm
            unit of c3d file
        rot : list ['x','y','z']
            swap axes of the c3d file
        """


        #if 'pg' in dir(self):
        # del self.pg
        # s, p, f, info = c3d.read_c3d(filename)
        self._s, self._p, self._f, info = c3d.ReadC3d(filename)

        if self._multi_subject_mocap:
            us = [us for us, s in enumerate(self._s) if self.name in s ]
            up = [up for up, p in enumerate(self._p) if self.name in p ]

            if len(us) == 0:
                raise AttributeError(self.name +' is not in the MOCAP file :' +filename)


            # in case of multiple body into the mocap file, 
            # mocap is restricted to nodes belonging to a single body.
            # the body is automatically selected by using the self.name
            # 

            self._f =self._f[:,up,:]
            self._s=[s for s in self._s if self.name in s ]
            self._p=[p for p in self._p if self.name in p ]
            



        self.mocapinfo = info

        self.filename = filename
        if nframes<>-1:
            self.nframes = nframes
        else:
            self.nframes = np.shape(self._f)[0]
        #
        # s : prefix
        # p : list of points name
        # f : nframe x npoints x 3
        #


        self.unit = unit
        if unit == 'cm':
            self._unit = 1e-2
        elif unit == 'mm':
            self._unit = 1e-3
        elif unit == 'm':
            self._unit = 1.
        else :
            raise AttributeError('unit'+unit + 'not recognized')
        # duration of the motion capture snapshot


        self._f=self._f*self._unit

        self.Tmocap = self.nframes / info['VideoFrameRate']

        # time base of the motion capture file (sec)
        self.time = np.linspace(0,self.Tmocap,self.nframes)

    def ccsfromc3d(self,config):
        """ Create ccs from C3D file

        Parameters
        ----------

        config : dictionnary

        """


        # dmn = dictionnary of mocap nodes position in self._p
        # for further ccs from marker creation

        self._dmn={n:un for un,n in enumerate(self._mocanodes)}
        self._ccs=np.empty((11,3,3,self.nframes))
        # T10 5 strn 7

        for k,v in config['ccs'].items():
            # clean bracket and coma
            vc = v.split('[')[1].split(']')[0].split(',')
            # get position in uc3d of marker
            uccs=map(lambda x: self._dmn[x],vc)


            # 1 vector carried by cylinder axis
            # 1.1 get cylinder number related to body part k
            upart = config['cylinder'][k]['i']
            # 1.2 get tail and head position in self.d
            kta = self.sl[upart,0].astype(int)
            khe = self.sl[upart,1].astype(int)
            # 1.3 create cylinder axis vector
            ca = self.d[:,kta,:]-self.d[:,khe,:]

            # 2 . create 2 extra vectors
            # 2.1 determine their positions
            # pccs = position of cylinder coordinates system (Nframe x Npts x 3)
            # determine associated vetors
            pccs = self._f[:,uccs,:]

            # vccs = vectors of cylinder coordinates system (3 x 2(Npt) x Nframe)
            vccs = self.d[:,kta,np.newaxis,:] - pccs[:,:,:].T
            # vccs = pccs[:,0,np.newaxis,:]-pccs[:,1:,:]
            # vccs = vectors of cylinder coordinates system (3 x 3(Npt) x Nframe)
            # import ipdb
            # ipdb.set_trace()

            vccs=np.concatenate((ca[:,np.newaxis,:],vccs),axis=1)
            self._ccs[self.dcyl[k],:,:,:]=geu.gram_schmidt(vccs)




            # check ccs continuity
            #~ W=self._ccs[self.dcyl[k],:,:,:]
            #~ W1=W[:,:,:-1]
            #~ W2=W[:,:,1:]
            #~ W1r=np.rollaxis(W1,2)
            #~ W2r=np.rollaxis(W2,2) 
            #~ W2ri=np.linalg.inv(W2r)
            #~ R=np.einsum('ijk,ikl->ijl',W1r,W2ri)
            #~ assert len(np.where(np.linalg.det(R)<1e-9)[0]) <1

    def cylfromc3d(self,centered = False):
        """ Create cylinders from C3D file

        Parameters
        ----------

        centered : boolean
        """
        #
        # motion capture data
        #
        # self.d : 3 x npoints x nframes
        #

        # number of points is determine by the ini file
        self.npoints = len(self.nodes_Id)

        # self.d = np.ndarray(shape=(3, self.npoints, self.nframes))

        #if self.d[2,:,:].max()>50:
        # extract only known nodes in nodes_Id
        self.d = np.zeros((3, self.npoints, self.nframes))
        for i in self.nodes_Id:
            # node name = 4 characters
            if not isinstance(self.nodes_Id[i],list) :
                idx = self._p.index(self._mocap_prefix + self.nodes_Id[i])
                self.d[:,i,:] = self._f[0:self.nframes, idx, :].T
            # perform center of mass of the nodes

            else :

                lnid = len(self.nodes_Id[i])
                for k in range(lnid):

                    nodename = self.nodes_Id[i][k].replace(' ','')

                    idx = self._p.index(self._mocap_prefix + nodename)
                    try:
                        tmp = tmp +self._f[0:self.nframes, idx, :].T
                    except:
                        tmp = self._f[0:self.nframes, idx, :].T
                self.d[:,i,:] = tmp / (1.*lnid)
                del tmp


        # f.T : 3 x npoints x nframe
        #
        # cm to meter conversion if required
        #


        self.d = self.d
        self.pg = np.sum(self.d,axis=1)/self.npoints
        self.pg[2,:] = 0
        #self.nodes_Id[15]='bottom'
        if centered:
            self.centered = False
            self.center()

        self.init_traj()




    def network(self):
        """ evaluate network topology and dynamics

        This function evaluates distance , velocity and acceleration of the
        radio network

        self.D2 : distances between radio nodes
        self.V2 : velocities between radio nodes
        self.A2 : accelerations between radio nodes

        """
        self.ddev = {}
        tdev = []
        for k in self.dev:
            self.ddev[k] = self.dev[k]['radiomarkname']
            tdev.append(self.dev[k]['uc3d'][0])
        tdev = np.array(tdev)

        Net = self._f[:,tdev,:]
        D = Net[:,:,np.newaxis,:]-Net[:,np.newaxis,:,:]
        V = (D[1:,:,:,:]-D[0:-1,:,:,:])/0.01
        A = (V[1:,:,:,:]-V[0:-1,:,:,:])/0.01
        Nt = D.shape[0]
        Nd = D.shape[1]
        D2 = np.sqrt(np.sum(D*D,axis=3))
        self.D2 = D2.reshape(Nt,Nd,Nd)
        V2 = np.sqrt(np.sum(V*V,axis=3))
        self.V2 = V2.reshape(Nt-1,Nd,Nd)
        A2 = np.sqrt(np.sum(A*A,axis=3))
        self.A2 = A2.reshape(Nt-2,Nd,Nd)



    def rdpdf(self):
        """ real device position dataframe
        """
        # dictionary of device dataframe
        df={}
        {df.update(
            {d:pd.DataFrame(
                columns=['dev_id','dev_x','dev_y','dev_z'],index=self.traj.index)})
            for d in self.dev.keys()}

        for d in self.dev:
            df[d]['dev_id'] = d

            try:
                df[d]['dev_x'] = np.mean(self._f[:len(self.time)-2,self.dev[d]['uc3d'],0],axis=1)
                df[d]['dev_y'] = np.mean(self._f[:len(self.time)-2,self.dev[d]['uc3d'],1],axis=1)
                df[d]['dev_z'] = np.mean(self._f[:len(self.time)-2,self.dev[d]['uc3d'],2],axis=1)
            except:
                df[d]['dev_x'] = self._f[:len(self.time)-2,self.dev[d]['uc3d'][0],0]
                df[d]['dev_y'] = self._f[:len(self.time)-2,self.dev[d]['uc3d'][0],1]
                df[d]['dev_z'] = self._f[:len(self.time)-2,self.dev[d]['uc3d'][0],2]
            # gather all devices in a single dataframe:
            addf = pd.DataFrame()
            for d in df:
                addf = pd.concat([addf,df[d]])
        addf = addf.sort_index()
        return addf

    def dpdf(self,tr=[],tunit='ns',poffset=False):
        """ device position dataframe
        return a dataframe with body and devices positions along the self.traj

        Parameters
        ----------

        tr : ndarray
            timerange

        Returns
        -------

        cdf: pd.DataFrame
            complete device data frame

        Example
        -------

        >>> from pylayers.mobility.ban.body import *
        >>> T = tr.Trajectories()
        >>> T.loadh5()
        >>> B=Body(traj=T[0])
        >>> cdf = B.dpdf()
        """

        if not isinstance(tr,np.ndarray):
            traj = self.traj
        else :
            traj = self.traj.copy()
            tstart = tr[0]
            tstop = tr[-1]
            tstep = tr[1]-tr[0]
            sf = traj.ts/tstep
            traj = traj.resample(sf = sf, tstart = tstart, tstop = tstop)
        # dictionary of device dataframe
        df={}
        {df.update(
            {d:pd.DataFrame(
                columns=['dev_id','dev_x','dev_y','dev_z'],index=traj.index)})
            for d in self.dev.keys()}


        dp=[]
        for it,t in enumerate(traj.time()):
            self.settopos(traj = traj,t=t,cs=True)
            dp.append(np.array(self.getdevp(df.keys())))
        dp =np.array(dp)
        for ud,d in enumerate(df.keys()):
            df[d]['dev_id']=d
            df[d].ix[:,['dev_x','dev_y','dev_z']]=dp[:,ud,:]



            # for ud,d in enumerate(df.keys()):
            #     df[d].ix[it,['dev_id']]=d
            #     df[d].ix[it,['dev_x','dev_y','dev_z']]=dp[ud]


        # gather all devices in a single dataframe:
        addf = pd.DataFrame()
        for d in df:
            addf = pd.concat([addf,df[d]])


        # join device dataframe with mobility data frame
        ddf = traj.join(addf)
        ddf['name'] = self.name
        # complete dataframe
        ddf['timestamp']= map(lambda x: str(x.hour).zfill(2) + ':' + str(x.minute).zfill(2) +  ':' + str(x.second).zfill(2) + '.' + str(x.microsecond).zfill(2)[:3],ddf.index)
        if tunit == 'ns':
            ddf['timestamp']= map(lambda x: x.microsecond*1e3+x.second*1e9+60*1e9*x.minute+3600*1e9*x.hour,ddf.index)

        if poffset:
            mx = min(min(ddf['x']),min(ddf['dev_x']))
            ddf['x']=ddf['x']-mx
            ddf['dev_x']=ddf['dev_x']-mx

            my = min(min(ddf['y']),min(ddf['dev_y']))
            ddf['y']=ddf['y']-my
            ddf['dev_y']=ddf['dev_y']-my

            mz = min(min(ddf['z']),min(ddf['dev_z']))
            ddf['z']=ddf['z']-mz
            ddf['dev_z']=ddf['dev_z']-mz

        return ddf

    def export_csv(self, unit = 'mm',df = [],_filename ='default.csv', col =['dev_id', 'dev_x', 'dev_y', 'dev_z', 'timestamp'],**kwargs):
        """
        """
        if _filename == 'default.csv':
            _filename = self.name + '.csv'
        filename =pyu.getlong(_filename,pstruc['DIRNETSAVE'])
        if isinstance(df,pd.DataFrame):
            ddf = df
        else :
            ddf = self.dpdf(**kwargs)

        ldf = ddf[col]
        ldf.rename(columns={'dev_id':'id',
                            'dev_x':'x',
                            'dev_y':'y',
                            'dev_z':'z'},inplace=True)

        if unit == 'm':
            _unit = 1.
        if unit == 'cm':
            _unit = 1e2
        if unit == 'mm':
            _unit = 1e3

        ldf.loc[:,'x']=ldf.loc[:,'x']*_unit
        ldf.loc[:,'y']=ldf.loc[:,'y']*_unit
        ldf.loc[:,'z']=ldf.loc[:,'z']*_unit


        ldf.to_csv(filename, sep = ' ',index=False)

        return ldf

    def init_traj(self):
        """ create trajectory object from given trajectory or mocap 
        """

        # speed vector of the gravity centernp.
        self.vg = self.pg[:,1:]-self.pg[:,0:-1]
        # duplicate last spped vector for size homogeneity
        self.vg = np.hstack((self.vg,self.vg[:,-1][:,np.newaxis]))
        # length of trajectory
        d = self.pg[0:-1,1:]-self.pg[0:-1,0:-1]
        # creates a trajectory associated to mocap file
        self.traj = tr.Trajectory()
        self.traj.generate(t=self.time,pt=self.pg.T,name=self.filename)
        self.smocap = np.cumsum(np.sqrt(np.sum(d*d,axis=0)))
        self.vmocap = self.smocap[-1]/self.Tmocap

    def center(self,force=False):
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

        if not self.centered or force:
            self.d = self.d - self.pg[:,np.newaxis,:]
            self.centered = True


    def posvel(self,traj,t):
        """ calculate position and velocity

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

        Notes
        -----

        This funtion takes as argument a trajectory which is a panda dataframe
        and a time value in the time scale of the trajectory.

        t value should of course be in the interval between trajecroty
        extremal times tmin and tmax.

        smax is the maximum distance covered in the whole motion capture
        sequence.

        sk is the distance covered from the begining of the trajectory until
        the trajectory time t.

        The ratio between those 2 distance is  rounded to the nearest integer.

        :math:`\delta = s_k -\lceil \frac{s_k}{s_{max}} s_{max}`

        k_f is the index of the topos motion capture into the MOCAP 


         _____________________________________________________
        |__________|__________|____________|___________|_____
                                     kf=2
        tmin                                                   tmax
        0         smax               sk

        """
        # t should be in the trajectory time range
        assert ((t>=traj.tmin) & (t<=traj.tmax)),'posvel: t not in trajectory time range'

        sk = traj.distance(t) # covered distance along trajectory at time t
        smax = self.smocap[-1]


        ks = int(np.floor(sk/smax)) # number of full MOCAP sequences of frames

        df = sk - ks*smax  # covered distance into the sequence

        kf = np.where(self.smocap>=df)[0][0]

        #tf = self.Tmocap/(1.0*self.nframes) # frame body time sampling period
        #timetraj = traj.time()
        #tt = timetraj[1]-timetraj[0]        # trajectory time sampling period

        kt = int(np.floor((t-traj.tmin)/traj.ts))        # trajectory time integer index
        # self.pg : 3 x Nframes
        # traj : Nptraj x 3 (t,x,y)


        #
        #  BODY SIDE
        #
        # vs  : speed vector along motion capture frame
        # vsn : unitary speed vector along motion capture frame
        #

        vs = self.pg[0:-1,kf] - self.pg[0:-1,kf-1]
        nvs = np.sqrt(np.dot(vs,vs))
        if nvs != 0 :
            vsn = vs/nvs
        else :
            vsn = vs
        wsn = np.array([vsn[1],-vsn[0]])

        #
        #
        #  TRAJECTORY SIDE  (Topos)
        #
        #
        # vt : speed vector along trajectory
        #

        vt = np.array([traj['vx'][kt],traj['vy'][kt]])
        nvt = np.sqrt(np.dot(vt,vt))
        if nvt != 0:
            vtn = vt/nvt
        else :
            vtn=vt
        wtn = np.array([vtn[1],-vtn[0]])

        # vt = traj[kt+1,1:] - traj[kt,1:]
        # vt = traj[kt+1,1:] - traj[kt,1:]

        return(kf,kt,vsn,wsn,vtn,wtn)

    def time2frame(self,t):
        return np.where(self.time<=15)[0][-1]

    def frame2time(self,f):
        return self.time[f]

    def settopos(self,traj=[],t=0,cs=True,treadmill=False,p0=np.array(([0.,0.]))):
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
        this method takes as argument a trajectory and a time value t in the
        trajectory time-scale.


        See Also
        --------

        pylayers.util.geomutil.affine

        """

        #
        # psa : origin source
        # psb = psa+vsn : a point in the direction of pedestrian motion
        #
        # pta : target translation
        # ptb = pta+vtn : a point in the direction of trajectory
        #
        # kt : trajectory integer index
        # kf : frame integer index


        if not isinstance(traj,tr.Trajectory):
            traj = self.traj


        kf,kt,vsn,wsn,vtn,wtn = self.posvel(traj,t)
        if self.mocap_loop:

            psa = np.array([0,0])
            psb = psa + vsn
            psc = psa + wsn
            if treadmill:
                pta=p0
            else: 
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

            #
            # TOPOS = A d + B     d == BODY at kf frame
            #
            self.topos = (np.dot(A,self.d[:,:,kf])+B)
        else :

            self.topos = self.d[:,:,kf]

        self.vtopos = np.hstack((vtn,np.array([0])))[:,np.newaxis]
        self.toposFrameId = kf

        # self.traj=traj

        kta = self.sl[:,0].astype(int)
        khe = self.sl[:,1].astype(int)
        self._pta = np.array([self.topos[0, kta], self.topos[1, kta], self.topos[2, kta]])
        self._phe = np.array([self.topos[0, khe], self.topos[1, khe], self.topos[2, khe]])
        # if asked for calculation of coordinates systems


        if cs:
            self.setcs(topos=True)


    def setcs(self, topos = True, frameId =0):
        """ set coordinates systems from a topos or a frame id

        Parameters
        ----------

        topos : boolean
                default : True
        frameId : int
                default 0

        See Also
        --------

        pylayers.mobility.ban.body.setccs()
        pylayers.mobility.ban.body.setdcs()
        pylayers.mobility.ban.body.setacs()

        """
        # calculate cylinder coordinate system 
        self.setccs(topos=topos,frameId=frameId)
        # calculate device coordinate system 
        self.setdcs(topos=topos,frameId=frameId)
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
            if self.dev[dev]['status'] == 'simulated':
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

            else :
                if 'toposFrameId'  in dir(self):
                    fId = self.toposFrameId
                else :
                    fId = frameId
                if len(self.dev[dev]['uc3d']) > 2:
                    # mp : marker pos
                    mp = self._f[fId,self.dev[dev]['uc3d'],:]
                    vm = np.diff(mp,axis=0)
                    vm = vm[:3]
                    mvm = np.sqrt(np.sum((vm*vm),axis=0))
                    T = vm/mvm

                    T[:,2]=np.cross(T[:,0],T[:,1])
                    T[:,1]=np.cross(T[:,0],T[:,2])
                    Tn = T/np.sqrt(np.sum(T*T,axis=0))
                    mp0 = mp[0]
                    # self.dcs[dev] = np.hstack((mp0[:,np.newaxis],T))
                else:
                    # associated cylinder
                    # if not self.dev[dev].has_key('asscyl'):
                    #     # find the closest cylinder to the device
                    #     c0=self.sl[:,0].astype(int)
                    #     c1=self.sl[:,1].astype(int)
                    #     pta = self.d[:,c0,0]
                    #     phe = self.d[:,c1,0]
                    #     import ipdb
                    #     ipdb.set_trace()
                    #     # vector tail head
                    #     th = phe - pta
                    #     thl =  np.sqrt(np.sum(th**2,axis=0))
                    #     # vector tail device
                    #     de = self._f[0,self.dev[dev]['uc3d'],:]
                    #     td = pta - de[0,:,np.newaxis]
                    #     tdl = np.sqrt(np.sum(td**2,axis=0))

                    #     do = [np.dot(th[:,i],td[:,i]) for i in range(th.shape[1])]
                    #     # distance matrix
                    #     d = abs(tdl*np.sin(do/(tdl*thl)))
                    #     # closest cylinder
                    #     md = np.where(min(d)==d)[0]
                    #     self.dev[dev]['asscyl']= md[0]

                    if not self.dev[dev].has_key('asscyl'):
                        # find the closest cylinder to the device
                        c0=self.sl[:,0].astype(int)
                        c1=self.sl[:,1].astype(int)
                        pta = self.d[:,c0,0]
                        phe = self.d[:,c1,0]
                        
                        de = self._f[0,self.dev[dev]['uc3d'][0],:][np.newaxis]
                        dtad = np.sqrt(np.sum((pta-de.T)**2,axis=0))
                        dhed = np.sqrt(np.sum((phe-de.T)**2,axis=0))
                        mta = np.min(dtad)
                        mhe = np.min(dhed)
                        # select the smallest distance as best candidate
                        if mta < mhe : 
                            um = np.where(dtad==mta)[0]
                        else :
                            um = np.where(dhed==mhe)[0]
                        self.dev[dev]['asscyl']= um[0]


                    mp0 = np.mean(self._f[fId,self.dev[dev]['uc3d'],:],axis=0)
                    Tn = self.ccs[self.dev[dev]['asscyl'],:,:]

                self.dcs[dev] = np.hstack((mp0[:,np.newaxis],Tn))


    def setacs(self):
        """ set antenna coordinate system (acs) from a topos or a set of frames

        """

        self.acs = {}
        for dev in self.dev.keys():
            if True:#self.dev[dev]['status'] == 'simulated':
                Rab = self.dev[dev]['T']
                U = self.dcs[dev]
                # extract only orthonormal basis
                Rbg = U[:,1:]
                self.acs[dev]  = np.dot(Rbg,Rab)
            else :
                self.acs[dev]  = self.dev[dev]['TTT']



    def c3d2traj(self):
        """ convert c3d file to trajectory
        """
        traj = tr.trajectory()
        return traj


    def getdevp(self,id=-1):
        """ get device position

        Parameters
        ----------

        id : str | list
            device id.
        frameId : int
            frameid
        t : float
            time

        Returns
        -------

        device position
        """

        # if frameId != []:
        #     self.setcs(topos=False,frameId=frameId)
        # elif t !=[]:
        #     fId,kt,vsn,wsn,vtn,wtn = self.posvel(self.traj,t)
        #     print fId
        #     self.setcs(topos=False,frameId=fId)
        # else:
        if not 'dcs' in dir(self):
            raise AttributeError('Body\'s dcs not yet set.\
                                  Use settopos() or precise a frameId or time')

        if id == -1:
            return np.array([self.dcs[i][:,0] for i in self.dcs.keys()])
        if isinstance(id,list):
            return np.array([self.dcs[i][:,0] for i in id])
        else :
            return self.dcs[id][:,0]


    def _rdposition(self):
        """ real devices positions
        """

        dp = {}
        dev = self.dev.keys()
        udev = np.array([self.dev[i]['uc3d'][0] for i in dev])
        p = self._f[:,udev,:]
        dist = np.sqrt(np.sum((p[:,:,np.newaxis,:]-p[:,np.newaxis,:,:])**2,axis=3))

        for ud,d in enumerate(dev) :
            dp[d]=p[:,ud,:]
        return dist,dp


    def getdevT(self,id=-1,t=[],frameId=[]):
        """ get device orientation

        Parameters
        ----------


        id : str | list
            device id.
        frameId : int
            frameid
        t : float
            time

        Returns
        -------

        device orientation
        """
        if not 'dcs' in dir(self):
            raise AttributeError('Body\'s dcs not yet set.\
                            Use settopos() or precise a frameId or time')

        if id == -1:
            return np.array([self.dcs[i][:,1:] for i in self.dcs.keys()])
        if isinstance(id,list):
            return np.array([self.dcs[i][:,1:] for i in id])
        else :
            return self.dcs[id][:,1:]


    def chronophoto(self,**kwargs):
        """ chronophotography of the body movement
            (a.k.a. position as a function of time)

        Parameters
        ----------

        tstart : float
            starting time in second
        tend : float
            ending time in second,
        tstep : float
            time step between tstart and tend
        sstep : float
            spatial step (distance between 2 instant)
        planes : list
            list of planes to be displayed ['xz','xy','yz']
        dev : bool
            show devices
        dev : bool
            show devices ids
        See Also
        --------

        http://en.wikipedia.org/wiki/Chronophotography

        """

        defaults = {'tstart':10,
                    'tend':40,
                    'tstep':5,
                    'sstep':2,
                    'planes':['xz','xy','yz'],
                    'figsize':(10,10),
                    'sharex':False
                    }

        fargs={}
        for k in defaults:
            if k not in kwargs:
                fargs[k] = defaults[k]
            else:
                fargs[k] = kwargs.pop(k)

        fstart=np.where(fargs['tstart']<=self.time)[0][0]
        fend=np.where(fargs['tend']<=self.time)[0][0]
        mocaptres = self.Tmocap/self.nframes
        step = int(fargs['tstep']/mocaptres)
        trange=np.arange(fargs['tstart'],fargs['tend'],fargs['tstep'])
        frange=range(fstart,fend,step)

        vstep=np.arange(0,len(frange))*fargs['sstep']

        fig,axs = plt.subplots(nrows =len(fargs['planes']),
                               ncols=1,sharex=fargs['sharex'],
                               figsize=fargs['figsize'])

        if not isinstance(axs,np.ndarray):
            axs=np.array([axs])
        for p,ax in enumerate(axs):
            for uf,f in enumerate(frange):
                fig,ax=self.show(color='b',plane=fargs['planes'][p],
                                 offset=vstep[uf], frameId=f,fig=fig,ax=ax, **kwargs)

            ax.set_aspect('auto')
            ax.set_ylabel('position (m)')
            ax.set_title('Plane ' + fargs['planes'][p])
            ax.set_xlabel('time (s)')
            ax.set_xlim(vstep[0]-fargs['sstep'] ,vstep[-1]+fargs['sstep'])
            if 'z' in fargs['planes'][p]:
                ax.set_ylim(0,2)
            else :
                ax.set_ylim(-2,2)
            tl = ax.get_xticklabels()
            labrange = np.round(10*np.linspace(fargs['tstart']-fargs['tstep'],fargs['tend']+fargs['tstep'],len(tl)))/10.
            ax.set_xticklabels(labrange)

        plt.tight_layout()
        return fig,ax

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

    @mlab.animate(delay=100)
    def anim(self):
        """ animate body

        Example
        -------

        >>> from pylayers.mobility.trajectory import *
        >>> from pylayers.mobility.ban.body import *
        >>> from pylayers.gis.layout import *
        >>> T=Trajectories()
        >>> T.loadh5()
        >>> L=Layout(T.Lfilename)
        >>> B = Body()
        >>> B.settopos(T[0],t=0,cs=True) 
        >>> L._show3()
        >>> B.anim(B)
        """
        self._show3()
        kta = self.sl[:,0].astype(int)
        khe = self.sl[:,1].astype(int)
        t=self.traj.time()

        anim = range(5000,self.nframes,10)

        ##init antennas
        if 'topos' in dir(self):
            Ant = {}
            for key in self.dcs.keys():
                Ant[key]=ant.Antenna(self.dev[key]['file'])
                if not hasattr(Ant[key],'SqG'):
                    Ant[key].Fsynth()
                Ant[key]._show3(po=self.dcs[key][:,0],
                               T=self.acs[key],
                               ilog=False,
                               minr=0.01,
                               maxr=0.2,
                               newfig=False,
                               title=False,
                               colorbar=False)
        while True:
            if 'topos' in dir(self):
                for k in anim:#range(len(t)):
                    self.settopos(t=t[k],cs=True)
                    # connections=zip(range(0,self.ncyl),range(self.ncyl,2*self.ncyl))
                    X=np.hstack((self._pta,self._phe))
                    # s = np.hstack((cylrad,cylrad))
                    self._mayapts.mlab_source.set(x=X[0,:], y=X[1,:], z=X[2,:])
                    for key in self.dcs.keys():
                        x, y, z ,k,scalar = Ant[key]._computemesh(po=self.dcs[key][:,0],
                                                   T=self.acs[key],
                                                   ilog=False,
                                                   minr=0.01,
                                                   maxr=0.2,
                                                   newfig=False,
                                                   title=False,
                                                   colorbar=False)
                        Ant[key]._mayamesh.mlab_source.set(x=x, y=y, z=z,scalar=scalar)
                    yield
            else:
                for k in anim:
                    pta =  np.array([self.d[0, kta, k], self.d[1, kta, k], self.d[2, kta, k]])
                    phe =  np.array([self.d[0, khe, k], self.d[1, khe, k], self.d[2, khe, k]])
                    # connections=zip(range(0,self.ncyl),range(self.ncyl,2*self.ncyl))
                    X=np.hstack((pta,phe))
                    # s = np.hstack((cylrad,cylrad))
                    self._mayapts.mlab_source.set(x=X[0,:], y=X[1,:], z=X[2,:])
                    yield

    @mlab.animate(delay=10)
    def animc3d(self):
        """ animate c3d file

        Example
        -------

        >>> from pylayers.mobility.trajectory import *
        >>> from pylayers.mobility.ban.body import *
        >>> B = Body()
        >>> B.animc3d(B)

        """
        self._plot3d(typ='c3d',text=False)
        # s, p, f, info = c3d.read_c3d(self.filename)
        s, p, f, info = c3d.ReadC3d(self.filename)
        f=f/10.
        fig = mlab.gcf()
        fig.scene.disable_render=False
        fig.scene.anti_aliasing_frames=0
        while True:
            for k in range(self.nframes):
                # s = np.hstack((cylrad,cylrad))
                self._mayapts.mlab_source.set(x=f[k,:,0],
                                              y=f[k,:,1],
                                              z=f[k,:,2])
                yield

    def _plot3d(self,**kwargs):
        """
            display points and their name for body or original C3D file

        Parameters
        ----------

        typ : str (body|c3d)
            choose points to be displayed (body or c3d)
        text: boolean
            diplay nodes name

        """

        defaults={'typ':'body',
                  'text':True,
                  'edge':False,
                  'ncolor':'b',
                  'ecolor':'b',
                  'iframe' : 0
                  }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        fig = mlab.gcf()

        if 'toposFrameId'  in dir(self):
            fId = self.toposFrameId
        else :
            fId = kwargs['iframe']

        cold = pyu.coldict()
        ncolhex = cold[kwargs['ncolor']]
        pt_color = tuple(pyu.rgb(ncolhex)/255.)
        ecolhex = cold[kwargs['ecolor']]
        ed_color = tuple(pyu.rgb(ecolhex)/255.)

        if kwargs['typ'] == 'c3d':
            # s, p, f, info = c3d.ReadC3d(self.filename)
            self._mayapts=mlab.points3d(self._f[fId,:,0],self._f[fId,:,1],self._f[fId,:,2],scale_factor=0.05,opacity=0.5)
            fig.children[-1].__setattr__('name',self.filename )
            if kwargs['text']:
                self._mayaptstxt=[mlab.text3d(self._f[fId,i,0],self._f[fId,i,1],self._f[fId,i,2],self._p[i].replace(self._s[0],''),
                                    scale=(50.*self._unit,50.*self._unit,50.*self._unit),
                                    color=(0,0,0)) for i in range(len(self._p))]

        else :
            kta = self.sl[:,0].astype(int)
            khe = self.sl[:,1].astype(int)
            cylrad = self.sl[:,2]
            if 'topos' in dir(self):
                pta =  np.array([self.topos[0, kta], self.topos[1, kta], self.topos[2, kta]])
                phe =  np.array([self.topos[0, khe], self.topos[1, khe], self.topos[2, khe]])
            else:
                pta =  np.array([self.d[0, kta, fId], self.d[1, kta, fId], self.d[2, kta, fId]])
                phe =  np.array([self.d[0, khe, fId], self.d[1, khe, fId], self.d[2, khe, fId]])

            X=np.hstack((pta,phe))
            s = np.hstack((cylrad,cylrad))
            self._mayapts = mlab.points3d(X[0,:],X[1,:], X[2,:], 
                                5*s ,
                                scale_factor=0.1,
                                resolution=10,
                                color =pt_color)
            if kwargs['edge']:
                connections=zip(range(0,self.ncyl),range(self.ncyl,2*self.ncyl))
                self._mayapts.mlab_source.dataset.lines = np.array(connections)
                tube = mlab.pipeline.tube(self._mayapts, tube_radius=0.005)
                mlab.pipeline.surface(tube,color=ed_color)
                fig.children[-1].__setattr__('name',self.name )
            if kwargs['text']:
                self._mayaptstxt=[mlab.text3d(X[0,i],X[1,i], X[2,i],self.idcyl[i],
                                scale=(50.*self._unit,50.*self._unit,50.*self._unit),
                                color=(0,0,0)) for i in range(self.ncyl)]

    def plot3d(self,iframe=0,topos=False,fig=[],ax=[],col='b'):
        """ scatter 3d plot

        Parameters
        ----------

        iframe : int
        topos : boolean
        fig :
        ax :
        col : string

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
                       self.d[2, :, iframe], color=col)
        else:
            ax.scatter(self.topos[0, :], self.topos[1, :],
                       self.topos[2, :], color=col)

        for k in range(self.sl.shape[0]):
            # e0 : tail node of cylinder segment
            e0 = self.sl[k, 0]
            # e1 : head node of cylinder segment
            e1 = self.sl[k, 1]
            if not topos:
                pA = self.d[:, e0, iframe].reshape(3,1)
                pB = self.d[:, e1, iframe].reshape(3,1)
            else:
                pA = self.topos[:, e0].reshape(3,1)
                pB = self.topos[:, e1].reshape(3,1)

            ax.plot(np.array([pA[0][0],pB[0][0]]),
                    np.array([pA[1][0],pB[1][0]]),
                    np.array([pA[2][0],pB[2][0]]), zdir='z',c=col)

        #ax.auto_scale_xyz([-2,2], [-2, 1], [-2, 2])
        ax.autoscale(enable=True)
        return(fig,ax)

    def _show3(self,**kwargs):
        """ mayavi visualization

        Parameters
        ----------



        name : boolean (False)
            display body name


        Cylinders 
        ---------
        cylinder : boolean (True)
            dispaly body cylinder
        widthfactor : float
            cylinder scaling factor (default 1.0)
        tube_sides : int (6)
            number of sides of polygone to approximate cylinder
        opacity : float (1)
            set body opacity
        ccs : boolean
            show ccs if True

        Devices 
        -------
        dcs : boolean (False)
            show dcs if True
        devlist: list
            list of devices to display 
        devtyp : list ([])
            list of device type ( when multiple RAT)
        devcolor : color ('g')
            device color
        devid : boolean 
            show device id
        devopacity : float (1)
            device opacity 
        devsize : float
            device size
        pattern : boolean (False)
            show pattern if True




        Config 
        -------
        iframe : int
            frame index (default 0 )
        k : frequency index
            select frequency index for displaying antenna pattern

        save : boolean (False)
            save _show3 into file
        returnfig': booleéan (False)
            return mlab.figure instance


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
                    'name':False,
                    'cylinder':True,
                    'widthfactor' : 1.,
                    'tube_sides' : 6,
                    'opacity':1,
                    'pattern':False,
                    'dev':False,
                    'devlist':[],
                    'ccs':False,
                    'dcs':False,
                    'devcolor':'green',
                    'devid':False,
                    'devopacity':1,
                    'devsize':5e1,
                    'devtyp':[],
                    'k':0,
                    'save':False,
                    'mocanodes':False,
                    'returnfig':False}

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        args = {}
        for k in kwargs:
            if k not in defaults:
                args[k] = kwargs[k]

        f = mlab.gcf()
        #visual.set_viewer(f)
        #f.scene.background=(1,1,1)


        fId = kwargs['iframe']


        cold = pyu.coldict()
        colhex = cold[self.color]
        body_color = tuple(pyu.rgb(colhex)/255.)

        kta = self.sl[:,0].astype(int)
        khe = self.sl[:,1].astype(int)
        cylrad = self.sl[:,2]
        if 'topos' in dir(self):
            X=np.hstack((self._pta,self._phe))
        else:
            pta =  np.array([self.d[0, kta, fId], self.d[1, kta, fId], self.d[2, kta, fId]])
            phe =  np.array([self.d[0, khe, fId], self.d[1, khe, fId], self.d[2, khe, fId]])
            X=np.hstack((pta,phe))

        if kwargs['cylinder']:
            connections=zip(range(0,self.ncyl),range(self.ncyl,2*self.ncyl))
            s = np.hstack((cylrad*kwargs['widthfactor'],cylrad*kwargs['widthfactor']))
            #pts = mlab.points3d(X[0,:],X[1,:], X[2,:], 5*s ,
                                                 # scale_factor=0.1, resolution=10)
            self._mayapts = mlab.pipeline.line_source(X[0,:],X[1,:], X[2,:], s ,
                                                 scale_factor=0.001, resolution=10)
            self._mayapts.mlab_source.dataset.lines = np.array(connections)
            tube = mlab.pipeline.tube(self._mayapts, tube_radius=0.05,tube_sides=kwargs['tube_sides'])
            tube.filter.radius_factor = 1.
            tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
            mlab.pipeline.surface(tube, color=body_color,opacity=kwargs['opacity'])
            f.children[-1].__setattr__('name',self.name )
                
        # ax = phe-pta
        # l = np.sqrt(np.sum((ax**2), axis=0))
        # cyl = [visual.Cylinder(pos=(pta[0, i],pta[1, i],pta[2, i]),
        #                        axis=(ax[0, i],ax[1, i],ax[2, i]), 
        #                        radius=cylrad[i]*kwargs['widthfactor'],
        #                        length=l[i], resolution=1) for i in range(self.ncyl)
        #                        ]
        # [mlab.pipeline.surface(cyl[i].polydata, color=body_color) 
        #  for i in range(self.ncyl)]
        # partnames = [self.name +' ' +self.idcyl[k] for k in range(self.ncyl)]
        # [f.children[k].__setattr__('name', partnames[k]+str(k))
        #  for k in range(self.ncyl)]

        if kwargs['mocanodes']:
            center = self.pg[:,fId]
            X=self._f[fId,:,:].T
            mlab.points3d(X[0,:],X[1,:], X[2,:], 
                          scale_factor=20*self._unit, 
                          resolution=10)
            [mlab.text3d(X[0,i],
                                 X[1,i],
                                 X[2,i],self._p[i].split(':')[1],
                                        scale=0.01,
                                        color=(1,0,0)) for i in range(len(self._p)) ]

        if kwargs['name']:
            uupper = np.where(X[2]==X[2].max())[0]
            self._mayaname = mlab.text3d(X[0,uupper][0],X[1,uupper][0],X[2,uupper][0],self.name,scale=0.05,color=(1,0,0))

        if kwargs['dev']:
            colhex = cold[kwargs['devcolor']]
            dev_color = tuple(pyu.rgb(colhex)/255.)

            if kwargs['devlist'] == []:
                dev = self.dev.keys()
            else :
                dev = [d  for d in self.dev if d in kwargs['devlist']]


            if 'dcs' in dir(self):
                X=np.array(self.getdevp(dev)).T
            else:
                udev = [self.dev[i]['uc3d'][0] for i in dev]

                center = self.pg[:,fId]
                X=self._f[fId,udev,:].T-center[:,np.newaxis]

            self._mayadev = mlab.points3d(X[0,:],X[1,:], X[2,:], 
                          scale_factor=kwargs['devsize']*self._unit, 
                          resolution=10, 
                          color = dev_color,
                          opacity=kwargs['devopacity'])
            nodename = self.dev.keys()
            
            if kwargs['devid']:

                if kwargs['devtyp']== []:
                    udt = np.arange(len(nodename))
                else:
                    devtyp = np.unique([self.dev[x]['name'] for x in dev])
                    ln =[filter(lambda x: d in x,nodename) for d in devtyp]
                    udev = [[nodename.index(n) for n in ln[i]] for i in range(len(ln))]

                for dt in kwargs['devtyp']:
                    udt = np.where(devtyp == dt)[0]

                    [mlab.text3d(X[0,i],
                                 X[1,i],
                                 X[2,i],nodename[i],
                                        scale=0.05,
                                        color=(1,0,0)) for i in range(len(nodename)) if i in udev[udt]]


        if kwargs['ccs']:
            # to be improved
            
            col=np.linspace(0,1,11)[:,np.newaxis]*np.ones((11,3))
            col[:,0]=0

            for k in range(len(self.ccs)):
                usl = self.sl[k]
                pt = (self.topos[:,usl[0].astype(int)]+self.topos[:,usl[1].astype(int)])/2.
                pt = pt+cylrad[k]*kwargs['widthfactor']*self.ccs[k, :, 0]
                pte = np.repeat(pt[:,np.newaxis],3,axis=1)
                ccs = mlab.quiver3d(pte[0], pte[1], pte[2],
                              self.ccs[k, 0], self.ccs[k, 1], self.ccs[k, 2],
                              scale_factor=2e2*self._unit,color=tuple(col[k]))



        # for k in range(self.ncyl):

        #     kta = int(self.sl[k,0])
        #     khe = int(self.sl[k,1])
        #     cylrad = self.sl[k,2]
        #     if kwargs['topos']:
        #         pta =  np.array([self.topos[0, kta], self.topos[1, kta], self.topos[2, kta]])
        #         phe =  np.array([self.topos[0, khe], self.topos[1, khe], self.topos[2, khe]])
        #     else:
        #         pta =  np.array([self.d[0, kta, fId], self.d[1, kta, fId], self.d[2, kta, fId]])
        #         phe =  np.array([self.d[0, khe, fId], self.d[1, khe, fId], self.d[2, khe, fId]])


            # cyl = visual.Cylinder(pos=(pta[0],pta[1],pta[2]),axis=(ax[0],ax[1],ax[2]), radius=cylrad*kwargs['widthfactor'],length=l)
            # mlab.pipeline.surface(cyl.polydata,color=body_color)
            # f.children[-1].name=self.name +' ' +self.idcyl[k]            

            #

            #
        if kwargs['dcs']:
            for key in self.dcs.keys():
                U = self.dcs[key]
                pt = U[:,0]
                pte  = np.repeat(pt[:,np.newaxis],3,axis=1)
                dcs = mlab.quiver3d(pte[0],pte[1],pte[2],self.dcs[key][0,1:],self.dcs[key][1,1:],self.dcs[key][2,1:],scale_factor=2e2*self._unit)


        if kwargs['pattern']:
            self.setacs()
            
            for key in self.dcs.keys():
                if not hasattr(self.dev[key]['ant'],'SqG'):
                    self.dev[key]['ant'].Fsynth()
                U = self.dcs[key]
                V = self.dev[key]['ant'].SqG[kwargs['k'],:,:]
                T = self.acs[key]

                self.dev[key]['ant']._show3(po=U[:,0],
                           T=T,
                           ilog=False,
                           minr=0.01,
                           maxr=0.2,
                           newfig=False,
                           title=False,
                           colorbar=False,
                           )
        if kwargs['save']:
            fig = mlab.gcf()
            mlab.savefig('Body.png',figure=fig)

        if kwargs['returnfig']:
            return fig


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
        offset = np.array()
            1,3



        """

        defaults = {'frameId' : 0,
                    'plane': 'yz',
                    'widthfactor' : 10,
                    'dev':False,
                    'devid':False,
                    'topos':False,
                    'offset':0}

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        offset = np.array([kwargs['offset'],0])[:,np.newaxis]
        args = {}
        for k in kwargs:
            if k not in defaults:
                args[k] = kwargs[k]

        if kwargs['plane'] == 'yz' or kwargs['plane'] == 'zy':
            ax1 = 1
            ax2 = 2
        elif kwargs['plane'] == 'xz' or kwargs['plane'] == 'zx':
            ax1 = 0
            ax2 = 2
        elif kwargs['plane'] == 'xy' or kwargs['plane'] == 'yx':
            ax1 = 0
            ax2 = 1
        else :
            raise AttributeError('Incorrect plane')

        fId = kwargs['frameId']

        for k in range(self.ncyl):

            kta = self.sl[k,0]
            khe = self.sl[k,1]
            cylrad = self.sl[k,2]
            if 'topos' in dir(self):
                pta =  np.array([self.topos[ax1, kta], self.topos[ax2, kta]])[:,np.newaxis]
                phe =  np.array([self.topos[ax1, khe], self.topos[ax2, khe]])[:,np.newaxis]
            else:
                pta =  np.array([self.d[ax1, kta, fId], self.d[ax2, kta, fId]])[:,np.newaxis]
                phe =  np.array([self.d[ax1, khe, fId], self.d[ax2, khe, fId]])[:,np.newaxis]

            fig,ax = plu.displot(pta+offset,phe+offset,linewidth = cylrad*kwargs['widthfactor'],**args)

        # display devices
        if kwargs['dev']:
            if 'dcs' in dir(self):
                pdev = np.array(self.getdevp(self.dev.keys()))
                ax.plot(pdev[:,ax1]+offset[0],pdev[:,ax2]+offset[1],'og')
            else :
                iudev = [(i,self.dev[i]['uc3d'][0]) for i in self.dev]
                udev = [i[1] for i in iudev]
                center = self.pg[:,fId]

                ax.plot(self._f[fId,udev,ax1]+offset[0]-center[ax1],
                        self._f[fId,udev,ax2]+offset[1]-center[ax2],'or')
                if kwargs['devid']:
                    for u in iudev:
                        ax.text(self._f[fId,u[1],ax1]+offset[0]-center[ax1],
                                self._f[fId,u[1],ax2]+offset[1]-center[ax2],u[0])



                # print center
                #ax.plot(self._f[fId,udev,1]+offset[0]-center[1],self._f[fId,udev,2]+offset[1]-center[0],'or')

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

        if self.mocapccs :
            if not topos:
                self.ccs=self._ccs[:,:,:,frameId]
            else:
                self.ccs=self._ccs[:,:,:,self.toposFrameId]
        else :
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
        """ intersect Body

        Parameters
        ----------

        A : np.array (3,)
        B : np.array (3,)
        topos : boolean
        frameId : 0
        cyl : list
            exclusion list

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
                    alpha = np.zeros(alpha.shape)
                if alpha > 1 :
                    alpha  = np.ones(alpha.shape)
                if beta < 0:
                    beta = np.zeros(beta.shape)
                if beta > 1:
                    beta = np.ones(beta.shape)
                dmin = np.sqrt(seg.dist (A,B,C,D,alpha,beta)[1])

                if dmin  < self.sl[k,2]:
                    intersect[k]=1
                    break

                # if 0 < alpha < 1 and 0 < beta < 1 :
                #     #print 'dmin = ', dmin
                #     #print 'r = ', self.sl[k,2]
                #     dAB = np.sqrt(sum((A-B)**2))
                #     if alpha <> 0:
                #         mu[k] =(dmin-self.sl[k,2])*np.sqrt(2/(lmd*dAB*abs(alpha)*abs(1-alpha)))

        return intersect

    def intersectBody3(self,A,B, topos = True, frameId = 0):
        """ intersect body new version

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

        loss_dB = 0 # in dB
        loss_lin  =1/10**(loss_dB/10.0)

        for k in [10]:        

            if topos  == True:
                kta  = int(self.sl[k,0])
                khe  = int(self.sl[k,1])
                C = self.topos[:,kta]
                D = self.topos[:,khe]
            else:
                kta  = self.sl[k,0]
                khe  = self.sl[k,1]
                C = self.d[:,kta,frameId]
                D = self.d[:,khe,frameId]

            alpha, beta,dmin = seg.dmin3d(A,B,C,D)
            if alpha < 0:
                alpha = np.zeros(alpha.shape)
            if alpha > 1 :
                alpha  = np.ones(alpha.shape)
            if beta < 0:
                beta = np.zeros(beta.shape)
            if beta > 1:
                beta = np.ones(beta.shape)
            dmin = np.sqrt(seg.dist (A,B,C,D,alpha,beta)[1])

            #~ if dmin  < self.sl[k,2]:
                #~ intersect=1
            lmd = 0.06 

            #~ if 0 < alpha < 1 and 0 < beta < 1 :
            loss1_dB = 0
            loss2_dB = 0 



            if dmin < self.sl[k,2]:

                """
                in this case intersection is True
                """
                #pdb.set_trace()
                dAB = np.sqrt(sum((A-B)**2))
                #nu1 =(self.sl[k,2]-dmin)*np.sqrt((2/lmd)*dAB*abs(alpha)*abs(1-alpha)))
                nu1 =(self.sl[k,2]-dmin)*np.sqrt(2/(lmd*dAB*abs(alpha)*abs(1-alpha)))*0.05
                nu2 =(dmin+self.sl[k,2])*np.sqrt(2/(lmd*dAB*abs(alpha)*abs(1-alpha)))*0.05

                if -0.7 < nu1 :

                    loss1_dB = 6.9 + 20*np.log10(np.sqrt((nu1-0.1)**2+1)+nu1-0.1)
                else :
                    loss1_dB = 0.0

                if -0.7 < nu2 :
                    loss2_dB = 6.9 + 20*np.log10(np.sqrt((nu2-0.1)**2+1)+nu2-0.1)
                else :
                    loss2_dB = 0.0

                loss_dB  =10*np.log10(10**(loss1_dB/10.0)+10**(loss2_dB/10.0))

                loss_lin  =1.0/10**(loss_dB/10.0)

        return loss_lin#, loss_dB, loss1_dB, loss2_dB

    def body_link(self, topos = True,frameId = 0):
        """ body link

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
        """ cylinder basis k

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
        """ cylinder antenna

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

    def _checkdevid(self):
            """ display 
            """
            for k in self.dev:
                print (k,self.dev[k]['uc3d']) 


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
    """ global trajectory

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
    """ change basis

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
    """ evaluate the distance between two points A and B

    Parameters
    ----------

    A : 2D point
    B : 2D point

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





class Cylinder(object):
    """
    """

    def __init__(self,name = 'Meriem_Cylindre:',
                _filemocap='/media/niamiot/DONNEES/svn2/measures/CORMORAN/RAW/12-06-2014/MOCAP/Nav_serie_006.c3d',
                unit='mm',
                color='white'):
        self.name = name
        self.loadC3D(_filemocap,unit=unit)
        self.init_traj()
        self.color=color
        self.settopos(t=0)

    def __repr__(self):
        st = ''

        st = "I am cylinder : " + self.name + '\n\n'

        if 'topos' not in dir(self):
            st = st+ '\nI am nowhere yet\n\n'
        else :
            st = st + '\n@ t=' +str(self.time[self.toposFrameId]) +' (frameID='+ str(self.toposFrameId) +'),\n'+'My centroid position is ' +str(self.pg[:2,self.toposFrameId])+"\n\n"
        
        st = st + '\n'

        return(st)

    def loadC3D(self, filename='07_01.c3d', nframes=-1 ,unit='cm'):
        """ load nframes of motion capture C3D file

        Parameters
        ----------

        filename : string
            file name
        nframes : int
            number of frames
        unit : str (mm|cm|mm
            unit of c3d file
        rot : list ['x','y','z']
            swap axes of the c3d file
        """


        #if 'pg' in dir(self):
        # del self.pg
        # s, p, f, info = c3d.read_c3d(filename)
        self._s, self._p, self._f, info = c3d.ReadC3d(filename)
        us = [us for us, s in enumerate(self._s) if self.name in s ]
        up = [up for up, p in enumerate(self._p) if self.name in p ]

        if len(us) == 0:
            raise AttributeError(self.name +' is not in the MOCAP file :' +filename)


            # in case of multiple body into the mocap file, 
            # mocap is restricted to nodes belonging to a single body.
            # the body is automatically selected by using the self.name
        # 

        self._f =self._f[:,up,:]
        self._s=[s for s in self._s if self.name in s ]
        self._p=[p for p in self._p if self.name in p ]
            



        self.mocapinfo = info

        if nframes<>-1:
            self.nframes = nframes
        else:
            self.nframes = np.shape(self._f)[0]
        #
        # s : prefix
        # p : list of points name
        # f : nframe x npoints x 3
        #


        self.unit = unit
        if unit == 'cm':
            self._unit = 1e-2
        elif unit == 'mm':
            self._unit = 1e-3
        elif unit == 'm':
            self._unit = 1.
        else :
            raise AttributeError('unit'+unit + 'not recognized')
        # duration of the motion capture snapshot


        self._f=self._f*self._unit
        self.d = self._f[0:self.nframes,:,:].T
        self.pg = np.mean(self.d,axis=1)
        self.pg[2,:]=0
        self.radius = np.sqrt(np.sum((self.d[:,0,:]-self.d[:,1,:])**2,axis=0))/2.
        heighestpos = np.max(self.d[2,:,:])
        self.topnode = np.where(self.d[2,:,:]==heighestpos)[0][0]

        self.Tmocap = self.nframes / info['VideoFrameRate']

        # time base of the motion capture file (sec)
        self.time = np.linspace(0,self.Tmocap,self.nframes)


    def settopos(self,t,**kwargs):

        ut = np.where(t<=self.time)[0][0]
        self.topos =self.d[:,:,ut]
        self.top = self.topos[:,self.topnode]
        self.bottom = copy.copy(self.top)
        self.bottom[2]=0
        self.toposradius =self.radius[ut]
        self.toposFrameId = ut

    def init_traj(self):
        """ create trajectory object from given trajectory or mocap 
        """

        # speed vector of the gravity centernp.
        self.vg = self.pg[:,1:]-self.pg[:,0:-1]
        # duplicate last spped vector for size homogeneity
        self.vg = np.hstack((self.vg,self.vg[:,-1][:,np.newaxis]))
        # length of trajectory
        d = self.pg[0:-1,1:]-self.pg[0:-1,0:-1]
        # creates a trajectory associated to mocap file
        self.traj = tr.Trajectory()
        self.traj.generate(t=self.time,pt=self.pg.T,name=self.name)
        self.smocap = np.cumsum(np.sqrt(np.sum(d*d,axis=0)))
        self.vmocap = self.smocap[-1]/self.Tmocap


    def _show3(self,**kwargs):
        """ mayavi visualization

        """
        defaults = {'iframe' : 0,
                    'name':False,
                    'widthfactor' : 1.,
                    'tube_sides' : 6,
                    'opacity':1,
                    'vecdir':True
                    }       

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        args = {}
        for k in kwargs:
            if k not in defaults:
                args[k] = kwargs[k]

        f = mlab.gcf()


        fId = kwargs['iframe']


        cold = pyu.coldict()
        colhex = cold[self.color]
        cyl_color = tuple(pyu.rgb(colhex)/255.)

       
        X=np.vstack((self.top,self.bottom))
       
        connections=(0,1)
    
        s = np.hstack((self.toposradius,self.toposradius))
        #pts = mlab.points3d(X[0,:],X[1,:], X[2,:], 5*s ,
                                             # scale_factor=0.1, resolution=10)
        self._mayapts = mlab.pipeline.line_source(X[:,0],X[:,1], X[:,2], s ,
                                             scale_factor=0.001, resolution=10)
        self._mayapts.mlab_source.dataset.lines = np.array([connections])
        tube = mlab.pipeline.tube(self._mayapts, tube_radius=0.05,tube_sides=kwargs['tube_sides'])
        tube.filter.radius_factor = 1.
        tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
        mlab.pipeline.surface(tube, color=cyl_color,opacity=kwargs['opacity'])
        f.children[-1].__setattr__('name',self.name )


        if kwargs['name']:
            self._mayaname = mlab.text3d(self.top[0],self.top[1],self.top[2],self.name,scale=0.05,color=(1,0,0))


        if kwargs['vecdir']:

            V = self.traj[['vx','vy','vz']].iloc[self.toposFrameId].values
            self._mayavdic =  mlab.quiver3d(self.top[0], self.top[1], self.top[2],
                              V[ 0], V[ 1], V[ 2],
                              scale_factor=2e2*self._unit)


    @mlab.animate(delay=100)
    def anim(self):
        """ animate cylinder

        Example
        -------

        >>> from pylayers.mobility.trajectory import *
        >>> from pylayers.mobility.ban.body import *
        >>> from pylayers.gis.layout import *
        >>> T=Trajectories()
        >>> T.loadh5()
        >>> L=Layout(T.Lfilename)
        >>> B = Body()
        >>> B.settopos(T[0],t=0,cs=True) 
        >>> L._show3()
        >>> B.anim(B)
        """
        self._show3()
        t=self.traj.time()

        anim = range(5000,self.nframes,10)

        
        while True:
            for k in anim:#range(len(t)):
                self.settopos(t=t[k],cs=True)
                # connections=zip(range(0,self.ncyl),range(self.ncyl,2*self.ncyl))
                X=np.vstack((self.top,self.bottom))
                # s = np.hstack((cylrad,cylrad))
                self._mayapts.mlab_source.set(x=X[:,0], y=X[:,1], z=X[:,2])
                yield