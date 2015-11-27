#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""
   Run simulation with full human trajectory

.. currentmodule:: pylayers.simul.simultraj

Simul class
===========

.. autosummary::
    :toctree: generated/

Run simulation and data exploitation
------------------------------------

.. autosummary::
    :toctree: generated/

    Simul.__init__
    Simul.run
    Simul._gen_net
    Simul.evaldeter
    Simul.evalstat
    Simul.show

Loading and Saving
------------------

.. autosummary::
    :toctree: generated/

    Simul._saveh5
    Simul._loadh5
    Simul.savepd
    Simul.loadpd

Extraction
----------

.. autosummary::
    :toctree: generated/

    Simul.get_link
    Simul.get_value

Miscellanous
------------

.. autosummary::
    :toctree: generated/

    Simul.setttime
    Simul.replace_data
    Simul.check_exist
    Simul.get_df_from_link
    Simul.update_pos

See Also
--------

pylayers.simul.link


"""
import doctest
import numpy as np
import copy
import matplotlib.pylab as plt
import pylayers.util.pyutil as pyu
import pylayers.signal.waveform as wvf
from pylayers.signal.device import Device

# Handle Layout
from pylayers.gis.layout import Layout
# Handle VectChannel and ScalChannel
from pylayers.antprop import antenna
from pylayers.network.network import Network
from pylayers.simul.link import *
from pylayers.measures.cormoran import *

# Handle directory hierarchy
from pylayers.util.project import *
# Handle UWB measurements
import pylayers.mobility.trajectory as tr
from pylayers.mobility.ban.body import *
from pylayers.antprop.statModel import *
import pandas as pd
import csv

class Simul(PyLayers):
    """
    Link oriented simulation

    A simulation requires :

        + A Layout
        + A Person
        + A Trajectory

    or a CorSer instance

    Members
    -------

    dpersons : dictionnary of persons (agent)
    dap : dictionnary of access points

    Methods
    -------

    load_simul : load configuration file
    load_Corser : load a Corser file
    _gen_net : generate network and asociated links
    show : show layout and network
    evaldeter : run simulation over time


    """

    def __init__(self, source ='simulnet_TA-Office.h5',verbose=False):
        """ object constructor

        Parameters
        ----------

        source : string
            h5 trajectory file default simulnet_TA-Office.h5

        verbose : boolean

        Notes
        -----

        The simultraj has a dataframe


        """


        # self.progress = -1  # simulation not loaded
        self.verbose = verbose
        self.cfield = []
        self.dpersons = {}

        self.dap = {}
        self.Nag = 0
        self.Nap = 0

        # source analysis
        if isinstance(source,str):
            self.filetraj = source
            self.load_simul(source)
            self.source = 'simul'
        elif 'pylayers' in source.__module__:
            self.filetraj = source._filename
            self.load_CorSer(source)
            cutoff=2
            self.source = 'CorSer'


        # generate the Network
        # The wireless standard and frequency is fixed in this function
        #
        self._gen_net()
        # initialize Stochastic Link
        self.SL = SLink()
        # initialize Deterministic Link
        self.DL = DLink(L=self.L,verbose=self.verbose)
        self.DL.cutoff=cutoff

        self.filename = 'simultraj_' + self.filetraj + '.h5'

        # data is a panda container which is initialized 
        #
        # We do not save all the simulation in a DataFRame anymore
        #
        #self.data = pd.DataFrame(columns=['id_a', 'id_b',
        #                                  'x_a', 'y_a', 'z_a',
        #                                  'x_b', 'y_b', 'z_b',
        #                                  'd', 'eng', 'typ',
        #                                  'wstd', 'fcghz',
        #                                  'fbminghz', 'fbmaxghz', 'fstep', 'aktk_id',
        #                                  'sig_id', 'ray_id', 'Ct_id', 'H_id'
        #                                  ])

        #self.data.index.name='t'
        self._filecsv = self.filename.split('.')[0] + '.csv'
        self.todo = {'OB': True,
                    'B2B': True,
                    'B2I': True,
                    'I2I': False}

        filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])

        if os.path.exists(filenameh5) :
            self.loadpd()
        self.settime(0.)
        # self._saveh5_init()


    def __repr__(self):


        s = 'Simul trajectories class\n'
        s = s + '------------------------\n'
        s = s +'\n'
        s = s + 'Used layout: ' + self.L.filename + '\n'
        s = s + 'Number of Agents: ' + str(self.Nag) + '\n'
        s = s + 'Number of Access Points: ' + str(self.Nap) + '\n'
        s = s + 'Link to be evaluated: ' + str(self.todo) + '\n'
        s = s + 'tmin: ' + str(self._tmin) + '\n'
        s = s + 'tmax: ' + str(self._tmax) + '\n'
        s = s +'\n'
        # network info
        s = s + 'self.N :\n'
        s = s + self.N.__repr__() + '\n'
        s = s + 'CURRENT TIME: ' + str(self.ctime) + '\n'




        return s

    def load_simul(self, source):
        """  load a simultraj configuration file

        Parameters
        ----------

        source : string
            name of simulation file to be loaded

        """
        self.filetraj = source
        if not os.path.isfile(source):
            raise AttributeError('Trajectory file'+source+'has not been found.\
             Please make sure you have run a simulnet simulation before runining simultraj.')

        # get the trajectory
        traj = tr.Trajectories()
        traj.loadh5(self.filetraj)

        # get the layout
        self.L = Layout(traj.Lfilename)

        # resample trajectory
        for ut, t in enumerate(traj):
            if t.typ == 'ag':
                person = Body(t.name + '.ini')
                tt = t.time()
                self.dpersons.update({t.name: person})
                self._tmin = tt[0]
                self._tmax = tt[-1]
                self.time = tt
            else:
                pos = np.array([t.x[0], t.y[0], t.z[0]])
                self.dap.update({t.ID: {'pos': pos,
                                        'ant': antenna.Antenna(),
                                        'name': t.name
                                        }
                                 })
        self.ctime = np.nan
        self.Nag = len(self.dpersons.keys())
        self.Nap = len(self.dap.keys())
        self.traj = traj

    def load_CorSer(self,source):
        """ load CorSer file for simulation

        Parameters
        ----------

        source :
            name of simulation file to be loaded

        """

        if isinstance(source.B,Body):
            B=[source.B]
        elif isinstance(source.B,list):
            B=source.B
        elif isinstance(source.B,dict):
            B=source.B.values()
        else:
            raise AttributeError('CorSer.B must be a list or a Body')

        self.L=source.L
        self.traj = tr.Trajectories()
        self.traj.Lfilename=self.L.filename

        for b in B:
            self.dpersons.update({b.name: b})
            self._tmin = b.time[0]
            self._tmax = b.time[-1]
            self.time = b.time
            self.traj.append(b.traj)

        for ap in source.din:
            techno,ID=ap.split(':')
            if techno == 'HKB':
                techno = 'hikob'
            if techno == 'TCR':
                techno = 'tcr'
            if techno == 'BS':
                techno = 'bespoon'


            self.dap.update({ap: {'pos': source.din[ap]['p'],
                                  'ant': source.din[ap]['ant'],
                                  'T': source.din[ap]['T'],
                                  'name': techno
                                        }
                                 })
        self.ctime = np.nan
        self.Nag = len(B)
        self.Nap = len(source.din)
        self.corser = source

    def _gen_net(self):
        """ generate Network and associated links

        Notes
        -----

        Create self.N : Network object

        See Also
        --------

        pylayers.network.network

        """

        #
        # Create Network
        #
        N = Network()
        #
        # get devices on bodies
        #
        # forall person
        #   forall device
        for p in self.dpersons:
            D = []
            for dev in self.dpersons[p].dev:
                aDev = Device(self.dpersons[p].dev[dev]['name'], ID = dev)
                D.append(aDev)

                D[-1].ant['A1']['name'] = self.dpersons[p].dev[dev]['file']
                D[-1].ant['antenna'] = self.dpersons[p].dev[dev]['ant']
            N.add_devices(D, grp=p)
        #
        # get access point devices
        #
        for ap in self.dap:
            D = Device(self.dap[ap]['name'], ID = ap)
            D.ant['antenna'] = self.dap[ap]['ant']
            N.add_devices(D, grp = 'ap', p = self.dap[ap]['pos'])
            N.update_orient(ap, self.dap[ap]['T'], now = 0.)
        # create Network
        #
        #    _get_wstd
        #    _get_grp
        #    _connect
        #    _init_PN
        #
        N.create()
        self.N = N

    def show(self):
        """ show actual simlulation configuration
        """
        fig, ax = self.L.showGs()
        fig, ax = self.N.show(fig=fig, ax=ax)
        return fig, ax

    def evaldeter(self, na, nb, wstd, fmod='force',nf=10,fGHz=[], **kwargs):
        """ deterministic evaluation of a link

        Parameters
        ----------

        na : string:
            node a id in self.N (Network)
        nb : string:
            node b id in self.N (Network)
        wstd : string:
            wireless standard used for commmunication between na and nb
        fmode : string ('center'|'band'|'force')
            mode of frequency evaluation
            center : single frequency (center frequency of a channel)
            band : nf points on the whole band
            force : takes directly fGHz
        nf : int:
            number of frequency points (if fmode = 'band')
        **kwargs : argument of DLink

        Returns
        -------

        (a, t )

        a : ndarray
            alpha_k
        t : ndarray
            tau_k

        See Also
        --------

        pylayers.simul.link.DLink

        """

        # todo in network :
        # take into consideration the postion and rotation of antenna and not device

        self.DL.Aa = self.N.node[na]['ant']['antenna']
        self.DL.a = self.N.node[na]['p']
        self.DL.Ta = self.N.node[na]['T']

        self.DL.Ab = self.N.node[nb]['ant']['antenna']
        self.DL.b = self.N.node[nb]['p']
        self.DL.Tb = self.N.node[nb]['T']

        #
        # The frequency band is chosen from the selected standard
        #  if fmode == 'center'
        #      only center frequency is calculated
        #
        #'
        if fmod == 'center':
            self.DL.fGHz = self.N.node[na]['wstd'][wstd]['fcghz']
        if fmod == 'band':
            fminGHz = self.N.node[na]['wstd'][wstd]['fbminghz']
            fmaxGHz = self.N.node[na]['wstd'][wstd]['fbmaxghz']
            self.DL.fGHz = np.linspace(fminGHz, fmaxGHz, nf)
        if fmod == 'force':
            assert len(fGHz)>0,"fGHz has not been defined"
            self.DL.fGHz = fGHz

        a, t = self.DL.eval(**kwargs)

        return a, t

    def evalstat(self, na, nb):
        """ statistical evaluation of a link

        Parameters
        ----------

        na : string:
            node a id in self.N (Netwrok)
        nb : string:
            node b id in self.N (Netwrok)

        Returns
        -------

        (a, t, eng)

        a : ndarray
            alpha_k
        t : ndarray
            tau_k
        eng : float
            engagement
        """

        pa = self.N.node[na]['p']
        pb = self.N.node[nb]['p']
        if self.source == 'simul':
            dida, name = na.split('_')
            didb, name = nb.split('_')
        elif self.source =='CorSer':
            bpa,absolutedida,dida,name,technoa = self.corser.devmapper(na)
            bpb,absolutedidb,didb,name,technob = self.corser.devmapper(nb)

        ak, tk, eng = self.SL.onbody(self.dpersons[name], dida, didb, pa, pb)

        return ak, tk, eng


    def settime(self,t):
        """ set current time
        """
        self.ctime = t
        self._traj=copy.copy(self.traj)
        self.update_pos(t)


    def run(self, **kwargs):
        """ run the link evaluation along a trajectory


        Parameters
        ----------

        OB: boolean
            perform on body statistical link evaluation
        B2B:  boolean
            perform body to body deterministic link evaluation
        B2I: boolean
            perform body to infrastructure deterministic link evaluation
        I2I:  boolean
            perform infrastructure to infrastructure deterministic link eval.
        links: dict
            dictionnary of link to be evaluated (key is wtsd and value is a list of links)
            (if [], all link are considered)
        wstd: list
            list of wstd to be evaluated
            (if [], all wstd are considered)
        t: np.array
            list of timestamp to be evaluated
            (if [], all timestamps are considered)
        tbr : boolean
            time in bit reverse order (tmin,tmax,N) Npoints=2**N
        replace_data: boolean (True)
            if True , reference id of all already simulated link will be erased
                and replace by new simulation id

        fGHz : np.array
            frequency in GHz


        Examples
        --------

            >>> from pylayers.simul.simultraj import *
            >>> from pylayers.measures.cormoran import *
            >>> C=CorSer()
            >>> S=Simul(C,verbose=True)
            >>> link={'ieee802154':[]}
            >>> link['ieee802154'].append(S.N.links['ieee802154'][0])
            >>> lt = [0,0.2,0.3,0.4,0.5]
            >>> S.run(links=link,t=lt)


        """
        defaults = {'OB': True,
                    'B2B': True,
                    'B2I': True,
                    'I2I': False,
                    'links': {},
                    'wstd': [],
                    't': np.array([]),
                    'btr':True,
                    'DLkwargs':{},
                    'replace_data':True,
                    'fmod':'force',
                    'fGHz':np.array([2.45])
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        DLkwargs = kwargs.pop('DLkwargs')
        links = kwargs.pop('links')
        wstd = kwargs.pop('wstd')
        OB = kwargs.pop('OB')
        B2B = kwargs.pop('B2B')
        B2I = kwargs.pop('B2I')
        I2I = kwargs.pop('I2I')
        fmod = kwargs.pop('fmod')
        self.fGHz = kwargs.pop('fGHz')

        self.todo.update({'OB':OB,'B2B':B2B,'B2I':B2I,'I2I':I2I})


        # Check link attribute

        if links == {}:
            links = self.N.links
        elif not isinstance(links, dict):
            raise AttributeError('links is {wstd:[list of links]}, see self.N.links')

        for k in links.keys():
            checkl = [l in self.N.links[k] for l in links[k]]
            if len(np.where(checkl==False)[0])>0:
            # if sum(checkl) != len(self.N.links):
                uwrong = np.where(np.array(checkl) is False)[0]
                raise AttributeError(str(np.array(links)[uwrong])
                                     + ' links does not exist in Network')

        wstd = links.keys()
        # # Check wstd attribute
        # if wstd == []:
        #     wstd = self.N.wstd.keys()
        # elif not isinstance(wstd, list):
        #     wstd = [wstd]

        checkw = [w in self.N.wstd.keys() for w in wstd]
        if sum(checkw) != len(wstd):
            uwrong = np.where(np.array(checkw) is False)[0]
            raise AttributeError(str(np.array(wstd)[uwrong])
                                 + ' wstd are not in Network')

        # force time attribute compliant

        if not isinstance(kwargs['t'],np.ndarray):
            if isinstance(kwargs['t'],list):
                lt = np.array(kwargs['t'])
            elif (isinstance(kwargs['t'], int)
                 or isinstance(kwargs['t'],float)):
                lt = np.array([kwargs['t']])
        else :
            lt = kwargs['t']

        #if len(lt) == 0:
        #    lt = self.time
        # check time attribute
        if kwargs['btr']:
            if (lt[0] < self._tmin) or\
               (lt[1] > self._tmax) :
                raise AttributeError('Requested time range not available')

        # self._traj is a copy of self.traj, which is affected by resampling.
        # it is only a temporary attribute for a given run
        # if len(lt) > 1:
        #     sf = 1/(1.*lt[1]-lt[0])
        #     self._traj = self.traj.resample(sf=sf, tstart=lt[0])

        # else:
        #     self._traj = self.traj.resample(sf=1.0, tstart=lt[0])
        #     self._traj.time()
        # self.time = self._traj.t
        # self._time = pd.to_datetime(self.time,unit='s')
        #
        # Nested Loops
        #
        #  time
        #    standard
        #      links
        #           evaldeter &| evalstat
        #
        #lt = self.get_sim_time(lt)
        #self._time=self.get_sim_time(lt)

        init = True
        if kwargs['btr']:
            tmin = lt[0]
            tmax = lt[1]
            Nt   = int(2**lt[2])
            ta   = np.linspace(tmin,tmax,Nt)
            it   = np.hstack((np.r_[0],np.r_[pyu.bitreverse(Nt,int(lt[2]))]))
            #trev = t[it]
        else:
            ta = kwargs['t']
            it = range(len(ta))

        ## Start to loop over time
        ##   ut : counter
        ##   t  : time value (s)
        #for ut, t in enumerate(lt):
        for ks,ut in enumerate(it):
            t  = ta[ut]
            self.ctime = t
            # update spatial configuration of the scene for time t
            self.update_pos(t)
            # print self.N.__repr__()
            ## Start to loop over available Wireless standard
            ##
            for w in wstd:
                ## Start to loop over the chosen links stored in links
                ##
                for na, nb, typ in links[w]:
                    # If type of link is valid (Body 2 Body,...)
                    #
                    if self.todo[typ]:
                        if self.verbose:
                            print '-'*30
                            print 'time:', t, '/',  lt[-1] ,' time idx:', ut,
                            '/',len(ta),'/',ks
                            print 'processing: ',na, ' <-> ', nb, 'wstd: ', w
                            print '-'*30
                        eng = 0
                        #
                        # Invoque link deterministic simulation 
                        #
                        #  node : na
                        #  node : nb
                        #  wstd : w
                        #
                        self.evaldeter(na, nb,
                                       w,
                                       applywav=False,
                                       fmod = fmod,
                                       fGHz = self.fGHz,
                                       **DLkwargs)
                        # if typ == 'OB':
                        #     self.evalstat(na, nb)
                        #     eng = self.SL.eng
                        #     L = self.DL + self.SL
                        #     self._ak = L.H.ak
                        #     self._tk = L.H.tk
                        # else :

                        # Get alphak an tauk
                        self._ak = self.DL.H.ak
                        self._tk = self.DL.H.tk
                        aktk_id = str(ut) + '_' + na + '_' + nb + '_' + w
                        # this is a dangerous way to proceed ! 
                        # the id as a finite number of characters
                        while len(aktk_id)<40:
                            aktk_id = aktk_id + ' '
                        df = pd.DataFrame({ 'id_a': na,
                                    'id_b': nb,
                                    'x_a': self.N.node[na]['p'][0],
                                    'y_a': self.N.node[na]['p'][1],
                                    'z_a': self.N.node[na]['p'][2],
                                    'x_b': self.N.node[nb]['p'][0],
                                    'y_b': self.N.node[nb]['p'][1],
                                    'z_b': self.N.node[nb]['p'][2],
                                    'd': self.N.edge[na][nb]['d'],
                                    'eng': eng,
                                    'typ': typ,
                                    'wstd': w,
                                    'fcghz': self.N.node[na]['wstd'][w]['fcghz'],
                                    'fbminghz': self.fGHz[0],
                                    'fbmaxghz': self.fGHz[-1],
                                    'nf': len(self.fGHz),
                                    'aktk_id':aktk_id,
                                    'sig_id': self.DL.dexist['sig']['grpname'],
                                    'ray_id': self.DL.dexist['ray']['grpname'],
                                    'Ct_id': self.DL.dexist['Ct']['grpname'],
                                    'H_id': self.DL.dexist['H']['grpname'],
                                                },columns=['id_a', 'id_b',
                                              'x_a', 'y_a', 'z_a',
                                              'x_b', 'y_b', 'z_b',
                                              'd', 'eng', 'typ',
                                              'wstd', 'fcghz',
                                              'fbminghz', 'fbmaxghz', 'fstep', 'aktk_id',
                                              'sig_id', 'ray_id', 'Ct_id', 'H_id'
                                              ],index= [t])  #self._time[ut]])

                        self.savepd(df)

    def replace_data(self, df):
        """check if a dataframe df already exists in self.data

        Parameters
        ----------
        df : pd.DataFrame

        Returns
        -------

        boolean
            True if already exists
            False otherwise

        """
        self.data[(self.data.index == df.index) &
                  (self.data['id_a'] == df['id_a'].values[0]) &
                  (self.data['id_b'] == df['id_b'].values[0]) &
                  (self.data['wstd'] == df['wstd'].values[0])]=df.values



    def check_exist(self, df):
        """check if a dataframe df already exists in self.data

        Parameters
        ----------
        df : pd.DataFrame

        Returns
        -------

        boolean
            True if already exists
            False otherwise

        """
        # check init case 
        if not len(self.data.index) == 0:

            ud = self.data[(self.data.index == df.index) & 
                           (self.data['id_a'] == df['id_a'].values[0]) & 
                           (self.data['id_b'] == df['id_b'].values[0]) & 
                           (self.data['wstd'] == df['wstd'].values[0])]

            if len(ud) == 0:
                return False
            else :
                return True
        else :
            return False


    def savepd(self,df):
        """ save data information of a simulation

        Parameters
        ----------

        df : one index data

        Notes
        -----


        """
        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        store = pd.HDFStore(filenameh5)
        #self.data=self.data.sort()
        store.append('df',df)
        store.close()

    def loadpd(self):
        """ load data from previous simulations
        """
        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        store = pd.HDFStore(filenameh5)
        #self.data = pd.read_hdf(filenameh5,'df')
        self.data = store.get('df')
        self.data.index.name='t'
        self.data = self.data.sort()

    def get_sim_time(self,t):
        """ retrieve closest time value in regard of passed t value in parameter
        """

        if not isinstance(t,list) and not isinstance(t,np.ndarray):
            return np.array([self.time[np.where(self.time <=t)[0][-1]]])
        else :
            return np.array([self.get_sim_time(tt) for tt in t])[:,0]


    def get_df_from_link(self,id_a,id_b,wstd=''):
        """ Return a restricted data frame for a specific link

            Parameters
            ----------

            id_a : str
                node id a
            id_b: str
                node id b
            wstd: str
                optionnal :wireslees standard
        """
        if wstd == '':
            return self.data[(self.data['id_a']==id_a) &
                             (self.data['id_b']==id_b)]
        else :
            return self.data[(self.data['id_a']==id_a) &
                             (self.data['id_b']==id_b) &
                             self.data['wstd']==wstd]


    def update_pos(self, t):
        """ update positions of devices and bodies for a given time index

        Parameters
        ----------

        t : int
            time value

        """

        # if a bodies are involved in simulation
        if ((self.todo['OB']) or (self.todo['B2B']) or (self.todo['B2I'])):
            nodeid = []
            pos = []
            devlist = []
            orient = []
            for up, person in enumerate(self.dpersons.values()):
                person.settopos(self._traj[up], t=t, cs=True)
                name = person.name
                dev = person.dev.keys()
                devlist.extend(dev)
                #nodeid.extend([n + '_' + name for n in dev])
                pos.extend([person.dcs[d][:, 0] for d in dev])
                orient.extend([person.acs[d] for d in dev])
            # TODO !!!!!!!!!!!!!!!!!!!!
            # in a future version , the network update must also update
            # antenna position in the device coordinate system
            self.N.update_pos(devlist, pos, now=t)
            self.N.update_orient(devlist, orient, now=t)
        self.N.update_dis()



    def get_value(self,**kwargs):
        """ retrieve output parameter at a specific time

        Parameters
        ----------

        typ : list
                list of parameters to be retrieved
                (R | C |H | ak | tk | rss )
        links: list
            dictionnary of link to be evaluated (key is wtsd and value is a list of links)
            (if [], all link are considered)
        t: int or np.array
            list of timestamp to be evaluated | singlr time instant

        Returns
        -------

        output: dict
                [link_key]['t']
                          ['ak']
                ...
        """



        # get time
        defaults = {'t': 0,
                    'typ':['ak'],
                    'links': {},
                    'wstd':[],
                    'angles':False
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        # allocate an empty dictionnary for wanted selected output
        output={}

        # manage time t can be a list or a float
        t = kwargs['t']
        t = self.get_sim_time(t)
        dt = self.time[1]-self.time[0]

        # manage links
        plinks = kwargs['links']
        links=[]
        if isinstance(plinks,dict):
            for l in plinks.keys():
                links.extend(plinks[l])

        if len(links) == 0:
            raise AttributeError('Please give valid links to get values')
        # output['t']=[]
        # output['time_to_simul']=[]
        # for each requested time step
        for tt in t :
            # for each requested links
            for link in links:
                linkname=link[0]+'-'+link[1]
                if not output.has_key(linkname):
                    output[linkname] = {}
                if not output[linkname].has_key('t'):
                    output[linkname]['t'] = []

                # restrict global dataframe self.data to the specific link
                df = self.get_df_from_link(link[0],link[1])
                # restrict global dataframe self.data to the specific z
                df = df[(df.index > tt-dt) & (df.index <= tt+dt)]

                if len(df) != 0:
                    output[linkname]['t'].append(tt)
                    if len(df)>1:
                        print 'Warning possible issue in self.get_value'
                    line = df.iloc[-1]
                    # # get info of the corresponding timestamp
                    # line = df[(df['id_a'] == link[0]) & (df['id_b'] == link[1])].iloc[-1]
                    # if len(line) == 0:
                    #     line = df[(df['id_b'] == link[0]) & (df['id_a'] == link[1])]
                    #     if len(line) == 0:
                    #         raise AttributeError('invalid link')

                    #retrieve correct position and orientation given the time
                    #self.update_pos(t=tt)
                    # antennas positions
                    #self.DL.a = self.N.node[link[0]]['p']
                    #self.DL.b = self.N.node[link[1]]['p']
                    # antennas orientation
                    #self.DL.Ta = self.N.node[link[0]]['T']
                    #self.DL.Tb = self.N.node[link[1]]['T']
                    # antennas object
                    #self.DL.Aa = self.N.node[link[0]]['ant']['antenna']
                    #self.DL.Ab = self.N.node[link[1]]['ant']['antenna']
                    # get the antenna index
                    #uAa_opt, uAa = self.DL.get_idx('A_map',self.DL.Aa._filename)
                    #uAb_opt, uAb = self.DL.get_idx('A_map',self.DL.Ab._filename)

                    if 'ak' in kwargs['typ'] or 'tk' in kwargs['typ'] or 'rss' in kwargs['typ']:
                        H_id = line['H_id'].decode('utf8')
                        # load the proper link
                        # parse index
                        lid = H_id.split('_')
                        #if (lid[5]==str(uAa))&(lid[6]==str(uAb)):
                        self.DL.load(self.DL.H,H_id)
                        if 'ak' in kwargs['typ']:
                            if not output[linkname].has_key('ak'):
                                output[linkname]['ak']=[]
                            output[linkname]['ak'].append(copy.deepcopy(self.DL.H.ak))
                        if 'tk' in kwargs['typ']:
                            if not output[linkname].has_key('tk'):
                                output[linkname]['tk']=[]
                            output[linkname]['tk'].append(copy.deepcopy(self.DL.H.tk))
                        if 'rss' in kwargs['typ']:
                            if not output[linkname].has_key('rss'):
                                output[linkname]['rss']=[]
                            output[linkname]['rss'].append(copy.deepcopy(self.DL.H.rssi()))

                    if 'R' in kwargs['typ']:
                        if not output[linkname].has_key('R'):
                            output[linkname]['R']=[]
                        ray_id = line['ray_id']
                        self.DL.load(self.DL.R,ray_id)
                        output[linkname]['R'].append(copy.deepcopy(self.DL.R))

                    if 'C' in kwargs['typ']:
                        if not output[linkname].has_key('C'):
                            output[linkname]['C']=[]
                        Ct_id = line['Ct_id']
                        self.DL.load(self.DL.C,Ct_id)

                        if kwargs['angles']:
                            self.DL.C.islocal=False
                            self.DL.C.locbas(Tt=self.DL.Ta, Tr=self.DL.Tb)
                        #T channel
                        output[linkname]['C'].append(copy.deepcopy(self.DL.C))


                    if 'H' in kwargs['typ']:
                        if not output[linkname].has_key('H'):
                            output[linkname]['H']=[]
                        H_id = line['H_id']
                        lid = H_id.split('_')
                        #if (lid[5]==str(uAa))&(lid[6]==str(uAb)):
                        self.DL.load(self.DL.H,H_id)
                        output[linkname]['H'].append(copy.deepcopy(self.DL.H))

                # if time value not found in dataframe
                else:
                    if not output[linkname].has_key('time_to_simul'):
                        output[linkname]['time_to_simul'] = []
                    output[linkname]['time_to_simul'].append(tt)


        for l in output.keys():
            if output[l].has_key('time_to_simul'):
                print 'link', l , 'require simulation for timestamps', output[l]['time_to_simul']


        return(output)


    def get_link(self,**kwargs):
        """ retrieve a Link specific time from a simultraj

        Parameters
        ----------

        typ : list
                list of parameters to be retrieved
                (ak | tk | R |C)
        links: list
            dictionnary of link to be evaluated (key is wtsd and value is a list of links)
            (if [], all link are considered)
        t: int or np.array
            list of timestamp to be evaluated | singlr time instant

        Returns
        -------

        DL : DLink

        Examples
        --------

        >>> from pylayers.simul.simultraj import *
        >>> from pylayers.measures.cormoran import *
        >>> C=CorSer(serie=6i,day=11)
        >>> S = Simul(C,verb ose=False)
        >>> DL = S.get_link(typ=['R','C','H'])
                ...
        """



        # get time
        defaults = {'t': 0,
                    'typ':['ak'],
                    'links': {},
                    'wstd':[],
                    'angles':False
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        output={}

        # manage time
        t = kwargs['t']
        t = self.get_sim_time(t)
        dt = self.time[1]-self.time[0]

        # manage links
        plinks = kwargs['links']
        links=[]
        if isinstance(plinks,dict):
            for l in plinks.keys():
                links.extend(plinks[l])

        if len(links) == 0:
            raise AttributeError('Please give valid links to get values')
        # output['t']=[]
        # output['time_to_simul']=[]
        # for each requested time step
        for tt in t :
            # for each requested links
            for link in links:
                linkname=link[0]+'-'+link[1]
                if not output.has_key(linkname):
                    output[linkname] = {}
                if not output[linkname].has_key('t'):
                    output[linkname]['t'] = []


                # restrict global dataframe self.data to the specific link
                df = self.get_df_from_link(link[0],link[1])
                # restrict global dataframe self.data to the specific z
                df = df[(df.index > tt-dt) & (df.index <= tt+dt)]

                if len(df) != 0:
                    output[linkname]['t'].append(tt)
                    if len(df)>1:
                        print 'Warning possible issue in self.get_link'
                    line = df.iloc[-1]
                    # # get info of the corresponding timestamp
                    # line = df[(df['id_a'] == link[0]) & (df['id_b'] == link[1])].iloc[-1]
                    # if len(line) == 0:
                    #     line = df[(df['id_b'] == link[0]) & (df['id_a'] == link[1])]
                    #     if len(line) == 0:
                    #         raise AttributeError('invalid link')

                    #retrieve correct position and orientation given the time
                    self.update_pos(t=tt)
                    self.DL.a = self.N.node[link[0]]['p']
                    self.DL.b = self.N.node[link[1]]['p']
                    self.DL.Ta = self.N.node[link[0]]['T']
                    self.DL.Tb = self.N.node[link[1]]['T']
                    #self.DL.Aa = self.N.node[link[0]]['ant']['antenna']
                    #self.DL.Ab = self.N.node[link[1]]['ant']['antenna']

                    #H_id = line['H_id'].decode('utf8')
                    #self.DL.load(self.DL.H,H_id)

                    if 'R' in kwargs['typ']:
                        ray_id = line['ray_id']
                        self.DL.load(self.DL.R,ray_id)

                    if 'C' in kwargs['typ']:
                        Ct_id = line['Ct_id']
                        self.DL.load(self.DL.C,Ct_id)

                        if kwargs['angles']:
                            self.DL.C.islocal=False
                            self.DL.C.locbas(Tt=self.DL.Ta, Tr=self.DL.Tb)

                    if 'H' in kwargs['typ']:
                        H_id = line['H_id']
                        self.DL.load(self.DL.H,H_id)

        return(self.DL)



    def _show3(self, **kwargs):
        """ 3D show using Mayavi

        Parameters
        ----------

        t: float
            time index
        link: list
            [id_a, id_b]
            id_a : node id a
            id_b : node id b
        'lay': bool
            show layout
        'net': bool
            show net
        'body': bool
            show bodies
        'rays': bool
            show rays
        """

        defaults = {'t': 0,
                    'link': [],
                    'wstd':[],
                    'lay': True,
                    'net': True,
                    'body': True,
                    'rays': True,
                    'ant': False

                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        link = kwargs['link']
        self.update_pos(kwargs['t'])

        if len(self.data) != 0:
            df = self.data[self.data.index == pd.to_datetime(kwargs['t'])]
            if len(df) != 0:
                raise AttributeError('invalid time')


            # default
            if link ==[]:
                line = df[df.index<=pd.to_datetime(0)]
                link = [line['id_a'].values[0],line['id_b'].values[0]]
            else :
                # get info of the corresponding timestamp
                line = df[(df['id_a'] == link[0]) & (df['id_b'] == link[1])]
            if len(line) == 0:
                line = df[(df['id_b'] == link[0]) & (df['id_a'] == link[1])]
                if len(line) == 0:
                    raise AttributeError('invalid link')
            rayid = line['ray_id'].values[0]


            self.DL.a = self.N.node[link[0]]['p']
            self.DL.b = self.N.node[link[1]]['p']
            self.DL.Ta = self.N.node[link[0]]['T']
            self.DL.Tb = self.N.node[link[1]]['T']
            self.DL.load(self.DL.R,rayid)





            self.DL._show3(newfig= False,
                           lay= kwargs['lay'],
                           rays= kwargs['rays'],
                           ant=False)
        else :
            self.DL._show3(newfig= False,
                           lay= True,
                           rays= False,
                           ant=False)
        if kwargs['net']:
            self.N._show3(newfig=False)
        if kwargs['body']:
            for p in self.dpersons:
                self.dpersons[p]._show3(newfig=False,
                                        topos=True,
                                        pattern=kwargs['ant'])

    # def _saveh5_init(self):
    #     """ initialization of the h5py file
    #     """
    #     filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
    #     import ipdb
    #     try:
    #         f5 = h5py.File(filenameh5, 'w')
    #         f5.create_dataset('time', shape=self.time.shape, data=self.time)
    #         f5.close()
    #     except:
    #         f5.close()
    #         raise NameError('simultra.saveinit: \
    #                         issue when writting h5py file')

    def _saveh5(self, ut, ida, idb, wstd):
        """ Save in h5py format

        Parameters
        ----------

        ut : int
            time index in self.time
        ida : string
            node a index
        idb : string
            node b index
        wstd : string
            wireless standard of used link

        Notes
        -----

        Dataset organisation:

        simultraj_<trajectory_filename.h5>.h5
            |
            |time
            |    ...
            |
            |/<tidx_ida_idb_wstd>/ |attrs
            |                      |a_k
            |                      |t_k


        Root dataset :
        time : array
            range of simulation time

        Group identifier :
            tidx : index in time dataset
            ida : node a index in Network
            idb : node b index in Network
            wstd : wireless standar of link interest


        Inside group:
            a_k : alpha_k values
            t_k : tau_k values

        See Also
        --------

        pylayers.simul.links

        """

        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        grpname = str(ut) + '_' + ida + '_' + idb + '_' + wstd
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            fh5 = h5py.File(filenameh5, 'a')
            if not grpname in fh5.keys():
                fh5.create_group(grpname)
                f = fh5[grpname]
                # for k in kwargs:
                #     f.attrs[k] = kwargs[k]

                f.create_dataset('alphak',
                                 shape=self._ak.shape,
                                 maxshape=(None),
                                 data=self._ak)
                f.create_dataset('tauk',
                                 shape=self._tk.shape,
                                 maxshape=(None),
                                 data=self._tk)
            else:
                pass#print grpname + ' already exists in ' + filenameh5


            fh5.close()
        except:
            fh5.close()
            raise NameError('Simultraj._saveh5: issue when writting h5py file')


    def _loadh5(self, grpname):
        """ Load in h5py format

        Parameters
        ----------

       grpname : string
            group name which can be found sin self.data aktk_idx column

        Returns
        -------
        (ak, tk, conf)

        ak : ndarray:
            alpha_k
        tk : ndarray:
            alpha_k
        """

        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            fh5 = h5py.File(filenameh5, 'r')
            if not grpname in fh5.keys():
                fh5.close()
                raise NameError(grpname + ' cannot be reached in ' + self.filename)
            f = fh5[grpname]
            # for k in f.attrs.keys():
            #     conf[k]=f.attrs[k]
            ak = f['alphak'][:]
            tk = f['tauk'][:]
            fh5.close()

            return ak, tk
        except:
            fh5.close()
            raise NameError('Simultraj._loadh5: issue when reading h5py file')


    def tocsv(self, ut, ida, idb, wstd,init=False):

        filecsv = pyu.getlong(self._filecsv,pstruc['DIRLNK'])

        with open(filecsv, 'a') as csvfile:
            fil = csv.writer(csvfile, delimiter=';',
                             quoting=csv.QUOTE_MINIMAL)
            if init:
                keys = self.data.iloc[-1].keys()
                data = [k for k in keys]
                data .append('ak')
                data .append('tk')
                fil.writerow(data)

            values = self.data.iloc[-1].values
            data = [v for v in values]
            sak = str(self._ak.tolist())
            stk = str(self._tk.tolist())
            data.append(sak)
            data.append(stk)
            fil.writerow(data)

if (__name__ == "__main__"):
    #plt.ion()
    doctest.testmod()
