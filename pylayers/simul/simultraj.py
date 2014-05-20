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
    Simul.gen_net
    Simul.evaldeter
    Simul.evalstat
    Simul.show

Internal configuration
----------------------

.. autosummary::
    :toctree: generated/

    Simul._saveh5
    Simul._loadh5

    
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
# Handle directory hierarchy
from pylayers.util.project import *
# Handle UWB measurements
import pylayers.mobility.trajectory as tr
from pylayers.mobility.ban.body import *
from pylayers.antprop.statModel import *
import pandas as pd
import csv

class Simul(object):

    """
    Link oriented simulation

    A simulation requires :

        A Layout
        A Person
        A Trajectory

    """

    def __init__(self, _filetraj='simulnet_TA-Office.h5',verbose=False):

        self.filetraj = _filetraj

        # self.progress = -1  # simulation not loaded
        self.verbose = verbose
        self.cfield = []
        self.dpersons = {}

        self.dap = {}
        self.Nag = 0
        self.Nap = 0
        self.load_config(_filetraj)
        self.gen_net()
        self.SL = SLink()
        self.DL = DLink(L=self.L,verbose=self.verbose)
        self.filename = 'simultraj_' + self.filetraj
        self.data = pd.DataFrame(columns=['id_a', 'id_b',
                                          'x_a', 'y_a', 'z_a',
                                          'x_b', 'y_b', 'z_b',
                                          'd', 'eng', 'typ',
                                          'wstd', 'fcghz',
                                          'fbminghz', 'fbmaxghz', 'fstep', 'aktk_id',
                                          'sig_id', 'ray_id', 'Ct_id', 'H_id'
                                          ])
        self.data.index.name='t'
        self._filecsv = self.filename.split('.')[0] + '.csv'
        self.todo = {'OB': True,
                    'B2B': True,
                    'B2I': True,
                    'I2I': False}
        filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])
        if os.path.exists(filenameh5) :
            self.loadpd()
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

    def load_config(self, _filetraj):
        """  load a simultraj configuration file

        Parameters
        ----------

        _filetraj : string
            name of simulation file to be loaded

        """
        self.filetraj = _filetraj


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


    def gen_net(self):
        """ generate Network and associated links

        Notes
        -----

        Create self.N : Network object

        See Also
        --------

        pylayers.network.network


        """

        N = Network()
        # get devices on bodies
        for p in self.dpersons:
            D = []
            for dev in self.dpersons[p].dev:
                D.append(
                    Device(self.dpersons[p].dev[dev]['name'], ID=dev + '_' + p))
            N.add_devices(D, grp=p)
        # get access point devices
        for ap in self.dap:
            D = Device(self.dap[ap]['name'], ID=ap)
            N.add_devices(D, grp='ap', p=self.dap[ap]['pos'])
        # create Network
        N.create()
        self.N = N

    def show(self):
        """ show actual simlulation configuration
        """
        fig, ax = self.L.showGs()
        fig, ax = self.N.show(fig=fig, ax=ax)
        return fig, ax

    def evaldeter(self, na, nb, wstd, fmode='band', nf=10):
        """ Deterministic evaluation of a link

        Parameters
        ----------

        na : string:
            node a id in self.N (Network)
        nb : string:
            node b id in self.N (Network)
        wstd : string:
            wireless standard used for commmunication between na and nb
        fmode : string ('center'|'band')
            mode of frequency evaluation
            center : only on the centered frequency
            band : on the whole band
        nf : int:
            number of frequency points (if fmode = 'band')

        Returns
        -------

        (a, t )

        a : ndarray
            alpha_k
        t : ndarray
            tau_k

        """
        self.DL.a = self.N.node[na]['p']
        self.DL.Ta = self.N.node[na]['T']
        self.DL.b = self.N.node[nb]['p']
        self.DL.Tb = self.N.node[nb]['T']
        if fmode == 'center':
            self.DL.fGHz = self.N.node[na]['wstd'][wstd]['fcghz']
        else:
            minb = self.N.node[na]['wstd'][wstd]['fbminghz']
            maxb = self.N.node[na]['wstd'][wstd]['fbmaxghz']
            self.DL.fGHz = np.linspace(minb, maxb, nf)
        a, t = self.DL.eval()

        return a, t

    def evalstat(self, na, nb):
        """ Statistical evaluation of a link

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
        dida, name = na.split('_')
        didb, name = nb.split('_')

        ak, tk, eng = self.SL.onbody(self.dpersons[name], dida, didb, pa, pb)

        return ak, tk, eng

    def run(self, **kwargs):
        """ Run the link evaluation along a trajectory


        Parameters
        ----------

        'OB': boolean
            perform on body statistical link evaluation
        'B2B':  boolean
            perform body to body deterministic link evaluation
        'B2I': boolean
            perform body to infrastructure deterministic link evaluation
        'I2I':  boolean
            perform infrastructure to infrastructure deterministic link eval.
        'llink': list
            list of link to be evaluated
            (if [], all link are considered)
        'wstd': list
            list of wstd to be evaluated
            (if [], all wstd are considered)
        't': np.array
            list of timestamp to be evaluated
            (if [], all timestamps are considered)


        """
        defaults = {'OB': True,
                    'B2B': True,
                    'B2I': True,
                    'I2I': False,
                    'llink': [],
                    'wstd': [],
                    't': np.array([]),
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        llink = kwargs.pop('llink')
        wstd = kwargs.pop('wstd')
        OB = kwargs.pop('OB')
        B2B = kwargs.pop('B2B')
        B2I = kwargs.pop('B2I')
        I2I = kwargs.pop('I2I')
        self.todo.update({'OB':OB,'B2B':B2B,'B2I':B2I,'I2I':I2I})

        # Check link attribute
        if llink == []:
            llink = self.N.links
        elif not isinstance(llink, list):
            llink = [llink]

        checkl = [l in self.N.links for l in llink]
        if sum(checkl) != len(self.N.links):
            uwrong = np.where(np.array(checkl) is False)[0]
            raise AttributeError(str(np.array(llink)[uwrong])
                                 + ' links does not exist in Network')

        # Check wstd attribute
        if wstd == []:
            wstd = self.N.wstd.keys()
        elif not isinstance(wstd, list):
            wstd = [wstd]

        checkw = [w in self.N.wstd.keys() for w in wstd]
        if sum(checkw) != len(self.N.wstd.keys()):
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
        
        if len(lt) == 0:
            lt = self.time
        # check time attribute
        if not lt[0] >= self._tmin and\
               lt[-1] <= self._tmax:
               raise AttributeError('Requested time range not available')

        # self._traj is a copy of self.traj, which is affected by resampling.
        # it is only a temporary attribute for a given run

        if len(lt) > 1:
            sf = 1/(1.*lt[1]-lt[0])
            self._traj = self.traj.resample(sf=sf, tstart=lt[0])

        else:

            self._traj = self.traj.resample(sf=1.0, tstart=lt[0])

            self._traj.time()
        
        self.time = self._traj.t
        self._time = pd.to_datetime(self.time,unit='s')

        #
        # Code
        #
       
        init = True
        for ut, t in enumerate(lt):
            self.ctime = t            
            self.update_pos(t)
            #print self.N.__repr__()
            print 'ut = ', ut, 't = ', t
            for w in wstd:                
                for na, nb, typ in llink[w]:
                    
                    if self.todo[typ]:
                        
                        if self.verbose:
                            print '-'*30
                            print 'time:', t, '/',  lt[-1] ,' time idx:', ut, '/',len(lt)
                            print 'processing: ',na, ' <-> ', nb, 'wstd: ', w
                            print '-'*30
                        eng = 0
                        self.evaldeter(na, nb, w)
                        if typ == 'OB':
                            #~ self.evalstat(na, nb)
                            #~ eng = self.SL.eng
                            L = self.DL #+ self.SL
                            self._ak = L.H.ak
                            self._tk = L.H.tk
                        else :
                            self._ak = self.DL.H.ak
                            self._tk = self.DL.H.tk
                        df = pd.DataFrame({\
                                    'id_a': na,
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
                                    'fbminghz': self.DL.fmin,
                                    'fbmaxghz': self.DL.fmax,
                                    'fstep': self.DL.fstep,
                                    'aktk_id': str(ut) + '_' + na + '_' + nb + '_' + w,
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
                                              ],index=[self._time[ut]])
                        #if not self.check_exist(df):
                        self.data = self.data.append(df)
                        # self._index = self._index + 1
                        # save csv
                        self.tocsv(ut, na, nb, w,init=init)
                        init=False

                        # save pandas self.data
                        #self.savepd()
                        # save ak tauk
                        self._saveh5(ut, na, nb, w)


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

        ud = self.data[(self.data.index == df.index) 
                        & (self.data['id_a'] == df['id_a'].values[0])
                        & (self.data['id_b'] == df['id_b'].values[0])
                        & (self.data['wstd'] == df['wstd'].values[0])
                        ]
        if len(ud) == 0:
            return False
        else :
            return True


    def savepd(self):
        """ save data information of a simulation
        """
        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        store = pd.HDFStore(filenameh5,'a')
        self.data=self.data.sort()
        store['df'] = self.data
        store.close()

    def loadpd(self):
        """ load data from previous simulations
        """
        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        self.data = pd.read_hdf(filenameh5,'df')

    def update_pos(self, t):
        ''' update positions of devices and bodies for a given time index

        Parameters
        ----------
        ut : int
            time value 
        '''

        # if a bodies are involved in simulation
        if ((self.todo['OB']) or (self.todo['B2B']) or (self.todo['B2I'])):
            nodeid = []
            pos = []
            orient = []
            for up, person in enumerate(self.dpersons.values()):   
                #person.settopos(self._traj[up], t=t, cs=True)
                person.settopos(self._traj[up], t=t, cs=True,treadmill = True, p0 = np.array([1.5,4.5]))
                name = person.name
                dev = person.dev.keys()
                nodeid.extend([n + '_' + name for n in dev])
                pos.extend([person.dcs[d][:, 0] for d in dev])
                orient.extend([person.acs[d] for d in dev])
            # in a future version , the network update must also update
            # antenna positon in the device coordinate system
            self.N.update_pos(nodeid, pos, now=t)
            self.N.update_orient(nodeid, orient, now=t)
        # TODO : to be moved on the network edges
        self.N.update_dis()

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

        df = self.data[self.data['t'] == kwargs['t']]
        if len(df) == 0:
            raise AttributeError('invalid time')

        # default
        if link ==[]:
            line = df[df.index==1]
            link = [line['id_a'].values[0],line['id_b'].values[0]]
        else :
            # get info of the corresponding timestamp
            line = df[(df['id_a'] == link[0]) & (df['id_b'] == link[1])]
        if len(line) == 0:
            line = df[(df['id_b'] == link[0]) & (df['id_a'] == link[1])]
            if len(line) == 0:
                raise AttributeError('invalid link')
        rayid = line['ray_id'].values[0]


        self.update_pos(kwargs['t'])
        self.DL.a = self.N.node[link[0]]['p']
        self.DL.b = self.N.node[link[1]]['p']
        self.DL.Ta = self.N.node[link[0]]['T']
        self.DL.Tb = self.N.node[link[1]]['T']
        self.DL.load(self.DL.R,rayid)

        self.DL._show3(newfig= False,
                       lay= kwargs['lay'],
                       rays= kwargs['rays'],
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

