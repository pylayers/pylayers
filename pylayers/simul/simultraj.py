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


    Simul.load_config
    Simul.load_config
    Simul._saveh5_init
    Simul._saveh5

"""
import doctest
import ConfigParser
import numpy as np
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

    def __init__(self, _filesimul='simultraj.ini',verbose=False):

        self.filesimul = _filesimul
        self.config = ConfigParser.ConfigParser()
        self.config.add_section("layout")
        self.config.add_section("person")
        self.config.add_section("trajectory")

        self.progress = -1  # simulation not loaded

        self.cfield = []
        fmin = 3
        fmax = 6
        fstep = 0.05
        nf = (fmax - fmin) / fstep
        self.fGHz = np.linspace(fmin, fmax, nf, endpoint=True)
        self.wav = wvf.Waveform()

        self.dpersons = {}
        self.dap = {}
        self.Nag = 0
        self.Nap = 0
        self.load_config(_filesimul)
        self.gen_net()
        self.SL = SLink()
        self.DL = DLink(L=self.L,verbose=verbose)
        self.filename = 'simultraj_' + self._trajname
        self.data = pd.DataFrame(columns=['t', 'id_a', 'id_b',
                                          'x_a', 'y_a', 'z_a',
                                          'x_b', 'y_b', 'z_b',
                                          'd', 'eng', 'typ',
                                          'wstd', 'fcghz',
                                          'fbminghz', 'fbmaxghz', 'fstep',
                                          'sig_id', 'ray_id', 'Ct_id', 'H_id'
                                          ])
        self.data.set_index('t')
        self._filecsv = self.filename.split('.')[0] + '.csv'
        # self._saveh5_init()


    def __repr__(self):

        try:
            s = 'Simul trajectories class\n'
            s = s + '------------------------\n\n'
            s = s + 'Used layout: ' + self.L.filename + '\n'
            s = s + 'Number of Agents: ' + str(self.Nag) + '\n'
            s = s + 'Number of Access Points: ' + str(self.Nap) + '\n'
            s = s + 'Number of Acces points:' + str(self.Nap) + '\n\n'

            s = s + 'Considered links:\n'
            s = s + '-----------------\n'

            if len(self.links) == 0:
                s = s + 'No link generated\n'
            else:
                s = s + "Body 2 body: " + str(self.B2B) + '\n'
                s = s + "On Body: " + str(self.OB) + '\n'
                s = s + "Body 2 Infrastructure " + str(self.B2I) + '\n'

        except:
            s = 'No trajectory loaded'

        return s

    def load_config(self, _filesimul):
        """  load a simultraj configuration file

        Parameters
        ----------

        _filesimul : string
            name of simulation file to be loaded
        """
        self.filesimul = _filesimul
        filesimul = pyu.getlong(self.filesimul, "ini")

        config = ConfigParser.ConfigParser()
        config.read(filesimul)
        sections = config.sections()
        di = {}

        for section in sections:
            di[section] = {}
            options = config.options(section)
            for option in options:
                di[section][option] = config.get(section, option)

        # get the layout
        self.L = Layout(di['layout']['layout'])
        # get the trajectory
        traj = tr.Trajectories()
        traj.loadh5(di['trajectory']['traj'])
        self._trajname = di['trajectory']['traj']
        # resample trajectory
        for ut, t in enumerate(traj):
            if t.typ == 'ag':
                t = t.resample(2)
                person = Body(t.name + '.ini')
                self.dpersons.update({t.name: person})
                self.time = t.time()

            else:
                pos = np.array([t.x[0], t.y[0], t.z[0]])
                self.dap.update({t.ID: {'pos': pos,
                                        'ant': antenna.Antenna(),
                                        'name': t.name
                                        }
                                 })
            traj[ut] = t
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
            node a id in self.N (Netwrok)
        nb : string:
            node b id in self.N (Netwrok)
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
        """ Run teh evaluation of link along a trajectory

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
        't': list
            list of timestamp to be evaluated 
            (if [], all timestamps are considered)


        """
        defaults = {'OB': True,
                    'B2B': False,
                    'B2I': True,
                    'I2I': False,
                    'llink': [],
                    'wstd': [],
                    't': [],
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
        todo = []
        if OB:
            todo.append('OB')
        if B2B:
            todo.append('B2B')
        if B2I:
            todo.append('B2I')
        if I2I:
            todo.append('I2I')

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

        # Check time attribute
        if kwargs['t'] == []:
            t = self.time
        else:
            if kwargs['t'][0] >= self.time[0] and\
               kwargs['t'][-1] <= self.time[-1]:
                t = self.time[kwargs['t']]
            else:
                raise AttributeError('Requested timestamp not available')

        #
        # Code
        #
        init=True
        for ut, it in enumerate(t):
            # if a bodies are involved in simulation
            if (('OB' in todo) or
                    ('B2B' in todo) or
                    ('B2I' in todo)):
                nodeid = []
                pos = []
                orient = []
                for up, person in enumerate(self.dpersons.values()):
                    person.settopos(self.traj[up], t=t[ut], cs=True)
                    name = person.name
                    dev = person.dev.keys()
                    nodeid.extend([n + '_' + name for n in dev])
                    pos.extend([person.dcs[d][:, 0] for d in dev])
                    orient.extend([person.acs[d] for d in dev])
                # in a future version , the network update must also update
                # antenna positon in the device coordinate system
                self.N.update_pos(nodeid, pos, now=it)
                self.N.update_orient(nodeid, orient, now=it)
            # TODO : to be moved on the network edges
            self.N.update_dis()

            for w in wstd:
                for na, nb, typ in llink[w]:
                    if typ in todo:
                        eng = 0
                        self.evaldeter(na, nb, w)
                        if typ == 'OB':
                            self.evalstat(na, nb)
                            eng = self.SL.eng
                            L = self.DL + self.SL
                            self._ak = L.H.ak
                            self._tk = L.H.tk
                        else : 
                            self._ak = self.DL.H.ak
                            self._tk = self.DL.H.tk

                        # kw = {'d': self._ddis[(na, nb)],
                        #       'eng': eng,
                        #       'typ': typ,
                        #       'fcghz': self.N.node[na]['wstd'][w]['fcghz'],
                        #       'fbminghz': self.DL.fmin,
                        #       'fbmaxghz': self.DL.fmax,
                        #       'fstep': self.DL.fstep,
                        #       'Link_sig_grpname': self.DL.dexist['sig']['grpname'],
                        #       'Link_ray_grpname': self.DL.dexist['ray']['grpname'],
                        #       'Link_Ct_grpname': self.DL.dexist['Ct']['grpname'],
                        #       'Link_H_grpname': self.DL.dexist['H']['grpname'],
                        #       }
                        # self._saveh5(ut, na, nb, w, **kw)

                    self.data = self.data.append(pd.DataFrame({\
                                't': ut,
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
                                'sig_id': self.DL.dexist['sig']['grpname'],
                                'ray_id': self.DL.dexist['ray']['grpname'],
                                'Ct_id': self.DL.dexist['Ct']['grpname'],
                                'H_id': self.DL.dexist['H']['grpname'],
                                            },columns=['t','id_a', 'id_b',
                                          'x_a', 'y_a', 'z_a',
                                          'x_b', 'y_b', 'z_b',
                                          'd', 'eng', 'typ',
                                          'wstd', 'fcghz',
                                          'fbminghz', 'fbmaxghz', 'fstep',
                                          'sig_id', 'ray_id', 'Ct_id', 'H_id'
                                          ],index=[it]))
   
                    self.tocsv(ut, na, nb, w,init=init)
                    init=False
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

    # def _saveh5(self, ut, ida, idb, wstd, **kwargs):
    #     """ Save in h5py format

    #     Parameters
    #     ----------

    #     ut : int
    #         time index in self.time
    #     ida : string
    #         node a index
    #     idb : string
    #         node b index
    #     wstd : string
    #         wireless standard of used link

    #     kwargs (attribute of dataset):

    #     d : ndarray
    #         distrance between nodes
    #     eng : float
    #         engagement
    #     typ : string
    #         type of link (OB|B2B|B2I|I2I)
    #     fcghz : float
    #         central frequency of evaluation
    #     fbminghz : float
    #         min frequency of bandwidth
    #     fbmaxghz : float
    #         max frequency of bandwidth
    #     fstep : int
    #         frequency step
    #     Link_sig_grpname : string
    #         DLink hdf5 identifier for sig
    #     Link_ray_grpname : string
    #         DLink hdf5 identifier for ray
    #     Link_Ct_grpname : string
    #         DLink hdf5 identifier for Ct
    #     Link_H_grpname : string
    #         DLink hdf5 identifier for H

    #     Notes
    #     -----

    #     Dataset organisation:

    #     simultraj_<trajectory_filename.h5>.h5
    #         |
    #         |time
    #         |    ...
    #         |
    #         |/<tidx_ida_idb_wstd>/ |attrs
    #         |                      |a_k
    #         |                      |t_k


    #     Root dataset :
    #     time : array
    #         range of simulation time

    #     Group identifier :
    #         tidx : index in time dataset
    #         ida : node a index in Network
    #         idb : node b index in Network
    #         wstd : wireless standar of link interest


    #     Inside group:
    #         a_k : alpha_k values
    #         t_k : tau_k values

    #     Attributes of group:
    #     'd': distance between nodes
    #     'eng': engagement
    #     'typ': type of link (OB,B2B,B2I,I2I),
    #     'fcghz': central frequency of evaluation
    #     'fbminghz': min freq of evalutaion
    #     'fbmaxghz': max freq of evalutaion
    #     'fstep' : numberof evaluation point for frequency
    #     'Link_sig_grpname': self.Li.dexist['sig']['grpname'],
    #     'Link_ray_grpname': self.Li.dexist['ray']['grpname'],
    #     'Link_Ct_grpname': self.Li.dexist['Ct']['grpname'],
    #     'Link_H_grpname': self.Li.dexist['H']['grpname'],

    #     See Also
    #     --------

    #     pylayers.simul.links

    #     """

    #     filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
    #     grpname = str(ut) + '_' + ida + '_' + idb + '_' + wstd
    #     # try/except to avoid loosing the h5 file if
    #     # read/write error
    #     try:
    #         fh5 = h5py.File(filenameh5, 'a')
    #         if not grpname in fh5.keys():
    #             fh5.create_group(grpname)
    #         else:
    #             print grpname + ' already exists in ' + filenameh5
    #         f = fh5[grpname]
    #         for k in kwargs:
    #             f.attrs[k] = kwargs[k]

    #         f.create_dataset('alphak',
    #                          shape=self._ak.shape,
    #                          maxshape=(None),
    #                          data=self._ak)
    #         f.create_dataset('tauk',
    #                          shape=self._tk.shape,
    #                          maxshape=(None),
    #                          data=self._tk)
    #         fh5.close()
    #     except:
    #         fh5.close()
    #         raise NameError('Simultraj._saveh5: issue when writting h5py file')


    # def _loadh5(self, t, ida, idb, wstd):
    #     """ Load in h5py format

    #     Parameters
    #     ----------

    #    ut : int
    #        time index in self.time
    #    ida : string
    #        node a index
    #    idb : string
    #        node b index
    #    wstd : string
    #        wireless standard of used link

    #     Returns
    #     -------
    #     (ak, tk, conf)

    #     ak : ndarray:
    #         alpha_k
    #     tk : ndarray:
    #         alpha_k
    #     conf : dict
    #         dictionnary containing configuration setup
    #     """

    #     filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
    #     grpname = str(ut) + '_' + ida + '_' + idb + '_' + wstd
    #     # try/except to avoid loosing the h5 file if
    #     # read/write error
    #     try:
    #         fh5 = h5py.File(filenameh5, 'r')
    #         if not grpname in fh5.keys():
    #             fh5.close()
    #             raise NameError(grpname + ' cannot be reached in ' + self.filename)
    #         f = fh5[grpname]
    #         conf={}
    #         for k in f.attrs.keys():
    #             conf[k]=f.attrs[k]
    #         ak = f['alphak']
    #         tk = f['tauk']
    #         fh5.close()

    #         return ak, tk, conf
    #     except:
    #         fh5.close()
    #         raise NameError('Simultraj._loadh5: issue when reading h5py file')


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
