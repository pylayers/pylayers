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
    Simul.run2
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
import pdb
import numpy as np
import matplotlib.pylab as plt
import pylayers.util.pyutil as pyu
import pylayers.signal.waveform as wvf
import pylayers.signal.bsignal as bs
from pylayers.signal.device import Device

# Handle Layout
from pylayers.gis.layout import Layout
# Handle VectChannel and ScalChannel
from pylayers.antprop import antenna, signature
from pylayers.network.network import Network
from pylayers.simul.link import *
# Handle directory hierarchy
from pylayers.util.project import *
# Handle UWB measurements
import pdb
import pylayers.mobility.trajectory as tr
from pylayers.mobility.ban.body import *
from pylayers.antprop.statModel import *


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
        self.Li = Link(L=self.L,verbose=verbose)
        self.filename = 'simultraj_' + self._trajname
        self._saveh5_init()

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

    def _saveh5_init(self):
        """ initialization of the h5py file
        """
        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        import ipdb
        try:
            f5 = h5py.File(filenameh5, 'w')
            f5.create_dataset('time', shape=self.time.shape, data=self.time)
            f5.close()
        except:
            f5.close()
            raise NameError('simultra.saveinit: \
                            issue when writting h5py file')

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
        self.Li.a = self.N.node[na]['p']
        self.Li.Ta = self.N.node[na]['T']
        self.Li.b = self.N.node[nb]['p']
        self.Li.Tb = self.N.node[nb]['T']
        if fmode == 'center':
            self.Li.fGHz = self.N.node[na]['wstd'][wstd]['fcghz']
        else:
            minb = self.N.node[na]['wstd'][wstd]['fbminghz']
            maxb = self.N.node[na]['wstd'][wstd]['fbmaxghz']
            self.Li.fGHz = np.linspace(minb, maxb, nf)

        a, t = self.Li.eval()

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

        (a, t )

        a : ndarray
            alpha_k
        t : ndarray
            tau_k
        """

        pa = self.N.node[na]['p']
        pb = self.N.node[nb]['p']
        dida, name = na.split('_')
        didb, name = nb.split('_')

        # inter to be replace by engaement
        inter = self.dpersons[name].intersectBody3(pa, pb, topos=True)

        condition = 'nlos'
        if inter == 1:
            condition = 'los'

        empA = self.dpersons[name].dev[dida]['cyl']
        empB = self.dpersons[name].dev[didb]['cyl']

        emp = empA
        if empA == 'trunkb':
            emp = empB
        if emp == 'forearml':
            emp = 'forearmr'
        if emp == 'trunku':
            condition = 'los'
        a, t = getchannel(
            emplacement=emp, condition=condition, intersection=inter)

        return a, t, inter

    def run2(self, **kwargs):
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
            self._ddis = self.N.update_dis()

            for w in wstd:
                for na, nb, typ in llink[w]:
                    if typ in todo:
                        eng = 0
                        if typ == 'OB':
                            _ak, _tk, eng = self.evalstat(na, nb)
                        else:
                            _ak, _tk = self.evaldeter(na, nb, w)
                            # fill panda dataframe 2D trajectory
                        self._ak = _ak
                        self._tk = _tk
                        kw = {'d': self._ddis[(na, nb)],
                              'eng': eng,
                              'typ': typ,
                              'fcghz': self.N.node[na]['wstd'][w]['fcghz'],
                              'fbminghz': self.Li.fmin,
                              'fbmaxghz': self.Li.fmax,
                              'fstep': self.Li.fstep,
                              'Link_sig_grpname': self.Li.dexist['sig']['grpname'],
                              'Link_ray_grpname': self.Li.dexist['ray']['grpname'],
                              'Link_Ct_grpname': self.Li.dexist['Ct']['grpname'],
                              'Link_H_grpname': self.Li.dexist['H']['grpname'],
                              }

                        self._saveh5(ut, na, nb, w, **kw)

    def _saveh5(self, ut, ida, idb, wstd, **kwargs):
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

        kwargs (attribute of dataset):

        d : ndarray
            distrance between nodes
        eng : float
            engagement
        typ : string
            type of link (OB|B2B|B2I|I2I)
        fcghz : float
            central frequency of evaluation
        fbminghz : float
            min frequency of bandwidth
        fbmaxghz : float
            max frequency of bandwidth
        fstep : int
            frequency step
        Link_sig_grpname : string
            Link hdf5 identifier for sig
        Link_ray_grpname : string
            Link hdf5 identifier for ray
        Link_Ct_grpname : string
            Link hdf5 identifier for Ct
        Link_H_grpname : string
            Link hdf5 identifier for H

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

        Attributes of group:
        'd': distance between nodes
        'eng': engagement
        'typ': type of link (OB,B2B,B2I,I2I),
        'fcghz': central frequency of evaluation
        'fbminghz': min freq of evalutaion
        'fbmaxghz': max freq of evalutaion
        'fstep' : numberof evaluation point for frequency
        'Link_sig_grpname': self.Li.dexist['sig']['grpname'],
        'Link_ray_grpname': self.Li.dexist['ray']['grpname'],
        'Link_Ct_grpname': self.Li.dexist['Ct']['grpname'],
        'Link_H_grpname': self.Li.dexist['H']['grpname'],

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
            else:
                print grpname + ' already exists in ' + filenameh5
            f = fh5[grpname]
            for k in kwargs:
                f.attrs[k] = kwargs[k]

            f.create_dataset('alphak',
                             shape=self._ak.shape,
                             maxshape=(None),
                             data=self._ak)
            f.create_dataset('tauk',
                             shape=self._tk.shape,
                             maxshape=(None),
                             data=self._tk)
            fh5.close()
        except:
            fh5.close()
            raise NameError('Simultraj._saveh5: issue when writting h5py file')


    def _loadh5(self, t, ida, idb, wstd):
        """ Load in h5py format

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

        Returns
        -------
        (ak, tk, conf)

        ak : ndarray:
            alpha_k
        tk : ndarray:
            alpha_k
        conf : dict
            dictionnary containing configuration setup
        """

        filenameh5 = pyu.getlong(self.filename, pstruc['DIRLNK'])
        grpname = str(ut) + '_' + ida + '_' + idb + '_' + wstd
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            fh5 = h5py.File(filenameh5, 'r')
            if not grpname in fh5.keys():
                fh5.close()
                raise NameError(grpname + ' cannot be reached in ' + self.filename)
            f = fh5[grpname]
            conf={}
            for k in f.attrs.keys():
                conf[k]=f.attrs[k]
            ak = f['alphak']
            tk = f['tauk']
            fh5.close()

            return ak, tk, conf
        except:
            fh5.close()
            raise NameError('Simultraj._loadh5: issue when reading h5py file')



    def gen_links(self, B2B=False, OB=True, B2I=False):
        """  DEPRECATED
        generate links
        Parameters
        ----------

        B2B : boolean
            Body to Body links
        OB  : boolean
            On-Body links
        B2I : boolean
            Body to infrastructure links



        """

        self.B2B = B2B
        self.OB = OB
        self.B2I = B2I
        lAP = self.dap.keys()
        lperson = self.dpersons.values()
        dlink = {}
        dlink['OB'] = []
        dlink['B2B'] = []
        dlink['B2I'] = []
        # pdb.set_trace()
        for person in lperson:
            if OB:
                for dev1 in person.dev:
                    if person.dev[dev1]['typ'] == 'dynamic':
                        for dev2 in person.dev:
                            if dev2 <> dev1:
                                dlink['OB'].append(
                                    [(person.name, dev1), (person.name, dev2)])
            if B2I:
                for ap in lAP:
                    # may be criteria on distance between person and AP
                    for dev in person.dev:
                        # if person.dev[dev1]['infra']:
                        # in the future port would be an antenna port for MIMO
                        # AP
                        dlink['B2I'].append([(ap, 'port'), (person.name, dev)])

            if B2B:
                for dev1 in person.dev:
                    # if person.dev[dev1]['tobody']:
                    for alter in lperson:
                        if alter.name <> person.name:
                            for dev2 in alter.dev:
                                # if alter.dev[dev2]['tobody']:
                                dlink['B2B'].append(
                                    [(person.name, dev1), (alter.name, dev2)])

        self.links = dlink

    def runB2I(self, llink=[], t=[]):
        """  DEPRECATED
        """

        if llink == []:
            llink = self.links['B2I']
        if t == []:
            time = self.traj[0].time()
        else:
            time = self.traj[0].time()[t]

        n_time = len(time)
        n_links = len(llink)
        Kmax = 100

        L = Link(L=self.L)

        for kt in range(0, n_time):

            #~ if time[kt]%int(time[kt]) == 0:
                #~ print 't = ', time[kt]
            print 't = ', time[kt]
            for kp, person in enumerate(self.dpersons.values()):
                person.settopos(self.traj[kp], t=time[kt], cs=True)

            for kl, l in enumerate(llink):

                # if first extremiti is an access point
                la = l[0]
                lb = l[1]
                # if a is an acces point
                import ipdb
                ipdb.set_trace()
                if la[0] in self.dap:
                    L.a = self.dap[la[0]]['pos']
                    L.Aa = self.dap[la[0]]['ant']
                    L.Ta = np.eye(3)

                else:
                    L.a = self.dpersons[la[0]].dcs[la[1]][:, 0]
                    L.Aa = antenna.Antenna(
                        self.dpersons[la[0]].dev[la[1]]['file'])
                    L.Ta = self.dpersons[la[0]].acs[la[1]]

                if lb[0] in self.dap:
                    L.b = self.dap[lb[0]]['pos']
                    L.Ab = self.dap[lb[0]]['ant']
                    L.Tb = np.eye(3)
                else:
                    L.b = self.dpersons[lb[0]].dcs[lb[1]][:, 0]
                    L.Ab = antenna.Antenna(
                        self.dpersons[lb[0]].dev[lb[1]]['file'])
                    L.Tb = self.dpersons[lb[0]].acs[lb[1]]

                L.eval()

                H = L.H
                H.applyFriis()
                ntraj = H.y.shape[0]
                alphak = np.zeros(shape=(Kmax))
                tauk = np.zeros(shape=(Kmax))
                if ntraj < Kmax:
                    # pdb.set_trace()
                    alphak[0:ntraj] = np.real(
                        np.sqrt(np.sum(H.y * np.conj(H.y), axis=1)) / len(self.fGHz))
                    tauk[0:ntraj] = H.tau0
                else:
                    alphak = np.real(
                        np.sqrt(np.sum(H.y * np.conj(H.y), axis=1)) / len(self.fGHz))[0:Kmax]
                    tauk = (H.tau0)[0:Kmax]
                    print ' warning ntraj > Kmax'
                tab = np.vstack((alphak, tauk)).T
                resultEnv[kt, kl,:,:] = tab
                print 'link = ', link, '  pr  = ', 10 * np.log10(sum((abs(alphak) ** 2)))

        return resultEnv  # , resultOb

    def run(self, llink=[], t=[]):
        """
        DEPRECATED

        Parameters
        ----------

        llink : list

        """

        if llink == []:
            llink = self.links
        if t == []:
            time = self.traj[0].time()
        else:
            time = self.traj[0].time()[t]

        n_time = len(time)
        n_links = len(llink)
        Kmax = 100
        resultEnv = np.zeros(shape=(n_time, n_links, Kmax, 2))
        resultOb = np.zeros(shape=(n_time, n_links, Kmax, 2))

        for kt in range(0, n_time):

            #~ if time[kt]%int(time[kt]) == 0:
                #~ print 't = ', time[kt]
            print 't = ', time[kt]
            for kp, person in enumerate(self.dpersons.values()):
                person.settopos(self.traj[kp], t=time[kt], cs=True)

            for kl in range(0, n_links):
                link = llink[kl]
                A = link[0]
                B = link[1]
                pA = self.dpersons[A[0]].dcs[A[1]][:, 0]
                TA = self.dpersons[A[0]].acs[A[1]]
                pB = self.dpersons[B[0]].dcs[B[1]][:, 0]
                TB = self.dpersons[B[0]].acs[B[1]]
                cylA = self.dpersons[A[0]].dev[A[1]]['cyl']
                cylB = self.dpersons[B[0]].dev[B[1]]['cyl']

                #interA = self.dpersons[A[0]].intersectBody3(pA,pB, topos = True)
                ##interB = self.dpersons[B[0]].intersectBody2(pA,pB, topos = True)

                #condition ='nlos'
                # if interA==1:
                    #condition = 'los'
                #empA = A[1]
                #empB = B[1]
                #emp  = empA
                # if empA == 'accelerometer':
                    #emp = empB
                # if emp == 'left_watch':
                    #emp = 'right_watch'
                    # ~ if condition == 'los':
                        # ~ condition = 'nlos'
                    # ~ else:
                        # ~ condition = 'los'
                # if emp == 'front_chest':
                    #condition = 'los'
                #alphakOb, taukOb = getchannel(emplacement = emp ,condition = condition, intersection = interA)
                #alphak =  np.zeros(shape=(Kmax))
                #tauk   =  np.zeros(shape=(Kmax))
                #ntraj = len(taukOb)
                # if ntraj < Kmax:
                    #alphak[0:ntraj] =  alphakOb[0:ntraj]
                    #tauk[0:ntraj] =  taukOb[0:ntraj]
                # else:
                    #alphak =  alphakOb[0:Kmax]
                    #tauk =  taukOb[0:Kmax]
                    # print ' warning ntraj > Kmax'
                #tab = np.vstack((alphak,tauk)).T
                #resultOb[kt,kl,:,:]= tab
                cycA = self.L.pt2cy(pt=pA)
                cycB = self.L.pt2cy(pt=pB)
                sig = signature.Signatures(self.L, cycA, cycB)
                #sig.run4(cutoff =1,algo='old')
                sig.run5(cutoff=3)
                tx_2D = pA[0:2]
                rx_2D = pB[0:2]
                r2d = sig.rays(tx_2D, rx_2D)
                r2d.pTx = pA
                r2d.pRx = pB
                r3d = r2d.to3D(self.L)
                r3d.locbas(self.L)
                r3d.fillinter(self.L)
                Cn = r3d.eval(fGHz=self.fGHz)
                Cn.locbas(Tt=TA, Tr=TB)
                AntA = antenna.Antenna(self.dpersons[A[0]].dev[A[1]]['file'])
                AntB = antenna.Antenna(self.dpersons[B[0]].dev[B[1]]['file'])
                H = Cn.prop2tran(a=AntA, b=AntB)
                H.applyFriis()
                ntraj = H.y.shape[0]
                alphak = np.zeros(shape=(Kmax))
                tauk = np.zeros(shape=(Kmax))
                if ntraj < Kmax:
                    # pdb.set_trace()
                    alphak[0:ntraj] = np.real(
                        np.sqrt(np.sum(H.y * np.conj(H.y), axis=1)) / len(self.fGHz))
                    tauk[0:ntraj] = H.tau0
                else:
                    alphak = np.real(
                        np.sqrt(np.sum(H.y * np.conj(H.y), axis=1)) / len(self.fGHz))[0:Kmax]
                    tauk = (H.tau0)[0:Kmax]
                    print ' warning ntraj > Kmax'
                tab = np.vstack((alphak, tauk)).T
                resultEnv[kt, kl,:,:] = tab
                print 'link = ', link, '  pr  = ', 10 * np.log10(sum((abs(alphak) ** 2)))

        return resultEnv  # , resultOb

    def runEnv(self, llink=[], t=[], show=False):
        """ DEPRECATED
        Parameters
        ----------

        llink : list

        """

        if llink == []:
            llink = self.links
        if t == []:
            time = self.traj[0].time()
        else:
            time = self.traj[0].time()[t]

        n_time = len(time)
        n_links = len(llink)
        Kmax = 100

        L = Link(L=self.L)

        resultEnv = np.zeros(shape=(n_time, n_links, Kmax, 2))
        taumin = 0
        taumax = 300
        taustep = 0.1
        x = np.arange(taumin, taumax, taustep)
        y = np.zeros(len(x))
        cira = bs.TUsignal(x, y)

        for kt in range(0, n_time):

            if time[kt] % int(time[kt]) == 0:
                print 't = ', time[kt]

            for kp, person in enumerate(self.dpersons.values()):
                person.settopos(self.traj[kp], t=time[kt], cs=True)

            for kl in range(0, n_links):
                link = llink[kl]
                import ipdb
                ipdb.set_trace()
                A = link[0]
                B = link[1]
                pA = self.dpersons[A[0]].dcs[A[1]][:, 0]
                TA = self.dpersons[A[0]].acs[A[1]]
                pB = self.dpersons[B[0]].dcs[B[1]][:, 0]
                TB = self.dpersons[B[0]].acs[B[1]]
                cylA = self.dpersons[A[0]].dev[A[1]]['cyl']
                cylB = self.dpersons[B[0]].dev[B[1]]['cyl']

                cycA = self.L.pt2cy(pt=pA)
                cycB = self.L.pt2cy(pt=pB)
                sig = signature.Signatures(self.L, cycA, cycB)
                sig.run4(cutoff=1, algo='old')
                tx_2D = pA[0:2]
                rx_2D = pB[0:2]
                r2d = sig.rays(tx_2D, rx_2D)
                r2d.pTx = pA
                r2d.pRx = pB
                r3d = r2d.to3D(self.L)
                r3d.locbas(self.L)
                r3d.fillinter(self.L)
                Cn = r3d.eval(fGHz=self.fGHz)
                Cn.locbas(Tt=TA, Tr=TB)
                AntA = antenna.Antenna(self.dpersons[A[0]].dev[A[1]]['file'])
                AntB = antenna.Antenna(self.dpersons[B[0]].dev[B[1]]['file'])
                H = Cn.prop2tran(a=AntA, b=AntB)
                H.applyFriis()
                ntraj = H.y.shape[0]
                alphak = np.zeros(shape=(Kmax))
                tauk = np.zeros(shape=(Kmax))

                if ntraj < Kmax:
                    # pdb.set_trace()
                    alphak[0:ntraj] = np.real(
                        np.sqrt(np.sum(H.y * np.conj(H.y), axis=1)) / len(self.fGHz))
                    tauk[0:ntraj] = H.tau0
                else:
                    alphak = np.real(
                        np.sqrt(np.sum(H.y * np.conj(H.y), axis=1)) / len(self.fGHz))[0:Kmax]
                    tauk = (H.tau0)[0:Kmax]
                    print ' warning ntraj > Kmax'
                tab = np.vstack((alphak, tauk)).T
                resultEnv[kt, kl,:,:] = tab
                if show:
                    print 'link  = ', link, '  pr env  = ', 10 * np.log10(sum((abs(alphak) ** 2)))
                cira.aggcir(alphak, tauk)
                self.chan = resultEnv
                self.cira = cira
        # return resultEnv, cira

    def runOb(self, llink=[], t=[], show=False):
        """DEPRECATED
        Parameters
        ----------

        llink : list

        """

        if llink == []:
            llink = self.links
        if t == []:
            time = self.traj[0].time()
        else:
            time = self.traj[0].time()[t]

        n_time = len(time)
        n_links = len(llink)
        Kmax = 100
        resultEnv = np.zeros(shape=(n_time, n_links, Kmax, 2))
        resultOb = np.zeros(shape=(n_time, n_links, Kmax, 2))

        #~ prev_lstate = list(np.zeros(shape = (n_links)))
        #~ prev_alphakOb = []
        #~ prev_taukOb = []

        for kt in range(0, n_time):

            for kp, person in enumerate(self.dpersons.values()):
                person.settopos(self.traj[kp], t=time[kt], cs=True)

            for kl in range(0, n_links):
                link = llink[kl]
                A = link[0]
                B = link[1]
                pA = self.dpersons[A[0]].dcs[A[1]][:, 0]
                TA = self.dpersons[A[0]].acs[A[1]]
                pB = self.dpersons[B[0]].dcs[B[1]][:, 0]
                TB = self.dpersons[B[0]].acs[B[1]]
                cylA = self.dpersons[A[0]].dev[A[1]]['cyl']
                cylB = self.dpersons[B[0]].dev[B[1]]['cyl']

                interA = self.dpersons[
                    A[0]].intersectBody3(pA, pB, topos=True)
                import ipdb
                ipdb.set_trace()
                #interB = self.dpersons[B[0]].intersectBody2(pA,pB, topos = True)

                condition = 'nlos'
                if interA == 1:
                    condition = 'los'
                devIdA = A[1]
                devIdB = B[1]
                empA = self.dpersons[A[0]].dev[devIdA]['cyl']
                empB = self.dpersons[B[0]].dev[devIdB]['cyl']

                emp = empA
                if empA == 'trunkb':
                    emp = empB
                if emp == 'forearml':
                    emp = 'forearmr'

                if emp == 'trunku':
                    condition = 'los'
                alphakOb, taukOb = getchannel(
                    emplacement=emp, condition=condition, intersection=interA)

                alphak = np.zeros(shape=(Kmax))
                tauk = np.zeros(shape=(Kmax))
                ntraj = len(taukOb)
                if ntraj < Kmax:
                    alphak[0:ntraj] = alphakOb[0:ntraj]
                    tauk[0:ntraj] = taukOb[0:ntraj]
                else:
                    alphak = alphakOb[0:Kmax]
                    tauk = taukOb[0:Kmax]
                    print ' warning ntraj > Kmax'

                tab = np.vstack((alphak, tauk)).T
                resultOb[kt, kl,:,:] = tab
                if show:
                    print 'link  = ', link, '  pr ob  = ', 10 * np.log10(sum((abs(alphak) ** 2)))

        return resultOb
