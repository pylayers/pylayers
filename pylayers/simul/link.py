#
# -*- coding: utf-8 -*-
#
from __future__ import print_function
r"""

.. currentmodule:: pylayers.simul.link

.. autosummary::
    :members:

"""

try:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi import mlab
except:
    print('Layout:Mayavi is not installed')

import doctest
import time
import numpy as np
import matplotlib.pylab as plt
import pylayers.signal.waveform as wvf
import pylayers.util.geomutil as geu
import logging

logger = logging.getLogger(__name__)

from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.simul.radionode import RadioNode
# Handle Layout
from pylayers.gis.layout import Layout
# Handle Antenna
from pylayers.antprop.antenna import Antenna

# Handle Signature
from pylayers.antprop.signature import Signatures,Signature
# Handle Rays
from pylayers.antprop.rays import Rays
# Handle VectChannel and ScalChannel
from pylayers.antprop.channel import Ctilde, Tchannel , AFPchannel
from pylayers.antprop.statModel import getchannel
import tqdm
import copy
import h5py
import pdb

logger = logging.getLogger(__name__)

class Link(PyLayers):
    """ Link class

    Members
    -------

    """

    def __init__(self):
        """ Link evaluation metaclass
        """
        PyLayers.__init__(self)


    def __add__(self,l):
        """ merge ak and tauk of 2 Links

        """
        L  = Link()
        tk = np.hstack((self.H.tk,l.H.tk))
        ak = np.hstack((self.H.ak,l.H.ak))
        us = np.argsort(tk)
        L.H.ak = ak[us]
        L.H.tk = tk[us]
        return L


class SLink(Link):

    def __init__(self):
        """ Statistical Link evaluation class
        """
        super(SLink, self).__init__()


    def onbody(self, B, dida, didb, a, b):
        """ Statistical evaluation of a on-body link

        Parameters
        ----------

        B : Body
            Body object on which devices are held
        dida: int
            device a id number on body
        didb: int
            device b id number on body
        a : nd array
            position of device a
        b : nd array
            position of device b

        Returns
        -------

        (ak, tk, eng )

        ak : ndarray
            alpha_k
        tk : ndarray
            tau_k


        See Also
        --------

        pylayers.mobility.ban.body

        """

        # inter to be replace by engaement
        eng = B.intersectBody3(a, b, topos=True)
        empa = B.dev[dida]['cyl']
        empb = B.dev[didb]['cyl']

        emp = empa
        if empa == 'trunkb':
            emp = empb
        if emp == 'forearml':
            emp = 'forearmr'

        self.H.ak, self.H.tk = getchannel(emplacement=emp, intersection=eng)
        self.eng = eng

        return self.H.ak, self.H.tk, self.eng



class DLink(Link):
    """ Deterministic Link Class

    Attributes
    ----------

        L : Layout
            Layout to be used
        Aa : Antenna
            Antenna of device dev_a
        Ab : Antenna
            Antenna of device dev_b
        a : np.ndarray (3,)
            position of a device dev_a
        b : np.ndarray (3,)
            position of a device dev_b
        ca : int
            cycle a number
        cb : int
            cycle b number
        Ta : np.ndarray (3,3)
            Rotation matrice of Antenna of device dev_a relative to global Layout scene
        Tb : np.ndarray (3,3)
            Rotation matrice of Antenna of device dev_b relative to global Layout scene
        fGHz : np.ndarray (Nf,)
            frequency range of Nf points used for evaluation of channel
        wav : Waveform
            Waveform to be applied on the channel
        save_idx : int
            number to identify the h5 file generated

    """

    def __init__(self, **kwargs):
        """ deterministic link evaluation


        Advanced (change only if you really know what you do !)

        save_opt : list (['sig','ray','Ct','H'])
            information to be saved in the Links h5 file. Should never be Modified !

        force_create : Boolean (False)
            force creating the h5py file (if already exist, will be erased)


        Notes
        -----

        All simulations are stored into a unique file in your <PyProject>/output directory
        using the following convention:

        Links_<save_idx>_<LayoutFilename>.h5

        where
            <save_idx> is an integer number to distinguish different links simulations
        and <LayoutFilename> is the Layout used for the link simulation.



        Dataset organisation:

        Links_<idx>_<Layout_name>.h5
            |
            |/sig/si_ID#0/
            |    /si_ID#1/
            |    ...
            |
            |/ray/ray_ID#0/
            |    /ray_ID#1/
            |    ...
            |
            |/Ct/Ct_ID#0/
            |   /Ct_ID#1/
            |    ...
            |
            |/H/H_ID#0/
            |  /H_ID#1/
            |    ...
            |
            |
            |p_map
            |c_map
            |f_map
            |A_map
            |T_map



        Roots Dataset :

        c_map : Cycles (Nc x 3)
        p_map : Positions (Np x 3)
        f_map : Frequency (Nf x 3)
        T_map : Rotation matrices (Nt x 3)
        A_map : Antenna name (Na x 3)

        Groups and subgroups:


            Signature identifier (si_ID#N):
                ca_cb_cutoff_th

            Ray identifier (ray_ID#N):
                cutoff_th_ua_ub

            Ctilde identifier (Ct_ID#N):
                ua_ub_uf

            H identifier (H_ID#N):
                ua_ub_uf_uTa_uTb_uAa_uAb

            with
            ca : cycle number of a
            cb : cycle number of b
            cutoff : signature.run cutoff
            th : signature.run threshold * 100
            ua : indice of a position in 'p_map' position dataset
            ub : indice of a position in 'p_map' position dataset
            uf : indice of freq position in 'f_map' frequency dataset
            uTa : indice of a position in 'T_map' Rotation dataset
            uTb : indice of a position in 'T_map' Rotation dataset
            uAa : indice of a position in 'A_map' Antenna name dataset
            uAb : indice of b position in 'A_map' Antenna name dataset



        Examples
        --------

            >>> from pylayers.simul.link import *
            >>> DL = DLink(L=Layout('DLR.lay'))
            >>> DL.eval()

        Notes
        -----

        When modifying the coordinates of the link it is important to modify
        the array and not one of its component. In that case the cycle will not
        be updated.


        """


        Link.__init__(self)

        defaults={ 'L': '',
                   'a': np.array(()),
                   'b': np.array(()),
                   'Aa': [],
                   'Ab': [],
                   'Ta': np.eye(3),
                   'Tb': np.eye(3),
                   'fGHz': np.array([2.4]),
                   'wav': wvf.Waveform(),
                   'cutoff': 3,
                   'threshold': 0.8,
                   'delay_excess_max_ns':500,
                   'save_opt': ['sig','ray2','ray','Ct','H'],
                   'save_idx':0,
                   'force_create':False,
                   'seed':0,
                   'graph':'tcvirw'
                }



        # self._ca = -1
        # self._cb = -1

        specset  = ['a','b','Aa','Ab','Ta','Tb','cutoff','L','fGHz','wav']

        # set default attribute

        for key, value in defaults.items():

            if key not in kwargs:
                if key in specset :
                    setattr(self,'_'+key,value)
                else :
                    setattr(self,key,value)
            else :
                if key in specset :
                    setattr(self,'_'+key,kwargs[key])

                else :
                    setattr(self,key,kwargs[key])

        # if self._L == '':
        #     raise AttributeError('Please specify a Layout')


        force = self.force_create
        delattr(self,'force_create')


        # dictionnary data exists
        self.dexist={'sig':{'exist':False,'grpname':''},
                     'ray':{'exist':False,'grpname':''},
                     'ray2':{'exist':False,'grpname':''},
                     'Ct':{'exist':False,'grpname':''},
                     'H':{'exist':False,'grpname':''}
                    }

        # The link frequency range depends on the antenna
        # self.fGHz = kwargs['fGHz']

        if self.Aa == []:
            self.Aa = Antenna(typ='Omni',fGHz=self.fGHz)
        if self.Ab == []:
            self.Ab = Antenna(typ='Omni',fGHz=self.fGHz)

        if isinstance(self._L,str):
            self._Lname = self._L
            self._L = Layout(self._Lname,bgraphs=True,bcheck=False)
        else:
            self._Lname = self._L._filename



        if self._Lname != '':

            self.filename = 'Links_' + str(self.save_idx) + '_' + self._Lname + '.h5'
            filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])
            # check if save file already exists
            if not os.path.exists(filenameh5) or force:
                print('Links save file for ' + self.L._filename + ' does not exist.')
                print('Creating file. You\'ll see this message only once per Layout')
                self.save_init(filenameh5)
            try:
                self.L.dumpr()
                print('Layout Graph loaded')
            except:
                print('This is the first time the Layout is used. Graphs have to be built. Please Wait')
                self.L.build(graph=self.graph)
                self.L.dumpw()

            #
            # In outdoor situation we delete transmission node involving
            # an indoor cycle at the exception of AIR
            #

            cindoor = [p for p in self.L.Gt.nodes() if self.L.Gt.node[p]['indoor']]

            if self._L.typ =='outdoor':
                u = self.L.Gi.node.keys()
                # lT : list of transmission interactions
                lT  =  [k for k in u if (len(k)==3)]
                # lTi : transmission connected at least to an indoor cycle
                lTi = [ k for k in lT if ((k[1]  in cindoor) or (k[2] in cindoor))]
                # lTiw : those which are walls (not those above buildings)
                lTiw = [ k for k in lTi if self.L.Gs.node[k[0]]['name']!='AIR' ]

                self.L.Gi.remove_nodes_from(lTiw)
                lE = list(self.L.Gi.edges())
                for k in range(len(lE)):
                    e = lE[k]
                    try:
                        output = self.L.Gi[e[0]][e[1]]['output']
                    except:
                        pdb.set_trace()
                    tbd = []
                    for l in output.keys():
                        if l in lTiw:
                            tbd.append(l)
                    for d in tbd : 
                        del output[d]

                    self.L.Gi[e[0]][e[1]]['output']=output
            #self.L.dumpw()
            #self.L.build()

            self.init_positions()



            ###########
            # init freq
            # TODO Check where it is used redundant with fGHz
            ###########
            #self.fmin  = self.fGHz[0]
            #self.fmax  = self.fGHz[-1]
            #self.fstep = self.fGHz[1]-self.fGHz[0]


            self.Si = Signatures(self.L,self.ca,self.cb,cutoff=self.cutoff)
            self.R = Rays(self.a,self.b)
            self.C = Ctilde()
            self.H = Tchannel()


    @property
    def Lname(self):
        return self._Lname

    @property
    def L(self):
        return self._L

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def ca(self):
        return self._ca

    @property
    def cb(self):
        return self._cb

    @property
    def Aa(self):
        return self._Aa

    @property
    def Ab(self):
        return self._Ab

    @property
    def Ta(self):
        return self._Ta

    @property
    def Tb(self):
        return self._Tb

    @property
    def fGHz(self):
        return self._fGHz

    @property
    def cutoff(self):
        return self._cutoff

    @property
    def delay_excess_max_ns(self):
        return self._delay_excess_max_ns

    @property
    def threshold(self):
        return self._threshold

    @property
    def wav(self):
        return self._wav

    @L.setter
    def L(self,L):
        # change layout and build/load
        plotfig=False
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            mlab.clf()
            plotfig=True
        if isinstance(L,str):
            self._L = Layout(L,bgraphs=False,bcheck=False)
            self._Lname = L
        elif isinstance(L,Layout):
            self._L = L
            self._Lname = L._filename

        self.reset_config()

        if plotfig:
            self._show3()

    @Lname.setter
    def Lname(self,Lname):
        # change layout and build/load
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            mlab.clf()
        self._L = Layout(Lname)
        self._Lname = Lname
        self.reset_config()


    @a.setter
    def a(self,position):
        if not self.L.ptin(position):
            # limit position in the visible region L.ax
            #if position[0] < self.L.ax[0]:
                #position[0] = self.L.ax[0]
            # if position[0] > self.L.ax[1]:
            #     position[0] = self.L.ax[1]
            # if position[1] < self.L.ax[2]:
            #     position[1] = self.L.ax[2]
            # if position[1] > self.L.ax[3]:
            #     position[1] = self.L.ax[3]
            raise NameError ('Warning : point a is not inside the Layout')
            # raise NameError ('Warning : point a is not inside the Layout')
        if not self.L.pt2cy(position) == self.ca:
            self.ca = self.L.pt2cy(position)
        self._a = position
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='a',delrays=True)


        if hasattr(self,'ca') and hasattr(self,'cb'):
            self._autocufoff()
            self.checkh5()

    @b.setter
    def b(self,position):
        if not self.L.ptin(position):
            # if position[0]<self.L.ax[0]:
            #     position[0]=self.L.ax[0]
            # if position[0]>self.L.ax[1]:
            #     position[0]=self.L.ax[1]
            # if position[1]<self.L.ax[2]:
            #     position[1]=self.L.ax[2]
            # if position[1]>self.L.ax[3]:
            #     position[1]=self.L.ax[3]
            raise NameError ('Warning : point b is not inside the Layout')
        if not self.L.pt2cy(position) == self.cb:
            self.cb = self.L.pt2cy(position)
        self._b = position
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='b',delrays=True)

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self._autocufoff()
            self.checkh5()

    @ca.setter
    def ca(self,cycle):
        if not cycle in self.L.Gt.nodes():
            raise NameError ('cycle ca is not inside Gt')

        self._ca = cycle
        self.a = self.L.cy2pt(cycle)

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()

    @cb.setter
    def cb(self,cycle):
        if not cycle in self.L.Gt.nodes():
            raise NameError ('cycle cb is not inside Gt')
        self._cb = cycle
        self.b = self.L.cy2pt(cycle)

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()

    @Aa.setter
    def Aa(self,Ant):

        if hasattr(self.Aa,'_mayamesh'):
            self.Aa._mayamesh.remove()
        # save rot
        rot = self.Ta
        self._Aa = Ant
        self.Ta = rot
        if Ant.fromfile:
            self.fGHz = Ant.fGHz
            print("Warning : frequency range modified by antenna Aa")
        else:
            self._Aa.fGHz = self.fGHz

        # self.initfreq()

        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='a')

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()


    @Ab.setter
    def Ab(self,Ant):
        if hasattr(self.Ab,'_mayamesh'):
           self.Ab._mayamesh.remove()

        #save rot
        rot = self.Tb
        self._Ab = Ant
        self.Tb = rot
        if Ant.fromfile:
           self.fGHz = Ant.fGHz
           print("Warning : frequency range modified by antenna Ab")
        else:
           self._Ab.fGHz = self.fGHz

        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='b')

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()


    @Ta.setter
    def Ta(self,orientation):
        self._Ta = orientation
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='a')

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()



    @Tb.setter
    def Tb(self,orientation):
        self._Tb = orientation
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='b')

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()

        # if self.dexist['Ct']['exist']:
        #     self.C.locbas(Ta=self.Ta, Tb=self.Tb)
        #     #T channel
        #     self.H = self.C.prop2tran(a=self.Aa,b=self.Ab,Friis=True)


    @cutoff.setter
    def cutoff(self,cutoff):

        co = max(cutoff,1)

        self._cutoff=co
        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()


    @threshold.setter
    def threshold(self,threshold):

        th = min(threshold,1.)
        th = max(threshold,0.)
        self._threshold= th

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()

    @delay_excess_max_ns.setter
    def delay_excess_max_ns(self,delay_excess_max_ns):

        delay_excess_max_ns = max(delay_excess_max_ns,0.)

        self._delay_excess_max_ns = delay_excess_max_ns

        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5


    @fGHz.setter
    def fGHz(self,freq):
        if not isinstance(freq,np.ndarray):
            freq=np.array([freq])

        diff_freq_a = (self.Aa.fGHz!=freq)
        diff_freq_b = (self.Ab.fGHz!=freq)

        if isinstance(diff_freq_a,bool):
            cond_a = diff_freq_a
        else: 
            cond_a = diff_freq_a.all()

        if  isinstance(diff_freq_b,bool):
            cond_b = diff_freq_b
        else: 
            cond_b = diff_freq_b.all()


        if (self.Aa.fromfile) & cond_a:
            print(" Antenna Aa frequency range is fixed, you cannot change frequency")
        elif (self.Ab.fromfile) & cond_b:
            print(" Antenna Ab frequency range is fixed,you cannot change frequency")
        else:
            self._fGHz = freq
            self.Aa.fGHz=self.fGHz
            self.Ab.fGHz=self.fGHz


        if hasattr(self,'ca') and hasattr(self,'cb'):
            self.checkh5()

        # if self.Aa.typ == 'Omni':
        #     self.Aa.fGHz = self.fGHz
        # if self.Ab.typ == 'Omni':
        #     self.Ab.fGHz = self.fGHz
        #if len(freq)>1:
        #    self.fmin = freq[0]
        #    self.fmax = freq[-1]
        #    self.fstep = freq[1]-freq[0]
        #else:
        #    self.fmin = freq
        #    self.fmax = freq
        #    self.step = 0

    @wav.setter
    def wav(self,waveform):
        self._wav = waveform
        if 'H' in dir(self):
            if len(self.H.taud[0])!=0:
                self.chanreal = self.H.get_cir(self.wav.sfg)


    def __repr__(self):

        if hasattr(self,'filename'):
            s = 'filename: ' + self.filename +'\n'
            #s = s + 'Layout file: ' + self.Lname + '\n'
            s = s + 'Node A   \n'
            s = s + '------  \n'
            s = s + '\tcoord : ' + str (self.a) + "  in cycle " + str(self.ca) + '\n'
            s = s + '\tantenna type: ' + str (self.Aa.typ) + '\n'
            if (self.Ta!=np.eye(3)).any():
                s = s + '\trotation matrice : \n ' + str (self.Ta) + '\n\n'
            s = s + 'Node B   \n'
            s = s + '------  \n'
            s = s + '\tcoord : ' + str (self.b) + " in cycle " + str(self.cb) + '\n'
            s = s + '\tantenna : ' + str (self.Ab.typ) + '\n'
            if (self.Ta!=np.eye(3)).any():
                s = s + 'rotation matrice : \n ' + str (self.Tb) + '\n\n'
            s = s + '---------------- \n'
            s = s + 'distance AB : ' + str("%6.3f" % np.sqrt(np.sum((self.a-self.b)**2))) + ' m \n'
            s = s + 'delay AB : ' + str("%6.3f" % (np.sqrt(np.sum((self.a-self.b)**2))/0.3)) + ' ns\n'
            rd2deg = 180/np.pi
            if not np.allclose(self.a,self.b):
                vsba = self.b-self.a
                a1 = geu.angledir(vsba[None,:])
                a2 = geu.angledir(-vsba[None,:])
                s = s + 'azimuth (A | B) : %.2f ' % (a1[0,1]*rd2deg) +' deg  | %.2f' % (a2[0,1]*rd2deg) + ' deg\n'
                s = s + 'elevation (A | B) : %.2f' % (a1[0,0]*rd2deg) + ' deg | %.2f ' % (a2[0,0]*rd2deg) + ' deg\n'
                s = s + 'tilt (A |  B) : '+str((a1[0,0]-np.pi/2)*rd2deg)+ ' deg  | '+ str((a2[0,0]-np.pi/2)*rd2deg)+ ' deg\n'
            #s = s + 'Frequency range :  \n'
            s = s + '------------- \n'
            Nf = len(self.fGHz)
            s = s + 'fGHz : %.2f, %.2f, %g ' %(self.fGHz[0],self.fGHz[-1],Nf) +'\n'
            if Nf>1:
                s = s + 'fstep (GHz) : ' + str(self.fGHz[1]-self.fGHz[0]) +'\n'
            d =  np.sqrt(np.sum((self.a-self.b)**2))
            if Nf>1:
                fcGHz = (self.fGHz[-1]+self.fGHz[0])/2.
            else:
                fcGHz = self.fGHz[0]
            L  = 32.4+20*np.log(d)+20*np.log10(fcGHz)
            s = s + '------------- \n'
            s = s + 'cutoff/threshold : %g / %.2f' %(self.cutoff, self.threshold)+'\n'
            s = s + 'max delay /dist: %.2f ns / %.2f m' %(self.delay_excess_max_ns,self.delay_excess_max_ns*0.3)+'\n'
            s = s + '-------------- \n'
            if hasattr(self,'Si'):
                s = s + '# Si : ' + str(len(self.Si))
            if hasattr(self,'r2d'):
                s = s + '\n# r2d : ' + str(len(self.r2d))
            if hasattr(self,'R'):
                s = s + '\n# R : ' + str(len(self.R))
            if hasattr(self,'C'):
                s = s + '\n# C.Ctt.y : ' + str(self.C.Ctt.y.shape)
                s = s + '\n# C.Ctp.y : ' + str(self.C.Ctp.y.shape)
                s = s + '\n# C.Cpt.y : ' + str(self.C.Cpt.y.shape)
                s = s + '\n# C.Cpp.y : ' + str(self.C.Cpp.y.shape)
        else:
            s = 'No Layout specified'
        return s

    def inforay(self,iray):
        """ provide full information about a specified ray

        Parameters
        ----------

        iray : int
            ray index

        """
        print("Ray : "+str(iray))
        if not self.R.evaluated:
            self.R.eval()

        PM = self.R.info(iray,ifGHz=0,matrix=1)
        print("Propagation Channel 2x2 (C):")
        self.C.inforay(iray)

        if self.C.islocal:
            # back to global frame
            self.C.locbas()
            dist = self.C.tauk[iray]*0.3
            C = dist*np.array([[self.C.Ctt.y[iray,0],self.C.Ctp.y[iray,0]],
                               [self.C.Cpt.y[iray,0],self.C.Cpt.y[iray,0]]] )
            b = np.allclose(PM,C)
            self.C.locbas()

    # def initfreq(self):
    #     """ Automatic freq determination from
    #         Antennas
    #     """
    #     #sf = self.fGHz[1]-self.fGHz[0]
    #     sf = 1e15
    #     if hasattr(self.Aa,'fGHz'):
    #         fa = self.Aa.fGHz
    #         if len(fa)==0:
    #             fa = np.array([2.4])
    #             self.Aa.fGHz = fa
    #             # raise AttributeError("Incompatible frequency range in Antenna. Consider change Dlink.fGHz")
    #             print "Incompatible frequency range in Antenna. WARNING  Dlink.fGHz changed to 2.4GHz"
    #         try:
    #             sa = fa[1]-fa[0]  # step frequecy
    #         except: #single frequency
    #             sa = fa[0]
    #         # step
    #         if len(self.fGHz)>0:
    #             minfa = max(min(fa),min(self.fGHz))
    #             maxfa = min(max(fa),max(self.fGHz))
    #         else:
    #             minfa = min(fa)
    #             maxfa = max(fa)
    #         sf = min(sa,sf)
    #         self.fGHz = np.arange(minfa,maxfa+sf,sf)

    #     elif hasattr(self.Ab,'fGHz'):
    #         fb = self.Ab.fGHz
    #         if len(fb)==0:
    #             # raise AttributeError("Incompatible frequency range in Antenna. Consider change Dlink.fGHz")
    #             fb = np.array([2.4])
    #             self.Ab.fGHz=fb
    #             # raise AttributeError("Incompatible frequency range in Antenna. Consider change Dlink.fGHz") 
    #             print "Incompatible frequency range in Antenna. WARNING  Dlink.fGHz changed to 2.4GHz"
    #         try:
    #             sb = fb[1]-fb[0] # step frequency 
    #         except:
    #             sb = fb[0]

    #         if len(self.fGHz)>0:
    #             minfb = max(min(self.fGHz),min(fb))
    #             maxfb = min(max(self.fGHz),max(fb))
    #         else:
    #             minfb = min(fb)
    #             maxfb = max(fb)

    #         sf = min(sf,sb)
    #         self.fGHz = np.arange(minfb,maxfb+sf,sf)
    #     else:
    #         self.fGHz = np.array([2.3,2.4,2.5])


    def init_positions(self,force=False):
        """ initialize random positions for a link

        Parameters
        ----------

        force : boolean 

        """
        ###########
        # init pos & cycles
        #
        # If a and b are not specified
        #  they are chosen as center of gravity of cycle 0
        #
        ###########
        nodes = self.L.Gt.nodes()
        #
        # pick the point outside building if Layout.indoor not activated 
        #
        if self.L.typ=='outdoor':
            nodes = [n for n in nodes if n!=0 and not self.L.Gt.node[n]['indoor']]
        elif self.L.typ=='indoor':
            nodes = [n for n in nodes if n!=0 and self.L.Gt.node[n]['indoor']] 
        else:
            nodes = [n for n in nodes if n!=0 ]

        # draw the link extremities randomly

        np.random.seed(self.seed)
        ia = np.random.randint(0,len(nodes))
        ib = np.random.randint(0,len(nodes))
        if len(self.a)==0 or force:
            self.ca = nodes[ia]
        else:
            if len(self.a) ==2:
                a=np.r_[self.a,1.0]
            else:
                a=self.a
            self.ca = self.L.pt2cy(a)
            self.a = a

        if len(self.b)==0 or force:
            self.cb = nodes[ib]
        else:
            if len(self.b) ==2:
                b=np.r_[self.b,1.0]
            else:
                b=self.b
            self.cb = self.L.pt2cy(b)
            self.b = b

    def reset_config(self):
        """ reset configuration when a new layout is loaded
        """
        try:
            self.L.dumpr()
        except:
            self.L.build()
            self.L.dumpw()

        # self.a = self.L.cy2pt(self.ca)
        # self.b = self.L.cy2pt(self.cb)

        # change h5py file if layout changed
        self.filename = 'Links_' + str(self.save_idx) + '_' + self._Lname + '.h5'
        filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])
        if not os.path.exists(filenameh5) :
            print('Links save file for ' + self.L._filename + ' does not exist.')
            print('It is beeing created. You\'ll see that message only once per Layout')
            self.save_init(filenameh5)


        self.ca = 1
        self.cb = 1
        self.init_positions(force=True)

        try:
            delattr(self,'Si')
        except:
            pass
        try:
            delattr(self,'R')
        except:
            pass
        try:
            delattr(self,'C')
        except:
            pass
        try:
            delattr(self,'H')
        except:
            pass


    def checkh5(self):
        """ check existence of previous simulations run with the same parameters.


        Returns
        -------

        update self.dexist dictionnary

        """
        # get identifier group name in h5py file
        self.get_grpname()
        # check if group name exists in the h5py file
        [self.check_grpname(k,self.dexist[k]['grpname'])   for k in self.save_opt]



    def save_init(self,filename_long):
        """ initialize the hdf5 file for link saving

        Parameters
        ----------

        filename_long : str
            complete path and filename

        'sig'    : Signatures
        'ray2'   : 2D rays
        'ray'    : 3D rays
        'Ct'     : Propagation channel
        'H'      : Transmission channel
        'p_map'  : points
        'c_map'  : cycles
        'f_map'  : frequency
        'A_map'  : antennas
        'T_map'  : rotation

        """


        f=h5py.File(filename_long,'w')
        # try/except to avoid loosing the h5 file if
        # read/write error

        try:

            f.create_group('sig')
            f.create_group('ray2')
            f.create_group('ray')
            f.create_group('Ct')
            f.create_group('H')
            # mapping point a
            f.create_dataset('p_map',shape=(0,3), maxshape=(None,3),dtype='float64')
            # mapping cycles
            f.create_dataset('c_map',shape=(0,3), maxshape=(None,3),dtype='int')
            # mapping (fmin,fmax,fstep)
            f.create_dataset('f_map',shape=(0,3), maxshape=(None,3),dtype='float64')
            # mapping Antenna name
            f.create_dataset('A_map',shape=(0,1), maxshape=(None,1),dtype="S10")
            # mapping rotation matrices Antenna
            f.create_dataset('T_map',shape=(0,3,3), maxshape=(None,3,3),dtype='float64')
            f.close()
        except:
            f.close()
            raise NameError('Links: issue when initializing h5py file')

    def stack(self,key,array):
        """ stack new array in h5py file
            for a given key (dataframe/group)

        Parameters
        ----------

        key : string

        array : np.ndarray

        Returns
        -------

        idx : int
            indice of last element of the array of key

        """
        try :
            lfilename=pyu.getlong(self.filename,pstruc['DIRLNK'])
            f=h5py.File(lfilename,'a')

            if type(array)==str:
                array = array.encode('utf-8')

            if key != 'T_map':
                sc = f[key].shape
                f[key].resize((sc[0]+1,sc[1]))
                f[key][-1,:] = array
            else:
                sc = f[key].shape
                f[key].resize((sc[0]+1,sc[1],sc[2]))
                f[key][-1,:,:] = array
            f.close()
            return np.array([sc[0]])
        except:
            f.close()
            raise NameError('Link stack: issue during stacking')

    def _delete(self,key,grpname):
        """ Delete a key and associated data into h5py file

        Parameters
        ----------

        key : string
            key of the h5py file
        grpname : string
            groupe name of the h5py file

        """
        lfilename = pyu.getlong(self.filename,pstruc['DIRLNK'])
        f = h5py.File(lfilename,'a')
        # try/except to avoid loosing the h5 file if
        # read/write error

        try:
            del f[key][grpname]
            # print 'delete ',key , ' in ', grpname
            f.close()
        except:
            f.close()
            raise NameError('Link._delete: issue when deleting in h5py file')



    def save(self,obj,key,grpname,force=False):
        """ Save a given object in the correct group

        Parameters
        ----------

        obj : Object
            (Signatures|Rays|Ctilde|Tchannel)
        key : string
            key of the h5py file
        gpname : string
            groupe name of the h5py file
        force : boolean or list
        """

        if not force :
            obj._saveh5(self.filename,grpname)
        # if save is forced, previous existing data are removed and
        # replaced by new ones.
        else :
            if self.dexist[key]['exist']:
                self._delete(key,grpname)

            #if type(grpname)==str:
            #    grpname.encode('utf-8')
            #print(key,grpname)
            #if key=='ray':
            #    pdb.set_trace()
            obj._saveh5(self.filename,grpname)


        logger.debug(str(obj.__class__).split('.')[-1] + ' from '+ grpname + ' saved')


    def load(self,obj,grpname,**kwargs):
        """ Load a given object in the correct grp

        Parameters
        ----------

        obj : Object
            (Signatures|Rays|Ctilde|Tchannel)
        grpname : string
            group name of the h5py file

        kwargs :
        layout for sig and rays


        """

        obj._loadh5(self.filename,grpname,**kwargs)
        logger.debug(str(obj.__class__).split('.')[-1] + ' from '+ grpname + ' loaded')


    def get_grpname(self):
        """ Determine the data group name for the given configuration

        Notes
        -----

        Update the key grpname of self.dexist[key] dictionnary,
        where key  = 'sig'|'ray'|'Ct'|'H'

        """

        ############
        # Signatures
        ############

        array = np.array(([self.ca,self.cb,self.cutoff]))
        ua_opt, ua = self.get_idx('c_map',array)
        th = str(int(np.round(self.threshold,decimals=2)*100))


        grpname = str(self.ca) + '_' +str(self.cb) + '_' + str(self.cutoff) + '_' + th
        self.dexist['sig']['grpname']=grpname



        ############
        # Rays
        #############

        # check existence of self.a in h5py file

        ua_opt, ua = self.get_idx('p_map',self.a)
        # check existence of self.b in h5py file
        ub_opt, ub = self.get_idx('p_map',self.b)
        # Write in h5py if no prior a-b link

        grpname = str(self.cutoff) + '_' + th + '_' + str(ua) + '_' +str(ub)
        self.dexist['ray2']['grpname'] = grpname
        self.dexist['ray']['grpname'] = grpname



        ############
        # Ctilde
        #############

        # check existence of frequency in h5py file
        #farray = np.array(([self.fmin,self.fmax,self.fstep]))
        Nf = len(self.fGHz)
        if Nf > 1:
            farray = np.array(([self.fGHz[0], self.fGHz[-1], self.fGHz[1]-self.fGHz[0]]))
        else:
            farray = np.array(([self.fGHz[0], self.fGHz[-1],0]))
        uf_opt, uf = self.get_idx('f_map', farray)

        grpname = str(ua) + '_' + str(ub) + '_' + str(uf)
        self.dexist['Ct']['grpname'] = grpname

        ############
        # H
        #############


        # check existence of Rot a (Ta) in h5py file
        uTa_opt, uTa = self.get_idx('T_map',self.Ta)
        # check existence of Rot b (Tb) in h5py file
        uTb_opt, uTb = self.get_idx('T_map',self.Tb)
        # check existence of Antenna a (Aa) in h5py file
        #uAa_opt, uAa = self.get_idx('A_map',self.Aa._filename)
        uAa_opt, uAa = self.get_idx('A_map',self.Aa.typ)
        # check existence of Antenna b (Ab) in h5py file
        uAb_opt, uAb = self.get_idx('A_map',self.Ab.typ)


        grpname = str(ua) + '_' + str(ub) + '_' + str(uf) + \
                  '_'  + str(uTa) + '_' + str(uTb) + \
                  '_'  + str(uAa) + '_' + str(uAb)

        self.dexist['H']['grpname'] = grpname


    def check_grpname(self,key,grpname):
        """Check if the key's data with a given groupname
            already exists in the h5py file

        Parameters
        ----------

        key: string
            key of the h5py group
        grpname : string
            groupe name of the h5py file

        Notes
        -----

        update the key grpname of self.dexist[key] dictionnary

        """
        try :
            lfilename=pyu.getlong(self.filename,pstruc['DIRLNK'])
            f = h5py.File(lfilename,'r')
            if grpname.encode('utf8') in f[key].keys():
                self.dexist[key]['exist'] = True
            else :
                self.dexist[key]['exist'] = False
            f.close()
        except:
            f.close()
            raise NameError('Link exist: issue during stacking')


    def get_idx(self,key,array,tol=1e-3):
        """ get the index of the requested array in the group key
            of the hdf5 file.
            If array doesn't exist, the hdf5file[key] array is stacked


        Parameters
        ----------

        key: string
            key of the h5py group
        array : np.ndarray
            array to check existence
        tol : np.float64
            tolerance (in meters for key == 'p_map')

        Returns
        -------

        (u_opt, u): tuple

        u : np.ndarray
            the index in the array of the file[key] group
        u_opt : string ('r'|'s')
            return 'r' if array has been read into h5py file
            return 's' if array has been stacked into the array of group key



        See Also:
        --------

        Links.array_exist

        """

        umap = self.array_exist(key,array,tol=tol)
        lu = len(umap)
        # if exists take the existing one
        # otherwise value is created
        if lu != 0:
            u = umap
            u_opt='r'
        else :
            u = self.stack(key,array)
            u_opt='s'
        return u_opt,u[0]


    def array_exist(self,key,array,tol=1e-3) :
        """ check an array key has already been stored in h5py file


        Parameters
        ----------

        key: string
            key of the h5py group
        array : np.ndarray
            array type to check existency
        tol : np.float64
            tolerance (in meter for key == 'p_map')

        Returns
        -------

        (ua)

        ua : np.ndarray
            the indice in the array of the file[key] group
            if the array is emtpy, value doesn't exist

        TODO
        ----

        Add a tolerance on the rotation angle (T_map)

        """

        lfilename=pyu.getlong(self.filename,pstruc['DIRLNK'])
        try :
            f=h5py.File(lfilename,'a')
            fa = f[key][...]
            f.close()
        except:
            f.close()
            raise NameError('Link check_exist: issue during reading')

        if key == 'c_map':
            eq = array == fa
            # sum eq = 3 means cy0,cy1 and cutoff are the same in fa and array
            ua = np.where(np.sum(eq,axis=1)==3)[0]

        elif key == 'p_map':

            da = np.sqrt(np.sum((array-fa)**2,axis=1))
            # indice points candidate in db for a

            ua = np.where(da<tol)[0]

        elif key == 'f_map':
            # import ipdb
            # ipdb.set_trace()
            #### fmin_h5 < fmin_rqst
            ufmi = fa[:,0]<=array[0]
            # old version
                # ufmi = np.where(fa[:,0]<=array[0])[0]
                # lufmi = len(ufmi)

            #### fmax_h5 > fmax_rqst
            ufma = fa[:,1]>=array[1]
            # old version
                # ufma = np.where(fa[:,1]>=array[1])[0]
                # lufma = len(ufma)

            ### fstep_h5 < fstep_rqst
            ufst = fa[:,2]<=array[2]
            # old version
                # ufst = np.where(fa[:,2]<=array[2])[0]
                # lufst = len(ufst)

            # if fmin, fmax or fstep
            #if (lufmi==0) and (lufma==0) and (lufst==0):
            if (not ufmi.any()) and (not ufma.any()):
                ua = np.array([])
            else:
                # find common lines of fmin and fmax and fstep
                ua = np.where(ufmi & ufma & ufst)[0]
                # ua = np.where(np.in1d(ufmi,ufma,ufst))[0]
                # # find common lines of fmin and fmax
                # ufmima = np.where(np.in1d(ufmi,ufma))[0]
                # # find common lines of fmin, fmax and fstep
                # ua = np.where(np.in1d(ufmima,ufst))[0]

        elif key == 'A_map':
            ua = np.where(fa == array.encode('utf-8'))[0]

        elif key == 'T_map':
            eq = array == fa
            seq = np.sum(np.sum(eq,axis=1), axis=1)
            ua = np.where(seq == 9)[0]
        else :
            raise NameError('Link.array_exist : invalid key')

        return ua

    def evalH(self,**kwargs):
        """ evaluate channel transfer function

        Notes
        -----

        This function modifies the orientation of the antenna at both sides
        via Ta and Tb 3x3 matrices and recalculates the channel transfer function
        for those new orientations.
        The self.H variable is updated

        """
        # Antenna Rotation
        self.C.locbas(Ta=self.Ta, Tb=self.Tb)
        # Transmission channel calculation
        H = self.C.prop2tran(a=self.Aa, b=self.Ab, Friis=True, debug=True)
        self.H = H

    def eval(self,**kwargs):
        """ link evaluation

        Parameters
        ----------

        applywav :boolean
         Apply waveform to H
        force : list
            Force the computation (['sig','ray2','ray,'Ct','H']) AND save (replace previous computations)
        alg : 1|'old'|'exp'|'exp2'
            version of run for signature
        si_progress: bollean ( False)
            display progression bar for signatures
        diffraction : boolean (False)
            takes into consideration diffraction points
        ra_number_mirror_cf : int
            rays.to3D number of ceil/floor reflexions
        ra_ceil_H: float, (default [])
            ceil height .
                If [] : Layout max ceil height
                If 0 : only floor reflection (outdoor case)
                If -1 : neither ceil nor floor reflection (2D case)
        ra_vectorized: boolean (True)
            if True used the (2015 new) vectorized approach to determine 2drays
        progressbar: str
            None: no progress bar
            python : progress bar in ipython


        Notes
        -----

        update self.ak and self.tk

        self.ak : ndarray
            alpha_k
        self.tk : ndarray
            tau_k


        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.simul.link import *
            >>> DL=DLink(L="defstr.lay")
            >>> DL.eval()
            >>> DL.show()
            >>> DL.plt_cir()
            >>> plt.show()


        See Also
        --------

        pylayers.antprop.signature
        pylayers.antprop.rays


        """

        defaults = {'applywav': False,
                   'si_progress': True,
                   'diffraction': True,
                   'ra_vectorized': True,
                   'ra_ceil_H': [],
                   'ra_number_mirror_cf': 1,
                   'force': True,
                   'bt': True,
                   'si_reverb': 4,
                   'nD': 2,
                   'nR': 10,
                   'nT': 10,
                   'debug': False,
                   'progressbar': None,
                   'rm_aw': True
                   }
        # check antenna frequency range compatibility
        if (self.Aa.fGHz!=self.Ab.fGHz).all():
            raise AttributeError("Antenna frequency range are not compatible")

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if 'delay_excess_max_ns' not in kwargs:
            kwargs['delay_excess_max_ns'] = self.delay_excess_max_ns
        else:
            self.delay_excess_max_ns = kwargs['delay_excess_max_ns']

        if 'cutoff' not in kwargs:
            kwargs['cutoff'] = self.cutoff
        else:
            self.cutoff=kwargs['cutoff']

        if 'threshold' not in kwargs:
            kwargs['threshold'] = self.threshold
        else:
            self.threshold=kwargs['threshold']

        if 'force' in kwargs:
            if not isinstance(kwargs['force'],list):
                if kwargs['force'] == True :
                    kwargs['force'] = ['sig','ray2','ray','Ct','H']
                else :
                    # Ct and H are not yet saved/loaded
                    # compliantly with the given configutain
                    # their are disabled here
                    kwargs['force'] = ['Ct','H']

        # self.checkh5()

        if isinstance(kwargs['progressbar'],str):
            if kwargs['progressbar'] =='notebook':
                pbar = tqdm.tqdm_notebook(total=100)
            elif kwargs['progressbar']=='python':
                pbar = tqdm.tqdm(total=100)
        elif isinstance(kwargs['progressbar'],tqdm.tqdm):
            pbar = kwargs['progressbar']



        ############
        # Signatures
        ############


        tic = time.time()
        Si = Signatures(self.L,
                        self.ca,
                        self.cb,
                        cutoff = kwargs['cutoff'],
                        threshold = kwargs['threshold'])

        if (self.dexist['sig']['exist'] and not ('sig' in kwargs['force'])):
            logger.info(" Load existing signatures from :%s",self.dexist['sig']['grpname'])
            self.load(Si, self.dexist['sig']['grpname'], L=self.L)
        else:
            logger.info(" Run signatures")

            self.nD = kwargs['nD']
            self.nT = kwargs['nT']
            self.nR = kwargs['nR']
            self.bt = kwargs['bt']

            Si.run(cutoff = self.cutoff,
                    diffraction = kwargs['diffraction'],
                    threshold = self.threshold,
                    delay_excess_max_ns = self.delay_excess_max_ns,
                    nD = self.nD,
                    nR = self.nR,
                    nT = self.nT,
                    progress = kwargs['si_progress'],
                    bt = self.bt)

            logger.info(" Save signature in %s ",self.dexist['sig']['grpname'])
            self.save(Si,'sig',self.dexist['sig']['grpname'],force = kwargs['force'])

        self.Si = Si

        toc = time.time()
        logger.info(" End signature in %d sec",toc-tic)

        try:
            pbar.update(20)
        except:
            pass


        ############
        # Rays
        ############

        logger.info(" Start Rays determination")
        tic = time.time()
        r2d = Rays(self.a,self.b)


        #############
        # get 2D rays
        #############
        if self.dexist['ray2']['exist'] and not ('ray2' in kwargs['force']):
            logger.info(" Load r2d from %s", self.dexist['ray2']['grpname'])
            self.load(r2d,self.dexist['ray2']['grpname'], L=self.L)
        else :
            # perform computation ...
            # ... with vectorized ray evaluation
            logger.debug(" a : (%d,%d,%d)", self.a[0], self.a[1], self.a[2])
            logger.debug(" b : (%d,%d,%d)", self.b[0], self.b[1], self.b[2])

            if kwargs['ra_vectorized']:
                logger.info(" Determine r2d vectorized version")
                r2d = Si.raysv(self.a,self.b)
            # ... or with original and slow approach ( to be removed in a near future)
            else :
                logger("Determine r2d non vectorized version")
                r2d = Si.rays(self.a,self.b)
            # save 2D rays

            logger.info(" save r2d in %s ",self.dexist['ray2']['grpname'])
            self.save(r2d,'ray2', self.dexist['ray2']['grpname'], force = kwargs['force'])

        self.r2d = r2d

        #############
        # get 3D rays
        #############

        R = Rays(self.a,self.b)
        R.is3D = True
        if self.dexist['ray']['exist'] and not ('ray' in kwargs['force']):
            logger.info(" Load r3d from %s", self.dexist['ray']['grpname'])
            self.load(R,self.dexist['ray']['grpname'], L=self.L)
        else :
            if kwargs['ra_ceil_H'] == []:
                if self.L.typ=='indoor':
                    ceilheight = self.L.maxheight
                else:
                    ceilheight = 0
            else:
                ceilheight = kwargs['ra_ceil_H']


            logger.info(" Run to3d H: %d, N: %d", ceilheight,kwargs['ra_number_mirror_cf'] )
            R = self.r2d.to3D(self.L, H=ceilheight, N=kwargs['ra_number_mirror_cf'])
            if kwargs['rm_aw']:
                R = R.remove_aw(self.L)


            logger.info(" Run R.locbas ")
            R.locbas(self.L)

            logger.info(" Run R.fillinter ")
            R.fillinter(self.L)

            # C = Ctilde()

            # C = R.eval(self.fGHz)

            # save 3D rays

            logger.info(" Save 3D rays in %s",self.dexist['ray']['grpname'])
            self.save(R, 'ray', self.dexist['ray']['grpname'], force = kwargs['force'])


        self.R = R
        toc = time.time()
        logger.info(" Stop rays %d",toc-tic)

        if self.R.nray == 0:
            raise NameError('No rays have been found. Try to re-run the simulation with a higher S.cutoff ')
        try:
            pbar.update(20)
        except:
            pass
        ############
        # Ctilde
        ############

        if self.dexist['Ct']['exist'] and not ('Ct' in kwargs['force']):
            C = Ctilde()
            self.load(C,self.dexist['Ct']['grpname'])
        else :
            #if not hasattr(R,'I'):
            # Ctilde...
            # Find an other criteria in order to decide if the R has
            # already been evaluated

            C = R.eval(self.fGHz)
            # ...save Ct
            self.save(C,'Ct',self.dexist['Ct']['grpname'],force = kwargs['force'])

        self.C = C

        try:
            pbar.update(20)
        except:
            pass
        ############
        # H
        ############

        H = Tchannel()

        if self.dexist['H']['exist'] and not ('H' in kwargs['force']):
            self.load(H,self.dexist['H']['grpname'])
        else :
            # Ctilde antenna
            C.locbas(Ta=self.Ta, Tb=self.Tb)
            #T channel
            H = C.prop2tran(a=self.Aa,b=self.Ab,Friis=True,debug=kwargs['debug'])
            self.save(H,'H',self.dexist['H']['grpname'],force = kwargs['force'])
        self.H = H
        try:
            pbar.update(20)
        except:
            pass

        if kwargs['applywav']:
            if self.H.isFriis:
                self.ir = self.H.get_cir(self.wav.sf)
            else:
                self.ir = self.H.get_cir(self.wav.sfg)
        try:
            pbar.update(20)
        except:
            pass
        self.checkh5()


    def adp(self,imax=1000):
        """ construct the angular delay profile

        Parameters
        ----------
        imax : int
        """
        self.adpmes = self.afpmes.toadp()
        self.adpmes.cut(imax=imax)
        self.adprt = self.afprt.toadp()
        self.adprt.cut(imax=imax)

    def afp(self,**kwargs):
        """ Evaluate angular frequency profile 

        Parameters
        ----------

        fGHz  : np.array 
            frequency range 
        az : azimuth angle (radian)  
        tilt : tilt angle (-pi/2<tilt<pi/2) 
        polar : string
        win : string 'rect' | 'hamming'
        _filemeas : string 
        _filecal : string 
        ang_offset : 
        BW : float 
            bandwidth
        ext : string 
            'txt' | 'mat'
        dirmeas : string 
            directory of the data in the project path 

        Notes
        -----

        If a measurement file is given the angular range is obtained from the measurement
        otherwise the variable az is used. 

        """
        defaults = {'fGHz':32.6,
                     'az': 0,
                     'tilt':0,
                     'polar':'V',
                     'win':'rect',
                     '_filemeas':'',
                     '_filecal':'',
                     'ang_offset' : 0.37,
                     'BW': 1.6,
                     'ext':'txt',
                     'dirmeas':'meas',
                     'refinement':False
                    }
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k] 

        fGHz = kwargs.pop('fGHz')
        az   = kwargs.pop('az')
        tilt = kwargs.pop('tilt')
        polar = kwargs.pop('polar')
        win = kwargs.pop('win')               # argument for loadmes
        _filemeas = kwargs.pop('_filemeas')   # argument for loadmes
        _filecal = kwargs.pop('_filecal')     # argument for loadmes
        dirmeas = kwargs.pop('dirmeas')      # argument for loadmes
        ang_offset = kwargs.pop('ang_offset') # argument for loadmes
        BW = kwargs.pop('BW')                 # argument for loadmes
        ext = kwargs.pop('ext')               # argument for loadmes
        refinement = kwargs.pop('refinement')  # argument for loadmes

        # read measurement if available
        if _filemeas!='':
            fcGHz = self.fGHz[0]
            self.afpmes = AFPchannel(tx=self.a,rx=self.b)
            self.afpmes.loadmes(_filemeas,
                                _filecal,
                                 fcGHz=fcGHz,
                                 BW=BW,
                                 win=win,
                                 ang_offset=ang_offset,
                                 ext=ext,
                                 dirmeas=dirmeas,
                                 refinement=refinement)
            az = self.afpmes.az
            #
            # afpmes.x
            # afpmes.y
            # afpmes.fcGHz
            # afpmes.az     measure angular range 
            # afpmes.azrt   ray tracing angular range  
        # create an empty AFP 
        # tx = a
        # rx = b 
        # angular range (a) : phi 
        #

        self.afprt = AFPchannel(tx=self.a,rx=self.b,az=az)
        for k,ph in enumerate(az.squeeze()):
            self.Tb = geu.MATP(self.Ab.sl,self.Ab.el,ph,tilt,polar)
            # self._update_show3(ant='b')
            # pdb.set_trace()
            self.evalH()
            E = self.H.energy()
            if k==0:
                self.dpadp = E[None,...]
            else:
                self.dpadp = np.concatenate((self.dpadp,E[None,...]),axis=0)

            if self.H.y.shape[3]!=1:
                S = np.sum(self.H.y*np.exp(-2*1j*np.pi*self.H.x[None,None,None,:]*self.H.taud[:,None,None,None]),axis=0)
            else:
                S = np.sum(self.H.y*np.exp(-2*1j*np.pi*fGHz*self.H.taud[:,None,None,None]),axis=0)

            try:
                self.afprt.y = np.vstack((self.afprt.y,np.squeeze(S)))
            except:
                self.afprt.y = np.squeeze(S)
        if self.H.y.shape[3]!=1:
            self.afprt.x = self.H.x
            self.afprt.fcGHz = self.afprt.x[len(self.afprt.x)/2]
        else:
            self.afprt.x = fGHz
            self.afprt.fcGHz = fGHz[len(fGHz)/2]

        self.dpdp = np.sum(self.dpadp,axis=0)
        self.dpap = np.sum(self.dpadp,axis=1)
        #return(afp)


    def select(self):
        fig,ax = self.show()
        self.cid = fig.canvas.mpl_connect('button_press_event',self.OnClick)
        return(fig,ax)

    def OnClick(self,event):
        x = event.xdata
        y = event.ydata
        if event.button==1:
            self.a=np.array([x,y,1.2])
            self.caf.set_offsets(np.array([[x,y]]))
            plt.draw()
        if event.button==3:
            self.b=np.array([x,y,1.2])
            self.cbf.set_offsets(np.array([[x,y]]))
            plt.draw()
        if event.button==2:
            self.eval()
            self._show3()
        print (x,y)

    def show(self,**kwargs):
        """ show the link

        Parameters
        ----------

        s   : int
            size of Tx/Rx marker in points
        ca  : string
            color of termination a (A)
        cb  : string
            color of termination b (B)
        alpha : float
            marker transparency (0 < alpha <1)
        axis : boolean
            display axis boolean (default True)
        figsize : tuple
            figure size if fig not specified default (20,10)
        fontsize : int
            default 20
        rays : boolean
            activation of rays vizalization (True)
        bsig : boolean
            activation of signature vizualization (False)
        bsave : boolean
            save in a file indexed by ix
        laddr : list
            list of signature addresses
        cmap : colormap
        radius : float
            radius in meters for layout vizualization
        labels : boolean
            enabling edge label (useful for signature identification)
        pol : string
            'tt','pp','tp','pt','co','cross',tot'
        col : string
            'cmap'
        width : float
        alpha : float
        dB    : boolean
            default False
        dyn : float
            dynamic in dB (def 70dB)

        Returns
        -------

        fig,ax

        Examples
        --------

        >>> from pylayers.simul.link import *
        >>> DL=Link()
        >>> DL.show(lr=-1,rays=True,dB=True,col='cmap',cmap = plt.cm.jet)
        >>> DL.show(laddr=[(6,2)],bsig=True)

        """
        defaults ={'s': 300,   # size points
                   'ca': '#6666ff', # color a
                   'cb': '#ff0000', # color b
                   'markera': "^", # tri_up
                   'markerb': "o", # circle
                   'alpha': 1,
                   'axis': True,
                   'lr': -1,
                   'ls': -1,
                   'fig': [],
                   'ax': [],
                   'figsize': (20,10),
                   'fontsize': 20,
                   'rays': True,
                   'bsig': False,
                   'bsave': False,
                   'laddr': [(1,0)],
                   'cmap': plt.cm.hot_r,
                   'pol': 'tot',
                   'color': 'k',
                   'linewidth': 1,
                   'alpha': 1,
                   'radius': -1,
                   'vmin': [],
                   'vmax': [],
                   'dB': True,
                   'labels': False,
                   'aw': False,
                   'dyn': 70,
                   'ix' : 0,
                   'vmin': -120,
                   'vmax': -40,
                   'bcolorbar': True}

        for key in defaults:
            if key not in kwargs:
                kwargs[key]=defaults[key]

        if kwargs['fig']==[]:
            fig = plt.figure(figsize=kwargs['figsize'])
        else:
            fig = kwargs['fig']

        if kwargs['ax']==[]:
            ax = fig.add_subplot(111)
        else:
            ax=kwargs['ax']

        #
        # Layout
        #

        fig, ax = self.L.showG('s',
                              nodes = False,
                              fig = fig,
                              ax = ax,
                              labels = kwargs['labels'],
                              aw = kwargs['aw'],
                              axis = kwargs['axis'])

        if kwargs['radius'] == -1:
            kwargs['radius'] = self.L.radius

        # background color
        #ax.set_facecolor('#cccccc')
        #
        # Plot Rays
        #
        if kwargs['rays'] and self.R.nray > 0:
            #ECtt,ECpp,ECtp,ECpt = self.C.energy()
            #if kwargs['pol']=='tt':
            #    val = ECtt
            #if kwargs['pol']=='pp':
            #    val = ECpp
            #if kwargs['pol']=='tp':
            #    val = ECtp
            #if kwargs['pol']=='pt':
            #    val = ECpt
            #if kwargs['pol']=='tot':
            #    val = ECtt+ECpp+ECpt+ECtp
            #if kwargs['pol']=='co':
            #    val = ECtt+ECpp
            #if kwargs['pol']=='cross':
            #"    val = ECtp+ECpt
            val = self.H.energy()[:,0,0]

            clm = kwargs['cmap']
            #
            # Select group of interactions
            #
            if ((type(kwargs['lr']) is list) or 
               (type(kwargs['lr']) is np.ndarray)):
                lr = kwargs['lr']
            else:
                if kwargs['lr']==-1:
                    lr  = np.arange(self.R.nray)
                else:
                    lr = [ int(kwargs['lr']) ]


            #
            #  Set the min and max of ray level
            #

            if kwargs['vmin'] == []:
                vmin = val.min()
                vmax = val.max()

                if kwargs['dB']:
                    vmin = 10*np.log10(vmin)
                    vmax = 10*np.log10(vmax)
            else:
                vmin = kwargs['vmin']
                vmax = kwargs['vmax']


            #
            # limitation of the vizualization zone around the center of the link
            #

            pm = (self.a + self.b)/2.
            R  = np.minimum(kwargs['radius'],1.5*self.L.radius)
            #ax.set_xlim(pm[0]-R,pm[0]+R)
            #ax.set_ylim(pm[1]-R,pm[1]+R)

            #
            # each ray ir from list lr has its own color
            #
            for ir  in lr:
                if kwargs['dB']:
                    valdB = np.array(10*np.log10(val[ir]))
                    valdB = np.maximum(vmin,valdB)
                    valdB = np.minimum(vmax,valdB)
                    #RayEnergy = max((10*np.log10(val[ir]/val.max())+kwargs['dyn']),0)/kwargs['dyn']
                    RayEnergy = (valdB-vmin)/(vmax-vmin)
                else:
                    valLin = val[ir]
                    valdB = np.maximum(vmin,valLin)
                    valdB = np.minimum(vmax,valLin)
                    RayEnergy = (valLin-vmin)/(vmax - vmin)

                if kwargs['color'] == 'cmap':
                    color = clm(RayEnergy)
                    #width = 10*RayEnergy
                    linewidth = kwargs['linewidth']
                    alpha = 1
                else:
                    color = kwargs['color']
                    linewidth = kwargs['linewidth']
                    alpha = kwargs['alpha']

                # plot ray (i,r)
                fig,ax = self.R.show(rlist = [ir],
                               color = color,
                               linewidth = 10*RayEnergy,
                               alpha = alpha,
                               fig = fig, ax = ax,
                               layout = False,
                               points = False,
                               bcolorbar = kwargs['bcolorbar'],
                               cmap = kwargs['cmap'],
                               vmin = vmin,
                               vmax = vmax )

            #if kwargs['color']=='cmap':
            #    sm = plt.cm.ScalarMappable(cmap=kwargs['cmap'], norm=plt.Normalize(vmin=kwargs['vmin'],
            #                                                  vmax=kwargs['vmax']))
            #    sm._A = []
            #    cb = plt.colorbar(sm)
            #    cb.ax.tick_params(labelsize=24)
            #    cb.set_label('Level (dB)', fontsize=24)
        #
        # Plot signature
        #

        if kwargs['bsig']:
            for addr in kwargs['laddr']:
                seq = self.Si[addr[0]][2*addr[1]:2*addr[1]+2,:]
                Si = Signature(seq)
                isvalid,r,u = Si.sig2ray(self.L,self.a[0:2],self.b[0:2])
                fig,ax = Si.show(self.L,self.a[0:2],self.b[0:2],fig=fig,ax=ax)

        #
        # Point A
        #
        self.caf = ax.scatter(self.a[0], self.a[1],
                              c = kwargs['ca'],
                              s = kwargs['s'],
                              marker = kwargs['markera'],
                              edgecolor='black',
                              facecolor=kwargs['ca'],
                              linewidth=2,
                              alpha =kwargs['alpha'],
                              zorder = 1000)

        #ax.text(self.a[0]+0.3,self.a[1]+0.3,'a',
        #        fontsize = kwargs['fontsize'], bbox=dict(facecolor='white',alpha=0.5))
        #
        # Point B
        #
        self.cbf = ax.scatter(self.b[0], self.b[1],
                              c = kwargs['cb'],
                              s = kwargs['s'],
                              marker = kwargs['markerb'],
                              edgecolor='black',
                              facecolor=kwargs['cb'],
                              linewidth=2,
                              alpha = kwargs['alpha'],
                              zorder = 1000)

        #ax.text(self.b[0]+0.3, self.b[1]+0.3, 'b',
        #        fontsize=kwargs['fontsize'],bbox=dict(facecolor='white',alpha=0.5))

        #
        # white scale
        #

        xe = 1
        ye = 3
        le = 1
        ax.plot(np.array([xe,xe+le]),np.array([ye,ye]),linewidth=4,color='k')
        ax.plot(np.array([xe,xe]),np.array([ye,ye+0.2]),linewidth=4,color='k')
        ax.plot(np.array([xe+le,xe+le]),np.array([ye,ye+0.2]),linewidth=4,color='k')
        ax.text(xe-0.1,ye-0.5,'1 meter',fontsize=18)
        #plt.axis('on')
        ax.tick_params(labelsize = 24)
        ax.set_xlabel('x meters',fontsize = 24)
        ax.set_ylabel('y meters',fontsize = 24)
        #plt.savefig('Link.eps')

        if kwargs['bsave']:
            plt.savefig('Link'+str(kwargs['ix'])+'.png')
            plt.close()

        plt.axis('auto')
        return fig,ax

    def _show3(self,rays=True, lay= True, ant= True, newfig= False, **kwargs):
        """ display the simulation scene using Mayavi

        Parameters
        ----------

        rays: boolean
        lay : boolean
        ant : boolean
        newfig : boolean (default : False)
        kwargs of Rays.show3()


        see also
        --------

        pylayers.gis.layout
        pylayers.antprop.antenna
        pylayers.antprop.rays

        Examples
        --------


            >>> from pylayers.simul.link import *
            >>> L=DLink()
            >>> L.eval()

        """

        if not newfig:
            self._maya_fig=mlab.gcf()
        else:
            self._maya_fig=mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))

        if 'scale' in kwargs:
            scale = kwargs.pop('scale')
        else:
            scale = 0.5

        if 'centered' in kwargs:
            centered = kwargs['centered']
        else :
            centered = False

        if centered:
            pg = np.zeros((3))
            pg[:2] = self.L.pg


        if centered:
            ptx = self.a-pg
            prx = self.b-pg
        else :
            ptx = self.a
            prx = self.b

        self._maya_fig.scene.disable_render = True


        if ant :
            Atx = self.Aa
            Arx = self.Ab
            Ttx = self.Ta
            Trx = self.Tb

            # evaluate antenna if required
            if not Atx.evaluated:
                Atx.eval()

            Atx._show3(T=Ttx.reshape(3,3),
                        po=ptx,
                        title=False,
                        bcolorbar=False,
                        bnewfig=False,
                        bcircle = False,
                        name = Atx._filename,
                        scale= scale,
                        binteract=False)

            if not Arx.evaluated:
                Arx.eval()

            Arx._show3(T=Trx.reshape(3,3),
                        po=prx,
                        title=False,
                        bcolorbar=False,
                        bnewfig=False,
                        bcircle = False,
                        name = Arx._filename,
                        scale= scale,
                        binteract=False)
        if lay:
            if self.L.typ == 'outdoor':
                show_ceil = False
                opacity = 1.
                ceil_opacity = 1.
            elif self.L.typ == 'indoor':
                show_ceil = True
                opacity = 0.7
                ceil_opacity = 0.

            self._maya_fig = self.L._show3(newfig=False,
                          opacity=opacity,
                          ceil_opacity=ceil_opacity,
                          show_ceil=show_ceil,
                          centered=centered,**kwargs)

        # mlab.text3d(self.a[0],self.a[1],self.a[2],'a',
        #             scale=1,
        #             color=(1,0,0))
        # mlab.text3d(self.b[0],self.b[1],self.b[2],'b',
        #             scale=1,
        #             color=(1,0,0))
        if rays :
            # check rays with energy
            # if hasattr(self,'H') and not kwargs.has_key('rlist'):
            #     urays = np.where(self.H.y!=0)[0]
            #     kwargs['rlist']=urays
            #     import ipdb
            #     ipdb.set_trace()

            if hasattr(self,'R'):
                if self.H.y.ndim>2:
                    ER = np.squeeze(self.H.energy())
                    kwargs['ER']=ER
                self.R._show3(L=[],**kwargs)

        fp = (self.a+self.b)/2.

        dab = np.sqrt(np.sum((self.a-self.b)**2))
        mlab.view(focalpoint=fp)#,distance=15*dab-55)
        self._maya_fig.scene.disable_render = False
        mlab.orientation_axes()
        mlab.show()
        return self._maya_fig
        #return(self._maya_fig)


    def _update_show3(self,ant='a',delrays=False):
        """
        """


        view=mlab.view()


        antenna = eval('self.A'+ant)
        rot = eval('self.T'+ant).reshape(3,3)
        pos = eval('self.'+ant)

        #if not antenna.full_evaluated:
        # if not antenna.full_evaluated:
        if not antenna.evaluated:
            antenna.eval()


        if hasattr(antenna,'_mayamesh'):
            # antenna.eval()
            x, y, z, k, scalar = antenna._computemesh(T=rot,po=pos,scale= 0.5)
            antenna._mayamesh.mlab_source.set(x=x,y=y,z=z,scalars=scalar)
        else:
            antenna._show3(T=rot,po=pos,
                title=False,
                bcolorbar=False,
                bcircle = False,
                bnewfig=False,
                scale= 0.5,
                name = antenna._filename,
                binteract=False)

        if delrays:
            import time
            for x in self._maya_fig.children[::-1]:
                if 'Rays' in x.name:
                    x.remove()
        mlab.view(view[0],view[1],view[2],view[3])
             # [x.remove() for x in self._maya_fig.children ]

        # # update wall opaccity
        # ds  =[i for i in self._maya_fig.children if self.L._filename in i.name][0]
        # a_in = self.L.Gt.node[self.ca]['indoor']
        # b_in = self.L.Gt.node[self.cb]['indoor']

        # if
        # if a_in or b_in:
        #     # indoor situation
        #     ds.children[0].children[0].actor.property.opacity=0.5
        # else:
        #     ds.children[0].children[0].actor.property.opacity=1.

    def plt_cir(self,**kwargs):
        """ plot link channel impulse response

        Parameters
        ----------

        BWGHz : Bandwidth
        Nf    : Number of frequency points
        fftshift : boolean
        rays : boolean
            display rays contributors
        fspl : boolean
            display free space path loss

        Returns
        -------

        fig,ax

        See Also
        --------

        pylayers.antprop.channel.Tchannel.getcir

        """

        defaults = {'fig' : [],
                    'ax' : [],
                    'BWGHz' :5,
                    'Nf' :1000,
                    'rays' :True,
                    'fspl' :True,
                    'vmin' :-120,
                    'vmax' : -40,
                    'taumin': 0,
                    'taumax': 160,
                    'bgrid':True,
                    'cmap':'jet',
                    'ix' : 0,
                    'bsave' : False,
                    'fontsize':18
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] == []:
            fig = plt.gcf()
        else:
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            ax = plt.gca()
        else:
            ax = kwargs['ax']

        fontsize = kwargs.pop('fontsize')
        vmin = kwargs.pop('vmin')
        vmax = kwargs.pop('vmax')
        taumin = kwargs.pop('taumin')
        taumax = kwargs.pop('taumax')
        #taumax = self.H.taud.max()
        BWGHz = kwargs['BWGHz']

        Nf = np.maximum(kwargs['Nf'],taumax*BWGHz).astype(int)
        # getcir is a Tchannel method

        self.ir = self.H.getcir(BWGHz=BWGHz, Nf=Nf)
        self.ir.plot(fig=fig, ax=ax, fontsize=fontsize)

        ax.set_ylim(vmin,vmax)


        delay = self.ir.x
        delay = delay[delay>0]
        dist = delay*0.3
        FSPL0 = -32.4- 20*np.log10(self.fGHz[0])-20*np.log10(dist)
        FSPLG = FSPL0 + self.Aa.GdBmax[0] + self.Ab.GdBmax[0]

        if kwargs['fspl']:
            # Free space path loss
            ax.plot(delay,FSPL0,linewidth=2,color='b',label='FSPL')
            # Free space path loss + gain
            ax.plot(delay,FSPLG,linewidth=3,color='k',label='FSPL+Gtmax+Grmax')

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(fontsize)

        if kwargs['bgrid']:
            ax.grid()

        if kwargs['rays']:
            # energy of each ray normaized between vmin(0) and vmax(1)
            ER = np.squeeze(self.H.energy())
            uER = ER.argsort()[::-1]
            ER = ER[uER]
            ERdB = 10*np.log10(ER)
            ERdB = np.minimum(ERdB,vmax)
            ERdB = np.maximum(ERdB,vmin)
            colors = (ERdB-vmin)/(vmax-vmin)
            #color_range = np.linspace( 0, 1., len(ER))#np.linspace( 0, np.pi, len(ER))
            # sort rays by increasing energy
            #colors = color_range[uER]
            # most important rays , it=0 ir=0 , if =0
            ax.scatter(self.H.taud[uER],ERdB,c=colors,s=200*colors,cmap=kwargs['cmap'],vmin=0,vmax=1)
            #ax.set_xlim([min(self.H.taud)-10,max(self.H.taud)+10])
            ax.set_xlim(taumin,taumax)

        ax.legend(fontsize=fontsize)
        if kwargs['bsave']:
            plt.savefig('cir'+str(kwargs['ix'])+'.png')
            plt.close()

        return fig,ax


    def plt_doa(self,**kwargs):
        """plot direction of arrival and departure

        Parameters
        ----------

        fig : plt.figure
        ax : plt.axis
        phi: tuple (-180, 180)
            phi angle
        normalize: bool
            energy normalized
        reverse : bool
            inverse theta and phi represenation
        polar : bool
            polar representation
        cmap: matplotlib.cmap
        mode: 'center' | 'mean' | 'in'
            see bsignal.energy
        s : float
            scatter dot size
        fontsize: float
        edgecolors: bool
        colorbar: bool
        title : bool

        See Also
        --------

        pylayers.antprop.channel.Tchannel.plotd

        """
        kwargs['d']='doa'
        return self.H.plotd(**kwargs)

    def plt_dod(self,**kwargs):
        """plot direction of arrival and departure

        Parameters
        ----------

        fig : plt.figure
        ax : plt.axis
        phi: tuple (-180, 180)
            phi angle
        normalize: bool
            energy normalized
        reverse : bool
            inverse theta and phi represenation
        polar : bool
            polar representation
        cmap: matplotlib.cmap
        mode: 'center' | 'mean' | 'in'
            see bsignal.energy
        s : float
            scatter dot size
        fontsize: float
        edgecolors: bool
        colorbar: bool
        title : bool

        See Also
        --------

        pylayers.antprop.channel.Tchannel.plotd

        """
        kwargs['d']='dod'
        return self.H.plotd(**kwargs)

    def plt_dspread(self,**kwargs):
        """ plot delay spread
        """
        defaults = { 'fig':[],
                     'ax':[]
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        if kwargs['fig'] == []:
            fig = plt.gcf()
        else:
            fig = kwargs['fig']
        if kwargs['ax'] == []:
            ax = plt.gca()
        else:
            ax = kwargs['ax']

        ax.hist(self.H.taud,bins=len(self.H.taud)/2)
        ax.set_xlim([0,max(self.H.taud)])
        return fig,ax

    def plt_aspread(self,**kwargs):
        """ plot angular spread
        """
        defaults = { 'fig':[],
                     'ax':[]
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        if kwargs['fig'] == []:
            fig = plt.gcf()
        else:
            fig = kwargs['fig']
        if kwargs['ax'] == []:
            ax = plt.gca()
        else:
            ax = kwargs['ax']

        ax.hist(self.H.doa[:,0],bins=len(self.H.doa[:,0])/2)
        ax.set_xlim([-np.pi,np.pi])
        return fig,ax

    def _autocufoff(self):
        """ automatically determine minimum cutoff

        See Also
        --------

        pylayers.antprop.loss.losst
        pylayers.gis.layout.angleonlink3

        """

        v = np.vectorize( lambda t:self.L.Gs.node[t]['name'])
        # determine incidence angles on segment crossing p1-p2 segment
        #data = L.angleonlink(p1,p2)
        if np.allclose(self.a,self.b):
            self.cutoff = 2
        else:
            data = self.L.angleonlink3(self.a,self.b)
            # as many slabs as segments and subsegments
            us  = data['s'] 
            if len(us) >0:
                sl = v(us)
                uus = np.where((sl != 'AIR') & (sl != '_AIR'))[0]
                self.cutoff = len(uus)
            else:
                self.cutoff = 2
        return self.cutoff

if (__name__ == "__main__"):
    #plt.ion()
    doctest.testmod()
