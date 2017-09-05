# -*- coding: utf-8 -*-
#
r"""

.. currentmodule:: pylayers.simul.link


=======
This module runs the electromagnetic simulation at the link level.

This module runs the electromagnetic simulation for a link.
A deterministic link has two termination points and an associated Layout
whereas a statistical link do not need any of those precursor object.

It stores simulated objects in `hdf5` format.


Link Class
===========

Link is a MetaClass, which derives from `Tchannel`.
Tchannel is a transmission channel i.e a radio channel
which includes both link termination antennas.

A common factor of both statistical (SLink) and deterministic channel (DLink)
is the exitence of :math:`\alpha_k` and :math:`\tau_k`

.. autosummary::
    :toctree: generated/

    Link.__add__


SLink Class
===========

Slink is for statistical links.

.. autosummary::
    :toctree: generated/

    SLink.onbody

DLink Class
===========

Dlink is for deterministic links

>>> from pylayers.simul.link import *
>>> L = DLink(verbose=False)
>>> aktk = L.eval()

DLink simulation
----------------

.. autosummary::
    :toctree: generated/

    DLink.eval
    DLink._show3


DLink init
----------

.. autosummary::
    :toctree: generated/

    DLink.__init__
    DLink.__repr__
    DLink.reset_config
    DLink.check_grpname


search in hdf5 file
-------------------

.. autosummary::
    :toctree: generated/

    DLink.checkh5
    DLink.array_exist
    DLink.get_grpname
    DLink.get_idx


Modify hdf5 file
----------------

.. autosummary::
    :toctree: generated/


    DLink.save_init
    DLink.save
    DLink.load
    DLink.stack
    DLink._delete




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


import h5py
import pdb



class Link(PyLayers):
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
        >>> L = DLink(verbose=False)
        >>> aktk = L.eval()


        """


        Link.__init__(self)

        defaults={ 'L':'',
                   'a':np.array(()),
                   'b':np.array(()),
                   'Aa':[],
                   'Ab':[],
                   'Ta':np.eye(3),
                   'Tb':np.eye(3),
                   'fGHz':np.array([2.4]),
                   'wav':wvf.Waveform(),
                   'cutoff':3,
                   'threshold':0.8,
                   'save_opt':['sig','ray2','ray','Ct','H'],
                   'save_idx':0,
                   'force_create':False,
                   'verbose':False,
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

        if self.Aa==[]:
            self.Aa=Antenna(typ='Omni',fGHz=self.fGHz)
        if self.Ab==[]:
            self.Ab=Antenna(typ='Omni',fGHz=self.fGHz)
        
        if isinstance(self._L,str):
            self._Lname = self._L
            self._L = Layout(self._Lname,bgraphs=True,bcheck=False)
        else:
            self._Lname = self._L._filename



        if self._Lname != '':

            self.filename = 'Links_' + str(self.save_idx) + '_' + self._Lname + '.h5'
            filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])
            # check if save file alreasdy exists
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
                # lTiw : those which are wall (not those above buildings) 
                lTiw = [ k for k in lTi if self.L.Gs.node[k[0]]['name']!='AIR' ]

                self.L.Gi.remove_nodes_from(lTiw)
                lE = self.L.Gi.edges()
                for k in range(len(lE)):
                    e = lE[k]
                    try:
                        output = self.L.Gi.edge[e[0]][e[1]]['output']
                    except:
                        pdb.set_trace()
                    for l in output.keys():
                        if l in lTiw:
                            del output[l]
                    self.L.Gi.edge[e[0]][e[1]]['output']=output
                
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
            if position[0]<self.L.ax[0]:
                position[0]=self.L.ax[0]
            if position[0]>self.L.ax[1]:
                position[0]=self.L.ax[1]
            if position[1]<self.L.ax[2]:
                position[1]=self.L.ax[2]
            if position[1]>self.L.ax[3]:
                position[1]=self.L.ax[3]
            # raise NameError ('Warning : point a is not inside the Layout')
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
            if position[0]<self.L.ax[0]:
                position[0]=self.L.ax[0]
            if position[0]>self.L.ax[1]:
                position[0]=self.L.ax[1]
            if position[1]<self.L.ax[2]:
                position[1]=self.L.ax[2]
            if position[1]>self.L.ax[3]:
                position[1]=self.L.ax[3]
            # raise NameError ('Warning : point b is not inside the Layout')
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
        """ __repr__
        """

        if hasattr(self,'filename'):
            s = 'filename: ' + self.filename +'\n'

            s = s + 'Link Parameters :\n'
            s = s + '------- --------\n'
            s = s + 'Layout : ' + self.Lname + '\n\n'
            s = s + 'Node a   \n'
            s = s + '------  \n'
            s = s + 'position : ' + str (self.a) + "  in cycle " + str(self.ca) + '\n'
            s = s + 'Antenna : ' + str (self.Aa.typ) + '\n'
            s = s + 'Rotation matrice : \n ' + str (self.Ta) + '\n\n'
            s = s + 'Node b   \n'
            s = s + '------  \n'
            s = s + 'position : ' + str (self.b) + " in cycle " + str(self.cb) + '\n'
            s = s + 'Antenna : ' + str (self.Ab.typ) + '\n'
            s = s + 'Rotation matrice : \n ' + str (self.Tb) + '\n\n'
            s = s + 'Link evaluation information : \n'
            s = s + '----------------------------- \n'
            s = s + 'distance : ' + str("%6.3f" % np.sqrt(np.sum((self.a-self.b)**2))) + ' m \n'
            s = s + 'delay : ' + str("%6.3f" % (np.sqrt(np.sum((self.a-self.b)**2))/0.3)) + ' ns\n'
            rd2deg = 180/np.pi
            if not np.allclose(self.a,self.b):
                vsba = self.b-self.a
                a1 = geu.angledir(vsba[None,:])
                a2 = geu.angledir(-vsba[None,:])
                s = s + 'azimuth (a | b ) : '+str(a1[0,1]*rd2deg)+' deg  | '+str(a2[0,1]*rd2deg)+ ' deg\n'
                s = s + 'elevation (a | b ) : '+str(a1[0,0]*rd2deg)+ ' deg |  '+str(a2[0,0]*rd2deg)+ ' deg\n'
                s = s + 'tilt (a |  b ) : '+str((a1[0,0]-np.pi/2)*rd2deg)+ ' deg  | '+ str((a2[0,0]-np.pi/2)*rd2deg)+ ' deg\n'
            #s = s + 'Frequency range :  \n'
            s = s + 'fmin (fGHz) : ' + str(self.fGHz[0]) +'\n'
            s = s + 'fmax (fGHz) : ' + str(self.fGHz[-1]) +'\n'
            Nf = len(self.fGHz)
            if Nf>1:
                s = s + 'fstep (fGHz) : ' + str(self.fGHz[1]-self.fGHz[0]) +'\n'
            else:
                s = s + 'fstep (fGHz) : ' + str(self.fGHz[0]-self.fGHz[0]) +'\n'
            s = s + 'Nf : ' + str(Nf) +'\n '
            d =  np.sqrt(np.sum((self.a-self.b)**2))
            if Nf>1:
                fcGHz = (self.fGHz[-1]+self.fGHz[0])/2.
            else:
                fcGHz = self.fGHz[0]
            L  = 32.4+20*np.log(d)+20*np.log10(fcGHz)
            s = s + 'Algorithm information : \n'
            s = s + '----------------------------- \n'
            s = s + 'cutoff : '+ str(self.cutoff)+'\n'
            s = s + 'threshold :'+ str(self.threshold)+'\n'
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
        print "Ray : "+str(iray)
        if not self.R.evaluated:
            self.R.eval()
        PM = self.R.info(iray,ifGHz=0,matrix=1)
        print "Propagation Channel 2x2 (C):"
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
            if key != 'T_map':
                sc = f[key].shape
                f[key].resize((sc[0]+1,sc[1]))
                f[key][-1,:]=array
            else:
                sc = f[key].shape
                f[key].resize((sc[0]+1,sc[1],sc[2]))
                f[key][-1,:,:]=array
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
              
            obj._saveh5(self.filename,grpname)

        if self.verbose :
            print(str(obj.__class__).split('.')[-1] + ' from '+ grpname + ' saved')


    def load(self,obj,grpname,**kwargs):
        """ Load a given object in the correct grp

        Parameters
        ----------

        obj : Object
            (Signatures|Rays|Ctilde|Tchannel)
        grpname : string
            groupe name of the h5py file

        kwargs : 
        layout for sig and rays

        Examples
        --------


        """



        obj._loadh5(self.filename,grpname,**kwargs)
        if self.verbose :
            print(str(obj.__class__).split('.')[-1] + ' from '+ grpname + ' loaded')


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
        self.dexist['ray2']['grpname']=grpname
        self.dexist['ray']['grpname']=grpname
        


        ############
        # Ctilde
        #############

        # check existence of frequency in h5py file
        #farray = np.array(([self.fmin,self.fmax,self.fstep]))
        Nf = len(self.fGHz)
        if Nf > 1:
            farray = np.array(([self.fGHz[0],self.fGHz[-1],self.fGHz[1]-self.fGHz[0]]))
        else:
            farray = np.array(([self.fGHz[0],self.fGHz[-1],0]))
        uf_opt, uf = self.get_idx('f_map',farray)

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
            if grpname.decode('utf8') in f[key].keys():
                self.dexist[key]['exist']=True
            else :
                self.dexist[key]['exist']=False
            f.close()
        except:
            f.close()
            raise NameError('Link exist: issue during stacking')


    def get_idx(self,key,array,tol=1e-3):
        """ try to get the index of the requested array in the group key
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
            # old version
                # ufmi = np.where(fa[:,0]<=array[0])[0]
                # lufmi = len(ufmi)

            #### fmax_h5 > fmax_rqst
            ufma = fa[:,1]>=array[1]
            # old version
                # ufma = np.where(fa[:,1]>=array[1])[0]
                # lufma = len(ufma)

            ### fstep_h5 < fstep_rqst
            ufst = fa[:,2]<=array[2]
            # old version
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
            ua = np.where(fa==array)[0]

        elif key == 'T_map':
            eq = array == fa
            seq = np.sum(np.sum(eq,axis=1),axis=1)
            ua = np.where(seq==9)[0]
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
        H = self.C.prop2tran(a=self.Aa,b=self.Ab,Friis=True,debug=True)
        self.H = H

    def eval(self,**kwargs):
        """ evaluate the link


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


    
        Notes    def eval(self,**kwargs):
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
            >>> L=DLink(verbose=False)
            >>> aktk = L.eval()


        See Also
        --------

        pylayers.antprop.signature
        pylayers.antprop.rays


        """



        defaults={ 'applywav':False,
                   'si_progress':True,
                   'diffraction':True,
                   'ra_vectorized':True,
                   'ra_ceil_H':[],
                   'ra_number_mirror_cf':1,
                   'force':[],
                   'bt':True,
                   'alg':1,
                   'si_reverb':4,
                   'nD':2,
                   'nR':10,
                   'nT':10,
                   'debug':False,
                   'verbose':[],
                   'progressbar':None,
                   }
        # check antenna frequency range compatibility
        if (self.Aa.fGHz!=self.Ab.fGHz).all():
            raise AttributeError("Antenna frequency range are not compatible")

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

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
                    # Ct and H are not yet saved/load 
                    # compliantly with the given configutain
                    # their are disabled here
                    kwargs['force'] = ['Ct','H']

        if kwargs['verbose'] != []:
            self.verbose=kwargs['verbose']

        
        # must be placed after all the init !!!!
        # if self.verbose :
        #     print("checkh5")
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
        
        if self.verbose :
            print("Start Signatures")
        tic = time.time()
        Si = Signatures(self.L,
                        self.ca,
                        self.cb,
                        cutoff=kwargs['cutoff'],
                        threshold = kwargs['threshold'])

        if (self.dexist['sig']['exist'] and not ('sig' in kwargs['force'])):
            self.load(Si,self.dexist['sig']['grpname'],L=self.L)
            if self.verbose :
                print("load signature")
        else :
            ## 1 is the default signature determination algorithm
            if kwargs['alg']==1:
                Si.run(cutoff = kwargs['cutoff'],
                        diffraction = kwargs['diffraction'],
                        threshold = kwargs['threshold'],
                        nD = kwargs['nD'],
                        nR = kwargs['nR'],
                        nT = kwargs['nT'],
                        progress = kwargs['si_progress'],
                        bt = kwargs['bt'])

                if self.verbose:
                    print("default algorithm")

            if kwargs['alg']=='exp':
                TMP = Si.run_exp(cutoff=kwargs['cutoff'],
                         cutoffbound=kwargs['si_reverb'])
                if self.verbose :
                    print("experimental (ex 2015)")

            if kwargs['alg']=='exp2':
                TMP = Si.run_exp2(cutoff=kwargs['cutoff'],
                        cutoffbound=kwargs['si_reverb'])
                if self.verbose :
                    print("algo exp2 ( ex 20152)")

        #Si.run6(diffraction=kwargs['diffraction'])
        # save sig
            
            self.save(Si,'sig',self.dexist['sig']['grpname'],force = kwargs['force'])

        self.Si = Si
        toc = time.time()
        if self.verbose :
            print("Stop signature",toc-tic)
        try:
            pbar.update(20)
        except: 
            pass



        ############
        # Rays
        ############

        if self.verbose :
            print("Start Rays")
        tic = time.time()
        r2d = Rays(self.a,self.b)
       

        #
        # get 2D rays 
        #
        if self.dexist['ray2']['exist'] and not ('ray2' in kwargs['force']):
            self.load(r2d,self.dexist['ray2']['grpname'],L=self.L)
        else :
            # perform computation ...
            # ... with vectorized ray evaluation 
            if kwargs['ra_vectorized']:
                r2d = Si.raysv(self.a,self.b)
            # ... or with original and slow approach ( to be removed in a near future)
            else :
                r2d = Si.rays(self.a,self.b)
            # save 2D rays
            self.save(r2d,'ray2',self.dexist['ray2']['grpname'],force = kwargs['force'])

        self.r2d = r2d

        #
        # get 3D rays 
        #
        R = Rays(self.a,self.b)
        R.is3D = True
        if self.dexist['ray']['exist'] and not ('ray' in kwargs['force']):
            self.load(R,self.dexist['ray']['grpname'],L=self.L)
        else :

            if kwargs['ra_ceil_H'] == []:
                if self.L.typ=='indoor':
                    ceilheight = self.L.maxheight
                else:
                    ceilheight = 0 
            else:
                ceilheight = kwargs['ra_ceil_H']


            R = self.r2d.to3D(self.L,H=ceilheight, N=kwargs['ra_number_mirror_cf'])

            R.locbas(self.L)
            

            R.fillinter(self.L)

            # C = Ctilde()

            # C = R.eval(self.fGHz)

            # save 3D rays 

            self.save(R,'ray',self.dexist['ray']['grpname'],force = kwargs['force'])

        
        self.R = R
        toc = time.time()
        if self.verbose :
            print("Stop rays",toc-tic)
        
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


    def afp(self,fGHz,vl,pl,az=0,tilt=0,polar='V'):
        """ Evaluate angular frequency profile 

        Parameters
        ----------

        fGHz  : np.array 
            frequency range 
        vl : np.array (,3)
            main beam direction in local frame
        pl : np.array (,3)
            polarisation vector in a plane perpendicular to vl 
        az : azimuth angle (radian)  
        tilt : tilt angle (-pi/2<tilt<pi/2) 
        polar : string

        """

        # create an empty AFP 
        # tx = a
        # rx = b 
        # angular range (a) : phi 
        #
        afp = AFPchannel(tx=self.a,rx=self.b,az=az)
        for ph in az:
            self.Tb = geu.MATP(vl,pl,ph,tilt,polar)
            # self._update_show3(ant='b')
            # pdb.set_trace()
            self.evalH()


            if self.H.y.shape[3]!=1:
                S = np.sum(self.H.y*np.exp(-2*1j*np.pi*self.H.x[None,None,None,:]*self.H.taud[:,None,None,None]),axis=0)
            else:
                S = np.sum(self.H.y*np.exp(-2*1j*np.pi*fGHz*self.H.taud[:,None,None,None]),axis=0)

            try:
                afp.y = np.vstack((afp.y,np.squeeze(S)))
            except:
                afp.y = np.squeeze(S)

        if self.H.y.shape[3]!=1:
            afp.x = self.H.x
        else:
            afp.x = fGHz



        return(afp)



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
        ca  : string
            color a
        cb  : string
            color b
        alpha : int
        axis : boolean 
            display axis boolean (default True)
        figsize : tuple
            (20,10)
        fontsize : int
            20
        rays : boolean
            False
        bsig : boolean 
            False    
        laddr : list 
            list of signature addresses 
        cmap : colormap
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
            dynamic in dB

        Examples
        --------

        >>> from pylayers.simul.link import *
        >>> DL=Link()
        >>> DL.show(lr=-1,rays=True,dB=True,col='cmap',cmap=plt.cm.jet)
        >>> DL.show(laddr=[(6,2)],bsig=True)

        """
        defaults ={'s':80,   # size points
                   'ca':'b', # color a 
                   'cb':'r', # color b 
                   'alpha':1,
                   'axis':True,
                   'lr':-1,
                   'ls':-1,
                   'figsize':(20,10),
                   'fontsize':20,
                   'rays':False,
                   'bsig':True,
                   'laddr':[(1,0)],
                   'cmap':plt.cm.hot,
                   'pol':'tot',
                   'col':'k',
                   'width':1,
                   'alpha':1,
                   'col':'k',
                   'dB':False,
                   'labels':False,
                   'aw':False,
                   'dyn':70}

        for key in defaults:
            if key not in kwargs:
                kwargs[key]=defaults[key]

        #
        # Layout
        #
        fig,ax = self.L.showG('s',nodes=False,figsize=kwargs['figsize'],labels=kwargs['labels'],aw=kwargs['aw'],axis=kwargs['axis'])
        #
        # Point A
        #
        self.caf = ax.scatter(self.a[0],
                   self.a[1], c=kwargs['ca'], s=kwargs['s'],
                   alpha=kwargs['alpha'])
        ax.text(self.a[0]+0.1,self.a[1]+0.1,'A',fontsize=kwargs['fontsize'])
        #
        # Point B
        #
        self.cbf = ax.scatter(self.b[0],
                   self.b[1], c=kwargs['cb'], s=kwargs['s'],
                   alpha=kwargs['alpha'])
        ax.text(self.b[0]-0.1,self.b[1]+0.1,'B',fontsize=kwargs['fontsize'])
        #
        # Plot Rays
        #
        if kwargs['rays']:
            ECtt,ECpp,ECtp,ECpt = self.C.energy()
            if kwargs['pol']=='tt':
                val = ECtt
            if kwargs['pol']=='pp':
                val = ECpp
            if kwargs['pol']=='tp':
                val = ECtp
            if kwargs['pol']=='pt':
                val = ECpt
            if kwargs['pol']=='tot':
                val = ECtt+ECpp+ECpt+ECtp
            if kwargs['pol']=='co':
                val = ECtt+ECpp
            if kwargs['pol']=='cross':
                val = ECtp+ECpt

            clm = kwargs['cmap']
            #
            # Select group of interactions
            #
            if kwargs['lr']==-1:
                lr  = np.arange(self.R.nray)
            else:
                lr = kwargs['lr']
            
            vmin = val.min()
            vmax = val.max() 
            if kwargs['dB']:
                vmin = 20*np.log10(vmin)
                vmax = 20*np.log10(vmax)

            for ir  in lr:
                if kwargs['dB']:
                    RayEnergy=max((20*np.log10(val[ir]/val.max())+kwargs['dyn']),0)/kwargs['dyn']
                else:
                    RayEnergy=val[ir]/val.max()

                if kwargs['col']=='cmap':
                    col = clm(RayEnergy)
                    width = 3*RayEnergy
                    alpha = 1
                else:
                    col = kwargs['col']
                    width = kwargs['width']
                    alpha = kwargs['alpha']
                
                # plot ray (i,r) 
                fig,ax = self.R.show(rlist=[ir],
                               colray=col,
                               widthray=width,
                               alpharay=alpha,
                               fig=fig,ax=ax,
                               layout=False,
                               points=False)
            if kwargs['col']=='cmap':
                sm = plt.cm.ScalarMappable(cmap=kwargs['cmap'], norm=plt.Normalize(vmin=vmin, vmax=vmax))
                sm._A = []
                plt.colorbar(sm)
        #
        # Plot Rays
        #
        if kwargs['bsig']:
            for addr in kwargs['laddr']: 
                seq = self.Si[addr[0]][2*addr[1]:2*addr[1]+2,:]
                Si = Signature(seq)
                fig,ax = Si.show(self.L,self.a[0:2],self.b[0:2],fig=fig,ax=ax)

        return fig,ax

    def _show3(self,rays=True, lay= True, ant= True, newfig= False, **kwargs):
        """ display the simulation scene using Mayavi
            using Mayavi


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
            >>> L=DLink(verbose=False)
            >>> L.eval()

        """

        if not newfig:
            self._maya_fig=mlab.gcf()
        else:
            self._maya_fig=mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))

        if 'centered' in kwargs:
            centered = kwargs['centered']
        else :
            centered = False

        if centered:
            pg =np.zeros((3))
            pg[:2]=self.L.pg


        if centered :
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
            try:
                Atx._show3(T=Ttx.reshape(3,3),
                           po=ptx,
                           title=False,
                           bcolorbar=False,
                           bnewfig=False,
                           bcircle = False,
                           name = Atx._filename,
                           scale= 0.5,
                           binteract=False)
            except:
                Atx.eval()
                Atx._show3(T=Ttx.reshape(3,3),
                            po=ptx,
                            title=False,
                            bcolorbar=False,
                            bnewfig=False,
                            bcircle = False,
                            name = Atx._filename,
                            scale= 0.5,
                            binteract=False)
            if not Arx.evaluated:
                Arx.eval()
            try:
                Arx._show3(T=Trx.reshape(3,3),
                            po=prx,
                            title=False,
                            bcolorbar=False,
                            bnewfig=False,
                            bcircle = False,
                            name = Arx._filename,
                            scale= 0.5,
                            binteract=False)
            except:
                Arx.eval()
                Arx._show3(T=Trx.reshape(3,3),
                            po=prx,
                            title=False,
                            bcolorbar=False,
                            bnewfig=False,
                            bcircle = False,
                            name = Arx._filename,
                            scale= 0.5,
                            binteract=False)
        if lay:
            # check if indoor/outdoor, outdoor or indoor situations
            # a_in = self.L.Gt.node[self.ca]['indoor']
            # b_in = self.L.Gt.node[self.cb]['indoor']

            # if (a_in) & (b_in):
            #     # indoor
            #     show_ceil=False
            #     opacity = 0.7
            #     ceil_opacity = 0.
            # elif ((not a_in) & (not b_in)):
            #     # outdoor
            #     show_ceil=True
            #     opacity = 1.
            #     ceil_opacity = 1.
            # else:
            #     # indoor/outdoor
            #     show_ceil=True
            #     opacity = 0.7
            #     ceil_opacity = 0.7

            if self.L.typ == 'outdoor':
                show_ceil=True
                opacity = 1.
                ceil_opacity = 1.
            elif self.L.typ == 'indoor':
                show_ceil=False
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
        return self._maya_fig
        #return(self._maya_fig)


    def _update_show3(self,ant='a',delrays=False):
        """
        """


        view=mlab.view()


        antenna = eval('self.A'+ant)
        rot = eval('self.T'+ant).reshape(3,3)
        pos = eval('self.'+ant)

        if not antenna.full_evaluated:
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

        # # update wall opaccity
        
        
        # ds  =[i for i in self._maya_fig.children if self.L._filename in i.name][0]
        # a_in = self.L.Gt.node[self.ca]['indoor']
        # b_in = self.L.Gt.node[self.cb]['indoor']

        # if 
        # if a_in or b_in:
        #     # indoor situation
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
        
        See Also
        --------

        pylayers.antprop.channel.Tchannel.getcir

        """

        defaults = {'fig':[],
                    'ax': [],
                     'BWGHz':5,
                    'Nf':1000,
                    'rays':True,
                    'fspl':True,
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

        # getcir is a Tchannel method 
        ir = self.H.getcir(BWGHz = kwargs['BWGHz'],Nf=kwargs['Nf'])
        ir.plot(fig=fig,ax=ax)
        delay = ir.x 
        dist = delay*0.3
        #FSPL = -32.4- 20*np.log10(self.fGHz[0])-20*np.log10(dist) + 19 + 2
        #if kwargs['fspl']:
        #    ax.plot(delay,FSPL,linewidth=2,color='b')


        if kwargs['rays'] : 
            ER = np.squeeze(self.H.energy())
            color_range = np.linspace( 0, 1., len(ER))#np.linspace( 0, np.pi, len(ER))
            uER = ER.argsort()[::-1]
            colors= color_range[uER]
            ax.scatter(self.H.taud[uER],20*np.log10(self.H.y[uER,0,0,0]),c=colors,cmap='hot')
            ax.set_xlim([min(self.H.taud)-10,max(self.H.taud)+10])


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
            us    = data['s'] 
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
