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
    print 'Layout:Mayavi is not installed'
import doctest
import time
import numpy as np
import matplotlib.pylab as plt
import pylayers.signal.waveform as wvf

from pylayers.util.project import *
import pylayers.util.pyutil as pyu
from pylayers.simul.radionode import RadioNode
# Handle Layout
from pylayers.gis.layout import Layout
# Handle Antenna
from pylayers.antprop.antenna import Antenna

# Handle Signauture
from pylayers.antprop.signature import Signatures
# Handle Rays
from pylayers.antprop.rays import Rays
# Handle VectChannel and ScalChannel
from pylayers.antprop.channel import Ctilde, Tchannel
from pylayers.antprop.statModel import getchannel
import tqdm


import h5py
import pdb



class Link(Tchannel):
    def __init__(self):
        """ Link evaluation metaclass
        """
        Tchannel.__init__(self)
    #    super(Link,self).__init__ ()


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
            forcecreating the h5py file (if already exist, will be erased)


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
                ca_cb_cutoff

            Ray identifier (ray_ID#N):
                cutoff_ua_ub

            Ctilde identifier (Ct_ID#N):
                ua_ub_uf

            H identifier (H_ID#N):
                ua_ub_uf_uTa_uTb_uAa_uAb

            with
            ca : cycle number of a
            cb : cycle number of b
            cutoff : signature.run cutoff
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

        defaults={ 'L':Layout('defstr.ini'),
                   'a':np.array(()),
                   'b':np.array(()),
                   'Aa':[],
                   'Ab':[],
                   'Ta':np.eye(3),
                   'Tb':np.eye(3),
                   'fGHz':[],
                   'wav':wvf.Waveform(),
                   'cutoff':3,
                   'save_opt':['sig','ray','Ct','H'],
                   'save_idx':0,
                   'force_create':False,
                   'verbose':False,
                   'graph':'tcvirw'
                }

        # self._ca = -1
        # self._cb = -1

        specset  = ['a','b','Aa','Ab','Ta','Tb','L','fGHz','wav']

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

    
        force=self.force_create
        delattr(self,'force_create')

       
        # The link frequency range depends on the antenna frequency range

        if self.fGHz == []:
            self.initfreq()
        else :
            pass

        
        if self.Aa==[]:
            self.Aa=Antenna(typ='Omni',fGHz=self.fGHz)
        if self.Ab==[]:
            self.Ab=Antenna(typ='Omni',fGHz=self.fGHz)
        


       
        self._Lname = self._L._filename
    

        self.filename = 'Links_' + str(self.save_idx) + '_' + self._Lname + '.h5'
        filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])
        # check if save file alreasdy exists
        if not os.path.exists(filenameh5) or force:
            print 'Links save file for ' + self.L._filename + ' does not exist.'
            print 'Creating file. You\'ll see this message only once per Layout'
            self.save_init(filenameh5)

        # dictionnary data exists
        self.dexist={'sig':{'exist':False,'grpname':''},
                     'ray':{'exist':False,'grpname':''},
                     'Ct':{'exist':False,'grpname':''},
                     'H':{'exist':False,'grpname':''}
                    }
       
        try:
            self.L.dumpr()
        except:
            print('This is the first time the Layout is used. Graphs have to be built. Please Wait')
            self.L.build(graph=self.graph)
            self.L.dumpw()
        #self.L.build()

        ###########
        # init pos & cycles
        #
        # If a and b are not specified
        #  they are chosen as center of gravity of cycle 0
        #
        ###########
        nodes = self.L.Gt.nodes()
        nodes = [n for n in nodes if n!=0 and self.L.Gt.node[n]['indoor']]
        if len(self.a)==0:
            self.ca = nodes[0]
            # self.a = self.L.cy2pt(self.ca)
        else:
            if len(kwargs['a']) ==2:
                a=np.r_[kwargs['a'],1.0]
            else:
                a=kwargs['a']
            self.ca = self.L.pt2cy(a)
            self.a = a

        if len(self.b)==0:
            if len(self.L.Gt.node)>2:
                self.cb = nodes[1]
            else:
                self.cb = nodes[0]
            self.b = self.L.cy2pt(self.cb)
        else:
            if len(kwargs['b']) ==2:
                b=np.r_[kwargs['b'],1.0]
            else:
                b=kwargs['b']
            self.cb = self.L.pt2cy(b)
            self.b = b


       
        ###########
        # init freq
        # TODO Check where it is used redocdundant with fGHz
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
            self._L = Layout(L)
            self._Lname = L
        elif isinstance(L,Layout):
            self._L = L
            self._Lname = L.filename

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
        

    @ca.setter
    def ca(self,cycle):
        if not cycle in self.L.Gt.nodes():
            raise NameError ('cycle ca is not inside Gt')

        self._ca = cycle
        self.a = self.L.cy2pt(cycle)

    @cb.setter
    def cb(self,cycle):
        if not cycle in self.L.Gt.nodes():
            raise NameError ('cycle cb is not inside Gt')
        self._cb = cycle
        self.b = self.L.cy2pt(cycle)

    @Aa.setter
    def Aa(self,Ant):

        if hasattr(self.Aa,'_mayamesh'):
            self.Aa._mayamesh.remove()

        # save rot
        rot = self.Ta
        self._Aa = Ant
        self.Ta = rot
        self.initfreq()

        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='a')


    @Ab.setter
    def Ab(self,Ant):
        if hasattr(self.Ab,'_mayamesh'):
           self.Ab._mayamesh.remove()

        #save rot
        rot = self.Tb
        self._Ab = Ant
        self.Tb = rot
        self.initfreq()

        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='b')


    @Ta.setter
    def Ta(self,orientation):
        self._Ta = orientation
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='a')

    @Tb.setter
    def Tb(self,orientation):
        self._Tb = orientation
        if hasattr(self,'_maya_fig') and self._maya_fig._is_running:
            self._update_show3(ant='b')

    @fGHz.setter
    def fGHz(self,freq):
        if not isinstance(freq,np.ndarray):
            freq=np.array([freq])
        self._fGHz = freq
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

        s = 'filename: ' + self.filename +'\n'

        s = s + 'Link Parameters :\n'
        s = s + '------- --------\n'
        s = s + 'Layout : ' + self.Lname + '\n\n'
        s = s + 'Node a   \n'
        s = s + '------  \n'
        s = s + 'position : ' + str (self.a) + '\n'
        s = s + 'Antenna : ' + str (self.Aa.typ) + '\n'
        s = s + 'Rotation matrice : \n ' + str (self.Ta) + '\n\n'
        s = s + 'Node b   \n'
        s = s + '------  \n'
        s = s + 'position : ' + str (self.b) + '\n'
        s = s + 'Antenna : ' + str (self.Ab.typ) + '\n'
        s = s + 'Rotation matrice : \n ' + str (self.Tb) + '\n\n'
        s = s + 'Link evaluation information : \n'
        s = s + '----------------------------- \n'
        s = s + 'distance : ' + str("%6.3f" % np.sqrt(np.sum((self.a-self.b)**2))) + ' m \n'
        s = s + 'delay : ' + str("%6.3f" % (np.sqrt(np.sum((self.a-self.b)**2))/0.3)) + ' ns\n'
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
        return s


    def initfreq(self):
        """ Automatic freq determination from
            Antennas
        """
        #sf = self.fGHz[1]-self.fGHz[0]
        sf = 1e15
        if hasattr(self.Aa,'fGHz'):
            fa = self.Aa.fGHz
            if len(fa)==0:
                fa = np.array([2.4])
                self.Aa.fGHz = fa
                # raise AttributeError("Incompatible frequency range in Antenna. Consider change Dlink.fGHz") 
                print "Incompatible frequency range in Antenna. WARNING  Dlink.fGHz changed to 2.4GHz"
            try:
                sa = fa[1]-fa[0]  # step frequecy 
            except: #single frequency
                sa = fa[0]
            # step
            if len(self.fGHz)>0:
                minfa = max(min(fa),min(self.fGHz))
                maxfa = min(max(fa),max(self.fGHz))
            else:
                minfa = min(fa)
                maxfa = max(fa)
            sf = min(sa,sf)
            self.fGHz = np.arange(minfa,maxfa+sf,sf)

        elif hasattr(self.Ab,'fGHz'):
            fb = self.Ab.fGHz
            if len(fb)==0:
                # raise AttributeError("Incompatible frequency range in Antenna. Consider change Dlink.fGHz")
                fb = np.array([2.4])
                self.Ab.fGHz=fb
                # raise AttributeError("Incompatible frequency range in Antenna. Consider change Dlink.fGHz") 
                print "Incompatible frequency range in Antenna. WARNING  Dlink.fGHz changed to 2.4GHz"
        
            try:
                sb = fb[1]-fb[0] # step frequency 
            except:
                sb = fb[0]

            if len(self.fGHz)>0:
                minfb = max(min(self.fGHz),min(fb))
                maxfb = min(max(self.fGHz),max(fb))
            else:
                minfb = min(fb)
                maxfb = max(fb)

            sf = min(sf,sb)
            self.fGHz = np.arange(minfb,maxfb+sf,sf)
        else:
            self.fGHz = np.array([2.3,2.4,2.5])


    def reset_config(self):
        """ reset configuration when a new layout is loaded
        """
        try:
            self.L.dumpr()
        except:
            self.L.build()
            self.L.dumpw()


        self.ca = 1
        self.cb = 1
        # self.a = self.L.cy2pt(self.ca)
        # self.b = self.L.cy2pt(self.cb)

        # change h5py file if layout changed
        self.filename = 'Links_' + str(self.save_idx) + '_' + self._Lname + '.h5'
        filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])
        if not os.path.exists(filenameh5) :
            print 'Links save file for ' + self.L.filename + ' does not exist.'
            print 'It is beeing created. You\'ll see that message only once per Layout'
            self.save_init(filenameh5)


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
        """ initialize save Link

        Parameters
        ----------

        filename_long : str
            complete path and filename


        """


        f=h5py.File(filename_long,'w')
        # try/except to avoid loosing the h5 file if
        # read/write error

        try:

            f.create_group('sig')
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
        lfilename=pyu.getlong(self.filename,pstruc['DIRLNK'])
        f=h5py.File(lfilename,'a')
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
            print str(obj.__class__).split('.')[-1] + ' from '+ grpname + ' saved'


    def load(self,obj,grpname):
        """ Load a given object in the correct grp

        Parameters
        ----------

        obj : Object
            (Signatures|Rays|Ctilde|Tchannel)
        grpname : string
            groupe name of the h5py file

        Examples
        --------


        """

        obj._loadh5(self.filename,grpname)
        if self.verbose :
            print str(obj.__class__).split('.')[-1] + ' from '+ grpname + ' loaded'


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

        grpname = str(self.ca) + '_' +str(self.cb) + '_' + str(self.cutoff)
        self.dexist['sig']['grpname']=grpname



        ############
        # Rays
        #############

        # check existence of self.a in h5py file

        ua_opt, ua = self.get_idx('p_map',self.a)
        # check existence of self.b in h5py file
        ub_opt, ub = self.get_idx('p_map',self.b)
        # Write in h5py if no prior a-b link

        grpname = str(self.cutoff) + '_' + str(ua) + '_' +str(ub)
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
        """ check if an array of a given key (h5py group)
            has already been stored into the h5py file


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
            # fmin_h5 < fmin_rqst
            ufmi = np.where(fa[:,0]<=array[0])[0]
            lufmi = len(ufmi)
            # fmax_h5 > fmax_rqst
            ufma = np.where(fa[:,1]>=array[1])[0]
            lufma = len(ufma)
            # fstep_h5 < fstep_rqst
            ufst = np.where(fa[:,2]<=array[2])[0]
            lufst = len(ufst)
            # if fmin, fmax or fstep
            if (lufmi==0) and (lufma==0) and (lufst==0):
                ua = np.array([])
            else :
                # find comon lines of fmin and fmax
                ufmima = np.where(np.in1d(ufmi,ufma))[0]
                # find comon lines of fmin, fmax and fstep
                ua = np.where(np.in1d(ufmima,ufst))[0]

        elif key == 'A_map':
            ua = np.where(fa==array)[0]

        elif key == 'T_map':
            eq = array == fa
            seq = np.sum(np.sum(eq,axis=1),axis=1)
            ua = np.where(seq==9)[0]
        else :
            raise NameError('Link.array_exist : invalid key')

        return ua


    def eval(self,**kwargs):
        """ evaluate the link


        Parameters
        ----------

        applywav :boolean
         Apply waveform to H
        force : list
            Force the computation (['sig','ray','Ct','H']) AND save (replace previous computations)
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


        Returns
        -------

        ak : ndarray
            alpha_k
        tk : ndarray
            tau_k

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
            >>> L = DLink(verbose=False)
            >>> aktk = L.eval()


        See Also
        --------

        pylayers.antprop.signature
        pylayers.antprop.rays

        Experimental
        ------------

        alg = 2015 | 20152 (best)
            vectorized signature research
        si_reverb : number of reverb in source/target cycle if alg=2015

        """



        defaults={ 'applywav':True,
                   'si_progress':False,
                   'diffraction':True,
                   'ra_vectorized':True,
                   'ra_ceil_H':[],
                   'ra_number_mirror_cf':1,
                   'force':[],
                   'alg':1,
                   'si_reverb':4,
                   'threshold':0.1,
                   'verbose':[],
                   'progressbar':None,
                   }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key]=value

        if 'cutoff' not in kwargs:
            kwargs['cutoff']=self.cutoff
        else:
            self.cutoff=kwargs['cutoff']

        if 'force' in kwargs:
            if not isinstance(kwargs['force'],list):
                if kwargs['force'] == True :
                    kwargs['force'] = ['sig','ray','Ct','H']
                else :
                    kwargs['force'] = []

        if kwargs['verbose'] != []:
            self.verbose=kwargs['verbose']


        #pdb.set_trace()
        # must be placed after all the init !!!!
        if self.verbose :
            print "checkh5"
        self.checkh5()

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
            print "Start Signatures"
        tic = time.time()
        Si = Signatures(self.L,self.ca,self.cb,cutoff=kwargs['cutoff'])

        if (self.dexist['sig']['exist'] and not ('sig' in kwargs['force'])):
            self.load(Si,self.dexist['sig']['grpname'])
            if self.verbose :
                print "load signature"
        else :
            if kwargs['alg']==1:
                Si.run(cutoff=kwargs['cutoff'],
                        diffraction=kwargs['diffraction'],
                        threshold=kwargs['threshold'],
                        progress=kwargs['si_progress'])
                if self.verbose :
                    print "default algorithm"

            if kwargs['alg']=='exp':
                TMP=Si.run_exp(cutoff=kwargs['cutoff'],
                        cutoffbound=kwargs['si_reverb'])
                if self.verbose :
                    print "experimental (ex 2015)"

            if kwargs['alg']=='exp2':
                TMP=Si.run_exp2(cutoff=kwargs['cutoff'],
                        cutoffbound=kwargs['si_reverb'])
                if self.verbose :
                    print "algo exp2 ( ex 20152)"

        #Si.run6(diffraction=kwargs['diffraction'])
        # save sig
            
            self.save(Si,'sig',self.dexist['sig']['grpname'],force = kwargs['force'])

        self.Si = Si
        toc = time.time()
        if self.verbose :
            print "Stop signature",toc-tic
        try:
            pbar.update(20)
        except: 
            pass



        ############
        # Rays
        ############

        if self.verbose :
            print "Start Rays"
        tic = time.time()
        R = Rays(self.a,self.b)

        if self.dexist['ray']['exist'] and not ('ray' in kwargs['force']):
            self.load(R,self.dexist['ray']['grpname'])

        else :

            # perform computation ...
            # ... with vetorized ray evaluation approach
            if kwargs['ra_vectorized']:
                r2d = Si.raysv(self.a,self.b)
            # ... or with original and slow approach ( to be removed in a near future)
            else :
                r2d = Si.rays(self.a,self.b)

            if kwargs['ra_ceil_H'] == []:
                ceilheight = self.L.maxheight
            else:
                ceilheight = kwargs['ra_ceil_H']

            R = r2d.to3D(self.L,H=ceilheight, N=kwargs['ra_number_mirror_cf'])

            R.locbas(self.L)
            # ...and save

            R.fillinter(self.L)

            C = Ctilde()
            C = R.eval(self.fGHz)
            self.save(R,'ray',self.dexist['ray']['grpname'],force = kwargs['force'])

        self.R = R
        toc = time.time()
        if self.verbose :
            print "Stop rays",toc-tic
        
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
            C=Ctilde()
            self.load(C,self.dexist['Ct']['grpname'])

        else :
            #if not hasattr(R,'I'):
            # Ctilde...
            # Find an other criteria in order to decide whether the R has
            # already been evaluated
            #pdb.set_trace()
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
            Cl=C.locbas(Tt=self.Ta, Tr=self.Tb)
            #T channel
            H = C.prop2tran(a=self.Aa,b=self.Ab,Friis=True,debug=True)
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
        return self.H.ak, self.H.tk

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
        figsize : tuple
            (20,10)
        fontsize : int
            20
        rays : boolean
            False
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
        >>> L=Link()
        >>> L.show(rays=True,dB=True)

        """

        defaults ={'s':80,
                   'ca':'b',
                   'cb':'r',
                   'alpha':1,
                   'i':-1,
                   'figsize':(20,10),
                   'fontsize':20,
                   'rays':False,
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
        fig,ax = self.L.showG('s',nodes=False,figsize=kwargs['figsize'],labels=kwargs['labels'],aw=kwargs['aw'])
        plt.axis('off')
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
        # Rays
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
            if kwargs['i']==-1:
                li  = self.R.keys()
            else:
                li = kwargs['i']

            for i  in li:
                lr = self.R[i]['rayidx']
                for r in range(len(lr)):
                    ir = lr[r]
                    try:
                        if kwargs['dB']:
                            RayEnergy=max((20*np.log10(val[ir]/val.max())+kwargs['dyn']),0)/kwargs['dyn']
                        else:
                            RayEnergy=val[ir]/val.max()
                    except:
                        pass
                    if kwargs['col']=='cmap':
                        col = clm(RayEnergy)
                        width = RayEnergy
                        alpha = 1
                        #alpha = RayEnergy
                    else:
                        col = kwargs['col']

                        width = kwargs['width']
                        alpha = kwargs['alpha']

                    fig,ax = self.R.show(i=i,r=r,
                                   colray=col,
                                   widthray=width,
                                   alpharay=alpha,
                                   fig=fig,ax=ax,
                                   layout=False,
                                   points=False)

        return fig,ax

    def _show3(self,rays=True, lay= True, ant= True, newfig= False, **kwargs):
        """ display the simulation scene using Mayavi
            using Mayavi


        Parameters
        ----------

        rays: Ray3d object :
            display the rays of the simulation
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
            >>> aktk = L.eval()

        """

        if not newfig:
            self._maya_fig=mlab.gcf()

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
                Atx._show3(T=Ttx.reshape(3,3),po=ptx,
                title=False,colorbar=False,newfig=False,interact=False)
            except:
                Atx.eval()
                Atx._show3(T=Ttx.reshape(3,3),po=ptx,
                title=False,colorbar=False,newfig=False,interact=False)
            if not Arx.evaluated:
                Arx.eval()
            try:
                Arx._show3(T=Trx.reshape(3,3),po=prx,
                title=False,colorbar=False,newfig=False,name = '',interact=False)
            except:
                Arx.eval()
                Arx._show3(T=Trx.reshape(3,3),po=prx,
                title=False,colorbar=False,newfig=False,name = '',interact=False)

        if lay:
            self.L._show3(newfig=False,opacity=0.7,centered=centered,**kwargs)

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
            try:
                self.R._show3(**kwargs)
            except:
                print 'Rays not computed yet'

        self._maya_fig.scene.disable_render = False


    def _update_show3(self,ant='a',delrays=False):
        """
        """

        antenna = eval('self.A'+ant)
        rot = eval('self.T'+ant).reshape(3,3)
        pos = eval('self.'+ant)

        if not antenna.full_evaluated:
            antenna.eval()


        if hasattr(antenna,'_mayamesh'):
            x, y, z, k, scalar = antenna._computemesh(T=rot,po=pos)
            antenna._mayamesh.mlab_source.set(x=x,y=y,z=z,scalars=scalar)
        else:
            antenna._show3(T=rot,po=pos,
                title=False,colorbar=False,newfig=False,name = '',interact=False)

        if delrays:
            import time
            for x in self._maya_fig.children[::-1]:
                if 'Rays' in x.name:
                    x.remove()
             # [x.remove() for x in self._maya_fig.children ]



if (__name__ == "__main__"):
    #plt.ion()
    doctest.testmod()
