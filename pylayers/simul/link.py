#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""

This module runs the electromagnetic simulation at the link level
It stores simulated objects in `hdf5` format.

DLink Class
==========

.. autosummary::
    :toctree: generated/


DLink simulation
---------------

.. autosummary::
    :toctree: generated/

    DLink.eval
    DLink._show3


DLink init
---------

.. autosummary::
    :toctree: generated/

    DLink.__init__
    DLink.__repr__
    DLink.reset_config
    DLink.fill_dexist


search in h5py file
-------------------

.. autosummary::
    :toctree: generated/

    DLink.checkh5
    DLink.array_exist
    DLink.get_grpname
    DLink.get_idx


Modify h5py file
----------------

.. autosummary::
    :toctree: generated/


    DLink.save_init
    DLink.save
    DLink.load
    DLink.stack
    DLink._delete




"""
import doctest
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
#from   Channel import *
import h5py
try:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi import mlab
except:
    print 'Layout:Mayavi is not installed'
import pdb



class Link(object):
    def __init__(self):
        """ Link evaluation class
        """
        pass

    def __add__(self,l):
        L = Link()
        tk = np.hstack((self.tk,l.tk))
        ak = np.hstack((self.ak,l.ak))
        us = np.argsort(tk)
        L.ak = ak[us]
        L.tk = tk[us]
        return L


class SLink(Link):

    def __init__(self):
        """ Statistical Link evaluation class
        """
        super(SLink,self).__init__()
        pass




class DLink(Link):

    def __init__(self, **kwargs):
        """ Deterministic Link evaluation

        Parameters
        ----------

        Basic

        L : Layout
            Layout to be used
        a : np.ndarray (3,)
            position of a device dev_a
        b : np.ndarray (3,)
            position of a device dev_b
        Aa : Antenna
            Antenna of device dev_a
        Ab : Antenna
            Antenna of device dev_b
        Ta : np.ndarray (3,3)
            Rotation matrice of Antenna of device dev_a relative to global Layout scene
        Tb : np.ndarray (3,3)
            Rotation matrice of Antenna of device dev_b relative to global Layout scene
        fGHz : np.ndarray (Nptf,)
            frequency range of Nptf poitns used for evaluation of channel in GHz
        wav : Waform
            Waveform to be applied on the channel
        save_idx : int
            number to differenciate the h5 file generated

        Advanced (change only if you really know what you do !)

        save_opt : list (['sig','ray','Ct','H'])
            information to be saved in the Links h5 file. Should never be Modified !

        force_create : Boolean (False)
            forcecreating the h5py file (if already exist, will be erased)


        Notes
        -----

        All simulation are stored into a unique file in your <PyProject>/output directory
        using the following convention:

        Links_<save_idx>_<LayoutFilename>.h5

        where
            <save_idx> is a integer number to be able to discriminate different links simulations
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
            uTb : indice of b position in 'T_map' Rotation dataset
            uAa : indice of a position in 'A_map' Antenna name dataset
            uAb : indice of b position in 'A_map' Antenna name dataset



        Examples
        --------

            >>> from pylayers.simul.link import *
            >>> L=DLink(verbose=False)
            >>> aktk = L.eval()


        """


        super(DLink,self).__init__()

        defaults={ 'L':Layout(),
                   'a':np.array(()),
                   'b':np.array(()),
                   'Aa':Antenna(),
                   'Ab':Antenna(),
                   'Ta':np.eye(3),
                   'Tb':np.eye(3),
                   'fGHz':np.linspace(2, 11, 181, endpoint=True),
                   'wav':wvf.Waveform(),
                   'cutoff':3,
                   'save_opt':['sig','ray','Ct','H'],
                   'save_idx':0,
                   'force_create':False,
                   'verbose':True
                }


        specset = ['a','b','Aa','Ab','Ta','Tb','L','fGHz']

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

        self._Lname = self._L.filename


        ###########
        # init ant
        ###########
        self.tx = RadioNode(name = '',
                            typ = 'tx',
                            _fileini = 'radiotx.ini',
                            _fileant = self.Aa._filename
                            )

        self.rx = RadioNode(name = '',
                            typ = 'rx',
                            _fileini = 'radiorx.ini',
                            _fileant = self.Ab._filename,
                            )


        ##############
        #### init save
        ###############
        self.filename = 'Links_' + str(self.save_idx) + '_' + self._Lname + '.h5' 
        filenameh5 = pyu.getlong(self.filename,pstruc['DIRLNK'])
        # check if save file alreasdy exists
        if not os.path.exists(filenameh5) or force:
            print 'Links save file for ' + self.L.filename + ' does not exist.'
            print 'It is beeing created. You\'ll see that message only once per Layout'
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
            print('This is the first time the Layout is used. Graphs have to be build. Please Wait')
            self.L.build()
            self.L.dumpw()


        ###########
        # init pos & cycles
        ###########

        self.ca = 0
        self.cb = 1
        self.a = self.L.cy2pt(self.ca)
        self.b = self.L.cy2pt(self.cb)


        ###########
        # init freq
        ###########
        self.fmin = self.fGHz[0]
        self.fmax = self.fGHz[-1]
        self.fstep = self.fGHz[1]-self.fGHz[0]

        # self.checkh5()

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

    @Lname.setter
    def Lname(self,Lname):
        # change layout and build/load
        self._L = Layout(Lname)
        self._Lname = Lname
        self.reset_config()

    @a.setter
    def a(self,position):
        if not self.L.ptin(position):
            raise NameError ('point a is not inside the Layout')
        self._a = position
        self.ca = self.L.pt2cy(position)
        self.tx.position = position

    @b.setter
    def b(self,position):
        if not self.L.ptin(position):
            raise NameError ('point b is not inside the Layout')
        self._b = position
        self.cb = self.L.pt2cy(position)
        self.rx.position = position

    @Aa.setter
    def Aa(self,Ant):
        self.tx = RadioNode(name = '',
                            typ = 'tx',
                            _fileini = 'radiotx.ini',
                            _fileant = Ant._filename
                            )
        self._Aa = Ant

    @Ab.setter
    def Ab(self,Ant):
        self.rx = RadioNode(name = '',
                            typ = 'rx',
                            _fileini = 'radiorx.ini',
                            _fileant = Ant._filename,
                            )
        self._Ab = Ant

    @Ta.setter
    def Ta(self,orientation):
        self._Ta = orientation
        self.tx.orientation = orientation

    @Tb.setter
    def Tb(self,orientation):
        self._Tb = orientation
        self.rx.orientation = orientation

    @fGHz.setter
    def fGHz(self,freq):
        self._fGHz = freq
        self.fmin = freq[0]
        self.fmax = freq[-1]
        self.fstep = freq[1]-freq[0]


    def __repr__(self):
        """ __repr__
        """

        s = 'filename: ' + self.filename +'\n'

        s = s + 'Actual Link considered is:\n'
        s = s + '--------------------------\n'
        s = s + 'Layout: ' + self.Lname + '\n\n'
        s = s + 'Device a:  \n'
        s = s + '---------  \n'
        s = s + 'position: ' + str (self.a) + '\n'
        s = s + 'Antenna: ' + str (self.Aa._filename) + '\n'
        s = s + 'Antenna rotation matrice : \n ' + str (self.Ta) + '\n\n'
        s = s + 'Device b:  \n'
        s = s + '---------  \n'
        s = s + 'position: ' + str (self.b) + '\n'
        s = s + 'Antenna: ' + str (self.Ab._filename) + '\n'
        s = s + 'Antenna rotation matrice : \n ' + str (self.Tb) + '\n\n'
        s = s + 'Link evaluation information : \n '
        s = s + '------------------ \n'
        s = s + 'distance: ' + str(np.sqrt(np.sum((self.a-self.b)**2))) + 'm \n'
        s = s + 'delay:' + str(np.sqrt(np.sum((self.a-self.b)**2))/0.3) + 'ns\n'
        s = s + 'Frequency range :  \n'
        s = s + 'fmin (fGHz):' + str(self.fGHz[0]) +'\n'
        s = s + 'fmax (fGHz):' + str(self.fGHz[-1]) +'\n'
        s = s + 'fstep (fGHz):' + str(self.fGHz[1]-self.fGHz[0]) +'\n '

        return s



    def reset_config(self):
        """ reset configuration when new layout loaded
        """
        try:
            self.L.dumpr()
        except:
            self.L.build()
            self.L.dumpw()


        self.ca = 0
        self.cb = 1
        self.a = self.L.cy2pt(self.ca)
        self.b = self.L.cy2pt(self.cb)

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
        """ check existence of previous simulation run with the same parameters.


        Returns
        -------

        update self.dexist dictionnary

        """
        # get identifier groupname in h5py file
        self.get_grpname()
        # check if grpnamee exist in the h5py file
        [self.fill_dexist(k,self.dexist[k]['grpname'])   for k in self.save_opt]



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
        -----------

        key : string

        array : np.ndarray

        Returns:
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
            return np.array([sc[0]+1])
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
            print 'delete ',key , ' in ', grpname
            f.close()
        except:
            f.close()
            raise NameError('Link._delete: issue when deleteting in h5py file')



    def save(self,obj,key,grpname,force=False):
        """ Save a given object in the correct group

        Parameters
        ----------

        obj : Object
            (Signatures|Rays|Ctilde|Tchannel)
        key : string
            key of the h5py file
        grpname : string
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


        """

        obj._loadh5(self.filename,grpname)
        if self.verbose :
            print str(obj.__class__).split('.')[-1] + ' from '+ grpname + ' loaded'


    def get_grpname(self):
        """ Determine if the group name of the data regarding the given configuration

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
        farray = np.array(([self.fmin,self.fmax,self.fstep]))
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
        uAa_opt, uAa = self.get_idx('A_map',self.Aa._filename)
        # check existence of Antenna b (Ab) in h5py file
        uAb_opt, uAb = self.get_idx('A_map',self.Ab._filename)


        grpname = str(ua) + '_' + str(ub) + '_' + str(uf) + \
                  '_'  + str(uTa) + '_' + str(uTb) + \
                  '_'  + str(uAa) + '_' + str(uAb)

        self.dexist['H']['grpname']=grpname


    def fill_dexist(self,key,grpname):
        """Check if the data of a key with a given groupname
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
            f=h5py.File(lfilename,'r')
            if grpname in f[key].keys():
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
            array type to check existency
        tol : np.float64
            tolerance (in meter for key == 'p_map')

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
        """ Evaluate the link


        Parameters
        ----------

        force_save : boolean
            Force the computation (even if obj already exists)
            AND save (replace previous computations)
        si.algo : str ('old'|'new')
            signature.run algo type
        ra.ceil_height_meter : int
            rays.to3D ceil height in mteres
        ra.number_mirror_cf : int
            rays.to3D number of ceil/floor reflexions


        Returns
        -------

        ak : ndarray
            alpha_k
        tk : ndarray
            tau_k


        update self.ak and self.tk

        self.ak : ndarray
            alpha_k
        self.tk : ndarray
            tau_k


        Example
        -------

        .. plot::
            :include-source:

            >>> from pylayers.simul.link import *
            >>> L=DLink(verbose=False)
            >>> iaktk = L.eval()


        See Also
        --------

        pylayers.antprop.signature
        pylayers.antprop.rays

        """

        defaults={ 'output':['sig','ray','Ct','H'],
                   'si.algo':'old',
                   'ra.ceil_height_meter':3,
                   'ra.number_mirror_cf':1,
                   'force':False,
                   }
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key]=value

        self.checkh5()


        ############
        # Signatures
        ############
        Si=Signatures(self.L,self.ca,self.cb,cutoff=self.cutoff)

        if self.dexist['sig']['exist'] and not kwargs['force']:
            self.load(Si,self.dexist['sig']['grpname'])

        else :
            Si.run5(cutoff=self.cutoff,algo=kwargs['si.algo'])
            # save sig
            self.save(Si,'sig',self.dexist['sig']['grpname'],force = kwargs['force'])

        self.Si = Si



        ############
        # Rays
        ############
        R = Rays(self.a,self.b)

        if self.dexist['ray']['exist'] and not kwargs['force']:
            self.load(R,self.dexist['ray']['grpname'])

        else :
            # perform computation...
            r2d = Si.rays(self.a,self.b)
            R = r2d.to3D(self.L,H=kwargs['ra.ceil_height_meter'], N=kwargs['ra.number_mirror_cf'])
            R.locbas(self.L)
            # ...and save
            self.save(R,'ray',self.dexist['ray']['grpname'],force = kwargs['force'])

        self.R = R

        if self.R.nray == 0:
            raise NameError('No ray have been founded. Try to re-run the simulation with a higher S.cutoff ')



        ############
        # Ctilde
        ############
        C=Ctilde()

        if self.dexist['Ct']['exist'] and not kwargs['force']:
            self.load(C,self.dexist['Ct']['grpname'])

        else :
            R.fillinter(self.L)
            # Ctilde...
            C=R.eval(self.fGHz)
            # ...save Ct
            self.save(C,'Ct',self.dexist['Ct']['grpname'],force = kwargs['force'])

        self.C = C

        ############
        # H
        ############
        H=Tchannel()

        if self.dexist['H']['exist'] and not kwargs['force']:
            self.load(H,self.dexist['H']['grpname'])


        else :
            # Ctilde antenna
            Cl=C.locbas(Tt=self.Ta, Tr=self.Tb)
            #T channel
            H=C.prop2tran(a=self.Aa,b=self.Ab)
            self.save(H,'H',self.dexist['H']['grpname'],force = kwargs['force'])

        self.H = H

        self.ak = H.ak
        self.tk = H.tk

        return self.ak, self.tk

    def _show3(self,rays=True,newfig = False,**kwargs):
        """ display the simulation scene using Mayavi

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

        if 'centered' in kwargs:
            centered = kwargs['centered']
        else :
            centered = False

        if centered:
            pg =np.zeros((3))
            pg[:2]=self.L.pg


        Atx = self.tx.A
        Arx = self.rx.A
        Ttx = self.tx.orientation
        Trx = self.rx.orientation

        if centered :
            ptx = self.tx.position-pg
            prx = self.rx.position-pg
        else :
            ptx = self.tx.position
            prx = self.rx.position

        # evaluate antenna if required
        if not Atx.evaluated:
            Atx.Fsynth()
        elif len(Atx.SqG.shape) == 2 :
            Atx.Fsynth()

        if not Arx.evaluated:
            Arx.Fsynth()
        elif len(Arx.SqG.shape) == 2 :
            Arx.Fsynth()

        self.L._show3(newfig=False,opacity=0.7,centered=centered)


        Atx._show3(T=Ttx.reshape(3,3),po=ptx,
            title=False,colorbar=False,newfig=False)
        Arx._show3(T=Trx.reshape(3,3),po=prx,
            title=False,colorbar=False,newfig=False,name = '')
        if rays :
            try:
                self.R._show3(**kwargs)
            except:
                print 'Rays not computed yet'

if (__name__ == "__main__"):
    #plt.ion()
    doctest.testmod()
