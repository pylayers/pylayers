#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""

This module runs the electromagnetic simulation at the link level
It stores simulated objects in `hdf5` format.

Link Class
==========

.. autosummary::
    :toctree: generated/


Link simulation
---------------

.. autosummary::
    :toctree: generated/

    Link.eval
    Link._show3


Link init
---------

.. autosummary::
    :toctree: generated/

    Link.__init__
    Link.__repr__
    Link.reset_config
    Link.fill_dexist


search in h5py file
-------------------

.. autosummary::
    :toctree: generated/

    Link.checkh5
    Link.array_exist
    Link.get_grpname
    Link.get_idx


Modify h5py file
----------------

.. autosummary::
    :toctree: generated/


    Link.save_init
    Link.save
    Link.load
    Link.stack
    Link._delete




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

    def __init__(self, **kwargs):
        """
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

        #.. plot::
        #    :include-source:
        #
        #    >>> from pylayers.simul.link import *
        #    >>> L=Link()
        #    >>> L.eval()
        #    >>> L._show3()


        """

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


        # if len(self.a) == 0 :
        #     self.a = self.L.cy2pt(0)

        # if len(self.b) == 0 :
        #     self.b = self.L.cy2pt(1)

        # if not self.L.ptin(self.a):
        #     raise NameError ('point a not inside the Layout')
        # if not self.L.ptin(self.b):
        #     raise NameError ('point b not inside the Layout')

        # self.ca = self.L.pt2cy(self.a)
        # self.cb = self.L.pt2cy(self.b)




        # self.fmin = self.fGHz[0]
        # self.fmax = self.fGHz[-1]
        # self.fstep = self.fGHz[1]-self.fGHz[0]

        # # update radio node

        # self.tx.position = self.a
        # self.rx.position = self.b
        # self.tx.orientation = self.Ta
        # self.rx.orientation = self.Tb
        # self.tx.A = self.Aa
        # self.rx.A = self.Ab

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

        a ,t

        a : ndarray
            alpha_k
        t : ndarray
            tau_k

        Friss is applyed on H 

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
        self.H.applyFriis()
        a = np.real(np.sqrt(np.sum(self.H.y * np.conj(self.H.y), axis=1))
                                                             / len(self.H.y))
        t = H.tau0

        return a, t

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

        #.. plot::
        #    :include-source:

        #    >>> from pylayers.simul.link import *
        #    >>> L=Link()
        #    >>> L.eval()
        #    >>> L._show3()

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
    # def save(self):
    #     """ save simulation file

    #     Simulation files are .ini files which are saved in a dedicated
    #     directory basename/ini in the Project tree

    #     """
    #     filesimul = pyu.getlong(self.filesimul, "ini")
    #     fd = open(filesimul, "w")
    #     # getting current spa file if any
    #     try:
    #         self.config.set("files", "tx", self.tx.fileini)
    #     except:
    #         pass
    #     try:
    #         self.config.set("files", "rx", self.rx.fileini)
    #     except:
    #         pass
    #     self.config.write(fd)
    #     fd.close()
    #     # save tx
    #     self.tx.save()
    #     # save rx
    #     self.rx.save()
    #     # save slab and mat file
    #     # --
    #     # self.slab.save(self.fileslabini)
    #     # self.slab.mat.save(self.filematini)
    #     # --
    #     # fix bug #189
    #     #   slab is a member of S.L not of S anymore
    #     #   mat is a member of S.L.sl not of S.slab
    #     try:
    #         self.L.sl.save(self.fileslabini)
    #     except:
    #         pass
    #     try:
    #         self.L.sl.mat.save(self.filematini)
    #     except:
    #         pass

#    def saveold(self):
#        """ save simulation file

#        """
#        filesimul = pyu.getlong(self.filesimul, "ini")
#        fd = open(filesimul, "w")
#        config = ConfigParser.ConfigParser()

#        #
#        # files section
#        #
#        #config.add_section("files")
#        self.config.set("files", "conf", self.fileconf)
#        self.config.set("files", "struc", self.filestr)
#        self.config.set("files", "slab", self.fileslab)
#        self.config.set("files", "mat", self.filemat)

#        try:
#            self.config.set("files", "tx", self.tx.filespa)
#        except:
#            pass
#        try:
#            self.config.set("files", "rx", self.rx.filespa)
#        except:
#            pass
#        try:
#            self.config.set("files", "txant", self.tx.fileant)
#        except:
#            pass
#        try:
#            self.config.set("files", "rxant", self.rx.fileant)
#        except:
#            pass
#        self.config.set("files", "patra", self.filepatra)
#        self.config.set("files", "palch", self.filepalch)
#        self.palch.save()
#        self.patra.save()

#        #
#        # tud section
#        #

#        self.config.set("tud", "purc", self.patud.purc)
#        self.config.set("tud", "nrmax", self.patud.nrmax)
#        self.config.set("tud", "num", self.patud.num)

#        #
#        # frequency section
#        #
#        self.config.set("frequency", "fghzmin", self.freq[0])
#        self.config.set("frequency", "fghzmax", self.freq[-1])
#        self.config.set("frequency", "Nf", len(self.freq))

#        #
#        # output section
#        #
#        #filelch exists
#        if self.progress > 0:
#            #config.add_section("output")
#            for k in range(len(self.filelch)):
#                _fileout = "out" + "???"
#                filename = self.filelch[k]
#                self.config.set("launch", str(k + 1), filename)

#        # filetra exists
#        for k in range(len(self.filelch)):
#            if self.progress > 1:
#                #self.config.add_section("trace")
#                for l in arange(len(self.filetra[k])):
#                    filename = self.filetra[k][l]
#                    self.config.set("trace", "rx" + str(l + 1), filename)

#        # .tang exists
#        # .rang exists
#        # .tud exists
#            if self.progress > 2:
#                #config.add_section("tud")
#                #config.add_section("tang")
#                #config.add_section("rang")

#                for l in arange(len(self.filetud[k])):
#                    ftud = self.filetud[k][l]
#                    self.config.set("tud", "rx" + str(l + 1), ftud)

#                for l in arange(len(self.filetang[k])):
#                    ftang = self.filetang[k][l]
#                    self.config.set("tang", "rx" + str(l + 1), ftang)

#                for l in arange(len(self.filerang[k])):
#                    frang = self.filerang[k][l]
#                    self.config.set("rang", "rx" + str(l + 1), frang)

#        # .field exist
#        # .tauk exist
#            if self.progress > 3:
#                #config.add_section("tauk")
#                #config.add_section("field")
#                for l in arange(len(self.filetud[k])):
#                    ftauk = self.filetud[k][l]
#                    self.config.set("tauk", "rx" + str(l + 1), ftauk)

#                for l in arange(len(self.filefield[k])):
#                    ffield = self.filefield[k][l]
#                    self.config.set("field", "rx" + str(l + 1), ffield)

#        self.config.write(fd)
#        fd.close()

#     def save_project(self):
#         """ save Simulation files in a zipfile

#         Simulation files are .ini files which are saved in a dedicated
#         directory basename/ini in the Project tree

#         """
#         root = Tkinter.Tk()
#         zipfileName = tkFileDialog.asksaveasfilename(parent=root,
#                             filetypes = [("zipped file","zip")] ,
#                             title="Save Project",
#                             )
#         pyu.zipd(basename,zipfileName)
#         root.withdraw()
#         print "Current project saved in", zipfileName

#     def load_project(self):
#         """ load Simulation files from a zipfile

#         Simulation files are .ini files which are saved in a dedicated
#         directory basename/ini in the Project tree

#         """
#         root = Tkinter.Tk()
#         zipfileName= tkFileDialog.askopenfile(parent=root,
#                                             mode='rb',
#                                             title='Choose a project')
#         dirname = tkFileDialog.askdirectory(parent=root,
#                                     initialdir=basename,
#                                     title='Please select a directory',
#                                     mustexist=0)
#         pyu.unzipd(dirname,zipfileName)
#         root.withdraw()
#         print "Current project loaded in", dirname

#     def choose(self):
#         """
#             Choose a simulation file in simuldir

#         """
#         import tkFileDialog as FD
#         fichsimul = FD.askopenfilename(filetypes=[("Fichiers simul ",
#                                                    "*.simul"),
#                                                   ("All", "*")],
#                                        title="Please choose a simulation file",
#                                        initialdir=simuldir)
#         self.filesimul = pyu.getshort(fichsimul)
#         self.load()

#     def load(self, _filesimul):
#         """ load a simulation configuration file

#          each transmiter simulation results in the creation of an .ini file
#          with the following sections
#          related to the results obtained for different receivers

#         Parameters
#         ----------

#         _filesimul   : file in the simul directory of the Project

#         """

#         self.filesimul = _filesimul
#         filesimul = pyu.getlong(self.filesimul, "ini")

#         self.config.read(filesimul)

#         sections = self.config.sections()
#         try:
#             _filetx = self.config.get("files", "tx")
#         except:
#             raise NameError('Error in section tx from '+ _filesimul)

#         try:
#             _filerx = self.config.get("files", "rx")
#         except:
#             raise NameError('Error in section rx from '+ _filesimul)


#         try:
#             _fileanttx = self.config.get("files", "txant")
#         except:
#             raise NameError('Error in section txant from '+ _filesimul)

#         try:
#            _fileantrx = self.config.get("files", "rxant")
#         except:
#             raise NameError('Error in section rxant from '+ _filesimul)

#         try:
#             self.filestr = self.config.get("files", "struc")
#             # force .str extension
#             f,e = self.filestr.split['.']
#             self.filestr = f+'.str'
#             print self.filestr
#         except:
#             raise NameError('Error in section struc from '+ _filesimul)

#         try:
#             self.tx = RadioNode(name = '',
#                                 typ = 'tx',
#                                 _fileini = _filetx,
#                                 _fileant = _fileanttx,
#                                 _filestr = self.filestr)

#             self.rx = RadioNode(name = '',
#                                 typ = 'rx',
#                                 _fileini = _filerx,
#                                 _fileant = _fileantrx,
#                                 _filestr = self.filestr)
#         except:
#             raise NameError('Error during Radionode load')
# #
# # Launching and Tracing parameters
# #

#         try:
#             self.palch = Palch(self.config.get("files", "palch"))
#             self.filepalch = self.config.get("files", "palch")
#         except:
#             raise NameError('Error in section palch from '+ _filesimul)

#         try:
#             self.patra = Patra(self.config.get("files", "patra"))
#             self.filepatra = self.config.get("files", "patra")
#         except:
#             raise NameError('Error in section patra from '+ _filesimul)
#         #_outfilename = "out"+self.ntx+".ini"
#         #self.outfilename = pyu.getlong(_outfilename,"simul")
# #  Load Simulation Mat File
# #
#         self.filematini = self.config.get("files", "mat")
#         self.mat = MatDB()
#         self.mat.load(self.filematini)
# #
# #  Load Simulation Slab File
# #
#         try:
#             self.fileslabini = self.config.get("files", "slab")
#             self.sl = SlabDB()
#             self.sl.mat = self.mat
#             self.sl.load(self.fileslabini)
#         except:
#             raise NameError('Slab load error')
# #
# # Load layout from .str or .str2 file
# #
#         try:
#             self.L = Layout(self.filestr,self.filematini, self.fileslabini)
# #            self.L.load(self.filestr)
#         except:
#             raise NameError('Layout load error')

#         try:
#             self.patud.nrmax = self.config.get("tud", "nrmax")
#             self.patud.num = self.config.get("tud", "num")
#             self.patud.purc = self.config.get("tud", "purc")
#         except:
#             pass

# #
# # Frequency base
# #
#         if "frequency" in sections:
#             try:
#                 self.fGHz = np.linspace(float(self.config.getfloat("frequency", "fghzmin")),
#                                         float(self.config.getfloat("frequency", "fghzmax")),
#                                         int(self.config.getint("frequency", "nf")),
#                                         endpoint=True)
#             except:
#                 raise NameError('Error in section frequency from '+ _filesimul)
#             # update .freq file in tud directory

#             filefreq = pyu.getlong(self.filefreq, pstruc['DIRTUD'])
#             fd = open(filefreq, "w")
#             chaine = self.config.get("frequency", "fghzmin") + ' ' + \
#                 self.config.get("frequency", "fghzmax") + ' ' + \
#                 self.config.get("frequency", "nf")
#             fd.write(chaine)
#             fd.close
# #
# # Simulation Progress
# #
#         self.output = {}
#         if "output" in sections:
#             for itx in self.config.options("output"):
#                 _filename  =  self.config.get("output", itx)
#                 self.dout[int(itx)] = _filename
#                 filename = pyu.getlong(_filename, "output")
#                 output = ConfigParser.ConfigParser()
#                 output.read(filename)
#                 secout = output.sections()
#                 self.dtra[int(itx)] = {}
#                 self.dtud[int(itx)] = {}
#                 self.dtang[int(itx)] = {}
#                 self.drang[int(itx)] = {}
#                 self.dtauk[int(itx)] = {}
#                 self.dfield[int(itx)] = {}
#                 self.dcir[int(itx)] = {}
#                 if "launch" in secout:
#                     self.progress = 1
#                     keys_launch = output.options("launch")
#                     for kl in keys_launch:
#                         self.dlch[int(kl)] = output.get("launch", kl)
#                 if "trace" in secout:
#                     self.progress = 2
#                     keys_tra = output.options("trace")
#                     for kt in keys_tra:
#                         self.dtra[int(itx)][int(kt)] = output.get("trace", kt)

#                 if "tang" in secout:
#                     self.progress = 3
#                     keys_tang = output.options("tang")
#                     for kt in keys_tang:
#                         self.dtang[int(itx)][int(kt)] = output.get("tang", kt)
#                         self.drang[int(itx)][int(kt)] = output.get("rang", kt)
#                         self.dtud[int(itx)][int(kt)] = output.get("tud", kt)
#                 if "field" in secout:
#                     self.progress = 4
#                     keys_field = output.options("field")
#                     for kt in keys_field:
#                         self.dfield[int(itx)][int(kt)] = output.get(
#                             "field", kt)
#                         self.dtauk[int(itx)][int(kt)] = output.get("tauk", kt)
#                 if "cir" in secout:
#                     self.progress = 5
#                     keys_cir = output.options("cir")
#                     for kt in keys_cir:
#                         self.dcir[int(itx)][int(kt)] = output.get("cir", kt)

#                 self.output[int(itx)] = output
#         #
#         # Waveform section
#         #
#         self.wav = wvf.Waveform()
#         self.wav.read(self.config)

#     def layout(self, _filestruc, _filematini='matDB.ini', _fileslabini='slabDB.ini'):
#         """ load a layout in the simulation oject

#         Parameters
#         ----------

#         _filestruc : string
#             short file name of the Layout object
#         _filematini   : string
#             short file name of the Mat object  (default matDB.ini)
#         _fileslab  : string
#             short file name of the Slab object (default slabDB.ini)

#         Examples
#         --------

#         >>> from pylayers.simul.simulem import *
#         >>> S = Simul()
#         >>> S.layout('defstr.str')

#         """
#         self.filestr = _filestruc
#         self.filematini = _filematini
#         self.fileslabini = _fileslabini

#         self.L = Layout(_filestruc,_filematini, _fileslabini)
#         # update config
#         self.config.set("files", "struc", self.filestr)
#         self.config.set("files", "slab", self.fileslabini)
#         self.config.set("files", "mat", self.filematini)
#         self.save()

#     def show(self, itx=[-1], irx=[-1], furniture=True, s=8, c='b', traj=False, num=False,fig=[],ax=[]):
#         """ show simulation

#             Parameters
#             -----------
#             itx        : list of tx indexes
#             irx        : list of rx indexes
#             furniture  : boolean for METALIC furniture display
#             s          : scale fir scatter plot  (default 8)
#             c          : color for scatter plot  (default 'b')
#             traj       : boolean  (def False)
#             num        : boolean  (def False)
#                 display a number


#             Examples
#             --------
#             >>> import matplotlib.pyplot as plt
#             >>> from pylayers.simul.simulem import *
#             >>> S = Simul()
#             >>> S.load('w1.ini')
#             >>> S.L.loadfur('FurW1.ini')
#             >>> S.show()
#             >>> plt.show()



#         """
#         if type(itx) == int:
#             itx = [itx]
#         if type(irx) == int:
#             irx = [irx]


#         if fig ==[]:
#             fig = plt.gcf()
#         if ax==[]:
#             ax = fig.gca()

#         #self.L.display['scaled']=False
#         fig,ax=self.L.showGs(fig=fig,ax=ax, show=False)
#         #
#         if furniture:
#             if 'lfur' in self.L.__dict__:
#                 for fur in self.L.lfur:
#                     if fur.Matname == 'METAL':
#                         fur.show(fig, ax)
#             else:
#                 print "Warning : no furniture file loaded"

#         if irx[0] == -1:
#             ax.scatter(self.rx.position[0,:],
#                        self.rx.position[1,:], c='b', s=s, alpha=0.5)
#             #ax.scatter(self.rx.position[0,0],self.rx.position[0,1],c='k',s=s,linewidth=0)
#             #ax.scatter(self.rx.position[1,0],self.rx.position[1,1],c='b',s=s,linewidth=0)
#             #ax.scatter(self.rx.position[2,0],self.rx.position[2,1],c='g',s=s,linewidth=0)
#             #ax.scatter(self.rx.position[3,0],self.rx.position[3,1],c='c',s=s,linewidth=0)
#         else:
#             for k in irx:
#                 ax.scatter(self.rx.position[0,k - 1],
#                            self.rx.position[1,k - 1], c='b', s=s, alpha=0.5)
#                 if num:
#                     ax.text(self.rx.position[0,k - 1],
#                             self.rx.position[1,k - 1],
#                             str(k), color='blue')

#         if itx[0] == -1:
#             ax.scatter(self.tx.position[0,:],
#                        self.tx.position[1,:], c='r', s=s)
#             if num:
#                 for k in range(302):
#                     ax.text(self.tx.position[0,k - 1],
#                             self.tx.position[1,k - 1],
#                             str(k), color='black')
#         else:
#             if traj:
#                 cpt = 1
#             for k in itx:
#                 ax.scatter(self.tx.position[0,k - 1],
#                            self.tx.position[1,k - 1],
#                            c=c, s=s, linewidth=0)
#                 if num:
#                     if traj:
#                         ax.text(self.tx.position[0,k - 1],
#                                 self.tx.position[1,k - 1],
#                                 str(cpt), color='black')
#                         cpt = cpt + 1
#                     else:
#                         ax.text(self.tx.position[0,k - 1],
#                                 self.tx.position[1,k - 1],
#                                 str(k), color='black')

#         return (fig,ax)
#         #for k in range(self.tx.N):
#         #    ax.text(self.tx.position[0,k],self.tx.position[1,k],str(k+1),color='black')
#         #    ax.scatter(self.tx.position[0,:],self.tx.position[0,:],

#         #for k in range(self.rx.N):
#         #   ax.text(self.rx.position[0,k],self.rx.position[1,k],str(k),color='black')

#     def PL(self, itx):
#         """ plot Path Loss

#         itx
#         """
#         td = []
#         tEa = []
#         tEo = []
#         for irx in self.dcir[itx].keys():
#             d = self.delay(itx, irx) * 0.3
#             cira, ciro = self.loadcir(itx, irx)
#             Ea = cira.energy()[0]
#             Eo = ciro.energy()[0]
#             td.append(d)
#             tEa.append(Ea)
#             tEo.append(Eo)

#         plt.semilogx(td, 10 * np.log10(tEa), 'xr')
#         plt.semilogx(td, 10 * np.log10(tEo), 'xb')
#         plt.show()
#         return td, tEa, tEo

#     def evalcir(self,cutoff=4,algo='new'):
#         """
#         Parameters
#         ----------

#         S
#         tx
#         rx
#         wav
#         cutoff

#         """

#         crxp =-1
#         ctxp =-1
#         tcir = {}
#         tx = self.tx.position
#         Ntx = len(tx[0])
#         rx = self.rx.position
#         Nrx = len(rx[0])

#         #for kt in range(1,Ntx-1):
#         #print kt+1
#         kt=0
#         tcir[kt] = {}
#         t = np.array([self.tx.position[0,kt],self.tx.position[1,kt],self.tx.position[2,kt]])
#         for kr in range(Nrx):
#             if (np.mod(kr,10)==0):
#                 print kr+1
#             r = np.array([self.rx.position[0,kr],self.rx.position[1,kr],self.rx.position[2,kr]])
#             ctx = self.L.pt2cy(t)
#             crx = self.L.pt2cy(r)
#             if (ctx<>ctxp)|(crx<>crxp):
#                 Si  = signature.Signatures(self.L,ctx,crx)
#                 ctxp = ctx
#                 crxp = crx
#                 Si.run4(cutoff=cutoff,algo=algo)
#             r2d = Si.rays(t,r)
#             #r2d.show(S.L)

#             r3d = r2d.to3D(self.L)
#             r3d.locbas(self.L)
#             r3d.fillinter(self.L) 
#             Ct  = r3d.eval(self.fGHz)
#             sca = Ct.prop2tran(self.tx.A,self.rx.A)
#             cir = sca.applywavB(self.wav.sfg)
#             tcir[kt][kr] = cir
#         return(tcir)

#     def loadcir(self, itx, irx):
#         """
#         Parameters
#         ----------
#         itx : Tx index
#         irx : Rx index

#         Returns
#         -------

#         cir(itx,irx)
#         """
#         _filecir = self.dcir[itx][irx] + '.mat'
#         ext = str(itx)
#         if len(ext) == 1:
#             ext = '00' + ext
#         if len(ext) == 2:
#             ext = '0' + ext

#         filecir = pyu.getlong(_filecir, pstruc['DIRCIR']+'/Tx' + ext)
#         D = spio.loadmat(filecir)

#         kxa = 'ta' + str(irx)
#         kya = 'cira' + str(irx)

#         kxo = 'to' + str(irx)
#         kyo = 'ciro' + str(irx)

#         cira = bs.TUsignal(D[kxa], D[kya][:, 0])
#         ciro = bs.TUsignal(D[kxo], D[kyo][:, 0])

#         return(cira, ciro)

#     def pltcir(self, itx=1, irx=1, mode='linear', noise=False, color='b',format='a',fig=[],ax=[]):
#         """ plot Channel Impulse Response

#         Parameters
#         ----------
#         itx : Tx index
#         irx : Rx index
#         mode : str
#             {'linear','dB'}
#             noise : boolean
#         color : string
#             default 'b'

#         >>> from pylayers.simul.simulem import *
#         >>> S = Simul()
#         >>> S.load('where2.ini')
#         >>> S.run(1,1)
#         >>> S.pltcir(1,1,mode='linear',noise=False,color='k')

#         """

#         if fig ==[]:
#             fig = plt.gcf()
#         #if ax==[]:
#         #    ax = fig.gca()

#         _filecir = self.dcir[itx][irx] + '.mat'
#         filecir = pyu.getlong(_filecir, pstruc['DIRCIR']+'/Tx' + str('%0.3d' % itx))
#         D = spio.loadmat(filecir)
#         ax = fig.add_subplot('211')

#         fig,ax=self.show(itx, irx,fig=fig,ax=ax)
#         ax=fig.add_subplot('212')
#         if 'a' in format :
#             kxa = 't'
#             kya = 'cir'
#             ta = D[kxa]
#             Tobs = ta[-1] - ta[0]
#             te = ta[1] - ta[0]
#             if noise:
#                 na = bs.Noise(Tobs + te, 1. / te)
#                 naf = na.gating(4.493, 0.5)
#             cira = bs.TUsignal(ta, D[kya][:, 0])


#         if 'o' in format:
#             kxo = 'to' + str(irx)
#             kyo = 'ciro' + str(irx)
#             to = D[kxo]
#             Tobs = to[-1] - to[0]
#             te = to[1] - to[0]
#             if noise:
#                 no = bs.Noise(Tobs + te, 1. / te)
#                 nof = no.gating(4.493, 0.5)
#             ciro = bs.TUsignal(to, D[kyo][:, 0])

#         if mode == 'linear':
#             #plt.plot(ta,naf.y,color='k',label='Noise')
#             plt.plot(ta, D[kya], label='Rx ' + str(irx), color=color)
#             plt.xlabel('Time (ns)')

#             '''if noise:
#                 naf.plot(col='k')
#             cira.plot(col=color)'''
#         else:
#             '''if noise:
#                 naf.plotdB(col='k')
#             cira.plotdB()'''
#             plt.plot(ta, 20 * np.log10(abs(D[kya])), label='Rx ' + str(irx), color=color)
#             plt.xlabel('Time (ns)')
# #        plt.legend()
#         plt.show()
#         #plt.savefig('Tx'+str(itx),format=pdf,dpi=300)

#     def scatter(self, itx, irx, values,
#                 cmap=plt.cm.gray,
#                 s=30,
#                 spl=221,
#                 title='',
#                 vaxis=((-30, 10, 2, 18)),
#                 vmin=0,
#                 vmax=1,
#                 colbool=False,
#                 cblabel='dB'):
#         """
#             Parameters
#             ----------
#             itx
#             irx
#             values
#             cmap
#             s
#             spl
#             title
#             vaxis
#             vmin
#             vmax
#             colbool
#             clabel

#         """
#         fig = plt.gcf()
#         ax = fig.add_subplot(spl)
#         xtx = self.tx.position[itx, 0]
#         ytx = self.tx.position[itx, 1]
#         xrx = self.rx.position[irx, 0]
#         yrx = self.rx.position[irx, 1]
#         self.L.display['title'] = title
#         self.L.showGs(ax)
#         for furk in siradel.siradel_furniture.keys():
#             fur = siradel.siradel_furniture[furk]
#             if fur.Matname == 'METAL':
#                 fur.show(fig, ax)

#         #self.show(furniture=True)
#         plt.axis(vaxis)
#         b1 = ax.scatter(xtx, ytx, s=s, c=values, cmap=cmap,
#                         linewidths=0, vmin=vmin, vmax=vmax)
#         ax.scatter(xrx, yrx, s=30, c='b', linewidths=0)
#         if colbool:
#             cb = colorbar(b1)
#             cb.set_label(cblabel, fontsize=14)
#         return(b1)

#     def info(self, itx=[], irx=[]):
#         """ display simulation information

#          Parameters
#          ----------
#          itx : Tx index
#          irx : Rx index

#         """
#         print self.filesimul
#         print '------------------------------------------'
#         try:
#             print "Layout Info : \n", self.L.info()
#         except:
#             print "provide a Layout in the simulation : S.L "
#             print ">>> S.layout(filename.str) "
#             print "or "
#             print ">>> S.layout(filename.str2 "
#             print "or "
#             print ">>> S.layout(filename.str,filematini,filematini) "
#             print "default files exists for filematini and fileslabini "

#             return
#         try:
#             print "Tx Info :\n", self.tx.info()
#         except:
#             print "provide a tx in the simulation : S.tx "
#             return
#         try:
#             print "Rx Info :\n", self.rx.info()
#         except:
#             print "provide a rx in the simulation : S.rx "
#             return

#         #    print "Tx : ",self.tx.points[itx]
#             print "Rx : ", self.rx.points[itx]
#             print "Delay (ns) :", self.delay(itx, irx)
#             print "Distance (m) :", 0.3 / self.delay(itx, irx)
#             print ""
#             if itx in self.dlch.keys():
#                 print "-----"
#                 print "Launching "
#                 print "-----"
#                 print " ", self.dlch[itx]
#             if irx in self.dtra[itx].keys():
#                 print "-----"
#                 print "Tracing "
#                 print "-----"
#                 print " ", self.dtra[itx][irx]
#                 gr = GrRay3D()
#                 gr.load(self.dtra[itx][irx], self.L)
#                 gr.info()
#             if irx in self.dtud[itx].keys():
#                 print "-----"
#                 print "Tud parameters "
#                 print "-----"
#                 print " ", self.dtud[itx][irx]
#                 print " ", self.dtang[itx][irx]
#                 print " ", self.drang[itx][irx]
#                 gt = GrRay3D.GrRayTud()
#                 gt.load(self.dtud[itx][irx],
#                         self.dtang[itx][irx],
#                         self.drang[itx][irx], self.sl)
#             if irx in self.dtauk[itx].keys():
#                 print self.dtauk[itx][irx]
#                 print self.dfield[itx][irx]
#                 VC = self.VC(itx, irx)
#             if irx in self.dcir[itx].keys():
#                 print self.dcir[itx][irx]

#     def info2(self):
#         for i, j in enumerate(self.__dict__.keys()):
#             print j, ':', self.__dict__.values()[i]

#     def filtray(self, itx, irx, tau0, tau1, col='b'):
#         """ filter rays

#         Parameters
#         ----------
#         itx :
#         irx :
#         tau0 :
#         tau1 :
#         col :

#         Display ray and nstr
#         """
#         gr = GrRay3D()
#         gr.load(self.dtra[itx][irx], self.L)
#         self.L.display['Thin'] = True
#         self.L.display['Node'] = False
#         self.L.display['NodeNum'] = False
#         self.L.display['EdgeNum'] = False
#         plt.axis('scaled')
#         #self.L.show(fig,ax,nodelist,seglist)
#         self.L.showGs()
#         delays = gr.delay()
#         rayset = np.nonzero((delays >= tau0) & (delays <= tau1))[0]
#         fig = plt.gcf()
#         ax = fig.get_axes()[0]
#         gr.show(ax, rayset, col=col, node=False)
#         plt.title('Tx' + str(itx) + '-Rx' + str(irx) + ' : ' + str(
#             tau0) + ' < tau < ' + str(tau1))

#     def showray(self, itx, irx, iray=np.array([]), fig=[], ax=[]):
#         """ show layout and rays for a radio link

#         Parameters
#         ----------
#         itx  : tx index
#         irx  : rx index
#         iray : list of rays to be displayed ndarray


#         """

#         #if ax==[]:
#         #    fig = plt.figure()
#         #    ax  = fig.add_subplot('111')

#         gr = GrRay3D()
#         gr.load(self.dtra[itx][irx], self.L)
#         if len(iray == 1):
#             ray = gr.ray3d[iray[0]]
#             nstr = ray.nstr[1:-1]
#             uneg = np.nonzero(nstr < 0)
#             upos = np.nonzero((nstr > 0) & (nstr <= self.L.Ne))
#             uceil = np.nonzero(nstr == self.L.Ne + 1)
#             ufloor = np.nonzero(nstr == self.L.Ne + 2)
#         #seglist  = nstr[upos[0]]-1
#         #nodelist = -nstr[uneg[0]]-1
#         #seglist2 = S.L.segpt(nodelist)
#         self.L.display['Thin'] = True
#         self.L.display['Node'] = False
#         self.L.display['NodeNum'] = False
#         self.L.display['EdgeNum'] = False
#         #self.L.show(fig,ax)
#         #print ray.nn
#         #print len(ray.nstr)
#         #print ray.nstr
#         #print nodelist
#         #print seglist
#         #print seglist2
#         #seglist = hstack((seglist,seglist2))
#         self.L.display['Node'] = False
#         self.L.display['Thin'] = False
#         self.L.display['NodeNum'] = True
#         self.L.display['EdgeNum'] = True
#         plt.axis('scaled')
#         #self.L.show(fig,ax,nodelist,seglist)
#         fig, ax = self.L.showGs(show=False)
#         gr.show(ax, iray, col='b', node=False)

#         if len(iray) == 1:
#             plt.title(str(nstr))
#         else:
#             plt.title('Tx' + str(itx) + '-Rx' + str(irx) +
#                       ' ' + str(min(iray)) + ' ' + str(max(iray)))

#     def show3l(self, itx, irx):
#         """ geomview display of a specific link

#         g = S.show3l(itx,irx)

#         Parameters
#         ----------
#         itx
#             transmitter index
#         irx
#             receiver index

#         """
#         filetra = self.dtra[itx][irx]
#         gr = GrRay3D()
#         gr.load(filetra, self.L)
#         gr.show3()

#         return(gr)



    # def show3(self,rays=[],**kwargs):
    #     """ geomview display of the simulation configuration

    #     Parameters
    #     ----------

    #     centered : boolean
    #         center the scene if True
    #     bdis  : boolean
    #         display local basis

    #     """
    #     try:
    #         self.tx.save()
    #     except:
    #         print('tx set is no defined')
    #     try:
    #         self.rx.save()
    #     except:
    #         print('rx set is no defined')
    #     _filename = self.filesimul.replace('.ini', '.off')
    #     filename = pyu.getlong(_filename, pstruc['DIRGEOM'])
    #     fo = open(filename, "w")
    #     fo.write("LIST\n")
    #     try:
    #         sttx = "{<" + self.tx.filegeom + "}\n"
    #     except:
    #         sttx = "\n"
    #     try:
    #         strx = "{<" + self.rx.filegeom + "}\n"
    #     except:
    #         strx = "\n"
    #     try:
    #         stst = "{<" + self.L.filegeom + "}\n"
    #     except:
    #         stst = "\n"
    #     fo.write(sttx)
    #     fo.write(strx)
    #     fo.write(stst)
    #     fo.write("{</usr/share/geomview/geom/xyz.vect}\n")
    #     if rays !=[]:
    #         kwargs['bdis']=False
    #         kwargs['L']=self.L
    #         kwargs['centered']=False
    #         fo.write("{<" + rays.show3(**kwargs) + "}")

    #     fo.close()


    #     command = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
    #     os.system(command)

    # def freq(self, GUI=False):
    #     """
    #         return the frequency base from the content of filefreq

    #     Parameters
    #     ----------
    #     GUI
    #         Boolean


    #     """
    #     filefreq = self.filefreq
    #     filefreq = pyu.getlong(filefreq, pstruc['DIRTUD'])
    #     fo = open(filefreq)
    #     l = fo.readline().split()
    #     fmin = eval(l[0])
    #     fmax = eval(l[1])
    #     Nf = eval(l[2])
    #     fo.close()
    #     if GUI:
    #         val = multenterbox('Enter frequency (GHz)', '',
    #                            ('start', 'stop', 'N'),
    #                            (str(fmin), str(fmax), str(Nf)))
    #         fmin = eval(val[0])
    #         fmax = eval(val[1])
    #         Nf = eval(val[2])
    #         fo = open(filefreq, 'w')
    #         data = val[0] + ' ' + val[1] + ' ' + val[2] + '\n'
    #         fo.write(data)
    #         fo.close()

    #     fGHz = np.linspace(fmin, fmax, Nf, endpoint=True)
    #     return(fGHz)

    # def getlaunch(self, k=1):
    #     """
    #      get the kth launch

    #     Parameters
    #     ----------
    #     k :  int
    #         launching index (default 1)

    #     Launching index starts at 1

    #     """
    #     if k in self.dlch.keys():
    #         filename = self.dlch[k]
    #         L = Launch()
    #         L.load(filename)
    #         return(L)
    #     else:
    #         print("Tx not available")

    # def gettra(self, itx, irx):
    #     """
    #      Gtra = S.gettra(itx,irx)
    #     """
    #     if (itx < self.tx.N) & (irx < self.rx.N):
    #         Gtra = GrRay3D()
    #         Gtra.load(self.filetra[itx][irx], self.indoor)
    #         return Gtra
    #     else:
    #         print "gettra warning : wrong tx or rx index"

    # def gettud(self, itx, irx):
    #     """
    #      Gtud = S.gettud(itx,irx)
    #     """
    #     if (itx < self.tx.N) & (irx < self.rx.N):
    #         Gtud = GrRayTud()
    #         Gtud.load(self.filetud[itx][irx], self.sl)
    #         return Gtud
    #     else:
    #         print "gettud warning : wrong tx or rx index"
    # #def gettud(self,k,l):
    # #def getfield(self,k,l):

    # def launching(self, itx=1,verbose=False):
    #     """ start the launching program and get the results files

    #     Parameters
    #     ----------
    #     itx : int
    #         transmiter index

    #     """

    #     filestr = os.path.splitext(self.filestr)[0] + '.str'

    #     if not os.path.exists(pyu.getlong(filestr,pstruc['DIRSTRUC'])):
    #         chaine = 'newstruc -str2 ' + filestr +'2 ' + filestr + ' -conf ' + basename +'/'+self.fileconf
    #         os.system(chaine)


    #     chaine = "launching -str  " + self.filestr + \
    #         " -slab " + self.fileslab + \
    #         " -palch " + self.filepalch + \
    #         " -spa " + self.tx.filespa + \
    #         " -conf " + basename + '/' + self.fileconf

    #     if verbose:
    #         print chaine

    #     self.claunching.append(chaine)
    #     os.system(chaine)
    #     aux = os.popen(chaine, "r")
    #     recup = aux.read()
    #     aux.close()
    #     if verbose:
    #         print recup

    #     self.recup = recup
    #     aux = recup.splitlines()
    #     len_aux = recup.count("\n")

    #     filelch = []
    #     for i in range(len_aux):
    #         if aux[i].find("filelchout") != -1:
    #             aux[i] = aux[i].replace('filelchout : ', '')
    #             fshort = pyu.getshort(aux[i])
    #             filelch.append(fshort)
    #             self.dlch[itx] = fshort

    #     self.filelch = filelch
    #     self.progress = 1
    #     self.nTx = len(self.filelch)
    #     self.dtra[itx] = {}
    #     self.dtud[itx] = {}
    #     self.dtang[itx] = {}
    #     self.drang[itx] = {}
    #     self.dtauk[itx] = {}
    #     self.dfield[itx] = {}
    #     self.dcir[itx] = {}
    #     # create a configuration file for Tx itx
    #     self.output[itx] = ConfigParser.ConfigParser()
    #     self.output[itx].add_section("launch")
    #     self.output[itx].add_section("trace")
    #     self.output[itx].add_section("tud")
    #     self.output[itx].add_section("rang")
    #     self.output[itx].add_section("tang")
    #     self.output[itx].add_section("field")
    #     self.output[itx].add_section("tauk")
    #     self.output[itx].set("launch", str(itx), self.dlch[itx])
    #     #
    #     # append output filename in section output
    #     #
    #     _outfilename = self.filesimul.replace('.ini', '') + str(itx) + ".ini"
    #     self.dout[itx] = _outfilename
    #     outfilename = pyu.getlong(_outfilename, pstruc['DIRLCH'])
    #     self.config.set("output", str(itx), self.dout[itx])

    #     fd = open(outfilename, "w")
    #     self.output[itx].write(fd)
    #     fd.close()

    #     filesimul = pyu.getlong(self.filesimul,'ini')
    #     fd = open(filesimul, "w")
    #     self.config.write(fd)
    #     fd.close()

    # def tracing(self, itx, irx,verbose=False):
    #     """ exec tracing

    #     Parameters
    #     ----------
    #     itx  : tx index
    #     irx  : rx index

    #     Notes
    #     -----
    #     This function should not be used for more than one rx.
    #     .. todo extend properly this function in order to handle properly
    #         multi-nodes in irx. The problem is to keep the association
    #         between the index number of the rx in the ini file and the
    #         rx in dtra.

    #     """
    #     #
    #     # Verify
    #     #
    #     if (self.progress >= 1):
    #         chaine = "tracing -lch " + self.dlch[itx] + \
    #             " -patra " + self.filepatra + \
    #             "  -spa " + self.rx.filespa + \
    #             " -conf " + basename + '/' + self.fileconf
    #         if verbose:
    #             print chaine
    #         self.ctracing.append(chaine)
    #         aux = os.popen(chaine, "r")
    #         recup = aux.read()
    #         #if verbose:
    #             #print recup
    #         aux.close()
    #         self.recup = recup
    #         aux = recup.splitlines()
    #         len_aux = recup.count("\n")
    #         #
    #         # set list of .tra files empty
    #         #filetra = []
    #         #for i in range(len_aux):
    #         #    if aux[i].find("filetraout") != -1:
    #         #        aux[i] = aux[i].replace('filetraout : ', '')
    #         #        filetra.append(pyu.getshort(aux[i]))
    #         #
    #         # Warning : this is a bad fix
    #         #
    #         #for filename in filetra:
    #         for i in range(len_aux):
    #             if aux[i].find("filetraout") != -1:
    #                 aux[i] = aux[i].replace('filetraout : ', '')
    #                 self.dtra[itx][irx] = pyu.getshort(aux[i])

    #         #if verbose:
    #             #print filetra
    #         #self.filetra.insert(ntx,filetra)
    #         self.progress = 2
    #     else:
    #         print "No launching available"

    #     _outfilename = self.config.get('output', str(itx))
    #     outfilename = pyu.getlong(_outfilename, pstruc['DIRLCH'])
    #     if irx in self.dtra[itx].keys():
    #         self.output[itx].set("trace", str(irx), self.dtra[itx][irx])
    #         fd = open(outfilename, "w")
    #         self.output[itx].write(fd)
    #         fd.close()

    # def tratotud(self, itx, irx,verbose=False):
    #     """ convert tracing in .tud

    #     Parameters
    #     ----------

    #      itx : integer
    #         transmitter index
    #      irx : integer
    #         receiver index

    #     """
    #     #
    #     # .. todo:: take value from the simulation class
    #     #
    #     nrmin = self.config.get("tud", "nrmax")
    #     num = self.config.get("tud", "num")
    #     purc = self.config.get("tud", "purc")
    #     chaine = "tratotud -tra " + self.dtra[itx][irx] + \
    #         " -min " + nrmin + \
    #         " -purc " + purc + \
    #         " -num " + num +  \
    #         " -conf " + basename + '/' + self.fileconf
    #     print 'DEBUG '+chaine
    #     self.ctratotud.append(chaine)
    #     if verbose:
    #         print chaine
    #     aux = os.popen(chaine, "r")
    #     os.system('echo $?')
    #     recup = aux.read()
    #     aux.close()
    #     aux = recup.splitlines()
    #     if verbose:
    #         print aux
    #     len_aux = recup.count("\n")
    #     for i in range(len_aux):
    #         if aux[i].find("filetudout") != -1:
    #             aux[i] = aux[i].replace('filetudout : ', '')
    #             filename = pyu.getshort(aux[i])
    #             # fix
    #             #filename = rename(filename,itx,irx,'output')
    #             if verbose:
    #                 print filename
    #             self.dtud[itx][irx] = filename
    #         elif aux[i].find("filetangout") != -1:
    #             aux[i] = aux[i].replace('filetangout : ', '')
    #             filename = pyu.getshort(aux[i])
    #             # fix
    #             #filename = rename(filename,itx,irx,'output')
    #             if verbose:
    #                 print filename
    #             self.dtang[itx][irx] = filename
    #         elif aux[i].find("filerangout") != -1:
    #             aux[i] = aux[i].replace('filerangout : ', '')
    #             filename = pyu.getshort(aux[i])
    #             # fix
    #             #filename = rename(filename,itx,irx,'output')
    #             if verbose:
    #                 print filename
    #             self.drang[itx][irx] = filename

    #     _outfilename = self.config.get('output', str(itx))
    #     outfilename = pyu.getlong(_outfilename, pstruc['DIRLCH'])
    #     if irx in self.dtud[itx].keys():
    #         self.output[itx].set("tud", str(irx), self.dtud[itx][irx])
    #     if irx in self.dtang[itx].keys():
    #         self.output[itx].set("tang", str(irx), self.dtang[itx][irx])
    #     if irx in self.drang[itx].keys():
    #         self.output[itx].set("rang", str(irx), self.drang[itx][irx])
    #     fd = open(outfilename, "w")
    #     try:
    #         self.output[itx].write(fd)
    #     except:
    #         raise NameError('error writing output ini file')
    #     fd.close()


    # def field(self, itx, irx,verbose=False):
    #     """ field calculation for Tx Rx using evalfield command

    #     Parameters
    #     ----------

    #     ntx : integer
    #              launching index
    #     nrx : integer
    #              tracing index
    #     verbose : Boolean
          

    #     """
    #     chaine = "evalfield -tud " + self.dtud[itx][irx] + \
    #              " -slab " + self.fileslab + \
    #              " -mat " + self.filemat + \
    #              " -freq " + self.filefreq + \
    #              " -conf " + basename + '/' + self.fileconf

    #     self.cfield.append(chaine)
    #     if verbose:
    #         print chaine
    #     os.system(chaine)
    #     aux = os.popen(chaine, "r")
    #     recup = aux.read()
    #     aux.close()
    #     aux = recup.splitlines()
    #     len_aux = recup.count("\n")
    #     for i in range(len_aux):
    #         if aux[i].find("filefieldout") != -1:
    #             aux[i] = aux[i].replace('filefieldout : ', '')
    #             filename = pyu.getshort(aux[i]).replace(' ', '')
    #             # fix
    #             # filename = rename(filename,itx,irx,'output')
    #             self.dfield[itx][irx] = filename
    #         elif aux[i].find("filetaukout") != -1:
    #             aux[i] = aux[i].replace('filetaukout : ', '')
    #             filename = pyu.getshort(aux[i]).replace(' ', '')
    #             # fix
    #             # filename = rename(filename,itx,irx,'output')
    #             self.dtauk[itx][irx] = filename
    #     _outfilename = self.config.get('output', str(itx))
    #     outfilename = pyu.getlong(_outfilename, pstruc['DIRLCH'])
    #     if irx in self.dtauk[itx].keys():
    #         self.output[itx].set("tauk", str(irx), self.dtauk[itx][irx])
    #     if irx in self.dfield[itx].keys():
    #         self.output[itx].set("field", str(irx), self.dfield[itx][irx])
    #     fd = open(outfilename, "w")
    #     try:
    #         self.output[itx].write(fd)
    #     except:
    #         raise NameError('error writing output ini file')
    #     fd.close()


    # def run2(self, link, cirforce=True,verbose=False,cutoff=4):
    #     """ run the simulation for 1 tx and a set of rx

    #         Parameters
    #         ----------

    #         itx      : tx index
    #         srx      : list of rx index
    #         cirforce : boolean

    #         Warnings
    #         --------

    #         index point start with 1

    #         Example
    #         -------

    #         >>> from pylayers.simul.simulem import *
    #         >>> itx = 1
    #         >>> srx = [1,2,3]
    #         >>> S   = Simul()
    #         >>> S.load('where2.ini')
    #         >>> out = S.run(itx,srx)


    #     """

    #     prefix = self.filesimul.replace('.ini', '')
    #     #
    #     #
    #     #
    #     lsig = []
    #     for k,il  in enumerate(link):
    #         tx = self.tx.points[il[0]]
    #         rx = self.rx.points[il[1]]

    #         ctx = S.L.pt2cy(tx)
    #         crx = S.L.pt2cy(rx)

    #         _filecir = prefix +'-cir-'+str(k)+'-'+str(link)+'-'+str((ctx,crx))
    #         D = {}
    #         D['Tx'] = tx
    #         D['Rx'] = rx

    #         if (ctx,crx) not in lsig:
    #             Si  = signature.Signatures(S.L,ctx,crx)
    #             #
    #             # Change the run number depending on
    #             # the algorithm used for signature determination
    #             #
    #             Si.run4(cutoff=cutoff)
    #             # keep track and save signature
    #             _filesir = prefix + '-sig-'+str((ctx,crx))
    #             fd = open(filesig,'w')
    #             pickle.dump(Si,filesig)
    #             fd.close()
    #             lsig.appeng((ctx,crx))
    #             Si.dump(S.L,(ctx,crx))



    #         r2d = Si.rays(tx,rx)
    #         r2d.show(S.L)

    #         r3d = r2d.to3D()
    #         r3d.locbas(S.L)
    #         r3d.fillinter(S.L)

    #         Ct  = r3d.eval(S.freq)
    #         sco = Ct.prop2tran(a='theta',b='phi')
    #         sca = Ct.prop2tran(a=S.tx.A,b=S.rx.A)

    #         ciro = sco.applywavB(self.wav.sfg)
    #         cira = sca.applywavB(self.wav.sfg)

    #         D['to'] = ciro.x
    #         D['ciro'] = ciro.y
    #         D['t'] = cira.x
    #         D['cir'] = cira.y

    #         filename = pyu.getlong(_filename, cirdir)
    #         spio.savemat(filename, D)

    # def run(self, itx, srx=[], cirforce=True,verbose=False):
    #     """ run the simulation for 1 tx and a set of rx

    #         Parameters
    #         ----------

    #         itx      : tx index
    #         srx      : list of rx index
    #         cirforce : boolean

    #         Warnings
    #         --------

    #         index point start with 1

    #         Example
    #         -------

    #         >>> from pylayers.simul.simulem import *
    #         >>> itx = 1
    #         >>> srx = [1,2,3]
    #         >>> S   = Simul()
    #         >>> S.load('where2.ini')
    #         >>> out = S.run(itx,srx)


    #     """

    #     self.checkh5()
    #     #t0 = time.clock()
    #     if type(srx) == int:
    #         srx = [srx]

    #     if srx == []:
    #         srx = self.rx.points.keys()

    #     if itx not in self.dlch.keys():
    #         # launching itx does not exist
    #         self.tx.filespa = 'tx' + str(itx) + '.spa'
    #         point = self.tx.points[itx]
    #         spafile(self.tx.filespa, point, pstruc['DIRLCH'])
    #         if verbose:
    #             print "---------------"
    #             print "Start Launching Tx : " + str(itx)
    #             print "---------------"
    #         self.launching(itx)

    #     #
    #     # Loop over a set of rx
    #     #
    #     #pdb.set_trace()
    #     for irx in srx:
    #         tracingabort = False
    #         if irx not in self.dtra[itx].keys():
    #             self.rx.filespa = 'rx' + str(irx) + '.spa'
    #             point = self.rx.points[irx]
    #             spafile(self.rx.filespa, point, pstruc['DIRTRA'])
    #             if verbose:
    #                 print "--------------------"
    #                 print "Start tracing  Rx : " + str(irx)
    #                 print "--------------------"
    #             tracingabort = self.tracing(itx, irx,verbose)

    #         if not tracingabort:
    #             if irx not in self.dtud[itx].keys():
    #                 if verbose:
    #                     print "---------------"
    #                     print "Start tratotud ", irx
    #                     print "---------------"

    #                 self.tratotud(itx, irx,verbose)
    #             if irx not in self.dfield[itx].keys():
    #                 if verbose:
    #                     print "---------------"
    #                     print "Start field  ", irx
    #                     print "---------------"
    #                 self.field(itx, irx,verbose)

    #             if ((irx not in self.dcir[itx].keys()) | cirforce):
    #                 if verbose:
    #                     print "---------------"
    #                     print "Start cir      ", irx
    #                     print "---------------"
    #                 if "waveform" in self.config.sections():
    #                     par = self.config.items("waveform")
    #                     self.wparam = {}
    #                     for k in range(len(par)):
    #                         key = par[k][0]
    #                         val = par[k][1]
    #                         if key == "band":
    #                             self.wparam[key] = float(val)
    #                         if key == "fcGHz":
    #                             self.wparam[key] = float(val)
    #                         if key == "feGHz":
    #                             self.wparam[key] = float(val)
    #                         if key == "threshdB":
    #                             self.wparam[key] = float(val)
    #                         if key == "twns":
    #                             self.wparam[key] = float(val)
    #                         if key == "typ":
    #                             self.wparam[key] = val

    #                     self.wav = wvf.Waveform(**self.wparam)
    #                     alpha = np.sqrt(1. / 30.0)
    #                     print "run debug ",itx,irx
    #                     self.cir([itx], [irx],
    #                              store_level=16 + 8 + 4 + 2 + 1, alpha=alpha)
    #                 else:
    #                     raise("Error no waveform in the config file ")
    #                     return(False)

    #     return(True)
    # def gt(self, itx, irx):
    #     """ gtud
    #     """
    #     gt = GrRayTud()
    #     filetud = self.dtud[itx][irx]
    #     filetang = self.dtang[itx][irx]
    #     filerang = self.drang[itx][irx]
    #     gt.load(filetud, filetang, filerang, self.sl)
    #     return(gt)

    # def delay(self, itx, irx):
    #     """
    #         calculate LOS link delay

    #         Parameters
    #         ----------

    #         itx
    #         irx

    #     """
    #     tx = self.tx.points[itx]
    #     rx = self.rx.points[irx]
    #     df = tx - rx
    #     dist = np.sqrt(np.dot(df, df))
    #     return(dist / 0.3)

    # def gr(self, itx, irx):
    #     """
    #        return a  cluster or fays from link itx-irx
    #     """
    #     gr = GrRay3D()
    #     gr.load(self.dtra[itx][irx], self.L)
    #     return(gr)

    # def VC(self, itx, irx):
    #     """
    #         return Vect Channel for link itx irx
    #     """

    #     VCl = channelc.VectChannel(self, itx, irx, False)
    #     return(VCl)

    # def cir(self, itx, irx, store_level=0, alpha=1.0, ext='', rep=pstruc['DIRCIR'],format='a'):
    #     """
    #     Calculate a set of channel impulse responses

    #     Parameters
    #     ----------
    #         itx : transmitter index iterable set
    #         irx : receiver index iterable set
    #         wav : applied waveform (st,sf,sfg)
    #             this waveform includes the gamma factor
    #         store_level : binary mask
    #                          bit  0 : save CVC
    #                          bit  1 : save CVCO
    #                          bit  2 : save CSCA
    #                          bit  3 : save CIRo
    #                          bit  4 : save CIRa
    #         alpha : normalization factor
    #         ext :
    #         rep :
    #         format : string
    #             a : with antenna
    #             o : omnidirectionnal
    #      Notes
    #      -----
    #      A factor ::math`\\sqrt{1}{30}` is necessary when applying the antenna

    #      Examples
    #      --------

    #     """
    #     if type(itx) == int:
    #         itx = [itx]

    #     if type(irx) == int:
    #         irx = [irx]

    #     if (store_level) & 1 == 1:
    #         self.CVC = []
    #     if (store_level) & 2 == 2:
    #         self.CSCO = []
    #     if (store_level) & 4 == 4:
    #         self.CSCA = []
    #     if (store_level) & 8 == 8:
    #         self.CIRo = []
    #     if (store_level) & 16 == 16:
    #         self.CIRa = []
    #     racine = self.filesimul.replace('.ini', '') + 'cir-'
    #     for l in itx:
    #         # create cir entry in outputTx if required
    #         _outfilename = self.config.get('output', str(l))
    #         outfilename = pyu.getlong(_outfilename, pstruc['DIRLCH'])
    #         if "cir" not in  self.output[l].sections():
    #             self.output[l].add_section("cir")

    #         CVC = []
    #         CSCO = []
    #         CSCA = []
    #         CIRo = []
    #         CIRa = []
    #         for k in irx:
    #             D = {}
    #             D['Tx'] = self.tx.points[l]
    #             if ext == '':
    #                 if (self.tx.name == '') and (self.rx.name == ''):
    #                     txrx = 'tx' + str('%0.3d' % l) + '-rx' + str('%0.3d' % k)
    #                 else :
    #                     txrx = self.tx.name +'-' + self.rx.name + str('-p%0.3d' % k)
    #             else:
    #                 txrx = ext

    #             _filename = racine + txrx
    #             self.dcir[l][k] = _filename
    #             if self.tx.name == '':
    #                 rep = rep + '/Tx' + str('%0.3d' % l)
    #             else :
    #                 rep = rep +'/' + self.tx.name
    #             if not os.path.isdir(basename+'/'+rep):
    #                 try:
    #                     os.mkdir(basename+'/'+rep)
    #                 except:
    #                     raise NameError(basename+'/'+rep)


    #             filename = pyu.getlong(_filename, rep)
    #             VCl = channelc.VectChannel(self, l, k, False)
    #             CVC.append(VCl)
    #             if not VCl.fail:
    #                 #SCO = VCl.vec2scal()
    #                 SCO = VCl.prop2tran(a='theta',b='theta')
    #                 CSCO.append(SCO)
    #                 #SCA = VCl.vec2scalA(self.tx.A, self.rx.A, alpha=alpha)
    #                 #
    #                 #  Apply the apha factor on waveform
    #                 #
    #                 SCA = VCl.prop2tran(a=self.tx.A,b=self.rx.A)
    #                 CSCA.append(SCA)
    #                 ciro = SCO.applywavB(self.wav.sfg)
    #                 CIRo.append(ciro)
    #                 cira = SCA.applywavB(self.wav.sfg)
    #                 CIRa.append(cira)
    #                 D['Rx' + str(k)] = self.rx.points[int(k)]
    #                 if 'o' in format:
    #                     D['to'] = ciro.x
    #                     D['ciro'] = ciro.y
    #                 if 'a' in format:
    #                     D['t'] = cira.x
    #                     D['cir'] = cira.y
    #                 spio.savemat(filename, D)
    #                 self.output[l].set("cir", str(k), self.dcir[l][k])
    #                 fd = open(outfilename, "w")
    #                 self.output[l].write(fd)
    #                 fd.close()
    #             else:
    #                 CSCO.append([])
    #                 CSCA.append([])
    #                 CIRo.append([])
    #                 CIRa.append([])

    #         if (store_level) & 1 == 1:
    #             self.CVC.append(CVC)
    #             self.CSCO.append(CSCO)
    #         if (store_level) & 4 == 4:
    #             self.CSCA.append(CSCA)
    #         if (store_level) & 8 == 8:
    #             self.CIRo.append(CIRo)
    #         if (store_level) & 16 == 16:
    #             self.CIRa.append(CIRa)


if (__name__ == "__main__"):
    #plt.ion()
    doctest.testmod()
