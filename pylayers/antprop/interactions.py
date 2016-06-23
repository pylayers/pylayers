#!/usr/bin/python
# -*- coding: utf8 -*-
"""

.. currentmodule:: pylayers.antprop.interactions

Inters Class
============

.. autosummary::
    :toctree: generated/

    Inter.__init__
    Inter.__repr__
    Inter.create_dusl
    Inter.sinsout
    Inter.stack

Interactions Class
==================


.. autosummary::
    :toctree: generated/

    Interactions.__init__
    Interactions.add
    Interactions.addi
    Interactions.eval

IntB Class
===========

.. autosummary::
    :toctree: generated/

    IntB.__init__
    IntB.__repr__
    IntB.eval

IntR Class
===========

.. autosummary::
    :toctree: generated/

    IntR.__init__
    IntR.__repr__
    IntR.eval

IntT Class
===========

.. autosummary::
    :toctree: generated/

    IntT.__init__
    IntT.__repr__
    IntT.eval

IntD Class
===========

.. autosummary::
    :toctree: generated/

    IntD.__init__
    IntD.__repr__
    IntD.eval
"""
import pdb
import os
import pdb
import glob
import doctest
import ConfigParser
import networkx as nx
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import struct as stru
import pylayers.util.geomutil as geu
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.antprop.slab import *
from pylayers.antprop.diffRT import *



class Inter(PyLayers):
    """ Interactions

    Meta class of interactions ( Interactions, IntB/IntL/IntT/intR/intD)

    Attributes
    ----------

        typ : int
            type of interaction
            1 : D
            2 : R
            3 : T
            0 : Tx or Rx
           -1 : B

        data: np.array
            data for the interaction
        idx:
            idx number of the interaction between 0 and (ray number * inter number)
        fGHz : np.array
            frequency range in GHz
        nf : int
            number of step into freq range
        olf : np.array
            np.ones((nf)) used for broadcasting


    """

    def __init__(self, 
            typ=0, 
            data=np.array(()), 
            idx=[],
            _filemat='matDB.ini',
            _fileslab='slabDB.ini',
            slab={}):
        """ Inter object constructor

        Parameters
        ----------

        typ : int
        data  : ndarray
        idx : list
        _filemat : string
        _fileslab : string
        slab : SlabDB

        """

        self.typ = typ
        self.data = data
        self.idx = idx
        if slab=={}:
            self.slab = SlabDB(filemat=_filemat, fileslab=_fileslab)
        else:
            self.slab = slab

        self.idx = []
        if idx != []:
            self.idx.append(idx)

        self.E = np.eye(2)

    def __repr__(self):
        if self.evaluated:
            s = self.T.__repr__()
            s = s + '\n' + self.R.__repr__()
            s = s + '\n'+ self.D.__repr__()
            return s
        else:
            return 'I not yet evaluated'

    def create_dusl(self,a):
        """ create dictionnary of used slab.

        Parameters
        ----------

        a : np.array of strings which contains ordered interactions
            ordered as in self.idx/self.data

        """
        for s in self.dusl:
            self.dusl[s]=np.where(a==s)[0]



    def sinsout(self):
        """ calculate sin sout of the interaction

        Notes
        -----

        typ

            1 : Diffraction
            2 : Reflection
            3 : Transmission

            si : self.data[:,1]
            so : self.data[:,2]

        typ = 0
            LOS

        typ = -1
           Basis

        """

        if self.typ in [2, 3]: #reflection & transmission
            self.si0 = self.data[:, 1]
            self.sout = self.data[:, 2]
        elif self.typ == 1: # diffraction
            self.si0 = self.data[:, 4]
            self.sout = self.data[:, 5]
        elif self.typ == 0: # loss
            self.sout = self.data[0]
        elif self.typ == -1: # B
            self.sout = np.zeros((len(self.data[:, 0])))

    def stack(self, data=np.array(()), idx=0, isdata=True):
        """ stack data and the associated idx

        Parameters
        ----------

        data : np.array()
            data to stack
        idx :
            index to stack
        isdata: bool
            False if you just want to stack idx (only used for intE class )

        Examples
        --------

        >>> from pylayers.antprop.rays import *
        >>> import numpy as np
        >>> I=Inter()
        >>> data = np.array(([3,4,5]))
        >>> idx = 0
        >>> I.stack(data,idx)
        >>> I.data
        array([3, 4, 5])
        >>> I.idx
        [0]
        >>> data = np.array(([3,4,5],[7,8,9]))
        >>> idx = [1,2]
        >>> I.stack(data,idx)
        >>> I.data
        array([[3, 4, 5],
               [3, 4, 5],
               [7, 8, 9]])
        >>> I.idx
        [0, 1, 2]
        """

        if isinstance(idx, int):
            try:
                if isdata:
                    self.data = np.vstack((self.data, data))
                self.idx.append(idx)
            except:
                if self.idx == []:
                    if isdata:
                        self.data = data
                    self.idx = [idx]
                else:
                    raise NameError('Issue in Inter.stack')

        elif isinstance(idx, list) or isinstance(idx, np.ndarray):
            try:
                self.data = np.vstack((self.data,data))
            except:
                self.data=data
            self.idx.extend(idx)

#            for ii, idx in enumerate(idx):
#                if isdata:
#                    try:
#                        self.data = np.vstack((self.data, data[ii]))
#                    except:
#                        self.data = data[ii]
#                self.idx.append(idx)


class Interactions(Inter,dict):
    """ Interaction parameters

        gather all type of interactions (IntB/L/R/T)

        Methods
        -------

        add(self,li): add a list of basis interactions
        addi(self,i): add a single interaction
        eval(self) : evaluate all the interactions added thanks to self.add or self.addi
                     and create the self.I which gather all thoses interactions

        5 following types of interactions

        B : local basis transformation matrix (unitary)
        L : LOS case
        R : Reflection
        T : Transmission
        D : Diffraction


    """

    def __init__(self,slab={}):
        """ object constructor
        """
        Inter.__init__(self,slab=slab)
        self['B'] = []
        self['L'] = []
        self[''] = []
        self['T'] = []
        self['D'] = []
        self.evaluated = False
        self.nimax = 0

    def add(self, li):
        """ add a list of interactions

        Parameters
        ----------

        li : list
            list of interactions

        """
        # determine the total number of interactions
        for i in li:
            if i.idx != []:
                self.nimax = max(self.nimax,max((i.idx)))+1

        for i in li:
            self.addi(i)

    def addi(self, i):
        """ add interactions into Interactions class

        Parameters
        ----------

        i : Inter object

        """

        if not isinstance(self.typ, np.ndarray):
            self.typ = np.zeros((self.nimax), dtype=str)
        if i.typ == -1:
            self.B = i
            self['B'] = i.idx
            self.typ[i.idx] = 'B'
        if i.typ == 0:
            self.L = i
            self['L'] = i.idx
            self.typ[i.idx] = 'L'
        if i.typ == 1:
            self.D = i
            self['D'] = i.idx
            self.typ[i.idx] = 'D'
        if i.typ == 2:
            self.R = i
            self['R'] = i.idx
            self.typ[i.idx] = 'R'
        if i.typ == 3:
            self.T = i
            self['T'] = i.idx
            self.typ[i.idx] = 'T'

    def eval(self,fGHz=np.array([2.4])):
        """ evaluate all the interactions

        Parameters
        ----------

        fGHz : np.array()

        Notes
        -----

        self.I : np.shape(self.I) = (self.nf,self.nimax,2,2)
            with self.nf :  number of frequences
            self.nimax : the total number of interactions ( of all rays)
        self.sout :
            distance from one interaction to the next one
        self.si0 :
            distance from the previous interaction to the one
        self.alpha :
            alpha as described in JFL Thesis
        self.gamma :
            !! gamma**2 !!! (squared included) as described

        """

        # Initialize the global I matrix which gathers all interactions
        # into a single np.array

        # f x i x 2 x 2

        self.fGHz = fGHz
        self.nf = len(fGHz)

        self.I = np.zeros((self.nf, self.nimax, 2, 2), dtype=complex)
        self.sout = np.zeros((self.nimax))
        self.si0 = np.zeros((self.nimax))
        self.alpha = np.ones((self.nimax), dtype=complex)
        self.gamma = np.ones((self.nimax), dtype=complex)

        # evaluate B and fill I 
        #OUT DATED , B MDA are stored outside of I
        # try:
        #     self.I[:, self.B.idx, :, :] = self.B.eval(fGHz=fGHz)
        #     self.sout[self.B.idx] = self.B.sout
        #     self.si0[self.B.idx] = self.B.si0

        # except:
        #     print 'Warning : No B interaction Evaluated'

        # evaluate L and fill I
        # OUT DATED , Los interaction is managed oustside of I
        # try:
        #     self.I[:, self.L.idx, :, :] = self.L.eval(fGHz=fGHz)
        #     self.sout[self.L.idx] = self.L.sout
        #     self.si0[self.L.idx] = self.L.si0

        # except:
        #     print 'Warning : No L interaction Evaluated'

        # evaluate R and fill I
        pdb.set_trace()
        if len(self.R.data)!=0:
            #try:
            self.I[:, self.R.idx, :, :] = self.R.eval(fGHz=fGHz)
            self.sout[self.R.idx] = self.R.sout
            self.si0[self.R.idx] = self.R.si0
            self.alpha[self.R.idx] = self.R.alpha
            self.gamma[self.R.idx] = self.R.gamma
            #except:
            #    print Warning('Warning Interaction.eval: No R interaction Evaluated,\ whereas Reflection rays found')
            #    pdb.set_trace()
        # evaluate T and fill I
        if len(self.T.data)!=0:
            #try:
            self.I[:, self.T.idx, :, :] = self.T.eval(fGHz=fGHz)
            self.sout[self.T.idx] = self.T.sout
            self.si0[self.T.idx] = self.T.si0
            self.alpha[self.T.idx] = self.T.alpha
            self.gamma[self.T.idx] = self.T.gamma

            #except:
            #     print Warning('Warning Interaction.eval: No T interaction Evaluated,\ whereas Transmission rays found')
            #    pdb.set_trace()
        # evaluate D and fill I

        if len(self.D.data)!=0:
            #try:
            self.I[:, self.D.idx, :, :] = self.D.eval(fGHz=fGHz)
            self.sout[self.D.idx] = self.D.sout
            self.si0[self.D.idx] = self.D.si0

                # self.alpha[self.D.idx] = self.D.alpha
                # self.gamma[self.D.idx] = self.D.gamma
            #except:
            #    print Warning('Warning Interaction.eval: No D interaction Evaluated,\ whereas Diffraction rays found')
            #    pdb.set_trace()

        self.evaluated = True



class IntB(Inter):
    """ Local Basis interaction class

        Basis interactions

        Attributes
        ----------

        data : np.array:
            WARNING np.shape(data) = (ninter x 4)
            the input matrix 2x2 is reshaped as 1x 4
        idx : list
            index of the corresponding ray and interaction

        Methods
        -------

        eval : evaluation of B interaction

        Notes
        -----

        The interaction object is np.array with shape (nf,ninter 2, 2)

    """
    def __init__(self, data=np.array(()), idx=[],slab={}):
        Inter.__init__(self, data=data, idx=idx, typ=-1,slab=slab)

    def __repr__(self):
        s = 'number of B basis :' + str(np.shape(self.data)[0])
        return s    

    def eval(self,fGHz=np.array([2.4])):
        """ evaluation of B interactions

        Parameters
        ----------

        fGHz : np.array()nn

            frequency range

        Returns
        -------

        self.data

        Examples
        --------

        >>> from pylayers.antprop.rays import *
        >>> M = np.eye(2).reshape(4)
        >>> B = IntB(M,0)
        >>> B.data
        array([ 1.,  0.,  0.,  1.])
        >>> B.stack(M,1)
        >>> B.data
        array([[ 1.,  0.,  0.,  1.],
               [ 1.,  0.,  0.,  1.]])
        >>> eB=B.eval()
        >>> nf = B.nf
        >>> ninter = len(B.idx)
        >>> np.shape(eB)
        (1, 2, 2, 2)

        """

        self.fGHz=fGHz
        self.nf=len(fGHz)


        self.sinsout()
        if len(self.data) != 0:
            lidx = len(self.idx)
            data = self.data.reshape(lidx, 2, 2)
            #return(self.olf[:, np.newaxis, np.newaxis, np.newaxis]*data[np.newaxis, :, :, :])
            return(np.ones((len(fGHz),1,1,1))*data[None, :, :, :])
        else:
            print 'no B interaction to evaluate'
            return(self.data[:, None, None, None])

#####Interaction Loss in not used for speed purpose
#the loss interaction is computed and added after the global computation
# class IntL(Inter):
#     """ Loss interaction

#     (nf ,ninter, 2, 2)

#     Attributes
#     ----------

#         data : np.array((ninter x [dist]))
#         idx : list
#             index of the corresponding ray and interaction

#     """
#     def __init__(self, data=np.array(()), idx=[]):
#         Inter.__init__(self, data=data, idx=idx, typ=0)

#     def __repr__(self):
#         s = Inter.__repr__(self)
#         return(s)

#     def eval(self,fGHz=np.array([2.4])):
#         """ evaluation of B interactions

#         Examples
#         --------

#         >>> from pylayers.antprop.rays import *
#         >>> d = np.array(([3]))
#         >>> L = IntL(d,0)
#         >>> L.data
#         array([3])
#         >>> L.stack(d+3,1)
#         >>> L.data
#         array([[3],
#                [6]])
#         >>> eL=L.eval()
#         >>> ninter = len(L.idx)
#         >>> nf = L.nf
#         >>> np.shape(eL)
#         (1, 1, 2, 2)

#         """


#         self.fGHz=fGHz
#         self.nf=len(fGHz)

#         self.sinsout()

#         if len(self.data != 0):
#             try:
#                 np.shape(self.data)[1]
#             except:
#                 # it means that self.data is not a matrix but a number
#                 self.data = self.data.reshape(1, 1)
#             ## Friis formula
# #            self.data=self.data[:,0]
# #            dis = (0.3 / (4*np.pi*self.data[np.newaxis,:]*self.f[:,np.newaxis]))
# #            return(dis[:,:,np.newaxis,np.newaxis]*self.E[np.newaxis,np.newaxis,:,:])
# #            dis = self.data
#             #return(self.olf[:, np.newaxis, np.newaxis, np.newaxis]*self.E[np.newaxis, np.newaxis, :, :])
#             return(np.ones((len(fGHz),1,1,1))*self.E[None, None, :, :])
#         else:
#             print 'no L interaction to evaluate'
#             return(self.data[:, None, None, None])


class IntR(Inter):
    """ Reflexion interaction class


    """
    def __init__(self, data=np.array(()), idx=[],slab={}):
#        self.theta = data[0]
#        self.si = data[1]
#        self.sr = data[2]
        Inter.__init__(self, data=data, idx=idx, typ=2,slab=slab)
        ## index for used slab
        self.uslidx = 0
        # dictionnary of used slab key = slab value = index of self.idx
        # WARNING The index of this dictionnary referes to the idex of self.idx
        # not to the global indx
        self.dusl = {}
#       self.dusl=dict.fromkeys(self.slab,np.array((),dtype=int))
        self.alpha = [1]
        self.gamma = [1]

    def __repr__(self):
        s = 'number of R interaction :' + str(np.shape(self.data)[0])
        return s    

    def eval(self,fGHz=np.array([2.4])):
        """ evaluation of reflexion interactions

        Parameters
        ----------

        fGHz : np.array (,Nf)


        Returns
        -------

        self.A  : evaluated interaction


        Examples
        --------

        >>> from pylayers.antprop.rays import *

        >>> # generate input data
        >>> theta1 = 0.1
        >>> theta2 = 0.4
        >>> si01 = 4
        >>> si02 = 0.6
        >>> sir1 = 3.15
        >>> sir2 = 3.4
        >>> data1=np.array((theta1,si01,sir1))
        >>> data2=np.array((theta2,si02,sir2))

        >>> # store input data to Instance
        >>> R = IntR(data1,idx=0)
        >>> R.data
        array([ 0.1 ,  4.  ,  3.15])
        >>> R.stack(data2,idx=1)
        >>> R.uslidx=1
        >>> R.dusl['WOOD']=[0,1]

        >>> # evaluation parameters (normally read from config.ini)
        >>> R.f = np.array([  2.,  11.])
        >>> R.nf = len(R.f)
        >>> R.olf = np.ones((R.nf))

        >>> # evaluation
        >>> eR=R.eval()

        >>> # examples
        >>> ninter = len(R.idx)
        >>> np.shape(eR)
        (1, 2, 2, 2)
        >>> eR[0,0,0]
        array([-0.0413822-0.12150975j,  0.0000000+0.j        ])

        Notes
        -----

        data = np.array((ninter x [theta,si,st]))

        """
        self.sinsout()

        self.fGHz=fGHz
        self.nf=len(fGHz)

        # A : f ri 2 2

        self.A = np.zeros((self.nf, len(self.idx), 2, 2), dtype=complex)

        if np.shape(self.data)[0]!=len(self.idx):
            self.data=self.data.T


        if len(self.data) != 0:
            mapp = []
            # loop on all type of materials used for reflexion
            for m in self.dusl.keys():
                # used theta of the given slab
                ut = self.data[self.dusl[m], 0]
                if not ut.size == 0:
                    # find the index of angles which satisfied the data
                    if m not in self.slab:
                        m = m.lower()
                    self.slab[m].ev(fGHz=fGHz, theta=ut, RT='R')
                    try:
                        R = np.concatenate((R, self.slab[m].R), axis=1)
                        mapp.extend(self.dusl[m])
                    except:
                        R = self.slab[m].R
                        mapp.extend(self.dusl[m])

            # replace in correct order the reflexion coeff
            self.A[:, np.array((mapp)), :, :] = R
            self.alpha = np.array(self.alpha*len(self.idx), dtype=complex)
            self.gamma = np.array(self.gamma*len(self.idx), dtype=complex)
            return(self.A)

        else:
            self.A = self.data[:, None, None, None]
            # print 'no R interaction to evaluate'
            return(self.A)


class IntT(Inter):
    """ Transmission  interaction class

        Attributes
        ----------

        data = np.array(( i x [theta,si,st]))

    """

    def __init__(self, data=np.array(()), idx=[],slab={}):

        Inter.__init__(self, data=data, idx=idx, typ=3,slab=slab)
        ## index for used slab
        self.uslidx = 0
        # dictionnary of used slab key = slab value = index
        self.dusl = {}
#        self.dusl=dict.fromkeys(self.slab,np.array((),dtype=int))
        self.alpha = []
        self.gamma = []


    def __repr__(self):
        s = 'number of T interaction :' + str(np.shape(self.data)[0])
        return(s)

    def eval(self,fGHz=np.array([2.4])):
        """ evaluate transmission

        Examples
        --------

        >>> from pylayers.antprop.rays import *

        >>> # generate input data
        >>> theta1 = 0.1
        >>> theta2 = 0.4
        >>> si01 = 4
        >>> si02 = 0.6
        >>> sir1 = 3.15
        >>> sir2 = 3.4
        >>> data1=np.array((theta1,si01,sir1))
        >>> data2=np.array((theta2,si02,sir2))

        >>> # store input data to Instance
        >>> T = IntT(data1,idx=0)
        >>> T.data
        array([ 0.1 ,  4.  ,  3.15])
        >>> T.stack(data2,idx=1)
        >>> T.uslidx=1
        >>> T.dusl['AIR']=[0,1]

        >>> # evaluation parameters (normally read from config.ini)
        >>> T.f = np.array([  2.,  11.])
        >>> T.nf = len(T.f)
        >>> T.olf = np.ones((T.nf))

        >>> # evaluation
        >>> eT=T.eval()

        >>> # examples
        >>> ninter = len(T.idx)
        >>> np.shape(eT)
        (1, 2, 2, 2)
        >>> eT[0,0]
        array([[ 1.+0.j,  0.+0.j],
               [ 0.+0.j,  1.+0.j]])
        """

        self.sinsout()

        self.fGHz=fGHz
        self.nf=len(fGHz)

        self.A = np.zeros((self.nf, len(self.idx), 2, 2), dtype=complex)
        self.alpha = np.zeros((len(self.idx)), dtype=complex)
        self.gamma = np.zeros((len(self.idx)), dtype=complex)
        self.sm = np.zeros((len(self.idx)), dtype=complex)

        if np.shape(self.data)[0]!=len(self.idx):
            self.data=self.data.T

        if len(self.data) != 0:
            mapp = []
            for m in self.dusl.keys():
                # used theta of the given slab
                ut = self.data[self.dusl[m], 0]
                if ut.size != 0:
                    if m not in self.slab:
                        m = m.lower()
                    # get alpha and gamma for divergence factor
                    if len(self.slab[m]['lmat']) > 1:
                        print 'Warning : IntR class implemented for mat with only 1 layer '
                    a = 1./np.sqrt(np.array(([self.slab[m]['lmat'][0]['epr']])) \
                               * np.ones(len(ut), dtype=complex))
                    g = (1.-np.sin(ut)**2)/(1.-a*np.sin(ut)**2)
                    try:
                        alpha = np.concatenate((alpha, a), axis=0)
                        gamma = np.concatenate((gamma, g), axis=0)
                    except:
    #                    print 'Warning : alpha or gamma fail'
                        alpha = a
                        gamma = g

                    # find the index of angles which satisfied the data
                    self.slab[m].ev(fGHz=fGHz, theta=ut, RT='T', compensate=True)

                    try:
                        T = np.concatenate((T, self.slab[m].T), axis=1)
                        mapp.extend(self.dusl[m])
                    except:
                        T = self.slab[m].T
                        mapp.extend(self.dusl[m])
            # replace in proper order the Transmission coeff
            self.A[:, np.array((mapp)), :, :] = T
            self.alpha[np.array((mapp))] = alpha
            self.gamma[np.array((mapp))] = gamma
#            self.sm[np.array((mapp))]=sm
            return(self.A)

        else:
            #print 'no T interaction to evaluate'
            return(self.data[:, None, None, None])


class IntD(Inter):
    """ diffraction interaction class
    """
    def __init__(self, data=np.array(()), idx=[],fGHz=np.array([2.4]),slab={}):
        Inter.__init__(self, data=data, idx=idx, typ=1,slab=slab)
        self.dusl = {}
    def __repr__(self):
        s = 'number of D interaction :' + str(np.shape(self.data)[0])
        return s

    def eval(self,fGHz=np.array([2.4])):
        """ evaluate diffraction interaction

        Parameters
        ----------

        fGHz : np.array

        """
        self.fGHz=fGHz
        self.nf=len(fGHz)
        self.A = np.zeros((self.nf, len(self.idx), 2, 2), dtype=complex)

        if len(self.data) != 0 :
            self.phi0 = self.data[:,0]
            self.phi = self.data[:,1]
            self.beta = self.data[:,2]
            self.N = self.data[:,3]
            self.sinsout()
            D = np.zeros([self.nf, len(self.phi), 2, 2], dtype=complex)
            mapp=[]
            for m in self.dusl.keys():
                idx = self.dusl[m]
                mats = m.split('_')
                mat0name = self.slab.di[eval(mats[0])]
                matNname = self.slab.di[eval(mats[1])]
                mat0 = self.slab[mat0name]['lmat'][0]
                matN = self.slab[matNname]['lmat'][0]
                # from IPython.core.debugger import Tracer
                # Tracer()()
                # Ds,Dh = diff(self.fGHz,self.phi0[idx],self.phi[idx],self.si0[idx],self.sout[idx],self.N[idx],mat0,matN,beta=self.beta[idx])
                Ds,Dh = diff(self.fGHz,self.phi0[idx],self.phi[idx],self.si0[idx],self.sout[idx],self.N[idx],mat0,matN,mode='tab',beta=self.beta[idx])
                D[:,idx,0,0]=-Dh
                D[:,idx,1,1]=Ds
                # from IPython.core.debugger import Tracer
                # Tracer()()
                mapp.extend(self.dusl[m])
                # import ipdb
                # ipdb.set_trace()
                # try:
                #     D = np.concatenate((D, self.slab[m].T), axis=1)
                #     mapp.extend(self.dusl[m])
                # except:
                #     D = np.empty([self.nf, len(idx), 2, 2], dtype=complex)
                #     D[:,:,0,0]=Ds
                #     D[:,:,1,1]=Dh
                #     mapp.extend(self.dusl[m])
            self.A[:, np.array((mapp)), :, :] = D[:,mapp,:,:]

            return(self.A)
        else :
            self.A = self.data[:, None, None, None]
            return(self.A)
            print 'not yet implemented'
        # else:
        #     self.A = self.data[:, None, None, None]
        #     print 'no D interaction to evaluate'
        #     return(self.A)

if (__name__ == "__main__"):
    plt.ion()
    print "testing pylayers/antprop/interactions.py"
    doctest.testmod()
