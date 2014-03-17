# -*- coding:Utf-8 -*-
"""

This module handles Scalar and Vector Spherical Harmonics in PyLayers

Class VectorCoeff
==================

.. autosummary::
    :toctree: generated/

    VectorCoeff.__init__
    VectorCoeff.inits1

class SSHCoeff
==============

.. autosummary::
    :toctree: generated/

     SSHCoeff.__init__
     SSHCoeff.s2tos3
     SSHCoeff.sets3


class SCoeff
============

.. autosummary::
    :toctree: generated/

     SCoeff.__init__
     SCoeff.__repr__
     SCoeff.inits2
     SCoeff.inits3
     SCoeff.delete
     SCoeff.put
     SCoeff.delete3
     SCoeff.put3
     SCoeff.s3tos2
     SCoeff.plot
     SCoeff.show

class VCoeff
============

.. autosummary::
    :toctree: generated/

     VCoeff.__init__
     VCoeff.__repr__
     VCoeff.inits1
     VCoeff.inits2
     VCoeff.inits3
     VCoeff.s1tos2
     VCoeff.delete
     VCoeff.put
     VCoeff.delete3
     VCoeff.put3
     VCoeff.s3tos2
     VCoeff.plot
     VCoeff.show

class VSHCoeff
==============

.. autosummary::
    :toctree: generated/

     VSHCoeff.__init__
     VSHCoeff.__repr__
     VSHCoeff.plot
     VSHCoeff.show
     VSHCoeff.s1tos2
     VSHCoeff.s2tos3_new
     VSHCoeff.s2tos3
     VSHCoeff.s3tos2
     VSHCoeff.strip3
     VSHCoeff.ens3
     VSHCoeff.drag3
     VSHCoeff.put3


Utility Functions
=================

.. autosummary::
    :toctree: generated/

     indexssh
     indexvsh
     index_vsh
     AFLegendre3
     AFLegendre2
     AFLegendre
     VW2
     VW
     VW0
     plotVW


"""
import doctest
import subprocess
import os
import re
import sys
import pdb
import numpy as np
import scipy as sp
import scipy.special as special
from scipy import io
import matplotlib.pylab as plt
from scipy.misc import factorial
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.util.plotutil import *
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import axes3d
from scipy import sparse
from matplotlib import rc
from matplotlib import cm

def indexssh(L,mirror=True):
    """ create [l,m] indexation from Lmax

    Parameters
    ----------

    L : maximum order
    mirror : boolean 
        if True the output contains negative m indices

    Returns
    -------

    t : np.array
        [l,m]   Ncoeff x 2

    Examples
    --------

        >>> from pylayers.antprop.spharm import *
        >>> indexssh(2)
        array([[ 0.,  0.],
               [ 1.,  0.],
               [ 2.,  0.],
               [ 1.,  1.],
               [ 2.,  1.],
               [ 2.,  2.],
               [ 1., -1.],
               [ 2., -1.],
               [ 2., -2.]])

    """

    for k in range(L+1):
        l = np.arange(k,L+1)
        m = k*np.ones(L+1-k)
        v = np.vstack((l,m)).T
        try:
            t = np.vstack((t,v))
        except:
            t = v

    if mirror:
        u =  t[L+1:,:]
        v =  np.vstack((u[:,0],-u[:,1])).T
        #v = v[::-1,:]
        t = np.vstack((t,v))

    return t

def indexvsh(L):
    """ indexvsh(L)

    Parameters
    ----------

    L : int
      degree max


    Returns
    -------

    t : ndarray ( (L+1)(L+2)/2 ,  2 )
       tab for indexing the upper triangle

    Examples
    --------

    >>> from pylayers.antprop.antenna import *
    >>> indexvsh(3)
    array([[1, 0],
           [1, 1],
           [2, 0],
           [2, 1],
           [2, 2],
           [3, 0],
           [3, 1],
           [3, 2],
           [3, 3]])

    """
    Kmax = (L + 1) * (L + 2) / 2
    #k = np.arange(Kmax) # old version with 0,0
    k = np.arange(1,Kmax)
    #k = np.arange(Kmax)
    #k = np.arange(Kmax)
    #k = np.arange(1,Kmax)
    l = np.ceil((-1 + np.sqrt(1 + 8 * (k + 1))) / 2) - 1
    m = k - l * (l + 1) / 2
    u = np.vstack((l, m)).T
    t = u.astype(int)
    return(t)

def index_vsh(L, M):
    """ vector spherical harmonics indexing

    Parameters
    ----------

    L : int
        degree max   sum(1..L)   L points
    M : int
        order max    sum(0..M)   M+1 points

    M <=L

    ind[0] = n
    ind[1] = m

    Notes
    -----

    This function is more generic than indexvsh because it allows to have M<>L

    See Also
    --------

    indexvsh

    """
    if M > L:
        print "indexvsh error M>L"

    Kmax1 = (M + 1) * (M + 2) / 2
    #k = np.arange(Kmax1)
    k = np.arange(1,Kmax1)
    l = np.ceil((-1 + np.sqrt(1 + 8 * (k + 1))) / 2) - 1
    m = k - l * (l + 1) / 2
    if (M < L):
        l1 = np.outer(np.arange(L - M) + M + 1, np.ones(M + 1)).ravel()
        m1 = np.outer(np.ones(L - M), np.arange(M + 1)).ravel()
        l = np.hstack((l, l1))
        m = np.hstack((m, m1))

    u = np.vstack((l, m)).T
    t = u.astype(int)
    return(t)

class VectorCoeff(object):

    def __init__(self, typ, fmin=0.6, fmax=6, data=np.array([]),
                 ind=np.array([]) ):
        """
        Parameters
        ----------

        typ :
        fmin :
            min frequency GHz
        fmax :
        data : np.array
        ind : np.array

        Notes 
        -----

        .. warning::
            seems deprecated

        """
        self.s1 = np.array([])
        self.s4 = np.array([])
        self.s3 = np.array([])

        self.fmin = fmin
        self.fmax = fmax

        if typ == 's1':
            self.inits1(data,ind)

    def inits1(self, data, ind):
        """

        Parameters
        ----------

        data :
        ind :

        """
        sh = np.shape(data)
        self.s1 = data
        self.ind_s1 = ind
        self.Nf = sh[0]

class SSHCoeff(object):
    def __init__(self, Cx,Cy,Cz):
        """

        Parameters
        ----------

        Cx : SCoeff
        Cy : SCoeff
        Cz : SCoeff 

        """

        self.Cx = Cx
        self.Cy = Cy
        self.Cz = Cz



    def s2tos3(self, threshold=1e-5):
        """
        convert scalar spherical coefficients from shape 2 to shape 3

        Parameters
        ----------

        threshold : float
            default 1e-20

        Energy thresholded coefficients

        """

        Ex = np.sum(np.abs(self.Cx.s2) ** 2, axis=0) # integrates energy over freq axis = 0 
        Ey = np.sum(np.abs(self.Cy.s2) ** 2, axis=0)
        Ez = np.sum(np.abs(self.Cz.s2) ** 2, axis=0)


        E = Ex + Ey + Ez

        ind = np.nonzero(E > (E.max() * threshold))[0]

        self.Cx.ind3 = self.Cx.ind2[ind]
        self.Cx.s3 = self.Cx.s2[:, ind]
        self.Cx.k2 = ind

        self.Cy.ind3 = self.Cy.ind2[ind]
        self.Cy.s3 = self.Cy.s2[:, ind]
        self.Cy.k2 = ind

        self.Cz.ind3 = self.Cz.ind2[ind]
        self.Cz.s3 = self.Cz.s2[:, ind]
        self.Cz.k2 = ind

    def sets3(self,Cx,Cy,Cz):
        """

        Parameters
        ----------

        Cx : SCoeff
        Cy : SCoeff
        Cz : SCoeff 

        """

        self.Cx.ind3 = Cx.ind3 
        self.Cx.s3 = Cx.s3
        self.Cx.k2 = Cx.k2

        self.Cy.ind3 = Cy.ind3
        self.Cy.s3 = Cy.s3
        self.Cy.k2 = Cy.k2

        self.Cz.ind3 = Cz.ind3
        self.Cz.s3 = Cz.s3
        self.Cz.k2 = Cz.k2


class SCoeff(object):
    """ Scalar Spherical Harmonics Coefficient

    d = np.array [Nf,N+1,M+1]

    Attributes
    ----------

    s2  shape 2   np.array [ Nf x (N+1)*(M+1)   ]
    s3  shape 3   np.array [ Nf x K     ]
    ind [ K x 2]

    """

    def __init__(self, typ='s2', fmin=0.6, fmax=6,lmax=20,  data=np.array([]),
                 ind=np.array([]), k=np.array([])):

    #~ def __init__(self, **kwargs):
        """ init VCoeff

         Parameters
         ----------
         typ : string
            's2' | 's3'
         fmin : float
         fmax : float
         data : ndarray
         ind  : ndarray
         k    : ndarray

         s2 , s3 containers are created
        """

        #~ defaults = { 'typ': 's2',
                    #~ 'fmin' : 0.6,
                    #~ 'fmax' : 6,
                    #~ 'lmax' : 20,
                    #~ 'data' : [],
                    #~ 'ind' : [],
                    #~ 'k'   : [] }
#~ 
        #~ for key, value in defaults.items():
            #~ if key not in kwargs:
                #~ kwargs[key] = value



        self.fmin = fmin
        self.fmax = fmax
        self.lmax = lmax

        if typ == 's2':
            self.s2 = np.array([])
            self.inits2(data,ind)
        if typ == 's3':
            self.s3 = np.array([])
            self.inits3(data, ind, k)

    def __repr__(self):

        st = "Nf   : " +  str(self.Nf) + "\n"
        st = st +  "fmin (GHz) : "+  str(self.fmin) + "\n"
        st = st +  "fmax (GHz) : "+  str(self.fmax) + "\n"

        if 's2' in self.__dict__.keys():
            sh2 = np.shape(self.s2)
            if sh2[0] != 0:
                st = st + "NCoeff s2  : " + str(len(self.ind2))+ "\n"

        if 's3' in self.__dict__.keys():
            sh3 = np.shape(self.s3)
            if sh3[0] != 0:
                st = st + "Ncoeff s3 : " + str(len(self.ind3))+ "\n"



        return(st)

    def inits2(self, data,ind):
        """ initialize shape 2 format

        Parameters
        ----------

        data : shape 2 data

        """

        sh = np.shape(data)
        # first axis is frequency
        self.Nf = sh[0]
        # second axis is the maximum number of coeff

        self.s2 = data

        #self.ind2 = indexssh(lmax)
        self.ind2 = ind

    def inits3(self, data, ind, k):
        """ initialize shape 3 format

        Parameters
        ----------

        data  : shape 3 data
        ind   : shape 3 indexing
        k     : k

        """

        sh = np.shape(data)
        self.Nf = sh[0]
        self.s3 = data
        self.ind3 = ind
        self.k2 = k

    def delete(self, ind, typ):
        """ delete coeff

        Parameters
        ----------

        ind   : int
        typ   : int
                2  shape 2  (Nf , N*M   )
                3  shape 3  (Nf , K )  T ( K x 2 )
        """

        if typ == 2:
            ind2 = self.ind2[ind]
            s2 = self.s2[:, ind]

            a = delete(self.ind2, ind, axis=0)
            b = delete(self.s2, ind, axis=1)
            self.ind2 = a
            self.s2 = b

        if typ == 3:

            ind3 = self.ind3[ind]
            k2 = self.k2[ind]
            s3 = self.s3[:, ind]

            a = delete(self.ind3, ind, axis=0)
            b = delete(self.k2, ind)
            c = delete(self.s3, ind, axis=1)

            self.ind3 = a
            self.k2 = b
            self.s3 = c

    def put(self, typ):
        """ recover last deleted coeff

        Parameters
        ----------

        typ : int
                2 : shape 2  (Nf , N*M   )
                3 : shape 3  (Nf , K )  T ( K x 2 )
        """

        if typ == 2:

            file_ind = pyu.getlong("outfile_i2.txt", pstruc['DIRANT'])
            aux = load(file_ind)
            ind = aux[0]
            ind2 = np.array([aux[1], aux[2]])

            file_s2 = pyu.getlong("outfile_s2.txt", pstruc['DIRANT'])
            s2 = load(file_s2)

            self.s2p = s2

            a = insert(self.ind2, ind, ind2, axis=0)
            b = insert(self.s2, ind, s2, axis=1)

            self.ind2 = a
            self.s2 = b

        if typ == 3:

            file_ind = pyu.getlong("outfile_i3.txt", pstruc['DIRANT'])
            aux = load(file_ind)
            ind = aux[0]
            ind3 = np.array([aux[1], aux[2]])
            k2 = aux[3]

            file_s3 = pyu.getlong("outfile_s3.txt", pstruc['DIRANT'])
            s3 = load(file_s3)

            a = insert(self.ind3, ind, ind3, axis=0)
            b = insert(self.k2, ind, k2)
            c = insert(self.s3, ind, s3[0], axis=1)

            self.ind3 = a
            self.k2 = b
            self.s3 = c

            os.remove(file_ind)
            os.remove(file_s3)

    def delete3(self, ind):
        """ delete3(self,ind): delete coeff.s3

        Parameters
        ----------
        ind :

        """
        a = delete(self.ind3, ind, axis=0)
        b = delete(self.k2, ind)
        c = delete(self.s3, ind, axis=1)
        self.ind3 = a
        self.k2 = b
        self.s3 = c

    def put3(self, i, i3):
        """ put3

        Parameters
        ----------
        i  :
        i3 :
        """

        k2 = i3[0] * (i3[0] + 1) / 2 + i3[1]
        ind3 = self.ind2[k2]
        s3 = self.s2[:, k2]

        a = insert(self.ind3, i, ind3, axis=0)
        b = insert(self.k2, i, k2)
        c = insert(self.s3, i, s3, axis=1)

        self.ind3 = a
        self.k2 = b
        self.s3 = c

    def s3tos2(self):
        """ transform shape3 to shape 2

        s2  shape 2   array [ Nf x (L+1)*(M+1) ]
        s3  shape 3   array [ Nf x K     ] ind [ K x 2]

        Notes
        -----

        The shape of s2 is (Lmax+1)*(Lmax+2)/2

        k2  : is the list of conserved indices in shape 3
        ind3 : np.array (K3, 2) are the conserved (l,m) indices 

        ind3 and k2 have one common dimension

        """
        # retrieve Nf and Lmax to build a void s2 structure
        Nf   = np.shape(self.s3)[0]
        Lmax = max(self.ind3[:,0])
        K2   = (Lmax+1)*(Lmax+2)/2
        self.s2 = np.zeros((Nf,K2),dtype=complex)

        # fill s2 with s3 at proper coefficient location
        self.s2[:,self.k2] = self.s3
        self.N2 = Lmax
        self.M2 = Lmax
        self.ind2 = indexvsh(Lmax)

    def plot(self,typ='s3',title='',xl=False,yl=False,log=False,stem=True,color='b'):
        """
        Parameters
        ----------

        typ : string
            's3'
        title 
        xl 
        yl 
        log
        stem: boolean
        color
        """
        if typ=='s3':
            indices = self.ind3
            tl = indices[:,0]
            C =[]
            for l in np.unique(tl):
                k = np.where(tl==l)
                a = np.real(np.sum(self.s3[:,k]*np.conj(self.s3[:,k])))
                C.append(a)
            C = np.real(np.array(C))
            Cs = np.sqrt(C)
            if log:
                Cs = 20*log10(Cs)
            if stem:
                plt.stem(np.unique(tl),Cs,markerfmt=color+'o')
            else:
                plt.plot(np.unique(tl),Cs,color=color)
            #plt.axis([0,max(tl),0,5])
            plt.title(title)
            if xl:
                plt.xlabel('degree l')
            if yl:
                plt.ylabel('Integrated Module of coeff')

    def show(self,
             typ='s1',
             k = 0,
             N = -1,
             M = -1,
             kmax = 1000,
             seuildb = 50,
             titre = 'SHC',
             xl = True,
             yl = True,
             fontsize=14,
             dB = True,
             cmap = plt.cm.hot_r,
             anim = True):
        """ show coeff

        Parameters
        ----------

        typ :  string
            default ('s1')
            's1'  shape 1  (Nf , N , M )
            's2'  shape 2  (Nf , N*M   )
            's3'  shape 3  (Nf , K )  T ( K x 2 )

        k  : integer
            frequency index default 0

        N, M = maximal value for degree, mode respectively
        (not to be defined if 's2' or 's3')

        """

        fa = np.linspace(self.fmin, self.fmax, self.Nf)
        if typ == 's1':
            if N == -1:
                N = self.N1
            if M == -1:
                M = self.M1
            Mg, Ng = plt.meshgrid(np.arange(M), np.arange(N))
            if anim:
                fig = plt.gcf()
                ax = fig.gca()
                v = np.abs(self.s1[k, 0:N, 0:M])
                if dB:
                    v = 20 * np.log10(v)
                p = plt.scatter(Mg, Ng, c=v, s=30, cmap=cmap,
                            linewidth=0, vmin=-seuildb, vmax=0)
                cb = plt.colorbar()
                cb.set_labe('Level dB')
                plt.draw()
            else:
                v = np.abs(self.s1[k, 0:N, 0:M])
                if dB:
                    vdB = 20 * np.log10(v + 1e-15)
                    plt.scatter(Mg, Ng, c=vdB, s=30, cmap=cmap, linewidth=0,
                                vmin=-seuildb, vmax=0)
                    plt.title(titre)
                    cb = plt.colorbar()
                    cb.set_labe('Level dB')
                else:
                    plt.scatter(Mg, Ng, c=v, s=30, cmap=cmap, linewidth=0)
                    plt.title(titre)
                    cb = plt.colorbar()
                    cb.set_labe('Level (linear scale)')

                if xl:
                    plt.xlabel('m', fontsize=fontsize)
                if yl:
                    plt.ylabel('n', fontsize=fontsize)

        if typ == 's2':
            if np.shape(self.s2)[1] <= 1:
                plt.plot(fa, 10 * np.log10(abs(self.s2[:, 0])))
            else:
                K = np.shape(self.s2)[1]
            kmax = min(kmax,K)
            db = 20 * np.log10(abs(self.s2[:, 0:kmax] + 1e-15))
            col = 1 - (db > -seuildb) * (db + seuildb) / seuildb
            #
            #gray
            #
            #pcolor(np.arange(K+1)[0:kmax],self.fa,col,cmap=cm.gray_r,vmin=0.0,vmax=1.0)
            #
            #color
            #
            plt.pcolor(np.arange(K + 1)[0:kmax], fa, col, cmap=plt.cm.hot, vmin=0.0, vmax=1.0)
            if xl:
                plt.xlabel('index', fontsize=fontsize)
            if yl:
                plt.ylabel('Frequency (GHz)', fontsize=fontsize)

        if typ == 's3':
            if np.shape(self.s3)[1] <= 1:
                plt.plot(fa, 10 * np.log10(abs(self.s3[:, 0])))
            else:
                K = np.shape(self.s3)[1]

            kmax = min(kmax,K)
            db = 20 * np.log10(abs(self.s3[:, 0:kmax] + 1e-15))
            col = 1 - (db > -seuildb) * (db + seuildb) / seuildb
            plt.pcolor(np.arange(K + 1)[0:kmax], fa, col,
                   cmap=plt.cm.hot, vmin=0.0, vmax=1.0)
            if xl:
                plt.xlabel('index', fontsize=fontsize)
            if yl:
                plt.ylabel('Frequency (GHz)', fontsize=fontsize)

                #echelle=[str(0), str(-10), str(-20), str(-30), str(-40), str(-50)]
        if (typ == 's2') | (typ =='s3') :

            echelle = [str(0), str(-seuildb + 40), str(-seuildb + 30),
                       str(-seuildb + 20), str(-seuildb + 10), str(-seuildb)]
            cbar = plt.colorbar(ticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
            cbar.ax.set_yticklabels(echelle)
            cbar.ax.set_ylim(1, 0)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            plt.title(titre, fontsize=fontsize + 2)

class VCoeff(object):
    """ Spherical Harmonics Coefficient

    d = np.array [Nf,N+1,M+1]

    Attributes
    ----------

    s1  shape 1   np.array [ Nf x (N+1) x (M+1) ]
    s2  shape 2   np.array [ Nf x (N+1)*(M+1)   ]
    s3  shape 3   np.array [ Nf x K     ]
    ind [ K x 2]

    """

    def __init__(self, typ, fmin=0.6, fmax=6, data=np.array([]),
                 ind=np.array([]), k=np.array([])):
        """ init VCoeff

         Parameters
         ----------
         typ : string
            's1' | 's2' | 's3'
         fmin : float
         fmax : float
         data : ndarray
         ind  : ndarray
         k    : ndarray

         s1, s2 , s3 containers are created
        """

        self.s1 = np.array([])
        self.s2 = np.array([])
        self.s3 = np.array([])
        self.fmin = fmin
        self.fmax = fmax

        if typ == 's1':
            self.inits1(data)
        if typ == 's2':
            self.inits2(data)
        if typ == 's3':
            self.inits3(data, ind, k)

    def __repr__(self):

        st = "Nf   : " +  str(self.Nf) + "\n"
        st = st +  "fmin (GHz) : "+  str(self.fmin) + "\n"
        st = st +  "fmax (GHz) : "+  str(self.fmax) + "\n"

        sh1 = np.shape(self.s1)
        sh2 = np.shape(self.s2)
        sh3 = np.shape(self.s3)

        if sh1[0] != 0:
            st =  "N1  : " + str(self.N1) + "\n"
            st = st + "M1  : " + str(self.M1)+ "\n"
            st = st + "Ncoeff s1 " + str(self.M1* self.N1)+ "\n"
        if sh2[0] != 0:
            st = st + "NCoeff s2  : " + str(len(self.ind2))+ "\n"
        if sh3[0] != 0:
            st = st + "Ncoeff s3 : " + str(len(self.ind3))+ "\n"

        return(st)

    def inits1(self, data):
        """ initialize shape 1 format

        Parameters
        ----------
        data  : shape 1 data

        """
        sh = np.shape(data)
        N = sh[1] - 1
        M = sh[2] - 1
        if M > N:
            print('VCoeff : M>N ')
            exit()
        else:
            self.s1 = data
            self.N1 = N
            self.M1 = M
        self.Nf = sh[0]

    def inits2(self, data):
        """ initialize shape 2 format

        Parameters
        ----------
         data : shape 2 data

        """
        sh = np.shape(data)
        self.Nf = sh[0]
        kmax = sh[1]
        nmax = np.ceil((-1 + np.sqrt(1 + 8 * (kmax + 1))) / 2) - 1
        t = indexvsh(nmax)
        N2 = t[:, 0].max() - 1
        M2 = t[:, 1].max() - 1
        self.s2 = data
        self.N2 = N2
        self.M2 = M2
        self.ind2 = index_vsh(N2, M2)

    def inits3(self, data, ind, k):
        """ initialize shape 3 format

        Parameters
        ----------
         data  : shape 3 data
         ind   : ishape 3 indexing
         k     : k

        """
        sh = np.shape(data)
        self.Nf = sh[0]
        self.s3 = data
        self.ind3 = ind
        self.k2 = k

    def s1tos2(self, N2=-1):
        """ convert shape 1 --> shape 2

        shape 1   array [ Nf , (L+1) , (M+1) ]
        shape 2   array [ Nf , (L+1) * (M+1) ]

        n = 0...N2
        m = 0...N2

        Parameters
        ----------

        N2 : int <= N1
            shape 1 has 3 axis - shape 2 has 2 axis
            by default all s1 coefficients are kept N2=-1 means N2=min(N1,M1) because M2 must be equal to N2

        See Also
        --------

        index_vsh

        """

        if N2 == -1:
            N2 = min(self.N1, self.M1)
        M2 = N2
        if (N2 <= self.N1):
            self.N2 = N2
            self.M2 = M2
            self.ind2 = index_vsh(N2, M2)
            self.s2 = self.s1[:, self.ind2[:, 0], self.ind2[:, 1]]
        else:
            print('error VCoeff s1tos2: N2>N1')

    def delete(self, ind, typ):
        """ delete coeff

        Parameters
        ----------

        ind   : int
        typ   : int
                2  shape 2  (Nf , N*M   )
                3  shape 3  (Nf , K )  T ( K x 2 )
        """

        if typ == 2:
            ind2 = self.ind2[ind]
            s2 = self.s2[:, ind]

            a = delete(self.ind2, ind, axis=0)
            b = delete(self.s2, ind, axis=1)
            self.ind2 = a
            self.s2 = b

        if typ == 3:

            ind3 = self.ind3[ind]
            k2 = self.k2[ind]
            s3 = self.s3[:, ind]

            a = delete(self.ind3, ind, axis=0)
            b = delete(self.k2, ind)
            c = delete(self.s3, ind, axis=1)
            self.ind3 = a
            self.k2 = b
            self.s3 = c

    def put(self, typ):
        """ recover last deleted coeff

        Parameters
        ----------
        typ : int
                2 : shape 2  (Nf , N*M   )
                3 : shape 3  (Nf , K )  T ( K x 2 )
        """

        if typ == 2:

            file_ind = pyu.getlong("outfile_i2.txt", pstruc['DIRANT'])
            aux = load(file_ind)
            ind = aux[0]
            ind2 = np.array([aux[1], aux[2]])

            file_s2 = pyu.getlong("outfile_s2.txt", pstruc['DIRANT'])
            s2 = load(file_s2)

            self.s2p = s2

            a = insert(self.ind2, ind, ind2, axis=0)
            b = insert(self.s2, ind, s2, axis=1)

            self.ind2 = a
            self.s2 = b

        if typ == 3:

            file_ind = pyu.getlong("outfile_i3.txt", pstruc['DIRANT'])
            aux = load(file_ind)
            ind = aux[0]
            ind3 = np.array([aux[1], aux[2]])
            k2 = aux[3]

            file_s3 = pyu.getlong("outfile_s3.txt", pstruc['DIRANT'])
            s3 = load(file_s3)

            a = insert(self.ind3, ind, ind3, axis=0)
            b = insert(self.k2, ind, k2)
            c = insert(self.s3, ind, s3[0], axis=1)

            self.ind3 = a
            self.k2 = b
            self.s3 = c

            os.remove(file_ind)
            os.remove(file_s3)

    def delete3(self, ind):
        """ delete3(self,ind): delete coeff.s3

        Parameters
        ----------
        ind :

        """
        a = delete(self.ind3, ind, axis=0)
        b = delete(self.k2, ind)
        c = delete(self.s3, ind, axis=1)
        self.ind3 = a
        self.k2 = b
        self.s3 = c

    def put3(self, i, i3):
        """ put3

        Parameters
        ----------
        i  :
        i3 :
        """

        k2 = i3[0] * (i3[0] + 1) / 2 + i3[1]
        ind3 = self.ind2[k2]
        s3 = self.s2[:, k2]

        a = insert(self.ind3, i, ind3, axis=0)
        b = insert(self.k2, i, k2)
        c = insert(self.s3, i, s3, axis=1)

        self.ind3 = a
        self.k2 = b
        self.s3 = c

    def s3tos2(self):
        """ transform shape3 to shape 2

        s2  shape 2   array [ Nf x (L+1)*(M+1) ]
        s3  shape 3   array [ Nf x K     ] ind [ K x 2]

        Notes
        -----

        The shape of s2 is (Lmax+1)*(Lmax+2)/2

        k2  : is the list of conserved indices in shape 3
        ind3 : np.array (K3, 2) are the conserved (l,m) indices 

        ind3 and k2 have one common dimension

        """
        # retrieve Nf and Lmax to build a void s2 structure
        Nf   = np.shape(self.s3)[0]
        Lmax = max(self.ind3[:,0])
        K2   = (Lmax+1)*(Lmax+2)/2
        self.s2 = np.zeros((Nf,K2),dtype=complex)

        # fill s2 with s3 at proper coefficient location
        self.s2[:,self.k2] = self.s3
        self.N2 = Lmax
        self.M2 = Lmax
        self.ind2 = indexvsh(Lmax)

    def plot(self,typ='s3',title='',xl=False,yl=False,log=False,stem=True,color='b'):
        """
        """
        if typ=='s3':
            indices = self.ind3
            tl = indices[:,0]
            C =[]
            for l in np.unique(tl):
                k = np.where(tl==l)
                a = np.real(np.sum(self.s3[:,k]*np.conj(self.s3[:,k])))
                C.append(a)
            C = np.real(np.array(C))
            Cs = np.sqrt(C)
            if log:
                Cs = 20*log10(Cs)
            if stem:
                plt.stem(np.unique(tl),Cs,markerfmt=color+'o')
            else:
                plt.plot(np.unique(tl),Cs,color=color)
            #plt.axis([0,max(tl),0,5])
            plt.title(title)
            if xl:
                plt.xlabel('degree l')
            if yl:
                plt.ylabel('Integrated Module of coeff')

    def show(self,
             typ='s1',
             k = 0,
             N = -1,
             M = -1,
             kmax = 1000,
             seuildb = 50,
             titre = 'SHC',
             xl = True,
             yl = True,
             fontsize=14,
             dB = True,
             cmap = plt.cm.hot_r,
             anim = True):
        """ show coeff

        Parameters
        ----------
        typ :  string 
            default ('s1')
            's1'  shape 1  (Nf , N , M )
            's2'  shape 2  (Nf , N*M   )
            's3'  shape 3  (Nf , K )  T ( K x 2 )

        k  : integer 
            frequency index default 0

        N, M = maximal value for degree, mode respectively
        (not to be defined if 's2' or 's3')

        """

        fa = np.linspace(self.fmin, self.fmax, self.Nf)
        if typ == 's1':
            if N == -1:
                N = self.N1
            if M == -1:
                M = self.M1
            Mg, Ng = plt.meshgrid(np.arange(M), np.arange(N))
            if anim:
                fig = plt.gcf()
                ax = fig.gca()
                v = np.abs(self.s1[k, 0:N, 0:M])
                if dB:
                    v = 20 * np.log10(v)
                p = plt.scatter(Mg, Ng, c=v, s=30, cmap=cmap,
                            linewidth=0, vmin=-seuildb, vmax=0)
                plt.colorbar()
                plt.draw()
            else:
                v = np.abs(self.s1[k, 0:N, 0:M])
                if dB:
                    vdB = 20 * np.log10(v + 1e-15)
                    plt.scatter(Mg, Ng, c=vdB, s=30, cmap=cmap, linewidth=0,
                                vmin=-seuildb, vmax=0)
                    plt.title(titre)
                    plt.colorbar()
                else:
                    plt.scatter(Mg, Ng, c=v, s=30, cmap=cmap, linewidth=0)
                    plt.title(titre)
                    plt.colorbar()

                if xl:
                    plt.xlabel('m', fontsize=fontsize)
                if yl:
                    plt.ylabel('n', fontsize=fontsize)

        if typ == 's2':
            if np.shape(self.s2)[1] <= 1:
                plt.plot(fa, 10 * np.log10(abs(self.s2[:, 0])))
            else:
                K = np.shape(self.s2)[1]
            
            kmax = min(kmax,K)
            db = 20 * np.log10(abs(self.s2[:, 0:kmax] + 1e-15))
            col = 1 - (db > -seuildb) * (db + seuildb) / seuildb
            #
            #gray
            #
            #pcolor(np.arange(K+1)[0:kmax],self.fa,col,cmap=cm.gray_r,vmin=0.0,vmax=1.0)
            #
            #color
            #
            plt.pcolor(np.arange(K + 1)[0:kmax], fa, col, cmap=plt.cm.hot, vmin=0.0, vmax=1.0)
            if xl:
                plt.xlabel('index', fontsize=fontsize)
            if yl:
                plt.ylabel('Frequency (GHz)', fontsize=fontsize)

        if typ == 's3':
            if np.shape(self.s3)[1] <= 1:
                plt.plot(fa, 10 * np.log10(abs(self.s3[:, 0])))
            else:
                K = np.shape(self.s3)[1]

            kmax = min(kmax,K)
            db = 20 * np.log10(abs(self.s3[:, 0:kmax] + 1e-15))
            col = 1 - (db > -seuildb) * (db + seuildb) / seuildb
            plt.pcolor(np.arange(K + 1)[0:kmax], fa, col,
                   cmap=plt.cm.hot, vmin=0.0, vmax=1.0)
            if xl:
                plt.xlabel('index', fontsize=fontsize)
            if yl:
                plt.ylabel('Frequency (GHz)', fontsize=fontsize)

                #echelle=[str(0), str(-10), str(-20), str(-30), str(-40), str(-50)]
        if (typ == 's2') | (typ =='s3') :

            echelle = [str(0), str(-seuildb + 40), str(-seuildb + 30), 
                       str(-seuildb + 20), str(-seuildb + 10), str(-seuildb)]
            cbar = plt.colorbar(ticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
            cbar.ax.set_yticklabels(echelle)
            cbar.ax.set_ylim(1, 0)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            plt.title(titre, fontsize=fontsize + 2)

class VSHCoeff(object):
    """ Vector Spherical Harmonics Coefficients class


    Attributes
    ----------

    Bi
    Br
    Ci
    Cr


    Notes
    ------

    Br = VCoeff(br)
    Bi = VCoeff(bi)
    Cr = VCoeff(cr)
    Ci = VCoeff(ci)

    C  = VSHCoeff(Br,Bi,Cr,Ci)

    """
    def __init__(self, Br, Bi, Cr, Ci):
        """ Init  VSHCoeff

        Parameters
        ----------

        Br
        Bi
        Cr
        Ci
        """

        self.Br = Br
        self.Bi = Bi
        self.Cr = Cr
        self.Ci = Ci

    def __repr__(self):
        """ 
        """
        st = "Br"+'\n'
        st = st + "-------------"+'\n'
        st = st + self.Br.__repr__()+'\n'
        st = st + "Bi"+'\n'
        st = st + "-------------"+'\n'
        st = st + self.Bi.__repr__()+'\n'
        st = st + "Cr"+'\n'
        st = st + "-------------"+'\n'
        st = st + self.Cr.__repr__()+'\n'
        st = st + "Ci"+'\n'
        st = st + "-------------"+'\n'
        st = st + self.Ci.__repr__()
        return(st) 

    def plot(self,typ='s3',titre='titre',log=False,stem=True,subp=True):
        """ plot coeff

        Parameters
        ----------

        typ
        titre
        log
        stem
        subp

        """
        fa = np.linspace(self.Br.fmin,self.Br.fmax,self.Br.Nf)
        st = titre+'  shape : '+typ
        plt.suptitle(st,fontsize=14)
        if subp:
            plt.subplot(221)
            titre = '$\sum_f \sum_m |Br_{l}^{(m)}(f)|$'
            self.Br.plot(typ=typ,title=titre, yl=True,color='r',stem=stem,log=log)
        else:
            self.Br.plot(typ=typ,color='r',stem=stem,log=log)
        if subp:
            plt.subplot(222)
            titre = '$\sum_f \sum_m |Bi_{l}^{(m)}(f)|$'
            self.Bi.plot(typ=typ,title=titrei,color='m',stem=stem,log=log)
        else:
            self.Bi.plot(typ=typ,color='m',stem=stem,log=log)
        if subp:
            plt.subplot(223)
            titre = '$\sum_f \sum_m |Cr_{l}^{(m)}(f)|$'
            self.Cr.plot(typ=typ,title=titre, xl=True, yl=True,color='b',stem=stem,log=log)
        else:
            self.Cr.plot(typ=typ,color='b',stem=stem,log=log)
        if subp:
            plt.subplot(224)
            titre = '$\sum_f \sum_m |Ci_{l}^{(m)}(f)|$'
            self.Ci.plot(typ=typ, title = titre, xl=True,color='c',stem=stem,log=log)
        else:
            self.Ci.plot(typ=typ,xl=True,yl=True,color='c',stem=stem,log=log)
        if not subp:
            plt.legend(('$\sum_f \sum_m |Br_{l}^{(m)}(f)|$',
                        '$\sum_f \sum_m |Bi_{l}^{(m)}(f)|$',
                        '$\sum_f \sum_m |Cr_{l}^{(m)}(f)|$',
                        '$\sum_f \sum_m |Ci_{l}^{(m)}(f)|$'))

    def show(self, typ='s1', k=1, N=-1, M=-1, kmax = 1000, seuildb=50,
             animate=False,titre=''):
        """ show VSH coeff

        Parameters
        ----------

        typ : str
            {'s1','s2','s3'}
        k  : int
            frequency index
        kmax : int
            maximum of the unfolded coefficient axes
        N  : int
        M  : int
        seuildB  : float
        animate : boolean
                default False
        """
        if not animate:
            fa = np.linspace(self.Br.fmin,self.Br.fmax,self.Br.Nf)
            st = titre+'  shape : '+typ
            plt.suptitle(st,fontsize=14)
            plt.subplot(221)
            titre = '$|Br_{n}^{(m)}|$'
            self.Br.show(typ=typ,titre=titre, xl=False, yl=True)
            plt.subplot(222)
            titre = '$|Bi_{n}^{(m)}|$'
            self.Bi.show(typ=typ,titre=titre, xl=False, yl=False)
            plt.subplot(223)
            titre = '$|Cr_{n}^{(m)}|$'
            self.Cr.show(typ=typ,titre=titre, xl=True, yl=True)
            plt.subplot(224)
            titre = '$|Ci_{n}^{(m)}|$'
            self.Ci.show(typ=typ, titre = titre, xl=True, yl=False)
        else:
            for k in np.arange(self.Br.Nf):
                plt.subplot(221)
                titre = '$|Br_{n}^{(m)}|$'
                self.Br.show(typ, titre = titre, xl=False, yl=True)
                plt.subplot(222)
                titre = '$|Bi_{n}^{(m)}|$'
                self.Bi.show(typ, titre = titre, xl=False, yl=False)
                plt.subplot(223)
                titre = '$|Cr_{n}^{(m)}|$'
                self.Cr.show(typ, titre = titre , xl=True, yl=True)
                plt.subplot(224)
                titre = '$|Ci_{n}^{(m)}|$'
                self.Ci.show(typ, titre = titre , xl=True, yl=False)
    #    show()

    def s1tos2(self, N2=-1):
        """ convert shape 1 to shape 2

        shape 1   array [ Nf x (L+1) x (M+1) ]
        shape 2   array [ Nf x (L+1)*(M+1)   ]

        Parameters
        ----------

        N2 : max level
            default (-1 means all values)

        """
        self.Bi.s1tos2(N2)
        self.Br.s1tos2(N2)
        self.Ci.s1tos2(N2)
        self.Cr.s1tos2(N2)

    def s2tos3_new(self, k):
        """ convert vector spherical coefficient from shape 2 to shape 3

        Parameters
        ----------

        k : number of coeff

        """

        EBr = np.sum(np.abs(self.Br.s2) ** 2, axis=0)
        EBi = np.sum(np.abs(self.Bi.s2) ** 2, axis=0)
        ECr = np.sum(np.abs(self.Cr.s2) ** 2, axis=0)
        ECi = np.sum(np.abs(self.Ci.s2) ** 2, axis=0)

        E  = EBr + EBi + ECr + ECi

        ib = np.argsort(E)[::-1]
        
        print self.Br.ind2[ib[k-1]]
        print self.Cr.ind2[ib[k-1]]
        print self.Ci.ind2[ib[k-1]]
        print self.Bi.ind2[ib[k-1]]
        #ind = np.nonzero(E > (E.max() * threshold))[0]
        self.Br.ind3 = self.Br.ind2[ib[range(k)]]
        self.Br.s3 = self.Br.s2[:, ib[range(k)]]
        self.Br.k2 = ib[range(k)]

        self.Bi.ind3 = self.Bi.ind2[ib[range(k)]]
        self.Bi.s3 = self.Bi.s2[:, ib[range(k)]]
        self.Bi.k2 = ib[range(k)]

        self.Cr.ind3 = self.Cr.ind2[ib[range(k)]]
        self.Cr.s3 = self.Cr.s2[:, ib[range(k)]]
        self.Cr.k2 = ib[range(k)]

        self.Ci.ind3 = self.Ci.ind2[ib[range(k)]]
        self.Ci.s3 = self.Ci.s2[:, ib[range(k)]]
        self.Ci.k2 = ib[range(k)]
        return E[ib[k-1]]
    
    def s2tos3(self, threshold=1e-5):
        """ convert vector spherical coefficients from shape 2 to shape 3

        Parameters
        ----------

        threshold : float
            default 1e-20

        Energy thresholded coefficients

        """

        EBr = np.sum(np.abs(self.Br.s2) ** 2, axis=0) # integrates energy over freq axis = 0 
        EBi = np.sum(np.abs(self.Bi.s2) ** 2, axis=0)
        ECr = np.sum(np.abs(self.Cr.s2) ** 2, axis=0)
        ECi = np.sum(np.abs(self.Ci.s2) ** 2, axis=0)

        E = EBr + EBi + ECr + ECi

        ind = np.nonzero(E > (E.max() * threshold))[0]

        self.Br.ind3 = self.Br.ind2[ind]
        self.Br.s3 = self.Br.s2[:, ind]
        self.Br.k2 = ind

        self.Bi.ind3 = self.Bi.ind2[ind]
        self.Bi.s3 = self.Bi.s2[:, ind]
        self.Bi.k2 = ind

        self.Cr.ind3 = self.Cr.ind2[ind]
        self.Cr.s3 = self.Cr.s2[:, ind]
        self.Cr.k2 = ind

        self.Ci.ind3 = self.Ci.ind2[ind]
        self.Ci.s3 = self.Ci.s2[:, ind]
        self.Ci.k2 = ind



    def s3tos2(self):
        """
        s3tos2
        """
        self.Br.s3tos2()
        self.Bi.s3tos2()
        self.Cr.s3tos2()
        self.Ci.s3tos2()

    def strip3(self):
        """ Thresholded coefficient conversion

        The s3 minimmum energy coefficient is deleted

        Returns
        -------
           ind
           ind3
        """
        EBr = sum(abs(self.Br.s3) ** 2, axis=0)
        EBi = sum(abs(self.Bi.s3) ** 2, axis=0)
        ECr = sum(abs(self.Cr.s3) ** 2, axis=0)
        ECi = sum(abs(self.Ci.s3) ** 2, axis=0)

        E = EBr + EBi + ECr + ECi

        Emin = min(E)
        ind  = find(E == Emin)
        ind3 = self.Br.ind3[ind]

        self.Br.delete3(ind)
        self.Bi.delete3(ind)
        self.Cr.delete3(ind)
        self.Ci.delete3(ind)

        return ind, ind3

    def ens3(self):
        """ return sorted energy values from minimal to maximal value

        Returns
        -------
        Es
            sorted energy values
        u
            index
        """
        EBr = np.sum(np.abs(self.Br.s3) ** 2, axis=0)
        EBi = np.sum(np.abs(self.Bi.s3) ** 2, axis=0)
        ECr = np.sum(np.abs(self.Cr.s3) ** 2, axis=0)
        ECi = np.sum(np.abs(self.Ci.s3) ** 2, axis=0)

        E = EBr + EBi + ECr + ECi
        u = np.argsort(E)
        Es = E[u]
        return(Es,u)

    def drag3(self, Emin):
        """ Thresholded coefficient conversion

        Parameters
        ----------
        Emin : Minimum energy

        """
        EBr = sum(abs(self.Br.s3) ** 2, axis=0)
        EBi = sum(abs(self.Bi.s3) ** 2, axis=0)
        ECr = sum(abs(self.Cr.s3) ** 2, axis=0)
        ECi = sum(abs(self.Ci.s3) ** 2, axis=0)
        E = EBr + EBi + ECr + ECi

        ind = find(E == Emin)

        ind3 = self.Br.ind3[ind]

        self.Br.delete3(ind)
        self.Bi.delete3(ind)
        self.Cr.delete3(ind)
        self.Ci.delete3(ind)

        return ind, ind3

    def put3(self, i, i3):
        """
        """
        self.Br.put3(i, i3)
        self.Bi.put3(i, i3)
        self.Cr.put3(i, i3)
        self.Ci.put3(i, i3)

def AFLegendre3(L, M, x):
    """ calculate Pmm1l and Pmp1l

    Parameters
    ----------
        L : int
            max order  (theta)   (also called l or level )
        M : int
            max degree (phi)
        x : np.array
            function argument

    Returns
    -------

    Pmm1l : ndarray (Nx , L , M )
        :math:`\\bar{P}_{l}^{(m-1)}(x)`

    Pmp1l : ndarray (Nx , L , M )
        :math:`\\bar{P}_{l}^{(m+1)}(x)`

    Notes
    -----

    This function returns :
        .. math::

            \\bar{P}_{l}^{(m-1)}(x)

            \\bar{P}_{l}^{(m+1)}(x)

     Where

        .. math::

            P_l^{(m)}(x)= \\sqrt{ \\frac{2}{2 l+1} \\frac{(l+m)!}{(l-m)!} } \\bar{P}_{l}^{(m)}(x)

    Examples
    --------


    >>> Pmm1l,Pmp1l = AFLegendre3(5,4,np.array([0,1]))

    Notes
    -----

    L has to be greater or equal than M

    See Also
    --------

    VW

    """
    PML = []
    nx = len(x)

    if M < L:
        MM = np.arange(M + 2).reshape(M+2,1,1)
        LL = np.arange(L + 1).reshape(1,L+1,1)
    else:
        MM = np.arange(M + 1).reshape(M+1,1,1)
        LL = np.arange(L + 1).reshape(1,L+1,1)

    x  = x.reshape(1,1,nx)

    #
    # Warning : this is a dangerous factorial ratio
    # surprinsingly it works well
    #
    C1 = np.sqrt((LL + 0.5) * factorial(LL - MM) / factorial(LL + MM))
    Pml = special.lpmv(MM,LL,x)*C1

    Pml = np.swapaxes(Pml,0,2)
    Pml = np.swapaxes(Pml,1,2)

    if M < L:
        Pmp1l = Pml[:, 1::1, :]
    else:
        Pmp1l = np.zeros((nx, M + 1, L + 1))
        Pmp1l[:, 0:-1, :] = Pml[:, 1::1, :]

    Pmm1l = np.zeros((nx, M + 1, L + 1))
    if M < L:
        Pmm1l[:, 1::1, :] = Pml[:, 0:-2, :]
    else:
        Pmm1l[:, 1::1, :] = Pml[:, 0:-1, :]
        Pmm1l[:, 0, :] = -Pml[:, 1, :]

    return Pmm1l, Pmp1l

def AFLegendre2(L, M, x):
    """ calculate Pmm1l and Pmp1l

    Parameters
    ----------
        L : int
            max order  (theta)   (also called l or level )
        M : int
            max degree (phi)
        x : np.array
            function argument

    Returns
    -------

    Pmm1l : ndarray (Nx , L , M )
        :math:`\\bar{P}_{l}^{(m-1)}(x)`

    Pmp1l : ndarray (Nx , L , M )
        :math:`\\bar{P}_{l}^{(m+1)}(x)`

    Notes
    -----

    This function returns :
        .. math::

            \\bar{P}_{l}^{(m-1)}(x)

            \\bar{P}_{l}^{(m+1)}(x)

     Where

        .. math::

            P_l^{(m)}(x)= \\sqrt{ \\frac{2}{2 l+1} \\frac{(l+m)!}{(l-m)!} } \\bar{P}_{l}^{(m)}(x)

    Examples
    --------


    >>> Pmm1l,Pmp1l = AFLegendre2(5,4,np.array([0,1]))

    Notes
    -----

    L has to be greater or equal than M

    See Also
    --------

    VW

    """
    PML = []
    nx = len(x)
    if M < L:
        MM = np.expand_dims(np.arange(M + 2),1)
        LL = np.expand_dims(np.arange(L + 1),0)
    else:
        MM = np.expand_dims(np.arange(M + 1),1)
        LL = np.expand_dims(np.arange(L + 1),0)
    #
    # Warning : this is a dangerous factorial ratio
    # surprinsingly it works well
    #
    C1 = np.sqrt((LL + 0.5) * factorial(LL - MM) / factorial(LL + MM))
    for i in range(nx):
        if M < L:
            pml = special.lpmn(M + 1, L, x[i])[0]
        else:
            pml = special.lpmn(M, L, x[i])[0]
        pml = pml * C1
        PML.append(pml)

    Pml = np.array(PML)
    if M < L:
        Pmp1l = Pml[:, 1::1, :]
    else:
        Pmp1l = np.zeros((nx, M + 1, L + 1))
        Pmp1l[:, 0:-1, :] = Pml[:, 1::1, :]

    Pmm1l = np.zeros((nx, M + 1, L + 1))
    if M < L:
        Pmm1l[:, 1::1, :] = Pml[:, 0:-2, :]
    else:
        Pmm1l[:, 1::1, :] = Pml[:, 0:-1, :]
        Pmm1l[:, 0, :] = -Pml[:, 1, :]

    return Pmm1l, Pmp1l

def AFLegendre(N, M, x):
    """ calculate Pmm1n and Pmp1n

    Parameters
    ----------

    N : int
        max order  (theta)   (also called l or level )
    M : int
        max degree (phi)
    x : np.array
            function argument

    Returns
    -------

    Pmm1l :  ndarray ( Ndir, M , L ) 

        :math:`\\bar{P}_{n}^{(m-1)}(x)`
    Pmp1l :  ndarray ( Ndir, M , L )
        :math:`\\bar{P}_{n}^{(m+1)}(x)`

    Notes
    -----

    This function returns :
        .. math::

            \\bar{P}_{l}^{(m-1)}(x)

            \\bar{P}_{l}^{(m+1)}(x)

     Where

        .. math::

            P_l^{(m)}(x)= \\sqrt{ \\frac{2}{2 l+1} \\frac{(l+m)!}{(l-m)!} } \\bar{P}_{l}^{(m)}(x)

    Examples
    --------

    >>> Pmm1n,Pmp1n = AFLegendre(5,4,np.array([0,1]))

    See Also
    --------

    VW

    """
    PMN = []
    nx = len(x)
    if M < N:
        MM = np.outer(np.arange(M + 2), np.ones(N + 1))
        NN = np.outer(np.ones(M + 2), np.arange(N + 1))
    else:
        MM = np.outer(np.arange(M + 1), np.ones(N + 1))
        NN = np.outer(np.ones(M + 1), np.arange(N + 1))
    #
    # Warning : this is a dangerous factorial ratio
    # surprinsingly it works well
    #
    C1 = np.sqrt((NN + 0.5) * factorial(NN - MM) / factorial(NN + MM))
    del MM
    del NN
    for i in range(nx):
        if M < N:
            pmn = special.lpmn(M + 1, N, x[i])[0]
        else:
            pmn = special.lpmn(M, N, x[i])[0]
        pmn = pmn * C1
        PMN.append(pmn)

    Pmn = np.array(PMN)
    if M < N:
        Pmp1n = Pmn[:, 1::1, :]
    else:
        Pmp1n = np.zeros((nx, M + 1, N + 1))
        Pmp1n[:, 0:-1, :] = Pmn[:, 1::1, :]

    Pmm1n = np.zeros((nx, M + 1, N + 1))
    if M < N:
        Pmm1n[:, 1::1, :] = Pmn[:, 0:-2, :]
    else:
        Pmm1n[:, 1::1, :] = Pmn[:, 0:-1, :]
        Pmm1n[:, 0, :] = -Pmn[:, 1, :]

    return Pmm1n, Pmp1n

def VW2(l, m, x, phi, Pmm1l, Pmp1l):
    """ evaluate vector Spherical Harmonics basis functions

    Parameters
    ----------
    l    : ndarray (1 x K)
        level
    m    : ndarray (1 x K)
        mode
    x    :  ndarray (1 x Nray)

    phi   : np.array (1 x Nray)

    Pmm1l : Legendre Polynomial

    Pmp1l : Legendre Polynomial

    Returns
    -------

    V  : ndarray (Nray , L, M)
    W  : ndarray (Nray , L, M)

    See Also
    --------

    AFLegendre

    Nx x M x L

    Examples
    --------

    """

    K   = len(l)
    Nr  = len(x)
    l   = l.reshape(1,K)
    m   = m.reshape(1,K)
    phi = phi.reshape(Nr,1)
    x   = x.reshape(Nr,1)

    t1 = np.sqrt((l + m) * (l - m + 1))
    t2 = np.sqrt((l - m) * (l + m + 1))

    Ephi = np.exp(1j*m*phi)

    Y1 = (t1 * Pmm1l[:,m,l] + t2 * Pmp1l[:,m,l]).reshape(Nr,K)
    Y2 = (t1 * Pmm1l[:,m,l] - t2 * Pmp1l[:,m,l]).reshape(Nr,K)

    W = Y1 * (-1.0) ** l / (2 * x * np.sqrt(l * (l + 1))) * Ephi
    W[np.isinf(W) | np.isnan(W)] = 0
    V = Y2 * (-1.0) ** l / (2 * np.sqrt(l * (l + 1))) * Ephi
    V[np.isinf(V) | np.isnan(V)] = 0
    return V, W

def VW(l, m, theta ,phi):
    """ evaluate vector Spherical Harmonics basis functions

    Parameters
    ----------
    l    : ndarray (1 x K)
        level
    m    : ndarray (1 x K)
        mode
    theta : np.array (1 x Nray)

    phi   : np.array (1 x Nray)


    Returns
    -------

    V  : ndarray (Nray , L, M)
    W  : ndarray (Nray , L, M)

    See Also
    --------

    AFLegendre

    Nray x M x L

    Examples
    --------
   
        >>> a = 1

    """

    if type(l) == float:
        l = np.array([l])
    if type(m) == float:
        m = np.array([m])

    L = np.max(l)
    M = np.max(m)
 
    # dirty fix
    index = np.where(abs(theta-np.pi/2)<1e-5)[0]
    if len(index)>0:
        theta[index]=np.pi/2-0.01
    x = -np.cos(theta)

    # The - sign is necessary to get the good reconstruction
    #     deduced from observation
    #     May be it comes from a different definition of theta in SPHEREPACK

    #Pmm1l, Pmp1l = AFLegendre(L, M, x)

    Pmm1l, Pmp1l = AFLegendre(L, L, x)

    K   = len(l)
    Nr  = len(x)

    l   = l.reshape(1,K)
    m   = m.reshape(1,K)
    phi = phi.reshape(Nr,1)
    x   = x.reshape(Nr,1)

    t1 = np.sqrt((l + m) * (l - m + 1))
    t2 = np.sqrt((l - m) * (l + m + 1))

    Ephi = np.exp(1j*m*phi)

    Y1 = (t1 * Pmm1l[:,m,l] + t2 * Pmp1l[:,m,l]).reshape(Nr,K)
    Y2 = (t1 * Pmm1l[:,m,l] - t2 * Pmp1l[:,m,l]).reshape(Nr,K)

    T =  (-1.0) ** l / (2 * np.sqrt(l * (l + 1))) * Ephi
    #W = Y1 * (-1.0) ** l / (2 * x * np.sqrt(l * (l + 1))) * Ephi
    #V = Y2 * (-1.0) ** l / (2 * np.sqrt(l * (l + 1))) * Ephi
    W = Y1 * T / x
    V = Y2 * T
    #
    # dirty fix
    #
    #W[np.isinf(W) | np.isnan(W)] = 0
    #V[np.isinf(V) | np.isnan(V)] = 0

    return V, W

def VW0(n, m, x, phi, Pmm1n, Pmp1n):
    """ evaluate vector Spherical Harmonics basis functions

    Parameters
    ----------
    n    : int
        level
    m    : int
        mode
    x    :  np.array
        function argument
    phi   : np.array
    Pmm1n : Legendre Polynomial
    Pmp1n : Legendre Polynomial


    Returns
    -------

    V
    W


    Examples
    --------

    >>> from pylayers.antprop.antenna import *
    >>> N = 2
    
    See Also
    --------

    AFLegendre


    """
    t1 = np.outer(np.ones(len(x)), np.sqrt((n + m) * (n - m + 1)))
    t2 = np.outer(np.ones(len(x)), np.sqrt((n - m) * (n + m + 1)))
    Y1 = t1 * Pmm1n[:, m, n] + t2 * Pmp1n[:, m, n]
    Y2 = t1 * Pmm1n[:, m, n] - t2 * Pmp1n[:, m, n]

    Mphi = np.outer(phi, m)
    Ephi = np.exp(1j * Mphi)
    del Mphi
    Y1 = t1 * Pmm1n[:, m, n] + t2 * Pmp1n[:, m, n]
    Y2 = t1 * Pmm1n[:, m, n] - t2 * Pmp1n[:, m, n]
    del t1
    del t2
    W = Y1 * np.outer(1.0 / x, (-1.0) ** n / (2 * np.sqrt(n * (n + 1)))) * Ephi
    W[np.isinf(W) | np.isnan(W)] = 0
    del Y1
    V = Y2 * np.outer( np.ones(len(x)), (-1.0) ** n / (2 * np.sqrt(n * (n + 1)))) * Ephi
    V[np.isinf(V) | np.isnan(V)] = 0
    del Y2
    return V, W

def plotVW(l, m, theta, phi, sf=False):
    """ plot VSH transform vsh basis in 3D plot

        (V in fig1 and W in fig2)

    Parameters
    ----------

    n,m   : integer values (m<=n)
    theta : ndarray
    phi   : ndarray
    sf    : boolean
        if sf : plotted figures are saved in a *.png file
        else  : plotted figures aren't saved

    Examples
    --------

    .. plot::
        :include-source:

        >>> from pylayers.antprop.spharm import *
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> n=5
        >>> m=3
        >>> theta = np.linspace(0,np.pi,30)
        >>> phi   = np.linspace(0,2*np.pi,60)
        >>> plotVW(n,m,theta,phi)
        >>> plt.show()

    """
    # calculate v and w
    if m <= l:
        theta[np.where(theta == np.pi / 2)[0]] = np.pi / 2 +  1e-10  # .. todo :: not clean
        x = -np.cos(theta)

        Pmm1l, Pmp1l = AFLegendre(l, m, x)

        t1 = np.sqrt((l + m) * (l - m + 1))
        t2 = np.sqrt((l - m) * (l + m + 1))
        y1 = t1 * Pmm1l[:, m, l] - t2 * Pmp1l[:, m, l]
        y2 = t1 * Pmm1l[:, m, l] + t2 * Pmp1l[:, m, l]

        Ephi = np.exp(1j * m * phi)
        cphi = np.cos(m * phi)
        if m == 0:
            sphi = 1e-10
        else:
            sphi = np.sin(m * phi)

        ny = len(y1)
        ne = len(Ephi)
        vy = np.ones(ny)
        ve = np.ones(ne)
        Y1 = np.outer(y1, ve)
        Y2 = np.outer(y2, ve)
        EPh = np.outer(vy, Ephi)

        const = (-1.0) ** l / (2 * np.sqrt(l * (l + 1)))
        V = const * Y1 * EPh
        #V[np.isinf(V)|isnan(V)]=0
        Vcos = cphi * V
        Vsin = sphi * V

        if m == 0:
            #W=np.zeros((len(theta),len(phi)))
            W = np.ones((len(theta), len(phi))) * 1e-10
        else:
            Waux = Y2 * EPh
            x1 = 1.0 / x
            W = np.outer(x1, const) * Waux

        Wcos = cphi * W
        Wsin = sphi * W

        #figdirV='/home/rburghel/Bureau/bases_decomposition_VW/base_V_Vsin_Vcos/'
        figdirV = './'
        ext1 = '.pdf'
        ext2 = '.eps'
        ext3 = '.png'
        slm = ' l = '+str(l)+' m = '+str(m)
        fig1 = plt.figure()
        pol3D(fig1,abs(V),theta,phi,title='$|V|$'+slm)
        fig2 = plt.figure()
        pol3D(fig2,abs(Vcos),theta,phi,title='$\Re V$'+slm)
        fig3 = plt.figure()
        pol3D(fig3,abs(Vsin),theta,phi,title='$\Im V$'+slm)
        fig4 = plt.figure()
        pol3D(fig4,abs(W),theta,phi,title='$|W|$'+slm)
        fig5 = plt.figure()
        pol3D(fig5,abs(Wcos),theta,phi,title='$\Re W'+slm)
        fig6 = plt.figure()
        pol3D(fig6,abs(Wsin),theta,phi,title='$\Im W$'+slm)
        plt.show()

    else:
        print "Error: m>n!!!"


if (__name__ == "__main__"):
    doctest.testmod()
