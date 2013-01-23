# -*- coding:Utf-8 -*-
"""
This module handles antennas in pylayers

To instantiate an antenna object : 

.. python::

    A = Antenna(typ,_filename,directory,nf,ntheta,nphi)

typ indicates the antenna file format to read 

Examples
--------
    >>> from pylayers.antprop.antenna import *
    >>> A = Antenna('mat','S1R1.mat','ant/UWBAN/Matfile')

The antenna can be represented in various formats

.vsh2
.vsh3



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
#from spharm import Spharmt,getspecindx
from pylayers.util.project import *
from sphere import spherepack, Wrapec, mathtogeo

from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import axes3d
from scipy import sparse
from matplotlib import rc
from matplotlib import cm


def indexvsh(L):
    """ indexvsh(L)

    Parameters
    ----------
         L : degree max


    Returns
    -------
        t : ndarray ( (L+1)(L+2)/2 ,  2 )
            tab for indexing the upper triangle
    Examples
    --------

        >>> from pylayers.antprop.antenna import *
        >>> indexvsh(3)
        array([[0, 0],
               [1, 0],
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
    #k = np.arange(Kmax)
    k = np.arange(Kmax)
    l = np.ceil((-1 + np.sqrt(1 + 8 * (k + 1))) / 2) - 1
    m = k - l * (l + 1) / 2
    u = np.vstack((l, m)).T
    t = u.astype(int)
    return(t)


def index_vsh(L, M):
    """ vector sperical harmonics indexing

    Parameters
    ----------
    L : int
        degree max   sum(1..L)   L points
    M : int
        order max    sum(0..M)   M+1 points

    M <=L

    ind[0] = n
    ind[1] = m

    """
    if M > L:
        print "indexvsh error M>L"

    Kmax1 = (M + 1) * (M + 2) / 2
    k = np.arange(Kmax1)
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


def geom_pattern(theta, phi, E, f, p, minr, maxr, racine, ilog=False):
    """ export antenna pattern in geomview format

    Parameters
    ----------
    theta : np.array (1 x Ntheta)
    phi   : np.array (1 x Nphi)
    E     : np.array complex  (Ntheta,Nphi)
    f     : frequency
    po    : origin (1x3)
    minr  : radius of minimum
    maxr  : radius of maximum
    ilog  : True (log) False (linear)


    """
    Nt = len(theta)
    Np = len(phi)

    if ilog:
        R = 10 * np.log10(abs(E))
    else:
        R = abs(E)

    Th = np.outer(theta, np.ones(Np))
    Ph = np.outer(np.ones(Nt), phi)

    T = (R - minr) / (maxr - minr)
    Ry = 5 + 2 * T
    x = Ry * np.sin(Th) * np.cos(Ph) + p[0]
    y = Ry * np.sin(Th) * np.sin(Ph) + p[1]
    z = Ry * np.cos(Th) + p[2]

    Npoints = Nt * Np
    Nfaces = (Nt - 1) * Np
    Nedge = 0
    #
    # Colormap
    #
    colmap = get_cmap()
    Ncol = colmap.N
    cmap = colmap(np.arange(Ncol))
    g = round((R - minr) * (Ncol - 1) / (maxr - minr))

    _filename = racine + str(1000 + f)[1:] + '.off'
    filename = pyu.getlong(_filename, pstruc['DIRGEOM'])
    fd = open(filename, 'w')
    fd.write('COFF\n')
    chaine = str(Npoints) + ' ' + str(Nfaces) + ' ' + str(Nedge) + '\n'
    fd.write(chaine)

    for ii in range(Nt):
        for jj in range(Np):
            cpos = str(x[ii, jj]) + ' ' + str(y[ii, jj]) + ' ' + str(z[ii, jj])
            cpos = cpos.replace(',', '.')
            ik = g[ii, jj]
            ccol = str(cmap[ik, 0]) + ' ' + str(cmap[ik, 1]) + \
                ' ' + str(cmap[ik, 2])
            ccol = ccol.replace(',', '.')
            fd.write(cpos + ' ' + ccol + ' 0.8\n')

    for ii in range(Nt - 1):
        for jj in range(Np):
            p1 = ii * Np + jj
            p2 = ii * Np + mod(jj + 1, Np)
            p3 = (ii + 1) * Np + jj
            p4 = (ii + 1) * Np + mod(jj + 1, Np)
            chaine = '4 ' + str(p1) + ' ' + str(p2) + ' ' + \
                str(p4) + ' ' + str(p3) + ' 0.5\n'
            fd.write(chaine)

    fd.close()


class SHCoeff(object):
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
        """ init SHCoeff

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
            print('SHCoeff : M>N ')
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

        shape 1   array [ Nf x (N+1) x (M+1) ]
        shape 2   array [ Nf x (N+1)*(M+1)   ]

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
            print('error SHCoeff s1tos2: N2>N1')

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

    def info(self):
        """ info about SHCoeff
        """

        print "Nf   : ", self.Nf
        print "fmin (GHz) : ", self.fmin
        print "fmax (GHz) : ", self.fmax

        sh1 = np.shape(self.s1)
        sh2 = np.shape(self.s2)
        sh3 = np.shape(self.s3)

        if sh1[0] != 0:
            print "N1  : ", self.N1
            print "M1  : ", self.M1
            print "Ncoeff s1 ", self.M1 * self.N1

        if sh2[0] != 0:
            print "NCoeff s2  : ", len(self.ind2)

        if sh3[0] != 0:
            print "Ncoeff s3 : ", len(self.ind3)

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
                plt.xlabel('index', fontsize=26)
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

    Br = SHCoeff(br)
    Bi = SHCoeff(bi)
    Cr = SHCoeff(cr)
    Ci = SHCoeff(ci)
    C  = VSHCoeff(Br,Bi,Cr,Ci)
    """
    def __init__(self, Br, Bi, Cr, Ci):
        """
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

    def info(self):
        """ VSH information
        """
        print "Br"
        print "-------------"
        self.Br.info()
        print "Bi"
        print "-------------"
        self.Bi.info()
        print "Cr"
        print "-------------"
        self.Cr.info()
        print "Ci"
        print "-------------"
        self.Ci.info()

    def show(self, typ='s1', k=1, N=-1, M=-1, kmax = 1000, seuildb=50, animate=False):
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
        plt.figure()
        if not animate:
            plt.subplot(221)
            titre = '$|Br_{n}^{(m)}|$'
            self.Br.show(typ, k, N, M, kmax, seuildb, titre, xl=False, yl=True)
            plt.subplot(222)
            titre = '$|Bi_{n}^{(m)}|$'
            self.Bi.show(typ, k, N, M, kmax, seuildb, titre, xl=False, yl=False)
            plt.subplot(223)
            titre = '$|Cr_{n}^{(m)}|$'
            self.Cr.show(typ, k, N, M, kmax, seuildb, titre, xl=True, yl=True)
            plt.subplot(224)
            titre = '$|Ci_{n}^{(m)}|$'
            self.Ci.show(typ, k, N, M, kmax, seuildb, titre, xl=True, yl=False)
        else:
            for k in np.arange(self.Br.Nf):
                plt.subplot(221)
                titre = '$|Br_{n}^{(m)}|$'
                self.Br.show(typ, k, N, M, kmax, seuildb, titre, xl=False, yl=True)
                plt.subplot(222)
                titre = '$|Bi_{n}^{(m)}|$'
                self.Bi.show(typ, k, N, M, kmax, seuildb, titre, xl=False, yl=False)
                plt.subplot(223)
                titre = '$|Cr_{n}^{(m)}|$'
                self.Cr.show(typ, k, N, M, kmax, seuildb, titre, xl=True, yl=True)
                plt.subplot(224)
                titre = '$|Ci_{n}^{(m)}|$'
                self.Ci.show(typ, k, N, M, kmax, seuildb, titre, xl=True, yl=False)
    #    show()

    def s1tos2(self, N2=-1):
        """ convert shape 1 to shape 2

        Parameters
        ----------
        N2 : max level
            default (-1 means all values)

        """
        self.Bi.s1tos2(N2)
        self.Br.s1tos2(N2)
        self.Ci.s1tos2(N2)
        self.Cr.s1tos2(N2)

    def s2tos3(self, threshold=1e-20):
        """ convert vector spherical coefficient from shape 2 to shape 3

        Parameters
        ----------

        threshold : float
            default 1e-20


        Energy thresholded coefficients
        This thresholding doesn't work

        """

        EBr = np.sum(np.abs(self.Br.s2) ** 2, axis=0)
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


class Antenna(object):
    """ Class Antenna

    Attributes
    ----------

    name   : Antenna name

    Nf     : number of frequency
    Nt     : number of theta
    Np     : number of phi

    Ftheta : Normalized Ftheta    (nf,ntheta,nphi)
    Fphi   : Normalized Fphi      (nf,ntheta,nphi)
    SqG    : Square root of gain  (nf,ntheta,nphi)

    Ttheta : theta                (nf,ntheta,nphi)
    Tphi   : phi                  (nf,ntheta,nphi)

    theta  : theta base            1 x ntheta
    phi    : phi base              1 x phi

    C      : VSH Coefficients

    Methods
    -------

    info()          : Display information about antenna
    vsh()           : calculates Vector Spherical Harmonics
    show3_geom()    : Geomview diagram
    show3()         : 3D diagram plotting using matplotlib toolkit
    Fsynth2         : Antenna F synthesis from coeff in s2

    Antenna trx file can be stored in various order
        natural : HFSS
        bcp     : near filed chamber
    It is important when initializing an antenna object to be aware of the typ of trx file


    .trx (ASCII Vectorial antenna Pattern)

    F   Phi   Theta  Fphi  Ftheta

    """
    def __init__(self, typ, _filename='', directory="ant", nf=104, ntheta=181, nphi=90):
        """

        Parameters
        ----------
        typ  : str
            type of file to read antenna {'mat','vsh2','vsh3','trx','trx1'}
        _filename : str
            antenna file name
        directory : str
            subdirectory of the current project where to find the antenna file
            the file is seek in the $PyProject/ant directory
        nf    :  int
            number of frequency (default 104)
        ntheta:   int
            number of theta (default 181)
        nph : int
            number of phi (default 90)

        Ant = Antenna(typ,_filename)

        if typ=='vsh3':
            Ant.loadvsh3()
        if typ=='vsh2':
            Ant.loadvsh2()
        if typ=='trx':
            Ant.loadtrx()


        """
        self.typ = typ
        self._filename = _filename
        if typ == 'vsh3':
            self.loadvsh3()
        if typ == 'vsh2':
            self.loadvsh2()
        if typ == 'trx':
            self.load_trx(directory, nf, ntheta, nphi)
        if typ == 'trx1':
            self.loadtrx(directory)
        if typ == 'mat':
            self.loadmat(directory)

    def loadmat(self, directory="ant"):
        """ load an antenna stored in a mat file

        Parameters
        ----------
        directory : str , optional
            default 'ant'

        Examples
        --------

            Read an Antenna file in UWBAN directory and plot a polar plot

        .. plot::
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> from pylayers.antprop.antenna import *
            >>> A = Antenna('mat','S1R1.mat','ant/UWBAN/Matfile')
            >>> pol1 = plt.polar(A.phi,abs(A.Ftheta[10,45,:]),'b')
            >>> pol2 = plt.polar(A.phi,abs(A.Ftheta[20,45,:]),'r')
            >>> pol3 = plt.polar(A.phi,abs(A.Ftheta[30,45,:]),'g')
            >>> txt = plt.title('S1R1 antenna : st loadmat')
            >>> plt.show()


        """
        _filemat = self._filename
        filemat = pyu.getlong(_filemat, directory)
        d = io.loadmat(filemat, squeeze_me=True, struct_as_record=False)
        ext = _filemat.replace('.mat', '')
        d = d[ext]
        #
        #
        #
        self.typ = 'mat'
        self.Date = str(d.Date)
        self.Notes = str(d.Notes)
        self.PhotoFile = str(d.PhotoFile)
        self.Serie = eval(str(d.Serie))
        self.Run = eval(str(d.Run))
        self.DataFile = str(d.DataFile)
        self.StartTime = str(d.StartTime)
        self.AntennaName = str(d.AntennaName)

        self.fa = d.freq/1.e9
        self.theta = d.theta
        self.phi = d.phi
        self.Ftheta = d.Ftheta
        self.Fphi = d.Fphi
        self.Nt = len(self.theta)
        self.Np = len(self.phi)
        self.Nf = len(self.fa)

    def load_trx(self, directory="ant", nf=104, ntheta=181, nphi=90, ncol=6):
        """ load a trx file

        Parameters
        ----------
        directory : str
            directory where to find the file

        """
        _filetrx = self._filename
        filename = pyu.getlong(_filetrx, directory)
        if ncol == 6:
            pattern = """^.*\t.*\t.*\t.*\t.*\t.*\t.*$"""
        else:
            pattern = """^.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*$"""
        fd = open(filename, 'r')
        d = fd.read().split('\r\n')
        fd.close()
        k = 0
        #while ((re.search(pattern1,d[k]) is None ) & (re.search(pattern2,d[k]) is None )):
        while re.search(pattern, d[k]) is None:
            k = k + 1

        d = d[k:]
        N = len(d)
        del d[N - 1]
        r = '\t'.join(d)
        r.replace(' ', '')
        d = np.array(r.split()).astype('float')

        #
        # TODO Parsing the header
        #
        #nf     = 104
        #nphi   = 90
        #ntheta = 181

        N = nf * nphi * ntheta
        d = d.reshape(N, 7)

        F = d[:, 0]
        PHI = d[:, 1]
        THETA = d[:, 2]
        Fphi = d[:, 3] + d[:, 4] * 1j
        Ftheta = d[:, 5] + d[:, 6] * 1j

        self.Fphi = Fphi.reshape((nf, nphi, ntheta))
        self.Ftheta = Ftheta.reshape((nf, nphi, ntheta))
        Ttheta = THETA.reshape((nf, nphi, ntheta))
        Tphi = PHI.reshape((nf, nphi, ntheta))
        Tf = F.reshape((nf, nphi, ntheta))

        self.Fphi = self.Fphi.swapaxes(1, 2)
        self.Ftheta = self.Ftheta.swapaxes(1, 2)
        Ttheta = Ttheta.swapaxes(1, 2)
        Tphi = Tphi.swapaxes(1, 2)
        Tf = Tf.swapaxes(1, 2)

        self.fa = Tf[:, 0, 0]
        self.theta = Ttheta[0, :, 0]
        #self.phi     = Tphi[0,0,:]

        #
        # Temporaire
        #
        A1 = self.Fphi[:, 90:181, :]
        A2 = self.Fphi[:, 0:91, :]
        self.Fphi = np.concatenate((A1, A2[:, ::-1, :]), axis=2)
        A1 = self.Ftheta[:, 90:181, :]
        A2 = self.Ftheta[:, 0:91, :]
        self.Ftheta = np.concatenate((A1, A2[:, ::-1, :]), axis=2)
        self.theta = np.linspace(0, np.pi, 91)
        self.phi = np.linspace(0, 2 * np.pi, 180, endpoint=False)
        self.Nt = 91
        self.Np = 180
        self.Nf = 104

    def errel(self, lmax, kf, dsf, typ='s3'):
        """ calculates error between antenna pattern and reference pattern

        This function works for a single frequency point

        Parameters
        ----------

        lmax : integer
            maximum order
        kf  : integer
            frequency index
        dsf : down sampling factor

        Returns
        -------

        errelTh : float
            relative error on :math:`F_{\\theta}`
        errelPh : float
            relative error on :math:`F_{\phi}`
        errel   : float

        Notes
        -----

        .. math::

            \epsilon_r^{\\theta} =
            \\frac{|F_{\\theta}(\\theta,\phi)-\hat{F}_{\\theta}(\\theta)(\phi)|^2}
                 {|F_{\\theta}(\\theta,\phi)|^2}

            \epsilon_r^{\phi} =
            \\frac{|F_{\phi}(\\theta,\phi)-\hat{F}_{\phi}(\\theta)(\phi)|^2}
                 {|F_{\\theta}(\\theta,\phi)|^2}


        """
        #
        # Retrieve angular bases from the down sampling factor dsf
        #
        theta = self.theta[::dsf]
        phi = self.phi[::dsf]
        Nt = len(theta)
        Np = len(phi)

        Th = np.kron(theta, np.ones(Np))
        Ph = np.kron(np.ones(Nt), phi)

        if typ =='s1':
            Fth, Fph = self.Fsynth1(Th, Ph)
        if typ =='s2':
            Fth, Fph = self.Fsynth2b(Th, Ph)
        if typ =='s3':
            Fth, Fph = self.Fsynth3(Th, Ph)

        FTh = Fth.reshape(self.Nf, Nt, Np)
        FPh = Fph.reshape(self.Nf, Nt, Np)
        #
        #  Jacobian
        #
        #st    = outer(sin(theta),ones(len(phi)))
        st = np.sin(theta).reshape((len(theta), 1))
        #
        # Construct difference between reference and reconstructed
        #
        dTh = (FTh[kf, :, :] - self.Ftheta[kf, ::dsf, ::dsf])
        dPh = (FPh[kf, :, :] - self.Fphi[kf, ::dsf, ::dsf])
        #
        # squaring  + Jacobian
        #
        dTh2 = np.real(dTh * np.conj(dTh)) * st
        dPh2 = np.real(dPh * np.conj(dPh)) * st

        vTh2 = np.real(self.Ftheta[kf, ::dsf, ::dsf] * np.conj(
            self.Ftheta[kf, ::dsf, ::dsf])) * st
        vPh2 = np.real(self.Fphi[kf, ::dsf, ::dsf] * np.conj(
            self.Fphi[kf, ::dsf, ::dsf])) * st

        mvTh2 = np.sum(vTh2)
        mvPh2 = np.sum(vPh2)

        errTh = np.sum(dTh2)
        errPh = np.sum(dPh2)

        errelTh = errTh / mvTh2
        errelPh = errPh / mvPh2
        errel = (errTh + errPh) / (mvTh2 + mvPh2)

        return(errelTh, errelPh, errel)

    def loadtrx(self, directory):
        """
        load trx file

        self._filename: short name of the antenna file

        the file is seek in the $PyProject/ant directory


        """

        _filetrx = self._filename
        _headtrx = 'header_' + _filetrx
        _headtrx = _headtrx.replace('trx', 'txt')
        headtrx = pyu.getlong(_headtrx, directory)
        filename = pyu.getlong(_filetrx, directory)
        #
        # Header
        #
        foh = open(headtrx)
        ligh = foh.read()
        foh.close()
        nf = eval(ligh.split()[2])
        nphi = eval(ligh.split()[5])
        ntheta = eval(ligh.split()[8])
        try:
            tau = eval(ligh.split()[9])  # tau : delay (ns)
        except:
            tau = 0

        #
        # Data
        #
        fi = open(filename)
        d = np.array(fi.read().split())
        N = len(d)
        M = N / 7
        d = d.reshape(M, 7)
        d = d.astype('float')

        f = d[:, 0]
        if f[0] == 0:
            print "error : frequency cannot be zero"
        if (f[0] > 200):
            f = f / 1.0e9

        phi = d[:, 1]
        theta = d[:, 2]
        #
        # type : refers to the way the angular values are stored in the file
        # Detection of file type
        #
        # Bcp
        #    f  phi theta
        #    2    1    0
        # Natural
        #    f  phi theta
        #    2    0    1

        dphi = abs(phi[0] - phi[1])
        dtheta = abs(theta[0] - theta[1])

        if (dphi == 0) & (dtheta != 0):
            typ = 'bcp'
        if (dtheta == 0) & (dphi != 0):
            typ = 'natural'

        self.typ = typ
        Fphi = d[:, 3] + d[:, 4] * 1j
        Ftheta = d[:, 5] + d[:, 6] * 1j
        #
        # Normalization
        #
        G = np.real(Fphi * np.conj(Fphi) + Ftheta * np.conj(Ftheta))
        SqG = np.sqrt(G)
        #Fphi    = Fphi/SqG
        #Ftheta  = Ftheta/SqG
        #Fphi    = Fphi
        #Ftheta  = Ftheta
        #
        # Reshaping
        #
        if typ == 'natural':
            self.Fphi = Fphi.reshape((nf, ntheta, nphi))
            self.Ftheta = Ftheta.reshape((nf, ntheta, nphi))
            self.SqG = SqG.reshape((nf, ntheta, nphi))
            self.Ttheta = theta.reshape((nf, ntheta, nphi))
            self.Tphi = phi.reshape((nf, ntheta, nphi))
            self.Tf = f.reshape((nf, ntheta, nphi))
        if typ == 'bcp':
            self.Fphi = Fphi.reshape((nf, nphi, ntheta))
            self.Ftheta = Ftheta.reshape((nf, nphi, ntheta))
            self.SqG = SqG.reshape((nf, nphi, ntheta))
            self.Ttheta = theta.reshape((nf, nphi, ntheta))
            self.Tphi = phi.reshape((nf, nphi, ntheta))
            self.Tf = f.reshape((nf, nphi, ntheta))
        #
        # Force natural order (f,theta,phi)
        # This is not the order of the satimo nfc which is  (f,phi,theta)
        #
        #  .. todo::  May be it is not necessary to load Ttheta and Tphi
        #

            self.Fphi = self.Fphi.swapaxes(1, 2)
            self.Ftheta = self.Ftheta.swapaxes(1, 2)
            self.SqG = self.SqG.swapaxes(1, 2)
            self.Ttheta = self.Ttheta.swapaxes(1, 2)
            self.Tphi = self.Tphi.swapaxes(1, 2)
            self.Tf = self.Tf.swapaxes(1, 2)

        self.fa = self.Tf[:, 0, 0]
        self.theta = self.Ttheta[0, :, 0]
        self.phi = self.Tphi[0, 0, :]

        self.Nf = nf
        self.Nt = ntheta
        self.Np = nphi
        self.tau = tau

    def checkpole(self, kf=0):
        """ display the reconstructed field on pole for integrity verification

        Parameters
        ----------
        kf : int
            frequency index default 0

        """
        Ft0 = self.Ftheta[kf, 0, :]
        Fp0 = self.Fphi[kf, 0, :]
        Ftp = self.Ftheta[kf, -1, :]
        Fpp = self.Fphi[kf, -1, :]
        phi = self.phi

        Ex0 = Ft0 * np.cos(phi) - Fp0 * np.sin(phi)
        Ey0 = Ft0 * np.sin(phi) + Fp0 * np.cos(phi)
        Exp = Ftp * np.cos(phi) - Fpp * np.sin(phi)
        Eyp = Ftp * np.sin(phi) + Fpp * np.cos(phi)

        plt.subplot(4, 2, 1)
        plt.plot(phi, np.real(Ex0))
        plt.subplot(4, 2, 2)
        plt.plot(phi, np.imag(Ex0))
        plt.subplot(4, 2, 3)
        plt.plot(phi, np.real(Ey0))
        plt.subplot(4, 2, 4)
        plt.plot(phi, np.imag(Ey0))
        plt.subplot(4, 2, 5)
        plt.plot(phi, np.real(Exp))
        plt.subplot(4, 2, 6)
        plt.plot(phi, np.imag(Exp))
        plt.subplot(4, 2, 7)
        plt.plot(phi, np.real(Eyp))
        plt.subplot(4, 2, 8)
        plt.plot(phi, np.imag(Eyp))

    def info(self):
        """ gives info about antenna object

           >>> A1 = Antenna('trx1','defant.trx')
           >>> A2 = Antenna('vsh3','defant.vsh3')
           >>> A3 = Antenna('mat','S1R1.mat','ant/UWBAN/Matfile')

        """
        print self._filename
        print "type : ", self.typ
        if self.typ == 'mat':
            print self.DataFile
            print self.AntennaName
            print self.Date
            print self.StartTime
            print self.Notes
            print self.Serie
            print self.Run
            print "Nb theta (lat) :", self.Nt
            print "Nb phi (lon) :", self.Np
        print "--------------------------"
        print "fmin (GHz) :", self.fa[0]
        print "fmax (GHz) :", self.fa[-1]
        print "Nf   :", len(self.fa)
        try:
            self.C.info()
        except:
            print "No vsh coefficient calculated yet"

    def help(self):
        """
        help on available functions
        """
        print "Antenna Class : List of available function"
        print "------------------------------------------"
        print "A.init(filename)"
        print "A.info()"
        print "A.demo()"
        print "A.help()"
        print "A.show3_geom(k=0,typ='Gain'|'Ftheta'|Fphi',mode='linear'|'not implemented',silent=(True)|False)"
        print "A.pol3d(k=0,R=1,St=1(Downsampling along theta),Sp=1(idem),silent=(False)[True)"
        print "A.mse(Fth,Fph,N=0) "
        print "A.elec_delay()"
        print "A.vsh() : Vector spherical harmonics transform"
        print "A.Fsynth2(theta,phi) : Calulate diagram from shape 2 coeff"
        print "A.Fsynth3(theta,phi) : Calulate diagram from shape 3 coeff"
        print "movie_vsh(mode='linear')"

    def polar(self, k=[0], it=0, ip=-1, dyn=6, GmaxdB=20, alpha=0.1):
        """ polar plot

            Parameters
            ----------
            k : list of int
                frequency index  (default 0)
            it  : int
                theta index      (default 0)
            ip  : int
                phi index        (default -1)
            GmaxdB :
                Max Gain (dB)
            dyn    :
                dynamic number of 5dB step
            alpha  : float
                default 0.1

            Examples
            --------

            .. plot::
                :include-source:

                >>> import matplotlib.pyplot as plt
                >>> from pylayers.antprop.antenna import *
                >>> A = Antenna('trx1','defant.trx')
                >>> A.polar(k=[0,10,50])
                >>> plt.show()

        """

        rc('grid', color='#316931', linewidth=1, linestyle='-')
        rc('xtick', labelsize=15)
        rc('ytick', labelsize=15)
        DyndB = dyn * 5
        GmindB = GmaxdB - DyndB
        # force square figure and square axes looks better for polar, IMO
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='#d5de9c')
        t1 = np.arange(5, DyndB + 5, 5)
        t2 = np.arange(GmindB + 5, GmaxdB + 5, 5)
        rline1, rtext1 = plt.rgrids(t1, t2)
        a1 = np.arange(0, 360, 30)
        a2 = [90, 60, 30, 0, 330, 300, 270, 240, 210, 180, 150, 120]
        rline2, rtext2 = plt.thetagrids(a1, a2)

        col = ['k', 'r', 'g', 'b', 'm', 'c', 'y']
        cpt = 0
        for ik in k:
            chaine = 'f = ' + str(self.fa[ik]) + ' GHz'
            if it == -1:
                itheta = np.arange(self.Nt)
                iphi1 = ip
                Np = self.Np
                if mod(Np, 2) == 0:
                    iphi2 = mod(ip + Np / 2, Np)
                else:
                    iphi2 = mod(ip + (Np - 1) / 2, Np)

                u1 = nonzero((self.theta <= np.pi / 2) & (self.theta >= 0))
                u2 = np.arange(self.Nt)
                u3 = nonzero((self.theta <= np.pi) & (
                    self.theta > np.pi / 2))
                r1 = -GmindB + 20 * np.log10(
                    alpha * self.SqG[ik, u1[0], iphi1])
                #r1  = self.SqG[k,u1[0],iphi1]
                negr1 = nonzero(r1 < 0)
                r1[negr1[0]] = 0
                r2 = -GmindB + 20 * np.log10(alpha * self.SqG[ik, u2, iphi2])
                #r2  = self.SqG[k,u2,iphi2]
                negr2 = nonzero(r2 < 0)
                r2[negr2[0]] = 0
                r3 = -GmindB + 20 * np.log10(
                    alpha * self.SqG[ik, u3[0], iphi1])
                #r3  = self.SqG[k,u3[0],iphi1]
                negr3 = nonzero(r3 < 0)
                r3[negr3[0]] = 0
                r = np.hstack((r1[::-1], r2, r3[::-1], r1[-1]))

                angle = np.linspace(0, 2 * np.pi, len(r), endpoint=True)
            else:
                iphi = np.arange(self.Np)
                itheta = it
                angle = self.phi[iphi]
                r = -GmindB + 20 * np.log10(self.SqG[ik, itheta, iphi])

            ax.plot(angle, r, color=col[cpt], lw=2, label=chaine)
            cpt = cpt + 1
        ax.legend()

    def show3_geom(self, k=0, typ='Gain', mode='linear', silent=False):
        """ show3 geomview

        Parameters
        ----------

        k : frequency index
        typ   = 'Gain' | 'Ftheta' | 'Fphi'
        mode   = 'linear'| 'not implemented'
        silent = True    | False
        """

        nt = self.Nt
        np = self.Np
        th = self.theta
        ph = self.phi
        f = self.fa[k]

        if typ == 'Gain':
            V = self.SqG[k, :, :]
        if typ == 'Ftheta':
            V = self.Ftheta[k, :, :]
        if typ == 'Fphi':
            V = self.Fphi[k, :, :]

        minr = abs(V).min()
        maxr = abs(V).max()

        po = np.array([0, 0, 0])

        geom_pattern(th, ph, V, k, po, minr, maxr, typ)

        _filename = typ + str(k) + '.off'
        filename = pyu.getlong(_filename, pstruc['DIRGEOM'])

        if not silent:
            chaine = "geomview -nopanel -b 1 1 1 " + filename + \
                " 2>/dev/null &"
            os.system(chaine)

    def show3(self, k=0, typ='Gain', col=True):
        """
        show3(self,k=0,typ='Gain',col=True)

        k : frequency index

        typ = 'Gain'
             = 'Ftheta'
             = 'Fphi'

        if col  -> color coded plot3D
        else    -> simple plot3D
        """

        nt = self.Nt
        np = self.Np
        th = self.theta
        ph = self.phi

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)

        if typ == 'Gain':
            V = self.SqG[k, :, :]
        if typ == 'Ftheta':
            V = self.Ftheta[k, :, :]
        if typ == 'Fphi':
            V = self.Fphi[k, :, :]

        vt = np.ones(nt)
        vp = np.ones(np)
        Th = np.outer(th, vp)
        Ph = np.outer(vt, ph)

        X = abs(V) * np.cos(Ph) * np.sin(Th)
        Y = abs(V) * np.sin(Ph) * np.sin(Th)
        Z = abs(V) * np.cos(Th)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        if col:
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        else:
            ax.plot3D(ravel(X), ravel(Y), ravel(Z))

    def pol3d(self, k=0, R=50, St=1, Sp=1, silent=False):
        """ Display polarisation diagram  in 3D

           pol3d(k=0,R=1,St=1,Sp=1,silent=False):

           Parameters
           ----------
           k  : int
               frequency index
           R  : float
                Radius of the sphere
           St : int
               Downsampling factor along theta
           Sp : int
               Downsampling factor along phi
           silent : Boolean
                (If True the file is created and not displayed')

           The file created is named : Polar{ifreq}.list
           it is placed in the /geom directory of the project

        """
        _filename = 'Polar' + str(10000 + k)[1:] + '.list'
        filename = pyu.getlong(_filename, pstruc['DIRGEOM'])
        fd = open(filename, "w")
        fd.write("LIST\n")
        Nt = self.Nt
        Np = self.Np
        N = 10
        plth = np.arange(0, Nt, St)
        plph = np.arange(0, Np, Sp)
        for m in plph:
            for n in plth:
                theta = self.theta[n]
                #print "m,theta= :",m,theta*180/np.pi
                phi = self.phi[m]
                #print "n,phi=:",n,phi*180/np.pi
                B = vec_sph(theta, phi)
                p = R * np.array((np.cos(phi) * np.sin(theta),
                                  np.sin(phi) * np.sin(theta),
                                  np.cos(theta)))
                fd.write('{\n')
                ellipse(fd, p, B[0, :], B[1, :], self.Ftheta[
                    k, n, m], self.Fphi[k, n, m], N)
                fd.write('}\n')
        fd.close()
        if not silent:
            chaine = "geomview " + filename + " 2>/dev/null &"
            os.system(chaine)

    def mse(self, Fth, Fph, N=0):
        """ mean square error

        Parameters
        ----------
        Fth  : np.array
        Fph  : np.array
        N    : int

        Notes
        -----

        Calculate the relative mean square error between original pattern A.Ftheta , A.Fphi and the
        pattern given as argument of the function  Fth , Fph

        The mse is evaluated on both polarization and normalized over the energy of each
        original pattern.

        The function returns the maximum between those two errors

        N is a parameter which allows to suppress value at the pole for the calculation of the error
        if N=0 all values are kept else   N < n < Nt - N

        """

        sh = np.shape(self.Ftheta)
        Nf = sh[0]
        Nt = sh[1]
        Np = sh[2]
        # plage de theta (exclusion du pole)
        pt = np.arange(N, Nt - N, 1)

        Fthr = Fth.reshape(sh)
        Fphr = Fph.reshape(sh)

        Gr = np.real(Fphr * np.conj(Fphr) + Fthr * np.conj(Fthr))
        SqGr = np.sqrt(Gr)

        Fthr = Fthr[:, pt, :].ravel()
        Fphr = Fphr[:, pt, :].ravel()
        SqGr = SqGr[:, pt, :].ravel()

        Ftho = self.Ftheta[:, pt, :].ravel()
        Fpho = self.Fphi[:, pt, :].ravel()
        SqGo = self.SqG[:, pt, :].ravel()

        Etho = np.sqrt(np.dot(np.conj(Ftho), Ftho))
        Epho = np.sqrt(np.dot(np.conj(Fpho), Fpho))
        Eo = np.sqrt(np.dot(np.conj(Ftho), Ftho) + np.dot(np.conj(Fpho), Fpho))

        errth = Ftho - Fthr
        errph = Fpho - Fphr
        Err = np.real(np.sqrt(np.dot(np.conj(errth), errth) + np.dot(np.conj(errph), errph)))

        Errth = np.real(np.sqrt(np.dot(np.conj(errth), errth)))
        Errph = np.real(np.sqrt(np.dot(np.conj(errph), errph)))

        #Errth_rel = Errth/Etho
        #Errph_rel = Errph/Epho
        Errth_rel = Errth / Eo
        Errph_rel = Errph / Eo
        Err_rel = Err / Eo

        return Err_rel, Errth_rel, Errph_rel

    def getdelay(self,freq,delayCandidates = np.arange(-10,10,0.001)):
        """ getelectrical delay

        Parameters
        ----------
        delayCandidates : ndarray
            default np.arange(-10,10,0.001)

        Returns
        -------
        electricalDelay  : float

        Author : Troels Pedersen (Aalborg University)
                 B.Uguen
        """
        maxPowerInd  = np.unravel_index(np.argmax(abs(self.Ftheta)),np.shape(self.Ftheta))
        electricalDelay  = delayCandidates[np.argmax(abs(
            np.dot(self.Ftheta[:,maxPowerInd[1],maxPowerInd[2]]
                ,np.exp(2j*np.pi*freq.reshape(len(freq),1)
                  *delayCandidates.reshape(1,len(delayCandidates))))
                ))]
        return(electricalDelay)

    def elec_delay(self):
        """ apply an electrical delay

         Notes
         -----
         This function apply an electrical delay math::`\\exp{2 j \\pi f \\tau)`
         on the phase of diagram math::``F_{\\theta}`` and math::`F_{\\phi}`

        """
        tau = self.tau
        Ftheta = self.Ftheta
        Fphi = self.Fphi
        sh = np.shape(Ftheta)
        e = np.exp(2 * pi * 1j * self.fa * tau)
        E = np.outer(e, ones(sh[1] * sh[2]))

        Fth = Ftheta.reshape(sh[0], sh[1] * sh[2])
        EFth = Fth * E
        self.Ftheta = EFth.reshape(sh[0], sh[1], sh[2])

        Fph = Fphi.reshape(sh[0], sh[1] * sh[2])
        EFph = Fph * E
        self.Fphi = EFph.reshape(sh[0], sh[1], sh[2])

    def vshd(self, dsf=1):
        """

        Parameters
        ----------
        dsf :  int
            down sampling factor  'default 1'

        Summary
        -------

        This function calculates the Vector Spherical Harmonics coefficients
        It makes use of the spherepack function vha

            m : phi    longitude
            n : theta  latitude

        Antenna pattern are stored       (f theta phi)
        Coeff are stored with this order (f , n , m )

        The vsh coefficient are organized differently
        should be better for compression along frequency axis


        """

        th = self.theta[::dsf]
        ph = self.phi[::dsf]

        nth = len(th)
        nph = len(ph)
        nf = self.Nf

        if (nph % 2) == 1:
            mdab = min(nth, (nph + 1) / 2)
        else:
            mdab = min(nth, nph / 2)

        ndab = nth

        Br = 1j * np.zeros((nf, ndab, mdab))
        Bi = 1j * np.zeros((nf, ndab, mdab))
        Cr = 1j * np.zeros((nf, ndab, mdab))
        Ci = 1j * np.zeros((nf, ndab, mdab))

        gridComp = Wrapec()
        wvha, lvha = gridComp.vhai(nth, nph)

        for k in range(nf):
            #
            # Real part
            #
            Fpr = self.Fphi[k][::dsf, ::dsf].real
            Ftr = self.Ftheta[k][::dsf, ::dsf].real
            #
            # Fpr     Ntheta,Nphi
            #
            brr, bir, crr, cir = gridComp.vha(nth, nph, 1,
                                              lvha, wvha,
                                              np.transpose(Fpr),
                                              np.transpose(Ftr))
            #
            # Imaginary part
            #
            Fpi = self.Fphi[k][::dsf, ::dsf].imag
            Fti = self.Ftheta[k][::dsf, ::dsf].imag
            bri, bii, cri, cii = gridComp.vha(nth, nph, 1,
                                              lvha, wvha,
                                              np.transpose(Fpi),
                                              np.transpose(Fti))

            Br[k, :, :] = brr + 1j * bri
            Bi[k, :, :] = bir + 1j * bii
            Cr[k, :, :] = crr + 1j * cri
            Ci[k, :, :] = cir + 1j * cii

        #
        # m=0 row is multiplied by 0.5
        #

        Br[:, :, 0] = 0.5 * Br[:, :, 0]
        Bi[:, :, 0] = 0.5 * Bi[:, :, 0]
        Cr[:, :, 0] = 0.5 * Cr[:, :, 0]
        Ci[:, :, 0] = 0.5 * Ci[:, :, 0]

        #print "self.fa[0] = ",self.fa[0]
        #print "self.fa[-1] = ",self.fa[-1]

        Br = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Br)
        Bi = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Bi)
        Cr = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Cr)
        Ci = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Ci)

        self.C = VSHCoeff(Br, Bi, Cr, Ci)

    def vsh(self):
        """ calculates the Vector Spherical Harmonics coefficients

        Summary
        -------

        It makes use of the spherepack function vha

            m : phi    longitude
            n : theta  latitude

        Antenna pattern are stored       (f theta phi)
        Coeff are stored with this order (f , n , m )

        The vsh coefficient are organized differently
        should be better for compression along frequency axis

        """

        nt = self.Nt
        np = self.Np
        nf = self.Nf

        if np % 2:
            mdab = min(nt, (np + 1) / 2)
        else:
            mdab = min(nt, np / 2)

        ndab = nt

        Br = 1j * np.zeros((nf, ndab, mdab))
        Bi = 1j * np.zeros((nf, ndab, mdab))
        Cr = 1j * np.zeros((nf, ndab, mdab))
        Ci = 1j * np.zeros((nf, ndab, mdab))

        gridComp = Wrapec()
        wvha, lvha = gridComp.vhai(nt, np)

        for k in range(nf):
            #
            # Real part
            #
            Fpr = self.Fphi[k].real
            Ftr = self.Ftheta[k].real
            #
            # Fpr     Ntheta,Nphi
            #
            brr, bir, crr, cir = gridComp.vha(nt, np, 1, lvha,
                                              wvha, transpose(Fpr), transpose(Ftr))
            #
            # Imaginary part
            #
            Fpi = self.Fphi[k].imag
            Fti = self.Ftheta[k].imag
            bri, bii, cri, cii = gridComp.vha(nt, np, 1, lvha,
                                              wvha, transpose(Fpi), transpose(Fti))

            Br[k, :, :] = brr + 1j * bri
            Bi[k, :, :] = bir + 1j * bii
            Cr[k, :, :] = crr + 1j * cri
            Ci[k, :, :] = cir + 1j * cii

        #
        # m=0 row is multiplied by 0.5
        #

        Br[:, :, 0] = 0.5 * Br[:, :, 0]
        Bi[:, :, 0] = 0.5 * Bi[:, :, 0]
        Cr[:, :, 0] = 0.5 * Cr[:, :, 0]
        Ci[:, :, 0] = 0.5 * Ci[:, :, 0]

        #print "self.fa[0] = ",self.fa[0]
        #print "self.fa[-1] = ",self.fa[-1]

        Br = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Br)
        Bi = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Bi)
        Cr = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Cr)
        Ci = SHCoeff(typ='s1', fmin=self.fa[0], fmax=self.fa[-1], data=Ci)

        self.C = VSHCoeff(Br, Bi, Cr, Ci)

    def demo(self):
        """ display few commands for executing little demo
        """
        print "A.C.s1tos2(30)"
        print "Fth , Fph = A.Fsynth2(th,ph)"
        print "FTh = Fth.reshape(A.Nf,A.Nt,A.Np)"
        print "FPh = Fph.reshape(A.Nf,A.Nt,A.Np)"
        print "compdiag(20,A,A.theta,A.phi,FTh,FPh) "

    #def Fsynth1(self, theta, phi, k=0):
    def Fsynth1(self, theta, phi):
        """ calculate complex antenna pattern  from VSH Coefficients (shape 1)

        Parameters
        ----------
        theta  : ndarray (1xNdir)
        phi    : ndarray (1xNdir)
        k      : int
            frequency index

        """

        nray = len(theta)

        #Br = self.C.Br.s1[k, :, :]
        #Bi = self.C.Bi.s1[k, :, :]
        #Cr = self.C.Cr.s1[k, :, :]
        #Ci = self.C.Ci.s1[k, :, :]

        Br = self.C.Br.s1[:, :, :]
        Bi = self.C.Bi.s1[:, :, :]
        Cr = self.C.Cr.s1[:, :, :]
        Ci = self.C.Ci.s1[:, :, :]

        N = self.C.Br.N1
        M = self.C.Br.M1
        #print "N,M",N,M
        #
        # The - sign is necessary to get the good reconstruction
        #     deduced from observation
        #     May be it comes from a different definition of theta in SPHEREPACK
        x = -np.cos(theta)
        Pmm1n, Pmp1n = AFLegendre3(N, M, x)
        ind = index_vsh(N, M)
        n = ind[:, 0]
        m = ind[:, 1]
        V, W = VW(n, m, x, phi, Pmm1n, Pmp1n)
        #
        # broadcasting along frequency axis
        #
        V = np.expand_dims(V,0)
        W = np.expand_dims(V,0)
        #
        #   k : frequency axis 
        #   l : coeff l 
        #   m  
        Fth = np.eisum('klm,kilm->ki',Br,np.real(V.T)) - \
              np.eisum('klm,kilm->ki',Bi,np.imag(V.T)) + \
              np.eisum('klm,kilm->ki',Ci,np.real(W.T)) + \
              np.eisum('klm,kilm->ki',Cr,np.imag(W.T)) 

        Fph = -np.eisum('klm,kilm->ki',Cr,np.real(V.T)) + \
              np.eisum('klm,kilm->ki',Ci,np.imag(V.T)) + \
              np.eisum('klm,kilm->ki',Bi,np.real(W.T)) + \
              np.eisum('klm,kilm->ki',Br,np.imag(W.T))

        #Fth = np.dot(Br, np.real(V.T)) - \
        #    np.dot(Bi, np.imag(V.T)) + \
        #    np.dot(Ci, np.real(W.T)) + \
        #    np.dot(Cr, np.imag(W.T))

        #Fph = -np.dot(Cr, np.real(V.T)) + \
        #    np.dot(Ci, np.imag(V.T)) + \
        #    np.dot(Bi, np.real(W.T)) + \
        #    np.dot(Br, np.imag(W.T))

        return Fth, Fph

    def Fsynth2b(self, theta, phi):
        """  pattern synthesis from shape 2 vsh coefficients

        Parameters
        ----------
        theta
        phi

        Notes
        -----

        Calculate complex antenna pattern from VSH Coefficients (shape 2)
        for the specified directions (theta,phi)
        theta and phi arrays needs to have the same size

        """

        Br = self.C.Br.s2
        Bi = self.C.Bi.s2
        Cr = self.C.Cr.s2
        Ci = self.C.Ci.s2

        N = self.C.Br.N2
        M = self.C.Br.M2

        #print "N,M",N,M
        #
        # The - sign is necessary to get the good reconstruction
        #     deduced from observation
        #     May be it comes from a different definition of theta in SPHEREPACK

        x = -np.cos(theta)

        Pmm1n, Pmp1n = AFLegendre3(N, M, x)
        ind = index_vsh(N, M)

        n = ind[:, 0]
        m = ind[:, 1]

        V, W = VW2(n, m, x, phi, Pmm1n, Pmp1n)


        Fth = np.dot(Br, np.real(V.T)) - np.dot(Bi, np.imag(V.T)) + \
              np.dot(Ci, np.real(W.T)) + np.dot(Cr, np.imag(W.T))
        Fph = -np.dot(Cr, np.real(V.T)) + np.dot(Ci, np.imag(V.T)) + \
              np.dot(Bi, np.real(W.T)) + np.dot(Br, np.imag(W.T))

        return Fth, Fph


    def Fsynth2(self, theta, phi):
        """  pattern synthesis from shape 2 vsh coeff

        Parameters
        ----------
        theta
        phi

        Notes
        -----

        Calculate complex antenna pattern from VSH Coefficients (shape 2)
        for the specified directions (theta,phi)
        theta and phi arrays needs to have the same size

        """

        Br = self.C.Br.s2
        Bi = self.C.Bi.s2
        Cr = self.C.Cr.s2
        Ci = self.C.Ci.s2

        N = self.C.Br.N2
        M = self.C.Br.M2

        #print "N,M",N,M
        #
        # The - sign is necessary to get the good reconstruction
        #     deduced from observation
        #     May be it comes from a different definition of theta in SPHEREPACK
        x = -np.cos(theta)

        Pmm1n, Pmp1n = AFLegendre3(N, M, x)
        ind = index_vsh(N, M)

        n = ind[:, 0]
        m = ind[:, 1]

        V, W = VW(n, m, x, phi, Pmm1n, Pmp1n)


        Fth = np.dot(Br, np.real(V.T)) - np.dot(Bi, np.imag(V.T)) + \
            np.dot(Ci, np.real(W.T)) + np.dot(Cr, np.imag(W.T))
        Fph = -np.dot(Cr, np.real(V.T)) + np.dot(Ci, np.imag(V.T)) + \
            np.dot(Bi, np.real(W.T)) + np.dot(Br, np.imag(W.T))

        return Fth, Fph


    def Fsynth3(self, theta, phi):
        """ synthesis of a complex antenna pattern from VSH coefficients (shape 3)

        Let Ndir be the number of directions

        Parameters
        ----------

        theta : ndarray (1xNdir)
        phi   : ndarray (1xNdir)

        Returns
        -------

        Fth   : ndarray (1xNdir)
        Fph   : ndarray (1xNdir)

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.antprop.antenna import *
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> A = Antenna('vsh3','defant.vsh3')
            >>> theta = np.linspace(0,np.pi,70)
            >>> phi = np.linspace(0,2*np.pi,180)
            >>> th = np.kron(theta,np.ones(len(phi)))
            >>> ph = np.kron(np.ones(len(theta)),phi)
            >>> Fth,Fph = A.Fsynth3(th,ph)

        """

        nray = len(theta)

        Br = self.C.Br.s3
        nBr = self.C.Br.ind3[:, 0]
        mBr = self.C.Br.ind3[:, 1]

        Bi = self.C.Bi.s3
        Cr = self.C.Cr.s3
        Ci = self.C.Ci.s3

        M = mBr.max()
        N = nBr.max()

        x = -np.cos(theta)

        Pmm1n, Pmp1n = AFLegendre3(N, M, x)

        V, W = VW(nBr, mBr, x, phi, Pmm1n, Pmp1n)

        Fth = np.dot(Br, np.real(V.T)) - \
            np.dot(Bi, np.imag(V.T)) + \
            np.dot(Ci, np.real(W.T)) + \
            np.dot(Cr, np.imag(W.T))

        Fph = -np.dot(Cr, np.real(V.T)) + \
            np.dot(Ci, np.imag(V.T)) + \
            np.dot(Bi, np.real(W.T)) + \
            np.dot(Br, np.imag(W.T))

        return Fth, Fph

    def movie_vsh(self, mode='linear'):
        """
            movie_vsh
        """

        Brmin = abs(self.C.Br[:, 0:20, 0:20]).min()
        Brmax = abs(self.C.Br[:, 0:20, 0:20]).max()
        Bimin = abs(self.C.Bi[:, 0:20, 0:20]).min()
        Bimax = abs(self.C.Bi[:, 0:20, 0:20]).max()

        Crmin = abs(self.C.Cr[:, 0:20, 0:20]).min()
        Crmax = abs(self.C.Cr[:, 0:20, 0:20]).max()
        Cimin = abs(self.C.Ci[:, 0:20, 0:20]).min()
        Cimax = abs(self.C.Ci[:, 0:20, 0:20]).max()

        print Brmin, Brmax, Bimin, Bimax, Crmin, Crmax, Cimin, Cimax

        for k in range(self.nf):
            plt.figure()
            stf = ' f=' + str(self.fa[k]) + ' GHz'
            subplot(221)
            pcolor(abs(self.C.Br.s1[k, 0:20, 0:20]),
                   vmin=Brmin, vmax=Brmax, edgecolors='k')
            #xlabel('m',fontsize=12)
            ylabel('n', fontsize=12)
            title('$|Br_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            subplot(222)
            pcolor(abs(self.C.Bi.s1[k, 0:20, 0:20]),
                   vmin=Bimin, vmax=Bimax, edgecolors='k')
            #xlabel('m',fontsize=12)
            ylabel('n', fontsize=12)
            title('$|Bi_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            subplot(223)
            pcolor(abs(self.C.Cr.s1[k, 0:20, 0:20]),
                   vmin=Crmin, vmax=Crmax, edgecolors='k')
            xlabel('m', fontsize=12)
            #ylabel('n',fontsize=12)
            title('$|Cr_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            subplot(224)
            pcolor(abs(self.C.Ci.s1[k, 0:20, 0:20]),
                   vmin=Cimin, vmax=Cimax, edgecolors='k')
            xlabel('m', fontsize=12)
            #ylabel('n',fontsize=12)
            title('$|Ci_{n}^{(m)}|$' + stf, fontsize=10)
            colorbar()
            filename = str('%03d' % k) + '.png'
            savefig(filename, dpi=100)
            clf()

        command = ('mencoder',
                   'mf://*.png',
                   '-mf',
                   'type=png:w=800:h=600:fps=1',
                   '-ovc',
                   'lavc',
                   '-lavcopts',
                   'vcodec=mpeg4',
                   '-oac',
                   'copy',
                   '-o',
                   'vshcoeff.avi')
        subprocess.check_call(command)

    def minsh3(self, emax=0.05):
        """ creates vsh3 with significant coeff until given relative reconstruction error

        Parameters
        ----------
        emax : float
            error default 0.05

        Summary
        -------

        Create antenna's vsh3 file which only contains
        the significant vsh coefficients in shape 3,
        in order to obtain a reconstruction maximal error =  emax

        This function requires a reading of .trx file before being executed

        """

        th = np.kron(self.theta, np.ones(self.Np))
        ph = np.kron(np.ones(self.Nt), self.phi)

        Fth3, Fph3 = self.Fsynth3(th, ph)
        Err = self.mse(Fth3, Fph3, 0)

        Enc = self.C.ens3()
        n = len(Enc)
        pos = 0
        while (pos < n) & (Err[0] < emax):

            Emin = Enc[pos]
            d = self.C.drag3(Emin)
            Fth3, Fph3 = self.Fsynth3(th, ph)
            Err = self.mse(Fth3, Fph3, 0)

            if Err[0] >= emax:
                i = d[0][0]
                i3 = d[1][0]
                self.C.put3(i, i3)
                Fth3, Fph3 = self.Fsynth3(th, ph)
                Err = self.mse(Fth3, Fph3, 0)

            pos = pos + 1

    def savevsh3(self):
        """
        A.savevsh3()

        Create a .vsh3 antenna file

        """

        # create vsh3 file

        _filevsh3 = self._filename.replace('.trx', '.vsh3')
        filevsh3 = pyu.getlong(_filevsh3, pstruc['DIRANT'])

        #filevsh3 = pyu.getlong(self._filename,'ant')

        if os.path.isfile(filevsh3):
            print filevsh3, ' already exist'
        else:
            print 'create ', filevsh3, ' file'

            coeff = {}
            coeff['fmin'] = self.fa[0]
            coeff['fmax'] = self.fa[-1]
            coeff['Br.ind'] = self.C.Br.ind3
            coeff['Bi.ind'] = self.C.Bi.ind3
            coeff['Cr.ind'] = self.C.Cr.ind3
            coeff['Ci.ind'] = self.C.Ci.ind3
            coeff['Br.k'] = self.C.Br.k2
            coeff['Bi.k'] = self.C.Bi.k2
            coeff['Cr.k'] = self.C.Cr.k2
            coeff['Ci.k'] = self.C.Ci.k2
            coeff['Br.s3'] = self.C.Br.s3
            coeff['Bi.s3'] = self.C.Bi.s3
            coeff['Cr.s3'] = self.C.Cr.s3
            coeff['Ci.s3'] = self.C.Ci.s3

            io.savemat(filevsh3, coeff, appendmat=False)

    def loadvsh3(self):
        """ Load antenna's vsh3 file

            it contains a thesholded version of vsh coefficients in shape 3
        """

        _filevsh3 = self._filename
        filevsh3 = pyu.getlong(_filevsh3, pstruc['DIRANT'])

        if os.path.isfile(filevsh3):
            coeff = io.loadmat(filevsh3, appendmat=False)
            #
            # This test is to fix a problem with 2 different
            # behavior of io.loadmat
            #
            if type(coeff['fmin']) == float:
                fmin = coeff['fmin']
                fmax = coeff['fmax']
            else:
                fmin = coeff['fmin'][0][0]
                fmax = coeff['fmax'][0][0]
            # .. Warning    
            # Warning modification take only one dimension for k 
            # if the .vsh3 format evolve it may not work anymore 
            #
            Br = SHCoeff('s3', fmin, fmax, coeff['Br.s3'],
                         coeff['Br.ind'], coeff['Br.k'][0])
            Bi = SHCoeff('s3', fmin, fmax, coeff['Bi.s3'],
                         coeff['Bi.ind'], coeff['Bi.k'][0])
            Cr = SHCoeff('s3', fmin, fmax, coeff['Cr.s3'],
                         coeff['Cr.ind'], coeff['Cr.k'][0])
            Ci = SHCoeff('s3', fmin, fmax, coeff['Ci.s3'],
                         coeff['Ci.ind'], coeff['Ci.k'][0])
            self.C = VSHCoeff(Br, Bi, Cr, Ci)
            Nf = np.shape(Br.s3)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filevsh3, ' does not exist'

    def savevsh2(self):
        """

        Create a .vsh2 antenna file

        """

        # create vsh2 file

        _filevsh2 = self._filename.replace('.trx', '.vsh2')
        filevsh2 = pyu.getlong(_filevsh2, pstruc['DIRANT'])

        if os.path.isfile(filevsh2):
            print filevsh2, ' already exist'
        else:
            print 'create ', filevsh2, ' file'

            coeff = {}
            coeff['fmin'] = self.fa[0]
            coeff['fmax'] = self.fa[-1]
            coeff['Br.ind'] = self.C.Br.ind2
            coeff['Bi.ind'] = self.C.Bi.ind2
            coeff['Cr.ind'] = self.C.Cr.ind2
            coeff['Ci.ind'] = self.C.Ci.ind2
            coeff['Br.s2'] = self.C.Br.s2
            coeff['Bi.s2'] = self.C.Bi.s2
            coeff['Cr.s2'] = self.C.Cr.s2
            coeff['Ci.s2'] = self.C.Ci.s2

            io.savemat(filevsh2, coeff, appendmat=False)

    def loadvsh2(self):
        """
            A.loadvsh2()

        Load antenna's vsh2 file which only contains
        the vsh coefficients in shape 2
        """

        _filevsh2 = self._filename
        filevsh2 = pyu.getlong(_filevsh2, pstruc['DIRANT'])

        if os.path.isfile(filevsh2):
            coeff = io.loadmat(filevsh2, appendmat=False)
            #
            # This test is to fix a problem with 2 different
            # behavior of io.loadmat
            #
            if type(coeff['fmin']) == float:
                fmin = coeff['fmin']
                fmax = coeff['fmax']
            else:
                fmin = coeff['fmin'][0][0]
                fmax = coeff['fmax'][0][0]
            Br = SHCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Br.s2'], ind=coeff['Br.ind'])
            Bi = SHCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Bi.s2'], ind=coeff['Bi.ind'])
            Cr = SHCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Cr.s2'], ind=coeff['Cr.ind'])
            Ci = SHCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Ci.s2'], ind=coeff['Ci.ind'])
            self.C = VSHCoeff(Br, Bi, Cr, Ci)
            Nf = np.shape(Br.s2)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filevsh2, ' does not exist'

    def loadvsh3_old(self):
        """ Load antenna vsh coefficients in shape 3

        """

        _filevsh3 = self._filename
        filevsh3 = getlong(_filevsh3, pstruc['DIRANT'])
        fmin = 2.
        fmax = 8.
        if os.path.isfile(filevsh3):
            coeff = io.loadmat(filevsh3, appendmat=False)
            Br = SHCoeff('s3', fmin, fmax, coeff['Br.s3'],
                         coeff['Br.ind'], coeff['Br.k'])
            Bi = SHCoeff('s3', fmin, fmax, coeff['Bi.s3'],
                         coeff['Bi.ind'], coeff['Bi.k'])
            Cr = SHCoeff('s3', fmin, fmax, coeff['Cr.s3'],
                         coeff['Cr.ind'], coeff['Cr.k'])
            Ci = SHCoeff('s3', fmin, fmax, coeff['Ci.s3'],
                         coeff['Ci.ind'], coeff['Ci.k'])
            self.C = VSHCoeff(Br, Bi, Cr, Ci)
            self.fa = np.linspace(fmin, fmax, 121)
        else:
            print _filevsh3, ' does not exist'

    def pol2cart(self, ith):
        """
        Conversion FTheta, FPhi to Fx,Fy,Fz for theta=ith

        Parameters
        ----------
        ith : theta index
        """
        Fth = self.Ftheta[:, ith, :]
        Fph = self.Fphi[:, ith, :]
        th = self.theta[ith]
        ph = self.phi

        Fx = Fth * np.cos(th) * np.cos(ph) - Fph * np.sin(ph)
        Fy = Fth * np.cos(th) * np.sin(ph) + Fph * np.cos(ph)
        Fz = (-1) * Fth * np.sin(th)

        return(Fx, Fy, Fz)

    def cart2pol(self, Fx, Fy, Fz, ith):
        """
        Conversion Fx,Fy,Fz vers Ftheta, Fphi for theta=ith
        Parameters
        ----------
        Fx
        Fy
        Fz
        ith : theta index

        """
        th = self.theta[ith]
        ph = self.phi

        Fth = Fx * np.cos(th) * np.cos(
            ph) + Fy * np.cos(th) * np.sin(ph) - Fz * np.sin(th)
        Fph = -Fx * np.sin(ph) + Fy * np.cos(th)

        SqG = np.sqrt(np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth)))
        self.SqG[:, ith, :] = SqG
        self.Ftheta[:, ith, :] = Fth
        self.Fphi[:, ith, :] = Fph


def forcesympol(A):
    """
        calculate A.Ftheta and A.Fphi in order to obtain A.SqG continuous at poles
        (theta=0,pi) for each frequency values

        Parameters
        ----------
        A : Antenna object

    """
    (Fx0, Fy0, Fz0) = A.pol2cart(0)
    (Fxp, Fyp, Fzp) = A.pol2cart(-1)

    for i in range(A.Nf):
        aux0 = mean(Fx0[i, :])
        auxp = mean(Fxp[i, :])
        Fx0[i, :] = aux0
        Fxp[i, :] = auxp

        aux0 = mean(Fy0[i, :])
        auxp = mean(Fyp[i, :])
        Fy0[i, :] = aux0
        Fyp[i, :] = auxp

        aux0 = mean(Fz0[i, :])
        auxp = mean(Fzp[i, :])
        Fz0[i, :] = aux0
        Fzp[i, :] = auxp

    A.cart2pol(Fx0, Fy0, Fz0, 0)
    A.cart2pol(Fxp, Fyp, Fzp, -1)

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

def VW3(l, m, theta ,phi ):
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

    """

    if type(l) == float:
        l = np.array([l])
    if type(m) == float:
        m = np.array([m])

    L = np.max(l)
    M = np.max(m)

    x = -np.cos(theta)

    # The - sign is necessary to get the good reconstruction
    #     deduced from observation
    #     May be it comes from a different definition of theta in SPHEREPACK

    Pmm1l, Pmp1l = AFLegendre(L, M, x)

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

    W[np.isinf(W) | np.isnan(W)] = 0
    V[np.isinf(V) | np.isnan(V)] = 0

    return V, W


def VW(n, m, x, phi, Pmm1n, Pmp1n):
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


    Example
    -------

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


def plotVW(n, m, theta, phi, sf=False):
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

        >>> from pylayers.antprop.antenna import *
        >>> import  matplotlib.pyplot as plt
        >>> import numpy as np
        >>> n=5
        >>> m=3
        >>> theta = np.linspace(0,np.pi,30)
        >>> phi   = np.linspace(0,2*np.pi,60)
        >>> plotVW(n,m,theta,phi)

    """
    # calculate v and w
    if m <= n:
        theta[np.where(theta == np.pi / 2)[0]] = np.pi / 2 + \
            1e-10  # .. todo :: not clean
        x = -np.cos(theta)
        Pmm1n, Pmp1n = AFLegendre(n, m, x)

        t1 = np.sqrt((n + m) * (n - m + 1))
        t2 = np.sqrt((n - m) * (n + m + 1))
        y1 = t1 * Pmm1n[:, m, n] - t2 * Pmp1n[:, m, n]
        y2 = t1 * Pmm1n[:, m, n] + t2 * Pmp1n[:, m, n]

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

        const = (-1.0) ** n / (2 * np.sqrt(n * (n + 1)))
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

        # plot V and W
        Ntheta = np.size(theta)
        vt = np.ones(Ntheta)
        Nphi = np.size(phi)
        vp = np.ones(Nphi)
        Phi = np.outer(vt, phi)
        Theta = np.outer(theta, vp)

        #figdirV='/home/rburghel/Bureau/bases_decomposition_VW/base_V_Vsin_Vcos/'
        figdirV = './'
        ext1 = '.pdf'
        ext2 = '.eps'
        ext3 = '.png'

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(V) * np.cos(Phi) * np.sin(Theta)
        Y = abs(V) * np.sin(Phi) * np.sin(Theta)
        Z = abs(V) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirV + 'V' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Vcos) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Vcos) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Vcos) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])

        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirV + 'Vcos' + str(n) + str(m) + '.jpg'
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Vsin) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Vsin) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Vsin) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirV + 'Vsin' + str(n) + str(m) + '.jpg'
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        #figdirW='/home/rburghel/Bureau/bases_decomposition_VW/base_W_Wsin_Wcos/'
        figdirW = './'

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(W) * np.cos(Phi) * np.sin(Theta)
        Y = abs(W) * np.sin(Phi) * np.sin(Theta)
        Z = abs(W) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirW + 'W' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Wcos) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Wcos) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Wcos) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirW + 'Wcos' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Wsin) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Wsin) * np.sin(Phi) * np.sin(Theta)

        fig = plt.figure()
        ax = axes3d.Axes3D(fig)
        X = abs(Wsin) * np.cos(Phi) * np.sin(Theta)
        Y = abs(Wsin) * np.sin(Phi) * np.sin(Theta)
        Z = abs(Wsin) * np.cos(Theta)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
        ax.set_xlim3d([-1, 1])
        ax.set_ylim3d([-1, 1])
        ax.set_zlim3d([-1, 1])
        if sf:
            sz = fig.get_size_inches()
            fig.set_size_inches(sz * 1.8)
            figname = figdirW + 'Wsin' + str(n) + str(m)
            fig.savefig(figname + ext1, orientation='portrait')
            fig.savefig(figname + ext2, orientation='portrait')
            fig.savefig(figname + ext3, orientation='portrait')

        plt.show()

    else:
        print "Error: m>n!!!"


def compdiag(k, A, th, ph, Fthr, Fphr, typ='modulus', lang='english', fontsize=18):
    """

    Comparison between Antenna diagram and reconstructed diagram

    Parameters
    ----------
    k : frequency index
    A : Antenna

    ph : phi base      (1 x Np)
    th : theta base    (1 x Nt)

    Fthr : Fth output of Fsynth Nf x (Ntheta*Tphi)
    Fphr : Fth output of Fsynth Nf x (Ntheta*Tphi)

    lang = 'french'
         = 'english'
    """

    Nf = np.shape(Fthr)[0]

    #Fthr = Fthr.reshape(Nf,len(th),len(ph))
    #Fphr = Fphr.reshape(Nf,len(th),len(ph))

    plt.figure()
    rc('text', usetex=True)
    Ftho = A.Ftheta
    Fpho = A.Fphi

    # limites module Fthr, Ftho, Fphr, Fpho
    maxTr = abs(Fthr[k, :, :]).max()
    maxTo = abs(Ftho[k, :, :]).max()
    MmT = max(maxTr, maxTo)

    minTr = abs(Fthr[k, :, :]).min()
    minTo = abs(Ftho[k, :, :]).min()
    mmT = min(minTr, minTo)

    maxPr = abs(Fphr[k, :, :]).max()
    maxPo = abs(Fpho[k, :, :]).max()
    MmP = max(maxPr, maxPo)

    minPr = abs(Fphr[k, :, :]).min()
    minPo = abs(Fpho[k, :, :]).min()
    mmP = min(minPr, minPo)

    # limites real Fthr, Ftho, Fphr, Fpho
    maxTrr = np.real(Fthr[k, :, :]).max()
    maxTor = np.real(Ftho[k, :, :]).max()
    MrT = max(maxTrr, maxTor)

    minTrr = np.real(Fthr[k, :, :]).min()
    minTor = np.real(Ftho[k, :, :]).min()
    mrT = min(minTrr, minTor)

    maxPrr = np.real(Fphr[k, :, :]).max()
    maxPor = np.real(Fpho[k, :, :]).max()
    MrP = max(maxPrr, maxPor)

    minPrr = np.real(Fphr[k, :, :]).min()
    minPor = np.real(Fpho[k, :, :]).min()
    mrP = min(minPrr, minPor)

    # limites real Fthr, Ftho, Fphr, Fpho
    maxTri = np.imag(Fthr[k, :, :]).max()
    maxToi = np.imag(Ftho[k, :, :]).max()
    MiT = max(maxTri, maxToi)

    minTri = np.imag(Fthr[k, :, :]).min()
    minToi = np.imag(Ftho[k, :, :]).min()
    miT = min(minTri, minToi)

    maxPri = np.imag(Fphr[k, :, :]).max()
    maxPoi = np.imag(Fpho[k, :, :]).max()
    MiP = max(maxPri, maxPoi)

    minPri = np.imag(Fphr[k, :, :]).min()
    minPoi = np.imag(Fpho[k, :, :]).min()
    miP = min(minPri, minPoi)

    # limithes arg Fth,Fph
    maxATr = np.angle(Fthr[k, :, :]).max()
    maxATo = np.angle(Ftho[k, :, :]).max()
    maT = max(maxATr, maxATo)
    minATr = np.angle(Fthr[k, :, :]).min()
    minATo = np.angle(Ftho[k, :, :]).min()
    maT0 = min(minATr, minATo)

    maxAPr = np.angle(Fphr[k, :, :]).max()
    maxAPo = np.angle(Fpho[k, :, :]).max()
    maP = max(maxAPr, maxAPo)
    minAPr = np.angle(Fphr[k, :, :]).min()
    minAPo = np.angle(Fpho[k, :, :]).min()
    maP0 = min(minAPr, minAPo)

    ax = plt.axes([0, 0, 360, 180])
    rtd = 180 / np.pi

    plt.subplot(221)
    if typ == 'modulus':
    #
    #cmap=cm.jet
        #pcolor(A.phi*rtd,A.theta*rtd,abs(Ftho[k,:,:]),vmin=0,vmax=mmT)
            #
    #cmap= gray
    #pcolor(A.phi*rtd,A.theta*rtd,abs(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
            #
    #cmap=cm.hot
        plt.pcolor(A.phi * rtd, A.theta * rtd, abs(Ftho[k, :, :]),
                   cmap=cm.hot_r, vmin=mmT, vmax=MmT)
        plt.title(r'$|F_{\theta}|$ original', fontsize=fontsize)

    if typ == 'real':
        #pcolor(A.phi*rtd,A.theta*rtd,real(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.real(Ftho[k, :, :]),
                   cmap=cm.hot_r, vmin=mrT, vmax=MrT)
        title(r'Re ($F_{\theta}$) original', fontsize=fontsize)
    if typ == 'imag':
        #pcolor(A.phi*rtd,A.theta*rtd,imag(Ftho[k,:,:]),cmap=cm.gray_r,vmin=0,vmax=mmT)
        pcolor(A.phi * rtd, A.theta * rtd, np.imag(Ftho[k, :, :]),
               cmap=cm.hot_r, vmin=miT, vmax=MiT)
        title(r'Im ($F_{\theta}$) original', fontsize=fontsize)
    if typ == 'phase':
        #pcolor(A.phi*rtd,A.theta*rtd,angle(Ftho[k,:,:]),cmap=cm.gray_r,vmin=maT0,vmax=maT)
        plt.pcolor(A.phi * rtd, A.theta * rtd, angle(Ftho[k, :, :]),
                   cmap=cm.hot_r, vmin=maT0, vmax=maT)
        if lang == 'french':
            plt.title(r'Arg ($F_{\theta}$) original', fontsize=fontsize)
        else:
            plt.title(r'Ang ($F_{\theta}$) original', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.ylabel(r'$\theta$ (deg)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)

    plt.subplot(222)
    if typ == 'modulus':
        plt.pcolor(A.phi * rtd, A.theta * rtd, abs(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=mmP, vmax=MmP)
        plt.title('$|F_{\phi}|$ original', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.real(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=mrP, vmax=MrP)
        plt.title('Re ($F_{\phi}$) original', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(A.phi * rtd, A.theta * rtd, np.imag(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=miP, vmax=MiP)
        plt.title('Im ($F_{\phi}$) original', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, angle(Fpho[k, :, :]),
                   cmap=cm.hot_r, vmin=maP0, vmax=maP)
        if lang == 'french':
            plt.title('Arg ($F_{\phi}$) original', fontsize=fontsize)
        else:
            plt.title('Ang ($F_{\phi}$) original', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)

    plt.subplot(223)
    if typ == 'modulus':
        plt.pcolor(ph * rtd, th * rtd, abs(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=mmT, vmax=MmT)
        if lang == 'french':
            plt.title(r'$|F_{\theta}|$ reconstruit', fontsize=fontsize)
        else:
            plt.title(r'$|F_{\theta}|$ reconstructed', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(ph * rtd, th * rtd, np.real(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=mrT, vmax=MrT)
        if lang == 'french':
            title(r'Re ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            title(r'Re ($F_{\theta}$) reconstructed', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(ph * rtd, th * rtd, np.imag(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=miT, vmax=MiT)
        if lang == 'french':
            plt.title(r'Im ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            plt.title(r'Im ($F_{\theta}$) reconstructed', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, angle(Fthr[k, :, :]),
                   cmap=cm.hot_r, vmin=maT0, vmax=maT)
        if lang == 'french':
            plt.title(r'Arg ($F_{\theta}$) reconstruit', fontsize=fontsize)
        else:
            plt.title(r'Ang ($F_{\theta}$) reconstructed', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.xlabel(r'$\phi$ (deg)', fontsize=fontsize)
    plt.ylabel(r'$\theta$ (deg)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)

    plt.subplot(224)
    if typ == 'modulus':
        plt.pcolor(ph * rtd, th * rtd, abs(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=mmP, vmax=MmP)
        if lang == 'french':
            plt.title('$|F_{\phi}|$ reconstruit', fontsize=fontsize)
        else:
            plt.title('$|F_{\phi}|$ reconstructed', fontsize=fontsize)
    if typ == 'real':
        plt.pcolor(ph * rtd, th * rtd, np.real(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=mrP, vmax=MrP)
        if lang == 'french':
            plt.title('Re ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Re ($F_{\phi}$) reconstructed', fontsize=fontsize)
    if typ == 'imag':
        plt.pcolor(ph * rtd, th * rtd, np.imag(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=miP, vmax=MiP)
        if lang == 'french':
            plt.title('Im ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Im ($F_{\phi}$) reconstructed', fontsize=fontsize)
    if typ == 'phase':
        plt.pcolor(A.phi * rtd, A.theta * rtd, angle(Fphr[k, :, :]),
                   cmap=cm.hot_r, vmin=maP0, vmax=maP)
        if lang == 'french':
            plt.title('Arg ($F_{\phi}$) reconstruit', fontsize=fontsize)
        else:
            plt.title('Ang ($F_{\phi}$) reconstructed', fontsize=fontsize)
    plt.axis([0, 360, 0, 180])
    plt.xlabel(r'$\phi$ (deg)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(fontsize)


def show3D(F, theta, phi, k, col=True):
    """ show 3D matplotlib diagram

    Parameters
    ----------

    F      : ndarray (Nf,Nt,Np)
    theta  : ndarray (1xNt)
        angle
    phi    : ndarray (1xNp)
        angle
    theta  : ndarray (Nt)
    k      : int
        frequency index
    col    : boolean
        if col  -> color coded plot3D
        if col == False -> simple plot3D

    Examples
    --------

    .. plot::
        :include-source:

        >>> import matplotlib.pyplot as plt
        >>> from pylayers.antprop.antenna import *
        >>> ifreq = 0
        >>> A     = Antenna('vsh3','defant.vsh3')
        >>> A.Nt  = 30
        >>> A.Np  = 60
        >>> A.Nf  = len(A.fa)
        >>> theta = np.linspace(0,np.pi,A.Nt)
        >>> phi   = np.linspace(0,2*np.pi,A.Np)
        >>> th    = np.kron(theta,np.ones(A.Np))
        >>> ph    = np.kron(np.ones(A.Nt),phi)
        >>> Fth3,Fph3 = A.Fsynth3(th,ph)
        >>> FTh3 = Fth3.reshape(A.Nf,A.Nt,A.Np)
        >>> FPh3 = Fph3.reshape(A.Nf,A.Nt,A.Np)
        >>> show3D(FTh3,theta,phi,ifreq)
        >>> txt = plt.title('show3D example')
        >>> plt.show()


    Warnings
    --------
            len(theta) must be equal with shape(F)[1]
            len(phi) must be equal with shape(F)[2]
    """
    nth = len(theta)
    nph = len(phi)

    if k >= np.shape(F)[0]:
        print 'Error: frequency index k not in F defined interval'
    if nth != np.shape(F)[1]:
        print 'Error: shape mistmatch between theta and F'

    if nph != np.shape(F)[2]:
        print 'Error: shape mistmatch between phi and F'

    fig = plt.figure()
    ax = axes3d.Axes3D(fig)

    V = F[k, :, :]

    vt = np.ones(nth)
    vp = np.ones(nph)
    Th = np.outer(theta, vp)
    Ph = np.outer(vt, phi)

    X = abs(V) * np.cos(Ph) * np.sin(Th)
    Y = abs(V) * np.sin(Ph) * np.sin(Th)
    Z = abs(V) * np.cos(Th)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    if (col):
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot_r)
    else:
        ax.plot3D(np.ravel(X), np.ravel(Y), np.ravel(Z))


if (__name__ == "__main__"):
    doctest.testmod()
