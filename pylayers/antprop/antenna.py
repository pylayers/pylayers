# -*- coding:Utf-8 -*-
"""
This module handles antennas in pylayers

To instantiate an antenna object :

.. python::

    A = Antenna(_filename,directory,nf,ntheta,nphi)

typ indicates the antenna file format to read

Examples
--------

    >>> from pylayers.antprop.antenna import *
    >>> A = Antenna('S1R1.mat','ant/UWBAN/Matfile')

The antenna can be represented in various formats

.vsh2
.vsh3



"""
import doctest
import os
import re
import sys
import pdb
import numpy as np
import scipy as sp
import scipy.special as special
from scipy import io
import matplotlib.pylab as plt
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
from pylayers.util.project import *
from pylayers.antprop.spharm import *
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import axes3d
#from scipy import sparse
from matplotlib import rc
from matplotlib import cm # colormaps
from pylayers.antprop.antssh import *

class Antenna(object):

    def Fsynth3(self, theta = [], phi=[], pattern=True):
        """ synthesis of a complex antenna pattern from VSH coefficients (shape 3)

        Ndir is the number of directions

        Parameters
        ----------

        theta : ndarray (1xNdir if not pattern)  (1xNtheta if pattern)
        phi   : ndarray (1xNdir if not pattter)  (1xNphi if pattern)
        
        pattern : boolean
            if True theta and phi are reorganized for building the pattern

        Returns
        -------

        if pattern:
            Fth   : ndarray (Ntheta x Nphi)
            Fph   : ndarray (Ntheta x Nphi)
        else:
            Fth   : ndarray (1 x Ndir)
            Fph   : ndarray (1 x Ndir)

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.antprop.antenna import *
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> A = Antenna('defant.vsh3')
            >>> theta = np.linspace(0,np.pi,70)
            >>> phi = np.linspace(0,2*np.pi,180)
            >>> F = A.Fsynth3(theta,phi,pattern=True)

        All Br,Cr,Bi,Ci have the same (l,m) index in order to evaluate only
        once the V,W function

        """
        typ = self._filename.split('.')[1]
        Nf = len(self.fa)
        if theta==[]:
            theta=np.linspace(0,np.pi,47)
        if phi == []:
            phi= np.linspace(0,2*np.pi,47)

        Nt = len(theta)
        Np = len(phi)

        if pattern:
                self.theta = theta[:,np.newaxis]
                self.phi = phi[np.newaxis,:] 
                theta = np.kron(theta, np.ones(Np))
                phi = np.kron(np.ones(Nt),phi)
                         
        
        if typ =='vsh3':        
            

            nray = len(theta)

            Br  = self.C.Br.s3
            lBr = self.C.Br.ind3[:, 0]
            mBr = self.C.Br.ind3[:, 1]

            Bi  = self.C.Bi.s3

            Cr  = self.C.Cr.s3

            Ci  = self.C.Ci.s3

            L = lBr.max()
            M = mBr.max()

            # vector spherical harmonics basis functions

            V, W = VW(lBr, mBr, theta, phi)
            Fth = np.dot(Br, np.real(V.T)) - \
                np.dot(Bi, np.imag(V.T)) + \
                np.dot(Ci, np.real(W.T)) + \
                np.dot(Cr, np.imag(W.T))

            Fph = -np.dot(Cr, np.real(V.T)) + \
                np.dot(Ci, np.imag(V.T)) + \
                np.dot(Bi, np.real(W.T)) + \
                np.dot(Br, np.imag(W.T))
            
            if pattern:
                
                Fth = Fth.reshape(Nf, Nt, Np)
                Fph = Fph.reshape(Nf, Nt, Np)
                self.Fphi = Fph
                self.Ftheta = Fth

                G = np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth))
                self.SqG = np.sqrt(G)

            
        if typ == 'sh3':
            cx = self.S.Cx.s3
            cy = self.S.Cy.s3
            cz = self.S.Cz.s3
            lmax = self.S.Cx.lmax
            Y ,indx = SSHFunc2(lmax, theta,phi)
            k = self.S.Cx.k2[:,0]
            if pattern :
                    
                Ex = np.dot(cx,Y[k])
                Ey = np.dot(cy,Y[k])
                Ez = np.dot(cz,Y[k])
                Fth,Fph = CartToSphere (theta, phi, Ex, Ey,Ez, bfreq = True, pattern = True ) 
                Fth = Fth.reshape(Nf,Nt,Np)
                Fph = Fph.reshape(Nf,Nt,Np)
            else:
                     
                Ex = np.dot(cx,Y[k])
                Ey = np.dot(cy,Y[k])
                Ez = np.dot(cz,Y[k])
                Fth,Fph = CartToSphere (theta, phi, Ex, Ey,Ez, bfreq = True, pattern = False)       
            
            self.Fphi = Fph
            self.Ftheta = Fth
            G = np.real(Fph * np.conj(Fph) + Fth * np.conj(Fth))
            self.SqG = np.sqrt(G)
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

        #th = np.kron(self.theta, np.ones(self.Np))
        #ph = np.kron(np.ones(self.Nt), self.phi)

        Fth3, Fph3 = self.Fsynth3(self.theta, self.phi,pattern=True)
        Err = self.mse(Fth3, Fph3, 0)

        Enc = self.C.ens3()
        n = len(Enc)
        pos = 0

        while (pos < n) & (Err[0] < emax):

            Emin = Enc[pos]
            d = self.C.drag3(Emin)
            Fth3, Fph3 = self.Fsynth3(self.theta, self.phi,pattern=True)
            Err = self.mse(Fth3, Fph3, 0)

            if Err[0] >= emax:
                i = d[0][0]
                i3 = d[1][0]
                self.C.put3(i, i3)
                Fth3, Fph3 = self.Fsynth3(self.theta,self.phi,pattern=True)
                Err = self.mse(Fth3, Fph3, 0)

            pos = pos + 1

    def savevsh3(self):
        """ save antenna in vsh3 format

        Create a .vsh3 antenna file


        """

        # create vsh3 file

        _filevsh3 = os.path.splitext(self._filename)[0]+'.vsh3'
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

    def savesh2(self):
        """
        Create a .sh2 antenna file

        """

        # create sh2 file
        typ = self._filename.split('.')[1]
        self.typ = typ
    
        _filesh2 = self._filename.replace('.'+ typ, '.sh2')
        filesh2 = pyu.getlong(_filesh2, pstruc['DIRANT'])
        if os.path.isfile(filesh2):
            print filesh2, ' already exist'
        else:
            print 'create ', filesh2, ' file'
            coeff = {}
            coeff['fmin'] = self.fa[0]
            coeff['fmax'] = self.fa[-1]

            
            coeff['Cx.ind'] = self.S.Cx.ind2
            coeff['Cy.ind'] = self.S.Cy.ind2
            coeff['Cz.ind'] = self.S.Cz.ind2
            coeff['Cx.lmax']= self.S.Cx.lmax           
            coeff['Cy.lmax']= self.S.Cy.lmax           
            coeff['Cz.lmax']= self.S.Cz.lmax           

            coeff['Cx.s2'] = self.S.Cx.s2
            coeff['Cy.s2'] = self.S.Cy.s2
            coeff['Cz.s2'] = self.S.Cz.s2
            io.savemat(filesh2, coeff, appendmat=False)
            
    def savesh3(self):
        
        """ 
        save antenna in sh3 format

        Create a .sh3 antenna file


        """
        
        
        # create sh3 file
        typ = self._filename.split('.')[1]
        self.typ = typ    
        _filesh3 = self._filename.replace('.'+ typ, '.sh3')
        filesh3 = pyu.getlong(_filesh3, pstruc['DIRANT'])
        if os.path.isfile(filesh3):
            print filesh3, ' already exist'

       
        else:
            print 'create ', filesh3, ' file'

            coeff = {}
            coeff['fmin'] = self.fa[0]
            coeff['fmax'] = self.fa[-1]
            coeff['Cx.ind'] = self.S.Cx.ind3
            coeff['Cy.ind'] = self.S.Cy.ind3
            coeff['Cz.ind'] = self.S.Cz.ind3
            
            coeff['Cx.k'] = self.S.Cx.k2
            coeff['Cy.k'] = self.S.Cy.k2
            coeff['Cz.k'] = self.S.Cz.k2
            

            coeff['Cx.lmax']= self.S.Cx.lmax           
            coeff['Cy.lmax']= self.S.Cy.lmax           
            coeff['Cz.lmax']= self.S.Cz.lmax 

            coeff['Cx.s3'] = self.S.Cx.s3
            coeff['Cy.s3'] = self.S.Cy.s3
            coeff['Cz.s3'] = self.S.Cz.s3
         

            io.savemat(filesh3, coeff, appendmat=False)

    def loadvsh3(self):
        """ Load antenna's vsh3 file

            vsh3 file contains a thesholded version of vsh coefficients in shape 3

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
            Br = VCoeff('s3', fmin, fmax, coeff['Br.s3'],
                         coeff['Br.ind'], coeff['Br.k'][0])
            Bi = VCoeff('s3', fmin, fmax, coeff['Bi.s3'],
                         coeff['Bi.ind'], coeff['Bi.k'][0])
            Cr = VCoeff('s3', fmin, fmax, coeff['Cr.s3'],
                         coeff['Cr.ind'], coeff['Cr.k'][0])
            Ci = VCoeff('s3', fmin, fmax, coeff['Ci.s3'],
                         coeff['Ci.ind'], coeff['Ci.k'][0])
            self.C = VSHCoeff(Br, Bi, Cr, Ci)
            Nf = np.shape(Br.s3)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filevsh3, ' does not exist'

    def loadsh3(self):
        """ Load antenna's sh3 file

            sh3 file contains a thesholded version of ssh coefficients in shape 3

        """
        _filesh3 = self._filename.split('.')[0]+'.sh3'
        filesh3 = pyu.getlong(_filesh3, pstruc['DIRANT'])

        if os.path.isfile(filesh3):
            coeff = io.loadmat(filesh3, appendmat=False)

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
            # Warning modification takes only one dimension for k 
            # if the .sh3 format evolve it may not work anymore 
            #

                      
            if type(coeff['Cx.lmax']) == float:
                lmax = coeff['Cx.lmax']
            else:
                lmax = coeff['Cx.lmax'][0][0]
            Cx = SCoeff(typ = 's3',
                        fmin = fmin ,
                        fmax = fmax , 
                        lmax = lmax,
                        data = coeff['Cx.s3'],
                        ind =  coeff['Cx.ind'],
                        k =  coeff['Cx.k'])

                        
            Cy = SCoeff(typ= 's3', 
                        fmin = fmin ,
                        fmax = fmax , 
                        lmax = lmax,
                        data = coeff['Cy.s3'],
                        ind =  coeff['Cy.ind'],
                        k =  coeff['Cy.k'])

                        
                         
            Cz = SCoeff(typ = 's3', 
                        fmin = fmin ,
                        fmax = fmax , 
                        data = coeff['Cz.s3'],

                        lmax = lmax,
                        ind =  coeff['Cz.ind'],
                        k =  coeff['Cz.k'])

            
            if not 'S' in self.__dict__.keys():
                self.S = SSHCoeff(Cx, Cy,Cz)
            else:
                self.S.sets3(Cx,Cy,Cz)
                    
            Nf = np.shape(Cx.s3)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filesh3, ' does not exist'
            
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
	

            
    def loadsh2(self):
        
        _filesh2 = self._filename.split('.')[0]+'.sh2'
        filesh2 = pyu.getlong(_filesh2, pstruc['DIRANT'])
        
        if os.path.isfile(filesh2):
            coeff = io.loadmat(filesh2, appendmat=False)
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

                
            if type(coeff['Cx.lmax']) == float:
                lmax = coeff['Cx.lmax']
            else:
                lmax = coeff['Cx.lmax'][0][0]
                
                
            Cx = SCoeff(typ='s2', 
                        fmin=fmin, 
                        fmax=fmax,
                        lmax = lmax,
                        data=coeff['Cx.s2'], 
                        ind=coeff['Cx.ind'])
                        
            Cy = SCoeff(typ='s2', 
                        fmin=fmin, 
                        fmax=fmax,
                        lmax = lmax,
                        data=coeff['Cy.s2'], 
                        ind=coeff['Cy.ind'])
            Cz = SCoeff(typ='s2', 
                        fmin=fmin, 
                        fmax=fmax,
                        lmax = lmax,
                        data=coeff['Cz.s2'], 
                        ind=coeff['Cz.ind'])
                         

            self.S = SSHCoeff(Cx, Cy,Cz)
            Nf = np.shape(Cx.s2)[0]
            self.fa = np.linspace(fmin, fmax, Nf)
        else:
            print _filesh2, ' does not exist'
        
            

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

            Br = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Br.s2'], ind=coeff['Br.ind'])
            Bi = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Bi.s2'], ind=coeff['Bi.ind'])
            Cr = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
                         data=coeff['Cr.s2'], ind=coeff['Cr.ind'])
            Ci = VCoeff(typ='s2', fmin=fmin, fmax=fmax,
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
            Br = VCoeff('s3', fmin, fmax, coeff['Br.s3'],
                         coeff['Br.ind'], coeff['Br.k'])
            Bi = VCoeff('s3', fmin, fmax, coeff['Bi.s3'],
                         coeff['Bi.ind'], coeff['Bi.k'])
            Cr = VCoeff('s3', fmin, fmax, coeff['Cr.s3'],
                         coeff['Cr.ind'], coeff['Cr.k'])
            Ci = VCoeff('s3', fmin, fmax, coeff['Ci.s3'],
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
        >>> A     = Antenna('defant.vsh3')
        >>> A.Nt  = 30
        >>> A.Np  = 60
        >>> A.Nf  = len(A.fa)
        >>> theta = np.linspace(0,np.pi,A.Nt)
        >>> phi   = np.linspace(0,2*np.pi,A.Np)
        >>> Fth3,Fph3 = A.Fsynth3(theta,phi)
        >>> FTh3 = Fth3.reshape(A.Nf,A.Nt,A.Np)
        >>> FPh3 = Fph3.reshape(A.Nf,A.Nt,A.Np)
        >>> show3D(FTh3,theta,phi,ifreq)
        >>> txt = plt.title('show3D example')


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
    plt.ion()
    doctest.testmod()
