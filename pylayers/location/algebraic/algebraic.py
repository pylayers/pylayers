# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#Foobar is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Foobar is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Bernard UGUEN          : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################

import numpy as np
import scipy as sp
from scipy import optimize
import numpy.linalg as nplg
import cvxmod as cvxm
import cvxopt as cvxo
import matplotlib.pylab as plt
import string

class algloc(object):
    """
    This class regroups all the algebraic localization sceanrios and
    techniques
    Attributes
    ----------
        RN : An Array that defines the Radio nodes implied in
           localization (coordiantes in meters)
           shape(RN)= (2 or 3,RNnum)
        c : speed of light
          float
    """

    def __init__(self, RN):
        if len(np.shape(RN))==2:
            self.RN = RN
        else:
            raise ValueError("inputs must be of shape (n,p), p=2 or 3")
        self.c = 3.e08

    def info(self):
        """
        Dispaly scenario information
        """
        print "Reference Radio Nodes : ", self.RN
        print "Speed of light : ", self.c

    def show(self):
        """
        plots the refernces nodes (scenario)

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> import matplotlib.pyplot as plt
            >>> RN = np.array([[0., 0., 10., 10.],[0., 10., 10., 0.]])
            >>> S = algloc(RN)
            >>> S.show() 
        """
        plt.plot(self.RN[0,:], self.RN[1,:],'D')
        plt.xlablel('X')
        plt.ylabel('Y')
        plt.show()

    def dist(self,x,y,ax):
        """
        calculates distance between two arrays along a given axis
        Parameters
        ----------
            x : numpy.ndarray
            y : numpy.ndarray
            ax : integer (0,1)

        Returns
        -------
            d : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> x = np.array([[0., 0., 10., 10.],[0., 10., 10., 0.]])
            >>> y = np.array([[5.],[5.]])
            >>> S = algloc(x)
            >>> ax = 0
            >>> S.dist(x,y,ax) 
        """
        d = np.sqrt(np.sum((x-y)**2, axis=ax))
        return d

        
class Toa(algloc):
    """
    This class contains:
    1- a set of RadioNodes (RN),
    2- a set of TOAs measurements (TOA) with associated std (TOA_std)

    This class manages the TOA-based localization techniques.
    Attributes 
    -----------

        RN : An Array that defines the Radio nodes implied in
           localization (coordiantes in meters)
           shape(RN)= (2 or 3,RNnum)
        TOA : A measurement vector of TOA associated to RN
            (TOA values in seconds)
            shape(TOA)= (RNnum,1)
        TOA_std : Associated std of TOA (std in seconds)
                shape(TOA_std)= (RNnum,1)
        c : speed of light
          float 

    """
    def __init__(self, RN, TOA, TOA_std):
        algloc.__init__(self, RN)
        if len(np.shape(TOA))==2:
            if np.shape(TOA)[0]==1:
                self.TOA = TOA[0,:]
            elif np.shape(TOA)[1]==1:
                self.TOA = TOA[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TOA))==1:
            self.TOA = TOA
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(TOA_std))==2:
            if np.shape(TOA_std)[0]==1:
                self.TOA_std = TOA_std[0,:]
            elif np.shape(TOA_std)[1]==1:
                self.TOA_std = TOA_std[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TOA_std))==1:
            self.TOA_std = TOA_std
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")

    def info(self):
        """
        Dispaly scenario information
        """
        print "Reference Radio Nodes : ", self.RN
        print "Measured TOA : ", self.TOA
        print "std of Measured TOA : ", self.TOA_std
        print "Light speed : ", self.c

    def ls_locate(self):
        """
        This method applies least squares (LS) approximation on TOA to
        get position P.    
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.ls_locate()
        """
        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0))
        k1 = RN2[1:] - RN2[0:1]
        ROA = self.c * self.TOA
        ROA2 = ROA * ROA
        k2 = ROA2[0:1] - ROA2[1:]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:].T - self.RN[:, 0]
        # Apply LS operator
        P = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
        return P.reshape(np.shape(self.RN[:,0:1]))

    def tls_locate(self):
        """
        This method applies total least squares (TLS) approximation on
        TOA to get position P.  
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.tls_locate()
        """

        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0))
        k1 = RN2[1:] - RN2[0:1]
        ROA = self.c * self.TOA
        ROA2 = ROA * ROA
        k2 = ROA2[0:1] - ROA2[1:]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:].T - self.RN[:, 0]
        A2 = np.dot(np.transpose(A), A)
        [U, S, V] = nplg.svd(A2)
        J = 1/S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= self.c * max(self.TOA_std):
            f = f + 1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0
        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        P = 0.5 * np.dot(A2i, np.dot(A.T, K))
        return P.reshape(np.shape(self.RN[:,0:1]))

    def wls_locate(self):
        """
        This method applies weighted least squares (WLS) approximation
        on TOA to get position P. 
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.wls_locate()
        """
        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0))
        k1 = RN2[1:] - RN2[0:1]
        ROA = self.c * self.TOA
        ROA2 = ROA * ROA
        k2 = ROA2[0:1] - ROA2[1:]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:].T - self.RN[:, 0]
        # Construct the Covariance Matrix
        ROA_std = self.c * self.TOA_std
        C = np.diag((ROA_std[1:]) ** 2)
        # Apply LS operator
        P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),A))), \
            np.dot(np.dot(A.T, nplg.inv(C)), K))
        return P.reshape(np.shape(self.RN[:,0:1]))

    def twls_locate(self):
        """
        This method applies total weighted least squares (TWLS)
        approximation on TOA to get position P.  
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.twls_locate()
        """

        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0))
        k1 = RN2[1:] - RN2[0:1]
        ROA = self.c * self.TOA
        ROA2 = ROA * ROA
        k2 = ROA2[0:1] - ROA2[1:]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:].T - self.RN[:, 0]
        # Construct the Covariance Matrix
        ROA_std = self.c * self.TOA_std
        C = np.diag((ROA_std[1:]) ** 2)
        A2 = np.dot(A.T, np.dot(nplg.inv(C), A))
        [U, S, V] = nplg.svd(A2)
        J = 1/S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= self.c * max(self.TOA_std):
            f = f + 1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0
        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        P = 0.5 * np.dot(A2i, np.dot(np.dot(A.T, nplg.inv(C)), K))
        return P.reshape(np.shape(self.RN[:,0:1]))

    def ts_locate(self, P0, Niter):
        """
        This method applies taylor series (TS) method on TOA
        to get position P.
        Parameters
        ----------
            P0 : numpy.ndarray
            Niter : integer 
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> Niter = 10
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.ts_locate(BN0, Niter)
        """

        P = P0
        ROA = self.c * self.TOA
        ROA_std = self.c * self.TOA_std
        for i in np.arange(Niter):
            # Construct the matrix A (see theory)
            A = ((P - self.RN) / np.sqrt( \
                np.sum((P - self.RN) ** 2, axis=0))).T
            # Construct the Covariance Matrix
            C = np.diag((ROA_std[:]) ** 2)
            # Construct the vector D (see theory)
            D = ROA - np.sqrt((np.sum((P - self.RN) ** 2, axis=0)))
            # construct the vector Delta (see theory)
            Delta = np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv( \
                C), A))), np.dot(np.dot(A.T, nplg.inv(C)), D))
            # update P
            P = P + Delta.reshape(np.shape(P))
        return P

    def sdp_locate(self):
        """
        This method applies semidefinite programming (SDP) on TOA
        to get position P.  
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.sdp_locate()
        """
        RN = cvxm.matrix(self.RN)
        TOA = cvxm.matrix(self.TOA)
        ROA = self.c * TOA
        ROA_std = self.c * self.TOA_std
        ROA_std = cvxm.matrix(ROA_std)
        mtoa, ntoa = cvxm.size(RN)
        Im = cvxm.eye(mtoa)
        Y = cvxm.optvar('Y', mtoa + 1, mtoa + 1)
        t = cvxm.optvar('t', ntoa, 1)
        prob = cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
        prob.constr.append(Y >= 0)
        prob.constr.append(Y[mtoa, mtoa] == 1)
        for i in range(ntoa):
            X0 = cvxm.matrix([[Im, -cvxm.transpose(RN[:, i])], \
                 [-RN[:, i], cvxm.transpose(RN[:, i]) * RN[:, i]]])
            prob.constr.append(-t[i] < (cvxm.trace(
                X0 * Y) - ROA[i] ** 2) * (1 / ROA_std[i]))
            prob.constr.append(t[i] > (cvxm.trace(
                X0 * Y) - ROA[i] ** 2) * (1 / ROA_std[i]))
        prob.solve(quiet=True)
        Pval = Y.value
        P = Pval[:mtoa, -1]
        return P

    def ml_function(self, P):
        """
        This defines the ML function to be minimized
        if ML estimator is used
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            ML : numpy.ndarray

        Notes
        -----
        Attention ! This function does not return a position

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.ml_function(BN)
        """
        RNnum = np.shape(self.RN)[1]
        ROA = self.c * self.TOA
        ROA_std = self.c * self.TOA_std
        # construct the ML function to be minimized
        RNmP = self.RN - np.outer(P, np.ones(RNnum))
        mRNmP = np.sqrt(np.diag(np.dot(RNmP.T,RNmP)))
        tk = (ROA-mRNmP) ** 2
        uk = tk/(2*ROA_std**2)+np.log(np.sqrt(2*np.pi)*ROA_std)
        ML = uk.sum(axis=0)
        return(ML)

    def ml_locate(self, P0):
        """
        This method applies maximum likelihood (ML) method on TOA
        to get position P.
        Parameters
        ----------
            P0 : numpy.ndarray  
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.ml_locate(BN0)
        """
        P = optimize.fmin(self.ml_function, P0, args=(),disp=0)
        return P.reshape(np.shape(P0))

    def fim(self, P):
        """
        Compute the fisher information matrix in P for the given
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            FIM : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.fim(BN)
        """
        ROA_std = self.c * self.TOA_std
        FIM = np.zeros((np.shape(self.RN)[0],np.shape(self.RN)[0]))
        for i in range(np.shape(self.RN)[1]):
            FIM += np.dot((P-self.RN[:,i:i+1]),(P-self.RN[:,i:i+1]).T)/\
            ((ROA_std[i]**2)*np.dot((P-self.RN[:,i:i+1]).T,\
            (P-self.RN[:,i:i+1])))
        return FIM

    def crb(self, P):
        """
        This method compute the cramer rao bound (CRB) at position P.
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            CRB : float

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> S = algloc(RN)
            >>> d = S.dist(RN,BN,0) # actual distances
            >>> TOF = d/S.c # actual TOAs
            >>> err_std = 1.0
            >>> TOA_std = err_std/S.c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
            >>> T = TOA(RN, TOA, TOA_std)
            >>> T.crb(BN)
        """
        FIM = self.fim(P)
        CRB = np.sqrt(np.trace(nplg.inv(FIM)))
        return CRB

class Rss(algloc):
    """
    This class contains:
    1- a set of RadioNodes (RN),
    2- a set of RSSs measurements (RSS) with associated std of shadowing
    (RSS_std) and associated propagation constant (RSS_np)

    This class manages the RSS-based localization techniques. The path
    loss model used here is log normal shadowing:
    RSS = PL0 + 10*np*log10(d/d0)+X(0,RSS_std)

    Attributes
    ----------

            RN : An Array that defines the Radio nodes implied in
               localization (coordiantes in meters)
               shape(RN)= (2 or 3,RNnum)
            RSS : A measurement vector of RSS associated to RN
                (RSS values in dB)
                shape(RSS)= (RNnum,1)
            d0 : the reference distance (usually equal to 1 meter)
               float
            RSS_std : std of shadowing (std in dB)
                    shape(RSS_std)= (RNnum,1)
            RSS_np : propagation constant
                   shape(RSS_np)= (RNnum,1)
            PL0 : RSS at d0
            Rest : range estimator ('mode', 'median', 'mean')
    """
    def __init__(self, RN, RSS, d0, RSS_std, RSS_np, PL0, Rest):
        algloc.__init__(self, RN)
        if len(np.shape(RSS))==2:
            if np.shape(RSS)[0]==1:
                self.RSS = RSS[0,:]
            elif np.shape(RSS)[1]==1:
                self.RSS = RSS[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(RSS))==1:
            self.RSS = RSS
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(RSS_std))==2:
            if np.shape(RSS_std)[0]==1:
                self.RSS_std = RSS_std[0,:]
            elif np.shape(RSS_std)[1]==1:
                self.RSS_std = RSS_std[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(RSS_std))==1:
            self.RSS_std = RSS_std
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(RSS_np))==2:
            if np.shape(RSS_np)[0]==1:
                self.RSS_np = RSS_np[0,:]
            elif np.shape(RSS_np)[1]==1:
                self.RSS_np = RSS_np[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(RSS_np))==1:
            self.RSS_np = RSS_np
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")

        if len(np.shape(PL0))==2:
            if np.shape(PL0)[0]==1:
                self.PL0 = PL0[0,:]
            elif np.shape(PL0)[1]==1:
                self.PL0 = PL0[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(PL0))==1:
            self.PL0 = PL0
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")

        self.d0 = d0
        self.Rest = Rest

    def info(self):
        """
        Dispaly scenario information
        """
        print "Reference Radio Nodes:\n", self.RN
        print "Measured RSS:\n", self.RSS
        print "References distances:\n", self.d0
        print "Propagation constants:\n", self.RSS_np
        print "std of Measured RSS shadowing:\n", self.RSS_std
        print "RSS at d0:\n", self.RSS_std
        print "Range estimator:\n", self.Rest

    def get_pl0(self, lamda):
        """
        Compute the thoretic path loss PL0 at distance d0
        Parameters
        ----------
            lamda : float
        Returns
        -------
            PL0 : float

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> R.get_plo(lamda)
        """
        PL0 = 20 * np.log10(4 * np.pi * self.d0 / lamda)
        return PL0

    def get_plmean(self, P):
        """
        Compute the averaged path loss at position P
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            PLmean : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> R.get_plmean(BN)
        """
        RNmP = self.dist(self.RN,P,0)
        PLmean = self.PL0 - 10 * self.RSS_np * np.log10(RNmP / self.d0)
        return PLmean

    def get_pl(self, P):
        """
        Compute the path loss PL at position P
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            PL : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> R.get_pl(BN)
        """

        X = self.RSS_std * np.random.randn(np.shape(self.PL0)[0])
        PL = self.get_plmean(P) + X
        return PL

    def get_range(self):
        """
        Compute the range using the "Rest" estimator
        Parameters
        ----------
            Rest : string
        Returns
        -------
            Range : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> R.get_range()
        """
        S = (np.log(10)/10)*self.RSS_std/self.RSS_np
        M = (np.log(10)/10)*(self.PL0-self.RSS)/self.RSS_np+\
            np.log(self.d0)

        
        if string.lower(self.Rest) == 'mode':
            Range = np.exp(M-S**2)
        elif string.lower(self.Rest) == 'median':
            Range = np.exp(M)
        elif string.lower(self.Rest) == 'mean':
            Range = np.exp(M+0.5*S**2)
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return Range

    def get_range_std(self):
        """
        Compute the range standard deviation using the "Rest" estimator
        Parameters
        ----------
            Rest : string
        Returns
        -------
            Range_std : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> R.get_range_std()
        """
        S = (np.log(10)/10)*self.RSS_std/self.RSS_np
        M = (np.log(10)/10)*(self.PL0-self.RSS)/self.RSS_np+\
             np.log(self.d0)
        if string.lower(self.Rest) == 'mode':
            Range_std = np.sqrt((np.exp(2*M-2*S**2))*(1-np.exp(-S**2)))
        elif string.lower(self.Rest) == 'median':
            Range_std = np.sqrt((np.exp(2*M+S**2))*(np.exp(S**2)-1))
        elif string.lower(self.Rest) == 'mean':
            Range_std = np.sqrt((np.exp(2*M+3*S**2))*(np.exp(S**2)-1))
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return Range_std

    def ls_locate(self):
        """
        This applies LS approximation on RSS based ranges to get
        position P.
        Parameters
        ----------
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.ls_locate()
        """
        shRN = np.shape(self.RN)
        RNnum = shRN[1]
        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0)).reshape(RNnum, 1)
        k1 = RN2[1:RNnum, :] - RN2[0, 0]
        ROA = self.get_range()
        ROA2 = (ROA * ROA).reshape(RNnum, 1)
        k2 = ROA2[0, 0] - ROA2[1:RNnum, :]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:RNnum].T - self.RN[:, 0].reshape(1, shRN[0])
        # Apply LS operator
        P = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
        return P

    def tls_locate(self):
        """
        This applies TLS approximation on RSS based ranges to get
        position P.
        Parameters
        ----------
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.tls_locate()
        """
        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0))
        k1 = RN2[1:] - RN2[0:1]
        ROA = self.get_range()
        ROA2 = ROA * ROA
        k2 = ROA2[0:1] - ROA2[1:]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:].T - self.RN[:, 0]
        A2 = np.dot(np.transpose(A), A)
        [U, S, V] = nplg.svd(A2)
        J = 1 / S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= max(self.get_range_std()):
            f = f + 1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        P = 0.5 * np.dot(A2i, np.dot(A.T, K))
        # Return the estimated position
        return P.reshape(np.shape(self.RN[:,0:1]))

    def wls_locate(self):
        """
        This applies WLS approximation on RSS based ranges to get
        position P.
        Parameters
        ----------
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.wls_locate()
        """
        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0))
        k1 = RN2[1:] - RN2[0:1]
        ROA = self.get_range()
        ROA_std = self.get_range_std()
        ROA2 = ROA * ROA
        k2 = ROA2[0:1] - ROA2[1:]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:].T - self.RN[:, 0]
        # Construct the Covariance Matrix
        C = np.diag((ROA_std[1:]) ** 2)
        # Apply LS operator
        P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))),\
            np.dot(np.dot(A.T, nplg.inv(C)), K))
        return P.reshape(np.shape(self.RN[:,0:1]))

    def twls_locate(self):
        """
        This applies TWLS approximation on RSS based ranges to get
        position P.
        Parameters
        ----------
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.twls_locate()
        """
        # Construct the vector K (see theory)
        RN2 = (np.sum(self.RN * self.RN, axis=0))
        k1 = RN2[1:] - RN2[0:1]
        ROA = self.get_range()
        ROA_std = self.get_range_std()
        ROA2 = ROA * ROA
        k2 = ROA2[0:1] - ROA2[1:]
        K = k1 + k2
        # Construct the matrix A (see theory)
        A = self.RN[:, 1:].T - self.RN[:, 0]
        # Construct the Covariance Matrix
        C = np.diag((ROA_std[1:]) ** 2)
        A2 = np.dot(A.T, np.dot(nplg.inv(C), A))
        [U, S, V] = nplg.svd(A2)
        J = 1 / S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= max(ROA_std):
            f = f + 1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        P = 0.5 * np.dot(A2i, np.dot(np.dot(A.T, nplg.inv(C)), K))
        return P.reshape(np.shape(self.RN[:,0:1]))

    def ts_locate(self, P0, Niter):
        """
        This method applies taylor series (TS) method on RSS
        to get position P.
        Parameters
        ----------
            Rest : string
            P0 : numpy.ndarray
            Niter : integer 
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> Niter = 10
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.ts_locate(BN0, Niter)
        """

        P = P0
        ROA = self.get_range()
        ROA_std = self.get_range_std()
        for i in np.arange(Niter):
            # Construct the matrix A (see theory)
            A = ((P - self.RN) / np.sqrt( \
                np.sum((P - self.RN) ** 2, axis=0))).T
            # Construct the Covariance Matrix
            C = np.diag((ROA_std[:]) ** 2)
            # Construct the vector D (see theory)
            D = ROA - np.sqrt((np.sum((P - self.RN) ** 2, axis=0)))
            # construct the vector Delta (see theory)
            Delta = np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv( \
                C), A))), np.dot(np.dot(A.T, nplg.inv(C)), D))
            # update P
            P = P + Delta.reshape(np.shape(P))
        return P

    def sdp_locate(self):
        """
        This method applies semidefinite programming (SDP) on RSS to get
         position P.
        Parameters
        ----------
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.sdp_locate()
        """
        ROA = self.get_range()
        RN = cvxm.matrix(self.RN)
        RSS = cvxm.matrix(self.RSS)
        RSS_np = cvxm.matrix(self.RSS_np)
        RSS_std = cvxm.matrix(self.RSS_std)
        PL0 = cvxm.matrix(self.PL0)
        ROA = cvxm.matrix(ROA)
        mrss, nrss = cvxm.size(RN)
        Si = np.array([(1/self.d0**2)*10**((RSS[0,0]-PL0[0,0])/(5.0*\
            RSS_np[0,0])),(1/self.d0**2)*10**((RSS[1,0]-PL0[1,0])/\
            (5.0*RSS_np[1,0])),(1/self.d0**2)*10**((RSS[2,0]-PL0[2,0])\
            /(5.0*RSS_np[2,0])),(1/self.d0**2)*10**((RSS[3,0]-PL0[3,0])\
            /(5.0*RSS_np[3,0]))])
        Im = cvxm.eye(mrss)
        Y = cvxm.optvar('Y', mrss + 1, mrss + 1)
        t = cvxm.optvar('t', nrss, 1)
        prob = cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
        prob.constr.append(Y >= 0)
        prob.constr.append(Y[mrss, mrss] == 1)
        for i in range(nrss):
            X0 = cvxm.matrix([[Im,-cvxm.transpose(RN[:,i])],\
                 [-RN[:,i],cvxm.transpose(RN[:,i])*RN[:, i]]])
            prob.constr.append(-RSS_std[i,0]*t[i]<Si[i]*\
                               cvxm.trace(X0*Y)-1)
            prob.constr.append(RSS_std[i,0]*t[i]>Si[i]*\
                               cvxm.trace(X0*Y)-1)
        prob.solve(quiet=True)
        Pval = Y.value
        P = Pval[:mrss, -1]
        return P

    def iml_function(self, P):
        """
        This defines the ML function (computed using RSS based ranges
        obtained using 'Rest' estimator) to be minimized if ML estimator
        (iml_locate) is used
        Parameters
        ----------
            P : numpy.ndarray
            Rest : string
        Returns
        -------
            ML : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.iml_function(BN)
        """
        RNnum = np.shape(self.RN)[1]
        ROA = self.get_range()
        ROA_std = self.get_range_std()
        # construct the ML function to be minimized
        RNmP = self.RN - np.outer(P, np.ones(RNnum))
        mRNmP = np.sqrt(np.diag(np.dot(RNmP.T,RNmP)))
        tk = (ROA-mRNmP) ** 2
        uk = tk/(2*ROA_std**2)+np.log(np.sqrt(2*np.pi)*ROA_std)
        ML = uk.sum(axis=0)
        return(ML)

    def iml_locate(self, P0):
        """
        This method applies maximum likelihood (ML) method on RSS based
        ranges to get position P.
        Parameters
        ----------
            P0 : numpy.ndarray
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.iml_locate(BN0)
        """

        P = optimize.fmin(self.iml_function, P0, args=(),disp=0)
        return P.reshape(np.shape(P0))

    def dml_function(self, P):
        """
        This defines the ML function (computed using RSS) to be
        minimized if ML estimator (dml_locate) is used
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            ML : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.dml_function(BN)
        """
        RNnum = np.shape(self.RN)[1]
        S = (np.log(10)/10)*self.RSS_std/self.RSS_np
        M = (np.log(10)/10)*(self.PL0-self.RSS)/self.RSS_np+\
            np.log(self.d0)
        # construct the ML function to be minimized
        RNmP = self.RN-np.outer(P, np.ones(RNnum))
        mRNmP = np.sqrt(np.diag(np.dot(RNmP.T,RNmP)))
        tk = (M-S**2-np.log(mRNmP))**2
        uk = tk/(2*S**2)
        ML = uk.sum(axis=0)
        return(ML)

    def dml_locate(self, P0):
        """
        This method applies direct maximum likelihood (ML) method on RSS
        to get position P.
        Parameters
        ----------
            P0 : numpy.ndarray
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.iml_locate(BN0)
        """

        P = optimize.fmin(self.dml_function, P0, args=(), disp=0)
        return P.reshape(np.shape(P0))

    def fim(self, P):
        """
        Compute the fisher information matrix in P for the given
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            FIM : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.fim(BN)

        """
        S = (np.log(10)/10)*self.RSS_std/self.RSS_np
        FIM = np.zeros((np.shape(self.RN)[0],np.shape(self.RN)[0]))
        for i in range(np.shape(self.RN)[1]):
            FIM += np.dot((P-self.RN[:,i:i+1]),(P-self.RN[:,i:i+1]).T)/\
                 ((S[0]**2)*(np.dot((P-self.RN[:,i:i+1]).T,\
                 (P-self.RN[:,i:i+1])))**2)
        return FIM

    def crb(self, P):
        """
        This method compute the cramer rao bound (CRB) at position P.
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            CRB : float

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> nRN = 4
            >>> dim = 2 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,nRN)
            >>> RSS_std = 4.34 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> RSS = 110*sp.rand(nRN)
            >>> PL0 = -34.7*np.ones(nRN)
            >>> lamda = 5e9
            >>> d0 = 1.
            >>> Rest = 'mode' # RSS based ranging estimator
            >>> R = Rss(RN, RSS, d0, RSS_std, RSS_np, PL0, Rest)
            >>> PL = R.get_pl(BN)
            >>> R.RSS = PL
            >>> R.crb(BN)
        """

        FIM = self.fim(P)
        CRB = np.sqrt(np.trace(nplg.inv(FIM)))
        return CRB

class Tdoa(algloc):
    """
    A TDoALocation contains:
    1- 2 sets of RadioNodes (RN and RNr)
    2- a set of TDoAs measurements (TDOA) with associated std (TDOA_std)

    This class manages the TDOA-based localization techniques.
    MEMBERS:

            RN : An Array that defines the set of radio nodes implied in
                 localization (coordiantes in meters)
                 shape(RN1)= (2 or 3,RNnum)
            RNr : An Array that defines the set of reference radio nodes
                  with whom which TDOA are computed
                  shape(RNr)= (2 or 3,RNnum)
            TDOA : A measurement vector of TDoA associated to RN
                   (TDoA values in seconds)
                   shape(TDoA)= (RNnum,1)
            TDOA_std : Associated std of TDoA (std in seconds)
                       shape(TDoAstd)= (RNnum,1)
    """

    def __init__(self, RN, RNr, TDOA, TDOA_std):
        algloc.__init__(self, RN)
        if len(np.shape(RNr))==2:
            self.RNr = RNr
        else:
            raise ValueError("inputs must be of shape (n,p), p=2 or 3")
        if len(np.shape(TDOA))==2:
            if np.shape(TDOA)[0]==1:
                self.TDOA = TDOA[0,:]
            elif np.shape(TDOA)[1]==1:
                self.TDOA = TDOA[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TDOA))==1:
            self.TDOA = TDOA
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(TDOA_std))==2:
            if np.shape(TDOA_std)[0]==1:
                self.TDOA_std = TDOA_std[0,:]
            elif np.shape(TDOA_std)[1]==1:
                self.TDOA_std = TDOA_std[:,1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TDOA_std))==1:
            self.TDOA_std = TDOA_std
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")


    def ls_locate(self):
        """
        This applies LS approximation on TDOA to get position P.
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.ls_locate()
        """
        shRN = np.shape(self.RN)
        # Construct the vector K (see theory)
        k1 = (np.sum((self.RN-self.RNr) * (self.RN-self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA2 = RDOA*RDOA
        k2 = RDOA2
        K = k1-k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T-self.RNr.T, \
        RDOA.reshape(np.shape(self.TDOA)[0],1)))
        # Apply LS operator
        Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
        P = Pr[:shRN[0]].reshape(shRN[0],1)
        return P

    def tls_locate(self):
        """
        This applies TLS approximation on TDOA to get position P.
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.tls_locate()
        """
        shRN = np.shape(self.RN)
        # Construct the vector K (see theory)
        k1 = (np.sum((self.RN-self.RNr) * (self.RN-self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA2 = RDOA*RDOA
        k2 = RDOA2
        K = k1-k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T-self.RNr.T, \
        RDOA.reshape(np.shape(self.TDOA)[0],1)))
        A2 = np.dot(A.T, A)
        [U, S, V] = nplg.svd(A2)
        J = 1/S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= self.c * max(self.TDOA_std):
            f = f+1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        # Apply LS operator
        Pr = 0.5 * np.dot(A2i, np.dot(A.T, K))
        P = Pr[:shRN[0]].reshape(shRN[0],1)
        return P

    def wls_locate(self):
        """
        This applies WLS approximation on TDOA to get position P.
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.wls_locate()
        """
        shRN = np.shape(self.RN)
        # Construct the vector K (see theory)
        k1 = (np.sum((self.RN-self.RNr) * (self.RN-self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA_std = self.c * self.TDOA_std
        RDOA2 = RDOA*RDOA
        k2 = RDOA2
        K = k1-k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T-self.RNr.T, \
        RDOA.reshape(np.shape(self.TDOA)[0],1)))
        # Construct the Covariance Matrix
        C = np.diag(RDOA_std[:] ** 2)
        # Apply LS operator
        Pr = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),A))),\
             np.dot(np.dot(A.T, nplg.inv(C)), K))
        P = Pr[:shRN[0]].reshape(shRN[0],1)
        return P

    def twls_locate(self):
        """
        This applies TWLS approximation on TDOA to get position P.
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.twls_locate()
        """
        shRN = np.shape(self.RN)
        # Construct the vector K (see theory)
        k1 = (np.sum((self.RN-self.RNr) * (self.RN-self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA_std = self.c * self.TDOA_std
        RDOA2 = RDOA*RDOA
        k2 = RDOA2
        K = k1-k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T-self.RNr.T, \
        RDOA.reshape(np.shape(self.TDOA)[0],1)))
        # Construct the Covariance Matrix
        C = np.diag(RDOA_std[:] ** 2)
        # Apply LS operator
        A2 = np.dot(A.T,np.dot(nplg.inv(C), A))
        [U, S, V] = nplg.svd(A2)
        J = 1/S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= self.c * max(self.TDOA_std):
            f = f+1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        Pr = 0.5*np.dot(A2i, np.dot(np.dot(A.T, nplg.inv(C)), K))
        P = Pr[:shRN[0]].reshape(shRN[0],1)
        return P

    def ts_locate(self, P0, Niter):
        """
        This method applies taylor series (TS) method on TDOA
        to get position P.
        Parameters
        ----------
            P0 : numpy.ndarray
            Niter : integer 
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> Niter = 10
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.ts_locate(BN0, Niter)
        """

        P = P0
        RDOA = self.c * self.TDOA
        RDOA_std = self.c * self.TDOA_std
        for i in np.arange(Niter):
            # Construct the matrix A (see theory)
            A = ((P-self.RN)/np.sqrt(np.sum((P-self.RN)**2,axis=0))).T-\
                ((P-self.RNr)/np.sqrt(np.sum((P-self.RNr)**2,axis=0))).T
            # Construct the Covariance Matrix
            C = np.diag((RDOA_std[:]) ** 2)
            # Construct the vector D (see theory)
            D = RDOA - (np.sqrt((np.sum((P-self.RN)**2, axis=0)))-\
                       np.sqrt((np.sum((P-self.RNr)**2, axis=0))))
            # construct the vector Delta (see theory)
            Delta = np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv( \
                C), A))), np.dot(np.dot(A.T, nplg.inv(C)), D))
            # update P
            P = P + Delta.reshape(np.shape(P))
        return P

    def ml_function(self, P):
        """
        This defines the ML function to be minimized
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            ML : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.ml_function(BN)
        """

        shRN = np.shape(self.RN)
        RNnum = shRN[1]
        RDOA = self.c*self.TDOA
        RDOA_std = self.c*self.TDOA_std

        # construct the ML function to be minimized
        RN1mP = self.RN - np.outer(P, np.ones(RNnum))
        mRN1mP = np.sqrt(np.diag(np.dot(RN1mP.T, RN1mP)))
        RN2mP = self.RNr - np.outer(P, np.ones(RNnum))
        mRN2mP = np.sqrt(np.diag(np.dot(RN2mP.T, RN2mP)))
        rRDOA = mRN1mP - mRN2mP

        tk = (RDOA - rRDOA) ** 2 / (2 * RDOA_std ** 2)
        ML = tk.sum(axis=0)

        return(ML)

    def ml_locate(self, P0):
        """
        This method applies maximum likelihood (ML) method on TDOA
        to get position P.
        Parameters
        ----------
            P0 : numpy.ndarray  
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.ml_locate(BN0)
        """
        P = optimize.fmin(self.ml_function, P0, args=(),disp=0)
        return P.reshape(np.shape(P0))

    def fim(self, P):
        """
        Compute the fisher information matrix in P for the given
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            FIM : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.fim(BN)
        """
        RDOA_std = self.c * self.TDOA_std
        FIM = np.zeros((np.shape(self.RN)[0],np.shape(self.RN)[0]))
        for i in range(np.shape(self.RN)[1]):
            PmRN = np.sqrt(np.dot((P-self.RN[:,i:i+1]).T,\
                   P-self.RN[:,i:i+1]))
            PmRNr = np.sqrt(np.dot((P-self.RNr[:,i:i+1]).T,\
                    P-self.RNr[:,i:i+1]))
            FIM += (1/RDOA_std[i]**2)*np.dot((P-self.RN[:,i:i+1])/PmRN-\
                (P-self.RNr[:,i:i+1])/PmRNr,((P-self.RN[:,i:i+1])/PmRN-\
                (P-self.RNr[:,i:i+1])/PmRNr).T)
        return FIM

    def crb(self, P):
        """
        This method compute the cramer rao bound (CRB) at position P.
        Parameters
        ----------
            P : numpy.ndarray
        Returns
        -------
            CRB : float

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> BN = L*sp.rand(dim,1)
            >>> RN = L*sp.rand(dim,4)
            >>> RNr = np.zeros((dim,nRN))
            >>> S = algloc(RN)
            >>> dr = S.dist(RNr,BN,0)
            >>> TDOF = (d-dr)/S.c # actual TDOA
            >>> TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.crb(BN)
        """
        FIM = self.fim(P)
        CRB = np.sqrt(np.trace(nplg.inv(FIM)))
        return CRB

    
