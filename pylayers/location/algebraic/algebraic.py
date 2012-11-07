# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#PYLAYERS is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#PYLAYERS is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with PYLAYERS.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Bernard UGUEN          : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import numpy as np
import scipy as sp
from scipy import optimize
import numpy.linalg as nplg
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import string
import doctest

class algloc(object):
    """ regroups algebraic localization scenarios and techniques

    Attributes
    ----------
        RN : An Array that defines the Radio nodes implied in
           localization (coordinates in meters)
           shape(RN)= (2 or 3,RNnum)
        c : speed of light
          float
    """

    def __init__(self, RN):
        if len(np.shape(RN)) == 2:
            self.RN = RN
        else:
            raise ValueError("inputs must be of shape (n,p), p=2 or 3")
        self.c = 3.e08

    def info(self):
        """ display scenario information

        """
        print "Reference Radio Nodes : ", self.RN
        print "Speed of light : ", self.c

    def show(self):
        """ plots the reference nodes (scenario)

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
        plt.plot(self.RN[0, :], self.RN[1, :], 'D')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.axis('scaled')
        plt.show()

    def dist(self, x, y, ax):
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
        d = np.sqrt(np.sum((x - y) ** 2, axis=ax))
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
        if len(np.shape(TOA)) == 2:
            if np.shape(TOA)[0] == 1:
                self.TOA = TOA[0, :]
            elif np.shape(TOA)[1] == 1:
                self.TOA = TOA[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TOA)) == 1:
            self.TOA = TOA
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(TOA_std)) == 2:
            if np.shape(TOA_std)[0] == 1:
                self.TOA_std = TOA_std[0, :]
            elif np.shape(TOA_std)[1] == 1:
                self.TOA_std = TOA_std[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TOA_std)) == 1:
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
        return P.reshape(np.shape(self.RN[:, 0:1]))

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
        J = 1 / S
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
        return P.reshape(np.shape(self.RN[:, 0:1]))

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
        P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))),
                         np.dot(np.dot(A.T, nplg.inv(C)), K))
        return P.reshape(np.shape(self.RN[:, 0:1]))

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
        J = 1 / S
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
        return P.reshape(np.shape(self.RN[:, 0:1]))

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
            A = ((P - self.RN) / np.sqrt(
                np.sum((P - self.RN) ** 2, axis=0))).T
            # Construct the Covariance Matrix
            C = np.diag((ROA_std[:]) ** 2)
            # Construct the vector D (see theory)
            D = ROA - np.sqrt((np.sum((P - self.RN) ** 2, axis=0)))
            # construct the vector Delta (see theory)
            Delta = np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(
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
            X0 = cvxm.matrix([[Im, -cvxm.transpose(RN[:, i])],
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
        mRNmP = np.sqrt(np.diag(np.dot(RNmP.T, RNmP)))
        tk = (ROA - mRNmP) ** 2
        uk = tk / (2 * ROA_std ** 2) + np.log(np.sqrt(2 * np.pi) * ROA_std)
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
        P = optimize.fmin(self.ml_function, P0, args=(), disp=0)
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
        FIM = np.zeros((np.shape(self.RN)[0], np.shape(self.RN)[0]))
        for i in range(np.shape(self.RN)[1]):
            FIM += np.dot((P - self.RN[:, i:i + 1]), (P - self.RN[:, i:i + 1]).T) /\
                ((ROA_std[i] ** 2) * np.dot((P - self.RN[:, i:i + 1]).T,
                                            (P - self.RN[:, i:i + 1])))
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
        if len(np.shape(RSS)) == 2:
            if np.shape(RSS)[0] == 1:
                self.RSS = RSS[0, :]
            elif np.shape(RSS)[1] == 1:
                self.RSS = RSS[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(RSS)) == 1:
            self.RSS = RSS
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(RSS_std)) == 2:
            if np.shape(RSS_std)[0] == 1:
                self.RSS_std = RSS_std[0, :]
            elif np.shape(RSS_std)[1] == 1:
                self.RSS_std = RSS_std[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(RSS_std)) == 1:
            self.RSS_std = RSS_std
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(RSS_np)) == 2:
            if np.shape(RSS_np)[0] == 1:
                self.RSS_np = RSS_np[0, :]
            elif np.shape(RSS_np)[1] == 1:
                self.RSS_np = RSS_np[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(RSS_np)) == 1:
            self.RSS_np = RSS_np
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")

        if len(np.shape(PL0)) == 2:
            if np.shape(PL0)[0] == 1:
                self.PL0 = PL0[0, :]
            elif np.shape(PL0)[1] == 1:
                self.PL0 = PL0[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(PL0)) == 1:
            self.PL0 = PL0
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")

        self.d0 = d0
        self.Rest = Rest

    def info(self):
        """
        Display scenario information
        """
        print "Reference Radio Nodes:\n", self.RN
        print "Measured RSS:\n", self.RSS
        print "References distances:\n", self.d0
        print "Propagation constants:\n", self.RSS_np
        print "std of Measured RSS shadowing:\n", self.RSS_std
        print "RSS at d0:\n", self.RSS_std
        print "Range estimator:\n", self.Rest

    def get_pl0(self, lamda):
        """ computes the theoretic path loss PL0 at distance d0

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
        RNmP = self.dist(self.RN, P, 0)
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
        S = (np.log(10) / 10) * self.RSS_std / self.RSS_np
        M = (np.log(10) / 10) * (self.PL0 - self.RSS) / self.RSS_np +\
            np.log(self.d0)

        if string.lower(self.Rest) == 'mode':
            Range = np.exp(M - S ** 2)
        elif string.lower(self.Rest) == 'median':
            Range = np.exp(M)
        elif string.lower(self.Rest) == 'mean':
            Range = np.exp(M + 0.5 * S ** 2)
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
        S = (np.log(10) / 10) * self.RSS_std / self.RSS_np
        M = (np.log(10) / 10) * (self.PL0 - self.RSS) / self.RSS_np +\
            np.log(self.d0)
        if string.lower(self.Rest) == 'mode':
            Range_std = np.sqrt((np.exp(
                2 * M - 2 * S ** 2)) * (1 - np.exp(-S ** 2)))
        elif string.lower(self.Rest) == 'median':
            Range_std = np.sqrt(
                (np.exp(2 * M + S ** 2)) * (np.exp(S ** 2) - 1))
        elif string.lower(self.Rest) == 'mean':
            Range_std = np.sqrt((np.exp(
                2 * M + 3 * S ** 2)) * (np.exp(S ** 2) - 1))
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
        return P.reshape(np.shape(self.RN[:, 0:1]))

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
        P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))),
                         np.dot(np.dot(A.T, nplg.inv(C)), K))
        return P.reshape(np.shape(self.RN[:, 0:1]))

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
        return P.reshape(np.shape(self.RN[:, 0:1]))

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
            A = ((P - self.RN) / np.sqrt(
                np.sum((P - self.RN) ** 2, axis=0))).T
            # Construct the Covariance Matrix
            C = np.diag((ROA_std[:]) ** 2)
            # Construct the vector D (see theory)
            D = ROA - np.sqrt((np.sum((P - self.RN) ** 2, axis=0)))
            # construct the vector Delta (see theory)
            Delta = np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(
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
        Si = np.array([(1 / self.d0 ** 2) * 10 ** ((RSS[0, 0] - PL0[0, 0]) / (5.0 * RSS_np[0, 0])), 
                       (1 / self.d0 ** 2) * 10 ** ((RSS[1, 0] - PL0[1, 0]) / (5.0 * RSS_np[1, 0])),
                       (1 / self.d0 ** 2) * 10 ** ((RSS[2, 0] - PL0[2, 0]) / (5.0 * RSS_np[2, 0])), 
                       (1 / self.d0 ** 2) * 10 ** ((RSS[3, 0] - PL0[3, 0]) / (5.0 * RSS_np[3, 0]))])

        Im = cvxm.eye(mrss)
        Y = cvxm.optvar('Y', mrss + 1, mrss + 1)
        t = cvxm.optvar('t', nrss, 1)
        prob = cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
        prob.constr.append(Y >= 0)
        prob.constr.append(Y[mrss, mrss] == 1)
        for i in range(nrss):
            X0 = cvxm.matrix([[Im, -cvxm.transpose(RN[:, i])],
                              [-RN[:, i], cvxm.transpose(RN[:, i]) * RN[:, i]]])
            prob.constr.append(-RSS_std[i, 0] * t[i] < Si[i] *
                               cvxm.trace(X0 * Y) - 1)
            prob.constr.append(RSS_std[i, 0] * t[i] > Si[i] *
                               cvxm.trace(X0 * Y) - 1)
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
        mRNmP = np.sqrt(np.diag(np.dot(RNmP.T, RNmP)))
        tk = (ROA - mRNmP) ** 2
        uk = tk / (2 * ROA_std ** 2) + np.log(np.sqrt(2 * np.pi) * ROA_std)
        ML = uk.sum(axis=0)
        return(ML)
=======
        nodes : dictionnary
        ldp : dictionnary
        c : speed of light
          float
        
    """

    def __init__(self, nodes, ldp):
        self.nodes = nodes
        self.ldp = ldp
        self.c = 0.3
>>>>>>> mlaaraie-master

    def info(self):
        """
        Display scenario information
        """
        print "Speed of light : ", self.c
        print "Nodes : ", self.nodes
        print "Location dependant parameters : ", self.ldp

<<<<<<< HEAD
        P = optimize.fmin(self.iml_function, P0, args=(), disp=0)
        return P.reshape(np.shape(P0))

    def dml_function(self, P):
=======
    def plot(self, rss, toa, tdoa):
>>>>>>> mlaaraie-master
        """
        Plot sceanrio
        Parameters
        ----------
<<<<<<< HEAD
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
            >>> ML = R.dml_function(BN)
        """
        RNnum = np.shape(self.RN)[1]
        S = (np.log(10) / 10) * self.RSS_std / self.RSS_np
        M = (np.log(10) / 10) * (self.PL0 - self.RSS) / self.RSS_np +\
            np.log(self.d0)
        # construct the ML function to be minimized
        RNmP = self.RN - np.outer(P, np.ones(RNnum))
        mRNmP = np.sqrt(np.diag(np.dot(RNmP.T, RNmP)))
        tk = (M - S ** 2 - np.log(mRNmP)) ** 2
        uk = tk / (2 * S ** 2)
        ML = uk.sum(axis=0)
        return(ML)
=======
            rss : boolean
            toa : boolean
            tdoa : boolean
        """
>>>>>>> mlaaraie-master

        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            fig = plt.figure()
            ax = Axes3D(fig)
            BN = self.nodes['BN']
            try:
                ax.plot(BN[0,:],BN[1,:],BN[2,:],'r*',zdir='z',\
                        label='blind node')
            except:
                plt.plot(BN[0,:],BN[1,:],'r*', label='blind node')
            if rss!=0:
                RN_RSS = self.nodes['RN_RSS']
                try:
                    ax.plot(RN_RSS[0,:],RN_RSS[1,:],RN_RSS[2,:],'ro',\
                            zdir='z', label='RSS node')
                except:
                    plt.plot(RN_RSS[0,:],RN_RSS[1,:],'ro',\
                            label='RSS node')
            if toa!=0:
                RN_TOA = self.nodes['RN_TOA']
                try:
                    ax.plot(RN_TOA[0,:],RN_TOA[1,:],RN_TOA[2,:],'gs',\
                            zdir='z', label='TOA node')
                except:
                    plt.plot(RN_TOA[0,:],RN_TOA[1,:], 'gs', \
                            label='TOA node')

            if toa!=0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                try:
                    ax.plot(RN_TDOA[0,:],RN_TDOA[1,:],RN_TDOA[2,:],'bD',\
                            zdir='z', label='TDOA node')
                    ax.plot(RNr_TDOA[0,:],RNr_TDOA[1,:],RNr_TDOA[2,:],\
                            'kD',zdir='z', label='Ref TDOA node')
                except:
                    plt.plot(RN_TDOA[0,:],RN_TDOA[1,:], 'bD', \
                            label='TDOA node')
                    plt.plot(RNr_TDOA[0,:],RNr_TDOA[1,:], 'kD', \
                            label='Ref TDOA node')

    def show(self, rss, toa, tdoa):
        """
        Plot sceanrio
        Parameters
        ----------
            rss : boolean
            toa : boolean
            tdoa : boolean

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> import numpy as np
            >>> from pylayers.location.algebraic.algebraic import *
            >>> from pylayers.util.geomutil import dist
            >>> nRN = 4
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> c = 0.3
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN_TOA = L*sp.rand(dim,nRN)
            >>> RN_RSS = L*sp.rand(dim,nRN)
            >>> RN_TDOA = L*sp.rand(dim,nRN)
            
            >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
            >>> TOF = d_TOA/c # actual TOA
            >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std

            >>> RSS_std = 0.001 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> PL0 = 34.7*np.ones(nRN)
            >>> d0 = 1.
            >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
            >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
            >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
            >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
            >>> d = dist(RN_TDOA,BN,0)
            >>> dr = dist(RNr_TDOA,BN,0)
            >>> TDOF = (d-dr)/c # actual TDOA
            >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std

            >>> nodes={}
            >>> nodes['BN']= BN
            >>> nodes['RN_RSS']= RN_RSS
            >>> nodes['RN_TOA']= RN_TOA
            >>> nodes['RN_TDOA']= RN_TDOA
            >>> nodes['RNr_TDOA']= RNr_TDOA

            >>> ldp={}
            >>> ldp['RSS'] = RSS
            >>> ldp['RSS_std'] = RSS_std
            >>> ldp['RSS_np'] = RSS_np
            >>> ldp['d0'] = d0
            >>> ldp['PL0'] = PL0
            >>> ldp['TOA'] = TOA
            >>> ldp['TOA_std'] = TOA_std
            >>> ldp['TDOA'] = TDOA
            >>> ldp['TDOA_std'] = TDOA_std
            
            >>> S = algloc(nodes, ldp)
            >>> S.show(1,1,1)
        """
        self.plot(rss, toa, tdoa)
        plt.legend(numpoints=1)
        plt.show()
        

    def get_range(self, Rest='mode'):
        """
        Compute the range from RSS using the "Rest" estimator
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
            >>> from pylayers.location.algebraic.algebraic import *
            >>> from pylayers.util.geomutil import dist
            >>> nRN = 4
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> c = 0.3
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN_TOA = L*sp.rand(dim,nRN)
            >>> RN_RSS = L*sp.rand(dim,nRN)
            >>> RN_TDOA = L*sp.rand(dim,nRN)
            
            >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
            >>> TOF = d_TOA/c # actual TOA
            >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std

            >>> RSS_std = 0.001 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> PL0 = 34.7*np.ones(nRN)
            >>> d0 = 1.
            >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
            >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
            >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
            >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
            >>> d = dist(RN_TDOA,BN,0)
            >>> dr = dist(RNr_TDOA,BN,0)
            >>> TDOF = (d-dr)/c # actual TDOA
            >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std

<<<<<<< HEAD
        """
        S = (np.log(10) / 10) * self.RSS_std / self.RSS_np
        FIM = np.zeros((np.shape(self.RN)[0], np.shape(self.RN)[0]))
        for i in range(np.shape(self.RN)[1]):
            FIM += np.dot((P - self.RN[:, i:i + 1]), (P - self.RN[:, i:i + 1]).T) /\
                ((S[0] ** 2) * (np.dot((P - self.RN[:, i:i + 1]).T,
                 (P - self.RN[:, i:i + 1]))) ** 2)
        return FIM
=======
            >>> nodes={}
            >>> nodes['BN']= BN
            >>> nodes['RN_RSS']= RN_RSS
            >>> nodes['RN_TOA']= RN_TOA
            >>> nodes['RN_TDOA']= RN_TDOA
            >>> nodes['RNr_TDOA']= RNr_TDOA

            >>> ldp={}
            >>> ldp['RSS'] = RSS
            >>> ldp['RSS_std'] = RSS_std
            >>> ldp['RSS_np'] = RSS_np
            >>> ldp['d0'] = d0
            >>> ldp['PL0'] = PL0
            >>> ldp['TOA'] = TOA
            >>> ldp['TOA_std'] = TOA_std
            >>> ldp['TDOA'] = TDOA
            >>> ldp['TDOA_std'] = TDOA_std
            
            >>> S = algloc(nodes, ldp)
            >>> r_mode = S.get_range('mode')
            >>> r_median = S.get_range('median')
            >>> r_mean = S.get_range('mean')
        """
        RSS = self.ldp['RSS']
        RSS_std = self.ldp['RSS_std']
        RSS_np = self.ldp['RSS_np']
        d0 = self.ldp['d0']
        PL0 = self.ldp['PL0']
        S = (np.log(10)/10)*RSS_std/RSS_np
        M = (np.log(10)/10)*(PL0-RSS)/RSS_np+np.log(d0)
        if string.lower(Rest) == 'mode':
            Range = np.exp(M-S**2)
        elif string.lower(Rest) == 'median':
            Range = np.exp(M)
        elif string.lower(Rest) == 'mean':
            Range = np.exp(M+0.5*S**2)
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return Range
>>>>>>> mlaaraie-master

    def get_range_std(self, Rest='mode'):
        """
<<<<<<< HEAD
        This method compute the Cramer Rao bound (CRB) at position P.

=======
        Compute the RSS range standard deviation using the "Rest" \
        estimator
>>>>>>> mlaaraie-master
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
            >>> from pylayers.location.algebraic.algebraic import *
            >>> from pylayers.util.geomutil import dist
            >>> nRN = 4
<<<<<<< HEAD
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
            >>> CRB = R.crb(BN)
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

    Attributes 
    ----------

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
        if len(np.shape(RNr)) == 2:
            self.RNr = RNr
        else:
            raise ValueError("inputs must be of shape (n,p), p=2 or 3")
        if len(np.shape(TDOA)) == 2:
            if np.shape(TDOA)[0] == 1:
                self.TDOA = TDOA[0, :]
            elif np.shape(TDOA)[1] == 1:
                self.TDOA = TDOA[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TDOA)) == 1:
            self.TDOA = TDOA
        else:
            raise ValueError("inputs shape must be (1,n) or (n,)")
        if len(np.shape(TDOA_std)) == 2:
            if np.shape(TDOA_std)[0] == 1:
                self.TDOA_std = TDOA_std[0, :]
            elif np.shape(TDOA_std)[1] == 1:
                self.TDOA_std = TDOA_std[:, 1]
            else:
                raise ValueError("inputs shape must be (1,n) or (n,)")
        elif len(np.shape(TDOA_std)) == 1:
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
=======
>>>>>>> mlaaraie-master
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> c = 0.3
            >>> BN = L*sp.rand(dim,1)
<<<<<<< HEAD
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
        k1 = (np.sum((self.RN - self.RNr) * (self.RN - self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA2 = RDOA * RDOA
        k2 = RDOA2
        K = k1 - k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T - self.RNr.T,
                       RDOA.reshape(np.shape(self.TDOA)[0], 1)))
        # Apply LS operator
        Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
        P = Pr[:shRN[0]].reshape(shRN[0], 1)
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
=======
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN_TOA = L*sp.rand(dim,nRN)
            >>> RN_RSS = L*sp.rand(dim,nRN)
            >>> RN_TDOA = L*sp.rand(dim,nRN)
            
            >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
            >>> TOF = d_TOA/c # actual TOA
            >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
>>>>>>> mlaaraie-master

            >>> RSS_std = 0.001 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> PL0 = 34.7*np.ones(nRN)
            >>> d0 = 1.
            >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
            >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
            >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
            >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
            >>> d = dist(RN_TDOA,BN,0)
            >>> dr = dist(RNr_TDOA,BN,0)
            >>> TDOF = (d-dr)/c # actual TDOA
            >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
<<<<<<< HEAD
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.tls_locate()
        """
        shRN = np.shape(self.RN)
        # Construct the vector K (see theory)
        k1 = (np.sum((self.RN - self.RNr) * (self.RN - self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA2 = RDOA * RDOA
        k2 = RDOA2
        K = k1 - k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T - self.RNr.T,
                       RDOA.reshape(np.shape(self.TDOA)[0], 1)))
        A2 = np.dot(A.T, A)
        [U, S, V] = nplg.svd(A2)
        J = 1 / S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= self.c * max(self.TDOA_std):
            f = f + 1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        # Apply LS operator
        Pr = 0.5 * np.dot(A2i, np.dot(A.T, K))
        P = Pr[:shRN[0]].reshape(shRN[0], 1)
        return P
=======
>>>>>>> mlaaraie-master

            >>> nodes={}
            >>> nodes['BN']= BN
            >>> nodes['RN_RSS']= RN_RSS
            >>> nodes['RN_TOA']= RN_TOA
            >>> nodes['RN_TDOA']= RN_TDOA
            >>> nodes['RNr_TDOA']= RNr_TDOA

            >>> ldp={}
            >>> ldp['RSS'] = RSS
            >>> ldp['RSS_std'] = RSS_std
            >>> ldp['RSS_np'] = RSS_np
            >>> ldp['d0'] = d0
            >>> ldp['PL0'] = PL0
            >>> ldp['TOA'] = TOA
            >>> ldp['TOA_std'] = TOA_std
            >>> ldp['TDOA'] = TDOA
            >>> ldp['TDOA_std'] = TDOA_std
            
            >>> S = algloc(nodes, ldp)
            >>> rs_mode = S.get_range_std('mode')
            >>> rs_median = S.get_range_std('median')
            >>> rs_mean = S.get_range_std('mean')
        """
        RSS = self.ldp['RSS']
        RSS_std = self.ldp['RSS_std']
        RSS_np = self.ldp['RSS_np']
        d0 = self.ldp['d0']
        PL0 = self.ldp['PL0']
        S = (np.log(10)/10)*RSS_std/RSS_np
        M = (np.log(10)/10)*(PL0-RSS)/RSS_np+np.log(d0)
        if string.lower(Rest) == 'mode':
            Range_std = np.sqrt((np.exp(2*M-2*S**2))*(1-np.exp(-S**2)))
        elif string.lower(Rest) == 'median':
            Range_std = np.sqrt((np.exp(2*M+S**2))*(np.exp(S**2)-1))
        elif string.lower(Rest) == 'mean':
            Range_std = np.sqrt((np.exp(2*M+3*S**2))*(np.exp(S**2)-1))
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return Range_std
        
    def ls_locate(self, rss, toa, tdoa, Rest):
        """
        This method applies least squares (LS) approximation to get
        position P.
        Parameters
        ----------
            rss : boolean
            toa : boolean
            tdoa : boolean
            Rest : string
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> import numpy as np
            >>> from pylayers.location.algebraic.algebraic import *
            >>> from pylayers.util.geomutil import dist
            >>> nRN = 4
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> c = 0.3
            >>> BN = L*sp.rand(dim,1)
<<<<<<< HEAD
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
        k1 = (np.sum((self.RN - self.RNr) * (self.RN - self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA_std = self.c * self.TDOA_std
        RDOA2 = RDOA * RDOA
        k2 = RDOA2
        K = k1 - k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T - self.RNr.T,
                       RDOA.reshape(np.shape(self.TDOA)[0], 1)))
        # Construct the Covariance Matrix
        C = np.diag(RDOA_std[:] ** 2)
        # Apply LS operator
        Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))),
                          np.dot(np.dot(A.T, nplg.inv(C)), K))
        P = Pr[:shRN[0]].reshape(shRN[0], 1)
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
=======
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN_TOA = L*sp.rand(dim,nRN)
            >>> RN_RSS = L*sp.rand(dim,nRN)
            >>> RN_TDOA = L*sp.rand(dim,nRN)
            
            >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
            >>> TOF = d_TOA/c # actual TOA
            >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std
>>>>>>> mlaaraie-master

            >>> RSS_std = 0.001 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> PL0 = 34.7*np.ones(nRN)
            >>> d0 = 1.
            >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
            >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
            >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
            >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
            >>> d = dist(RN_TDOA,BN,0)
            >>> dr = dist(RNr_TDOA,BN,0)
            >>> TDOF = (d-dr)/c # actual TDOA
            >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
<<<<<<< HEAD
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.twls_locate()
        """
        shRN = np.shape(self.RN)
        # Construct the vector K (see theory)
        k1 = (np.sum((self.RN - self.RNr) * (self.RN - self.RNr), axis=0))
        RDOA = self.c * self.TDOA
        RDOA_std = self.c * self.TDOA_std
        RDOA2 = RDOA * RDOA
        k2 = RDOA2
        K = k1 - k2
        # Construct the matrix A (see theory)
        A = np.hstack((self.RN.T - self.RNr.T,
                       RDOA.reshape(np.shape(self.TDOA)[0], 1)))
        # Construct the Covariance Matrix
        C = np.diag(RDOA_std[:] ** 2)
        # Apply LS operator
        A2 = np.dot(A.T, np.dot(nplg.inv(C), A))
        [U, S, V] = nplg.svd(A2)
        J = 1 / S
        rA = np.rank(A)
        m, n = np.shape(A)
        f = 0
        if np.log10(nplg.cond(A)) >= self.c * max(self.TDOA_std):
            f = f + 1
            for i in range(n - rA):
                u = np.where(J == max(J))
                J[u] = 0

        A2i = np.dot(np.dot(V.T, np.diag(J)), U.T)
        Pr = 0.5 * np.dot(A2i, np.dot(np.dot(A.T, nplg.inv(C)), K))
        P = Pr[:shRN[0]].reshape(shRN[0], 1)
        return P
=======

            >>> nodes={}
            >>> nodes['BN']= BN
            >>> nodes['RN_RSS']= RN_RSS
            >>> nodes['RN_TOA']= RN_TOA
            >>> nodes['RN_TDOA']= RN_TDOA
            >>> nodes['RNr_TDOA']= RNr_TDOA

            >>> ldp={}
            >>> ldp['RSS'] = RSS
            >>> ldp['RSS_std'] = RSS_std
            >>> ldp['RSS_np'] = RSS_np
            >>> ldp['d0'] = d0
            >>> ldp['PL0'] = PL0
            >>> ldp['TOA'] = TOA
            >>> ldp['TOA_std'] = TOA_std
            >>> ldp['TDOA'] = TDOA
            >>> ldp['TDOA_std'] = TDOA_std
            
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.ls_locate(1, 0, 0, 'mode')
            >>> P_toa = S.ls_locate(0, 1, 0, 'mode')
            >>> P_tdoa = S.ls_locate(0, 0, 1, 'mode')
            >>> P_rsstoa = S.ls_locate(1, 1, 0, 'mode')
            >>> P_rsstdoa = S.ls_locate(1, 0, 1, 'mode')
            >>> P_toatdoa = S.ls_locate(0, 1, 1, 'mode')
            >>> P_rsstoatdoa = S.ls_locate(1, 1, 1, 'mode')
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K = k1-k2
                # Construct the matrix A (see theory)
                A = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                # Construct the vector K (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.c * TOA
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = P.reshape(np.shape(RN_TOA[:,0:1]))
                
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                # Construct the vector K (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.get_range(Rest)
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TOA.reshape(np.shape(K_RSS)[0],1)))
                A = np.vstack((A_RSS, A_TOA))
                P = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif rss==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2 = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1 + k2
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1_TOA = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2_TOA = RDOA2
                K_TDOA = k1_TOA-k2_TOA
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_TOA.reshape(np.shape(K_TOA)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_TOA, A_TDOA))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)

            elif toa == 0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_RSS, A_TDOA))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
>>>>>>> mlaaraie-master

            else:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2 = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1 + k2
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1_TOA = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA2 = RDOA*RDOA
                k2_TOA = RDOA2
                K_TDOA = k1_TOA-k2_TOA
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((np.vstack((K_RSS.reshape(np.shape(K_RSS)\
                    [0],1), K_TOA.reshape(np.shape(K_TOA)[0],1))),\
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((np.vstack((A_RSS, A_TOA)),A_TDOA))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,A)), np.dot(A.T,K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            return P


    def wls_locate(self, rss, toa, tdoa, Rest):
        """
        This method applies weighted least squares (WLS) approximation
        to get position P.
        Parameters
        ----------
<<<<<<< HEAD
            P0 : numpy.ndarray
            Niter : integer
=======
            rss : boolean
            toa : boolean
            tdoa : boolean
            Rest : string
>>>>>>> mlaaraie-master
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> import numpy as np
            >>> from pylayers.location.algebraic.algebraic import *
            >>> from pylayers.util.geomutil import dist
            >>> nRN = 4
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> c = 0.3
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN_TOA = L*sp.rand(dim,nRN)
            >>> RN_RSS = L*sp.rand(dim,nRN)
            >>> RN_TDOA = L*sp.rand(dim,nRN)
            
            >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
            >>> TOF = d_TOA/c # actual TOA
            >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std

            >>> RSS_std = 0.001 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> PL0 = 34.7*np.ones(nRN)
            >>> d0 = 1.
            >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
            >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
            >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
            >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
            >>> d = dist(RN_TDOA,BN,0)
            >>> dr = dist(RNr_TDOA,BN,0)
            >>> TDOF = (d-dr)/c # actual TDOA
            >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std

<<<<<<< HEAD
        P = P0
        RDOA = self.c * self.TDOA
        RDOA_std = self.c * self.TDOA_std
        for i in np.arange(Niter):
            # Construct the matrix A (see theory)
            A = ((P - self.RN) / np.sqrt(np.sum((P - self.RN) ** 2, axis=0))).T -\
                ((P - self.RNr) / np.sqrt(
                    np.sum((P - self.RNr) ** 2, axis=0))).T
            # Construct the Covariance Matrix
            C = np.diag((RDOA_std[:]) ** 2)
            # Construct the vector D (see theory)
            D = RDOA - (np.sqrt((np.sum((P - self.RN) ** 2, axis=0))) -
                        np.sqrt((np.sum((P - self.RNr) ** 2, axis=0))))
            # construct the vector Delta (see theory)
            Delta = np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(
                C), A))), np.dot(np.dot(A.T, nplg.inv(C)), D))
            # update P
            P = P + Delta.reshape(np.shape(P))
        return P
=======
            >>> nodes={}
            >>> nodes['BN']= BN
            >>> nodes['RN_RSS']= RN_RSS
            >>> nodes['RN_TOA']= RN_TOA
            >>> nodes['RN_TDOA']= RN_TDOA
            >>> nodes['RNr_TDOA']= RNr_TDOA

            >>> ldp={}
            >>> ldp['RSS'] = RSS
            >>> ldp['RSS_std'] = RSS_std
            >>> ldp['RSS_np'] = RSS_np
            >>> ldp['d0'] = d0
            >>> ldp['PL0'] = PL0
            >>> ldp['TOA'] = TOA
            >>> ldp['TOA_std'] = TOA_std
            >>> ldp['TDOA'] = TDOA
            >>> ldp['TDOA_std'] = TDOA_std
            
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.wls_locate(1, 0, 0, 'mode')
            >>> P_toa = S.wls_locate(0, 1, 0, 'mode')
            >>> P_tdoa = S.wls_locate(0, 0, 1, 'mode')
            >>> P_rsstoa = S.wls_locate(1, 1, 0, 'mode')
            >>> P_rsstdoa = S.wls_locate(1, 0, 1, 'mode')
            >>> P_toatdoa = S.wls_locate(0, 1, 1, 'mode')
            >>> P_rsstoatdoa = S.wls_locate(1, 1, 1, 'mode')
            
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K = k1-k2
                # Construct the matrix A (see theory)
                A = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                RDOA.reshape(np.shape(TDOA)[0],1)))
                # Construct the Covariance Matrix
                C = np.diag(RDOA_std[:] ** 2)
                # Apply LS operator
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C)\
                ,A))),np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                # Construct the vector K (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.c * TOA
                ROA_std = self.c * TOA_std
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Construct the Covariance Matrix
                C = np.diag(ROA_std[1:] ** 2)
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),\
                A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = P.reshape(np.shape(RN_TOA[:,0:1]))
                
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                # Construct the vector K (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1 = RN2[1:] - RN2[0:1]
                ROA = self.get_range(Rest)
                ROA_std = self.get_range_std(Rest)
                ROA2 = ROA * ROA
                k2 = ROA2[0:1] - ROA2[1:]
                K = k1 + k2
                # Construct the matrix A (see theory)
                A = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Construct the Covariance Matrix
                C = np.diag((ROA_std[1:]) ** 2)
                # Apply LS operator
                P = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),\
                A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA_RSS_std = self.get_range_std(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = RN_RSS[:, 1:].T - RN_RSS[:, 0]
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA_TOA_std = self.c * TOA_std
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = RN_TOA[:, 1:].T - RN_TOA[:, 0]
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TOA.reshape(np.shape(K_RSS)[0],1)))
                A = np.vstack((A_RSS, A_TOA))
                C = np.diag(np.hstack((ROA_RSS_std[1:]**2,\
                    ROA_TOA_std[1:]**2)))
                P = 0.5*np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C),\
                A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = P.reshape(np.shape(RN_RSS[:,0:1]))

            elif rss==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA_TOA_std = self.c * TOA_std
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_TOA.reshape(np.shape(K_TOA)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_TOA, A_TDOA))
                C = np.diag(np.hstack((ROA_TOA_std[1:]**2,\
                    RDOA_std[:]**2)))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,np.dot(nplg.inv(C),\
                    A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)

            elif toa == 0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA_RSS_std = self.get_range_std(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((K_RSS.reshape(np.shape(K_RSS)[0],1), \
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((A_RSS, A_TDOA))
                C = np.diag(np.hstack((ROA_RSS_std[1:]**2,\
                    RDOA_std[:]**2)))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,np.dot(nplg.inv(C),\
                    A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
>>>>>>> mlaaraie-master

            else:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                # Construct the vector K_RSS (see theory)
                RN2 = np.sum(RN_RSS * RN_RSS, axis=0)
                k1_RSS = RN2[1:] - RN2[0:1]
                ROA_RSS = self.get_range(Rest)
                ROA_RSS_std = self.get_range_std(Rest)
                ROA2_RSS = ROA_RSS * ROA_RSS
                k2_RSS = ROA2_RSS[0:1] - ROA2_RSS[1:]
                K_RSS = k1_RSS + k2_RSS
                # Construct the matrix A_RSS (see theory)
                A_RSS = np.hstack((RN_RSS[:, 1:].T - RN_RSS[:, 0],\
                        np.zeros((np.shape(RN_RSS)[1]-1,1))))
                # Construct the vector K_TOA (see theory)
                RN2 = (np.sum(RN_TOA * RN_TOA, axis=0))
                k1_TOA = RN2[1:] - RN2[0:1]
                ROA_TOA = self.c * TOA
                ROA_TOA_std = self.c * TOA_std
                ROA2_TOA = ROA_TOA * ROA_TOA
                k2_TOA = ROA2_TOA[0:1] - ROA2_TOA[1:]
                K_TOA = k1_TOA + k2_TOA
                # Construct the matrix A_TOA (see theory)
                A_TOA = np.hstack((RN_TOA[:, 1:].T - RN_TOA[:, 0],\
                        np.zeros((np.shape(RN_TOA)[1]-1,1))))
                # Construct the vector K_TDOA (see theory)
                k1 = (np.sum((RN_TDOA-RNr_TDOA) * (RN_TDOA-RNr_TDOA), \
                      axis=0))
                RDOA = self.c * TDOA
                RDOA_std = self.c * TDOA_std
                RDOA2 = RDOA*RDOA
                k2 = RDOA2
                K_TDOA = k1-k2
                # Construct the matrix A (see theory)
                A_TDOA = np.hstack((RN_TDOA.T-RNr_TDOA.T, \
                        RDOA.reshape(np.shape(TDOA)[0],1)))
                # Apply LS operator
                K = np.vstack((np.vstack((K_RSS.reshape(np.shape(K_RSS)\
                    [0],1), K_TOA.reshape(np.shape(K_TOA)[0],1))),\
                    K_TDOA.reshape(np.shape(K_TDOA)[0],1)))
                A = np.vstack((np.vstack((A_RSS, A_TOA)),A_TDOA))
                C = np.diag(np.hstack((np.hstack((ROA_RSS_std[1:]**2,\
                    ROA_TOA_std[1:]**2)),RDOA_std[:]**2)))
                Pr = 0.5*np.dot(nplg.inv(np.dot(A.T,np.dot(nplg.inv(C),\
                    A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:shRN[0]].reshape(shRN[0],1)
                
            return P

    def ml_function(self, P, rss, toa, tdoa):
        """
        This defines the ML function to be minimized if ML estimator
        is used
        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean
        Returns
        -------
            ML : numpy.ndarray
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                RNnum = shRN[1]
                RDOA = self.c*TDOA
                RDOA_std = self.c*TDOA_std
                # construct the ML function to be minimized
                RN1mP = RN_TDOA - np.outer(P, np.ones(RNnum))
                mRN1mP = np.sqrt(np.diag(np.dot(RN1mP.T, RN1mP)))
                RN2mP = RNr_TDOA - np.outer(P, np.ones(RNnum))
                mRN2mP = np.sqrt(np.diag(np.dot(RN2mP.T, RN2mP)))
                rRDOA = mRN1mP - mRN2mP
                tk = (RDOA - rRDOA) ** 2 / (2 * RDOA_std ** 2)
                ML = tk.sum(axis=0)
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RNnum = np.shape(RN_TOA)[1]
                ROA = self.c * TOA
                ROA_std = self.c * TOA_std
                # construct the ML function to be minimized
                RNmP = RN_TOA - np.outer(P, np.ones(RNnum))
                mRNmP = np.sqrt(np.diag(np.dot(RNmP.T,RNmP)))
                tk = (ROA-mRNmP) ** 2
                uk = tk/(2*ROA_std**2)+np.log(np.sqrt(2*np.pi)*ROA_std)
                ML = uk.sum(axis=0)
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RNnum = np.shape(RN_RSS)[1]
                S = (np.log(10)/10)*RSS_std/RSS_np
                M = (np.log(10)/10)*(PL0-RSS)/RSS_np+np.log(d0)
                # construct the ML function to be minimized
                RNmP = RN_RSS-np.outer(P, np.ones(RNnum))
                mRNmP = np.sqrt(np.diag(np.dot(RNmP.T,RNmP)))
                tk = (M-S**2-np.log(mRNmP))**2
                uk = tk/(2*S**2)
                ML = uk.sum(axis=0)
            elif rss==0:
                ML = self.ml_function(P,0,1,0)+self.ml_function(P,0,0,1)
            elif toa==0:
                ML = self.ml_function(P,1,0,0)+self.ml_function(P,0,0,1)
            elif tdoa==0:
                ML = self.ml_function(P,1,0,0)+self.ml_function(P,0,1,0)
            else:
                ML =self.ml_function(P,1,0,0)+self.ml_function(P,0,1,0)\
                    + self.ml_function(P,0,0,1)

<<<<<<< HEAD
        shRN = np.shape(self.RN)
        RNnum = shRN[1]
        RDOA = self.c * self.TDOA
        RDOA_std = self.c * self.TDOA_std

        # construct the ML function to be minimized
        RN1mP = self.RN - np.outer(P, np.ones(RNnum))
        mRN1mP = np.sqrt(np.diag(np.dot(RN1mP.T, RN1mP)))
        RN2mP = self.RNr - np.outer(P, np.ones(RNnum))
        mRN2mP = np.sqrt(np.diag(np.dot(RN2mP.T, RN2mP)))
        rRDOA = mRN1mP - mRN2mP

        tk = (RDOA - rRDOA) ** 2 / (2 * RDOA_std ** 2)
        ML = tk.sum(axis=0)

        return(ML)
=======
            return ML
>>>>>>> mlaaraie-master

    def ml_locate(self, P0, rss, toa, tdoa):
        """
        This applies ML estimator to compute position P
        Parameters
        ----------
            P0 : numpy.ndarray
<<<<<<< HEAD
=======
            rss : boolean
            toa : boolean
            tdoa : boolean
>>>>>>> mlaaraie-master
        Returns
        -------
            P : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> import numpy as np
            >>> from pylayers.location.algebraic.algebraic import *
            >>> from pylayers.util.geomutil import dist
            >>> nRN = 4
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> c = 0.3
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN_TOA = L*sp.rand(dim,nRN)
            >>> RN_RSS = L*sp.rand(dim,nRN)
            >>> RN_TDOA = L*sp.rand(dim,nRN)
            
            >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
            >>> TOF = d_TOA/c # actual TOA
            >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std

            >>> RSS_std = 0.001 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> PL0 = 34.7*np.ones(nRN)
            >>> d0 = 1.
            >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
            >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
            >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
            >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
            >>> d = dist(RN_TDOA,BN,0)
            >>> dr = dist(RNr_TDOA,BN,0)
            >>> TDOF = (d-dr)/c # actual TDOA
            >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std
<<<<<<< HEAD
            >>> D = Tdoa(RN, RNr, TDOA, TDOA_std)
            >>> D.ml_locate(BN0)
        """
        P = optimize.fmin(self.ml_function, P0, args=(), disp=0)
        return P.reshape(np.shape(P0))
=======
>>>>>>> mlaaraie-master

            >>> nodes={}
            >>> nodes['BN']= BN
            >>> nodes['RN_RSS']= RN_RSS
            >>> nodes['RN_TOA']= RN_TOA
            >>> nodes['RN_TDOA']= RN_TDOA
            >>> nodes['RNr_TDOA']= RNr_TDOA

            >>> ldp={}
            >>> ldp['RSS'] = RSS
            >>> ldp['RSS_std'] = RSS_std
            >>> ldp['RSS_np'] = RSS_np
            >>> ldp['d0'] = d0
            >>> ldp['PL0'] = PL0
            >>> ldp['TOA'] = TOA
            >>> ldp['TOA_std'] = TOA_std
            >>> ldp['TDOA'] = TDOA
            >>> ldp['TDOA_std'] = TDOA_std
            
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.ml_locate(BN0, 1, 0, 0)
            >>> P_toa = S.ml_locate(BN0, 0, 1, 0)
            >>> P_tdoa = S.ml_locate(BN0, 0, 0, 1)
            >>> P_rsstoa = S.ml_locate(BN0, 1, 1, 0)
            >>> P_rsstdoa = S.ml_locate(BN0, 1, 0, 1)
            >>> P_toatdoa = S.ml_locate(BN0, 0, 1, 1)
            >>> P_rsstoatdoa = S.ml_locate(BN0, 1, 1, 1)
        """
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            P = optimize.fmin(self.ml_function,P0,\
                args=(rss,toa,tdoa),disp=0)
            return P.reshape(np.shape(P0))

    def fim(self, P, rss, toa, tdoa):
        """
        Compute the fisher information matrix in P for the given
        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean
        Returns
        -------
            FIM : numpy.ndarray
        """
<<<<<<< HEAD
        RDOA_std = self.c * self.TDOA_std
        FIM = np.zeros((np.shape(self.RN)[0], np.shape(self.RN)[0]))
        for i in range(np.shape(self.RN)[1]):
            PmRN = np.sqrt(np.dot((P - self.RN[:, i:i + 1]).T,
                                  P - self.RN[:, i:i + 1]))
            PmRNr = np.sqrt(np.dot((P - self.RNr[:, i:i + 1]).T,
                                   P - self.RNr[:, i:i + 1]))
            FIM += (1 / RDOA_std[i] ** 2) * np.dot((P - self.RN[:, i:i + 1]) / PmRN -
                                                   (P - self.RNr[:, i:i + 1]) / PmRNr, ((P - self.RN[:, i:i + 1]) / PmRN -
                                                                                        (P - self.RNr[:, i:i + 1]) / PmRNr).T)
        return FIM
=======
        if rss==0 and toa==0 and tdoa==0:
            raise ValueError("inputs missed")
        else:
            if rss==0 and toa==0:
                RN_TDOA = self.nodes['RN_TDOA']
                RNr_TDOA = self.nodes['RNr_TDOA']
                TDOA = self.ldp['TDOA']
                TDOA_std = self.ldp['TDOA_std']
                shRN = np.shape(RN_TDOA)
                RNnum = shRN[1]
                RDOA = self.c*TDOA
                RDOA_std = self.c*TDOA_std
                FIM = np.zeros((np.shape(RN_TDOA)[0],\
                        np.shape(RN_TDOA)[0]))
                for i in range(np.shape(RN_TDOA)[1]):
                    PmRN = np.sqrt(np.dot((P-RN_TDOA[:,i:i+1]).T,\
                        P-RN_TDOA[:,i:i+1]))
                    PmRNr = np.sqrt(np.dot((P-RNr_TDOA[:,i:i+1]).T,\
                        P-RNr_TDOA[:,i:i+1]))
                    FIM += (1/RDOA_std[i]**2)*\
                            np.dot((P-RN_TDOA[:,i:i+1])/PmRN-\
                            (P-RNr_TDOA[:,i:i+1])/PmRNr,\
                            ((P-RN_TDOA[:,i:i+1])/PmRN-\
                            (P-RNr_TDOA[:,i:i+1])/PmRNr).T)
            elif rss==0 and tdoa==0:
                RN_TOA = self.nodes['RN_TOA']
                TOA = self.ldp['TOA']
                TOA_std = self.ldp['TOA_std']
                RNnum = np.shape(RN_TOA)[1]
                ROA = self.c * TOA
                ROA_std = self.c * TOA_std
                FIM = np.zeros((np.shape(RN_TOA)[0],np.shape(RN_TOA)[0]))
                for i in range(np.shape(RN_TOA)[1]):
                    FIM += np.dot((P-RN_TOA[:,i:i+1]),\
                        (P-RN_TOA[:,i:i+1]).T)/((ROA_std[i]**2)*\
                        np.dot((P-RN_TOA[:,i:i+1]).T,\
                        (P-RN_TOA[:,i:i+1])))
            elif toa==0 and tdoa==0:
                RN_RSS = self.nodes['RN_RSS']
                RSS = self.ldp['RSS']
                RSS_std = self.ldp['RSS_std']
                RSS_np = self.ldp['RSS_np']
                d0 = self.ldp['d0']
                PL0 = self.ldp['PL0']
                RNnum = np.shape(RN_RSS)[1]
                S = (np.log(10)/10)*RSS_std/RSS_np
                FIM = np.zeros((np.shape(RN_RSS)[0],\
                    np.shape(RN_RSS)[0]))
                for i in range(np.shape(RN_RSS)[1]):
                    FIM += np.dot((P-RN_RSS[:,i:i+1]),\
                        (P-RN_RSS[:,i:i+1]).T)/((S[0]**2)*\
                        (np.dot((P-RN_RSS[:,i:i+1]).T,\
                        (P-RN_RSS[:,i:i+1])))**2)
            elif rss==0:
                FIM = self.fim(P,0,1,0)+self.fim(P,0,0,1)
            elif toa==0:
                FIM = self.fim(P,1,0,0)+self.fim(P,0,0,1)
            elif tdoa==0:
                FIM = self.fim(P,1,0,0)+self.fim(P,0,1,0)
            else:
                FIM = self.fim(P,1,0,0)+self.fim(P,0,1,0)+\
                    self.fim(P,0,0,1)
>>>>>>> mlaaraie-master

            return FIM

    def crb(self, P, rss, toa, tdoa):
        """
        This method compute the cramer rao bound (CRB) at position P.
        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean
        Returns
        -------
            CRB : float

        Examples
        --------
        .. plot::
            :include-source:
            
            >>> import numpy as np
            >>> from pylayers.location.algebraic.algebraic import *
            >>> from pylayers.util.geomutil import dist
            >>> nRN = 4
            >>> dim = 3 # 2 for 2D, 3 for 3D
            >>> L = 20.
            >>> c = 0.3
            >>> BN = L*sp.rand(dim,1)
            >>> BN0 = L*sp.rand(dim,1)
            >>> RN_TOA = L*sp.rand(dim,nRN)
            >>> RN_RSS = L*sp.rand(dim,nRN)
            >>> RN_TDOA = L*sp.rand(dim,nRN)
            
            >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
            >>> TOF = d_TOA/c # actual TOA
            >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
            >>> TOA = TOF + TOA_std

            >>> RSS_std = 0.001 * np.ones(nRN)
            >>> RSS_np = 2.645 * np.ones(nRN)
            >>> PL0 = 34.7*np.ones(nRN)
            >>> d0 = 1.
            >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
            >>> X = RSS_std * np.random.randn(np.shape(PL0)[0])
            >>> RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    
            >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
            >>> d = dist(RN_TDOA,BN,0)
            >>> dr = dist(RNr_TDOA,BN,0)
            >>> TDOF = (d-dr)/c # actual TDOA
            >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
            >>> TDOA = TDOF + TDOA_std

            >>> nodes={}
            >>> nodes['BN']= BN
            >>> nodes['RN_RSS']= RN_RSS
            >>> nodes['RN_TOA']= RN_TOA
            >>> nodes['RN_TDOA']= RN_TDOA
            >>> nodes['RNr_TDOA']= RNr_TDOA

            >>> ldp={}
            >>> ldp['RSS'] = RSS
            >>> ldp['RSS_std'] = RSS_std
            >>> ldp['RSS_np'] = RSS_np
            >>> ldp['d0'] = d0
            >>> ldp['PL0'] = PL0
            >>> ldp['TOA'] = TOA
            >>> ldp['TOA_std'] = TOA_std
            >>> ldp['TDOA'] = TDOA
            >>> ldp['TDOA_std'] = TDOA_std
            
            >>> S = algloc(nodes, ldp)
            >>> crb_rss = S.crb(BN, 1, 0, 0)
            >>> crb_toa = S.crb(BN, 0, 1, 0)
            >>> crb_tdoa = S.crb(BN, 0, 0, 1)
            >>> crb_rsstoa = S.crb(BN, 1, 1, 0)
            >>> crb_rsstdoa = S.crb(BN, 1, 0, 1)
            >>> crb_toatdoa = S.crb(BN, 0, 1, 1)
            >>> crb_rsstoatdoa = S.crb(BN, 1, 1, 1)
        """

        FIM = self.fim(P, rss, toa, tdoa)
        CRB = np.sqrt(np.trace(nplg.inv(FIM)))
        return CRB

if __name__ == "__main__":
    doctest.testmod()
<<<<<<< HEAD
=======

        
                
                
  
>>>>>>> mlaaraie-master
