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
#Bernard UGUEN          : buguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import os
from numpy import *
from scipy import *
from scipy import optimize
from numpy.linalg import *
import cvxmod as cvxm
import cvxopt as cvxo
from pylayers.location.algebraic.crb import *


class TDoALocation(object):
    """
    A TDoALocation contains:
    1- 2 sets of RadioNodes (RN1 and RN2) with associated position accuracies (RN1QoS and RN2QoS),
    2- a set of TDoAs measurements (TDoA) with associated STD (TDoAStd)

    This class manages the TDoA-based localization techniques.
    MEMBERS:

            RN1     : An Array that defines the set of first side Radio nodes implied in localization (coordiantes in meters)
                    : shape(RN1)= (2 or 3,RNnum)
            RN1QoS  : An Array that defines the precision of positions of RN1 (std in meters)
                    : shape(RN1QoS)= (2 or 3, RNnum)
            RN2     : An Array that defines the set of second side Radio nodes implied in localization (coordiantes in meters)
                    : shape(RN2)= (2 or 3,RNnum)
            RN2QoS  : An Array that defines the precision of positions of RN2 (std in meters)
                    : shape(RN2QoS)= (2 or 3, RNnum)
            TDoA    : A measurement vector of TDoA associated to RN (TDoA values in seconds)
                    : shape(TDoA)= (RNnum,1)
            TDoAStd : Associated STD of TDoA (std in seconds)
                    : shape(TDoAStd)= (RNnum,1)



    Provided Methods:

            info()                                          : Display information about the scenario

            LSTDoALocate(RN, TDoA)                          : Applies Least Square approximation and estimate position
            WLSTDoALocate(RN, TDoA, TDoAStd)                : Applies Weighted Least Square approximation and estimate position
            TSTDoALocation(P0, RN, TDoA, TDoAStd, Niter)    : Applies Taylor Series method and estimate position after Niter iterations

            TDoAOptimizer(RN, TDoA, TDoAStd)                : Defines the function to be optimized
            MLTDoALocate(P0, RN, TDoA, TDoAStd)             : Optimize TDoAOptimizer() and estimate Position (P0:initial guess)

            CRBTDoALocate(self, P, RN, TDoA, TDoAStd)       : Compute the CRB in P for the given scenario
    """
    """
    def __init__(self,RN1, RN2, TDoA, TDoAStd):
            self.RN1        = RN1
            self.RN2        = RN2
            self.TDoA       = TDoA
            self.TDoAStd    = TDoAStd
    """
    def __init__(self, RN1):
        self.RN1 = RN1

    def info(self):
        """
        Dispaly scenario information
        """
        print "First Reference Radio Nodes:\n", self.RN1
        print "Second Reference Radio Nodes:\n", self.RN2
        print "Measured TDoA:\n", self.TDoA
        print "STD of Measured TDoA:\n", self.TDoAStd

    def LSTDoALocate(self, RN1, RN2, TDoA):
        """
        This applies LS approximation on TDoA to get position P.
        Return P
        """
        shRN = shape(RN1)                                    # shape of RN
        RNnum = shRN[1]
            # Number of reference nodes
        c = 3e08                                          # Speed of light
        # Construct the vector K (see theory)
        k1 = (sum((RN1 - RN2) * (RN1 - RN2), axis=0)).reshape(
            RNnum, 1)    # first half of K

        RDoA = c * TDoA                                        # Range of arrival (meters)
        RDoA2 = (RDoA * RDoA).reshape(RNnum, 1)
        k2 = RDoA2                                         # second half of K

        K = k1 - k2

        # Construct the matrix A (see theory)
        A = hstack((RN1.T - RN2.T, RDoA))

        # Apply LS operator
        Pr = 0.5 * dot(linalg.inv(dot(A.T, A)), dot(A.T, K))
        P = Pr[:shRN[0], :]
        # Return the estimated position
        return P

    def TLSTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
        """
        This applies LS approximation on TDoA to get position P.
        Return P
        """
        shRN = shape(RN1)                                    # shape of RN
        RNnum = shRN[1]
            # Number of reference nodes
        c = 3e08                                          # Speed of light
        # Construct the vector K (see theory)
        k1 = (sum((RN1 - RN2) * (RN1 - RN2), axis=0)).reshape(
            RNnum, 1)    # first half of K

        RDoA = c * TDoA                                        # Range of arrival (meters)
        RDoA2 = (RDoA * RDoA).reshape(RNnum, 1)
        k2 = RDoA2                                         # second half of K

        K = k1 - k2

        # Construct the matrix A (see theory)
        A = hstack((RN1.T - RN2.T, RDoA))
        A2 = dot(transpose(A), A)

        [U, S, V] = svd(A2)
        J = 1 / S
        rA = rank(A)
        m, n = shape(A)
        f = 0
        if log10(cond(A)) >= c * max(TDoAStd):
            f = f + 1
            for i in range(n - rA):
                u = where(J == max(J))
                J[u] = 0

        A2i = dot(dot(V.T, diag(J)), U.T)
        # Apply LS operator
        Pr = 0.5 * dot(A2i, dot(A.T, K))
        P = Pr[:shRN[0], :]
        # Return the estimated position
        return P

    def WLSTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
        """
        This applies WLS approximation on TDoA assuming TDoAStd to get position P.
        Return P
        """
        shRN = shape(RN1)                                    # shape of RN
        RNnum = shRN[1]
            # Number of reference nodes
        c = 3e08
        RDoAStd = c * TDoAStd                                             # Speed of light
        # Construct the vector K (see theory)
        k1 = (sum((RN1 - RN2) * (RN1 - RN2), axis=0)).reshape(
            RNnum, 1)    # first half of K

        RDoA = c * TDoA                                        # Range of arrival (meters)
        RDoA2 = (RDoA * RDoA).reshape(RNnum, 1)
        k2 = RDoA2                                         # second half of K

        K = k1 - k2

        # Construct the matrix A (see theory)
        A = hstack((RN1.T - RN2.T, RDoA))

        # Construct the Covariance Matrix
        C = diag(RDoAStd[:, 0] ** 2)

        # Apply LS operator
        Pr = 0.5 * dot(linalg.inv(dot(A.T, dot(linalg.inv(C),
                                               A))), dot(dot(A.T, linalg.inv(C)), K))
        P = Pr[:shRN[0], :]
        # Return the estimated position
        return P

    def TWLSTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
        """
        This applies WLS approximation on TDoA assuming TDoAStd to get position P.
        Return P
        """
        shRN = shape(RN1)                                    # shape of RN
        RNnum = shRN[1]
            # Number of reference nodes
        c = 3e08
        RDoAStd = c * TDoAStd                                             # Speed of light
        # Construct the vector K (see theory)
        k1 = (sum((RN1 - RN2) * (RN1 - RN2), axis=0)).reshape(
            RNnum, 1)    # first half of K

        RDoA = c * TDoA                                        # Range of arrival (meters)
        RDoA2 = (RDoA * RDoA).reshape(RNnum, 1)
        k2 = RDoA2                                         # second half of K

        K = k1 - k2

        # Construct the matrix A (see theory)
        A = hstack((RN1.T - RN2.T, RDoA))

        # Construct the Covariance Matrix
        C = diag(RDoAStd[:, 0] ** 2)

        A2 = dot(A.T, dot(linalg.inv(C), A))

        [U, S, V] = svd(A2)
        J = 1 / S
        rA = rank(A)
        m, n = shape(A)
        f = 0
        if log10(cond(A)) >= c * max(TDoAStd):
            f = f + 1
            for i in range(n - rA):
                u = where(J == max(J))
                J[u] = 0

        A2i = dot(dot(V.T, diag(J)), U.T)
        # Apply LS operator
        Pr = 0.5 * dot(A2i, dot(dot(A.T, linalg.inv(C)), K))
        P = Pr[:shRN[0], :]
        # Return the estimated position
        return P

    def TSTDoALocation(self, P0, RN1, RN2, TDoA, TDoAStd, Niter):
        '''
        Applies Taylor Series method and estimate position after Niter iterations
        '''

        P = P0                                            # Initialisation of P as equal to intial guess P0
        shRN = shape(RN1)                                    # shape of RN
        RNnum = shRN[1]
            # Number of reference nodes
        c = 3e08                                          # Speed of light
        RDoA = c * TDoA
        RDoAStd = c * TDoAStd

        for i in arange(Niter):

            # Construct the matrix A (see theory)
            A = ((outer(P, ones(RNnum)) - RN1) / sqrt(sum((outer(P, ones(RNnum)) - RN1) ** 2, axis=0))).T - ((outer(P, ones(RNnum)) - RN2) / sqrt(sum((outer(P, ones(RNnum)) - RN2) ** 2, axis=0))).T

            # Construct the Covariance Matrix
            C = diag((RDoAStd[:, 0]) ** 2)

            # Construct the vector D (see theory)
            D = RDoA - (sqrt((sum((outer(P, ones(RNnum)) - RN1) ** 2, axis=0)).reshape(shape(RDoA))) - sqrt((sum((outer(P, ones(RNnum)) - RN2) ** 2, axis=0)).reshape(shape(RDoA))))

            # construct the vector Delta (see theory)
            Delta = dot(linalg.inv(dot(A.T, dot(linalg.inv(
                C), A))), dot(dot(A.T, linalg.inv(C)), D))

            # update P
            P = P + Delta

        # Return the estimated position
        return P

    def TDoAOptimizer(self, P, RN1, RN2, TDoA, TDoAStd):
        """
        This defines the ML function to be minimized
        """

        shRN = shape(RN1)                                    # shape of RN
        RNnum = shRN[1]
            # Number of reference nodes
        c = 3e08                                          # Speed of light
        RDoA = c * TDoA
        RDoAStd = c * TDoAStd

        # construct the ML function to be minimized
        RN1mP = RN1 - outer(P, ones(RNnum))
        mRN1mP = (sqrt(diag(dot(RN1mP.T, RN1mP)))).reshape(RNnum, 1)
        RN2mP = RN2 - outer(P, ones(RNnum))
        mRN2mP = (sqrt(diag(dot(RN2mP.T, RN2mP)))).reshape(RNnum, 1)
        rRDoA = mRN1mP - mRN2mP

        tk = (RDoA - rRDoA) ** 2 / (2 * RDoAStd ** 2)
        uk = tk[:, 0]  # *(RN1mP/mRN1mP[:,0]-RN2mP/mRN2mP[:,0])
        suk = uk.sum(axis=0)
        #msuk    = sqrt(dot(suk,suk.T))

        return(suk)

    def MLTDoALocate(self, P, P0, RN1, RN2, TDoA, TDoAStd):
        """
        Optimization Routine
        """

        P = optimize.fmin(self.TDoAOptimizer, P0, args=(RN1,
                                                        RN2, TDoA, TDoAStd), xtol=1e-10, ftol=1e-10)

        return P.reshape(shape(P0))

    def SDPTDoALocate(self, RN1, RN2, TDoA, TDoAStd):
        """
        Apply SDP approximation and localization
        """
        RN1 = cvxm.matrix(RN1)
        RN2 = cvxm.matrix(RN2)
        TDoA = cvxm.matrix(TDoA)
        c = 3e08
        RDoA = c * TDoA
        RDoAStd = cvxm.matrix(c * TDoAStd)
        mtdoa, ntdoa = cvxm.size(RN1)
        Im = cvxm.eye(mtdoa)
        Y = cvxm.optvar('Y', mtdoa + 1, mtdoa + 1)
        t = cvxm.optvar('t', ntdoa, 1)
        prob = cvxm.problem(cvxm.minimize(cvxm.norm2(t)))
        prob.constr.append(Y >= 0)
        prob.constr.append(Y[mtdoa, mtdoa] == 1)
        for i in range(ntdoa):
            X0 = cvxm.matrix([[Im, -cvxm.transpose(RN1[:, i])], [-RN1[
                :, i], cvxm.transpose(RN1[:, i]) * RN1[:, i]]])
            X1 = cvxm.matrix([[Im, -cvxm.transpose(RN2[:, i])], [-RN2[
                :, i], cvxm.transpose(RN2[:, i]) * RN2[:, i]]])
            prob.constr.append(-RDoAStd[i, 0] * t[i] < cvxm.trace(
                X0 * Y) + cvxm.trace(X1 * Y) - RDoA[i, 0] ** 2)
            prob.constr.append(RDoAStd[i, 0] * t[i] > cvxm.trace(
                X0 * Y) + cvxm.trace(X1 * Y) - RDoA[i, 0] ** 2)
        prob.solve()
        Pval = Y.value
        X_cvx = Pval[:2, -1]

        return X_cvx

    def CRBTDoALocate(self, P, RN1, RN2, TDoAStd):
        """
        Compute the CRB in P for the given scenario
        """

        '''c       = 3e08
        shP     = shape(P)
        shRN    = shape(RN1)
        RNnum   = shRN[1]
        RDoAStd = c*TDoAStd

        RN1mP   = outer(P,ones(RNnum))- RN1
        mRN1mP  = (sqrt(diag(dot(RN1mP.T,RN1mP)))).reshape(RNnum,1)
        RN2mP   = outer(P,ones(RNnum))- RN2
        mRN2mP  = (sqrt(diag(dot(RN2mP.T,RN2mP)))).reshape(RNnum,1)

        num     = sum(2/(RDoAStd[:,0]**2)*(1-sum((RN1mP/mRN1mP[:,0])*(RN2mP/mRN2mP[:,0]),axis=0)),axis=0)       # the numerator of the CRLB


        div1    = sum((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]-RN2mP/ \
            mRN2mP[:,0])**2,axis=1).reshape(shP)
        don1    = div1.prod(axis=0)[0]          # first term of the doniminator

        div2    = sum(prod((1/RDoAStd[:,0]**2)*(RN1mP/mRN1mP[:,0]- \
            RN2mP/mRN2mP[:,0]),axis=0),axis=0)
        don2    = div2**2                       # second term of the doniminator

        CRB     = num/(don1-don2)               # the CRB'''
        crlb = CRBLocation(RN1)
        CRB = crlb.CRB_TDOA_fim(P, RN1, RN2, TDoAStd)

        return sqrt(abs(CRB))
