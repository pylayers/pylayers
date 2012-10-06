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
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import os
from numpy import *
import numpy as np 
import scipy as sp
from scipy import *
from scipy import optimize
from string import *
from numpy.linalg import *


class CRBLocation(object):
    """
    A CRBLocation contains:
    1- a set of RadioNodes (RN) with associated position accuracies (RNQoS),
    2- a set of measurements (RSS, TOA, TDOA) with associated accuracies.

    This class manages the CRB techniques.
    MEMBERS:

            RN      : An np.array that defines the Radio nodes implied in localization (coordiantes in meters)
                    : shape(RN)= (2 or 3,RNnum)
            RNQoS   : An np.array that defines the precision of positions of RN (std in meters)
                    : shape(RNQoS)= (2 or 3, RNnum)
            param   : a set of parameters depending on the type of measurement [list].


    Provided Methods:

            info()                                                  : Display information about the scenario


    """

    def __init__(self, RN):
        self.RN = RN

    def info(self):
        """
        Dispaly scenario information
        """
        print "Reference Radio Nodes:\n", self.RN
        print "parameters:\n", self.param

    def FIM_RSS(self, P, RN_RSS, RSSnp, RSSStd):
        """
        This computes FIM in P of RSS positioning
        """

        shP = shape(P)
        shRN = shape(RN_RSS)
        RNnum = shRN[1]

        S = (log(10) / 10) * RSSStd / RSSnp

        RNmP = RN_RSS - np.outer(P, ones(RNnum))
        mRNmP = (sqrt(diag(dot(RNmP.T, RNmP))))

        j11 = sum(((1 + S[:, 0] ** 2) * RNmP[0, :]
                   ** 2) / ((S[:, 0] ** 2) * mRNmP ** 4), axis=0)
        j22 = sum(((1 + S[:, 0] ** 2) * RNmP[1, :]
                   ** 2) / ((S[:, 0] ** 2) * mRNmP ** 4), axis=0)
        j12 = j21 = sum(((1 + S[:, 0] ** 2) * (RNmP.prod(
            axis=0))) / ((S[:, 0] ** 2) * mRNmP ** 4), axis=0)

        FIM = np.array([[j11, j12], [j21, j22]])
        return FIM

    def FIM_TOA(self, P, RN_TOA, TOAStd):
        """
        Compute the FIM in P for the given scenario
        """

        c = 3e08
        shP = shape(P)
        shRN = shape(RN_TOA)
        RNnum = shRN[1]
        RoAStd = c * TOAStd

        num = sum(
            1 / (RoAStd ** 2), axis=0)[0]  # the numerator of the CRLB

        RNmP = RN_TOA - np.outer(P, ones(RNnum))
        mRNmP = (sqrt(diag(dot(RNmP.T, RNmP))))

        j11 = sum(RNmP[0, :] ** 2 / ((RoAStd[:, 0]
                                      ** 2) * mRNmP ** 2), axis=0)
        j22 = sum(RNmP[1, :] ** 2 / ((RoAStd[:, 0]
                                      ** 2) * mRNmP ** 2), axis=0)
        j12 = j21 = sum((RNmP.prod(
            axis=0)) / ((RoAStd[:, 0] ** 2) * mRNmP ** 2), axis=0)

        FIM = np.array([[j11, j12], [j21, j22]])
        return FIM

    def FIM_TDOA(self, P, RN1_TDOA, RN2_TDOA, TDOAStd):
        """
        Compute the FIM in P for the given scenario
        """

        c = 3e08
        shP = shape(P)
        shRN = shape(RN1_TDOA)
        RNnum = shRN[1]
        RDoAStd = c * TDOAStd

        RN1mP = np.outer(P, ones(RNnum)) - RN1_TDOA
        mRN1mP = (sqrt(diag(dot(RN1mP.T, RN1mP)))).reshape(RNnum, 1)
        RN2mP = np.outer(P, ones(RNnum)) - RN2_TDOA
        mRN2mP = (sqrt(diag(dot(RN2mP.T, RN2mP)))).reshape(RNnum, 1)

        j11 = sum((1 / RDoAStd[:, 0] ** 2) * (RN1mP /
                                              mRN1mP[:, 0] - RN2mP / mRN2mP[:, 0]) ** 2, axis=1)[0]
        j22 = sum((1 / RDoAStd[:, 0] ** 2) * (RN1mP /
                                              mRN1mP[:, 0] - RN2mP / mRN2mP[:, 0]) ** 2, axis=1)[1]
        j12 = j21 = sum(prod((1 / RDoAStd[:, 0] ** 2) * (RN1mP /
                                                         mRN1mP[:, 0] - RN2mP / mRN2mP[:, 0]), axis=0), axis=0)

        FIM = np.array([[j11, j12], [j21, j22]])
        return FIM

    def CRB_RSS_fim(self, P, RN_RSS, RSSnp, RSSStd):
        """
        This computes CRB of RSS positioning as the trace of inv of fim
        """
        FIM = self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd)

        return trace(inv(FIM))

    def CRB_TOA_fim(self, P, RN_TOA, TOAStd):
        """
        Compute the CRB in P for the given scenario  as the trace of inv of fim
        """
        FIM = self.FIM_TOA(P, RN_TOA, TOAStd)

        return trace(inv(FIM))

    def CRB_TDOA_fim(self, P, RN1_TDOA, RN2_TDOA, TDOAStd):
        """
        Compute the CRB in P for the given scenario as the trace of inv of fim
        """
        FIM = self.FIM_TDOA(P, RN1_TDOA, RN2_TDOA, TDOAStd)

        return trace(inv(FIM))

    def CRB_RSS_TOA_fim(self, P, RN_RSS, RN_TOA, RSSnp, RSSStd, TOAStd):
        """
        This computes CRB of RSS/TOA positioning as the trace of inv of fim
        """
        FIM = self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd) + \
            self.FIM_TOA(P, RN_TOA, TOAStd)

        return trace(inv(FIM))

    def CRB_RSS_TDOA_fim(self, P, RN_RSS, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TDOAStd):
        """
        This computes CRB of RSS/TDOA positioning as the trace of inv of fim
        """
        FIM = self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd) + \
            self.FIM_TDOA(P, RN1_TDOA, RN2_TDOA, TDOAStd)

        return trace(inv(FIM))

    def CRB_TOA_TDOA_fim(self, P, RN_TOA, RN1_TDOA, RN2_TDOA, TOAStd, TDOAStd):
        """
        This computes CRB of TOA/TDOA positioning as the trace of inv of fim
        """
        FIM = self.FIM_TOA(P, RN_TOA, TOAStd) + self.FIM_TDOA(P,
                                                              RN1_TDOA, RN2_TDOA, TDOAStd)

        return trace(inv(FIM))

    def CRB_RSS_TOA_TDOA_fim(self, P, RN_RSS, RN_TOA, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TOAStd, TDOAStd):
        """
        This computes CRB of RSS/TOA/TDOA positioning as the trace of inv of fim
        """
        FIM = self.FIM_RSS(P, RN_RSS, RSSnp, RSSStd) + self.FIM_TOA(P, RN_TOA,
                                                                    TOAStd) + self.FIM_TDOA(P, RN1_TDOA, RN2_TDOA, TDOAStd)

        return trace(inv(FIM))

    def MCRB_RSS_fim(self, L, RN_RSS, RSSnp, RSSStd):
        """
        This computes mean CRB of RSS positioning over the area of length L
        """

        delta = L / 50.0
        CRB = []
        for x in arange(0.001, L + 0.1, delta):

            for y in arange(0.001, L + 0.1, delta):

                P = np.array([[x], [y]])
                f1 = self.CRB_RSS_fim(P, RN_RSS, RSSnp, RSSStd)
                if isnan(f1) == 0:
                    CRB.append(f1)

        moy = mean(CRB)
        return moy

    def MCRB_TOA_fim(self, L, RN_TOA, TOAStd):
        """
        This computes mean CRB of TOA positioning over the area of length L
        """

        delta = L / 50.0
        CRB = []
        for x in arange(0.001, L + 0.1, delta):

            for y in arange(0.001, L + 0.1, delta):

                P = np.array([[x], [y]])
                f1 = self.CRB_TOA_fim(P, RN_TOA, TOAStd)
                if isnan(f1) == 0:
                    CRB.append(f1)

        moy = mean(CRB)
        return moy

    def MCRB_TDOA_fim(self, L, RN1_TDOA, RN2_TDOA, TDOAStd):
        """
        This computes mean CRB of TDOA positioning over the area of length L
        """

        delta = L / 50.0
        CRB = []
        for x in arange(0.001, L + 0.1, delta):

            for y in arange(0.001, L + 0.1, delta):

                P = np.array([[x], [y]])
                f1 = self.CRB_TDOA_fim(P, RN1_TDOA, RN2_TDOA, TDOAStd)
                if isnan(f1) == 0:
                    CRB.append(f1)

        moy = mean(CRB)
        return moy

    def MCRB_RSS_TOA_fim(self, L, RN_RSS, RN_TOA, RSSnp, RSSStd, TOAStd):
        """
        This computes mean CRB of RSS/TOA positioning over the area of length L
        """

        delta = L / 50.0
        CRB = []
        for x in arange(0.001, L + 0.1, delta):

            for y in arange(0.001, L + 0.1, delta):

                P = np.array([[x], [y]])
                f1 = self.CRB_RSS_TOA_fim(
                    P, RN_RSS, RN_TOA, RSSnp, RSSStd, TOAStd)
                if isnan(f1) == 0:
                    CRB.append(f1)

        moy = mean(CRB)
        return moy

    def MCRB_RSS_TDOA_fim(self, L, RN_RSS, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TDOAStd):
        """
        This computes mean CRB of RSS/TDOA positioning over the area of length L
        """

        delta = L / 50.0
        CRB = []
        for x in arange(0.001, L + 0.1, delta):

            for y in arange(0.001, L + 0.1, delta):

                P = np.array([[x], [y]])
                f1 = self.CRB_RSS_TDOA_fim(P, RN_RSS, RN1_TDOA,
                                           RN2_TDOA, RSSnp, RSSStd, TDOAStd)
                if isnan(f1) == 0:
                    CRB.append(f1)

        moy = mean(CRB)
        return moy

    def MCRB_TOA_TDOA_fim(self, L, RN_TOA, RN1_TDOA, RN2_TDOA, TOAStd, TDOAStd):
        """
        This computes mean CRB of TOA/TDOA positioning over the area of length L
        """

        delta = L / 50.0
        CRB = []
        for x in arange(0.001, L + 0.1, delta):

            for y in arange(0.001, L + 0.1, delta):

                P = np.array([[x], [y]])
                f1 = self.CRB_TOA_TDOA_fim(P, RN_TOA,
                                           RN1_TDOA, RN2_TDOA, TOAStd, TDOAStd)
                if isnan(f1) == 0:
                    CRB.append(f1)

        moy = mean(CRB)
        return moy

    def MCRB_RSS_TOA_TDOA_fim(self, L, RN_RSS, RN_TOA, RN1_TDOA, RN2_TDOA, RSSnp, RSSStd, TOAStd, TDOAStd):
        """
        This computes mean CRB of RSS/TOA/TDOA positioning over the area of length L
        """

        delta = L / 50.0
        CRB = []
        for x in arange(0.001, L + 0.1, delta):

            for y in arange(0.001, L + 0.1, delta):

                P = np.array([[x], [y]])
                f1 = self.CRB_RSS_TOA_TDOA_fim(P, RN_RSS, RN_TOA, RN1_TDOA,
                                               RN2_TDOA, RSSnp, RSSStd, TOAStd, TDOAStd)
                if isnan(f1) == 0:
                    CRB.append(f1)

        moy = mean(CRB)
        return moy

    def CRB_RSS(self, P, RN_RSS, RSSnp, RSSStd):
        """
        This computes CRB of RSS positioning
        """

        shP = shape(P)
        shRN = shape(RN_RSS)
        RNnum = shRN[1]

        S = (log(10) / 10) * RSSStd / RSSnp

        RNmP = RN_RSS - np.outer(P, ones(RNnum))
        mRNmP = (sqrt(diag(dot(RNmP.T, RNmP))))

        num = sum((1 + S[:, 0] ** 2) / ((S[:, 0] ** 2) *
                                        mRNmP ** 2), axis=0)      # the numerator of the CRLB

        div1 = sum(((1 + S[:, 0] ** 2) * RNmP ** 2) / ((S[:, 0]
                                                        ** 2) * mRNmP ** 4), axis=1).reshape(shP)
        don1 = div1.prod(axis=0)  # [0]         # first term of the doniminator

        div2 = sum(((1 + S[:, 0] ** 2) * (RNmP.prod(
            axis=0))) / ((S[:, 0] ** 2) * mRNmP ** 4), axis=0)
        don2 = div2 ** 2                       # second term of the doniminator

        CRB = num / (don1 - don2)               # the CRB

        return CRB

    def CRB_TOA(self, P, RN_TOA, TOAStd):
        """
        Compute the CRB in P for the given scenario
        """

        c = 3e08
        shP = shape(P)
        shRN = shape(RN_TOA)
        RNnum = shRN[1]
        RoAStd = c * TOAStd

        num = sum(
            1 / (RoAStd ** 2), axis=0)[0]  # the numerator of the CRLB

        RNmP = RN_TOA - np.outer(P, ones(RNnum))
        mRNmP = (sqrt(diag(dot(RNmP.T, RNmP))))

        div1 = sum(RNmP ** 2 / ((RoAStd[:, 0] ** 2) *
                                mRNmP ** 2), axis=1).reshape(shP)
        don1 = div1.prod(axis=0)[0]          # first term of the doniminator

        div2 = sum((RNmP.prod(
            axis=0)) / ((RoAStd[:, 0] ** 2) * mRNmP ** 2), axis=0)
        don2 = div2 ** 2                       # second term of the doniminator

        CRB = num / (don1 - don2)               # the CRB

        return CRB

    def CRB_TDOA(self, P, RN1_TDOA, RN2_TDOA, TDOAStd):
        """
        Compute the CRB in P for the given scenario

        Parameters
        ----------
        P
        RNA_TDOA
        RN2_TDOA
        TDOAStd

        """

        c = 3e08
        shP = shape(P)
        shRN = shape(RN1_TDOA)
        RNnum = shRN[1]
        RDoAStd = c * TDOAStd

        RN1mP = np.outer(P, ones(RNnum)) - RN1_TDOA
        mRN1mP = (sqrt(diag(dot(RN1mP.T, RN1mP)))).reshape(RNnum, 1)
        RN2mP = np.outer(P, ones(RNnum)) - RN2_TDOA
        mRN2mP = (sqrt(diag(dot(RN2mP.T, RN2mP)))).reshape(RNnum, 1)

        num = sum(2 / (RDoAStd[:, 0] ** 2) * 
                  (1 - sum((RN1mP / mRN1mP[:, 0]) *
                           (RN2mP / mRN2mP[:, 0]), axis=0)), axis=0)       # the numerator of the CRLB

        div1 = sum((1 / RDoAStd[:, 0] ** 2) * 
                   (RN1mP / mRN1mP[:, 0] - RN2mP / mRN2mP[:, 0]) ** 2, axis=1).reshape(shP)
        don1 = div1.prod(axis=0)[0]          # first term of the doniminator

        div2 = sum(prod((1 / RDoAStd[:, 0] ** 2) * (RN1mP /
                                                    mRN1mP[:, 0] - RN2mP / mRN2mP[:, 0]), axis=0), axis=0)
        don2 = div2 ** 2                       # second term of the doniminator

        CRB = num / (don1 - don2)               # the CRB

        return CRB

    def P_CRB_RSS(self, P0, RN_RSS, RSSnp, RSSStd):
        """
        This uses RSS CRB to compute the best position of additional reference node
        """

        #P0=np.array(([[0.0],[0.0]])
        P = optimize.fmin(self.CRB_RSS, P0, args=(RN_RSS, RSSnp,
                                                  RSSStd), xtol=1e-10, ftol=1e-10)

        return P

    def P_CRB_TOA(self, P0, RN_TOA, TOAStd):
        """
        This uses TOA CRB to compute the best position of additional reference node
        """

        #P0=np.array(([[0.0],[0.0]])
        P = optimize.fmin(self.CRB_TOA, P0, args=(RN_TOA, TOAStd),
                          xtol=1e-10, ftol=1e-10)

        return P
