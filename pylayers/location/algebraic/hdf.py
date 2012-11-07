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
import cvxmod as cvxm
import cvxopt as cvxo
import matplotlib.pylab as plt
import string


class algloc(object):
    """
    This class regroups all the algebraic localization sceanrios and
    techniques for both hybrid and non-hybrid scenarios using TOA, RSS,
    and TDOA
    Attributes
    ----------
        n_TOA : number of TOAs
        n_TDOA : number of TDOAs
        n_RSS : number of RSSs
        RN_TOA : An Array that defines the Radio nodes implied in
           TOA localization (coordiantes in meters)
        TOA : A measurement vector of TOA associated to RN
            (TOA values in seconds)
        TOA_std : Associated std of TOA (std in seconds)
        RN_RSS : An Array that defines the Radio nodes implied in
           RSS localization (coordiantes in meters)
        RSS : A measurement vector of RSS associated to RN
            (RSS values in dB)
        d0 : the reference distance (usually equal to 1 meter)
        RSS_std : std of shadowing (std in dB)
        RSS_np : propagation constant
        PL0 : RSS at d0
        Rest : range estimator ('mode', 'median', 'mean')
        RN_TDOA : An Array that defines the Radio nodes implied in
           TDOA localization (coordiantes in meters)
        RNr_TDOA : An Array that defines the set of reference radio
           nodes with whom TDOA are computed
        TDOA : A measurement vector of TDoA associated to RN
            (TDoA values in seconds)
        TDOA_std : Associated std of TDoA (std in seconds)
        c : speed of light
          float
    """

    def __init__(self, n_TOA, n_RSS, n_TDOA, RN_TOA, RN_RSS, ):
        if len(np.shape(RN))==2:
            self.RN = RN
        else:
            raise ValueError("inputs must be of shape (n,p), p=2 or 3")
        self.c = 3.e08
