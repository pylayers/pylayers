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
#Bernard UGUEN      : buguen@univ-rennes1.fr
#Mohamed LAARAIEDH  : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################


import numpy as np
import scipy as sp
from scipy import optimize
import numpy.linalg as nplg
import cvxmod as cvxm
import cvxopt as cvxo
import matplotlib.pylab as plt
from string import *
from pylayers.location.algebraic.algebraic import *


if __name__=="__main__":

    nRN = 4
    dim = 2 # 2 for 2D, 3 for 3D
    Niter = 10
    L = 20.
    BN = L*sp.rand(dim,1)
    BN0 = L*sp.rand(dim,1)
    RN_TOA = L*sp.rand(dim,nRN)
    RN_RSS = L*sp.rand(dim,nRN)
    RN_TDOA = L*sp.rand(dim,nRN)
    S = algloc()
    

    # test Toa
    print 'TOA-TOA-TOA'
    d_TOA = S.dist(RN_TOA,BN,0) # actual distances
    TOF = d_TOA/S.c # actual TOA
    TOA_std = 1./S.c*np.ones(np.shape(TOF))
    TOA = TOF + TOA_std
    T = Toa(RN_TOA, TOA, TOA_std)
    print 'Actual position '
    print BN
    print 'LS solution'
    print T.ls_locate()
    print 'TLS solution'
    print T.tls_locate()
    print 'WLS solution'
    print T.wls_locate()
    print 'TWLS solution'
    print T.twls_locate()
    print 'TS solution'
    print T.ts_locate(BN0, Niter)
    print 'SDP solution'
    print T.sdp_locate()
    print 'ML solution'
    print T.ml_locate(BN0)
    print 'CRB'
    print T.crb(BN)
    if TOA_std.all()==0.:
        assert (T.ls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.tls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.wls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.twls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.ts_locate(BN0, Niter).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.sdp_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.ml_locate(BN0).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"

    # test Rss
    print 'RSS-RSS-RSS'
    RSS_std = 4.34 * np.ones(nRN)
    RSS_np = 2.645 * np.ones(nRN)
    RSS = 110*sp.rand(nRN)
    PL0 = -34.7*np.ones(nRN)
    lamda = 5e9
    d0 = 1.
    Rest = 'mode' # RSS based ranging estimator
    R = Rss(RN_RSS, RSS, d0, RSS_std, RSS_np, PL0, Rest)
    PL0 = R.get_pl0(lamda)
    print "theoretic PL0"
    print PL0
    PLmean = R.get_plmean(BN)
    print "mean PL"
    print PLmean
    PL = R.get_pl(BN)
    print "PL"
    print PL
    R.RSS = PL
    Range = R.get_range()
    Range_std = R.get_range_std()
    print 'Estimated Ranges'
    print Range
    print 'STD of Estimated Ranges'
    print Range_std
    print 'LS solution'
    print R.ls_locate()
    print 'TLS solution'
    print R.tls_locate()
    print 'WLS solution'
    print R.wls_locate()
    print 'TWLS solution'
    print R.twls_locate()
    print 'TS solution'
    print R.ts_locate(BN0, Niter)
    print 'SDP solution'
    print R.sdp_locate()
    print 'indirect ML solution'
    print R.iml_locate(BN0)
    print 'direct ML solution'
    print R.dml_locate(BN0)
    print 'CRB'
    print R.crb(BN)

    if RSS_std.all()==0.:
        assert (R.ls_locate().all() == BN.all()),\
         "Without any error, actual and estimated positions are equal"
        assert (R.tls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (R.wls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (R.twls_locate().all() == BN.all()),\
         "Without any error, actual and estimated positions are equal"
        assert (R.ts_locate(BN0, Niter).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (R.sdp_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (R.iml_locate(BN0).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (R.dml_locate(BN0).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"

        
    # test Tdoa
    RNr = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
    d = S.dist(RN_TDOA,BN,0)
    dr = S.dist(RNr,BN,0)
    TDOF = (d-dr)/S.c # actual TDOA
    TDOA_std = 0.001/S.c*np.ones(np.shape(TDOF))
    TDOA = TDOF + TDOA_std
    D = Tdoa(RN_TDOA, RNr, TDOA, TDOA_std)
    print 'Actual position '
    print BN
    print 'LS solution'
    print D.ls_locate()
    print 'TLS solution'
    print D.tls_locate()
    print 'WLS solution'
    print D.wls_locate()
    print 'TWLS solution'
    print D.twls_locate()
    print 'TS solution'
    print D.ts_locate(BN0, Niter)
    print 'ML solution'
    print D.ml_locate(BN0)
    print 'CRB'
    print D.crb(BN)
    if TDOA_std.all()==0.:
        assert (D.ls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (D.tls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (D.wls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (D.twls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (D.ts_locate(BN0, Niter).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (D.ml_locate(BN0).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"

    # test Toarss
    print 'TOA/RSS'
    print np.shape(PL0)
    TR = Toarss(RN_TOA, TOA, TOA_std, RN_RSS, R.RSS, d0, RSS_std, \
    RSS_np, R.PL0, Rest)
    print 'Actual position '
    print BN
    print 'LS solution'
    print TR.ls_locate()
    print 'WLS solution'
    print TR.wls_locate()
    print 'ML solution'
    print TR.ml_locate(BN0)
    '''
    print 'CRB'
    print T.crb(BN)
    if TOA_std.all()==0.:
        assert (T.ls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.tls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.wls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.twls_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.ts_locate(BN0, Niter).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.sdp_locate().all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
        assert (T.ml_locate(BN0).all() == BN.all()), \
        "Without any error, actual and estimated positions are equal"
    '''
    
    
