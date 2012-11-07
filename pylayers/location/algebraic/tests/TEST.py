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
from pylayers.location.algebraic.algebraic import *
from pylayers.util.geomutil import dist


if __name__=="__main__":

    nRN = 4
    dim = 3 # 2 for 2D, 3 for 3D
    L = 20.
    c = 0.3
    BN = L*sp.rand(dim,1)
    BN0 = L*sp.rand(dim,1)
    RN_TOA = L*sp.rand(dim,nRN)
    RN_RSS = L*sp.rand(dim,nRN)
    RN_TDOA = L*sp.rand(dim,nRN)

    d_TOA = dist(RN_TOA,BN,0) # actual distances
    TOF = d_TOA/c # actual TOA
    TOA_std = 0.001/c*np.ones(np.shape(TOF))
    TOA = TOF + TOA_std

    RSS_std = 0.001 * np.ones(nRN)
    RSS_np = 2.645 * np.ones(nRN)
    PL0 = 34.7*np.ones(nRN)
    d0 = 1.
    d_RSS = dist(RN_RSS,BN,0) # actual distances
    X = RSS_std * np.random.randn(np.shape(PL0)[0])
    RSS = PL0-10*RSS_np*np.log10(d_RSS/d0)+X
    

    RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
    d = dist(RN_TDOA,BN,0)
    dr = dist(RNr_TDOA,BN,0)
    TDOF = (d-dr)/c # actual TDOA
    TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
    TDOA = TDOF + TDOA_std

    nodes={}
    nodes['BN']= BN
    nodes['RN_RSS']= RN_RSS
    nodes['RN_TOA']= RN_TOA
    nodes['RN_TDOA']= RN_TDOA
    nodes['RNr_TDOA']= RNr_TDOA

    ldp={}
    ldp['RSS'] = RSS
    ldp['RSS_std'] = RSS_std
    ldp['RSS_np'] = RSS_np
    ldp['d0'] = d0
    ldp['PL0'] = PL0
    ldp['TOA'] = TOA
    ldp['TOA_std'] = TOA_std
    ldp['TDOA'] = TDOA
    ldp['TDOA_std'] = TDOA_std
    S = algloc(nodes, ldp)

    S.show(1,1,1)
    
    print 'BN'
    print BN
    print '--------------------------------'

    # Test TOA
    print 'Test TOA'
    print '--------------------------------'
    print 'LS solution'
    print S.ls_locate(0,1,0,'mode')
    print 'WLS solution'
    print S.wls_locate(0,1,0,'mode')
    print 'ML solution'
    print S.ml_locate(BN0,0,1,0)
    print 'CRB'
    print S.crb(BN,0,1,0)

    # Test RSS
    print 'Test RSS'
    print '--------------------------------'
    print 'LS solution'
    print S.ls_locate(1,0,0,'mode')
    print 'WLS solution'
    print S.wls_locate(1,0,0,'mode')
    print 'ML solution'
    print S.ml_locate(BN0,1,0,0)
    print 'CRB'
    print S.crb(BN,1,0,0)

    # Test TDOA
    print 'Test TDOA'
    print '--------------------------------'
    print 'LS solution'
    print S.ls_locate(0,0,1,'mode')
    print 'WLS solution'
    print S.wls_locate(0,0,1,'mode')
    print 'ML solution'
    print S.ml_locate(BN0,0,0,1)
    print 'CRB'
    print S.crb(BN,0,0,1)

    # Test RSS/TOA
    print 'Test RSS/TOA'
    print '--------------------------------'
    print 'LS solution'
    print S.ls_locate(1,1,0,'mode')
    print 'WLS solution'
    print S.wls_locate(1,1,0,'mode')
    print 'ML solution'
    print S.ml_locate(BN0,1,1,0)
    print 'CRB'
    print S.crb(BN,1,1,0)

    # Test RSS/TDOA
    print 'Test RSS/TDOA'
    print '--------------------------------'
    print 'LS solution'
    print S.ls_locate(1,0,1,'mode')
    print 'WLS solution'
    print S.wls_locate(1,0,1,'mode')
    print 'ML solution'
    print S.ml_locate(BN0,1,0,1)
    print 'CRB'
    print S.crb(BN,1,0,1)

    # Test TOA/TDOA
    print 'Test TOA/TDOA'
    print '--------------------------------'
    print 'LS solution'
    print S.ls_locate(0,1,1,'mode')
    print 'WLS solution'
    print S.wls_locate(0,1,1,'mode')
    print 'ML solution'
    print S.ml_locate(BN0,0,1,1)
    print 'CRB'
    print S.crb(BN,0,1,1)

    # Test RSS/TOA/TDOA
    print 'Test RSS/TOA/TDOA'
    print '--------------------------------'
    print 'LS solution'
    print S.ls_locate(1,1,1,'mode')
    print 'WLS solution'
    print S.wls_locate(1,1,1,'mode')
    print 'ML solution'
    print S.ml_locate(BN0,1,1,1)
    print 'CRB'
    print S.crb(BN,1,1,1)
