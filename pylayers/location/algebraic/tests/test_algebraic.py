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


    nodes, ldp, BN0 = scenario()

    S = algloc(nodes,ldp)

    S.show(1,1,1)

    print 'BN'
    print nodes['BN']
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
    print S.crb(nodes['BN'],0,1,0)

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
    print S.crb(nodes['BN'],1,0,0)

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
    print S.crb(nodes['BN'],0,0,1)

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
    print S.crb(nodes['BN'],1,1,0)

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
    print S.crb(nodes['BN'],1,0,1)

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
    print S.crb(nodes['BN'],0,1,1)

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
    print S.crb(nodes['BN'],1,1,1)
