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
#Bernard UGUEN        : buguen@univ-rennes1.fr
#Mohamed LAARAIEDH    : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
from numpy import *
from scipy import *
from scipy import optimize
from numpy.linalg import *
import cvxmod as cvxm
import cvxopt as cvxo
from string import *
from HDFLocation import *

import os
import sys


sys.path.append('../Util') 
from Values_compute import *
from Model import *    


if __name__=="__main__":

    nbnod = 4

    L=20
    bn   = L*sp.rand(2)
    H1=5
    H2=5
    H3=5

    Ap  = np.zeros(2)
    MS  = np.zeros(2)
    BS  = np.zeros(2)

    ############# RSS #############
    Ap1 = np.array([0,0])
    Ap2 = np.array([L,0])
    Ap3 = np.array([0,L])
    Ap4 = np.array([L,L])

    Ap = np.vstack((Ap,Ap1))
    Ap = np.vstack((Ap,Ap2))
    Ap = np.vstack((Ap,Ap3))
    Ap = np.vstack((Ap,Ap4))
    Ap = np.delete(Ap,0,0)
    std_RSS = 5
    std_RSS_v=std_RSS*np.ones(4)

    VC=Values_compute()
    VC.nodes = Ap#[Ap1,Ap2,Ap3,Ap4]
    VC.bn    = bn
    VC.noise = True
    VC.ctype = 'RSS' 
    VC.std = std_RSS
    noise_vect,v_RSS = VC.compute()
    #####################

    ########### TOA ##########
    MS1 = np.array([L/2,L/4])
    MS2 = np.array([3*L/4,L/2])
    MS3 = np.array([L/4,L/2])
    MS4 = np.array([L/2,3*L/4])

    MS = np.vstack((MS,MS1))
    MS = np.vstack((MS,MS2))
    MS = np.vstack((MS,MS3))
    MS = np.vstack((MS,MS4))
    MS = np.delete(MS,0,0)



    std_TOA = 5
    std_TOA_v=std_TOA*np.ones(4)

    VC=Values_compute()
    VC.nodes = [MS1,MS2,MS3,MS4]
    VC.bn    = bn
    VC.noise = True
    VC.ctype = 'TOA' 
    VC.std = std_TOA
    noise_vect,v_TOA = VC.compute()
    #####################



    ######## TDOA #############

    BS1 = np.array([L/2,-L])
    BS2 = np.array([2*L,L/2])
    BS3 = np.array([L/2,2*L])
    BS4 = np.array([-L,L/2])

    MS = np.vstack((MS,MS1))
    MS = np.vstack((MS,MS2))
    MS = np.vstack((MS,MS3))
    MS = np.vstack((MS,MS4))
     MS = np.delete(MS,0,0)



    std_TDOA = 5
    std_TDOA_v=std_TDOA*np.ones(3)

    VC=Values_compute()
    VC.nodes = Ap#[Ap1,Ap2,Ap3,Ap4]
    VC.bn    = bn
    VC.noise = True
    VC.ctype = 'TDOA' 
    VC.std = std_TDOA
    noise_vect,v_TDOA = VC.compute()

    ###########GLUE PART##########
    v_ToA= np.zeros(1)

    for i in range (len(v_TOA)):
            v_ToA=np.vstack((v_ToA,v_TOA[i]))
    v_ToA = np.delete(v_ToA,0,0)

    v_TDoA= np.zeros(1)
    for i in range (len(v_TDOA)):
            v_TDoA=np.vstack((v_TDoA,v_TDOA[i]))
    v_TDoA = np.delete(v_TDoA,0,0)


    v_Rss= np.zeros(1)
    for i in range (len(v_RSS)):
            v_Rss=np.vstack((v_Rss,v_RSS[i]))
    v_Rss = np.delete(v_Rss,0,0)

    ############

    M=Model()

    hdf=HDFLocation(Ap,MS,BS)

    hdf.TWLSHDFLocate(Ap.T,MS.T,BS.T, v_ToA, std_TOA_v, v_TDoA, std_TDOA_v, M.PL0, 1.0, v_Rss, 2.4*np.ones(4,1) , std_RSS_v, 'mean')
