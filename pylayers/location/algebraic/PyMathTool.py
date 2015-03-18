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
import os
from numpy import *
from scipy import *
from numpy.linalg import inv

#################################
def dist(a=array([]),b=array([])):
    """
    Compute euclidean distance between 2 points given by 2 arrays.

    """

    n1=len(a)
    n2=len(b)
    d=0.0
    if (n1==n2):
        d2     = (a - b)*(a - b)
        d    = sum(d2,axis=0)
        return (sqrt(d)[0])
    else :
        print("ERROR: Coordinates are not in the same base")

##################################
def lsop(H=array([])):

    """
    Compute the least square operator 
    
    """
    return (dot(inv(dot(transpose(H),H)),transpose(H)))
