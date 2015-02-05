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
import matplotlib.pyplot as plt

#################################
def cdf(x,colsym="",lab="",lw=1):
    """
    Plot the cumulative density function

    """

    x  = sort(x)
    n  = len(x)
    x2 = repeat(x, 2)
    y2 = hstack([0.0, repeat(arange(1,n) / float(n), 2), 1.0])
    plt.plot(x2,y2,colsym,label=lab,linewidth=lw)

def histo(x, n, fc,norm,xlab,ylab):
    plt.hist(x, n, facecolor=fc,normed=norm)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
