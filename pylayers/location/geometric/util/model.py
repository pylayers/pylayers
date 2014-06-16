# -*- coding:Utf-8 -*-
#
# This file is part of RGPA.

# Foobar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
# authors :
# Nicolas AMIOT		: nicolas.amiot@univ-rennes1.fr
# Bernard UGUEN		: bernard.uguen@univ-rennes1.fr
# Mohamed LAARAIEDH	: mohamed.laaraiedh@univ-rennes1.fr
#
import numpy as np
import scipy as sp


class Model(object):
    """ RSS model

    Attributes
    ----------

    f : frequency in GHz
    n : path loss exponent
    method : used model

    """

    def __init__(self, f=3.0, n=2.64, d0=1.0, method='OneSlope', d=np.array([]), PL=np.array([])):
        self.f = f
        self.n = n
        self.d0 = d0
        self.PL0_c()
        self.method = method

    def info(self):
        print 'frequency (f in GHz) : ', self.f
        print 'path loss exponent (n) : ', self.n
        print 'PL0 : ', self.PL0

    def PL0_c(self, Gt=0, Gr=0):
        """
        PL0_c : PL0 compute
        f  : frequency GHz
        Gt : transmitting antenna gain dB (default 0 dB)
        Gr : receiving antenna gain dB (default 0 dB)
        """
        Gt = 10 ** (Gt / 10.)
        Gr = 10 ** (Gr / 10.)
        ld = 0.3 / self.f
        self.PL0 = -20 * np.log10(ld / (4.0 * np.pi * self.d0))

    def OneSlope(self, d):
        """
        OneSlope model : give Power Level from distance  with OneSlope method

        f : frequency  GHz
        n : path loss exponent
        D : Distance array

        """
        self.PL0_c()
        PL = self.PL0 + 10 * self.n * np.log10(d)
        return(PL)

    def iOneSlope(self, PL):
        """
        inverse OneSlope model : give distance from Power Level with OneSlope method

        f : frequency  GHz
        n : path loss exponent
        D : Distance array

        """
        self.PL0_c()
        d = 10 ** ((PL - self.PL0) / (10 * self.n))
        return(d)

    def get_PL(self, d):
        """
        Get Power Level from a given distance

        d : distance
        """

        met = self.method
        if met == 'OneSlope':
            PL = self.OneSlope(d)

        # ajouter ici les differentes methodes
#		elif mod =='AUTRE'
#			self.AuTRE()

        return(PL)

    def get_d(self, PL):
        """
        Get  distance from a given Power Level

        d : distance
        """
        met = self.method
        if met == 'OneSlope':
            d = self.iOneSlope(PL)

        # ajouter ici les differentes methodes
#		elif mod =='AUTRE'
#			self.AuTRE()

        return(d)
