# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of Network and RGPA.

#PyLayers is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#PyLayers is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with PyLayers.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Nicolas AMIOT        : nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN        : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH    : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import numpy as np
import scipy as sp
from pylayers.util.project import *
import pdb
import doctest


class PLSmodel(PyLayers):
    """ Path Loss Shadowing model

    Attributes
    ----------

    f : float
        frequency in GHz
    rssnp : float
        path loss exponent
    d0    : float
        PL0 distance
    sigrss : float
        shadowing variance

    method : used model



    """
    def __init__(self,f=3.0,rssnp=2.64,d0=1.0,sigrss=3.0,method='mode'):
        self.f  = f
        self.d0  = d0
        self.rssnp  = rssnp
        self.sigrss = sigrss
        self.getPL0()
        self.method = method
        self.param = dict(f=self.f,
                      d0=self.d0,
                      rssnp=self.rssnp,
                      sigrss=self.sigrss,
                      PL0=self.PL0,
                      method=self.method)


    def __repr__(self):
        st =  'frequency (GHz) : ' +str(self.f)+'\n'
        st = st +  'path loss exponent (n) : '+str(self.rssnp)+'\n'
        st = st +  'PL0 (dB): '+str(self.PL0)
        return(st)

    def getPL0(self,Gt=0,Gr=0):
        """ get Path Loss at reference distance d0

        Parameters
        ----------

        PL0_c : PL0 compute
        Gt : transmitting antenna gain dB (default 0 dB)
        Gr : receiving antenna gain dB (default 0 dB)

        Examples
        --------

        >>> from pylayers.network.model import *
        >>> plm = PLSmodel()
        >>> plm.getPL0()

        """

        Gt  = 10**(Gt/10.)
        Gr  = 10**(Gr/10.)
        ld  = 0.3/self.f

        self.PL0 = -20*np.log10(ld/(4.0*np.pi*self.d0))


    def OneSlope(self,r):
        """ OneSlope model : give Power Level from distance  with OneSlope method

        Parameters
        ----------

        r : range (meters)

        Returns
        -------

        PL : float
            path loss values

        """
        try:
            PL = self.PL0+10*self.np*np.log10(r)
        except:
            self.getPL0()
            PL = self.PL0+10*self.np*np.log10(r)
        return(PL)

    def iOneSlope(self,PL):
        """ goes from PL to estimated distance

        inverse OneSlope model : give distance from Power Level with OneSlope method

        Parameters
        ----------

        PL :
            path loss in dB

        Returns
        -------

        r : range array

        """

        try :
            r = 10**((PL-self.PL0)/(10*self.rssnp))
        except:
            self.getPL0()
            r = 10**((PL-self.PL0)/(10*self.rssnp))

        return(r)

    def getPLmean(self, d):
        """ compute PL mean

        Notes
        -----

            $$\bar{PL}=PL_0 - 10 n_p \log_{10}{\frac{d}{d_0}}$$

        """

        PLmean = self.PL0-10*self.rssnp*np.log10(d/self.d0)

        return  PLmean


    def getPL(self,r,RSSStd):
        """ Get Power Level from a given distance

        Parameters
        ----------

        r : range
        RSSStd : range standard deviation

        Examples
        --------

        >>> M = PLSmodel(f=0.3,rssnp=2.64,d0=1,sigrss=3,method='mode')
        >>> PL =  M.getPL(16,1)

        """

        if self.method =='OneSlope':
            PL=self.OneSlope(r)

        elif self.method == 'mode' or self.method == 'median' or self.method == 'mean':
            PLmean          = self.getPLmean(r)
            try:
                shPLmean        = np.shape(PLmean)
                Xrand           = RSSStd*sp.randn(shPLmean[0])
            except:
                Xrand           = RSSStd*sp.randn()
            PL        = PLmean+Xrand

        else :
            raise NameError('Pathloss method name')

        return(PL)



    def getRange(self,RSS,RSSStd):
        """ Get  distance from a given Power Level

        Parameters
        ----------

        RSS :
        R

        Returns
        -------

        r :
        """
        if self.method =='OneSlope':
            r    = self.iOneSlope(RSS)

        elif self.method == 'mode':
            S    = -(np.log(10)/10)* RSSStd/self.rssnp                    # STD of ranges distribution
            M    = (np.log(10)/10)*(self.PL0-RSS)/self.rssnp + np.log(self.d0)        # Mean of ranges distribution
            r    =  np.exp(M-S**2)

        elif self.method == 'median':
            S    = -(np.log(10)/10)* RSSStd/self.rssnp                    # STD of ranges distribution
            M    = (np.log(10)/10)*(self.PL0-RSS)/self.rssnp + np.log(self.d0)        # Mean of ranges distribution
            r    = np.exp(M)

        elif self.method == 'mean':
            S    = -(np.log(10)/10)* RSSStd/self.rssnp                    # STD of ranges distribution
            M    = (np.log(10)/10)*(self.PL0-RSS)/self.rssnp + np.log(self.d0)        # Mean of ranges distribution
            r    =  np.exp(M+0.5*S**2)

        else :
            raise NameError('invalid Pathloss method name for range computation')

        return(r)



    def getRangeStd(self, RSS, RSSStd):
        """Compute Ranges std associated to "Rest" estimator

        Parameters
        ----------

        RSS :
        RSSStd :


        """

        if self.method == 'mode':
            S       = -(np.log(10)/10)* RSSStd/self.rssnp                                    # STD of ranges distribution
            M       = (np.log(10)/10)*(self.PL0-RSS)/self.rssnp + np.log(self.d0)             # Mean of ranges distribution
            r    =  np.sqrt((np.exp(2*M-2*S**2))*(-np.exp(-S**2)+1))

        elif self.method == 'median':
            S       = -(np.log(10)/10)* RSSStd/self.rssnp                                    # STD of ranges distribution
            M       = (np.log(10)/10)*(self.PL0-RSS)/self.rssnp + np.log(self.d0)             # Mean of ranges distribution
            r    = np.sqrt((np.exp(2*M+S**2))*(np.exp(S**2)-1))

        elif self.method == 'mean':
            S       = -(np.log(10)/10)* RSSStd/self.rssnp                                    # STD of ranges distribution
            M       = (np.log(10)/10)*(self.PL0-RSS)/self.rssnp + np.log(self.d0)             # Mean of ranges distribution
            r    =  np.sqrt((np.exp(2*M+3*S**2))*(np.exp(S**2)-1))

        else :
            raise NameError('invalid Pathloss method name STD range computation')

        return (r)


if __name__=='__main__':
    doctest.testmod()

