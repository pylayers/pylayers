# -*- coding:Utf-8 -*-
r"""

.. currentmodule:: pylayers.location.algebraic.algebraic

.. autosummary::
    :toctree: generated

"""
#####################################################################
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
#Mohamed LAARAIEDH
#Bernard Uguen
#####################################################################
import numpy as np
import doctest
import scipy as sp
from scipy import optimize
import numpy.linalg as nplg
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from pylayers.util.geomutil import dist
from pylayers.antprop.loss import PL0
import string
import pdb


class algloc(object):
    """
<<<<<<< HEAD
    This class gathers algebraic localization algorithms

    Attributes
    ----------

    nodes : dictionnary
    ldp : dictionnary
=======
    This class regroups algebraic localization algorithms
        All nodes positions are gatherred in nodes dictionnary
        All location dependent parameters are gathered in ldp dictionnary
>>>>>>> 95d88122037bbe1e5d26beff4fd79e7448eb3dd8

    Methods
    --------
    ls_locate : perform least square position evaluation
    wls_locate : perform weighted least square position evaluation
    ml_locate : perform maximum likelihood position evaluation

    plot : plot scenario
    show : plot scenario
    get_range
    """

    def __init__(self, **kwargs):
        """

        Parameters
        ----------

        dim ( dimension 2 or 3 )
        N : number of nodes

        nodes : dict
            the keys of nodes dict :


        nodes['BN'] : blind node position (dim x N)
        an_rss= RSS nodes position (dim x N)
        an_toa= TOA nodes position (dim x N)
        nodes['RN_TDOA']= TDOA nodes position (dim x N)
        nodes['RNr_TDOA']= TDOA reference nodes position (dim x N)
            you need to repeat along N the reference node ( to be broacasted)

        ldp : dict
            dictionnary containing the location depend parameters
            (between parenthesis is the shape of the value)

        ldp['RSS'] = received power value in dB ( N )
        ldp['RSS_std'] = received power standard ddeviation dB ( N )
        ldp['RSS_np'] = path loss exponent ( N )
        ldp['d0'] = reference distance (integer)
        ldp['PL0'] = path loss @ d0 ( N )
        ldp['TOA'] = Time of arrival in nano seconds (N)
        ldp['TOA_std'] = standard ddeviation in ns ( N )
        ldp['TDOA'] = Time difference of arrival in ns ( N )
        ldp['TDOA_std'] = Time difference of arrival standard dev. in ns ( N )


    Examples
    --------

    .. plot::
        :include-source:


        >>> from pylayers.location.algebraic.algebraic import *
        >>> import scipy as sp
        >>> from pylayers.util.geomutil import dist

        >>> nodes,ldp = scenario()
        >>> L = algloc()
        >>> L.nodes=nodes
        >>> L.ldp=ldp

        >>> print "least square locate RSS + TDOA"
        >>> L.ls_locate(rss=True,toa= False,tdoa = True)
        """


        defaults = { 'an_toa': [],
                     'an_tdoa': [],
                     'an_rss': [],
                     'toa': [],
                     'tdoa': [],
                     'rss': [],
                     'toa_std': [],
                     'tdoa_std': [],
                     'rss_std': [],
                     'rss_np': [],
                     'PL0':[],
                     'd0': 1.,
                     }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]


        self.an_toa = kwargs['an_toa']
        self.an_tdoa = kwargs['an_tdoa']
        self.an_rss = kwargs['an_rss']

        self.toa = kwargs['toa']
        self.tdoa = kwargs['tdoa']
        self.rss = kwargs['rss']

        self.toa_std = kwargs['toa_std']
        self.tdoa_std = kwargs['tdoa_std']
        self.rss_std = kwargs['rss_std']

        self.rss_np = kwargs['rss_np']
        self.d0 = kwargs['d0']

        self.c = 0.2997924583

        # available ldp
        self._av_ldp = []


        # TOA check
        toa_check = np.array([(self.an_toa != defaults['an_toa']), (self.toa != defaults['toa']), (self.toa_std != defaults['toa_std'])])
        if np.alltrue(toa_check) or np.alltrue(~toa_check):
            self._av_ldp.append('toa')
        else:
            raise AttributeError('missing info for TOA')




        


        if self.toa != []:
            self._av_ldp.append('toa')
        if self.tdoa != []:
            self._av_ldp.append('tdoa')
        if self.rss != []:
            self._av_ldp.append('rss')

        # used ldp initilalised at available ldp
        self.used_ldp = self._av_ldp


    def __repr__(self):
        s = 'Available ldp:\n'
        s = s + str(self._av_ldp)[1:-1]
        s = s + '\nUsed ldp for Loc.:\n'
        if 'tdoa' in self.used_ldp:
            s = s + 'tdoa:\n'
            s = s + str(self.an_tdoa)
        if 'toa' in self.used_ldp:
            s = s + '\ntoa:\n'
            s = s + str(self.an_toa)
        if 'rss' in self.used_ldp:
            s = s + '\nrss:\n'
            s = s + str(self.an_rss)

        return(s)


    def plot(self):
        """ plot scenario

        Parameters
        ----------

        rss : boolean
            False
        toa : boolean
            True
        tdoa : boolean
            False

        """

        if self.ldp == []:
            raise ValueError("inputs missed")
        else:
            fig = plt.figure()
            ax = Axes3D(fig)
            # bn = self.bn
            # try:
            #     ax.plot(bn[0, :], bn[1, :], bn[2, :], 'r*', zdir='z', label='blind node')
            # except:
            #     plt.plot(bn[0, :], bn[1, :], 'r*', label='blind node')
            if 'rss' in self.used_ldp:
                rn = self.an_rss
                try:
                    ax.plot(rn[0, :], rn[1, :], rn[2, :], 'ro', zdir='z', label='RSS node')
                except:
                    plt.plot(rn[0, :], rn[1, :], 'ro', label='RSS node')
            if 'toa' in self.used_ldp:
                rn = self.an_toa
                try:
                    ax.plot(rn[0, :], rn[1, :], rn[2, :], 'gs', zdir='z', label='TOA node')
                except:
                    plt.plot(rn[0, :], rn[1, :], 'gs', label='TOA node')

            if 'tdoa' in self.used_ldp:
                rn= self.an_tdoa
                rnr = self.nodes['RNr_TDOA']
                try:
                    ax.plot(rn[0, :], rn[1, :], rn[2, :], 'bD', zdir='z', label='TDOA node')
                    ax.plot(rnr[0, :], rnr[1, :], rnr[2, :], 'kD', zdir='z', label='Ref TDOA node')
                except:
                    plt.plot(rn[0, :], rn[1, :], 'bD', label='TDOA node')
                    plt.plot(rnr[0, :], rnr[1, :], 'kD', label='Ref TDOA node')

        return(fig,ax)

    def show(self):
        """ show scenario

        Parameters
        ----------

        rss : boolean
            False
        toa : boolean
            True
        tdoa : boolean
            False

        Examples
        --------

        .. plot::
            :include-source:

            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> S.show(1,1,1)
            >>> plt.show()
        """

        self.plot()
        plt.legend(numpoints=1)
        plt.show()

    def get_range(self, Rest='mode'):
        """ Compute the range from RSS using the "Rest" estimator

        Parameters
        ----------

         Rest : string

        Returns
        -------

        rg : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> r_mode = S.get_range('mode')
            >>> r_median = S.get_range('median')
            >>> r_mean = S.get_range('mean')

        """

        rss_db = self.rss
        rss_std = self.rss_std
        rss_np = self.rss_np
        d0 = self.d0
        pl0 = self.PL0
        s = (np.log(10) / 10) * rss_std / rss_np
        m = (np.log(10) / 10) * (pl0 - rss_db) / rss_np + np.log(d0)
        if string.lower(Rest) == 'mode':
            rg = np.exp(m - s**2)
        elif string.lower(Rest) == 'median':
            rg = np.exp(m)
        elif string.lower(Rest) == 'mean':
            rg = np.exp(m + 0.5*s**2)
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return rg

    def get_range_std(self, Rest='mode'):
        """
        Compute the RSS range standard deviation using the "Rest" \
        estimator

        Parameters
        ----------
            Rest : string

        Returns
        -------
            rg_std : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> rs_mode = S.get_range_std('mode')
            >>> rs_median = S.get_range_std('median')
            >>> rs_mean = S.get_range_std('mean')
        """

        rss_db = self.rss
        rss_std = self.rss_std
        rss_np = self.rss_np
        d0 = self.d0
        pl0 = self.PL0
        s = (np.log(10) / 10) * rss_std / rss_np
        m = (np.log(10) / 10) * (pl0 - rss_db) / rss_np + np.log(d0)

        if string.lower(Rest) == 'mode':
            rg_std = np.sqrt((np.exp(2*m - 2*s**2))*(1 - np.exp(-s**2)))
        elif string.lower(Rest) == 'median':
            rg_std = np.sqrt(
                (np.exp(2*m + s**2)) * (np.exp(s**2) - 1))
        elif string.lower(Rest) == 'mean':
            rg_std = np.sqrt((np.exp(2*m + 3*s**2))*(np.exp(s**2) - 1))
        else:
            raise ValueError(Rest + ": no such ranging estimator")
        return rg_std

    def ls_locate(self, Rest='mode'):
        """
        This method applies least squares (LS) approximation to get
        position P.

        Parameters
        ----------


        Rest : string

        Returns
        -------

        P : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.ls_locate(1, 0, 0, 'mode')
            >>> P_toa = S.ls_locate(0, 1, 0, 'mode')
            >>> P_tdoa = S.ls_locate(0, 0, 1, 'mode')
            >>> P_rsstoa = S.ls_locate(1, 1, 0, 'mode')
            >>> P_rsstdoa = S.ls_locate(1, 0, 1, 'mode')
            >>> P_toatdoa = S.ls_locate(0, 1, 1, 'mode')
            >>> P_rsstoatdoa = S.ls_locate(1, 1, 1, 'mode')

        """
        if len(self.used_ldp) == 0:
            raise ValueError("inputs missed")
        else:
            # only TDOA
            if 'tdoa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                sh = np.shape(rn_tdoa)

                if sh[1] >= sh[0]:
                    # Construct the vector K (see theory)
                    # Overcoming Singularities In TDOA Based Location
                    # Estimation Using Total Least Square
                    #
                    # (xn-xref)**2+(yn-yref)**2 - r(n,ref)**2
                    #
                    k1 = (np.sum(rn_tdoa*rn_tdoa,axis=0) - np.sum(rnr_tdoa * rnr_tdoa, axis=0))
                    drg = self.c * tdoa_ns
                    k2 = drg*drg
                    K = k1 - k2
                    # Construct the matrix A (see theory)
                    A = np.hstack((rn_tdoa.T-rnr_tdoa.T,drg.reshape(np.shape(tdoa_ns)[0],1)))
                    # Apply LS operator
                    #Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                    Pr = 0.5 * np.dot(nplg.pinv(A),K)
                    P = Pr[:sh[0]].reshape(sh[0], 1)
                else:
                    raise ValueError("Data are not sufficient to perform localization")

            # only TOA
            elif 'toa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_toa = self.an_toa
                toa_ns = self.toa
                sh = np.shape(rn_toa)
                if sh[1] > sh[0]:
                    # Construct the vector K (see theory)
                    rn2 = (np.sum(rn_toa * rn_toa, axis=0))
                    k1 = rn2[1:] - rn2[0:1]
                    rg = self.c * toa_ns
                    rg2 = rg * rg
                    k2 = rg2[0:1] - rg2[1:]
                    K = k1 + k2
                    # Construct the matrix A (see theory)
                    A = rn_toa[:, 1:].T - rn_toa[:, 0]
                    # Apply LS operator
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_toa[:, 0:1]))
                else:
                    raise ValueError("Data are not sufficient to perform localization")

            # only Rss
            elif 'rss' in self.used_ldp and len(self.used_ldp) == 1:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                sh = np.shape(rn_rss)
                if sh[1] > sh[0]:
                    # Construct the vector K (see theory)
                    rn2 = np.sum(rn_rss * rn_rss, axis=0)
                    k1 = rn2[1:] - rn2[0:1]
                    rg = self.get_range(Rest)
                    rg2 = rg * rg
                    k2 = rg2[0:1] - rg2[1:]
                    K = k1 + k2
                    # Construct the matrix A (see theory)
                    A = rn_rss[:, 1:].T - rn_rss[:, 0]
                    # Apply LS operator
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))

                else:
                    raise ValueError("Data are not sufficient to perform localization")

            # TOA + RSS
            elif 'toa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_toa = self.an_toa
                toa_ns = self.toa
                sh1 = np.shape(rn_rss)
                sh2 = np.shape(rn_toa)
                if sh1[1] > 1 and sh1[1]+sh2[1] > sh1[0]:
                    # Construct the vector K_rss (see theory)
                    rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                    k1_rss = rn2_rss[1:] - rn2_rss[0:1]
                    rg_rss = self.get_range(Rest)
                    rg2_rss = rg_rss * rg_rss
                    k2_rss = rg2_rss[0:1] - rg2_rss[1:]
                    K_rss = k1_rss + k2_rss
                    # Construct the matrix A_rss (see theory)
                    A_rss = rn_rss[:, 1:].T - rn_rss[:, 0]
                    # Construct the vector K_toa (see theory)
                    rn2_toa = (np.sum(rn_toa * rn_toa, axis=0))
                    k1_toa = rn2_toa[:] - rn2_rss[0:1]
                    rg_toa = self.c * toa_ns
                    rg2_toa = rg_toa * rg_toa
                    k2_toa = rg2_rss[0:1] - rg2_toa[:]
                    K_toa = k1_toa + k2_toa
                    # Construct the matrix A_toa (see theory)
                    A_toa = rn_toa[:, :].T - rn_rss[:, 0]
                    # Apply LS operator
                    sh3 = np.shape(K_rss)[0]
                    sh4 = np.shape(K_toa)[0]
                    K = np.vstack((K_rss.reshape(sh3, 1),K_toa.reshape(sh4, 1)))
                    A = np.vstack((A_rss, A_toa))
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                elif sh2[1] > 1 and sh1[1]+sh2[1] > sh1[0]:
                    # Construct the vector K_toa (see theory)
                    rn2_toa = (np.sum(rn_toa * rn_toa, axis=0))
                    k1_toa = rn2_toa[1:] - rn2_toa[0:1]
                    rg_toa = self.c * toa_ns
                    rg2_toa = rg_toa * rg_toa
                    k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                    K_toa = k1_toa + k2_toa
                    # Construct the matrix A_toa (see theory)
                    A_toa = rn_toa[:, 1:].T - rn_toa[:, 0]
                    # Construct the vector K_rss (see theory)
                    rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                    k1_rss = rn2_rss[:] - rn2_toa[0:1]
                    rg_rss = self.get_range(Rest)
                    rg2_rss = rg_rss * rg_rss
                    k2_rss = rg2_toa[0:1] - rg2_rss[:]
                    K_rss = k1_rss + k2_rss
                    # Construct the matrix A_rss (see theory)
                    A_rss = rn_rss[:, :].T - rn_toa[:, 0]
                    # Apply LS operator
                    sh3 = np.shape(K_rss)[0]
                    sh4 = np.shape(K_toa)[0]
                    K = np.vstack((K_rss.reshape(sh3, 1),K_toa.reshape(sh4, 1)))
                    A = np.vstack((A_rss, A_toa))
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                else:
                    raise ValueError("Data are not sufficient to perform localization")
                    
            # TDOA + TOA
            elif 'tdoa' in self.used_ldp and \
                 'toa' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                rn_toa = self.an_toa
                toa_ns = self.toa
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                sh1 = np.shape(rn_toa)
                sh2 = np.shape(rn_tdoa)
                # Construct the vector K_Ttoa (see theory)
                rn2 = (np.sum(rn_toa * rn_toa, axis=0))
                k1_toa = rn2[1:] - rn2[0:1]
                rg_toa = self.c * toa_ns
                rg2_toa = rg_toa * rg_toa
                k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                K_toa = k1_toa + k2_toa
                # Construct the matrix A_toa (see theory)
                A_toa = np.hstack((rn_toa[:, 1:].T - rn_toa[:, 0],
                                   np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) * (rn_tdoa - rnr_tdoa),axis=0))
                drg = self.c * tdoa_ns
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack((rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_toa)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_toa.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_toa, A_tdoa))
                Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                P = Pr[:sh2[0]].reshape(sh2[0], 1)

            # TDOA + RSS
            elif 'tdoa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                sh1 = np.shape(rn_rss)
                sh2 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2 = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2[1:] - rn2[0:1]
                rg = self.get_range(Rest)
                rg2 = rg * rg
                k2_rss = rg2[0:1] - rg2[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack((rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) * (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                k2_tdoa = drg * drg
                K_tdoa= k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack((rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_rss)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_rss.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_rss, A_tdoa))
                Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                P = Pr[:sh2[0]].reshape(sh2[0], 1)

            # TDOA + TOA + RSS
            else:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_toa = self.an_toa
                toa_ns = self.toa
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                sh1 = np.shape(rn_toa)
                sh2 = np.shape(rn_rss)
                sh3 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2_rss[1:] - rn2_rss[0:1]
                rg_rss = self.get_range(Rest)
                rg2_rss = rg_rss * rg_rss
                k2_rss = rg2_rss[0:1] - rg2_rss[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack((rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh2[1] - 1, 1))))
                # Construct the vector K_toa (see theory)
                rn2_toa = (np.sum(rn_toa * rn_toa, axis=0))
                k1_toa = rn2_toa[1:] - rn2_toa[0:1]
                rg_toa = self.c * toa_ns
                rg2_toa = rg_toa * rg_toa
                k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                K_toa = k1_toa + k2_toa
                # Construct the matrix A_toa (see theory)
                A_toa = np.hstack((rn_toa[:, 1:].T - rn_toa[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) * (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg2 = drg * drg
                k2_tdoa = drg2
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack((rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh3[1], 1)))
                # Apply LS operator
                sh4 = np.shape(K_rss)[0]
                sh5 = np.shape(K_toa)[0]
                sh6 = np.shape(K_tdoa)[0]
                K = np.vstack((np.vstack((K_rss.reshape(sh4, 1), K_toa.reshape(sh5, 1))),K_tdoa.reshape(sh6, 1)))
                A = np.vstack((np.vstack((A_rss, A_toa)), A_tdoa))
                Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, A)), np.dot(A.T, K))
                P = Pr[:sh3[0]].reshape(sh3[0], 1)

            return P

    def wls_locate(self, Rest='mode'):
        """
        This method applies weighted least squares (WLS) approximation
        to get position P.

        Parameters
        ----------
            rss : boolean
            toa : boolean
            tdoa : boolean
            Rest : string

        Returns
        -------
            P : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.wls_locate(1, 0, 0, 'mode')
            >>> P_toa = S.wls_locate(0, 1, 0, 'mode')
            >>> P_tdoa = S.wls_locate(0, 0, 1, 'mode')
            >>> P_rsstoa = S.wls_locate(1, 1, 0, 'mode')
            >>> P_rsstdoa = S.wls_locate(1, 0, 1, 'mode')
            >>> P_toatdoa = S.wls_locate(0, 1, 1, 'mode')
            >>> P_rsstoatdoa = S.wls_locate(1, 1, 1, 'mode')

        """
        if len(self.used_ldp) == 0:
            raise ValueError("inputs missed")
        else:
            # only TDOA
            if 'tdoa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.ldp['TDOA_std']
                sh = np.shape(rn_tdoa)
                if sh[1] >= sh[0]:
                    # Construct the vector K (see theory)
                    k1 = (np.sum((rn_tdoa - rnr_tdoa) * (rn_tdoa - rnr_tdoa), axis=0))
                    drg = self.c * tdoa_ns
                    drg_std = self.c * tdoa_std
                    k2 = drg * drg
                    K = k1 - k2
                    # Construct the matrix A (see theory)
                    A = np.hstack((rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh[1], 1)))
                    # Construct the Covariance Matrix
                    C = np.diag(drg_std[:] ** 2)
                    # Apply LS operator
                    Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                    P = Pr[:sh[0]].reshape(sh[0], 1)
                else:
                    raise ValueError("Data are not sufficient to perform localization")
                

            # only TOA
            elif 'toa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_toa = self.an_toa
                toa_ns = self.toa
                toa_std = self.toa_std
                sh = np.shape(rn_toa)
                if sh[1] > sh[0]:
                    # Construct the vector K (see theory)
                    rn2 = (np.sum(rn_toa * rn_toa, axis=0))
                    k1 = rn2[1:] - rn2[0:1]
                    rg = self.c * toa_ns
                    rg_std = self.c * toa_std
                    rg2 = rg * rg
                    k2 = rg2[0:1] - rg2[1:]
                    K = k1 + k2
                    # Construct the matrix A (see theory)
                    A = rn_toa[:, 1:].T - rn_toa[:, 0]
                    # Construct the Covariance Matrix
                    C = np.diag(rg_std[1:] ** 2)
                    # Apply LS operator
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                    P = P.reshape(np.shape(rn_toa[:, 0:1]))
                else:
                    raise ValueError("Data are not sufficient to perform localization")
                

            # only Rss
            elif 'rss' in self.used_ldp and len(self.used_ldp) == 1:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                sh = np.shape(rn_rss)
                if sh[1] > sh[0]:
                    # Construct the vector K (see theory)
                    rn2 = np.sum(rn_rss * rn_rss, axis=0)
                    k1 = rn2[1:] - rn2[0:1]
                    rg = self.get_range(Rest)
                    rg_std = self.get_range_std(Rest)
                    rg2 = rg * rg
                    k2 = rg2[0:1] - rg2[1:]
                    K = k1 + k2
                    # Construct the matrix A (see theory)
                    A = rn_rss[:, 1:].T - rn_rss[:, 0]
                    # Construct the Covariance Matrix
                    C = np.diag((rg_std[1:]) ** 2)
                    # Apply LS operator
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                else:
                    raise ValueError("Data are not sufficient to perform localization")


            # TOA + RSS
            elif 'toa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_toa = self.an_toa
                toa_ns = self.toa
                toa_std = self.toa_std
                sh1 = np.shape(rn_rss)
                sh2 = np.shape(rn_toa)
                if sh1[1] > 1 and sh1[1]+sh2[1] > sh1[0]:
                    # Construct the vector K_rss (see theory)
                    rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                    k1_rss = rn2_rss[1:] - rn2_rss[0:1]
                    rg_rss = self.get_range(Rest)
                    rg_rss_std = self.get_range_std(Rest)
                    rg2_rss = rg_rss * rg_rss
                    k2_rss = rg2_rss[0:1] - rg2_rss[1:]
                    K_rss = k1_rss + k2_rss
                    # Construct the matrix A_rss (see theory)
                    A_rss = rn_rss[:, 1:].T - rn_rss[:, 0]
                    # Construct the vector K_toa (see theory)
                    rn2_toa = (np.sum(rn_toa * rn_toa, axis=0))
                    k1_toa = rn2_toa[:] - rn2_rss[0:1]
                    rg_toa = self.c * toa_ns
                    rg_toa_std = self.c * toa_std
                    rg2_toa = rg_toa * rg_toa
                    k2_toa = rg2_rss[0:1] - rg2_toa[:]
                    K_toa = k1_toa + k2_toa
                    # Construct the matrix A_toa (see theory)
                    A_toa = rn_toa[:, :].T - rn_rss[:, 0]
                    # Apply LS operator
                    sh3 = np.shape(K_toa)[0]
                    sh4 = np.shape(K_rss)[0]
                    K = np.vstack((K_rss.reshape(sh4, 1),
                                   K_toa.reshape(sh3, 1)))
                    A = np.vstack((A_rss, A_toa))
                    C = np.diag(np.hstack((rg_rss_std[1:] ** 2,
                                           rg_toa_std[:] ** 2)))
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                elif sh2[1] > 1 and sh1[1]+sh2[1] > sh1[0]:
                    # Construct the vector K_toa (see theory)
                    rn2_toa = (np.sum(rn_toa * rn_toa, axis=0))
                    k1_toa = rn2_toa[1:] - rn2_toa[0:1]
                    rg_toa = self.c * toa_ns
                    rg_toa_std = self.c * toa_std
                    rg2_toa = rg_toa * rg_toa
                    k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                    K_toa = k1_toa + k2_toa
                    # Construct the matrix A_toa (see theory)
                    A_toa = rn_toa[:, 1:].T - rn_toa[:, 0]
                    # Construct the vector K_rss (see theory)
                    rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                    k1_rss = rn2_rss[:] - rn2_toa[0:1]
                    rg_rss = self.get_range(Rest)
                    rg_rss_std = self.get_range_std(Rest)
                    rg2_rss = rg_rss * rg_rss
                    k2_rss = rg2_toa[0:1] - rg2_rss[:]
                    K_rss = k1_rss + k2_rss
                    # Construct the matrix A_rss (see theory)
                    A_rss = rn_rss[:, :].T - rn_toa[:, 0]
                    # Apply LS operator
                    sh3 = np.shape(K_toa)[0]
                    sh4 = np.shape(K_rss)[0]
                    K = np.vstack((K_rss.reshape(sh4, 1),
                                   K_toa.reshape(sh3, 1)))
                    A = np.vstack((A_rss, A_toa))
                    C = np.diag(np.hstack((rg_rss_std[:] ** 2,
                                           rg_toa_std[1:] ** 2)))
                    P = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                else:
                    raise ValueError("Data are not sufficient to perform localization")


            # TDOA + TOA
            elif 'tdoa' in self.used_ldp and \
                 'toa' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                rn_toa = self.an_toa
                toa_ns = self.toa
                toa_std = self.toa_std
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.ldp['TDOA_std']
                sh1 = np.shape(rn_toa)
                sh2 = np.shape(rn_tdoa)
                # Construct the vector K_toa (see theory)
                rntoa2 = (np.sum(rn_toa * rn_toa, axis=0))
                k1_toa = rntoa2[1:] - rntoa2[0:1]
                rg_toa = self.c * toa_ns
                rg_toa_std = self.c * toa_std
                rg2_toa = rg_toa * rg_toa
                k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                K_toa = k1_toa + k2_toa
                # Construct the matrix A_toa (see theory)
                A_toa = np.hstack((rn_toa[:, 1:].T - rn_toa[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) * (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack((rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_toa)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_toa.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_toa, A_tdoa))
                C = np.diag(np.hstack((rg_toa_std[1:] ** 2, drg_std[:] ** 2)))
                Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:sh2[0]].reshape(sh2[0], 1)

            # TDOA + RSS
            elif 'tdoa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.ldp['TDOA_std']
                sh1 = np.shape(rn_rss)
                sh2 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2 = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2[1:] - rn2[0:1]
                rg = self.get_range(Rest)
                rg_std = self.get_range_std(Rest)
                rg2 = rg * rg
                k2_rss = rg2[0:1] - rg2[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack((rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) * (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack((rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_rss)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_rss.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_rss, A_tdoa))
                C = np.diag(np.hstack((rg_std[1:] ** 2, drg_std[:] ** 2)))
                Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:sh2[0]].reshape(sh2[0], 1)

            else:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_toa = self.an_toa
                toa_ns = self.toa
                toa_std = self.toa_std
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.ldp['TDOA_std']
                sh1 = np.shape(rn_toa)
                sh2 = np.shape(rn_rss)
                sh3 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2 = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2[1:] - rn2[0:1]
                rg_rss = self.get_range(Rest)
                rg_rss_std = self.get_range_std(Rest)
                rg2_rss = rg_rss * rg_rss
                k2_rss = rg2_rss[0:1] - rg2_rss[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack((rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_toa (see theory)
                rntoa2 = (np.sum(rn_toa * rn_toa, axis=0))
                k1_toa = rntoa2[1:] - rntoa2[0:1]
                rg_toa = self.c * toa_ns
                rg_toa_std = self.c * toa_std
                rg2_toa = rg_toa * rg_toa
                k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                K_toa = k1_toa + k2_toa
                # Construct the matrix A_toa (see theory)
                A_toa = np.hstack((rn_toa[:, 1:].T - rn_toa[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) * (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A (see theory)
                A_tdoa = np.hstack((rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh4 = np.shape(K_toa)[0]
                sh5 = np.shape(K_rss)[0]
                sh6 = np.shape(K_tdoa)[0]
                K = np.vstack((np.vstack((K_rss.reshape(sh5, 1), K_toa.reshape(sh4, 1))), K_tdoa.reshape(sh6, 1)))
                A = np.vstack((np.vstack((A_rss, A_toa)), A_tdoa))
                C = np.diag(np.hstack((np.hstack((rg_rss_std[1:] ** 2, rg_toa_std[1:] ** 2)), drg_std[:] ** 2)))
                Pr = 0.5 * np.dot(nplg.inv(np.dot(A.T, np.dot(nplg.inv(C), A))), np.dot(np.dot(A.T, nplg.inv(C)), K))
                P = Pr[:sh3[0]].reshape(sh3[0], 1)

            return P

    def ml_function(self, P):
        """
        This defines the ML function to be minimized if ML estimator
        is used

        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean

        Returns
        -------
            ML : numpy.ndarray
        """
        if len(self.used_ldp) == 0:
            raise ValueError("inputs missed")
        else:
            # only TDOA
            if 'tdoa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.ldp['TDOA_std']
                sh1 = np.shape(rn_tdoa)[1]
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                # construct the ML function to be minimized
                rnmp = rn_tdoa - np.outer(P, np.ones(sh1))
                mrnmp = np.sqrt(np.diag(np.dot(rnmp.T, rnmp)))
                rnrmp = rnr_tdoa - np.outer(P, np.ones(sh1))
                mrnrmp = np.sqrt(np.diag(np.dot(rnrmp.T, rnrmp)))
                rdrg = mrnmp - mrnrmp
                dd = (drg - rdrg) ** 2 / (2 * drg_std ** 2)
                ML = dd.sum(axis=0)

            # only TOA
            elif 'toa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_toa = self.an_toa
                toa_ns = self.toa
                toa_std = self.toa_std
                sh1 = np.shape(rn_toa)[1]
                rg = self.c * toa_ns
                rg_std = self.c * toa_std
                # construct the ML function to be minimized
                rnmp = rn_toa - np.outer(P, np.ones(sh1))
                mrnmp = np.sqrt(np.diag(np.dot(rnmp.T, rnmp)))
                dd = (rg - mrnmp) ** 2
                uu = dd / (2 * rg_std ** 2) + np.log(np.sqrt(2 * np.pi) * rg_std)
                ML = uu.sum(axis=0)

            # only Rss
            elif 'rss' in self.used_ldp and len(self.used_ldp) == 1:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                sh1 = np.shape(rn_rss)[1]
                S = (np.log(10) / 10) * rss_std / rss_np
                M = (np.log(10) / 10) * (pl0 - rss_db) / rss_np + np.log(d0)
                # construct the ML function to be minimized
                rnmp = rn_rss - np.outer(P, np.ones(sh1))
                mrnmp = np.sqrt(np.diag(np.dot(rnmp.T, rnmp)))
                dd = (M - S ** 2 - np.log(mrnmp)) ** 2
                uu = dd / (2 * S ** 2)
                ML = uu.sum(axis=0)

            # TDOA + TOA
            elif 'tdoa' in self.used_ldp and \
                 'toa' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                ML = self.ml_function(P, 0, 1, 0) +\
                     self.ml_function(P, 0, 0, 1)

            # TDOA + RSS
            elif 'tdoa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                ML = self.ml_function(P, 1, 0, 0) +\
                     self.ml_function(P, 0, 0, 1)


            # TOA + RSS
            elif 'toa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                ML = self.ml_function(P, 1, 0, 0) +\
                     self.ml_function(P, 0, 1, 0)
            # TDOA + TOA + RSS
            else:
                ML = self.ml_function(P, 1, 0, 0) +\
                     self.ml_function(P, 0, 1, 0) +\
                     self.ml_function(P, 0, 0, 1)

            return ML

    def ml_locate(self, P0):
        """
        This applies ML estimator to compute position P

        Parameters
        ----------
            P0 : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean

        Returns
        -------
            P : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> P_rss = S.ml_locate(BN0, 1, 0, 0)
            >>> P_toa = S.ml_locate(BN0, 0, 1, 0)
            >>> P_tdoa = S.ml_locate(BN0, 0, 0, 1)
            >>> P_rsstoa = S.ml_locate(BN0, 1, 1, 0)
            >>> P_rsstdoa = S.ml_locate(BN0, 1, 0, 1)
            >>> P_toatdoa = S.ml_locate(BN0, 0, 1, 1)
            >>> P_rsstoatdoa = S.ml_locate(BN0, 1, 1, 1)
        """
        if len(self.used_ldp) == 0:
            raise ValueError("inputs missed")
        else:
            P = optimize.fmin(self.ml_function, P0, disp=0)
            return P.reshape(np.shape(P0))

    def fim(self, P):
        """
        Compute the fisher information matrix in P for the given
        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean

        Returns
        -------
            FIM : numpy.ndarray

        """
        if len(self.used_ldp) == 0:
            raise ValueError("inputs missed")

        else:
            # only TDOA
            if 'tdoa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_tdoa = self.an_tdoa
                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.ldp['TDOA_std']
                sh1 = np.shape(rn_tdoa)
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                FIM = np.zeros((sh1[0],sh1[0]))
                for i in range(sh1[1]):
                    f1 = P - rn_tdoa[:, i:i + 1]
                    f2 = P - rnr_tdoa[:, i:i + 1]
                    pmrn = np.sqrt(np.dot(f1.T,f1))
                    pmrnr = np.sqrt(np.dot(f2.T,f2))
                    FIM += (1 / drg_std[i] ** 2) * np.dot(f1 / pmrn - f2 / pmrnr, (f1 / pmrn - f2 / pmrnr).T)

            # only TOA
            elif 'toa' in self.used_ldp and len(self.used_ldp) == 1:
                rn_toa = self.an_toa
                toa_ns = self.toa
                toa_std = self.toa_std
                sh1 = np.shape(rn_toa)
                rg = self.c * toa_ns
                rg_std = self.c * toa_std
                FIM = np.zeros((sh1[0], sh1[0]))
                for i in range(sh1[1]):
                    f1 = P - rn_toa[:, i:i + 1]
                    FIM += np.dot(f1 , f1.T) / ((rg_std[i] ** 2) *np.dot(f1.T , f1))
                    FIM += np.dot(f1 , f1.T) / ((rg_ramerstd[i] ** 2) *np.dot(f1.T , f1))

            # only Rss
            elif 'rss' in self.used_ldp and len(self.used_ldp) == 1:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                sh1 = np.shape(rn_rss)
                d0 = self.d0
                pl0 = self.PL0
                S = (np.log(10) / 10) * rss_std / rss_np
                FIM = np.zeros((sh1[0],sh1[0]))
                for i in range(sh1[1]):
                    f1 = P - rn_rss[:, i:i + 1]
                    FIM += np.dot(f1, f1.T) / ((S[0] ** 2) * (np.dot( f1.T, f1)) ** 2)
            # TDOA + TOA
            elif 'tdoa' in self.used_ldp and \
                 'toa' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                FIM = self.fim(P, 0, 1, 0) + self.fim(P, 0, 0, 1)


            # TDOA + RSS
            elif 'tdoa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                FIM = self.fim(P, 1, 0, 0) + self.fim(P, 0, 0, 1)

            # TOA + RSS
            elif 'toa' in self.used_ldp and \
                 'rss' in self.used_ldp and \
                  len(self.used_ldp) == 2:
                FIM = self.fim(P, 1, 0, 0) + self.fim(P, 0, 1, 0)

            # TDOA + TOA + RSS
            else:
                FIM = self.fim(P, 1, 0, 0) + self.fim(P, 0, 1, 0) +\
                    self.fim(P, 0, 0, 1)

            return FIM

    def crb(self, P):
        """
        This method compute the cramer rao bound (CRB) at position P.

        Parameters
        ----------
            P : numpy.ndarray
            rss : boolean
            toa : boolean
            tdoa : boolean

        Returns
        -------
            CRB : float

        Examples
        --------
        .. plot::
            :include-source:

            >>> nodes, ldp, BN0 = scenario()
            >>> S = algloc(nodes, ldp)
            >>> crb_rss = S.crb(BN, 1, 0, 0)
            >>> crb_toa = S.crb(BN, 0, 1, 0)
            >>> crb_tdoa = S.crb(BN, 0, 0, 1)
            >>> crb_rsstoa = S.crb(BN, 1, 1, 0)
            >>> crb_rsstdoa = S.crb(BN, 1, 0, 1)
            >>> crb_toatdoa = S.crb(BN, 0, 1, 1)
            >>> crb_rsstoatdoa = S.crb(BN, 1, 1, 1)
        """

        FIM = self.fim(P)
        CRB = np.sqrt(np.trace(nplg.inv(FIM)))
        return CRB


def scenario():
    """
    This method is not a class member, it defines a sample scenario.
    Defines a sample scenario for testing
    
    Returns
    -------
        nodes : dictionnary 
            'BN' 
            'RN_RSS'
            'RN_TDOA'
            'RN_TOA'
            'RNr_TDOA'
        ldp   : all ldps of bn0
            'PL0'
            'RSS'
            'RSS_np'
            'RSS_std'
            'TDOA'
            'TDOA_std'
            'TOA'
            'TOA_std'
        bn0   : The true position of blind node


    Examples
    --------

    .. plot::
        :include-source:


        >>> from pylayers.location.algebraic import *
        >>> import scipy as sp
        >>> from pylayers.util.geomutil import dist
        >>> nRN = 4
        >>> dim = 3 # 2 for 2D, 3 for 3D
        >>> L = 20.
        >>> c = 0.3
        >>> BN = L*sp.rand(dim,1)
        >>> BN0 = L*sp.rand(dim,1)
        >>> RN_TOA = L*sp.rand(dim,nRN)
        >>> RN_RSS = L*sp.rand(dim,nRN)
        >>> RN_TDOA = L*sp.rand(dim,nRN)

        >>> d_TOA = dist(RN_TOA,BN,0) # actual distances
        >>> TOF = d_TOA/c # actual TOA
        >>> TOA_std = 0.001/c*np.ones(np.shape(TOF))
        >>> TOA = TOF + TOA_std

        >>> rss_std = 0.001 * np.ones(nRN)
        >>> rss_np = 2.645 * np.ones(nRN)
        >>> pl0 = 34.7*np.ones(nRN)
        >>> d0 = 1.
        >>> d_RSS = dist(RN_RSS,BN,0) # actual distances
        >>> X = rss_std * np.random.randn(np.shape(pl0)[0])
        >>> rss_db = pl0-10*rss_np*np.log10(d_RSS/d0)+X

        >>> RNr_TDOA = np.zeros((dim,nRN))#L*sp.rand(dim,nRN)
        >>> d = dist(RN_TDOA,BN,0)
        >>> dr = dist(RNr_TDOA,BN,0)
        >>> TDOF = (d-dr)/c # actual TDOA
        >>> TDOA_std = 0.001/c*np.ones(np.shape(TDOF))
        >>> TDOA = TDOF + TDOA_std

        >>> nodes={}
        >>> nodes['BN']= BN
        >>> an_rss= RN_RSS
        >>> an_toa= RN_TOA
        >>> nodes['RN_TDOA']= RN_TDOA
        >>> nodes['RNr_TDOA']= RNr_TDOA

        >>> ldp={}
        >>> ldp['RSS'] = rss_db
        >>> ldp['RSS_std'] = rss_std
        >>> ldp['RSS_np'] = rss_np
        >>> ldp['d0'] = d0
        >>> ldp['PL0'] = pl0
        >>> ldp['TOA'] = TOA
        >>> ldp['TOA_std'] = TOA_std
        >>> ldp['TDOA'] = TDOA
        >>> ldp['TDOA_std'] = TDOA_std

        >>> print 'Nodes'
        >>> print nodes
        >>> print 'LDPs'
        >>> print ldp
        >>> print 'BN0:initial guess for ML estimator'
        >>> print BN0
        """
    nrn = 4
    dim = 3  # 2 for 2D, 3 for 3D
    L = 20.
    c = 0.3
    bn = L * sp.rand(dim, 1)
    bn0 = L * sp.rand(dim, 1)
    rn_toa = L * sp.rand(dim, nrn)
    rn_rss = L * sp.rand(dim, nrn)
    rn_tdoa = L * sp.rand(dim, nrn)

    d_toa = dist(rn_toa, bn, 0)  # actual distances
    tof = d_toa / c  # actual TOA
    toa_std = 0.001 / c * np.ones(np.shape(tof))
    toa_ns = tof + toa_std

    rss_std = 0.001 * np.ones(nrn)
    rss_np = 2.645 * np.ones(nrn)
    pl0 = 34.7 * np.ones(nrn)
    d0 = 1.
    d_rss = dist(rn_rss, bn, 0)  # actual distances
    x = rss_std * np.random.randn(np.shape(pl0)[0])
    rss_db = pl0 - 10 * rss_np * np.log10(d_rss / d0) + x

    rnr_tdoa = np.zeros((dim, nrn))  # L*sp.rand(dim,nRN)
    d = dist(rn_tdoa, bn, 0)
    dr = dist(rnr_tdoa, bn, 0)
    tdof = (d - dr) / c  # actual TDOA
    tdoa_std = 0.001 / c * np.ones(np.shape(tdof))
    tdoa_ns = tdof + tdoa_std
    
    nodes = {}
    nodes['BN'] = bn
    an_rss = rn_rss
    an_toa = rn_toa
    nodes['RN_TDOA'] = rn_tdoa
    nodes['RNr_TDOA'] = rnr_tdoa

    ldp = {}
    ldp['RSS'] = rss_db
    ldp['RSS_std'] = rss_std
    ldp['RSS_np'] = rss_np
    ldp['d0'] = d0
    ldp['PL0'] = pl0
    ldp['TOA'] = toa_ns
    ldp['TOA_std'] = toa_std
    ldp['TDOA'] = tdoa_ns
    ldp['TDOA_std'] = tdoa_std

    return nodes, ldp, bn0


if __name__ == "__main__":
    doctest.testmod()
