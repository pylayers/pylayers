# -*- coding:Utf-8 -*-
r"""

.. currentmodule:: pylayers.location.algebraic.algebraic

.. autosummary::
    :toctree: generated

"""

# authors :
# Mohamed LAARAIEDH
# Nicolas AMIOT
# Bernard Uguen

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


class Algloc(object):
    """
    This class gathers algebraic localization algorithms

    Attributes
    ----------


    Methods
    --------
    ls_locate : perform least square position evaluation
    wls_locate : perform weighted least square position evaluation
    ml_locate : perform maximum likelihood position evaluation

    plot : plot scenario
    """

    def __init__(self, **kwargs):
        """

        Parameters
        ----------
        an_toa: np.ndarray(shape=(3,0))
            anchor nodes for toa
        an_tdoa: np.ndarray(shape=(3,0))
            anchor nodes for tdoa
        an_rss: np.ndarray(shape=(3,0))
            anchor nodes for rss
        toa: np.ndarray(shape=(0))
            toa in ns 
        tdoa: np.ndarray(shape=(0))
            tdoa in ns 
        rss: np.ndarray(shape=(0))
            rss value in dB
        toa_std: np.ndarray(shape=(0))
            toa standard deviation in ns 
        tdoa_std: np.ndarray(shape=(0))
            tdoa standard deviation in ns 
        rss_std: np.ndarray(shape=(0))
            rss standard deviation in dB
        rss_np: np.ndarray(shape=(0))
            pathloss exponent
        PL0: np.ndarray(shape=(0))
            loss at d0
        d0: 1.
        Rest: str
            'mode'(default)|'median'|'mean'
            range estimation mode

        used_ldp : list 
            dictionnary to determien which ldp are used for localization etaimtion

        bnGT : np.ndarray(shape=(3,0))
            blind node position

    Examples
    --------

    .. plot::
        :include-source:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from pylayers.location.algebraic.algebraic import Algloc
        >>> from pylayers.location.observables import Observables
        >>> an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
        >>> an = an.T
        >>> bn = np.array([1, 1, 2.])
        >>> O_toa = Observables(an=an, bn=bn, mode='toa')
        >>> A = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                       toa_std=O_toa.noise_model['std'], bnGT=bn)
        >>> print('Blind node :' + str(bn) + '\n')
        >>> print('TOA')
        >>> print('---')
        >>> print('LS')
        >>> print(A.ls_locate())
        >>> print('WLS')
        >>> print(A.wls_locate())
        >>> print('ML')
        >>> print(A.ml_locate())
        >>> A.show()
        >>> plt.show()
        """

        defaults = {'an_toa': np.ndarray(shape=(3, 0)),
                    'an_tdoa': np.ndarray(shape=(3, 0)),
                    'an_rss': np.ndarray(shape=(3, 0)),
                    'toa': np.ndarray(shape=(0)),
                    'tdoa': np.ndarray(shape=(0)),
                    'tdoa_ref': 0,
                    'rss': np.ndarray(shape=(0)),
                    'toa_std': np.ndarray(shape=(0)),
                    'tdoa_std': np.ndarray(shape=(0)),
                    'rss_std': np.ndarray(shape=(0)),
                    'rss_np': np.ndarray(shape=(0)),
                    'PL0': np.ndarray(shape=(0)),
                    'd0': 1.,
                    'Rest': 'mode',
                    'bnGT': np.ndarray(shape=(3, 0))
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]
            if not isinstance(kwargs[k], np.ndarray):
                if not isinstance(kwargs[k], str):
                    kwargs[k] = np.array([kwargs[k]])
            setattr(self, k, kwargs[k])
        self.c = 0.2997924583

        # available ldp
        self._av_ldp = []

        if self.bnGT.ndim > 2:
            raise AttributeError('blind node \'bn\' shape must be (3 x 1)')
        if self.bnGT.ndim < 2:
            self.bnGT = self.bnGT.reshape((self.bnGT.shape[0], 1))
        if self.bnGT.shape[0] == 2:
            self.bnGT = np.vstack(
                (self.bnGT, np.zeros(self.bnGT.shape[1])))
        elif self.bnGT.shape[0] != 3:
            raise AttributeError('Blind node first dimension reserved to space\
                                  (x,y,z) coordinates')

        ###########
        # TOA check
        ###########

        if self.an_toa.ndim > 2:
            raise AttributeError('Anchors \'an\' shape must be (3 x Na)')
        if self.an_toa.ndim < 2:
            self.an_toa = self.an_toa.reshape((self.an_toa.shape[0], 1))
        if self.an_toa.shape[0] == 2:
            self.an_toa = np.vstack(
                (self.an_toa, np.zeros(self.an_toa.shape[1])))
        elif self.an_toa.shape[0] != 3:
            raise AttributeError('Anchors TOA first dimension reserved to space\
                                  (x,y,z) coordinates')

        self.Ntoa = self.an_toa.shape[1]

        if len(self.toa.shape) > 1:
            try:
                self.toa = self.toa.reshape(self.Ntoa)
            except:
                raise AttributeError('Wrong shape for toa')
        if len(self.toa_std.shape) > 1:
            try:
                self.toa_std = self.toa_std.reshape(self.Ntoa)
            except:
                raise AttributeError('Wrong shape for toa_std')

        if (self.an_toa.shape[1] != self.toa.shape[0]):
            raise AttributeError('toa shape mishmatch or missing')
        elif (self.an_toa.shape[1] != self.toa_std.shape[0]):
            if len(self.toa_std) == 1:
                self.toa_std = np.array(self.toa_std) * np.ones(self.Ntoa)
            else:
                raise AttributeError('toa_std shape mishmatch or missing')

        if self.Ntoa > 0:
            self._av_ldp.append('toa')

        ###########
        # TDOA check
        ###########

        if self.an_tdoa.ndim > 2:
            raise AttributeError('Anchors \'an\' shape must be (3 x Na)')
        if self.an_tdoa.ndim < 2:
            self.an_tdoa = self.an_tdoa.reshape((self.an_tdoa.shape[0], 1))
        if self.an_tdoa.shape[0] == 2:
            self.an_tdoa = np.vstack(
                (self.an_tdoa, np.zeros(self.an_tdoa.shape[1])))
        elif self.an_tdoa.shape[0] != 3:
            raise AttributeError('Anchors TDOA first dimension reserved to space\
                                  (x,y,z) coordinates')

        self.Ntdoa = self.an_tdoa.shape[1]

        if len(self.tdoa.shape) > 1:
            try:
                self.tdoa = self.tdoa.reshape(self.Ntdoa - 1)
            except:
                raise AttributeError('Wrong shape for tdoa')
        if len(self.tdoa_std.shape) > 1:
            try:
                self.tdoa_std = self.tdoa_std.reshape(self.Ntdoa - 1)
            except:
                raise AttributeError('Wrong shape for tdoa_std')

        if (self.an_tdoa.shape[1] != self.tdoa.shape[0] + 1) and\
                self.an_tdoa.shape[1] != 0:
            raise AttributeError('tdoa shape mishmatch or missing')
        elif (self.an_tdoa.shape[1] != self.tdoa_std.shape[0] + 1) and\
                self.an_tdoa.shape[1] != 0:
            if len(self.tdoa_std) == 1:
                self.tdoa_std = np.array(self.tdoa_std) * np.ones(self.Ntdoa)
            else:
                raise AttributeError('tdoa_std shape mishmatch or missing')

        if self.Ntdoa > 0:
            self._av_ldp.append('tdoa')

        # RSS

        if self.an_rss.ndim > 2:
            raise AttributeError('Anchors \'an\' shape must be (3 x Na)')
        if self.an_rss.ndim < 2:
            self.an_rss = self.an_rss.reshape((self.an_rss.shape[0], 1))
        if self.an_rss.shape[0] == 2:
            self.an_rss = np.vstack(
                (self.an_rss, np.zeros(self.an_rss.shape[1])))
        elif self.an_rss.shape[0] != 3:
            raise AttributeError('Anchors RSS first dimension reserved to space\
                                  (x,y,z) coordinates')

        self.Nrss = self.an_rss.shape[1]

        if self.rss.ndim > 1:
            try:
                self.rss = self.rss.reshape(self.Nrss)
            except:
                raise AttributeError('Wrong shape for rss')

        if self.rss_std.ndim > 1:
            try:
                self.rss_std = self.rss_std.reshape(self.Nrss)
            except:
                raise AttributeError('Wrong shape for rss_std')

        if not (self.an_rss.shape[1] == self.rss.shape[0]):
            raise AttributeError('rss shape mishmatch or missing')
        elif (self.an_rss.shape[1] != self.rss_std.shape[0]):
            if len(self.rss_std) == 1:
                self.rss_std = np.array(self.rss_std) * np.ones(self.Nrss)
            else:
                raise AttributeError('rss_std shape mishmatch or missing')
        elif (self.an_rss.shape[1] != self.rss_np.shape[0]):
            if len(self.rss_np) == 1:
                self.rss_np = np.array(self.rss_np) * np.ones(self.Nrss)
            else:
                raise AttributeError('rss_np shape mishmatch or missing')
        elif (self.an_rss.shape[1] != self.PL0.shape[0]):
            if len(self.PL0) == 1:
                self.PL0 = np.array(self.PL0) * np.ones(self.Nrss)
            else:
                raise AttributeError('PL0 shape mishmatch or missing')

        if np.alltrue(self.rss < 0):
            self.rss = -self.rss

        if self.Nrss > 0:
            self._av_ldp.append('rss')

        # used ldp initilalised at available ldp
        self.used_ldp = {'toa': False, 'tdoa': False, 'rss': False}
        self.used_ldp.update({v: True for v in self._av_ldp})

    def __repr__(self):
        s = 'Available ldp:\n'
        s = s + str(self._av_ldp)
        s = s + '\nUsed ldp for Loc.:\n'
        s = s + str(self.used_ldp) + '\n'
        if self.used_ldp['tdoa']:
            s = s + '\ntdoa:\n'
            s = s + str(self.an_tdoa)
        if self.used_ldp['toa']:
            s = s + '\ntoa:\n'
            s = s + str(self.an_toa)
        if self.used_ldp['rss']:
            s = s + '\nrss:\n'
            s = s + str(self.an_rss)

        return(s)

    def __add__(self, A):

        nd = {}

        nd['an_toa'] = np.hstack((self.an_toa, A.an_toa))
        nd['toa'] = np.hstack((self.toa, A.toa))
        nd['toa_std'] = np.hstack((self.toa_std, A.toa_std))


        if 'tdoa' in self._av_ldp and 'tdoa' in A._av_ldp:
            if not np.alltrue(self.an_tdoa[self.tdoa_ref] == A.an_tdoa[A.tdoa_ref]):
                raise AttributeError(
                    'tdoa reference nodes are note the same => cannot add these 2 algloc Algloc objects')
            ext_an = np.delete(A.an_tdoa,A.tdoa_ref,1)

        else :
            ext_an = A.an_tdoa

        nd['an_tdoa'] = np.hstack((self.an_tdoa, ext_an))
        nd['tdoa'] = np.hstack((self.tdoa, A.tdoa))
        nd['tdoa_std'] = np.hstack((self.tdoa_std, A.tdoa_std))
        nd['tdoa_ref'] = self.tdoa_ref

        nd['an_rss'] = np.hstack((self.an_rss, A.an_rss))
        nd['rss'] = np.hstack((self.rss, A.rss))
        nd['rss_std'] = np.hstack((self.rss_std, A.rss_std))
        nd['rss_np'] = np.hstack((self.rss_np, A.rss_np))
        nd['PL0'] = np.hstack((self.PL0, A.PL0))
        nd['Rest'] = self.Rest
        nd['d0'] = self.d0
        nd['bnGT'] = self.bnGT

        return Algloc(**nd)


    def show(self, **kwargs):
        self.plot(**kwargs)

    def plot(self, **kwargs):
        """ plot scenario

        """

        defaults = {'fig': [],
                    'ax': [],
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        if kwargs['fig'] == []:
            fig = plt.figure()
        else:
            fig = kwargs['fig']
        if kwargs['ax'] == []:
            ax = Axes3D(fig)
        else:
            ax = kwargs['ax']

        ax.plot(self.bnGT[0, :], self.bnGT[1, :], self.bnGT[
                2, :], 'r*', zdir='z', label='blind node')

        if self.used_ldp['rss']:
            ax.plot(self.an_rss[0, :], self.an_rss[1, :], self.an_rss[2, :],
                    'ks', zdir='z', label='RSS node')

        if self.used_ldp['toa']:
            ax.plot(self.an_toa[0, :], self.an_toa[1, :], self.an_toa[2, :],
                    'gs', zdir='z', label='TOA node')

        if self.used_ldp['tdoa']:

            ax.plot(self.an_tdoa[0, :], self.an_tdoa[1, :], self.an_tdoa[2, :],
                    'bs', zdir='z', label='TDOA node')
            ax.plot(self.an_tdoa[0, self.tdoa_ref], self.an_tdoa[1, self.tdoa_ref], self.an_tdoa[2, self.tdoa_ref],
                    'ks', zdir='z', label='Ref TDOA node')

        if hasattr(self, 'pe_ls'):
            X = np.concatenate([self.pe_ls, self.bnGT], axis=1)

            ax.plot(self.pe_ls[0, :], self.pe_ls[1, :], self.pe_ls[2, :],
                    'gD', zdir='z', label='Pe LS')
            ax.plot(X[0, :], X[1, :], X[2, :], linewidth=0.5, color='k')

        if hasattr(self, 'pe_wls'):
            X = np.concatenate([self.pe_wls, self.bnGT], axis=1)
            ax.plot(self.pe_wls[0, :], self.pe_wls[1, :], self.pe_wls[2, :],
                    'cD', zdir='z', label='Pe WLS')
            ax.plot(X[0, :], X[1, :], X[2, :], linewidth=0.5, color='k')

        if hasattr(self, 'pe_ml'):
            X = np.concatenate([self.pe_ml, self.bnGT], axis=1)
            ax.plot(self.pe_ml[0, :], self.pe_ml[1, :], self.pe_ml[2, :],
                    'mD', zdir='z', label='Pe ML')
            ax.plot(X[0, :], X[1, :], X[2, :], linewidth=0.5, color='k')

        ax.legend()

        return(fig, ax)

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

    def get_range(self):
        """ Compute the range from RSS using the "self.Rest" estimator

        Returns
        -------

        rg : numpy.ndarray



        """

        rss_db = self.rss
        rss_std = self.rss_std
        rss_np = self.rss_np
        d0 = self.d0
        pl0 = self.PL0
        s = (np.log(10) / 10) * rss_std / rss_np
        m = (np.log(10) / 10) * (pl0 - rss_db) / rss_np + np.log(d0)
        if string.lower(self.Rest) == 'mode':
            rg = np.exp(m - s**2)
        elif string.lower(self.Rest) == 'median':
            rg = np.exp(m)
        elif string.lower(self.Rest) == 'mean':
            rg = np.exp(m + 0.5 * s**2)
        else:
            raise ValueError(self.Rest + ": no such ranging estimator")
        return rg

    def get_range_std(self):
        """
        Compute the RSS range standard deviation using the "self.Rest" \
        estimator


        """

        rss_db = self.rss
        rss_std = self.rss_std
        rss_np = self.rss_np
        d0 = self.d0
        pl0 = self.PL0
        s = (np.log(10) / 10) * rss_std / rss_np
        m = (np.log(10) / 10) * (pl0 - rss_db) / rss_np + np.log(d0)

        if string.lower(self.Rest) == 'mode':
            rg_std = np.sqrt((np.exp(2 * m - 2 * s**2)) * (1 - np.exp(-s**2)))
        elif string.lower(self.Rest) == 'median':
            rg_std = np.sqrt(
                (np.exp(2 * m + s**2)) * (np.exp(s**2) - 1))
        elif string.lower(self.Rest) == 'mean':
            rg_std = np.sqrt((np.exp(2 * m + 3 * s**2)) * (np.exp(s**2) - 1))
        else:
            raise ValueError(self.Rest + ": no such ranging estimator")
        return rg_std

    def ls_locate(self):
        """
        This method applies least squares (LS) approximation to get
        position P.


        Returns
        -------

        pe_ls : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from pylayers.location.algebraic.algebraic import Algloc
        >>> from pylayers.location.observables import Observables
        >>> an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
        >>> an = an.T
        >>> bn = np.array([1, 1, 2.])
        >>> O_toa = Observables(an=an, bn=bn, mode='toa')
        >>> A = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                       toa_std=O_toa.noise_model['std'], bnGT=bn)
        >>> A.ls_locate()

        """
        if np.sum(self.used_ldp.values()) == 0:
            raise ValueError("inputs missed")
        else:
            # only TDOA
            if self.used_ldp['tdoa'] and\
               not self.used_ldp['toa'] and\
               not self.used_ldp['rss']:
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
                # rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                sh = np.shape(rn_tdoa)
                if sh[1] >= sh[0]:
                    # Construct the vector K (see theory)
                    # Overcoming Singularities In TDOA Based Location
                    # Estimation Using Total Least Square
                    #
                    # (xn-xref)**2+(yn-yref)**2 - r(n,ref)**2
                    #
                    k1 = (np.sum(rn_tdoa * rn_tdoa, axis=0) -
                          np.sum(rnr_tdoa * rnr_tdoa, axis=0))
                    drg = self.c * tdoa_ns
                    k2 = drg * drg
                    K = k1 - k2
                    # Construct the matrix A (see theory)
                    A = np.hstack(
                        (rn_tdoa.T - rnr_tdoa.T, drg.reshape(np.shape(tdoa_ns)[0], 1)))
                    # Apply LS operator
                    #Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
                    Pr = 0.5 * np.dot(nplg.pinv(A), K)
                    P = Pr[:sh[0]].reshape(sh[0], 1)
                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # only TOA
            elif self.used_ldp['toa'] and\
                    not self.used_ldp['tdoa'] and\
                    not self.used_ldp['rss']:
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
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_toa[:, 0:1]))
                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # only Rss
            elif self.used_ldp['rss'] and\
                    not self.used_ldp['tdoa'] and\
                    not self.used_ldp['toa']:
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
                    rg = self.get_range()
                    rg2 = rg * rg
                    k2 = rg2[0:1] - rg2[1:]
                    K = k1 + k2
                    # Construct the matrix A (see theory)
                    A = rn_rss[:, 1:].T - rn_rss[:, 0]
                    # Apply LS operator
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))

                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # TOA + RSS
            elif self.used_ldp['toa'] and \
                    self.used_ldp['rss'] and \
                    not self.used_ldp['tdoa']:
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
                if sh1[1] > 1 and sh1[1] + sh2[1] > sh1[0]:
                    # Construct the vector K_rss (see theory)
                    rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                    k1_rss = rn2_rss[1:] - rn2_rss[0:1]
                    rg_rss = self.get_range()
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
                    K = np.vstack(
                        (K_rss.reshape(sh3, 1), K_toa.reshape(sh4, 1)))
                    A = np.vstack((A_rss, A_toa))
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                elif sh2[1] > 1 and sh1[1] + sh2[1] > sh1[0]:
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
                    rg_rss = self.get_range()
                    rg2_rss = rg_rss * rg_rss
                    k2_rss = rg2_toa[0:1] - rg2_rss[:]
                    K_rss = k1_rss + k2_rss
                    # Construct the matrix A_rss (see theory)
                    A_rss = rn_rss[:, :].T - rn_toa[:, 0]
                    # Apply LS operator
                    sh3 = np.shape(K_rss)[0]
                    sh4 = np.shape(K_toa)[0]
                    K = np.vstack(
                        (K_rss.reshape(sh3, 1), K_toa.reshape(sh4, 1)))
                    A = np.vstack((A_rss, A_toa))
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # TDOA + TOA
            elif self.used_ldp['tdoa'] and \
                    self.used_ldp['toa'] and \
                    not self.used_ldp['rss']:
                rn_toa = self.an_toa
                toa_ns = self.toa
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
                #rnr_tdoa = self.nodes['RNr_TDOA']
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
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) *
                                  (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack(
                    (rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_toa)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_toa.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_toa, A_tdoa))
                Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
                P = Pr[:sh2[0]].reshape(sh2[0], 1)

            # TDOA + RSS
            elif self.used_ldp['tdoa'] and \
                    self.used_ldp['rss'] and \
                    not self.used_ldp['toa']:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
#                rnr_tdoa = self.nodes['RNr_TDOA']

                tdoa_ns = self.tdoa
                sh1 = np.shape(rn_rss)
                sh2 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2 = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2[1:] - rn2[0:1]
                rg = self.get_range()
                rg2 = rg * rg
                k2_rss = rg2[0:1] - rg2[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack(
                    (rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) *
                                  (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack(
                    (rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_rss)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_rss.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_rss, A_tdoa))
                Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
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
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
#                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                sh1 = np.shape(rn_toa)
                sh2 = np.shape(rn_rss)
                sh3 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2_rss[1:] - rn2_rss[0:1]
                rg_rss = self.get_range()
                rg2_rss = rg_rss * rg_rss
                k2_rss = rg2_rss[0:1] - rg2_rss[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack(
                    (rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh2[1] - 1, 1))))
                # Construct the vector K_toa (see theory)
                rn2_toa = (np.sum(rn_toa * rn_toa, axis=0))
                k1_toa = rn2_toa[1:] - rn2_toa[0:1]
                rg_toa = self.c * toa_ns
                rg2_toa = rg_toa * rg_toa
                k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                K_toa = k1_toa + k2_toa
                # Construct the matrix A_toa (see theory)
                A_toa = np.hstack(
                    (rn_toa[:, 1:].T - rn_toa[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) *
                                  (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg2 = drg * drg
                k2_tdoa = drg2
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack(
                    (rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh3[1], 1)))
                # Apply LS operator
                sh4 = np.shape(K_rss)[0]
                sh5 = np.shape(K_toa)[0]
                sh6 = np.shape(K_tdoa)[0]
                K = np.vstack((np.vstack(
                    (K_rss.reshape(sh4, 1), K_toa.reshape(sh5, 1))), K_tdoa.reshape(sh6, 1)))
                A = np.vstack((np.vstack((A_rss, A_toa)), A_tdoa))
                Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, A)), np.dot(A.T, K))
                P = Pr[:sh3[0]].reshape(sh3[0], 1)

            self.pe_ls = P
            return P

    def wls_locate(self):
        """
        This method applies weighted least squares (WLS) approximation
        to get position P.

        Returns
        -------
            pe_wls : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from pylayers.location.algebraic.algebraic import Algloc
        >>> from pylayers.location.observables import Observables
        >>> an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
        >>> an = an.T
        >>> bn = np.array([1, 1, 2.])
        >>> O_toa = Observables(an=an, bn=bn, mode='toa')
        >>> A = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                       toa_std=O_toa.noise_model['std'], bnGT=bn)
        >>> A.wls_locate()

        """
        if np.sum(self.used_ldp.values()) == 0:
            raise ValueError("inputs missed")
        else:
            # only TDOA
            if self.used_ldp['tdoa'] and\
               not self.used_ldp['toa'] and\
               not self.used_ldp['rss']:
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
                #rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.tdoa_std
                sh = np.shape(rn_tdoa)
                if sh[1] >= sh[0]:
                    # Construct the vector K (see theory)
                    k1 = (np.sum((rn_tdoa - rnr_tdoa) *
                                 (rn_tdoa - rnr_tdoa), axis=0))
                    drg = self.c * tdoa_ns
                    drg_std = self.c * tdoa_std
                    k2 = drg * drg
                    K = k1 - k2
                    # Construct the matrix A (see theory)
                    A = np.hstack(
                        (rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh[1], 1)))
                    # Construct the Covariance Matrix
                    C = np.diag(drg_std[:] ** 2)
                    # Apply LS operator
                    Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                      np.dot(np.dot(A.T, nplg.pinv(C)), K))
                    P = Pr[:sh[0]].reshape(sh[0], 1)
                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # only TOA
            elif self.used_ldp['toa'] and\
                    not self.used_ldp['rss'] and\
                    not self.used_ldp['tdoa']:
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
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                     np.dot(np.dot(A.T, nplg.pinv(C)), K))
                    P = P.reshape(np.shape(rn_toa[:, 0:1]))
                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # only Rss
            elif self.used_ldp['rss'] and\
                    not self.used_ldp['toa'] and\
                    not self.used_ldp['tdoa']:
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
                    rg = self.get_range()
                    rg_std = self.get_range_std()
                    rg2 = rg * rg
                    k2 = rg2[0:1] - rg2[1:]
                    K = k1 + k2
                    # Construct the matrix A (see theory)
                    A = rn_rss[:, 1:].T - rn_rss[:, 0]
                    # Construct the Covariance Matrix
                    C = np.diag((rg_std[1:]) ** 2)
                    # Apply LS operator
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                     np.dot(np.dot(A.T, nplg.pinv(C)), K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # TOA + RSS
            elif self.used_ldp['toa'] and \
                    self.used_ldp['rss'] and \
                    not self.used_ldp['tdoa']:
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
                if sh1[1] > 1 and sh1[1] + sh2[1] > sh1[0]:
                    # Construct the vector K_rss (see theory)
                    rn2_rss = np.sum(rn_rss * rn_rss, axis=0)
                    k1_rss = rn2_rss[1:] - rn2_rss[0:1]
                    rg_rss = self.get_range()
                    rg_rss_std = self.get_range_std()
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
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                     np.dot(np.dot(A.T, nplg.pinv(C)), K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                elif sh2[1] > 1 and sh1[1] + sh2[1] > sh1[0]:
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
                    rg_rss = self.get_range()
                    rg_rss_std = self.get_range_std()
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
                    P = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                     np.dot(np.dot(A.T, nplg.pinv(C)), K))
                    P = P.reshape(np.shape(rn_rss[:, 0:1]))
                else:
                    raise ValueError(
                        "Data are not sufficient to perform localization")

            # TDOA + TOA
            elif self.used_ldp['tdoa'] and \
                    self.used_ldp['toa'] and \
                not self.used_ldp['rss']:
                rn_toa = self.an_toa
                toa_ns = self.toa
                toa_std = self.toa_std
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
#                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.tdoa_std
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
                A_toa = np.hstack(
                    (rn_toa[:, 1:].T - rn_toa[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) *
                                  (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack(
                    (rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_toa)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_toa.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_toa, A_tdoa))
                C = np.diag(np.hstack((rg_toa_std[1:] ** 2, drg_std[:] ** 2)))
                Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                  np.dot(np.dot(A.T, nplg.pinv(C)), K))
                P = Pr[:sh2[0]].reshape(sh2[0], 1)

            # TDOA + RSS
            elif self.used_ldp['tdoa'] and\
                    self.used_ldp['rss'] and\
                    not self.used_ldp['toa']:
                rn_rss = self.an_rss
                rss_db = self.rss
                rss_std = self.rss_std
                rss_np = self.rss_np
                d0 = self.d0
                pl0 = self.PL0
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
#                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.tdoa_std
                sh1 = np.shape(rn_rss)
                sh2 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2 = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2[1:] - rn2[0:1]
                rg = self.get_range()
                rg_std = self.get_range_std()
                rg2 = rg * rg
                k2_rss = rg2[0:1] - rg2[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack(
                    (rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) *
                                  (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A_tdoa (see theory)
                A_tdoa = np.hstack(
                    (rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh2[1], 1)))
                # Apply LS operator
                sh3 = np.shape(K_rss)[0]
                sh4 = np.shape(K_tdoa)[0]
                K = np.vstack((K_rss.reshape(sh3, 1), K_tdoa.reshape(sh4, 1)))
                A = np.vstack((A_rss, A_tdoa))
                C = np.diag(np.hstack((rg_std[1:] ** 2, drg_std[:] ** 2)))
                Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                  np.dot(np.dot(A.T, nplg.pinv(C)), K))
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
                rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
                rnr_tdoa = self.an_tdoa[
                    :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
#                rnr_tdoa = self.nodes['RNr_TDOA']
                tdoa_ns = self.tdoa
                tdoa_std = self.tdoa_std
                sh1 = np.shape(rn_toa)
                sh2 = np.shape(rn_rss)
                sh3 = np.shape(rn_tdoa)
                # Construct the vector K_rss (see theory)
                rn2 = np.sum(rn_rss * rn_rss, axis=0)
                k1_rss = rn2[1:] - rn2[0:1]
                rg_rss = self.get_range()
                rg_rss_std = self.get_range_std()
                rg2_rss = rg_rss * rg_rss
                k2_rss = rg2_rss[0:1] - rg2_rss[1:]
                K_rss = k1_rss + k2_rss
                # Construct the matrix A_rss (see theory)
                A_rss = np.hstack(
                    (rn_rss[:, 1:].T - rn_rss[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_toa (see theory)
                rntoa2 = (np.sum(rn_toa * rn_toa, axis=0))
                k1_toa = rntoa2[1:] - rntoa2[0:1]
                rg_toa = self.c * toa_ns
                rg_toa_std = self.c * toa_std
                rg2_toa = rg_toa * rg_toa
                k2_toa = rg2_toa[0:1] - rg2_toa[1:]
                K_toa = k1_toa + k2_toa
                # Construct the matrix A_toa (see theory)
                A_toa = np.hstack(
                    (rn_toa[:, 1:].T - rn_toa[:, 0], np.zeros((sh1[1] - 1, 1))))
                # Construct the vector K_tdoa (see theory)
                k1_tdoa = (np.sum((rn_tdoa - rnr_tdoa) *
                                  (rn_tdoa - rnr_tdoa), axis=0))
                drg = self.c * tdoa_ns
                drg_std = self.c * tdoa_std
                k2_tdoa = drg * drg
                K_tdoa = k1_tdoa - k2_tdoa
                # Construct the matrix A (see theory)
                A_tdoa = np.hstack(
                    (rn_tdoa.T - rnr_tdoa.T, drg.reshape(sh3[1], 1)))
                # Apply LS operator
                sh4 = np.shape(K_toa)[0]
                sh5 = np.shape(K_rss)[0]
                sh6 = np.shape(K_tdoa)[0]
                K = np.vstack((np.vstack(
                    (K_rss.reshape(sh5, 1), K_toa.reshape(sh4, 1))), K_tdoa.reshape(sh6, 1)))
                A = np.vstack((np.vstack((A_rss, A_toa)), A_tdoa))
                C = np.diag(np.hstack(
                    (np.hstack((rg_rss_std[1:] ** 2, rg_toa_std[1:] ** 2)), drg_std[:] ** 2)))
                Pr = 0.5 * np.dot(nplg.pinv(np.dot(A.T, np.dot(nplg.pinv(C), A))),
                                  np.dot(np.dot(A.T, nplg.pinv(C)), K))
                P = Pr[:sh3[0]].reshape(sh3[0], 1)

            self.pe_wls = P

            return P

    def _ml_function_TOA(self, P):
        """
        This defines the ML function to be minimized if ML estimator
        is used for TOA

        Parameters
        ----------
            pe_ml : numpy.ndarray
        """
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
        uu = dd / (2 * rg_std ** 2) + \
            np.log(np.sqrt(2 * np.pi) * rg_std)
        ML = uu.sum(axis=0)

        return ML

    def _ml_function_TDOA(self, P):
        """
        This defines the ML function to be minimized if ML estimator
        is used for TDOA

        Parameters
        ----------
            P : numpy.ndarray
        """
        rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
        rnr_tdoa = self.an_tdoa[
            :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
#                rnr_tdoa = self.nodes['RNr_TDOA']
        tdoa_ns = self.tdoa
        tdoa_std = self.tdoa_std
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

        return ML

    def _ml_function_RSS(self, P):
        """
        This defines the ML function to be minimized if ML estimator
        is used for RSS

        Parameters
        ----------
            P : numpy.ndarray
        """
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

        return ML

    def ml_function(self, P):
        """
        This defines the ML function to be minimized if ML estimator
        is used

        Parameters
        ----------
            P : numpy.ndarray
                gradient initial guess 
        Returns
        -------
            ML : numpy.ndarray
        """
        if np.sum(self.used_ldp.values()) == 0:
            raise ValueError("inputs missed")
        else:
            # only TDOA
            if self.used_ldp['tdoa'] and\
               not self.used_ldp['toa'] and\
               not self.used_ldp['rss']:
                ML = self._ml_function_TDOA(P)

            # only TOA
            elif self.used_ldp['toa'] and\
                    not self.used_ldp['tdoa'] and\
                    not self.used_ldp['rss']:
                ML = self._ml_function_TOA(P)

            # only Rss
            elif self.used_ldp['rss'] and\
                    not self.used_ldp['toa'] and\
                    not self.used_ldp['tdoa']:
                ML = self._ml_function_RSS(P)

            # TDOA + TOA
            elif self.used_ldp['tdoa'] and \
                    self.used_ldp['toa'] and \
                    not self.used_ldp['rss']:
                ML = self._ml_function_TDOA(P) +\
                    self._ml_function_TOA(P)

            # TDOA + RSS
            elif self.used_ldp['tdoa'] and \
                    self.used_ldp['rss'] and \
                    not self.used_ldp['toa']:
                ML = self._ml_function_TDOA(P) +\
                    self._ml_function_RSS(P)

            # TOA + RSS
            elif self.used_ldp['toa'] and \
                    self.used_ldp['rss'] and \
                    not self.used_ldp['tdoa']:
                ML = self._ml_function_TOA(P) +\
                    self._ml_function_RSS(P)

            # TDOA + TOA + RSS
            else:
                ML = self._ml_function_TDOA(P) +\
                    self._ml_function_TOA(P) +\
                    self._ml_function_RSS(P)

            return ML

    def ml_locate(self, P0=np.zeros((3, 1))):
        """
        This applies ML estimator to compute position P

        Parameters
        ----------
            P0 : numpy.ndarray
                initial guess of gradient

        Returns
        -------
            pe_ml : numpy.ndarray

        Examples
        --------

        .. plot::
            :include-source:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from pylayers.location.algebraic.algebraic import Algloc
        >>> from pylayers.location.observables import Observables
        >>> an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
        >>> an = an.T
        >>> bn = np.array([1, 1, 2.])
        >>> O_toa = Observables(an=an, bn=bn, mode='toa')
        >>> A = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                       toa_std=O_toa.noise_model['std'], bnGT=bn)
        >>> A.ml_locate()
        """
        if np.sum(self.used_ldp.values()) == 0:
            raise ValueError("inputs missed")
        else:
            P = optimize.fmin(self.ml_function, P0, disp=0)

            P = P.reshape(np.shape(P0))
            self.pe_ml = P

            return P


    def fim(self, P):
        """
        Compute the fisher information matrix in P for the given

        Parameters
        ----------
            P : numpy.ndarray

        Returns
        -------
            FIM : numpy.ndarray


        Examples
        --------

        .. plot::
            :include-source:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from pylayers.location.algebraic.algebraic import Algloc
        >>> from pylayers.location.observables import Observables
        >>> an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
        >>> an = an.T
        >>> bn = np.array([1, 1, 2.])
        >>> O_toa = Observables(an=an, bn=bn, mode='toa')
        >>> A = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                       toa_std=O_toa.noise_model['std'], bnGT=bn)
        >>> A.fim(A.bnGT)

        """
        if np.sum(self.used_ldp.values()) == 0:
            raise ValueError("inputs missed")

        else:
            # only TDOA
            if self.used_ldp['tdoa'] and np.sum(self.used_ldp.values()) == 1:
                FIM = self._fim_TDOA(P)

            # only TOA
            elif self.used_ldp['toa'] and np.sum(self.used_ldp.values()) == 1:
                FIM = self._fim_TOA(P)

            # only Rss
            elif self.used_ldp['rss'] and np.sum(self.used_ldp.values()) == 1:
                FIM = self._fim_RSS(P)

            # TDOA + TOA
            elif self.used_ldp['tdoa'] and \
                    self.used_ldp['toa'] and \
                    not self.used_ldp['rss']:
                FIM = self._fim_TDOA(P) + self._fim_TOA(P)

            # TDOA + RSS
            elif self.used_ldp['tdoa'] and \
                    self.used_ldp['rss'] and \
                    not self.used_ldp['toa']:
                FIM = self._fim_TDOA(P) + self._fim_RSS(P)

            # TOA + RSS
            elif self.used_ldp['toa'] and \
                    self.used_ldp['rss'] and \
                    not self.used_ldp['tdoa']:
                FIM = self._fim_TOA(P) + self._fim_RSS(P)

            # TDOA + TOA + RSS
            else:
                FIM = self._fim_TDOA(P) + self._fim_TOA(P) +\
                    self._fim_RSS(P)

            return FIM


    def _fim_TDOA(self,P):
        """
        Compute the fisher information matrix in P for the given
        for TDOA observables

        Parameters
        ----------
            P : numpy.ndarray

        Returns
        -------
            FIM : numpy.ndarray

        """

        rn_tdoa = np.delete(self.an_tdoa, self.tdoa_ref, 1)
        rnr_tdoa = self.an_tdoa[
            :, self.tdoa_ref] * np.ones((3, self.Ntdoa - 1))
    #                rnr_tdoa = self.nodes['RNr_TDOA']
        tdoa_ns = self.tdoa
        tdoa_std = self.tdoa_std
        sh1 = np.shape(rn_tdoa)
        drg = self.c * tdoa_ns
        drg_std = self.c * tdoa_std
        FIM = np.zeros((sh1[0], sh1[0]))
        for i in range(sh1[1]):
            f1 = P - rn_tdoa[:, i:i + 1]
            f2 = P - rnr_tdoa[:, i:i + 1]
            pmrn = np.sqrt(np.dot(f1.T, f1))
            pmrnr = np.sqrt(np.dot(f2.T, f2))
            FIM += (1 / drg_std[i] ** 2) * np.dot(f1 /
                                                  pmrn - f2 / pmrnr, (f1 / pmrn - f2 / pmrnr).T)
        return FIM

    def _fim_TOA(self,P):
        """
        Compute the fisher information matrix in P for the given
        for TOA observables

        Parameters
        ----------
            P : numpy.ndarray

        Returns
        -------
            FIM : numpy.ndarray

        """
        rn_toa = self.an_toa
        toa_ns = self.toa
        toa_std = self.toa_std
        sh1 = np.shape(rn_toa)
        rg = self.c * toa_ns
        rg_std = self.c * toa_std
        FIM = np.zeros((sh1[0], sh1[0]))
        for i in range(sh1[1]):
            f1 = P - rn_toa[:, i:i + 1]
            FIM += np.dot(f1, f1.T) / \
                ((rg_std[i] ** 2) * np.dot(f1.T, f1))
            FIM += np.dot(f1, f1.T) / \
                ((rg_std[i] ** 2) * np.dot(f1.T, f1))
        return FIM

    def _fim_RSS(self,P):
        """
        Compute the fisher information matrix in P for the given
        for RSS observables

        Parameters
        ----------
            P : numpy.ndarray

        Returns
        -------
            FIM : numpy.ndarray

        """
        rn_rss = self.an_rss
        rss_db = self.rss
        rss_std = self.rss_std
        rss_np = self.rss_np
        sh1 = np.shape(rn_rss)
        d0 = self.d0
        pl0 = self.PL0
        S = (np.log(10) / 10) * rss_std / rss_np
        FIM = np.zeros((sh1[0], sh1[0]))
        for i in range(sh1[1]):
            f1 = P - rn_rss[:, i:i + 1]
            FIM += np.dot(f1, f1.T) / \
                ((S[0] ** 2) * (np.dot(f1.T, f1)) ** 2)
        return FIM


    def crb(self, P):
        """
        This method compute the cramer rao bound (CRB) at position P.

        Parameters
        ----------
            P : numpy.ndarray

        Returns
        -------
            CRB : float

        Examples
        --------

        .. plot::
            :include-source:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from pylayers.location.algebraic.algebraic import Algloc
        >>> from pylayers.location.observables import Observables
        >>> an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
               [2, 3, -1], [1, 0, 5], [1, 4, 0]])
        >>> an = an.T
        >>> bn = np.array([1, 1, 2.])
        >>> O_toa = Observables(an=an, bn=bn, mode='toa')
        >>> A = Algloc(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                       toa_std=O_toa.noise_model['std'], bnGT=bn)
        >>> A.crb(A.bnGT)
        """

        FIM = self.fim(P)
        CRB = np.sqrt(np.trace(nplg.pinv(FIM)))
        return CRB



if __name__ == "__main__":
    doctest.testmod()
