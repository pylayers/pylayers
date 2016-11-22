# -*- coding:Utf-8 -*-
import numpy as np

from pylayers.location.observables import Observables
from pylayers.location.geometric.constraints.cla import CLA
from pylayers.location.geometric.constraints.rss import RSS
from pylayers.location.geometric.constraints.toa import TOA
from pylayers.location.geometric.constraints.tdoa import TDOA
from pylayers.location.geometric.constraints.exclude import Exclude
from pylayers.location.algebraic.algebraic import Algloc


class Localization(object):
    """ Handle localization engine of agents

    Attributes
    ----------

    args
    config
    cla
    algloc
    idx

    """

    def __init__(self, **kwargs):
        """
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

        self.alg = Algloc(**kwargs)
        kwargs.update(self.alg.__dict__)

        for k in kwargs:
            setattr(self, k, kwargs[k])

        self.geo = CLA()
        CST = 0
        for k in range(self.Ntoa):
            self.geo.append(
                TOA(id=CST, p=self.an_toa[:, k], value=self.toa[k], std=self.toa_std[k]))
            CST += 1
        for k in range(self.Ntdoa - 1):
            if k != self.tdoa_ref:
                an = self.an_tdoa[np.ix_([0, 1, 2], [self.tdoa_ref[0], k])].T
                self.geo.append(
                    TDOA(id=CST, p=an, value=self.tdoa[k], std=self.tdoa_std[k]))
                CST += 1
        for k in range(self.Nrss):
            model = {'rss_np': self.rss_np[
                k], 'PL0': self.PL0, 'd0': self.d0, 'Rest': self.Rest}
            self.geo.append(
                RSS(id=CST, p=self.an_rss[:, k], value=self.rss[k], std=self.rss_std[k], model=model))
            CST += 1

    # def _update_used_ldp(self):
    #     """ check and update  modification in used _ldp have been made
    #         and apply change
    #     """
    #     self.alg.used_ldp.update(self.used_ldp)

    #     usable = np.array(self.geo.usable)
    #     for ldp in self.used_ldp:
    #         ug = np.where(np.array(self.geo.type) == ldp.upper())[0]
    #         for u in ug:
    #             self.geo.c[u].usable = self.used_ldp[ldp]
    #     self.geo.update()
    #     import ipdb
    #     ipdb.set_trace()

    def locate(self, mode=['ls', 'wls', 'ml', 'geo', 'crb']):
        """ Perform localization with given mode

        Parameters
        ----------
            mode = str | list
                'ls': algebraic localization using Least Square
                'wls': algebraic localization using Weighted Least Square
                'ml': algebraic localization using Maximum Likelihood
                'geo'  : geometric localization based on RGPA
                'crb' : Carmer-Rao bound

        Examples
        --------

        >>> from pylayers.location.observables import Observables
        >>> from pylayers.location.localization import *
        >>> import matplotlib.pyplot as plt
        >>> O = Observables()
        >>> an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
                        [2, 3, -1], [1, 0, 5], [1, 4, 0]])
        >>> an = an.T
        >>> bn = np.array([1, 1, 2.])
        >>> O_toa = Observables(an=an, bn=bn, mode='toa')
        >>> O_tdoa = Observables(an=an, bn=bn, mode='tdoa')
        >>> O_rss = Observables(an=an, bn=bn, mode='rss')
        >>> L = Localization(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                          toa_std=O_toa.noise_model['std'],
                          an_tdoa=O_tdoa.an, tdoa=O_tdoa.drng,
                          tdoa_ref=O_tdoa.an_ref, tdoa_std=0.05,
                          an_rss=O_rss.an, rss=O_rss.rp, rss_std=O_rss.noise_model[
                          'std'], rss_np=2., PL0=40.04, d0=1., bnGT=bn)
        >>> L.locate()
        >>> L.show()
        >>> plt.show()

        """

        if isinstance(mode, str):
            mode = [mode]

        if 'ls' in mode:
            self.pe_ls = self.alg.locate('ls')
        if 'wls' in mode:
            self.pe_wls = self.alg.locate('wls')
        if 'ml' in mode:
            self.pe_ml = self.alg.locate('ml')
        if 'crb' in mode:
            if self.bnGT.shape[1] != 0:
                self.crb = self.alg.crb(self.bnGT)
            else:
                Warning('CRB not computed, because self.bnGT void')
        if 'geo' in mode:
            self.geo.compute()
            self.pe_geo = self.geo.pe.reshape(3, 1)

    def show(self, **kwargs):

        defaults = {'legend': True,
                    }
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        legend = kwargs['legend']
        kwargs['legend'] = False

        fig, ax = self.alg.show(**kwargs)
        if hasattr(self, 'pe_geo'):
            X = np.concatenate([self.pe_geo, self.bnGT], axis=1)
            ax.plot(self.pe_geo[0, :], self.pe_geo[
                    1, :], self.pe_geo[2, :], "rD",label='geometric')
            ax.plot(X[0, :], X[1, :], X[2, :], linewidth=0.5, color='k')
        if legend:
            ax.legend()
        return fig, ax


if (__name__ == "__main__"):
    pass
