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
        # for k in range(self.Nrss):
        #     model =
        #     self.geo.append(RSS(id=CST, p=self.an_rss[:,k], value=self.rss[k], std=self.rss_std[k]))
        #     CST += 1

    def locate(self, mode=['ls', 'wls', 'ml', 'geo', 'crb']):
        """ Perform localization with given mode

            mode = 'ls' | 'wls' | 'ml' | 'geo'  
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
            if self.bnGT.shape[1]!=0:
                self.crb = self.alg.crb(self.bnGT)
            else:
                Warning('CRB not computed, because self.bnGT void')
        if 'geo' in mode:
            self.geo.compute()
            self.pe_geo = self.geo.pe.reshape(3,1)


if (__name__ == "__main__"):
    O = Observables()
    an = np.array([[0, 1, 2.], [0, 3, 1], [2, 1, 3],
                   [2, 3, -1], [1, 0, 5], [1, 4, 0]])
    an = an.T

    bn = np.array([1, 1, 2.])

    O_toa = Observables(an=an, bn=bn, mode='toa')
    O_tdoa = Observables(an=an, bn=bn, mode='tdoa')

    L = Localization(an_toa=O_toa.an, toa=O_toa.rng + O_toa.noise,
                     toa_std=O_toa.noise_model['std'],
                     an_tdoa=O_tdoa.an, tdoa=O_tdoa.drng,
                     tdoa_ref=O_tdoa.an_ref, tdoa_std=0.05,
                     bnGT=bn)
