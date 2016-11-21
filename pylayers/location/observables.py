from __future__ import print_function
import numpy as np
import scipy as sp
from copy import copy
import pylayers.antprop.loss as plm
import matplotlib.pyplot as plt

class Observables(object):
    """ Generate observables for localization prupose
    """
    def __init__(self, an=10 * sp.rand(3, 5), bn=5 * sp.rand(3,4),mode = 'toa'):
        """
        Init

        Parameters
        ----------
        an : ndarray
            anchors node(3 x Na)
        bn : ndarray
            blind node (3 x Nb)
        mode : str
            "toa" : compute ranges
            "tdoa" : compute diff of ranges
            "rss" : compute recived power

        Notes
        -----
        self.Na : Number of anchor nodes (Na)
        self.Nb : Number of blind nodes  (Nb)

        self.dist : distances matrix (Nb x Na)

        self.rng : range matrix (Nb x Na)
        self.drng : difference of ranges matrix (Na x Nb x Na)

        self.rp : received power ( Nb x Na)
        self.rp_model : power model

        self.noise : noise matrix ()
        self.noise_model : noise matrix ()

        """

        if len(an.shape) > 2:
            raise AttributeError('Anchors \'an\' shape must be (3 x Na)')

        if len(bn.shape) > 2:
            raise AttributeError('Blind nodes \'bn\' shape must be (3 x Nb)')

        if len(an.shape) == 1:
            an = np.array([an])
            if an.shape[0] == 1:
                an = an.T
        if len(bn.shape) == 1:
            bn = np.array([bn])
            if bn.shape[0] == 1:
                bn = bn.T


        if an.shape[0] == 2:
            an = np.vstack((an, np.zeros(an.shape[1])))
        elif an.shape[0] != 3:
            raise AttributeError('Anchors an first dimension reserved to space\
                                  (x,y,z) coordinates')

        if bn.shape[0] == 2:
            bn = np.vstack((bn, np.zeros(bn.shape[1])))
        elif bn.shape[0] != 3:
            raise AttributeError('Blind nodes bn first dimension reserved to \
                                  space (x,y,z) coordinates')



        self._implemented_law = ['norm']


        self.an = an
        self.bn = bn

        self.Na = self.an.shape[1]
        self.Nb = self.bn.shape[1]

        self.mode = mode

        self.compute_distances()
        if mode.lower() == 'toa':
            self.compute_ranges()
        if mode.lower() == 'tdoa':
            self.an_ref = 0
            self.compute_diff_distances()
            self.compute_diff_ranges()
            self._change_an_ref()
        if mode.lower() == 'rss':
            self.compute_rpower()
        self.param_noise(param={})


    @property
    def an_ref(self):
        return self._an_ref

    @an_ref.setter
    def an_ref(self, value):
        if hasattr(self, 'an_ref'):
            self._change_an_ref(value)
        else:
            # first call from compute_power
            self._an_ref = value

    # @property
    # def rp_model(self):
    #     return self._rp_model

    # @rp_model.setter
    # def rp_model(self, value):
    #     if hasattr(self, 'rp_model'):
    #         if self._rp_model != value:
    #             self.compute_rpower(param=value)
    #     else:
    #         # first call from compute_power
    #         self._rp_model = value

    # @property
    # def noise_model(self):
    #     return self._noise_model

    # @noise_model.setter
    # def noise_model(self, value):
    #     if hasattr(self, 'noise_model'):
    #         if not isinstance(value['law'],list):
    #             if self.mode != 'tdoa':
    #                 value['law']=[value['law']]*self.Na
    #             else: 
    #                 value['law']=[value['law']]*self.Na-1
    #         idem = np.alltrue(np.array([i in self._implemented_law for i in value['law']]))
    #         # check 
    #         if not idem:
    #             raise AttributeError('A specified law of noise model is not yet implemented')
    #         else:
    #             self._noise_model = value
    #             self.generate_noise(param=value)
    #     else:
    #         # first call from compute_power
    #         self._noise_model = value


    def __repr__(self):

        s = 'self.Na : Number of anchor nodes (Na)'
        s = s + '\n' + 'self.Nb : Number of blind nodes  (Nb)'

        s = s + '\n' + 'self.dist : distances matrix (Nb x Na)'

        if self.mode =='toa':
            s = s + '\n' + 'self.rng : range matrix (Nb x Na)'

        if self.mode =='tdoa':
            s = s + '\n' + 'self.drng : difference of ranges matrix (Na x Nb x Na)'

        if self.mode =='rss':
            s = s + '\n' + 'self.rp : received power ( Nb x Na)'
            s = s + '\n' + 'self.rp_model : power model'

        s = s + '\n' + 'self.noise : noise matrix (Nb x Na)'
        s = s + '\n' + 'self.noise_model : noise param dict'

        s = s + '\n' + 'self.generate_noise_samples() : update self.noise'


        s = s + '\n\n' +str(self.Na) + ' Anchors:\n'
        if self.mode =='tdoa':
            s = s + 'reference node: ' + str(self.an_ref) + '\n'
        s = s + '--------------\n\n'

        for a in self.an:
            s = s + str(a) + "\n"


        s = s + '\n' + str(self.Nb) + ' Blind nodes:\n'
        s = s + '--------------\n\n'

        for b in self.bn:
            s = s + str(b) + "\n"

        s = s + '\n' +'noise paramuration:\n'
        s = s + '---------------------\n\n'

        for n in self.noise_model:
            s = s + n + ' : ' + str(self.noise_model[n ])+ "\n"

        s = s + '\n' +'noise realizations:\n'
        s = s + '-------------------\n\n'
        s = s + str(self.noise) + "\n"

        return s

    def compute_distances(self):
        """
            Compute ditance between all anchors an and all blind nodes bn

        Return
        ------

        self.dist : nd.array
            (Nb x Nc)

        """

        self.dist = np.sqrt(np.sum((self.an[:, None, :] - self.bn[:, :, None])**2, axis=0))

    def compute_ranges(self):
        """
        Compute range in nanoseconds between all anchors an
        and all blind nodes bn

        Return
        ------

        self.rng : ndarray
            range nd array(Nb x Na)
        """

        if not hasattr(self, 'dist'):
            self.compute_distances()
        self.rng = self.dist / 0.3

    def compute_diff_distances(self):
        """
        Compute difference of ditance

        Return
        ------

        self.ddist : ndarray
            difference of distances for each node as refernce:
            (Nb x Na-1)
        """
        ddist = np.ndarray(shape=(0, self.Nb, self.Na))
        for a in xrange(self.Na):
            diff = self.dist[:, :] - self.dist[:, a][:, None]
            ddist = np.vstack((ddist, diff[None, ...]))
        self._ddist = ddist

    def compute_diff_ranges(self):
        """
        Compute difference of ditance

        Return
        ------

        self.drng : ndarray
            difference of ranges in nanoseconds for each node as refernce:
            (Nb x Na-1)
        """
        if not hasattr(self,'ddist'):
            self.compute_diff_distances()
        self._drng = self._ddist/0.3

    def _change_an_ref(self, an_ref = 0):
        """ 
        Changing the reference node of TDOA

        Parameters
        ----------

        an_ref : int
            index of an column indicating the refeernce node
        """

        self.ddist = np.delete(self._ddist[an_ref,:,:],an_ref,1)
        self.drng = np.delete(self._drng[an_ref,:,:],an_ref,1)

    def compute_rpower(self,param={}):
        """
        Compute received power given a model

        Parameters
        ----------
        param : dict
            for Pathloss shadowing :
                param['model']='PL'
                param['d0'] : reference distance
                param['fGHz'] : frequency in GHz
                param['pl_exp'] : pathloss exponent

        """

        implemented_model=['PL']

        if param == {}:
            param['model'] = 'PL'
            param['d0'] = 1.
            param['fGHz'] = 2.4
            param['pl_exp'] = 2.
        else:
            if isinstance(param, dict):
                if param.has_key('model'):
                    if param['model'] in implemented_model:
                        pass
                    else:
                        raise AttributeError('model ' + str(param['model']) +
                                              ' is not yet implemented')
                else:
                    raise AttributeError('param dict has no \'model\' key')
            else:
                raise AttributeError('param must be a dict instance')

        if param['model'] == 'PL':
            self.rp = -plm.PL0(param['fGHz'], param['d0']) +\
                          10 * param['pl_exp'] * np.log10(self.dist / param['d0'])
            self.rp_model = param

    def param_noise(self, param={}):
        """

        Create a noise matrix for obserables

        Parameters
        ----------

        model : dict
            param['law'] = name of distrib
            'norm'

        Returns
        -------
        """

        if self.mode != 'tdoa':
            Na = self.Na
        else:
            Na = self.Na -1

        if param == {}:
            param['law'] = 'norm'
            param['mean'] = 0.
            param['std'] = 0.2
        else:
            if isinstance(param, dict):
                if param.has_key('law'):
                    if not isinstance(param['law'],list):
                        param['law']=[param['law']]*Na
                    for na in range(Na):
                        if param['law'][na] in self._implemented_law:
                            pass
                        else:
                            raise AttributeError('law ' + str(param['law'][na]) +
                                                  ' is not yet implemented')
                else:
                    raise AttributeError('param dict has no \'law\' key')
            else:
                raise AttributeError('param must be a dict instance')

        if not isinstance(param['law'],list):
            param['law']=[param['law']]*Na
        if not isinstance(param['mean'],list):
            param['mean']=[param['mean']]*Na
        if not isinstance(param['std'],list):
            param['std']=[param['std']]*Na

        self.law = [np.nan]*Na

        for na in range(Na):
            if param['law'][na] == 'norm':
                self.law[na] = sp.stats.norm(loc=param['mean'][na], scale=param['std'][na])


        self.noise_model = param
        del param
        self.generate_noise_samples()

    def generate_noise_samples(self):
        """
        Generate new noise samples relying on the law paramuration 
        setup in self.paramure_noise
        """
        if self.mode != 'tdoa':
            Na = self.Na
        else:
            Na = self.Na -1

        self.noise = np.ndarray((self.Nb,Na))

        for na in range(Na):
            self.noise[:,na] = self.law[na].rvs((self.Nb))

    def show(self, **kwargs):
        """
            Show scene
        """
        defaults = {'fig': [],
                    'ax': []
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] == []:
            fig = plt.figure(figsize=(5, 5))
        else:
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']


        ax.scatter(self.an[0,:],self.an[1,:], c='k', marker='o', label='Anchors')
        ax.scatter(self.bn[0,:],self.bn[1,:], c='r', marker='o', label='Blind nodes')
        plt.legend()

        return fig, ax

if (__name__ == "__main__"):
    import copy
    O = Observables()