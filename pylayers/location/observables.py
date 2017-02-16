from __future__ import print_function
import numpy as np
import scipy as sp
import pylayers.antprop.loss as plm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Observables(object):
    """ Generate observables for localization prupose
    """
    def __init__(self, an=10 * sp.rand(3, 5), bn=5 * sp.rand(3,4)):
        """
        Init

        Parameters
        ----------
        an : ndarray
            anchors node(3 x Na)
        bn : ndarray
            blind node (3 x Nb)


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

        if len(an.shape) == 1:
            an = np.array([an])
        if len(bn.shape) == 1:
            bn = np.array([bn])



        self.an = an
        self.bn = bn

        self.Na = self.an.shape[1]
        self.Nb = self.bn.shape[1]

        self.compute_distances()
        self.compute_diff_distances()
        self.compute_ranges()
        self.compute_diff_ranges()
        self.compute_rpower()

        self.config_noise()

    @property
    def rp_model(self):
        return self._rp_model

    @rp_model.setter
    def rp_model(self, value):
        if hasattr(self, 'rp_model'):
            if self._rp_model != value:
                self.compute_rpower(config=value)
        else:
            # first call from compute_power
            self._rp_model = value

    @property
    def noise_model(self):
        return self._noise_model

    @noise_model.setter
    def noise_model(self, value):
        if hasattr(self, 'noise_model'):
            idem = np.alltrue(np.array([i in value.items() for i in self.noise_model.items()]))
            # check if any key/values of self.noise != new assigned value
            if not idem:
                self._noise_model = value
                self.create_noise(config=value)
        else:
            # first call from compute_power
            self._noise_model = value


    def __repr__(self):

        s = str(self.Na) + ' Anchors:\n'
        s = s + '--------\n\n'

        for a in self.an:
            s = s + str(a) + "\n"

        s = s + '\n' + str(self.Nb) + ' Blind nodes:\n'
        s = s + '------------\n\n'

        for b in self.bn:
            s = s + str(b) + "\n"

        s = s + '\n'
        s = s + '\n' + 'self.Na : Number of anchor nodes (Na)'
        s = s + '\n' + 'self.Nb : Number of blind nodes  (Nb)'

        s = s + '\n' + 'self.dist : distances matrix (Nb x Na)'
        s = s + '\n' + 'self.rng : range matrix (Nb x Na)'
        s = s + '\n' + 'self.drng : difference of ranges matrix (Na x Nb x Na)'

        s = s + '\n' + 'self.rp : received power ( Nb x Na)'
        s = s + '\n' + 'self.rp_model : power model'

        s = s + '\n' + 'self.noise : noise matrix ()'
        s = s + '\n' + 'self.noise_model : noise matrix ()'

        s = s + '\n' + 'self.generate_noise_samples() : update self.noise'

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
            (Na x Nb x Na)
        """
        ddist = np.ndarray(shape=(0, self.Nb, self.Na))
        for a in xrange(self.Na):
            diff = self.dist[:, a][:, None] - self.dist[:, :]

            ddist = np.vstack((ddist, diff[None, ...]))
        self.ddist = ddist

    def compute_diff_ranges(self):
        """
        Compute difference of ditance

        Return
        ------

        self.drng : ndarray
            difference of ranges in nanoseconds for each node as refernce:
            (Na x Nb x Na)
        """
        if not hasattr(self,'ddist'):
            self.compute_diff_distances()
        self.drng = self.ddist/0.3

    def compute_rpower(self,config={}):
        """
        Compute received power given a model

        Parameters
        ----------
        config : dict
            for Pathloss shadowing :
                config['model']='PL'
                config['d0'] : reference distance
                config['fGHz'] : frequency in GHz
                config['pl_exp'] : pathloss exponent

        """

        implemented_model=['PL']

        if config == {}:
            config['model'] = 'PL'
            config['d0'] = 1.
            config['fGHz'] = 2.4
            config['pl_exp'] = 2.
        else:
            if isinstance(config, dict):
                if config.has_key('model'):
                    if config['model'] in implemented_model:
                        pass
                    else:
                        raise AttributeError('model ' + str(config['model']) +
                                              ' is not yet implemented')
                else:
                    raise AttributeError('config dict has no \'model\' key')
            else:
                raise AttributeError('config must be a dict instance')

        if config['model'] == 'PL':
            self.rp = -plm.PL0(config['fGHz'], config['d0']) +\
                          10 * config['pl_exp'] * np.log10(self.dist / config['d0'])
            self.rp_model = config

    def config_noise(self, config={}):
        """

        Create a noise matrix for obserables

        Parameters
        ----------

        model : dict
            config['law'] = name of distrib
            'norm'

        Returns
        -------
        """

        implemented_law = ['norm']

        if config == {}:
            config['law'] = 'norm'
            config['mean'] = 0.
            config['std'] = 2.
        else:
            if isinstance(config, dict):
                if config.has_key('law'):
                    if config['law'] in implemented_law:
                        pass
                    else:
                        raise AttributeError('law ' + str(config['law']) +
                                              'is not yet implemented')
                else:
                    raise AttributeError('config dict has no \'law\' key')
            else:
                raise AttributeError('config must be a dict instance')

        if config['law'] == 'norm':
            self.law = sp.stats.norm(loc=config['mean'], scale=config['std'])

        self.noise_model = config
        self.generate_noise_samples()

    def generate_noise_samples(self):
        """
        Generate new noise samples relying on the law configuration 
        setup in self.configure_noise
        """

        self.noise = self.law.rvs((self.Nb, self.Na))

    def show(self, **kwargs):
        """
            Show scene
        """
        defaults = {'fig': [],
                    'ax': [],
                    'mode3d':False ,
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] == []:
            fig = plt.figure(figsize=(5, 5))
        else:
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            if kwargs['d3']:
                ax = fig.add_subplot(111,projection='3d')
            else:
                ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']


        if kwargs['d3']:
            ax.scatter(self.an[0,:],self.an[1,:],self.an[2,:], c='k', marker='o', label='Anchors')
            ax.scatter(self.bn[0,:],self.bn[1,:],self.bn[2,:], c='r', marker='o', label='Blind nodes')
        else:
            ax.scatter(self.an[0,:],self.an[1,:], c='k', marker='o', label='Anchors')
            ax.scatter(self.bn[0,:],self.bn[1,:], c='r', marker='o', label='Blind nodes')

        plt.legend()

        return fig, ax

if (__name__ == "__main__"):

    O = Observables()
