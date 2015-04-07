import scipy.stats as st
import numpy as np


def getchannel(emplacement = 'trunku',intersection = 1):
    """ get channel

    Parameters
    ----------

    emplacement : 'trunku' | 'thighr' | 'forearm' | 'calfr'
    intersection : 1 = LOS 0 : NLOS

    Returns
    -------

    alphak : np.array
    tauk   : np.array

    Notes
    -----

    See : "Delay Dispersion of the On-Body Channel" (Raffaele Derrico, Laurent Ouvry)

    """


    pdp = {'trunku':{'los':{'gamma':25.9,'gamma0':3.37,'k':26,'lambda-1':2.13,'sigma':4.02},
                'trans':{'gamma':'none','gamma0':'none','k':'none','lambda-1':'none','sigma':'none'},
                'nlos':{'gamma':'none','gamma0':'none','k':'none','lambda-1':'none','sigma':'none'}},
       'thighr':{'los':{'gamma':24.9,'gamma0':2.1,'k':7,'lambda-1':3.13,'sigma':4.62},
                'trans':{'gamma':46.7,'gamma0':1.6,'k':19,'lambda-1':2.46,'sigma':4.58},
                'nlos':{'gamma':58.9,'gamma0':1.6,'k':30,'lambda-1':2.35,'sigma':4.32}},
       'forearmr':{'los':{'gamma':22.3,'gamma0':2.3,'k':5,'lambda-1':2.27,'sigma':4.55},
                'trans':{'gamma':45.5,'gamma0':2.7,'k':13,'lambda-1':2.69,'sigma':4.28},
                'nlos':{'gamma':63,'gamma0':2.5,'k':27,'lambda-1':2.39,'sigma':4.31}},
       'calfr':{'los':{'gamma':29,'gamma0':3.3,'k':10,'lambda-1':2.98,'sigma':4.57},
                'trans':{'gamma':40.5,'gamma0':2.2,'k':23,'lambda-1':2.37,'sigma':4.61},
                'nlos':{'gamma':55,'gamma0':1.4,'k':32,'lambda-1':2.21,'sigma':4.65}}
       }

    g0 = {'trunku':{'mu0':-52.62,'sigma0':4.35},
        'thighr':{'mu0':-63.30,'sigma0':2.31},
        'forearmr':{'mu0':-59.96,'sigma0':3.28},
        'calfr':{'mu0':-62.93,'sigma0':1.69},
        }

    condition = 'nlos'
    if intersection == 1:
        condition = 'los'

    if emplacement == 'trunku':
        condition = 'los'

    #number of paths
    K = pdp[emplacement][condition]['k']

    #delay
    Lambda = 1./pdp[emplacement][condition]['lambda-1']

    Tk = st.expon(0,Lambda)
    sampleTk = Tk.rvs(K)
    tauk = np.cumsum(sampleTk)

    #exponential decay
    gamma  = pdp[emplacement][condition]['gamma']

    #rician factor
    gamma0 = pdp[emplacement][condition]['gamma0']

    #path amplitude
    alpha_k_dB_m = gamma0 + 10*np.log10(np.exp(-tauk)/gamma)

    sigmas = pdp[emplacement][condition]['sigma']
    alpha_k_dB = st.norm(alpha_k_dB_m,sigmas)

    alphak = 10**(alpha_k_dB.rvs(size = len(tauk))/10)
    alphak = alphak/np.sqrt(sum(abs(alphak)**2))

    # mean channel gain
    mu0 = g0[emplacement]['mu0']-5
    sigma0 = g0[emplacement]['sigma0']
    GdB_dist = st.norm(mu0,sigma0)
    GdB = GdB_dist.rvs()
    G = 10**(GdB/10.0)
    alphak = np.sqrt(G)*alphak


    alphak = intersection*alphak

    return alphak,tauk







