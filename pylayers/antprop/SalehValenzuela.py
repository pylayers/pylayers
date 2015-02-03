import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import pylayers.signal.bsignal as bs

def SalehValenzuela(**kwargs):
    """ generic Saleh and Valenzuela Model

    Parameters
    ----------

    Lam : clusters Poisson Process parameter (ns)
    lam : rays Poisson Process parameter (ns)
    Gam : clusters exponential decay factor
    gam : rays exponential decay factor
    tauM : maximum delay


    """
    defaults = { 'Lam' : 10.,
                 'lam' : 5.,
                 'Gam' : 30.,
                 'gam' : 5. ,
                 'tauM': 1000.}

    for k in defaults:
        if k not in kwargs:
            kwargs[k]=defaults[k]

    Lam = kwargs['Lam']
    lam = kwargs['lam']
    Gam = kwargs['Gam']
    gam = kwargs['gam']
    tauM = kwargs['tauM']
    Nc = tauM/Lam
    Nr = tauM/lam

    p1 = st.poisson(Lam)
    p2 = st.poisson(lam)
    # cluster time of arrival
    tc = np.cumsum(p1.rvs(Nc))
    tauc = np.kron(tc,np.ones((1,Nr)))[0,:]
    # rays time of arrival
    taur = np.cumsum(p2.rvs((Nr,Nc)),axis=0).ravel()
    # exponential decays of cluster and rays
    etc = np.exp(-tauc/(1.0*Gam))
    etr = np.exp(-taur/(1.0*gam))
    et = etc*etr
    tau = tauc+taur
    # reordering in delay domain
    u = np.argsort(tau)
    taus = tau[u]
    ets = et[u]
    # limiting in delay domain
    v = np.where(taus<tauM)[0]
    taus = taus[v]
    ets = ets[v]
    SVir = bs.Bsignal(taus,ets)
    return(SVir)
