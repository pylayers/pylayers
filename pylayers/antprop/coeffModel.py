"""

Function SSH model
===================

.. autosummary::
   :toctree: generated/

   relative_error
   RepAzimuth1
   mode_energy
   mode_energy2
   level_energy
   modeMax
   lmreshape
   sshModel

"""
import doctest
import pdb
import numpy as np
import scipy as sp
import scipy.special as special
import matplotlib.pylab as plt
from numpy import zeros


def relative_error(Eth_original, Eph_original,Eth_model, Eph_model,theta, phi, dsf=1,kf=-1):
    """ calculate relative error betwee original and model
    Parameters
    ----------

    Eth_original : np.array
    Eph_original : np.array
    Eth_model    : np.array
    Eph_model    : np.array
    theta : np.array
    phi : np.phi
    dsf : int
        down sampling factor
    kf : int

    """

    st = np.sin(theta).reshape((len(theta), 1))
    #
    # Construct difference between reference and reconstructed
    #
    if kf<>-1:
        dTh = (Eth_model[kf, :, :] - Eth_original[kf, ::dsf, ::dsf])
        dPh = (Eph_model[kf, :, :] - Eph_original[kf, ::dsf, ::dsf])
        #
        # squaring  + Jacobian
        #
        dTh2 = np.real(dTh * np.conj(dTh)) * st
        dPh2 = np.real(dPh * np.conj(dPh)) * st

        vTh2 = np.real(Eth_original[kf, ::dsf, ::dsf] \
             * np.conj(Eth_original[kf, ::dsf, ::dsf])) * st
        vPh2 = np.real(Eph_original[kf, ::dsf, ::dsf] \
             * np.conj(Eph_original[kf, ::dsf, ::dsf])) * st

        mvTh2 = np.sum(vTh2)
        mvPh2 = np.sum(vPh2)

        errTh = np.sum(dTh2)
        errPh = np.sum(dPh2)
    else:
        dTh = (Eth_model[:, :, :] - Eth_original[:, ::dsf, ::dsf])
        dPh = (Eph_model[:, :, :] - Eph_original[:, ::dsf, ::dsf])
        #
        # squaring  + Jacobian
        #
        dTh2 = np.real(dTh * np.conj(dTh)) * st
        dPh2 = np.real(dPh * np.conj(dPh)) * st

        vTh2 = np.real(Eth_original[:, ::dsf, ::dsf] \
             * np.conj(Eth_original[:, ::dsf, ::dsf])) * st
        vPh2 = np.real(Eph_original[:, ::dsf, ::dsf] \
             * np.conj(Eph_original[:, ::dsf, ::dsf])) * st

        mvTh2 = np.sum(vTh2)
        mvPh2 = np.sum(vPh2)

        errTh = np.sum(dTh2)
        errPh = np.sum(dPh2)

    errelTh = (errTh / mvTh2)
    errelPh = (errPh / mvPh2)
    errel =( (errTh + errPh) / (mvTh2 + mvPh2))

    return(errelTh, errelPh, errel)


def RepAzimuth1 (Ec, theta, phi, th= np.pi/2,typ = 'Gain'):
    """
    Parameters
    ----------
    Ec 
    theta :
    phi
    th : 
    typ : string
        'Gain' 
    """

    pos_th = np.where(theta == th)[0][0]
    start = pos_th*len(phi)
    stop =  start + len(phi)
    if typ=='Gain':
        V   =  np.sqrt(np.real(Ec[0,:,start:stop]*np.conj(Ec[0,:,start:stop])+Ec[1,:,start:stop]*np.conj(Ec[1,:,start:stop])+Ec[2,:,start:stop]*np.conj(Ec[2,:,start:stop])))
    if typ=='Ex':
        V = np.abs(Ec[0,:,start:stop])
    if typ=='Ey':
        V = np.abs(Ec[1,:,start:stop])
        
    if typ=='Ez':
        V = np.abs(Ec[2,:,start:stop])
    
    VdB = 20*np.log10(V)
    VdBmin = -40
    VdB = VdB - VdBmin
    V = VdB
    #plt.polar(phi,V)
    #plt.title('theta = '+str(th)) 
    return V

def mode_energy(C,M,L =20, ifreq = 46):
    
    """
    shape C = (dim = 3,Ncoef = (1+L)**2)
    """
    Em = []
    Lc = (1+L)**2
    for m in range(M+1):
        im = m*(2*L+3-m)/2
        bind = (1+L)*(L+2)/2 + im-L-1 
        #if ifreq > 0:
        if m == 0:
            em  = np.sum(np.abs(C[:,ifreq,im:im+L-m+1])**2)
        else:
            em  = np.sum(np.abs(C[:,ifreq,im:im+L-m+1])**2) + np.sum(np.abs(C[:,ifreq,bind:bind + L-m+1])**2)
        Et = np.sum(np.abs(C[:,ifreq,:])**2)
       
        Em.append(em)    
    return  np.array(Em)/Et

def mode_energy2(A,m, ifreq=46, L= 20):
    
    cx = lmreshape(A.S.Cx.s2)
    cy = lmreshape(A.S.Cy.s2)
    cz = lmreshape(A.S.Cz.s2)
    if ifreq >0:
        em = np.sum(np.abs(cx[ifreq,:,L+m])**2+np.abs(cy[ifreq,:,L+m])**2+np.abs(cz[ifreq,:,L+m])**2)
        Et = np.sum(np.abs(cx[ifreq])**2+np.abs(cy[ifreq])**2+np.abs(cz[ifreq])**2)
    return em/Et

def level_energy(A,l, ifreq = 46,L=20):
    cx = lmreshape(A.S.Cx.s2)
    cy = lmreshape(A.S.Cy.s2)
    cz = lmreshape(A.S.Cz.s2)
    if ifreq >0:
        el = np.sum(np.abs(cx[ifreq,l,:])**2+np.abs(cy[ifreq,l,:])**2+np.abs(cz[ifreq,l,:])**2)
        Et = np.sum(np.abs(cx[ifreq])**2+np.abs(cy[ifreq])**2+np.abs(cz[ifreq])**2)
    return el/Et

def modeMax(coeff,L= 20, ifreq  = 46):
    
    Em_dB = 20*np.log10(mode_energy(C = coeff,M = L, ifreq=ifreq))
    
    max_mode = np.where(Em_dB <-20 )[0][0]-1
    return max_mode

def lmreshape(coeff,L= 20):
    
    sh = coeff.shape
    
    coeff_lm = zeros(shape = (sh[0],1+L, 1+2*L), dtype = complex )
    
    for m in range(0,1+L):        
        im = m*(2*L+3-m)/2
        coeff_lm[:,m:L+1,L+m] = coeff[:,im:im +L+1-m]
       
    for m in range(1,L):
        im = m*(2*L+3-m)/2
        bind = (1+L)*(L+2)/2 + im-L-1 
        coeff_lm[:,m:L+1,L-m]= coeff[:,bind: bind + L-m+1]
    
    return coeff_lm 



def sshModel(c,dmm, L = 20, ifreq = 46):
    """ modify the free space pattern of a body mounted antenna

    Parameters
    ----------
    c : SHCoeff
        Free space coefficient of the antenna
    dmm : float 
        distance in millimeters between antenna and dielectric cylinder
    L : 
    ifreq : frequency index
    model : dict 
        dictionnary with model parameters a,b,a',b',a",b"

    Notes 
    -----
    
    This model is presented in "Modelling of UWB Antenna Perturbed by Human Phantom in Spherical Harmonics Space"
    
    EUCAP 2014 
    
    See Also
    --------
    
    pylayers.antprop.spharm
    
    """
    defaults = { L : 20 ,
                 ifreq: 46, 
                 model : { a : 0.002,
                           b : 0.55,
                           ap : -0.002,
                           bp : 1.55,
                           app : -0.006,
                           bpp : 1.22
                         }
                }
    
    for k in defaults.keys():
        if k not in kwargs:
            kwargs[k]=defaults[k]
    
    mdl = kwargs['model']
    L = kwargs['L']
    ifreq =kwargs['ifreq']
            
    Lc = (1+L)**2
    sh = np.shape(c)
    cm = np.zeros(shape = sh , dtype = complex)
    m0 =  modeMax(c, ifreq= ifreq, L= 20)
    im0 = m0*(2*L+3-m0)/2
    
    # m0 corresponds to Mfs
    # M corresponds to Mp
    # alpha corresponds to phim 
    
    # formula (10) [mhed14]
    M = m0 + int(0.06*dmm) + 4
    #a0 = 0.002*dmm+0.55
    a0 = mdl['a']*dmm+mdl['b']
    #am = -0.002*dmm + 1.55
    am = mdl['ap']*dmm + mdl['bp']
    #alpha = -0.006*dmm+1.22
    alpha = mdl['app']*dmm+mdl['bpp']
    
    for m in range(0,m0):

        im = m*(2*L+3-m)/2
        if m == 0:
            dephm = 0
            cm[:,:,im: im+L+1-m] = a0*c[:,:,im: im+L+1-m]
        else:
            dephm = (m-m0)*alpha
            cm[:,:,im: im+L+1-m] = a0*c[:,:,im: im+L+1-m]*np.exp(1j*dephm)
            bind = (1+L)*(L+2)/2 + im-L-1
            cm[:,:,bind: bind + L-m+1] = ((-1)**m)*cm[:,:,im: im+L+1-m]


    for m in range(m0,M+1):

        dephm = (m-m0)*alpha
        if m == m0:

            im = m*(2*L+3-m)/2
            cm[:,:,im: im+L+1-m] = (am/(m-m0+1))*c[:,:,im0 : im0+L-m+1]*np.exp(1j*dephm)
            bind = (1+L)*(L+2)/2 + im -L-1
            cm[:,:,bind: bind + L-m+1] =  ((-1)**m)*(cm[:,:,im: im+L+1-m])

        else:

            im = m*(2*L+3-m)/2
            cm[:,:,im: im+L+1-m] = (am/(m-m0+1))*c[:,:,im0 : im0+L-m+1]*np.exp(1j*dephm)
            bind = (1+L)*(L+2)/2 + im -L-1
            cm[:,:,bind: bind + L-m+1] =  ((-1)**m)*(cm[:,:,im: im+L+1-m])

    cm[0:2] = c[0:2]
    return cm
    
    
def sshModel2(c,d, L = 20, ifreq = 46, alpha =0):
    """
    sshModel

    """


    Lc = (1+L)**2
    sh = np.shape(c)
    cm = np.zeros(shape = sh , dtype = complex)
    m0 =  modeMax(c, ifreq= ifreq, L= 20)
    im0 = m0*(2*L+3-m0)/2
    M = m0 + int(0.06*d) + 4
    a0 = 0.002*d+0.55
    am = -0.002*d + 1.55
    #alpha = -0.009*d+1.84 +np.pi/2

    for m in range(0,m0):

        im = m*(2*L+3-m)/2
        if m == 0:
            dephm = m*alpha
            cm[:,:,im: im+L+1-m] = a0*c[:,:,im: im+L+1-m]
        else:
            dephm = (m-m0)*alpha
            cm[:,:,im: im+L+1-m] = a0*c[:,:,im: im+L+1-m]*np.exp(1j*dephm)
            bind = (1+L)*(L+2)/2 + im-L-1
            cm[:,:,bind: bind + L-m+1] = ((-1)**m)*cm[:,:,im: im+L+1-m]


    for m in range(m0,M+1):

        dephm = m*alpha
        if m == m0:

            im = m*(2*L+3-m)/2
            cm[:,:,im: im+L+1-m] = (am/(m-m0+1))*c[:,:,im0 : im0+L-m+1]*np.exp(1j*dephm)
            bind = (1+L)*(L+2)/2 + im -L-1
            cm[:,:,bind: bind + L-m+1] =  ((-1)**m)*(cm[:,:,im: im+L+1-m])

        else:

            im = m*(2*L+3-m)/2
            cm[:,:,im: im+L+1-m] = (am/(m-m0+1))*c[:,:,im0 : im0+L-m+1]*np.exp(1j*dephm)
            bind = (1+L)*(L+2)/2 + im -L-1
            cm[:,:,bind: bind + L-m+1] =  ((-1)**m)*(cm[:,:,im: im+L+1-m])

    #cm[0:2] = c[0:2]
    return cm





if (__name__=="__main__"):
    doctest.testmod()


