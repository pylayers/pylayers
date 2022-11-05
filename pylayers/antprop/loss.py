# -*- coding:Utf-8 -*-
#from numpy import *
from __future__ import print_function
"""
.. currentmodule:: pylayers.antprop.loss

.. autosummary::
    :members:

"""
import doctest
import os
import glob
import doctest
import logging
import numpy as np
from scipy import io
import matplotlib.pylab as plt
import pylayers.gis.gisutil as gu
from pylayers.util.project import logger
import numpy.linalg as la
import pdb
import time
from numba import jit

def PL0(fGHz,GtdB=0,GrdB=0,R=1):
    """  Path Loss at frequency fGHZ @ R

    Parameters
    ----------

    fGHz: float
          frequency GHz
    GtdB: float
          transmitting antenna gain dB (default 0 dB)
    GrdB: float
          receiving antenna gain dB (default 0 dB)
    R   : float
        distance in m


    Returns
    -------

    PL0 : float
          path @ R

    Notes
    -----

    .. math:: PL_0 = -20 log_{10}(\\frac{\\lambda}{4\\pi}) - GtdB -GrdB

    Examples
    --------

    >>> fGHz  = 2.4
    >>> PL = PL0(fGHz)
    >>> assert (PL<41)&(PL>40),"something wrong"

    """

    if not isinstance(fGHz,np.ndarray):
        fGHz=np.array([fGHz])

    ld  = 0.3/fGHz
    PL0 = -20*np.log10(ld/(4.0*np.pi*R))-GtdB-GrdB

    return PL0

def Dgrid_points(points,Px):
    """ distance point to grid

    Parameters
    ----------

    points : np.array
             grid Np x 2 array

    Px     : np.array
             point 2 x 1  array

    Returns
    -------

    D: Euclidian distance matrix

    """

    Dx = points[:,0] - Px[0]
    Dy = points[:,1] - Px[1]

    D = np.sqrt( Dx*Dx + Dy*Dy )

    return(D)


def FMetisShad2(fGHz,r,D,sign=1):
    """ F Metis shadowing function

    Parameters
    ----------

    fGHz : np.array(Nf)
        frequency GHz
    r : np.array(Nseg,)
        distance between Tx and Rx
    D : np.array(Nseg,Nscreen)
        indirect distance between Tx and Rx (screen effect)
    sign : np.array(Nseg,Nscreen)
        == 1  : Shadowing NLOS situation
        ==-1  : No shadowing LOS situation


    Returns
    -------

    F : np.array(Nseg,Nscreen,Nf)

    Notes
    -----

    Provides an implementation of formula (6.6) in D1.4 of METIS project


    See Also
    --------

    LossMetisShadowing

    The arctan is approximating the enveloppe of the Fresnel integral.

    """
    lamda = 0.3/fGHz[None,None,:]
    F = np.arctan(sign[:,:,None]*np.pi/2.*(np.sqrt((np.pi/lamda)*(D[:,:,None]-r[:,None,None])))) / np.pi
    return(F)

def FMetisShad(fGHz,r,D,sign=1):
    """ F Metis shadowing function

    Parameters
    ----------

    fGHz : float
        frequency GHz
    r : float
        distance between Tx and Rx
    D : float
        indirect distance between Tx and Rx (screen effect)
    sign : int
        == 1  : Shadowing NLOS situation
        ==-1  : No shadowing LOS situation   


    Notes
    -----

    Provides an implementation of formula (6.6) in D1.4 of METIS project


    See Also
    --------

    LossMetisShadowing


    """
    lamda = 0.3/fGHz
    F = np.arctan(sign*np.pi/2.*(np.sqrt((np.pi/lamda)*(D-r)))) / np.pi
    return(F)

def LossMetisShadowing(fGHz,tx,rx,pg,uw,uh,w,h):
    """ Calculate the Loss from 

    Parameters
    ----------
    fGHz : float

    tx  : np.array (,3) of floats
        transmiter coordinates
    rx  : np.array (,3) of floats
        receiver coordinates
    pg  : np.array (,3) of floats
        center of gravity of the screen 
    uw  : np.array (,3) of floats 
        unitary vector along width dimension
    uh  : np.array (,3) of floats 
        unitary vector along height dimension
    w   : float 
        width in meters
    h   : float 
        height in meters 

    Returns
    -------

    Lsh : float 
        Loss in dB to add to the FS path Loss


    Notes
    -----

    This function provides an implementation of formula 6.5 of D1.4 deliverable of METIS project

    [Metis D1.4](Ahttps://www.metis2020.com/wp-content/uploads/METIS_D1.4_v3.pdf)

    # geometry parametric issue : find M in [tx-rx] defined as M = alpha*rx + (1-alpha)tx where alpha in [0-1].
    # if alpha = 0 then M = tx ; if alpha = 1 then M = rx.
    # Besides, M is defined as M = pg + beta*uw + gamma*uh then  alpha*rx + (1-alpha)tx = pg + beta*uw + gamma*uh
    # [rx-tx , -uw, -uh]*[alpha,beta,gamma].T = pg - tx <==> Ax = b solved by la.solve ; x[0]=alpha, x[1]=beta and

    TODO To be vectorized

    """



    rxtx = rx - tx # LOS distance

    # x[2]=gamma.
    A = np.vstack((rxtx,-uw,-uh)).T
    b = pg - tx
    x = la.solve(A,b)

    # condition of shadowing
    condseg = ((x[0]>1) or (x[0]<0))
    condw = ((x[1]>w/2.) or (x[1]<-w/2.))
    condh = ((x[2]>h/2.) or (x[2]<-h/2.))

    visi = condseg or condw or condh
    if visi:
        shad = -1
    else:
        shad = 1

    r = np.dot(rxtx,rxtx)**0.5
    w1 = pg + uw*w/2.
    w2 = pg - uw*w/2.
    h1 = pg + uh*h/2.
    h2 = pg - uh*h/2.


    Dtw1 = np.dot(tx-w1,tx-w1)**0.5
    Drw1 = np.dot(rx-w1,rx-w1)**0.5
    Dtw2 = np.dot(tx-w2,tx-w2)**0.5
    Drw2 = np.dot(rx-w2,rx-w2)**0.5
    Dth1 = np.dot(tx-h1,tx-h1)**0.5
    Drh1 = np.dot(rx-h1,rx-h1)**0.5
    Dth2 = np.dot(tx-h2,tx-h2)**0.5
    Drh2 = np.dot(rx-h2,rx-h2)**0.5

    D1w = Dtw1+Drw1
    D1h = Dth1+Drh1
    D2w = Dtw2+Drw2
    D2h = Dth2+Drh2


    if shad == 1:
        signw1 = 1
        signw2 = 1
        signh1 = 1
        signh2 = 1
    else:
        if condw:
            if D1w>D2w:
                signw1=1
                signw2=-1
            else:
                signw1=-1
                signw2=1
        else:
            signw1 = 1
            signw2 = 1
        
        if condh:
            if D1h>D2h:
                signh1=1
                signh2=-1
            else:
                signh1=-1
                signh2=1
        else:
            
            signh1 = 1
            signh2 = 1

    
            
    Fw1 = FMetisShad(fGHz,r,D1w,sign=signw1)
    Fh1 = FMetisShad(fGHz,r,D1h,sign=signh1)
    Fw2 = FMetisShad(fGHz,r,D2w,sign=signw2)
    Fh2 = FMetisShad(fGHz,r,D2h,sign=signh2)
    tmp = (Fh1+Fh2)*(Fw1+Fw2)
    Lsh = -20*np.log10(1-tmp)

    #return(Lsh,shad,tmp,Fw1,Fh1,Fw2,Fh2,condh,condw)
    return(Lsh)

def LossMetisShadowing2(fGHz,tx,rx,pg,uw,uh,w,h):
    """ Calculate the Loss from 

    Parameters
    ----------

    fGHz : np.array(,Nf)
        
    tx  : np.array (3,Nseg) of floats  
        transmiter coordinates 
    rx  : np.array (3,Nseg) of floats  
        receiver coordinates 
    pg  : np.array (3,Nscreen) of floats 
        center of gravity of the screen 
    uw  : np.array (3,Nscreen) of floats 
        unitary vector along width dimension
    uh  : np.array (3,Nscreen) of floats 
        unitary vector along height dimension
    w   : np.array (,Nscreen)
        width in meters
    h   : np.array (,Nscreen)
        height in meters 

    Returns
    -------

    Lsh : np.array (Nseg,Nscreen,Nf)
        Loss in dB to add to the FS path Loss


    Notes
    -----

    This function provides an implementation of formula 6.5 of D1.4 deliverable of METIS project

    [Metis D1.4](Ahttps://www.metis2020.com/wp-content/uploads/METIS_D1.4_v3.pdf)

    # geometry parametric issue : find M in [tx-rx] defined as M = alpha*rx + (1-alpha)tx where alpha in [0-1].
    # if alpha = 0 then M = tx ; if alpha = 1 then M = rx.
    # Besides, M is defined as M = pg + beta*uw + gamma*uh then  alpha*rx + (1-alpha)tx = pg + beta*uw + gamma*uh
    # [rx-tx , -uw, -uh]*[alpha,beta,gamma].T = pg - tx <==> Ax = b solved by la.solve ; x[0]=alpha, x[1]=beta and

    
    """

    Nseg = tx.shape[1]
    Nscreen = uw.shape[1]

    rxtx = rx - tx # (3,Nseg) LOS distance
   

    # A : (Nseg,Nscreen,3,3)
    # b : (Nseg,Nscreen,3)
    # rxtx.T (Nseg,3)
    # uw.T (Nscreen, 3)
    # uh.T (Nscreen, 3)

    U = rxtx.T[:,None,:,None]
    W = uw.T[None,:,:,None]
    H = uh.T[None,:,:,None]

    We = W + np.zeros(U.shape)
    He = H + np.zeros(U.shape)
    Ue = U + np.zeros(He.shape)

    A = np.concatenate((Ue,-We,-He),axis=3)
    #A = np.vstack((rxtx,-uw,-uh)).T 

    # pg.T Nscreen, 3
    # tx.T Nseg,3
    b = pg.T[None,:,:]-tx.T[:,None,:] 
    #b = pg - tx
    x = la.solve(A,b)
    
    
    # condition of shadowing
    condseg = ((x[:,:,0]>1) + (x[:,:,0]<0)) 
    condw = ((x[:,:,1]>w[None,:]/2.) + (x[:,:,1]<-w[None,:]/2.)) 
    condh = ((x[:,:,2]>h[None,:]/2.) + (x[:,:,2]<-h[None,:]/2.)) 
    
    visi = (condseg + condw + condh)%2

    
    # if visi:
    #     shad = -1
    # else:
    #     shad = 1
    #shad = - visi    
    
    r = np.sum(rxtx*rxtx,axis=0)**0.5


    w1 = pg + uw*w[None,:]/2.
    w2 = pg - uw*w[None,:]/2.
    h1 = pg + uh*h[None,:]/2.
    h2 = pg - uh*h[None,:]/2.


    
    Dtw1 = np.sum((tx[...,None]-w1[:,None,:])*(tx[...,None]-w1[:,None,:]),axis=0)**0.5
    Drw1 = np.sum((rx[...,None]-w1[:,None,:])*(rx[...,None]-w1[:,None,:]),axis=0)**0.5
    Dtw2 = np.sum((tx[...,None]-w2[:,None,:])*(tx[...,None]-w2[:,None,:]),axis=0)**0.5
    Drw2 = np.sum((rx[...,None]-w2[:,None,:])*(rx[...,None]-w2[:,None,:]),axis=0)**0.5

    Dth1 = np.sum((tx[...,None]-h1[:,None,:])*(tx[...,None]-h1[:,None,:]),axis=0)**0.5
    Drh1 = np.sum((rx[...,None]-h1[:,None,:])*(rx[...,None]-h1[:,None,:]),axis=0)**0.5
    Dth2 = np.sum((tx[...,None]-h2[:,None,:])*(tx[...,None]-h2[:,None,:]),axis=0)**0.5
    Drh2 = np.sum((rx[...,None]-h2[:,None,:])*(rx[...,None]-h2[:,None,:]),axis=0)**0.5

    # Drw1 = np.dot(rx-w1,rx-w1)**0.5
    # Dtw2 = np.dot(tx-w2,tx-w2)**0.5
    # Drw2 = np.dot(rx-w2,rx-w2)**0.5
    # Dth1 = np.dot(tx-h1,tx-h1)**0.5
    # Drh1 = np.dot(rx-h1,rx-h1)**0.5
    # Dth2 = np.dot(tx-h2,tx-h2)**0.5
    # Drh2 = np.dot(rx-h2,rx-h2)**0.5
    
    
    D1w = Dtw1+Drw1
    D1h = Dth1+Drh1
    D2w = Dtw2+Drw2
    D2h = Dth2+Drh2
    
    signw1 = np.ones((Nseg,Nscreen))
    signw2 = np.ones((Nseg,Nscreen))
    signh1 = np.ones((Nseg,Nscreen))
    signh2 = np.ones((Nseg,Nscreen))


    condw1 = (visi*condw*(D1w<=D2w)).astype(bool)
    condw2 = (visi*condw*(D1w>D2w)).astype(bool)
    signw1[condw1]=-1
    signw2[condw2]=-1
    condh1 = (visi*condh*(D1h<=D2h)).astype(bool)
    condh2 = (visi*condh*(D1h>D2h)).astype(bool)
    signh1[condh1]=-1
    signh2[condh2]=-1


    Fw1 = FMetisShad2(fGHz,r,D1w,sign=signw1)
    Fh1 = FMetisShad2(fGHz,r,D1h,sign=signh1)
    Fw2 = FMetisShad2(fGHz,r,D2w,sign=signw2)
    Fh2 = FMetisShad2(fGHz,r,D2h,sign=signh2)
    tmp = (Fh1+Fh2)*(Fw1+Fw2)
    Lsh = -20*np.log10(1-tmp)

    #return(Lsh,shad,tmp,Fw1,Fh1,Fw2,Fh2,condh,condw)
    return(Lsh)

def Dgrid_zone(zone,Px):
    """ Distance point to zone

    A zone is a quadrilateral zone.

    Parameters
    ----------

    zone : dictionnary
           xmin xmax Nx
           ymin ymax Ny

    Px : np.array
         point

    Returns
    -------

    D : np.array Nx x Ny
        Euclidian distance matrix

    Notes
    -----

    Build the distance matrix between Tx and points in the zone

    use broadcasting instead

    """

    rx = np.linspace(zone['xmin'],zone['xmax'],zone['Nx'])
    ry = np.linspace(zone['ymin'],zone['ymax'],zone['Ny'])

    R_x = np.outer(np.ones(len(ry)),rx)
    R_y = np.outer(ry,np.ones(len(rx)))

    Dx = R_x - Px[0]
    Dy = R_y - Px[1]
    D = np.sqrt(Dx*Dx+Dy*Dy)
    return (D)

def OneSlopeMdl(D,n,fGHz):
    """ one slope model

    Parameters
    ----------

    D   : np.array
          distance array
    n   : float
          path loss exponent
    fGHz : np.array
           frequency  GHz

    Returns
    -------

    PL : np.array
          path loss as a function of distance

    """

    PL = PL0(fGHz) + 10*n*np.log10(D)

    return(PL)

def cost231(pBS,pMS,hroof,phir,wr,fMHz,wb=20,dB=True,city='medium'):
    """ Walfish Ikegami model (COST 231)

    Parameters
    ----------

    pBS   : np.array (3xNlink)
    pMS   : np.array (3xNlink)
    hroof : np.array (1xNlink)
    phir  : np.array (1xNlink)
        degrees
    wr    : np.array (1xNlink)
    fMHz  : np.array (1xNf)
    wb    : float
        average building separation
    dB    : boolean

    Returns
    -------

    PathLoss  (Nlink,Nf)

    References
    ----------

    http://morse.colorado.edu/~tlen5510/text/classwebch3.html

    Examples
    --------

    .. plot::
        :include-source:

        >>> from pylayers.antprop.loss import *
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> # Number of links and BS and MS heights
        >>> Nlink = 100
        >>> hBS = 300
        >>> hMS = 1.5
        >>> # hroof and phir are drawn uniformily at random 
        >>> hroof = 40*np.random.rand(Nlink)
        >>> wr = 10*np.ones(Nlink)
        >>> phir = 90*np.random.rand(Nlink)
        >>> pMS = np.vstack((np.linspace(10,2500,Nlink),np.zeros(Nlink),hMS*np.ones(Nlink)))
        >>> pBS = np.vstack((np.zeros(Nlink),np.zeros(Nlink),hBS*np.ones(Nlink)))
        >>> # frequency range 
        >>> fMHz = np.linspace(700,1900,120)
        >>> pl = cost231(pBS,pMS,hroof,phir,wr,fMHz)
        >>> im = plt.imshow(pl,extent=(0,100,0.7,1.9))
        >>> cb = plt.colorbar()
        >>> cb.set_label('Loss (dB)')
        >>> plt.axis('tight')
        >>> plt.xlabel('Frequency (GHz)')
        >>> plt.ylabel('Link Number')
        >>> plt.title('100 WI Path Loss realizations ')
        >>> plt.show()

    """

    hBS = pBS[2,:][:,np.newaxis]
    hMS = pMS[2,:][:,np.newaxis]
    wr  = wr[:,np.newaxis]
    hroof = hroof[:,np.newaxis]
    phir = phir[:,np.newaxis]
    fMHz = fMHz[np.newaxis,:]
    dm  = np.sqrt(np.sum((pBS-pMS)*(pBS-pMS),axis=0))[:,np.newaxis]
    dkm = dm/1000.
    Nlink = len(dm)

    pl0 = 32.4 + 20*np.log10(dkm) + 20*np.log10(fMHz)

    delta_base = hBS-hroof

    u035  = np.where((phir>=0) & (phir<35))
    u3555 = np.where((phir>=35) & (phir<55))
    u5590 = np.where((phir>=55) & (phir<90))

    plori = np.zeros(Nlink)[:,np.newaxis]
    # path loss due to orientation w.r.t road
    plori[u035] = -10+0.354*phir[u035]
    plori[u3555] = 2.5+0.075*phir[u3555]
    plori[u5590] = 4.0-0.114*(phir[u5590]-55)

    # rooftop to street
    plrts = -16.9-10*np.log10(wr)+10*np.log10(fMHz)+20*np.log10(hroof-hMS)+plori
    uroofsupBS  = np.where(hBS>hroof)
    uroofinfBS  = np.where(hBS<=hroof)
    udistsup500 = np.where((hBS<=hroof)&(dkm>0.5))
    udistinf500 = np.where((hBS<=hroof)&(dkm<0.5))

    plbsh = np.zeros((Nlink,1))
    plbsh[uroofsupBS] = -18*np.log10(1+delta_base[uroofsupBS])

    ka  = 54*np.ones((Nlink,1))
    ka[udistsup500] = ka[udistsup500]-0.8*delta_base[udistsup500]
    ka[udistinf500] = ka[udistinf500]-0.8*delta_base[udistinf500]*dkm[udistinf500]/0.5

    kd  = 18*np.ones((Nlink,1))
    kd[uroofinfBS] = kd[uroofinfBS]-15*delta_base[uroofinfBS]/hroof[uroofinfBS]

    if city=='medium':
        kf = -4+0.7*(fMHz/925.-1)
    else:
        kf = -4+1.5*(fMHz/925.-1)

    plmsd = plbsh+ka+kd*np.log10(dkm)+kf*np.log10(fMHz)-9*np.log10(wb)

    pl = pl0
    padd = plmsd + plrts
    ulosspos = np.where(padd>0)[0]
    pl[ulosspos]=pl[ulosspos]+padd[ulosspos]

    if not dB:
        pl = 10**(-pl/20.)

    return(pl)

def cost259(pMS,pBS,fMHz):
    """ cost259 model

    Parameters
    ----------

    pMS : np.array (position of Mobile Station)
    pBS : np.array (position of Base station)
    fMHz : float

    """

    dm  = np.sqrt((pBS-pMS)*(pBS-pMS))
    lmbd = 300/fMHz
    pl = 10*2.6*np.log10(dm)+20*log10(4*np.pi/lmbd)

    if not dB:
        pl = 10**(-pl/20.);
    return(pl)

def hata(pMS,pBS,fGHz,hMS,hBS,typ):
    """ Hata Path loss model

    Parameters
    ----------

    pMS : np.array
        Mobile position (meters)
    pBS : np.array
        Base station position (meters)
    fGHz : np.array
    hMS : height mobile station (m)
    hBS : height base station (m)
    typ : 'small' or 'big'

    Returns
    -------

    L : Attenuation (dB)

    Examples
    --------

        >>> d = np.linspace(100,5000,120)
        >>> hBS = 30
        >>> hMS = 1.5
        >>> fGHz = 0.9
        >>> pMS = np.array([d,0,hMS])
        >>> pBS = np.array([0,0,hBS])
        >>> L = hata(pMS,pBS,fGHz,hMS,hBS,'small')

    Notes
    -----

    This model is valid until 1.5GHz, for higher frequency see

    COST231-Hata model

    References
    ----------

    OKUMURA (Y.), OHMORI (E.), KAWANO (T.)
    et FUKUA (K.).  Field strength and its varia-
    bility in UHF and VHF land-mobile radio ser-
    vice. Rev. Elec. Commun. Lab., vol. 16, n o 9,
    1968.

    HATA (M.).  Empirical formula for propaga-
    tion loss in land mobile radio services. IEEE
    Trans. Veh. Technol., vol. 29, pp. 317-325,
    Aug. 1980

    """
    dm  = np.sqrt(np.sum((pBS-pMS)*(pBS-pMS),axis=0))
    dkm = dm/1000
    fMHz = fGHz*1000
    if (typ=='small'):
       CH = (1.1*np.log10(fMHz)-0.7)*hMS-(1.56*np.log10(fMHz)-0.8)
    if (typ=='big'):
        if fGHz<0.2:
            CH = 8.29*(np.log10(1.54*hMS)**2)-1.1
        else:# valid until 1.5GHz
            CH = 3.2*(np.log10(11.75*hMS)**2)-4.97
    L = 69.55+26.16*np.log10(fMHz)-13.82*np.log10(hBS)+(44.9-6.55*np.log10(hBS))*np.log10(dkm)-CH

    return(L)

def cost2100(pMS,pBS,fGHz,nfloor=1,dB=True):
    """ cost 2100 model

    Parameters
    ----------

    pMS :
    pBS :
    fGHz : float
    nfloor : int
    dB : boolean

    """
    # distance (meters)
    dm  = np.sqrt((pBS-pMS)*(pBS-pMS))
    pl0 = 32.4+20*log10(dm)+20*np.log10(fGHz)
    pld = nfloor*30
    pl  = pl0+pld
    if not dB:
        pl = 10**(-pl/20.)
    return(pl)

def PL(fGHz,pts,p,n=2.0,dB=True,d0=1):
    """ calculate Free Space Path Loss

    Parameters
    ----------

    fGHz   : float
             frequency (GHz)
    pts    : np.array (2xNp)
             points
    p      : np.array (2x1) or (2xNp)
    n      : float
            path loss exponent (default = 2)

    dB : : boolean
        return result in dB


    Returns
    -------

    PL : np.array
         path loss w.r.t distance and frequency

    """
    shp = np.shape(p)
    # assert(shp[0]==2)

    D = np.sqrt(np.sum((pts-p)**2,axis=0))
    # f x grid x ap
    #PL = np.array([PL0(fGHz)])[:,np.newaxis] + 10*n*np.log10(D)[np.newaxis,:]

    PL = PL0(fGHz,d0)[:,np.newaxis] + 10*n*np.log10(D/d0)[np.newaxis,:]

    if not dB:
        PL=10**(-PL/10)

    return(PL)

def Lconcrete(fGHz):
    """ 3GPP above 6GHz concrete loss

    From Table 7.4.3.1 of 3GPP TR38.900
    Study on channel model for frequency spectrum above 6 GHz

    """

    return(5+4*fGHz)

def Lglass(fGHz):
    """ 3GPP above 6GHz glass loss

    From Table 7.4.3.1 of 3GPP TR38.900
    Study on channel model for frequency spectrum above 6 GHz

    """
    return(2 + 0.2* fGHz)

def LIRRglass(fGHz):
    """ 3GPP above 6GHz IRR (Infra Red Reflecting) glass loss

    From Table 7.4.3.1 of 3GPP TR38.900
    Study on channel model for frequency spectrum above 6 GHz

    """
    return(23 + 0.3* fGHz)

def Lwood(fGHz):
    """ 3GPP above 6GHz wood loss

    From Table 7.4.3.1 of 3GPP TR38.900
    Study on channel model for frequency spectrum above 6 GHz

    """
    return(4.85 + 0.12* fGHz)




def LossPenetration(fGHz, alpha = 0.3, typ='low'):
    if typ=='low':
        PLTW = 5 - 10*np.log10(alpha*10**(-Lglass(fGHz)/10)+(1-alpha)*10**(-Lconcrete(fGHz)))

    if typ=='high':
        PLTW = 5 - 10*np.log10((1-alpha)*10**(-LIRRglass(fGHz)/10)+alpha*10**(-Lconcrete(fGHz)))

    return(PLTW)


def Losst(L,fGHz,p1,p2,dB=True,bceilfloor=False):
    """  calculate Losses between links p1-p2

    Parameters
    ----------

    L   : Layout object
    fGHz : np.array
           frequency GHz
    p1 : source points
        (3 x Np1) array or (3,) array
    p2 : observation points
        (3 x Np2) array or (3,) array
    dB : boolean
    bceilfloor : boolean

    Examples
    --------

    .. plot::
        :include-source:

        >>> import matplotlib.pyplot as plt
        >>> from pylayers.measures.mesuwb import *
        >>> from pylayers.antprop.loss import *
        >>> S = Simul()
        >>> S.layout('WHERE1.lay')
        >>> fGHz = 4
        >>> Tx,Rx = ptw1()
        >>> Lwo,Lwp,Edo,Edp = Losst(S.L,fGHz,Tx.T,Rx[1,0:2],dB=True)
        >>> fig=plt.figure(figsize=(20,10))
        >>> fig,ax = S.L.showGs(fig=fig)
        >>> tit = plt.title('test Losst')
        >>> sc2 = ax.scatter(Rx[1,0],Rx[1,1],s=20,marker='x',c='k')
        >>> sc1 = ax.scatter(Tx[:,0],Tx[:,1],s=20,c=Lwo,linewidth=0)
        >>> cb = plt.colorbar(sc1)
        >>> cb.set_label('dB')
        >>> plt.show()

    See Also
    --------

    pylayers.antprop.coverage
    pylayers.slab.Interface.losst

    """

    if (type(fGHz)==float) | (type(fGHz)==int):
        fGHz=np.array([fGHz],dtype=float)

    sh1 = np.shape(p1)
    sh2 = np.shape(p2)

    if (len(sh1)>1) & (len(sh2)>1):
        Nlink = max(sh1[1],sh2[1])
    if (len(sh1)>1) & (len(sh2)<2):
        Nlink = sh1[1]
    if (len(sh1)<2) & (len(sh2)>1):
        Nlink = sh2[1]
    if (len(sh1)<2) & (len(sh2)<2):
        Nlink = 1

    # determine incidence angles on segment crossing p1-p2 segment
    #data = L.angleonlink(p1,p2)
    logger.debug('losst before angleonlink3')
    data = L.angleonlink3(p1,p2)

    # as many slabs as segments and subsegments
    us = data['s']
    slabs = np.array([ L.Gs.node[x]['name'] for x in us ])


    #slabs = L.sla[us]
    check = np.where(slabs=='')

    #
    # As segment numbering is not necessarily contiguous
    # there exist void string '' in slabs
    cslab = list(np.unique(slabs))
    if '' in cslab:
        cslab.remove('')
    if 'AIR' in cslab:
        cslab.remove('AIR')
    if '_AIR' in cslab:
        cslab.remove('_AIR')


    LossWallo = np.zeros((len(fGHz),Nlink))
    LossWallp = np.zeros((len(fGHz),Nlink))
    EdWallo = np.zeros((len(fGHz),Nlink))
    EdWallp = np.zeros((len(fGHz),Nlink))

    for slname in cslab:
        # u index of slabs of name slname
        # data['a'][u] angle
        # data['s'][u] segment number including subsegment
        u = np.nonzero(np.array(slabs)==slname)[0]
        #
        # calculate Loss for slab slname
        #
        lko,lkp  = L.sl[slname].losst(fGHz,data['a'][u])
        #
        # calculate Excess delay for slab slname
        #
        do , dp  = L.sl[slname].excess_grdelay(theta=data['a'][u])
        # data['i'][u] links number
        indexu = data['i'][u]
        # reduce to involved links
        involved_links, indices = np.unique(indexu,return_index=True)
        indicep = np.hstack((indices[1:],np.array([len(indexu)])))
        # range on involved links
        irange = np.arange(len(involved_links))
        #
        # sum contribution of slab of a same link
        #
        Wallo = np.array([ np.sum(lko[:,indices[x]:indicep[x]],axis=1) for x in irange ] ).T
        Wallp = np.array([ np.sum(lkp[:,indices[x]:indicep[x]],axis=1) for x in irange ] ).T

        Edo = np.array([np.sum(do[indices[x]:indicep[x]]) for x in irange]).T
        Edp = np.array([np.sum(dp[indices[x]:indicep[x]]) for x in irange]).T

        LossWallo[:,involved_links] = LossWallo[:,involved_links] + Wallo
        LossWallp[:,involved_links] = LossWallp[:,involved_links] + Wallp

        EdWallo[:,involved_links] = EdWallo[:,involved_links] + Edo
        EdWallp[:,involved_links] = EdWallp[:,involved_links] + Edp

    if bceilfloor:
    # Managing Ceil / Floor transmission
    # At that point there is only a single type of ceil and floor
    # it should be defined ideally as a specific entity
        #
        # TODO fix shape error p2 is not always (3 x N)
        #
        if (p1[2,:]> L.zceil).any() or (p2[2,:]> L.zceil).any():

            # WARNING : this test should be done individually
            if (p1[2]>p2[2]).all():
                v0 = p1
                v1 = p2
            else:
                v0 = p2
                v1 = p1

            uu = v0 - v1
            # 1 x N
            nu = np.sqrt(np.sum(uu * uu, axis=0))
            # 3 x N
            un = uu / nu[np.newaxis, :]
            dotp = np.einsum('ij,i->j',un,np.array([0,0,1]))
            alphas = np.arccos(dotp)

            #
            # calculate Loss for slab CEIL
            #
            lkco,lkcp  = L.sl['CEIL'].losst(fGHz,alphas)
            #
            # calculate Excess delay for slab CEIL
            #
            dco , dcp  = L.sl['CEIL'].excess_grdelay(theta=alphas)



            LossWallo = LossWallo + lkco
            LossWallp = LossWallp + lkcp

            EdWallo = EdWallo + dco
            EdWallp = EdWallp + dcp


        # check crossing floor
        if (p1[2,:]< L.zfloor).any() or (p2[2,:]< L.zfloor).any():


            # WARNING : this test should be done individually
            if (p1[2]>p2[2]).all():
                v0 = p1
                v1 = p2
            else:
                v0 = p2
                v1 = p1

            uu = v0 - v1
            # 1 x N
            nu = np.sqrt(np.sum(uu * uu, axis=0))
            # 3 x N
            un = uu / nu[np.newaxis, :]
            dotp = np.einsum('ij,i->j',un,np.array([0,0,1]))
            alphas = np.arccos(dotp)

            #
            # calculate Loss for slab CEIL
            #
            lkfo,lkfp  = L.sl['FLOOR'].losst(fGHz,alphas)
            #
            # calculate Excess delay for slab CEIL
            #
            dfo , dfp  = L.sl['FLOOR'].excess_grdelay(theta=alphas)

            LossWallo = LossWallo + lkfo
            LossWallp = LossWallp + lkfp

            EdWallo = EdWallo + dfo
            EdWallp = EdWallp + dfp

    if not dB:
        LossWallo = 10**(-LossWallo/10)
        LossWallp = 10**(-LossWallp/10)

    return(LossWallo,LossWallp,EdWallo,EdWallp)

def gaspl(d,fGHz,T,PhPa,wvden):
    """ attenuation due to atmospheric gases

    Parameters
    ----------
    d : np.array
        range (meters)
    fGHz : np.array
        frequency (GHz)
    T : float
        Temprature in degree Celcius
    PhPa : float
        Pressure in hPa
    wvden : float
        Water vapor density (g/m**3)

    Examples
    --------

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> T = 15
    >>> PhPa = 1013
    >>> wvden = 7.5
    >>> d = 1000
    >>> fGHz = np.linspace(1,1000,100)
    >>> L = gaspl(d,fGHz,T,PhPa,wvden)
    >>> plt.plot(fGHz,L)

    Notes
    -----

    This function implements the recommandation UIT-P676-10

    """

    affO2=np.array([
    [50.474214,0.975,9.651,6.690,0.0,2.566,6.850],
    [50.987745,2.529,8.653,7.170,0.0,2.246,6.800],
    [51.503360,6.193,7.709,7.640,0.0,1.947,6.729],
    [52.021429,14.320,6.819,8.110,0.0,1.667,6.640],
    [52.542418,31.240,5.983,8.580,0.0,1.388,6.526],
    [53.066934,64.290,5.201,9.060,0.0,1.349,6.206],
    [53.595775,124.600,4.474,9.550,0.0,2.227,5.085],
    [54.130025,227.300,3.800,9.960,0.0,3.170,3.750],
    [54.671180,389.700,3.182,10.370,0.0,3.558,2.654],
    [55.221384,627.100,2.618,10.890,0.0,2.560,2.952],
    [55.783815,945.300,2.109,11.340,0.0,-1.172,6.135],
    [56.264774,543.400,0.014,17.030,0.0,3.525,-0.978],
    [56.363399,1331.800,1.654,11.890,0.0,-2.378,6.547],
    [56.968211,1746.600,1.255,12.230,0.0,-3.545,6.451],
    [57.612486,2120.100,0.910,12.620,0.0,-5.416,6.056],
    [58.323877,2363.700,0.621,12.950,0.0,-1.932,0.436],
    [58.446588,1442.100,0.083,14.910,0.0,6.768,-1.273],
    [59.164204,2379.900,0.387,13.530,0.0,-6.561,2.309],
    [59.590983,2090.700,0.207,14.080,0.0,6.957,-0.776],
    [60.306056,2103.400,0.207,14.150,0.0,-6.395,0.699],
    [60.434778,2438.000,0.386,13.390,0.0,6.342,-2.825],
    [61.150562,2479.500,0.621,12.920,0.0,1.014,-0.584],
    [61.800158,2275.900,0.910,12.630,0.0,5.014,-6.619],
    [62.411220,1915.400,1.255,12.170,0.0,3.029,-6.759],
    [62.486253,1503.000,0.083,15.130,0.0,-4.499,0.844],
    [62.997984,1490.200,1.654,11.740,0.0,1.856,-6.675],
    [63.568526,1078.000,2.108,11.340,0.0,0.658,-6.139],
    [64.127775,728.700,2.617,10.880,0.0,-3.036,-2.895],
    [64.678910,461.300,3.181,10.380,0.0,-3.968,-2.590],
    [65.224078,274.000,3.800,9.960,0.0,-3.528,-3.680],
    [65.764779,153.000,4.473,9.550,0.0,-2.548,-5.002],
    [66.302096,80.400,5.200,9.060,0.0,-1.660,-6.091],
    [66.836834,39.800,5.982,8.580,0.0,-1.680,-6.393],
    [67.369601,18.560,6.818,8.110,0.0,-1.956,-6.475],
    [67.900868,8.172,7.708,7.640,0.0,-2.216,-6.545],
    [68.431006,3.397,8.652,7.170,0.0,-2.492,-6.600],
    [68.960312,1.334,9.650,6.690,0.0,-2.773,-6.650],
    [118.750334,940.300,0.010,16.640,0.0,-0.439,0.079],
    [368.498246,67.400,0.048,16.400,0.0,0.000,0.000],
    [424.763020,637.700,0.044,16.400,0.0,0.000,0.000],
    [487.249273,237.400,0.049,16.000,0.0,0.000,0.000],
    [715.392902,98.100,0.145,16.000,0.0,0.000,0.000],
    [773.839490,572.300,0.141,16.200,0.0,0.000,0.000],
    [834.145546,183.100,0.145,14.700,0.0,0.000,0.000]])
    ## spectroscopic data for H20
    ## f0 b1 b2 b3 b4 b5 b6
    affH2O=np.array([
    [22.235080,0.1130,2.143,28.11,0.69,4.800,1.00],
    [67.803960,0.0012,8.735,28.58,0.69,4.930,0.82],
    [119.995940,0.0008,8.356,29.48,0.70,4.780,0.79],
    [183.310091,2.4200,0.668,30.50,0.64,5.300,0.85,],
    [321.225644,0.0483,6.181,23.03,0.67,4.690,0.54],
    [325.152919,1.4990,1.540,27.83,0.68,4.850,0.74],
    [336.222601,0.0011,9.829,26.93,0.69,4.740,0.61],
    [380.197372,11.5200,1.048,28.73,0.54,5.380,0.89,],
    [390.134508,0.0046,7.350,21.52,0.63,4.810,0.55],
    [437.346667,0.0650,5.050,18.45,0.60,4.230,0.48,],
    [439.150812,0.9218,3.596,21.00,0.63,4.290,0.52,],
    [443.018295,0.1976,5.050,18.60,0.60,4.230,0.50],
    [448.001075,10.3200,1.405,26.32,0.66,4.840,0.67],
    [470.888947,0.3297,3.599,21.52,0.66,4.570,0.65,],
    [474.689127,1.2620,2.381,23.55,0.65,4.650,0.64],
    [488.491133,0.2520,2.853,26.02,0.69,5.040,0.72],
    [503.568532,0.0390,6.733,16.12,0.61,3.980,0.43],
    [504.482692,0.0130,6.733,16.12,0.61,4.010,0.45],
    [547.676440,9.7010,0.114,26.00,0.70,4.500,1.00],
    [552.020960,14.7700,0.114,26.00,0.70,4.500,1.00],
    [556.936002,487.4000,0.159,32.10,0.69,4.110,1.00],
    [620.700807,5.0120,2.200,24.38,0.71,4.680,0.68],
    [645.866155,0.0713,8.580,18.00,0.60,4.000,0.50],
    [658.005280,0.3022,7.820,32.10,0.69,4.140,1.00],
    [752.033227,239.6000,0.396,30.60,0.68,4.090,0.84],
    [841.053973,0.0140,8.180,15.90,0.33,5.760,0.45],
    [859.962313,0.1472,7.989,30.60,0.68,4.090,0.84],
    [899.306675,0.0605,7.917,29.85,0.68,4.530,0.90],
    [902.616173,0.0426,8.432,28.65,0.70,5.100,0.95],
    [906.207325,0.1876,5.111,24.08,0.70,4.700,0.53],
    [916.171582,8.3400,1.442,26.70,0.70,4.780,0.78],
    [923.118427,0.0869,10.220,29.00,0.70,5.000,0.80],
    [970.315022,8.9720,1.920,25.50,0.64,4.940,0.67],
    [987.926764,132.1000,0.258,29.85,0.68,4.550,0.90],
    [1780.000000,22,300.0000,0.952,176.20,0.50,30.500]])
    
    dkm = d/1000.
    TK = T + 273.15
    theta = 300./TK #3

    fO2 = affO2[:,0]
    a1 = affO2[:,1]
    a2 = affO2[:,2]
    a3 = affO2[:,3]
    a4 = affO2[:,4]
    a5 = affO2[:,5]
    a6 = affO2[:,6]
    fH2O = affH2O[:,0]
    b1 = affH2O[:,1]
    b2 = affH2O[:,2]
    b3 = affH2O[:,3]
    b4 = affH2O[:,4]
    b5 = affH2O[:,5]
    b6 = affH2O[:,6]
    e = wvden*TK/216.7  # 4

    SO2  = a1*1e-7*PhPa*(theta**3)*np.exp(a2*(1-theta))  # 3
    DO2  = a3*1e-4*(PhPa*(theta**(0.8-a4))+1.1*e*theta)  # 6a

    SH2O = b1*1e-1*e*(theta**(3.5))*np.exp(b2*(1-theta)) # 3
    DH2O = b3*1e-4*(PhPa*theta**b4+b5*e*theta**b6)  # 6a

    DO2_m = np.sqrt(DO2**2+2.25e-6) # 6b
    DH2O_m = 0.535*DH2O+np.sqrt(0.217*DH2O**2+(2.1316*1e-12*fH2O**2)/theta)
    deltaO2 = (a5+a6*theta)*1e-4*(PhPa+e)*theta**(0.8)
    #
    #  O2
    #
    uO2 = fO2[:,None]-fGHz[None,:]
    vO2 = fO2[:,None]+fGHz[None,:]
    n1O2 = DO2_m[:,None]-deltaO2[:,None]*uO2
    n2O2 = DO2_m[:,None]-deltaO2[:,None]*vO2
    d1O2 = uO2**2 + DO2_m[:,None]**2
    d2O2 = vO2**2 + DO2_m[:,None]**2
    FO2 = (fGHz[None,:]/fO2[:,None])*(n1O2/d1O2+n2O2/d2O2)
    UO2 = SO2[:,None]*FO2
    #
    # H2O
    #
    uH2O = fH2O[:,None]-fGHz[None,:]
    vH2O = fH2O[:,None]+fGHz[None,:]
    nH2O = DH2O_m[:,None]
    d1H2O = uH2O**2 + DH2O_m[:,None]**2
    d2H2O = vH2O**2 + DH2O_m[:,None]**2
    FH2O = (fGHz[None,:]/fH2O[:,None])*(nH2O/d1H2O+nH2O/d2H2O)
    UH2O = SH2O[:,None]*FH2O

    # Nsec  (8)
    
    dD = 5.6e-4*(PhPa+e)*theta**(0.8)
    t1 = 6.14e-5/(dD*(1.+(fGHz/dD)**2))
    t2 = 1.4e-12*PhPa*(theta**(1.5))/(1+1.9e-5*fGHz**(1.5))
    Nsec = fGHz*PhPa*(theta**2)*(t1+t2)  # 9 
    
    ulow  = np.where(fGHz<118.750343)[0]
    uhigh = np.where(fGHz>=118.750343)[0] 
    UO2low = UO2[:,ulow]
    UO2high = UO2[:,uhigh] 
    SO2low  = np.sum(UO2low,axis=0)
    SO2high = np.sum(UO2high[38:,:],axis=0)
    sSO2 = np.hstack((SO2low,SO2high))
    Npp = sSO2 + np.sum(UH2O,axis=0)+Nsec
    Npp = np.sum(UO2,axis=0) + np.sum(UH2O,axis=0)+Nsec
    gamma = 0.1820*fGHz*Npp 
    LgasdB = gamma*dkm 

    return(LgasdB)


def Loss0(S,rx,ry,f,p):
    """ calculate Loss through Layers for theta=0 deg

    Parameters
    ----------

    S  : Simulation object
    rx : extremity of link 
    ry : extremity of link 
    fGHz : float
        frequency GHz
    p :

    """
    Nx  = len(rx)
    Ny  = len(ry)
    Lw  = np.zeros((Nx,Ny))
    print(shape(Lw))
    i   = 0
    for x in rx:
        j   = 0
        for y in ry:
            Loss = 0
            pxy = np.array([x,y])
            seglist,theta = L.angleonlinkold(p,pxy)
            for k in seglist:
                name = L.name[k]
                lk = L.sl[name].loss0(f)
                Loss  = Loss + lk[0]
            Lw[i,j] = Loss
            j = j+1
        i = i+1
    return(Lw)

def Loss_diff(u):
    """ calculate Path Loss of the diffraction
    """
    if u < -0.7:
        Ld = 0
    elif u > 1.5:
        Ld = 13 + 20*np.log10(u)
    else:
        Ld = 6.9 + 20*np.log10(np.sqrt((u-0.1)**2+1)+u-0.1)

    return(Ld)

def calnu(h,d1,d2,fGHz):
    r""" Calculate the diffraction Fresnel parameter

    Parameters
    ----------

    h  : signed height w.r.t LOS (meter)
    d1 : distance 1 (meter)
    d2 : distance 2 (meter)
    fGHz  : frequency GHz

    Notes
    -----

    .. math::   \nu = h \sqrt{\frac{2}{\lambda} \frac{d_1+d_2}{d_1 d_2}}

    """

    ld = 0.3/fGHz
    nu = h*np.sqrt(2*(d1+d2)/(ld*d1*d2))

    return(nu)


def route(X, Y, Z, Ha, Hb, fGHz, K, method='deygout'):
    """ diffraction loss along a route


    Parameters
    ----------

    X : np.array (Nphi,Nr)
        cartesian coordinate grid
    Y : np.array (Nphi,Nr)
        cartesian coordinate grid
    Z : np.array (Nphi,Nr)
        height (meters)

    Ha : float
    Hb : float
    fGHz : np.array (,Nf)
        frequency in GHz
    method : 'deygout' | 'bullington'

    Returns
    -------

    L : Losses (dB)


    """
    Nphi, Nr = Z.shape

    if (type(fGHz) == float):
        fGHz = np.array([fGHz])

    Nf = len(fGHz)
    L = np.zeros((Nphi, Nf))
    L0 = np.zeros(Nf)
    # loop over azimut
    for ip in range(Nphi):
        x = X[ip, :]
        y = Y[ip, :]
        z = Z[ip, :]
        d = np.sqrt((x-x[0])**2+(y-y[0])**2)
        # effect of refraction in equivalent earth curvature
        dh = d*(d[::-1])/(2*K*6375e3)
        z = z + dh
        LOS = 32.4 + 20*np.log10(fGHz) + 20*np.log10(d[-1])
        z[0] = z[0] + Ha
        z[-1] = z[-1] + Hb
        if method == 'deygout':
            LDiff = deygout(d, z, fGHz, L0, 0)
        if method == 'bullington':
            LDiff, deq, heq = bullington(d, z, fGHz)
        L[ip, :] = LDiff+LOS
    return(L)


def cover(X, Y, Z, Ha, Hb, fGHz, K, method='deygout'):
    """ outdoor coverage on a region

    Parameters
    ----------

    X : np.array (Nphi,Nr)
        cartesian coordinate grid
    Y : np.array (Nphi,Nr)
        cartesian coordinate grid
    Z : np.array (Nphi,Nr)
        height (meters)

    Ha : float
    Hb : float
    fGHz : np.array (,Nf)
        frequency in GHz
    method : 'deygout' | 'bullington'

    Returns
    -------

    L : Losses (dB)


    """
    Nphi, Nr = Z.shape

    if (type(fGHz) == float):
        fGHz = np.array([fGHz])

    Nf = len(fGHz)
    L = np.zeros((Nphi, Nr, Nf))
    L0 = np.zeros(Nf)
    # loop over azimut
    for ip in range(Nphi):
        # loop over range
        # il : 2 ... Nr-2
        # uk : 0 ....Nr-1
        for il in np.arange(2, Nr-1):
            uk = np.arange(0, il+1)
            z = np.empty(len(uk))
            x = X[ip, uk]
            y = Y[ip, uk]
            z[uk] = Z[ip, uk]
            d = np.sqrt((x-x[0])**2+(y-y[0])**2)
            # effect of refraction in equivalent earth curvature
            dh = d*(d[::-1])/(2*K*6375e3)
            z = z + dh
            LOS = 32.4 + 20*np.log10(fGHz) + 20*np.log10(d[-1])
            z[0] = z[0] + Ha
            z[-1] = z[-1] + Hb
            if method == 'deygout':
                LDiff = deygout(d, z, fGHz, L0, 0)
            if method == 'bullington':
                LDiff, deq, heq = bullington(d, z, fGHz)
            L[ip, il, :] = LDiff[None, :]+LOS[None,:]
    return(L)

def cover2(X, Y, Z, Ha, Hb, fGHz, K, method='deygout'):
    """ outdoor coverage on a region

    Parameters
    ----------

    X : np.array (Nphi,Nr)
        cartesian coordinate grid
    Y : np.array (Nphi,Nr)
        cartesian coordinate grid
    Z : np.array (Nphi,Nr)
        height (meters)

    Ha : float
    Hb : float
    fGHz : np.array (,Nf)
        frequency in GHz
    method : 'deygout' | 'bullington'

    Returns
    -------

    L : Losses (dB)


    """
    Nphi, Nr = Z.shape

    if (type(fGHz) == float):
        fGHz = np.array([fGHz])

    Nf = len(fGHz)
    L = np.zeros((Nphi, Nr, Nf))
    L0 = np.zeros(Nf)
    # loop over azimut
    for il in np.arange(2, Nr-1):
        dmax = LOS = 32.4 + 20*np.log10(fGHz) + 20*np.log10(d[-1])

        for ip in range(Nphi):
        # loop over range
        # il : 2 ... Nr-2
        # uk : 0 ....Nr-1
            uk = np.arange(0, il+1)
            z = np.empty(len(uk))
            x = X[ip, uk]
            y = Y[ip, uk]
            z[uk] = Z[ip, uk]
            d = np.sqrt((x-x[0])**2+(y-y[0])**2)
            # effect of refraction in equivalent earth curvature
            dh = d*(d[::-1])/(2*K*6375e3)
            z = z + dh
            z[0] = z[0] + Ha
            z[-1] = z[-1] + Hb
            if method == 'deygout':
                LDiff = deygout(d, z, fGHz, L0, 0)
            if method == 'bullington':
                LDiff, deq, heq = bullington(d, z, fGHz)
            L[ip, il, :] = LDiff[None, :]+LOS[None,:]
    return(L)

def deygout(d, height, fGHz, L, depth):
    """ Deygout attenuation

    Parameters
    ----------

    d : np.array (,N)
        horizontal distance
    height : np.array (,N)
        height profile
    fGHz : np.array (,Nf)
        frequency GHz
    L : np.array (,Nf)
        Additional Loss
    depth : recursive depth

    Notes
    -----

    This function is recursive

    """
    lmbda = 0.3/fGHz
    L0 = np.zeros(len(fGHz))
    depth = depth+1
    N = len(height)
    if depth < 3:
        if N > 3:
            u = np.arange(N)/(N-1.0)  # float
            # l : straight line between termination (LOS)
            l = (height[0])*(1-u)+(height[-1])*u
            # h excludes termination points
            h = height[1:-1] - l[1:-1]
            # Fresnel parameter (engagement)
            nu = h[:, None] * np.sqrt((2/lmbda[None, :]) *
                             (1/d[1:-1, None]+1/(d[-1]-d[1:-1, None])))
            imax = np.unique(np.nanargmax(nu, axis=0))[0]
            numax = nu[imax, :]
        else:
            numax = -10*np.ones(len(fGHz))
        if (numax > -0.78).any():
            w = numax - 0.1
            L = L + np.maximum(6.9 + 20*np.log10(np.sqrt(w**2+1)+w), 0)
            # left link
            height1 = height[0:imax+2]
            d1 = d[0:imax+2]
            Ll = deygout(d1, height1, fGHz, L0, depth)
            # right link
            height2 = height[imax+1:]
            d2 = d[imax+1:]
            Lr = deygout(d2, height2, fGHz, L0, depth)
            # add losses
            L = L + Lr + Ll

    return(L)

def bullington(d, height, fGHz):
    """ edges attenuation with Bullington method

    Parameters
    ----------

    d : np.array
    height : np.array
        antenna height is includes in height[0] and height[-1]
    fGHz : np.array

    Returns
    -------

    L : np.array
        total loss

    """

    def recl(d, height):
        """ determine left interception point
        Parameters
        ----------
        d : np.array
        height : np.array

        """
        N = len(height)
        u = np.arange(N)/(N-1.)
        # l : straight line between termination (LOS)
        l = height[0]*(1-u)+(height[-1])*u
        h = height - l
        # imax : index of the maximum height offset
        imax = np.argmax(h)
        if imax>0:
            # hmax : maximum height offset
            hmax = h[imax]
            # parameterization from 0 to imax
            ul = np.arange(imax)/(imax-1.)
            # straight line
            dhl = h[0]*(1-ul) + hmax*ul
            # el : offset if <0 split again
            el = dhl - h[0:imax]
            if np.min(el) < 0:
                u, v = recl(d[0:imax+1], height[0:imax+1])
            else:
                u = d[0:imax+1]
                v = h[0:imax+1]
        else:
            u = d[0:1]
            v = d[0:1]
        return(u, v)
    #if min(er)<0:
    #    u,v = rec(d[imax-1:-1],dhl)
    #else:
    #er  = dhr - h[imax-1:-1]
    def recr(d, height):
        """ determine the right interception point
        """
        N = len(height)
        u = np.arange(N)/(N-1.)
        l = height[0]*(1-u)+(height[-1])*u
        h = height - l
        imax = np.argmax(h)
        hmax = h[imax]
        ur = np.arange(N-imax)/(N-imax-1.)
        dhr = hmax*(1-ur) + h[-1]*ur

        er = dhr - h[imax:]
        if np.min(er) < 0:
            u, v = recr(d[imax:],h[imax:])
        else:
            u = d[imax:]
            v = h[imax:]

        return(u,v)
        #if min(er)<0:
    #    u,v = rec(d[imax-1:-1],dhl)
    #else:
    #er  = dhr - h[imax-1:-1]

    lmbda = 0.3/fGHz
    u = np.arange(len(height))/(len(height)-1.)
    l = (height[0])*(1-u)+(height[-1])*u
    h = height - l
    if (h>0).any():
        ul, vl = recl(d, height)
        ur, vr = recr(d, height)
        idtx = len(ul)
        idrx = len(h) - len(ur)
        dtx = d[idtx]
        drx = d[-1]-d[idrx]
        htx = h[idtx-1]
        hrx = h[idrx]
        deq = (dtx*hrx)*d[-1]/(drx*htx+dtx*hrx)
        heq = deq*(htx/dtx)
    else:
        heq = -np.min(np.abs(h[1:-1]))
        ieq = np.where(h==heq)[0][0]
        deq = d[ieq]

    nu = heq*np.sqrt((2/lmbda)*(1/deq+1/(d[-1]-deq)))
    w = nu - 0.1
    L = np.maximum(6.9 + 20*np.log10(np.sqrt(w**2+1)+w), 0)
    return(L, deq, heq)

def two_rays_flatearth(fGHz, **kwargs):
    """
    Parameters
    ----------

    p0 : transmitter position
        (3 x Np1) array or (2,) array
    p1 : receiver position
        (3 x Np2) array or (2,) array


    OR :

    d : distance between Tx and Rx
        (Np1,)
    ht : Tx height
    hr : Rx height
        (Np1)
    GtdB : float (0)
        Transmitter Antenna Gain (dB)
    GrdB : float(0)
        Receiver Antenna Gain (dB)
    fGHz : float (2.4)
        frequency (GHz)
    gamma : complex (-1.+0.j)
        Reflexion coeff


    dB : boolean (True)
        return result in d


    Returns
    -------

    P :
        received power

    Examples
    --------
    .. plot::
        :include-source:

        >>> from pylayers.antprop.loss import *
        >>> NPT=10000
        >>> x=np.array([0,0,8])
        >>> x=x.reshape(3,1)
        >>> y = np.ones((3,NPT))
        >>> y[0,:]=0
        >>> y[1,:]=np.arange(NPT)
        >>> y[2,:]=2
        >>> g0=1
        >>> g1=1
        >>> fGHz=2.4
        >>> PL2R=two_rays_flatearth(p0=x,p1=y,fGHz=fGHz,GtdB=g0,GrdB=g1)
        >>> PL1R = PL(fGHz,x,y,2)
        >>> plt.semilogx(PL2R,label='two-ray model')
        >>> plt.semilogx(-PL1R[0,:],label='one slope model')
        >>> plt.axis([10,NPT,-150,-50])
        >>> plt.title('Loss 2-rays model vs one slope model')
        >>> plt.xlabel('distance (m)')
        >>> plt.ylabel('Loss Pr/Pt (dB)')
        >>> plt.legend()
        >>> plt.show()

        >>> d=np.arange(1,1000)
        >>> PL2Rd = two_rays_flatearth(d=d,ht=np.array([5]),hr=np.array([10]),fGHz=fGHz,GtdB=g0,GrdB=g1)
        >>> plt.semilogx(PL2Rd,label='two-ray model')
        >>> plt.semilogx(-PL1R[0,:],label='one slope model')
        >>> plt.axis([10,NPT,-150,-50])
        >>> plt.title('Loss 2-rays model vs one slope model')
        >>> plt.xlabel('distance (m)')
        >>> plt.ylabel('Loss Pr/Pt (dB)')
        >>> plt.legend()
        >>> plt.show()

    References
    ----------

    https://en.wikipedia.org/wiki/Two-ray_ground-reflection_model#As_a_case_of_log_distance_path_loss_model
    http://morse.colorado.edu/~tlen5510/text/classwebch3.html#x15-590003.3.3

    """


    defaults = { 'p0':np.array((0,0,10)),
                 'p1':np.array((0,10,10)),
                 'd':[],
                 'ht':10,
                 'hr':10,
                 'GtdB':0.,
                 'GrdB':0.,
                 'gamma': -1.+0.j,
                 'pol':'v',
                 'eps' :[],
                 'sig':0.,
                 'dB':True
               }

    for k in defaults:
       if k not in kwargs:
           kwargs[k]=defaults[k]

    GtdB=kwargs.pop('GtdB')
    GrdB=kwargs.pop('GrdB')

    Gt = 10**((1.*GtdB)/10.)
    Gr = 10**((1.*GrdB)/10.)
    gamma=kwargs.pop('gamma')
    pol=kwargs.pop('pol')
    eps=kwargs.pop('eps')
    sig=kwargs.pop('sig')


    if kwargs['d'] == []:
        p0=kwargs['p0']
        p1=kwargs['p1']
        assert p0.shape[0] == 3, 'p0 is not 3D'
        assert p1.shape[0] == 3, 'p1 is not 3D'


        if len(p0.shape) == 1:
            p0=p0.reshape(p0.shape[0],1)
        if len(p1.shape) == 1:
            p1=p1.reshape(p1.shape[0],1)

        p0=1.*p0
        p1=1.*p1

        ht = p0[2,:]
        hr = p1[2,:]
        dloss = np.sqrt(np.sum((p0-p1)**2,axis=0)) #l0
    else:

        dloss=kwargs['d']
        ht=kwargs['ht']
        hr=kwargs['hr']

    Gt = 10**((1.*Gt)/10.)
    Gr = 10**((1.*Gr)/10.)

    d0 = np.sqrt( dloss**2 - 1.*(ht-hr)**2 ) # d0
    dref = np.sqrt(d0**2+1.*(ht+hr)**2) #l0'


    if eps != []:
        psy = np.arcsin((ht+hr)/dref)
        er = eps - 60.j*sig*0.3/fGHz
        if pol == 'v':
            Z = (1./er)* np.sqrt(er-np.cos(psy)**2)
        elif pol == 'h':
            Z = np.sqrt(er-np.cos(psy)**2)

        gamma = (np.sin(psy)-Z)/((np.sin(psy)+Z))



    deltad = dref-dloss
    deltaphi = (2*np.pi*fGHz*deltad)/0.3
    E= (0.3/(4*np.pi*fGHz)) * (np.sqrt(Gt*Gr)/dloss + gamma * np.sqrt(Gr*Gr)*(np.exp(-1.j*deltaphi))/dref)
    P = abs(E)**2

    # import ipdb
    # ipdb.set_trace()
    if kwargs['dB'] :
        return 10*np.log10(P)
    else:
        return P


def lossref_compute(P,h0,h1,k=4/3.) :
    """
    compute loss and reflection rays on curved earth

    Parameters
    ----------

    P : float |list

        if len(P) == 1 => P is a distance
        if len(P) == 4 => P is a list of [lon0,lat0,lon1,lat1]

        where :
        lat0 : float |string
            latitude first point (decimal |deg min sec Direction)
        lat1 : float |string
            latitude second point (decimal |deg min sec Direction)
        lon0 : float |string
            longitude first point (decimal |deg min sec Direction)
        lon1 : float |string
            longitude second point (decimal |deg min sec Direction)
    h0 : float:
        height of 1st point
    h1 : float:
        height of 2nd point
    k : electromagnetic earth factor


    Returns
    -------

    dloss : float
        length of direct path (meter)
    dref : float
        length of reflective path (meter)
    psy : float
        Reflection angle

    References
    ----------
    B. R. Mahafza, Radar systems analysis and design using MATLAB, Third edition. Boca Raton; London: CRC/Taylor & Francis, chapter 8, 2013.

    """


    if isinstance(P,float) or isinstance(P,int) :
        #P is a distance
        r=P
        mode = 'dist'
    elif isinstance(P,np.ndarray) or isinstance(P,list):
        if len(P) == 1:
            #P is a distance
            r=P
            mode = 'dist'
        elif len(P) == 4:
            #P is a lonlat
            lat0=P[0]
            lon0=P[1]
            lat1=P[2]
            lon1=P[2]
            mode = 'lonlat'
        else :
            raise AttributeError('P must be a list [lat0,lon0,lat1,lon0] or a distance')
    else :
        raise AttributeError('Invalid P format ( list |ndarray )')


    # if h0<h1:
    #     h1,h0 = h0,h1

    r0 = 6371.e3 # earth radius
    re = k*r0 # telecom earth radius


    if mode == 'lonlat':
        # r = distance curvilignenp.arcsin((h1/R1)-R1/(2.*re)) entre TXetRX / geodesic
        r = gu.distance_on_earth(lat0, lon0, lat1, lon1)
    else :
        r=P

    r=1.*r

    # import ipdb
    # ipdb.set_trace()
    p = 2/(np.sqrt(3))*np.sqrt(re*(h0+h1)+(r**2/4.)) #eq 8.45
    eps = np.arcsin(2*re*r*(h1-h0)/p**3) # eq 8.46



    #distance of reflection on curved earth
    r1 = r/2 - p*np.sin(eps/3) #eq 8.44

    r2 = r -r1

    phi1 = r1/re #8.47
    phi2 = r2/re # 8.48

    R1 = np.sqrt(h0**2+4*re*(re+h0)*(np.sin(phi1/2))**2) # 8.51
    R2 = np.sqrt(h1**2+4*re*(re+h1)*(np.sin(phi2/2))**2) #8.52

    Rd = np.sqrt((h1-h0)**2+4*(re+h1)*(re+h0)*np.sin((phi1+phi2)/2.)**2) # 8.53
    # tangente angle on earth
    psy = np.arcsin((h1/R1)-R1/(2.*re)) #eq 8.55
    deltaR = 4*R1*R2*np.sin(psy)**2/(R1+R2+Rd)

    dloss = Rd
    dref = R1+R2

    return psy,dloss,dref

def two_rays_curvedearthold(P,h0,h1,fGHz=2.4,**kwargs):
    """


    Parameters
    ----------

    P : float |list

        if len(P) == 1 => P is a distance
        if len(P) == 4 => P is a list of [lon0,lat0,lon1,lat1]

        where :
        lat0 : float |string
            latitude first point (decimal |deg min sec Direction)
        lat1 : float |string
            latitude second point (decimal |deg min sec Direction)
        lon0 : float |string
            longitude first point (decimal |deg min sec Direction)
        lon1 : float |string
            longitude second point (decimal |deg min sec Direction)
    h0 : float:
        height of 1st point
    h1 : float:
        height of 2nd point
    fGHz : float
        frequency (GHz)


    k : float
        electromagnetic earth factor
    GtdB : float
        Transmitter Antenna Gain (dB)
    GrdB : float
        Receiver Antenna Gain (dB)
    gamma : complex (-1.+0.j)
        Reflexion coeff if eps and sig are not precised

    'pol': string ('v')
        polarization ('v'|'h')
    'eps' : float ([])
        lossless relative permittivity [],
    'sig': float (0.)
        conductivity


    dB : boolean (True)
        return result in dB


    Returns
    -------

    P :
        received power

    Examples
    --------

    .. plot::
        :include-source:

        >>> from pylayers.antprop.loss import *
        >>> import matplotlib.pyplot as plt
        >>> fGHz=2.4
        >>> p0=np.array(([0,0,20]))
        >>> p1=np.array(([0,1,20]))
        >>> p0=p0.reshape(3,1)
        >>> p1=p1.reshape(3,1)
        >>> TRF = [] #Two Ray model on flat earth
        >>> TRC = [] #Two Ray model on curved earth
        >>> PLoss=[]
        >>> for d in np.arange(1,10000,1):
        >>>     p1[1,:]=d
        >>>     TRF.append(two_rays_flatearth(p0[:,0],p1[:,0],fGHz,GtdB=0.,GrdB=0.,))
        >>>     TRC.append(two_rays_curvedearth(d,p0[2,:],p1[2,:],fGHz))
        >>>     PLoss.append(PL(fGHz, p0[:,0],p1[:,0], n=2.0, dB=True, d0=np.array([1])))
        >>> PLoss=np.array(PLoss)[:,0,0]
        >>> plt.semilogx(TRF,label='two-rays model flat earth')
        >>> plt.semilogx(TRC,label='two-rays model curved earth')
        >>> plt.semilogx(-PLoss,label='Path Loss')
        >>> plt.legend()
        >>> plt.show()

    """


    defaults = { 'GtdB':0.,
                 'GrdB':0.,
                 'k':4/3.,
                 'gamma': -1.+0.j,
                 'pol':'v',
                 'eps' :[],
                 'sig':0.,
                 'mode':'PL',
                 'dB':True
               }

    for k in defaults:
        if k not in kwargs:
            kwargs[k]=defaults[k]

    GtdB=kwargs.pop('GtdB')
    GrdB=kwargs.pop('GrdB')

    Gt = 10**((1.*GtdB)/10.)
    Gr = 10**((1.*GrdB)/10.)
    k=kwargs.pop('k')
    gamma=kwargs.pop('gamma')
    pol=kwargs.pop('pol')
    eps=kwargs.pop('eps')
    sig=kwargs.pop('sig')


    h0=1.*h0
    h1=1.*h1

    psy,dloss,dref = lossref_compute(P,h0,h1,k)

    if eps != []:
        er = eps  - 60.j*sig*0.3/fGHz
        if pol == 'v':
            Z= (1./er)* np.sqrt(er-np.cos(psy)**2)
        elif pol == 'h':
            Z= np.sqrt(er-np.cos(psy)**2)

        gamma = (np.sin(psy)-Z)/((np.sin(psy)+Z))

    deltad= dref-dloss
    deltaphi = (2*np.pi*fGHz*deltad)/0.3
    E= (0.3/(4*np.pi*fGHz) ) *(np.sqrt(Gt*Gr)/dloss + gamma * np.sqrt(Gr*Gr)*(np.exp(-1.j*deltaphi))/dref)
    P = abs(E)**2

    # import ipdb
    # ipdb.set_trace()
    if kwargs['dB'] :
        return 10*np.log10(P)
    else:
        return P


def visuPts(S,nu,nd,Pts,Values,fig=[],sp=[],vmin=0,vmax=-1,label=' ',tit='',size=25,colbar=True,xticks=False):
    """ visuPt  : Visualization of values a given points

    Parameters
    ----------

    S       : Simul
    nu      : useful Points
    nd      : Points deleted
    Pts     : Points coordinates
    Value
    """

    vx  = Pts[nu,0]
    vy  = Pts[nu,1]
    vxd = Pts[nd,0]
    vyd = Pts[nd,1]

    if vmax<vmin:
        #vmin = min(Values)
        vmin = -150
        vmax = max(Values)
    S.L.showGs()
    if xticks:
        for loc, spine in sp.spines.iteritems():
              if loc in ['left','bottom']:
                spine.set_position(('outward',10)) # outward by 10 points
                sp.yaxis.set_ticks_position('left')
                sp.xaxis.set_ticks_position('bottom')
              elif loc in ['right','top']:
                 spine.set_color('none') # don't draw spine
              else:
                raise ValueError('unknown spine location: %s'%loc)
    else:
        for loc, spine in sp.spines.iteritems():
              if loc in ['left']:
                spine.set_position(('outward',10)) # outward by 10 points
                sp.yaxis.set_ticks_position('left')
              elif loc in ['right','top','bottom']:
                 spine.set_color('none') # don't draw spine
                 sp.xaxis.set_ticks([])
              else:
                raise ValueError('unknown spine location: %s'%loc)



        # no xaxis ticks
        #ax.xaxis.set_ticks([])
    #sp.spines['left'].set_position('center')
        #sp.spines['left'].set_color('none')
    #sp.spines['right'].set_position('center')
        #sp.spines['right'].set_color('none')
    #sp.spines['bottom'].set_position('center')
    #sp.xaxis.set_ticks_position('bottom')
    #sp.yaxis.set_ticks_position('left')
        #sp.spines['bottom'].set_color('none')
    #sp.spines['top'].set_position('center')
        #sp.spines['top'].set_color('none')
    #
    # Rooms annotation
    #
    annotate('R 8',xy=(-19,14.1))
    annotate('R 9',xy=(-24.5,6.5))
    annotate('R 14',xy=(-20,6.5))
    annotate('R 10',xy=(-16.5,6.5))
    annotate('R 7',xy=(-10.5,14.1))
    annotate('R 11',xy=(-2.5,13.5))
    annotate('R 12',xy=(-8.7,6.5))
    annotate('R 13',xy=(-5.2,14.5))
    annotate('R 1',xy=(3.5,8))
    annotate('R 2',xy=(1.5,13.8))
    annotate('R 6',xy=(-3.6,6.5))


    n=scatter(vx,vy,marker='o',c=Values,s=size,vmin=vmin,vmax=vmax)
    n.set_edgecolor('face')
#    m=scatter(vxd,vyd,marker='o',c='k',s=size)
#    m.set_edgecolor('face')
    axis('scaled')
    title(tit)
    ylabel('meters')
    if xticks:
        xlabel('meters')
    if colbar:
        cbar=colorbar(orientation='vertical')
        cbar.set_label(label)

def cdf(x,colsym="",lab="",lw=4):
    """ plot the cumulative density function

    Parameters
    ----------

    x : np.array()
    colsym : string
    lab : string
    lw : int
        linewidth

    Examples
    --------

    .. plot::
        :include-source:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt 
        >>> x = np.random.randn(10000)
        >>> cdf(x)
        >>> plt.show()

    """
    rcParams['legend.fontsize']=20
    rcParams['font.size']=20

    x  = np.sort(x)
    n  = len(x)
    x2 = np.repeat(x, 2)
    y2 = np.hstack([0.0, repeat(np.arange(1,n) / float(n), 2), 1.0])
    plt.plot(x2,y2,colsym,label=lab,linewidth=lw)
    plt.grid('on')
    plt.legend(loc=2)
    plt.xlabel('Ranging Error[m]')
    plt.ylabel('Cumulative Probability')

if __name__=="__main__":
    doctest.testmod()
