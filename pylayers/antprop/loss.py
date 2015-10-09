# -*- coding:Utf-8 -*-
#from numpy import *
"""

.. currentmodule:: pylayers.antprop.loss

Loss module
===========

This module implements path loss models for various situations.


.. autosummary::
    :toctree: generated/

    PL0
    PL
    Losst
    Loss0_v2
    Loss0
    Loss_diff
    OneSlopeMdl
    cost231
    cost259
    cost2100
    Dgrid_points
    Dgrid_zone
    calnu
    visuPts
    cdf

"""

import doctest
import logging
import numpy as np
from scipy import io
import matplotlib.pylab as plt
import pylayers.simul.simulem
import pylayers.measures.mesuwb
import pylayers.gis.gisutil as gu

import pdb

def PL0(fGHz,GtdB=0,GrdB=0):
    """  Path Loss at frequency fGHZ @ 1m 

    Parameters
    ----------

    fGHz: float
          frequency GHz
    GtdB: float
          transmitting antenna gain dB (default 0 dB)
    GrdB: float
          receiving antenna gain dB (default 0 dB)

    Returns
    -------

    PL0 : float
          path @ 1m

    Notes
    -----

    .. math:: PL_0 = -20 log_{10}(\\frac{\\lambda}{4\\pi}) - GtdB -GrdB

    Examples
    --------

    >>> fGHz  = 2.4
    >>> PL = PL0(fGHz)
    >>> assert (PL<41)&(PL>40),"something wrong"

    """

    ld  = 0.3/fGHz
    PL0 = -20*np.log10(ld/(4.0*np.pi))-GtdB-GrdB

    return PL0

def Dgrid_points(points,Px):
    """ distance point to grid

    Parameters
    ----------

    points : np.array
             grid Np x 2 array

    Px     : np.array
             point 2 x 1  array

    """

    Dx = points[:,0] - Px[0]
    Dy = points[:,1] - Px[1]

    D = np.sqrt( Dx*Dx + Dy*Dy )

    return(D)

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

    Build the distance matrix between Tx and points in the zone

    Notes
    -----

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


    Notes
    -----

    http://

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

    Returns
    -------

    L : Attenuation (dB)


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
    dm  = np.sqrt((pBS-pMS)*(pBS-pMS))
    if (typ=='small'):
       CH = (1.1*np.log10(fGHz*1000)-0.7)*hMS-(1.56*np.log10(fGHz*1000)-0.8)
    if (typ=='big'):
        if fGHz<0.2:
            CH = 8.29*(np.log10(1.54*h_MS)**2)-1.1
        else:# valid until 1.5GHz
            CH = 3.2*(np.log10(11.75*h_MS)**2)-4.97

    L = 69.55+26.16*np.log10(fGHz*1000)-13.82*np.log10(hBS)+(44.9-6.55*np.log10(hBS))*np.log10(dm/1000.)-CH

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

def PL(fGHz,pts,p,n=2.0,dB=True):
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

    PL = np.array([PL0(fGHz)])[:,np.newaxis] + 10*n*np.log10(D)[np.newaxis,:]

    if not dB:
        PL=10**(-PL/10)

    return(PL)

def Losst(L,fGHz,p1,p2,dB=True):
    """  calculate Losses between links p1 p2

    Parameters
    ----------

    L   : Layout object
    fGHz : np.array
           frequency GHz
    p1 : source points
        (2 x Np1) array or (2,) array
    p2 : observation point
        (2 x Np2) array or (2,) array
    dB : boolean

    Examples
    --------

    .. plot::
        :include-source:

        >>> import matplotlib.pyplot as plt
        >>> from pylayers.simul.simulem import *
        >>> from pylayers.measures.mesuwb import *
        >>> from pylayers.antprop.loss import *
        >>> S = Simul()
        >>> S.layout('WHERE1.ini')
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
    data = L.angleonlink(p1,p2)

    # as many slabs as segments
    slabs = L.sla[data['s']]

    cslab = np.unique(slabs)

    LossWallo = np.zeros((len(fGHz),Nlink))
    LossWallp = np.zeros((len(fGHz),Nlink))
    EdWallo = np.zeros((len(fGHz),Nlink))
    EdWallp = np.zeros((len(fGHz),Nlink))

    for slname in cslab:
        # u index of slabs of name slname
        # data['a'][u] angle
        # data['s'][u] segment number including subsegment
        u = np.nonzero(slabs==slname)[0]
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
        Wallo = np.array(map(lambda x: np.sum(lko[:,indices[x]:indicep[x]],axis=1),irange)).T
        Wallp = np.array(map(lambda x: np.sum(lkp[:,indices[x]:indicep[x]],axis=1),irange)).T

        Edo = np.array(map(lambda x: np.sum(do[indices[x]:indicep[x]]),irange)).T
        Edp = np.array(map(lambda x: np.sum(dp[indices[x]:indicep[x]]),irange)).T

        LossWallo[:,involved_links] = LossWallo[:,involved_links] + Wallo
        LossWallp[:,involved_links] = LossWallp[:,involved_links] + Wallp

        EdWallo[:,involved_links] = EdWallo[:,involved_links] + Edo
        EdWallp[:,involved_links] = EdWallp[:,involved_links] + Edp

    if not dB:
        LossWallo = 10**(-LossWallo/10)
        LossWallp = 10**(-LossWallp/10)

    return(LossWallo,LossWallp,EdWallo,EdWallp)

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
    print shape(Lw)
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

    ld  = 0.3/fGHz
    nu  = h*np.sqrt(2*(d1+d2)/(ld*d1*d2))

    return(nu)

def showfurniture(fig,ax):
    """ show furniture (not the good module) 
    """
    #R1_A.show(fig,ax)
    #R1_B1.show(fig,ax)
    #R1_B2.show(fig,ax)
    #R2_A.show(fig,ax)
    #R6_A.show(fig,ax)
    #R6_B.show(fig,ax)
    #R6_C.show(fig,ax)
    R6_D.show(fig,ax)
    R6_E.show(fig,ax)
    R6_F.show(fig,ax)
    R6_G.show(fig,ax)
    R6_HA.show(fig,ax)
    R6_HB.show(fig,ax)
    R6_IA.show(fig,ax)
    R6_IB.show(fig,ax)
    R13_A.show(fig,ax)
    R13_B.show(fig,ax)
    R13_C.show(fig,ax)
    R12_A1.show(fig,ax)
    R12_A2.show(fig,ax)
    #R12_B.show(fig,ax)
    #R12_C1.show(fig,ax)
    #R12_C2.show(fig,ax)
    R12_D.show(fig,ax)
    #R7_A1.show(fig,ax)
    #R7_A2.show(fig,ax)
    #R7_B.show(fig,ax)
    R7_C1.show(fig,ax)
    R7_C2.show(fig,ax)
    R7_D.show(fig,ax)
    #R10_A.show(fig,ax)
    #R10_B.show(fig,ax)
    R10_C.show(fig,ax)
    R10_D.show(fig,ax)
    #R11_A.show(fig,ax)
    #R11_B.show(fig,ax)
    R11_C.show(fig,ax)
    R11_D1.show(fig,ax)
    R11_D2.show(fig,ax)
    R11_E1.show(fig,ax)
    R11_E2.show(fig,ax)
    R11_F.show(fig,ax)
    #R8_A.show(fig,ax)
    #R8_B.show(fig,ax)
    #R9_A1.show(fig,ax)
    #R9_A2.show(fig,ax)
    R9_B1.show(fig,ax)
    R9_B2.show(fig,ax)
    R9_B3.show(fig,ax)
    #R9_C.show(fig,ax)
    R9_D.show(fig,ax)
    R9_E1.show(fig,ax)
    R9_E2.show(fig,ax)
    R9_F.show(fig,ax)
    R9_G.show(fig,ax)
    axis('scaled')


def two_rays_flatearth(p0,p1,Gt,Gr,fGHz,gamma= -1.+0.j,mode='PL',dB=True):
    """
    Parameters
    ----------

    p1 : transmitter position
        (3 x Np1) array or (2,) array
    p2 : receiver position
        (3 x Np2) array or (2,) array

    Gt : float (0) 
        Transmitter Antenna Gain (dB)
    Gr : float(0)
        Receiver Antenna Gain (dB)
    fGHz : float (2.4)
        frequency (GHz)
    gamma : complex (-1.+0.j)
        Reflexion coeff

    mode : PL | E (default : PL)
        return Energy (E) or Path loss/power loss (PL)
    dB : boolean (True)
        return result in d


    Returns
    -------

    PL : 
        path loss w.r.t distance and frequency


    E  : 
        energy loss w.r.t distance and frequency



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
        >>> PL2R=two_rays_flatearth(x,y,g0,g1,fGHz)
        >>> PL1R = PL(fGHz,x,y,2)
        >>> plt.semilogx(PL2R,label='two-ray model')
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

    assert p0.shape[0] == 3, 'p0 is not 3D'
    assert p1.shape[0] == 3, 'p1 is not 3D'


    if len(p0.shape) == 1:
        p0=p0.reshape(p0.shape[0],1)
    if len(p1.shape) == 1:
        p1=p1.reshape(p1.shape[0],1)


    Gt = 10**((1.*Gt)/10.)
    Gr = 10**((1.*Gr)/10.)

    ht = p0[2,:]
    hr = p1[2,:]

    dloss = np.sqrt(np.sum((p0-p1)**2,axis=0)) # l0
    d0 = np.sqrt( dloss**2 - (ht-hr)**2 ) # d0
    dref = np.sqrt(d0**2+(ht+hr)**2) #l0'


    lossterm = np.exp(2.j*np.pi*dloss*fGHz/0.3) / (1.*dloss)
    refterm = np.exp(2.j*np.pi*dref*fGHz/0.3) / (1.*dref)


    E = np.sqrt(Gt*Gr)*0.3/(4*np.pi*fGHz)* (lossterm + gamma*refterm)

    if mode == 'E':
        return E
    if mode == 'PL':
        if dB :
            return 20*np.log10(E)
        else:
            return E**2


def lossref_compute(P,h0,h1,k=4/3.) :
    """
    compute loss and reflection rays on curved earth

    Parameters
    ----------

    P : float | list 

        if len(P) == 1 => P is a distance
        if len(P) == 4 => P is a list of [lon0,lat0,lon1,lat1]

        where :
        lat0 : float | string
            latitude first point (decimal | deg min sec Direction)
        lat1 : float | string
            latitude second point (decimal | deg min sec Direction)
        lon0 : float | string
            longitude first point (decimal | deg min sec Direction)
        lon1 : float | string
            longitude second point (decimal | deg min sec Direction)
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
    B. R. Mahafza, Radar systems analysis and design using MATLAB, Third edition. Boca Raton ; London: CRC/Taylor & Francis, chapter 8, 2013.

    """


    if isinstance(P,float) or isinstance(P,int) :
        # P is a distance
        r=P
        mode = 'dist'
    elif isinstance(P,np.ndarray) or isinstance(P,list):
        if len(P) == 1:
            # P is a distance
            r=P
            mode = 'dist'
        elif len(P) == 4:
            # P is a lonlat
            lat0=P[0]
            lon0=P[1]
            lat1=P[2]
            lon1=P[2]
            mode = 'lonlat'
        else :
            raise AttributeError('P must be a list [lat0,lon0,lat1,lon0] or a distance')
    else :
        raise AttributeError('Invalid P format ( list | ndarray )')


    # if h0<h1:
    #     h1,h0 = h0,h1

    r0 = 6371e3 # earth radius
    re = k*r0 # telecom earth radius


    if mode == 'lonlat':
        # r = distance curviligne entre TXetRX / geodesic
        r = gu.distance_on_earth(lat0, lon0, lat1, lon1)
    else :
        r=P

    # import ipdb
    # ipdb.set_trace()
    p = 2/(np.sqrt(3))*np.sqrt(re*(h0+h1)+(r**2/4.)) # eq 8.45
    eps = np.arcsin(2*re*r*(h1-h0)/p**3) # eq 8.46



    # distance of reflection on curved earth
    r1 = r/2 - p*np.sin(eps/3) # eq 8.44

    r2 = r -r1

    phi1 = r1/re # 8.47
    phi2 = r2/re # 8.48

    R1 = np.sqrt(h0**2+4*re*(re+h0)*(np.sin(phi1/2))**2) # 8.51
    R2 = np.sqrt(h1**2+4*re*(re+h1)*(np.sin(phi2/2))**2) #8.52

    Rd = np.sqrt((h1-h0)**2+4*(re+h1)*(re+h0)*np.sin((phi1+phi2)/2.)**2) # 8.53

    # tangente angle on earth
    psy = np.arcsin((h1/R1)-R1/(2*re)) # eq 8.55
    deltaR = 4*R1*R2*np.sin(psy)**2/(R1+R2+Rd)

    dloss = Rd
    dref = R1+R2

    return psy,dloss,dref

def two_ray_curvedearth(P,h0,h1,fGHz=2.4,**kwargs):
    """


    Parameters
    ----------

    P : float | list 

        if len(P) == 1 => P is a distance
        if len(P) == 4 => P is a list of [lon0,lat0,lon1,lat1]

        where :
        lat0 : float | string
            latitude first point (decimal | deg min sec Direction)
        lat1 : float | string
            latitude second point (decimal | deg min sec Direction)
        lon0 : float | string
            longitude first point (decimal | deg min sec Direction)
        lon1 : float | string
            longitude second point (decimal | deg min sec Direction)
    h0 : float:
        height of 1st point 
    h1 : float:
        height of 2nd point 
    fGHz : float
        frequency (GHz)


    k : float
        electromagnetic earth factor
    Gt : float
        Transmitter Antenna Gain (dB)
    Gr : float
        Receiver Antenna Gain (dB)
    gamma : complex (-1.+0.j)
        Reflexion coeff if eps and sig are not precised

    'pol': string ('v')
        polarization ('v'|'h')
    'eps' : float ([])
        lossless relative permittivity [],
    'sig': float (0.)
        conductivity 


    mode : PL | E (default : PL)
        return Energy (E) or Path loss/power loss (PL)
    dB : boolean (True)
        return result in dB


    Example
    -------

    .. plot::
        :include-source:
        
            >>> from pylayers.antprop.loss import *
            >>> import matplotlib.pyplot as plt
            >>> fGHz=2.4
            >>> p0=np.array(([0,0,20]))
            >>> p1=np.array(([0,1,20]))
            >>> p0=p0.reshape(3,1)
            >>> p1=p1.reshape(3,1)
            >>> OR = [] # One Ray model
            >>> TRF = [] # Two Ray model on flat earth
            >>> TRC = [] # Two Ray model on curved earth
            >>> for d in np.arange(1,10000,1):
            >>>     p1[1,:]=d
            >>>     OR.append(-PL(fGHz,p0[:2,:],p1[:2,:],2)[0])
            >>>     TRF.append(two_rays_flatearth(p0[:,0],p1[:,0],0.,0.,fGHz))
            >>>     TRC.append(two_ray_curvedearth(d,p0[2,:],p1[2,:],fGHz))
            >>> plt.semilogx(TRF,label='two-ray model flat earth')
            >>> plt.semilogx(TRC,label='two-ray model curved earth')
            >>> plt.semilogx(OR,label='one-ray model')
            >>> plt.legend()
            >>> plt.show()

    """


    defaults = { 'Gt':0.,
                 'Gr':0.,
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

    Gt=kwargs.pop('Gt')
    Gr=kwargs.pop('Gr')

    Gt = 10**((1.*Gt)/10.)
    Gr = 10**((1.*Gr)/10.)
    k=kwargs.pop('k')
    gamma=kwargs.pop('gamma')
    pol=kwargs.pop('pol')
    eps=kwargs.pop('eps')
    sig=kwargs.pop('sig')


    psy,dloss,dref = lossref_compute(P,h0,h1,k)

    if eps != []:
        er = eps  - 60.j*sig*0.3/fGHz
        if pol == 'v':
            Z= (1./er)* np.sqrt(er-np.cos(psy)**2)
        elif pol == 'h':
            Z= np.sqrt(er-np.cos(psy)**2)

        gamma = (np.sin(psy)-Z)/((np.sin(psy)+Z))


    lossterm = np.exp(2.j*np.pi*dloss*fGHz/0.3) / (1.*dloss)
    refterm = np.exp(2.j*np.pi*dref*fGHz/0.3) / (1.*dref)

    E = np.sqrt(Gt*Gr)*0.3/(4*np.pi*fGHz)* (lossterm + gamma*refterm)

    if kwargs['mode'] == 'E':
        return E
    if kwargs['mode'] == 'PL':
        if kwargs['dB'] :
            return 20*np.log10(E)
        else:
            return E**2


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

    >>> import numpy as np

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
