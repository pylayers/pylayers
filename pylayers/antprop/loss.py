# -*- coding:Utf-8 -*-
#from numpy import *
"""

Loss module
===========

.. autosummary::
    :toctree: generated/

    PL0
    PL
    Losst
    Loss0_v2
    Loss_mur_the
    Loss0
    Loss_diff
    Loss_obstacle
    Loss_2obstacle
    OneSlopeMdl
    Dgrid_points
    Dgrid_zone
    calnu
    Carretosegment
    Intersection
    Dis
    Interline
    showfurniture
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

    >>> from pylayers.antprop.loss import *
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> Nlink = 100
    >>> hBS = 300
    >>> hMS = 1.5
    >>> hroof = 40*np.random.rand(Nlink)
    >>> wr = 10*np.ones(Nlink)
    >>> phir = 90*np.random.rand(Nlink)
    >>> pMS = np.vstack((np.linspace(10,2500,Nlink),np.zeros(Nlink),hMS*np.ones(Nlink)))
    >>> pBS = np.vstack((np.zeros(Nlink),np.zeros(Nlink),hBS*np.ones(Nlink)))
    >>> fMHz = np.linspace(700,1900,120)
    >>> pl = cost231(pBS,pMS,hroof,phir,wr,fMHz)
    >>> im = plt.imshow(pl)
    >>> cb = plt.colorbar()
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
    """
    """
    hBS = pBS[3,:]
    hMS = pMS[3,:]
    dm  = np.sqrt((pBS-pMS)*(pBS-pMS))
    lmbd = 300/fMHz
    pl = 10*2.6*np.log10(dm)+20*log10(4*np.pi/lmbd) 
    if not dB:
        pl = 10**(-pl/20.);
    return(pl)

def cost2100(pMS,pBS,fGHz,nfloor=1,dB=True):
    """
    """
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

    Returns
    -------

    PL : np.array
         path loss w.r.t distance and frequency

    """
    shp = np.shape(p)
    assert(shp[0]==2)

    D = np.sqrt(np.sum((pts-p)**2,axis=0))

    # f x grid x ap
    PL = PL0(fGHz)[:,np.newaxis] + 10*n*np.log10(D)[np.newaxis,:]

    if not dB:
        PL=10**(-PL/10)

    return(PL)

def Losst(L,fGHz,p1,p2,dB=True):
    """  Calculate Loss between links p1  p2

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
#    i = 0
#    for k in seglist:
#        if k != 0:
#            try:
#                # TODO use z to determine ss_name
#                name = L.Gs.node[k]['ss_name'][-2]
#            except:
#                name = L.Gs.node[k]['name']
#            #if k in S.indoor.ce.keys():
#            #if k in S.L.ce.keys():
#            # nom du sous-segment  
#            #    indss = S.L.ce[k][0]
#            #    name  = S.L.sl.di[indss]
#            #    print name
#            #"else:  
#            # nom du segment   
#            #    name = S.L.Gs.node[k]['name'] 
#            the = theta[i]
#            # todo : comparison multiwall th=0 th=variable
#            # comparison mesurement
#            #the = 0
#
#            i   = i + 1
#            #
#            # Loss0 du slab
#            #
#            lko,lkp  = L.sl[name].losst(fGHz,the)
#            do , dp  = L.sl[name].excess_grdelay(theta=the)
#            edo = edo - np.mean(do)
#            edp = edp - np.mean(dp)
##           print lko
##           print lkp
#            Lo   = Lo + lko[0]
#            Lp   = Lp + lkp[0]
#    Lwo = np.hstack((Lwo,Lo))
#    Lwp = np.hstack((Lwp,Lp))
#    Edo = np.hstack((Edo,edo))
#    Edp = np.hstack((Edp,edp))

#    return(Lwo,Lwp,Edo,Edp)
def Loss0_v2(L,Pts,fGHz,p):
    """ calculate loss

    Parameters
    ----------

    L   : Layout
          Layout object

    Pts : observation grid 
        (Np x 2) array

    fGHz : np.array
           frequency 
    p  : point
        source points

    Returns
    -------

    Lwo : Losses in wall polarization o
    Lwp : Losses in wall polarization p
    Edo :  polarization o
    Edp :  polarization p

    Examples
    --------

    .. plot::
        :include-source:

#        >>> import matplotlib.pyplot as plt
#        >>> from pylayers.simul.simulem import *
#        >>> from pylayers.measures.mesuwb import *
#        >>> from pylayers.antprop.loss import *
#        >>> S = Simul()
#        >>> S.layout('Where1.ini')
#        >>> fGHz = 4
#        >>> Tx,Rx = ptw1()
#        >>> Lwo,Lwp,Edo,Edp = Loss0_v2(S.L,Tx,fGHz,Rx[1,0:2])
#        >>> fig,ax = S.L.showGs()
#        >>> tit = plt.title('test Loss0_v2')
#        >>> sc2 = ax.scatter(Rx[1,0],Rx[1,1],s=20,marker='x',c='k')
#        >>> sc1 = ax.scatter(Tx[:,0],Tx[:,1],s=Edo,c=Edo,linewidth=0)
#        >>> plt.show()

    Notes
    -----

    DEPRECATED : Use losst instead

    """

    logging.warning('DEPRECATED function')

    N   = np.shape(Pts)[0]
    Lwo = np.array([])
    Lwp = np.array([])
    Edo = np.array([])
    Edp = np.array([])

    for i in range(N):
        Lo = 0.0
        Lp = 0.0
        edo = 0.0
        edp = 0.0
        pi = Pts[i,:]
        seglist,theta = L.angleonlinkold(p,pi)
        i = 0
        for k in seglist:
            if k != 0:
                if k==124:
                    pdb.set_trace()
                try:
                    # TODO use z to determine ss_name
                    name = L.Gs.node[k]['ss_name'][0]
                except:
                    name = L.Gs.node[k]['name']
                #if k in S.indoor.ce.keys():
                #if k in S.L.ce.keys():
                # nom du sous-segment
                #    indss = S.L.ce[k][0]
                #    name  = S.L.sl.di[indss]
                #    print name
                #"else:
                # nom du segment
                #    name = S.L.Gs.node[k]['name']
                the = theta[i]
                # todo : comparison multiwall th=0 th=variable
                # comparison mesurement
                #the = 0

                i   = i + 1
                #
                # Loss0 du slab
                #
                lko,lkp  = L.sl[name].losst(fGHz,the)
                do , dp  = L.sl[name].excess_grdelay(theta=the)
                edo = edo - np.mean(do)
                edp = edp - np.mean(dp)
    #           print lko
    #           print lkp
                Lo   = Lo + lko[0]
                Lp   = Lp + lkp[0]
        Lwo = np.hstack((Lwo,Lo))
        Lwp = np.hstack((Lwp,Lp))
        Edo = np.hstack((Edo,edo))
        Edp = np.hstack((Edp,edp))

    return(Lwo,Lwp,Edo,Edp)

def Loss_mur_the(S,pi,f,p):
    """

    """
    seglist,theta = S.L.angleonlinkold(p,pi)
    return(seglist,theta)

def Loss0(S,rx,ry,f,p):
    """ Calculate Loss through Layers theta=0 deg

    Parameters
    ----------

    S : Simulation object
    rx:
    ry:
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
            L = 0
            pxy = np.array([x,y])
            seglist,theta = S.L.angleonlinkold(p,pxy)
            for k in seglist:
                name = S.L.name[k]
                lk = S.sl[name].loss0(f)
                L  = L + lk[0]
            Lw[i,j] = L
            j = j+1
        i = i+1
    return(Lw)

def Loss_diff(u):
    """ Calculate Path Loss of the diffraction
    """
    if u < -0.7:
        Ld = 0
    elif u > 1.5:
        Ld = 13 + 20*np.log10(u)
    else:
        Ld = 6.9 + 20*np.log10(np.sqrt((u-0.1)**2+1)+u-0.1)

    return(Ld)

def calnu(h,d1,d2,fGHz):
    """ Calculate the diffraction Fresnel parameter

    Parameters
    ----------

    h  : signed height w.r.t LOS (meter)
    d1 : distance 1 (meter)
    d2 : distance 2 (meter)
    fGHz  : frequency GHz

    Notes
    -----

    .. math::   \\nu = h \\sqrt{\\frac{2}{\\lambda} \\frac{d_1+d_2}{d_1 d_2}}

    """

    ld  = 0.3/fGHz
    nu  = h*np.sqrt(2*(d1+d2)/(ld*d1*d2))

    return(nu)

def Carretosegment(Mob):
    """ define 4 segment using the position of the rectangle

    """
    PC = Mob.position()

    seg1 = (PC[0][0],PC[0][1],PC[1][0],PC[1][1])
    seg2 = (PC[1][0],PC[1][1],PC[2][0],PC[2][1])
    seg3 = (PC[2][0],PC[2][1],PC[3][0],PC[3][1])
    seg4 = (PC[3][0],PC[3][1],PC[0][0],PC[0][1])

    return(seg1,seg2,seg3,seg4)  

def Intersection(x1,y1,x2,y2,x3,y3,x4,y4):
    """ check intersecton  of 2 segments 

    Parameters
    ----------

    segment 1 : 
        ({x1,y1),(x2,y2)}
    segment 2 : 
        ({x3,y3),(x4,y4)}
    Returns
    -------

    Notes
    -----
    This function is implemented better in GeomUtil intersect
    """
    if max(x1,x2) <  min(x3,x4):
        return False
    else:
    #    if x3 == x4:
            if abs(x3-x4) < 1e-9:
            #  y = A1*x+B1  y = x3
                A1 = (y2-y1)/(x2-x1)
                B1 = (x2*y1-x1*y2)/(x2-x1)
                Xa = x3
                Ya = x3*A1+B1
                if Ya > max(y3,y4) or Ya < min(y3,y4) or Ya > max(y1,y2) or Ya < min(y1,y2):
                    return False
                else:
                    return(Xa,Ya)
    #    elif y3 == y4:
            elif abs(y3-y4) < 1e-9:
            # y = A1*x+B1  x = y3
                A1 = (y2-y1)/(x2-x1)
                B1 = (x2*y1-x1*y2)/(x2-x1)
                Ya = y3
                Xa = (y3-B1)/A1
                if Xa > max(x3,x4) or Xa < min(x3,x4) or Xa > max(x1,x2) or Xa < min(x1,x2):
                    return False
                else:
                    return(Xa,Ya)  
            else:
        #  
        # y = A1*x+B1     y = A2*x+B2
        #
                A1 = (y2-y1)/(x2-x1)
                A2 = (y3-y4)/(x3-x4)
                B1 = (x2*y1-x1*y2)/(x2-x1)
                B2 = (x4*y3-x3*y4)/(x4-x3)

                if A1 == A2:
                    return False

            # intersect point (Xa,Ya)
                else:

                    Xa = (B2-B1)/(A1-A2)
                    Ya = Xa*A1+B1

                    if ((Xa < max(min(x1,x2), min(x3,x4))) or  (Xa > min(max(x1,x2), max(x3,x4)))):
                        return False
          
                    else:
                        return(Xa,Ya)

def Dis(x1,y1,x2,y2,x3,y3):
    """ Distance between a point and a line

    Parameters
    ----------
        point (x1,y1)
        line  {(x2,y2),(x3,y3)}

    """
    A = (y3-y2)/(x3-x2)
    B = (x3*y2-x2*y3)/(x3-x2)

    # h: distance
    # ax+by+c=0   h = |a*xo+b*yo+c|/np.sqrt(a*a+b*b)
    #

    h = abs(A*x1-y1+B)/np.sqrt(A*A+1)

    return h

def Interline(x1,y1,x2,y2,Obstacle):
    """
    Parameters
    ----------
    x1
    y1
    x2
    y2
    Obstacle

    Returns
    -------
    SS
    """
    SS =  np.array([])
    SEG = Carretosegment(Obstacle)
    for l in range(4):
        Pint = Intersection(x1,y1,x2,y2,SEG[l][0],SEG[l][1],SEG[l][2],SEG[l][3])
        if Pint ==False:
            print 'No intersection'
        else:
            SS = hstack((SS,l))
    return(SS)

def Loss_obstacle(SS,x1,y1,x2,y2,Obstacle):
    """
    Parameters
    ----------
    SS

        
    """
    LD =  np.array([])
    N = [0,1,2,3]
    K = [val for val in N if val not in SS]
    SEG = Carretosegment(Obstacle)
    height = Obstacle.height
    if len(SS)==0:
        LM = 0
    elif len(SS)==1:
        LM = 0
    else:
        if SS[0]+2 ==SS[1]:
            for n in range(2):
                P1  = [SEG[K[n]][0],SEG[K[n]][1]]
                P2  = [SEG[K[n]][2],SEG[K[n]][3]]
                h1  = Dis(P1[0],P1[1],x1,y1,x2,y2)
                h2  = Dis(P2[0],P2[1],x1,y1,x2,y2)
                d11 = np.sqrt((P1[0]-x1)*(P1[0]-x1)+(P1[1]-y1)*(P1[1]-y1)-h1*h1)
                d12 = np.sqrt((P1[0]-x2)*(P1[0]-x2)+(P1[1]-y2)*(P1[1]-y2)-h1*h1)
                d21 = np.sqrt((P2[0]-x1)*(P2[0]-x1)+(P2[1]-y1)*(P2[1]-y1)-h2*h2)
                d22 = np.sqrt((P2[0]-x2)*(P2[0]-x2)+(P2[1]-y2)*(P2[1]-y2)-h2*h2)
                v1 = calnu(h1,d11,d12,f)
                v2 = calnu(h2,d21,d22,f)
                v = max(v1,v2)
                Ld1 = Loss_diff(v)
              
                if v1 ==v:
                    p0 = P1
                    pm = P2
                else:
                    p0 = P2
                    pm = P1
                h0  = Dis(p0[0],p0[1],x1,y1,x2,y2)
                hm  = Dis(pm[0],pm[1],x1,y1,x2,y2)
                d01 = np.sqrt((p0[0]-x1)*(p0[0]-x1)+(p0[1]-y1)*(p0[1]-y1)-h0*h0)
                dm1 = np.sqrt((pm[0]-x1)*(pm[0]-x1)+(pm[1]-y1)*(pm[1]-y1)-hm*hm)
                if d01 < dm1:
                    # new line: p0 Rx
                    pt = [x2,y2]
                else:
                    # new line: p0 Tx
                    pt = [x1,y1]
          
                d0t =  np.sqrt((p0[0]-pt[0])*(p0[0]-pt[0])+(p0[1]-pt[1])*(p0[1]-pt[1])-h0*h0)
                dmt =  np.sqrt((pm[0]-pt[0])*(pm[0]-pt[0])+(pm[1]-pt[1])*(pm[1]-pt[1])-hm*hm)
                d0m = d0t - dmt
                hm2 = (dmt/d0t)*h0
                if hm2 > hm:
                    Ld2 = 0
                else:
                    hmt = hm - hm2
                    v0t = calnu(hmt,dmt,d0m,f)
                    Ld2 = Loss_diff(v)

                Ld = Ld1 + Ld2
                LD = hstack((LD,Ld))                       
            hh = height - 1.2
            Pint1 = Intersection(x1,y1,x2,y2,SEG[int(SS[0])][0],SEG[int(SS[0])][1],SEG[int(SS[0])][2],SEG[int(SS[0])][3])  
            Pint2 = Intersection(x1,y1,x2,y2,SEG[int(SS[1])][0],SEG[int(SS[1])][1],SEG[int(SS[1])][2],SEG[int(SS[1])][3])
            dis11 = np.sqrt((Pint1[0]-x1)*(Pint1[0]-x1)+(Pint1[1]-y1)*(Pint1[1]-y1)+hh*hh)
            dis12 = np.sqrt((Pint1[0]-x2)*(Pint1[0]-x2)+(Pint1[1]-y2)*(Pint1[1]-y2)+hh*hh)
            dis21 = np.sqrt((Pint2[0]-x1)*(Pint2[0]-x1)+(Pint2[1]-y1)*(Pint2[1]-y1)+hh*hh)
            dis22 = np.sqrt((Pint2[0]-x2)*(Pint2[0]-x2)+(Pint2[1]-y2)*(Pint2[1]-y2)+hh*hh)
            vh1 = calnu(hh,dis11,dis12,f)
            vh2 = calnu(hh,dis21,dis22,f)
            vh = max(vh1,vh2)
            Ldh = Loss_diff(vh)
            LD = hstack((LD,Ldh))
        #    consider the path througth from haut of the furniture  
        #    LM = -10*log10(10**((-LD[0]/10))+10**((-LD[1]/10)))
            LM = -10*log10(10**((-LD[0]/10))+10**((-LD[1]/10))+10**((-LD[2]/10)))
        else:  
            p1 = [SEG[int(SS[0])][0],SEG[int(SS[0])][1]]
            p2 = [SEG[int(SS[0])][2],SEG[int(SS[0])][3]]
            p3 = [SEG[int(SS[1])][0],SEG[int(SS[1])][1]]
            p4 = [SEG[int(SS[1])][2],SEG[int(SS[1])][3]]
            if p1==p4:
                Po = p1
            if p2==p3:
                Po = p2  
        #    Po = [val for val in SEG[S[0]] if val in SEG[S[1]]]
            h = Dis(Po[0],Po[1],x1,y1,x2,y2)
            d1 = np.sqrt((Po[0]-x1)*(Po[0]-x1)+(Po[1]-y1)*(Po[1]-y1)-h*h)                  
            d2 = np.sqrt((Po[0]-x2)*(Po[0]-x2)+(Po[1]-y2)*(Po[1]-y2)-h*h)
            v = calnu(h,d1,d2,f)
            Ld = Loss_diff(v)
            LD = hstack((LD,Ld))

            p1 = [SEG[K[0]][0],SEG[K[0]][1]]
            p2 = [SEG[K[0]][2],SEG[K[0]][3]]
            p3 = [SEG[K[1]][0],SEG[K[1]][1]]
            p4 = [SEG[K[1]][2],SEG[K[1]][3]]

            if p2 == p3:
                pa = p2
                pb = p1
                pc = p4
            if p1 == p4:
                pa = p1
                pb = p2
                pc = p3
  
            ha = Dis(pa[0],pa[1],x1,y1,x2,y2)
            hb = Dis(pb[0],pb[1],x1,y1,x2,y2)
            hc = Dis(pc[0],pc[1],x1,y1,x2,y2)
          
            da1 = np.sqrt((pa[0]-x1)*(pa[0]-x1)+(pa[1]-y1)*(pa[1]-y1)-ha*ha)                  
            da2 = np.sqrt((pa[0]-x2)*(pa[0]-x2)+(pa[1]-y2)*(pa[1]-y2)-ha*ha)

            db1 = np.sqrt((pb[0]-x1)*(pb[0]-x1)+(pb[1]-y1)*(pb[1]-y1)-hb*hb)
            db2 = np.sqrt((pb[0]-x2)*(pb[0]-x2)+(pb[1]-y2)*(pb[1]-y2)-hb*hb)

            dc1 = np.sqrt((pc[0]-x1)*(pc[0]-x1)+(pc[1]-y1)*(pc[1]-y1)-hc*hc)                  
            dc2 = np.sqrt((pc[0]-x2)*(pc[0]-x2)+(pc[1]-y2)*(pc[1]-y2)-hc*hc)

            va = calnu(ha,da1,da2,f)
            vb = calnu(hb,db1,db2,f)
            vc = calnu(hc,dc1,dc2,f)
            vmax = max(va,vb,vc)
            Ld1 = Loss_diff(vmax)

            if vmax == va:
                p0 = pa
                pm = pb
                pn = pc
                h0 = ha
                hm = hb
                hn = hc
            elif vmax ==vb:
                p0 = pb
                pm = pa
                pn = pc
                h0 = hb
                hm = ha
                hn = hc
            else:
                p0 = pc
                pm = pa
                pn = pb
                h0 = hc
                hm = ha
                hn = hb
            d01 = np.sqrt((p0[0]-x1)*(p0[0]-x1)+(p0[1]-y1)*(p0[1]-y1)-h0*h0)
            d02 = np.sqrt((p0[0]-x2)*(p0[0]-x2)+(p0[1]-y2)*(p0[1]-y2)-h0*h0)
            dm1 = np.sqrt((pm[0]-x1)*(pm[0]-x1)+(pm[1]-y1)*(pm[1]-y1)-hm*hm)
            dm2 = np.sqrt((pm[0]-x2)*(pm[0]-x2)+(pm[1]-y2)*(pm[1]-y2)-hm*hm)
            dn1 = np.sqrt((pn[0]-x1)*(pn[0]-x1)+(pn[1]-y1)*(pn[1]-y1)-hn*hn)
            dn2 = np.sqrt((pn[0]-x2)*(pn[0]-x2)+(pn[1]-y2)*(pn[1]-y2)-hn*hn)
            nl  = find(np.hstack((dm1,dn1))<d01)
            nr  = find(np.hstack((dm2,dn2))<d02)
            if len(nl) == 0:
                Ldl = 0
            else:
                pt = [x1,y1]
                VL =  np.array([])
                for i in range(len(nl)):
                    hh = (np.hstack((dm1,dn1))[nl[i]]/d01)*h0
                    hreal = np.hstack((hm,hn))[nl[i]]
                    if hreal < hh:
                        Ldl = 0
                        VL = hstack((VL,Ldl))
                    else:
                        hmt = hreal - hh
                        pp = [np.hstack((pm,pn))[2*nl[i]],np.hstack((pm,pn))[2*nl[i]+1]]
                        ldt = np.sqrt((pp[0]-x1)*(pp[0]-x1)+(pp[1]-y1)*(pp[1]-y1)-hreal*hreal)               
                        ld0 = d01 - ldt
                        vl = calnu(hmt,ld0,ldt,f)
                        Ldl = Loss_diff(vl)
                        VL = hstack((VL,Ldl))   
                Ldl = max(VL)

            if len(nr) == 0:
                Ldr = 0          
            else:
                pt = [x2,y2]
                VR =  np.array([])
                for i in range(len(nr)):
                    hh = (hstack((dm2,dn2))[nr[i]]/d02)*h0
                    hreal = hstack((hm,hn))[nr[i]]
                    if hreal <hh:
                        Ldr = 0
                        VR = hstack((VR,Ldr))
                    else:
                        hmt = hreal - hh
                        pp = [hstack((pm,pn))[2*nr[i]],hstack((pm,pn))[2*nr[i]+1]]              
                        rdt = np.sqrt((pp[0]-x2)*(pp[0]-x2)+(pp[1]-y2)*(pp[1]-y2)-hreal*hreal)
                        rd0 = d02 -rdt
                        vr = calnu(hmt,rd0,rdt,f)
                        Ldr = Loss_diff(vr)
                        VR = hstack((VR,Ldr))
                Ldr = max(VR)
      
            Ld = Ld1+Ldl+Ldr
            LD =hstack((LD,Ld))
          
            hh = height - 1.2
            Pint1 = Intersection(x1,y1,x2,y2,SEG[int(SS[0])][0],SEG[int(SS[0])][1],SEG[int(SS[0])][2],SEG[int(SS[0])][3])  
            Pint2 = Intersection(x1,y1,x2,y2,SEG[int(SS[1])][0],SEG[int(SS[1])][1],SEG[int(SS[1])][2],SEG[int(SS[1])][3])
            dis11 = np.sqrt((Pint1[0]-x1)*(Pint1[0]-x1)+(Pint1[1]-y1)*(Pint1[1]-y1)+hh*hh)
            dis12 = np.sqrt((Pint1[0]-x2)*(Pint1[0]-x2)+(Pint1[1]-y2)*(Pint1[1]-y2)+hh*hh)
            dis21 = np.sqrt((Pint2[0]-x1)*(Pint2[0]-x1)+(Pint2[1]-y1)*(Pint2[1]-y1)+hh*hh)
            dis22 = np.sqrt((Pint2[0]-x2)*(Pint2[0]-x2)+(Pint2[1]-y2)*(Pint2[1]-y2)+hh*hh)
            vh1 = calnu(hh,dis11,dis12,f)
            vh2 = calnu(hh,dis21,dis22,f)
            vh = max(vh1,vh2)
            Ldh = Loss_diff(vh)
            LD = hstack((LD,Ldh))
          
            LM = -10*log10(10**((-LD[0]/10))+10**((-LD[1]/10))+10**((-LD[2]/10)))  
            #LM = -10*log10(10**((-LD[0]/10))+10**((-LD[1]/10)))
          
    return(LM)                              

def Loss_2obstacle(x1,y1,x2,y2,Obstacle1,Obstacle2):
    """ Yu Lei function

    Parameters
    ----------

    x1
    y1
    x2
    y2
    Obstacle1
    Obstacle2

    """

    S1 = Interline(x1,y1,x2,y2,Obstacle1)
    S2 = Interline(x1,y1,x2,y2,Obstacle2)
    N = [0,1,2,3]
    K = [val for val in N if val not in S1]
    SEG = Carretosegment(Obstacle1)
    height1 = Obstacle1.height
    height2 = Obstacle2.height
    if S1[0]+2==S1[1]:
        p1 = [SEG[K[0]][0],SEG[K[0]][1]]
        p2 = [SEG[K[0]][2],SEG[K[0]][3]]
        p3 = [SEG[K[1]][0],SEG[K[1]][1]]
        p4 = [SEG[K[1]][2],SEG[K[1]][3]]
        Ga = [p1,p2]
        Gb = [p3,p4]
    else:
        p1 = [SEG[int(S1[0])][0],SEG[int(S1[0])][1]]
        p2 = [SEG[int(S1[0])][2],SEG[int(S1[0])][3]]
        p3 = [SEG[int(S1[1])][0],SEG[int(S1[1])][1]]
        p4 = [SEG[int(S1[1])][2],SEG[int(S1[1])][3]]

        if p1==p4:
            po = p1
        if p2==p3:
            po = p2
        Ga = [po]

        Pa = [SEG[K[0]][0],SEG[K[0]][1]]
        Pb = [SEG[K[0]][2],SEG[K[0]][3]]
        Pc = [SEG[K[1]][0],SEG[K[1]][1]]
        Pd = [SEG[K[1]][2],SEG[K[1]][3]]

        if Pa==Pd:
            pp = Pa
            Gb = [Pa,Pb,Pc]

        if Pb==Pc:
            pp = Pb
            Gb = [Pa,Pc,Pd]

    SEG = Carretosegment(Obstacle2)
    K = [val for val in N if val not in S2]

    if S2[0]+2==S2[1]:
        p1 = [SEG[K[0]][0],SEG[K[0]][1]]
        p2 = [SEG[K[0]][2],SEG[K[0]][3]]
        p3 = [SEG[K[1]][0],SEG[K[1]][1]]
        p4 = [SEG[K[1]][2],SEG[K[1]][3]]

        Ga2 = [p1,p2]
        Gb2 = [p3,p4]
    else: 
        p1 = [SEG[int(S2[0])][0],SEG[int(S2[0])][1]]
        p2 = [SEG[int(S2[0])][2],SEG[int(S2[0])][3]]
        p3 = [SEG[int(S2[1])][0],SEG[int(S2[1])][1]]
        p4 = [SEG[int(S2[1])][2],SEG[int(S2[1])][3]]

        if p1==p4:
                po = p1
        if p2==p3:
                po = p2
        Ga2 = [po]

        Pa = [SEG[K[0]][0],SEG[K[0]][1]]
        Pb = [SEG[K[0]][2],SEG[K[0]][3]]
        Pc = [SEG[K[1]][0],SEG[K[1]][1]]
        Pd = [SEG[K[1]][2],SEG[K[1]][3]]

        if Pa==Pd:
            pp = Pa
            Gb2 = [Pa,Pb,Pc]

        if Pb==Pc:
            pp = Pb
            Gb2 = [Pa,Pc,Pd]

    if Intersection(xt,yt,xr,yr,Ga[0][0],Ga[0][1],Ga2[0][0],Ga2[0][1]) == False:
        Gupper  = vstack((Ga,Ga2))
        Gbottom = vstack((Gb,Gb2))
    else:
        Gupper  = vstack((Ga,Gb2))
        Gbottom = vstack((Ga2,Gb))

    h_upper = Dis(Gupper[:,0],Gupper[:,1],x1,y1,x2,y2)
    h_bottom = Dis(Gbottom[:,0],Gbottom[:,1],x1,y1,x2,y2)
    d1_upper = np.sqrt((Gupper[:,0]-x1)*(Gupper[:,0]-x1)+(Gupper[:,1]-y1)*(Gupper[:,1]-y1)-h_upper[:]*h_upper[:])
    d2_upper = np.sqrt((Gupper[:,0]-x2)*(Gupper[:,0]-x2)+(Gupper[:,1]-y2)*(Gupper[:,1]-y2)-h_upper[:]*h_upper[:])
    v_upper = calnu(h_upper[:],d1_upper[:],d2_upper[:],f)
    v_max1=max(v_upper)
    nn = find(v_upper==v_max1)
    Ld = Loss_diff(v_max1)

    left = find(d1_upper<d1_upper[nn])
    right = find(d1_upper>d1_upper[nn])

    if len(left)!=0:
        Hstd = h_upper[nn][0]*d1_upper[left][:]/d1_upper[nn]
        hreal = h_upper[left]
        nl = find((Hstd-hreal)<0)
        if len(nl)!=0:
            height = (hreal - Hstd)[nl]
            dl1 = np.sqrt((Gupper[left][nl][:,0]-xt)*(Gupper[left][nl][:,0]-xt)+(Gupper[left][nl][:,1]-yt)*(Gupper[left][nl][:,1]-yt)-hreal[nl][:]*hreal[nl][:])
            dl2 = d1_upper[nn] - dl1
            vl = calnu(height[:],dl1[:],dl2[:],f)
            Ll = Loss_diff(max(vl))
        else:
            Ll = 0

    else:
         Ll = 0

    if len(right)!=0:
        Hstd = h_upper[nn][0]*d2_upper[right][:]/d2_upper[nn]
        hreal = h_upper[right]
        nr = find((Hstd-hreal)<0)
        if len(nr)!=0:
            height = (hreal - Hstd)[nr]
            dr2 = np.sqrt((Gupper[right][nr][:,0]-xr)*(Gupper[right][nr][:,0]-xr)+(Gupper[right][nr][:,1]-yr)*(Gupper[right][nr][:,1]-yr)-hreal[nr][:]*hreal[nr][:])  
            dr1 = d2_upper[nn] - dr2
            vr = calnu(height[:],dr1[:],dr2[:],f)
            Lr = Loss_diff(max(vr))
        else:
            Lr = 0
    else:
        Lr = 0

    Lupper = Ld + Ll + Lr
#    Lupper = -10*log10(10**((-Ld/10))+10**((-Ll/10))+10**((-Lr/10)))


    d1_bottom = np.sqrt((Gbottom[:,0]-x1)*(Gbottom[:,0]-x1)+(Gbottom[:,1]-y1)*(Gbottom[:,1]-y1)-h_bottom[:]*h_bottom[:])
    d2_bottom = np.sqrt((Gbottom[:,0]-x2)*(Gbottom[:,0]-x2)+(Gbottom[:,1]-y2)*(Gbottom[:,1]-y2)-h_bottom[:]*h_bottom[:])
    v_bottom = calnu(h_bottom[:],d1_bottom[:],d2_bottom[:],f)
    v_max2=max(v_bottom)
    nn2 = find(v_bottom==v_max2)
    Ld2 = Loss_diff(v_max2)

    left2= find(d1_bottom<d1_bottom[nn2])
    right2= find(d1_bottom>d1_bottom[nn2])

    if len(left2)!=0:
        Hstd = h_bottom[nn2][0]*d1_bottom[left2][:]/d1_bottom[nn2]
        hreal = h_bottom[left2]
        nl = find((Hstd-hreal)<0)
        if len(nl)!=0:
            height = (hreal - Hstd)[nl]
            dl1 = np.sqrt((Gbottom[left2][nl][:,0]-xt)*(Gbottom[left2][nl][:,0]-xt)+(Gbottom[left2][nl][:,1]-yt)*(Gbottom[left2][nl][:,1]-yt)-hreal[nl][:]*hreal[nl][:])
            dl2 = d1_bottom[nn2] - dl1
            vl = calnu(height[:],dl1[:],dl2[:],f)
            Ll = Loss_diff(max(vl))
        else:
            Ll = 0

    else:
         Ll = 0

    if len(right2)!=0:
        Hstd = h_bottom[nn2][0]*d2_bottom[right2][:]/d2_bottom[nn2]
        hreal = h_bottom[right2]
        nr = find((Hstd-hreal)<0)
        if len(nr)!=0:
            height = (hreal -Hstd)[nr]
            dr2 = np.sqrt((Gbottom[right2][nr][:,0]-xr)*(Gbottom[right2][nr][:,0]-xr)
                         +(Gbottom[right2][nr][:,1]-yr)*(Gbottom[right2][nr][:,1]-yr)
                          -hreal[nr][:]*hreal[nr][:])
            dr1 = d2_bottom[nn2] - dr2
            vr = calnu(height[:],dr1[:],dr2[:],f)
            Lr = Loss_diff(max(vr))
        else:
            Lr = 0
    else:
        Lr = 0

    Lbottom = Ld2 + Ll + Lr
#    Lbottom = -10*log10(10**((-Ld2/10))+10**((-Ll/10))+10**((-Lr/10)))


    hh1 = height1 - 1.2
    SEG1 = Carretosegment(Obstacle1)
    SEG2 = Carretosegment(Obstacle2)
    Pint1 = Intersection(x1,y1,x2,y2,SEG1[int(S1[0])][0],SEG1[int(S1[0])][1],SEG1[int(S1[0])][2],SEG1[int(S1[0])][3])  
    Pint2 = Intersection(x1,y1,x2,y2,SEG1[int(S1[1])][0],SEG1[int(S1[1])][1],SEG1[int(S1[1])][2],SEG1[int(S1[1])][3])
    #print Pint1
    #print Pint2
    dis11 = np.sqrt((Pint1[0]-x1)*(Pint1[0]-x1)+(Pint1[1]-y1)*(Pint1[1]-y1)+hh1*hh1)
    dis12 = np.sqrt((Pint1[0]-x2)*(Pint1[0]-x2)+(Pint1[1]-y2)*(Pint1[1]-y2)+hh1*hh1)
    dis21 = np.sqrt((Pint2[0]-x1)*(Pint2[0]-x1)+(Pint2[1]-y1)*(Pint2[1]-y1)+hh1*hh1)
    dis22 = np.sqrt((Pint2[0]-x2)*(Pint2[0]-x2)+(Pint2[1]-y2)*(Pint2[1]-y2)+hh1*hh1)
    vh1 = calnu(hh1,dis11,dis12,f)
    vh2 = calnu(hh1,dis21,dis22,f)

    hh2 = height2-1.2
    Pint3 = Intersection(x1,y1,x2,y2,SEG2[int(S2[0])][0],SEG2[int(S2[0])][1],SEG2[int(S2[0])][2],SEG2[int(S2[0])][3])  
    Pint4 = Intersection(x1,y1,x2,y2,SEG2[int(S2[1])][0],SEG2[int(S2[1])][1],SEG2[int(S2[1])][2],SEG2[int(S2[1])][3])
    dis31 = np.sqrt((Pint3[0]-x1)*(Pint3[0]-x1)+(Pint3[1]-y1)*(Pint3[1]-y1)+hh2*hh2)
    dis32 = np.sqrt((Pint3[0]-x2)*(Pint3[0]-x2)+(Pint3[1]-y2)*(Pint3[1]-y2)+hh2*hh2)
    dis41 = np.sqrt((Pint4[0]-x1)*(Pint4[0]-x1)+(Pint4[1]-y1)*(Pint4[1]-y1)+hh2*hh2)
    dis42 = np.sqrt((Pint4[0]-x2)*(Pint4[0]-x2)+(Pint4[1]-y2)*(Pint4[1]-y2)+hh2*hh2)
    vh3 = calnu(hh2,dis31,dis32,f)
    vh4 = calnu(hh2,dis41,dis42,f)


    vh = max(vh1,vh2,vh3,vh4)
    Lheight = Loss_diff(vh)


    Ltotal = -10*log10(10**((-Lupper/10))+10**((-Lbottom/10))+10**((-Lheight/10)))
#    Ltotal = -10*log10(10**((-Lupper/10))+10**((-Lbottom/10)))

    return(Ltotal)

def showfurniture(ax):
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


def visuPts(S,nu,nd,Pts,Values,fig=[],sp=[],vmin=0,vmax=-1,label=' ',tit='',size=25,colbar=True,xticks=False):
    """
    visuPt  : Visualization of values a given points

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
    """ 
        Plot the cumulative density function
    """
    rcParams['legend.fontsize']=20
    rcParams['font.size']=20

    x  = sort(x)
    n  = len(x)
    x2 = repeat(x, 2)
    y2 = hstack([0.0, repeat(arange(1,n) / float(n), 2), 1.0])
    plt.plot(x2,y2,colsym,label=lab,linewidth=lw)
#    plt.semilogx(x2,y2,colsym,label=lab,linewidth=lw)
    plt.grid('on')
    legend(loc=2)
    xlabel('Ranging Error[m]')
    ylabel('Cumulative Probability')

if __name__=="__main__":
    doctest.testmod()
