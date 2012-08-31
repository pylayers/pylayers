# -*- coding:Utf-8 -*-
#from numpy import *
import doctest
import numpy as np 
from scipy import io
import matplotlib.pylab as plt 
import pylayers.simul.simulem 
import pylayers.measures.mesuwb

def LOSS_furniture(Tx,Rx,furn):
    """
      Yu Lei 
      Pas utilisé
    """
    for i in range(1,5):
        rx = Rx[i]
        for j in range(1,377):
            tx = Tx[j]
            x1 = tx[0]
            y1 = tx[1]
            x2 = rx[0]  
            y2 = rx[1]          
            T = furn
            position = T.position()

def PL0(fGHz,GtdB=0,GrdB=0):
    """  Path Loss at frequency f @ 1m 

    Parameters
    ----------
    fGHz:
        frequency GHz
    GtdB:
        transmitting antenna gain dB (default 0 dB)
    GrdB:   
        receiving antenna gain dB (default 0 dB)

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
    """
    Dgrid_points(points,Px):
    points : Np x 2 array
    Px     : 2 x 1  array
    """
  
    Dx = points[:,0]-Px[0]
    Dy = points[:,1]-Px[1]
    D  = np.sqrt( Dx*Dx + Dy*Dy )

    return(D)

def Dgrid_zone(zone,Px):
    """
    Dgrid_zone(zone,Px):

    A zone is a quadrilateral zone.

    Zone dictionnary :
        xmin xmax Nx
        ymin ymax Ny

    Build the distance matrix between Tx and the zone points

    """

    rx = np.linspace(zone['xmin'],zone['xmax'],zone['Nx'])
    ry = np.linspace(zone['ymin'],zone['ymax'],zone['Ny'])
    R_x = np.outer(np.ones(len(ry)),rx)
    R_y = np.outer(ry,np.ones(len(rx)))
    Dx = R_x - Px[0]
    Dy = R_y - Px[1]
    D  = np.sqrt(Dx*Dx+Dy*Dy)
    return (D)

def OneSlopeMdl(D,n,fGHz):
    """
    Parameters
    ----------
    n 
        path loss exponent
    D 
        Distance array
    fGHz 
        frequency  GHz

    """
    PL = PL0(fGHz)+10*n*np.log10(D)
    return(PL)

def Loss0_v2(S,Pts,f,p):
    """
    Parameters
    ----------
    S
        Simulation object
    Pts
        (Np x 2) array
    f
        frequency (GHz)
    p
        source points

    Examples
    --------

    >>> import matplotlib.pyplot as plt 
    >>> import Simul
    >>> import MesCEA
    >>> import MultiWall
    >>> S = Simul.Simul()
    >>> S.layout('sircut.str','simul8.mat','simul8.slab')
    >>> fGHz = 4 
    >>> Tx,Rx = MesCEA.ptSiradel()
    >>> Lwo,Lwp = MultiWall.Loss0_v2(S,Tx,fGHz,Rx[1,0:2])
    >>> fig,ax = S.L.showGs()
    >>> tit = plt.title('test Loss0_v2')
    >>> sc2 = ax.scatter(Rx[1,0],Rx[1,1],s=20,marker='x',c='k')
    >>> sc1 = ax.scatter(Tx[:,0],Tx[:,1],s=Lwo,c=Lwo,linewidth=0)
    >>> plt.show()

    .. plot::

	 import matplotlib.pyplot as plt 
         import Simul
         import MesCEA
    	 import MultiWall
         S = Simul.Simul()
         S.layout('sircut.str','simul8.mat','simul8.slab')
         fGHz = 4 
         Tx,Rx = MesCEA.ptSiradel()
         Lwo,Lwp = MultiWall.Loss0_v2(S,Tx,fGHz,Rx[1,0:2])
         fig,ax = S.L.showGs()
         tit = plt.title('test Loss0_v2')
         sc2 = ax.scatter(Rx[1,0],Rx[1,1],s=20,marker='x',c='k')
         sc1 = ax.scatter(Tx[:,0],Tx[:,1],s=Lwo,c=Lwo,linewidth=0)
         plt.show()


    """
    N   = np.shape(Pts)[0]
    Lwo = np.array([])
    Lwp = np.array([])
    for i in range(N):
        Lo = 0.0
        Lp = 0.0
        pi = Pts[i,:]
        seglist,theta = S.L.angleonlink(p,pi)
        i = 0
        for k in seglist:
            try:
                name = S.L.Gs.node[k]['ss_name']
            except:
                name = S.L.Gs.node[k]['name']
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
            # idea paper : comparison multiwall th=0 th=variable
            # comparison mesurement
            #the = 0

            i   = i + 1
            #
            # Loss0 du slab
            #
            lko,lkp  = S.L.sl[name].losst(f,the)
#            print lko
#            print lkp
            Lo   = Lo + lko[0]
            Lp   = Lp + lkp[0]
        Lwo = np.hstack((Lwo,Lo))
        Lwp = np.hstack((Lwp,Lp))

    return(Lwo,Lwp)  
  

def Loss0_v2_separe(S,pi,f,p):
    """
    """
    # for calibrate the loss multiwall
    lwo   = np.array([])
    lwp   = np.array([])
    Theta = np.array([])
    seglist,theta = S.L.angleonlink(p,pi)
    i = 0
    for k in seglist:
        if k in S.L.ce.keys():
            indss = S.L.ce[k][0]
            name = S.L.sl.di[indss]
        else:
            name =S.L. name[k]
        the = theta[i]
        i = i + 1

        lko,lkp = S.sl[name].losst(f,the)
        if name == 'PARTITION':
            lwo = hstack((lwo,lko))
            lwp = hstack((lwp,lkp))
            Theta = hstack((Theta,the))
        else:
            lwo = lwo
            lwp = lwp
            Theta = Theta
    return(seglist,lwo,lwp,Theta)

def Loss_mur_the(S,pi,f,p):
    """

    """
    seglist,theta = S.L.angleonlink(p,pi)
    return(seglist,theta)

def Loss0(S,rx,ry,f,p):
    """ Calculate Loss through Layers theta=0 deg
    Parameters
    ----------
    S
        Simulation 
    f
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
            seglist,theta = S.L.angleonlink(p,pxy)
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

def Diffraction_parameter(h,d1,d2,f):
    """ Calculate the diffraction parameter

    Parameters
    ----------
    h 
        height 
    d1
        distance 1
    d2
        distance 2
    fGHz
        frequency GHz

    Notes
    -----
    .. math::   \\nu = h \\sqrt{2\\frac{d_1+d_2}{\\lambda d_1d_2}}
    """
    ld  = 0.3/f
    nu  = h*np.sqrt(2*(d1+d2)/(ld*d1*d2))

    return(nu)

def Carretosegment(Mob):
    """
    define 4 segment using the position of the rectangle
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
                v1 =Diffraction_parameter(h1,d11,d12,f)
                v2 =Diffraction_parameter(h2,d21,d22,f)
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
                    v0t = Diffraction_parameter(hmt,dmt,d0m,f)
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
            vh1 = Diffraction_parameter(hh,dis11,dis12,f)
            vh2 = Diffraction_parameter(hh,dis21,dis22,f)
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
            v = Diffraction_parameter(h,d1,d2,f)
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

            va = Diffraction_parameter(ha,da1,da2,f)
            vb = Diffraction_parameter(hb,db1,db2,f)
            vc = Diffraction_parameter(hc,dc1,dc2,f)
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
                        vl = Diffraction_parameter(hmt,ld0,ldt,f)
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
                        vr = Diffraction_parameter(hmt,rd0,rdt,f)
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
            vh1 = Diffraction_parameter(hh,dis11,dis12,f)
            vh2 = Diffraction_parameter(hh,dis21,dis22,f)
            vh = max(vh1,vh2)
            Ldh = Loss_diff(vh)
            LD = hstack((LD,Ldh))
          
            LM = -10*log10(10**((-LD[0]/10))+10**((-LD[1]/10))+10**((-LD[2]/10)))  
            #LM = -10*log10(10**((-LD[0]/10))+10**((-LD[1]/10)))
          
    return(LM)                              

def Loss_2obstacle(x1,y1,x2,y2,Obstacle1,Obstacle2):
    """ Yu Lei function 
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
    v_upper = Diffraction_parameter(h_upper[:],d1_upper[:],d2_upper[:],f)
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
            vl = Diffraction_parameter(height[:],dl1[:],dl2[:],f)
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
            vr = Diffraction_parameter(height[:],dr1[:],dr2[:],f)
            Lr = Loss_diff(max(vr))
        else:
            Lr = 0
    else:
        Lr = 0      

    Lupper = Ld + Ll + Lr
#    Lupper = -10*log10(10**((-Ld/10))+10**((-Ll/10))+10**((-Lr/10)))

              
    d1_bottom = np.sqrt((Gbottom[:,0]-x1)*(Gbottom[:,0]-x1)+(Gbottom[:,1]-y1)*(Gbottom[:,1]-y1)-h_bottom[:]*h_bottom[:])
    d2_bottom = np.sqrt((Gbottom[:,0]-x2)*(Gbottom[:,0]-x2)+(Gbottom[:,1]-y2)*(Gbottom[:,1]-y2)-h_bottom[:]*h_bottom[:])
    v_bottom = Diffraction_parameter(h_bottom[:],d1_bottom[:],d2_bottom[:],f)
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
            vl = Diffraction_parameter(height[:],dl1[:],dl2[:],f)
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
            dr2 = np.sqrt((Gbottom[right2][nr][:,0]-xr)*(Gbottom[right2][nr][:,0]-xr)+(Gbottom[right2][nr][:,1]-yr)*(Gbottom[right2][nr][:,1]-yr)-hreal[nr][:]*hreal[nr][:])  
            dr1 = d2_bottom[nn2] - dr2
            vr = Diffraction_parameter(height[:],dr1[:],dr2[:],f)
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
    vh1 = Diffraction_parameter(hh1,dis11,dis12,f)
    vh2 = Diffraction_parameter(hh1,dis21,dis22,f)
  
    hh2 = height2-1.2
    Pint3 = Intersection(x1,y1,x2,y2,SEG2[int(S2[0])][0],SEG2[int(S2[0])][1],SEG2[int(S2[0])][2],SEG2[int(S2[0])][3])  
    Pint4 = Intersection(x1,y1,x2,y2,SEG2[int(S2[1])][0],SEG2[int(S2[1])][1],SEG2[int(S2[1])][2],SEG2[int(S2[1])][3])
    dis31 = np.sqrt((Pint3[0]-x1)*(Pint3[0]-x1)+(Pint3[1]-y1)*(Pint3[1]-y1)+hh2*hh2)
    dis32 = np.sqrt((Pint3[0]-x2)*(Pint3[0]-x2)+(Pint3[1]-y2)*(Pint3[1]-y2)+hh2*hh2)
    dis41 = np.sqrt((Pint4[0]-x1)*(Pint4[0]-x1)+(Pint4[1]-y1)*(Pint4[1]-y1)+hh2*hh2)
    dis42 = np.sqrt((Pint4[0]-x2)*(Pint4[0]-x2)+(Pint4[1]-y2)*(Pint4[1]-y2)+hh2*hh2)
    vh3 = Diffraction_parameter(hh2,dis31,dis32,f)
    vh4 = Diffraction_parameter(hh2,dis41,dis42,f)
  

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
    # Siradel Rooms annotation
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
#
# Siradel environment materials
#
# Sigint Values after calibration with S2 (3.5GHz)
#
#     WALLS  (PARTITION)   8    - 0.2 j      0.1 m
#     WOODEN_DOOR          2.5  - 0.03 j     0.1 m
#     INSULATION           2    - 0.5 j      0.2 m
#     PLASTERBOARD_14CM    2.5  - 0.3 j      0.14
#     PLASTERBOARD_7CM     2.5  - 0.3 j      0.07
#     CONCRETE_15CM3D      5.5  - 0.25 j     0.15
#     CONCRETE_20CM3D      5.5  - 0.25 j     0.2
#     CONCRETE_7CM3D       5.5  - 0.25 j     0.07
#     CONCRETE_25CM3D      5.5  - 0.25 j     0.25
#     3D_WINDOW_GLASSES    4.5  - 0.03 j     0.05
#     WOODEN_TABLES_DESK2  2.5  - 0.03 j     0.01
#     FLOOR                4.5  - 0.03 j     0.5
#     TOP_CEILING          4.5  - 0.6 j      0.5
#     CEILING              2.5  - 0.3 j      0.005
#     Metallics            1    - 1e6 j      0.005
#     PLASTIC              4    - 0.4 j      0.1
#
    #R1_A=furniture('Table',array([1.33,5.35]),120,0.75,1.4,0.7,0.03,'WOOD')
    #R1_B1=furniture('Desk',array([4.05,7.5]),90,0.75,1.22,0.8,0.03,'WOOD')
    #R1_B2=furniture('Desk',array([4.05,7.1]),0,0.75,0.85,0.4,0.03,'WOOD')
  
    #R2_A=furniture('Block of Tables',array([0.5,13.5]),0,0.75,3.30,1.2,0.03,'WOOD')
  
    #R6_A=furniture('1.3m long table',array([-5.47,8]),90,0.75,1.30,0.65,0.03,'WOOD')
    #R6_B=furniture('1.6m long table',array([-5.32,6]),90,0.75,1.60,0.8,0.03,'WOOD')
    #R6_C=furniture('0.9m high cupboard',array([-0.7,5]),0,0.9,0.60,0.6,0.01,'WOOD')
#    R6_D=furniture('refrigerator', np.array([-0.7,5.6]),0,0.9,0.60,0.6,0.005,'METAL')
#    R6_E=furniture('Drinking fountain', np.array([-0.59,6.2]),0,1.22,0.4,0.35,0.005,'METAL')
#    R6_F=furniture('1.7m high cupboard', np.array([-0.05,7.7]),90,1.7,0.9,0.42,0.005,'METAL')
#    R6_G=furniture('0.55m high printer', np.array([-5.88,6]),0,0.55,0.53,0.5,0.005,'METAL')
#    R6_HA=furniture('copy machine', np.array([-5.76,5]),0,1.17,1.2,0.6,0.005,'METAL')
#    R6_HB=furniture('copy machine', np.array([-5.46,5]),0,1.32,0.9,0.6,0.005,'METAL')
#    R6_IA=furniture('high printer', np.array([-4.3,4.85]),0,0.5,0.75,0.57,0.005,'METAL')
#    R6_IB=furniture('high printer', np.array([-4.3,4.85]),0,1.06,1.0,0.57,0.005,'METAL')
#  
#    R13_A = furniture('1.33-high cupboard', np.array([-6.17,13.05]),0,1.33,0.64,0.4,0.005,'METAL')
#    R13_B = furniture('2m-high cupboard', np.array([-6,13.5]),0,2,1.6,0.45,0.005,'METAL')
#    R13_C = furniture('1.7-high cupboard', np.array([-5.7,15.3]),90,1.7,0.9,0.42,0.005,'METAL')
#  
#  
#    R12_A1= furniture('1.98m high cupboards', np.array([-7.54,9.88]),0,1.98,1.2,0.42,0.005,'METAL')
#    R12_A2= furniture('1.98m high cupboards', np.array([-6.16,6.2]),90,1.98,1.2,0.42,0.005,'METAL')
#    #R12_B= furniture('Table', np.array([-8,7.4]),0,0.75,1.1,1.1,0.03,'WOOD')
#    #R12_C1= furniture('Desk', np.array([-6.17,4.92]),90,0.75,1.22,0.8,0.03,'WOOD')
#    #R12_C2= furniture('Desk', np.array([-7.82,6.16]),0,1.22,0.85,0.40,0.03,'WOOD')
#    R12_D= furniture('1.02m-high cupbpard', np.array([-8.73,7.28]),90,1.02,1.2,0.42,0.005,'METAL')
#  
#    #R7_A1= furniture('Groups of 4 desk', np.array([-10.70,12.4]),90,0.75,3.32,2.82,0.03,'WOOD')
#    #R7_A2= furniture('Groups of 4 desk', np.array([-6.92,12.4]),90,0.75,3.32,2.82,0.03,'WOOD')
#    #R7_B = furniture('Block of 6 desks', np.array([-9.49,5.42]),90,0.75,4.82,3.22,0.03,'WOOD')
#    R7_C1= furniture('1.98m-high cupboards', np.array([-13.71,11.69]),90,1.98,1.20,0.42,0.005,'METAL')
#    R7_C2= furniture('1.98m-high cupboards', np.array([-13.71,8.7]),90,1.98,1.20,0.42,0.005,'METAL')
#    R7_D = furniture('Group of 1.02m-high cupboards', np.array([-13.71,4.98]),90,1.02,3.62,0.42,0.005,'METAL')
#      
#    #R10_A = furniture('1.80m-longtable', np.array([-16.7,5]),0,0.75,1.80,0.80,0.03,'WOOD')
#    #R10_B = furniture('1.60m-longtable', np.array([-14.1,5]),90,0.75,1.60,0.80,0.03,'WOOD')
#    R10_C = furniture('0.70m-high cupboards', np.array([-14.1,8.43]),90,0.80,0.90,0.42,0.005,'METAL')
#    R10_D = furniture('1.98m-high cupboards', np.array([-15.52,9.80]),0,1.98,1.20,0.42,0.005,'METAL')
#      
#    #R11_A = furniture('Desk', np.array([-20.1,6.45]),0,0.75,2.10,1.04,0.03,'WOOD')
#    #R11_B= furniture('Table', np.array([-19.3,8.6]),0,0.75,1.1,1.1,0.03,'WOOD')
#    R11_C= furniture('1.98m-high cupbpard', np.array([-17.15,6.6]),90,1.98,1.2,0.42,0.005,'METAL')
#    R11_D1 = furniture('1.02m-high cupbpard', np.array([-17.15,5.4]),90,1.02,1.2,0.42,0.005,'METAL')
#    R11_D2 = furniture('1.02m-high cupbpard', np.array([-20.68,6.15]),90,1.02,1.2,0.42,0.005,'METAL')
#    R11_E1= furniture('0.70m-high cupbpard', np.array([-17.15,4.6]),90,0.70,0.80,0.42,0.005,'METAL')
#    R11_E2= furniture('0.70m-high cupbpard', np.array([-20.68,5.35]),90,0.70,0.80,0.42,0.005,'METAL')
#    R11_F = furniture('1.01m-high cupbpard', np.array([-20.68,4.64]),90,1.01,0.61,0.43,0.01,'METAL')
#    #R8_A = furniture('Block of tabel', np.array([-21.4,12.6]),0,0.75,5.57,2.30,0.03,'WOOD')
#    #R8_B = furniture('Table', np.array([-21.96,14.38]),90,0.75,1.40,0.70,0.03,'WOOD')
#  
#    #R9_A1 = furniture('Block of desks', np.array([-26.75,12.55]),0,0.75,3.2,3.2,0.03,'WOOD')
#    #R9_A2 = furniture('Block of desks', np.array([-25,5.52]),0,0.75,3.2,3.2,0.03,'WOOD')
#    R9_B1 = furniture('1.98m-high cupboards', np.array([-22.74,12.3]),90,1.98,1.20,0.42,0.005,'METAL')
#    R9_B2 = furniture('1.98m-high cupboards', np.array([-21.15,8.81]),90,1.98,1.20,0.42,0.005,'METAL')
#    R9_B3 = furniture('1.98m-high cupboards', np.array([-21.15,4.7]),90,1.98,1.20,0.42,0.005,'METAL')
#    #R9_C = furniture('Table', np.array([-27.71,11.11]),0,0.75,1.8,0.8,0.03,'WOOD')
#    R9_D = furniture('0.70m-high cupboards',np.array([-23,4.6]),0,0.70,0.80,0.43,0.005,'METAL')
#    R9_E1 = furniture('1.98m-high cupboards',np.array([-25.78,4.7]),0,1.33,0.64,0.40,0.005,'METAL')
#    R9_E2 = furniture('1.98m-high cupboards',np.array([-25.78,9.02]),0,1.33,0.64,0.40,0.005,'METAL')
#    R9_F = furniture('1.04m-high cupboards',np.array([-25.6,9.42]),0,1.04,0.42,0.30,0.005,'METAL')
#    R9_G = furniture('Printer',np.array([-23.29,11.7]),0,0.90,0.55,0.40,0.005,'METAL')
#  
#    Fur1 = R6_E
#    Fur2 = R6_F
#    Fur3 = R6_HB
#    Fur4 = R13_A
#    Fur5 = R13_B
#    Fur6 = R13_C
#    Fur7 = R12_A1
#    Fur8 = R12_A2
#    Fur9 = R7_C1
#    Fur10 = R7_C2
#    Fur11 = R10_D
#    Fur12 = R11_C
#    Fur13 = R9_B1
#    Fur14 = R9_B2
#    Fur15 = R9_B3
#    Fur16  = R9_E1
#    Fur17 = R9_E2
#  
#    FL = [Fur1,Fur2,Fur3,Fur4,Fur5,Fur6,Fur7,Fur8,Fur9,Fur10,Fur11,Fur12,Fur13,Fur14,Fur15,Fur16,Fur17]

#    ###
#    ### Dictionary of  the room-point
#    ###
#    # .. todo:: automatiser cela 
#
#    Group = {}
#    Group['R1']=np.array(range(297,333,1))
#    Group['R2']=np.array(range(333,377,1))
#    Group['R6']=np.array(range(255,297,1))
#    Group['R7']=np.array(range(133,232,1))
#    Group['R8']=np.array(range(82,122,1))
#    Group['R9']=np.hstack((np.array(range(1,82,1)),np.array(range(122,133,1))))
#    Group['R10']=np.array([])
#    Group['R11']=np.array([])
#    Group['R12']=np.array([])
#    Group['R13']=np.array([])
#    Group['Couloir right']=np.array(range(232,255,1))
#  
#    S   = Simul.Simul()
#    S.layout('sircut.str','simul9.mat','simul9.slab')
#    S.load
#    Tx,Rx = ptSiradel()
#    #S.load('simul-siradel')
#    S.L.display['Node']=False
#    S.L.display['Scaled']=False
#    S.L.display['Thin']=True
#  
#
#    b    = S.L.boundary()
#    #zone = {}
#    #zone['xmin'] =b[0,0]
#    #zone['xmax'] =b[0,1]
#    #zone['ymin'] =b[1,0]
#    #zone['ymax'] =b[1,1]
#    #zone['Nx']   = 75
#    #zone['Ny']   = 25
#
#    #D=io.loadmat('M1-h1.mat')
#    #D = io.loadmat('newM1h1.mat')
#    D = io.loadmat('M1_essai_new.mat')
#    #
#    # Deal with Matlab file problem
#    #
#    #    Txs  = D['M1h1'][0][0].pTx
#    #    Emax = D['M1h1'][0][0].Emax
#    #    Etot = D['M1h1'][0][0].Etot
#    #    tau_rms =  D['M1h1'][0][0].tau_rms
#    #    tau_m =  D['M1h1'][0][0].tau_moy
#    #    toa_max =  D['M1h1'][0][0].toa_max
#    #    toa_th =  D['M1h1'][0][0].toa_th
#    #    toa_cum =  D['M1h1'][0][0].toa_cum
#    #    toa_new2 = D['M1h1'][0][0].toa_new2
#    #    tau0    =  D['M1h1'][0][0].distance/0.3
#    Txs  = D['M1']['pTx'][0][0]
#    Emax = D['M1']['Emax'][0][0]
#    Etot = D['M1']['Etot'][0][0]
#    Etau0 = D['M1']['Etau0'][0][0]
#    Efirst = D['M1']['Efirst'][0][0]
#    tau_rms =  D['M1']['tau_rms'][0][0]
#    tau_m =  D['M1']['tau_moy'][0][0]
#    toa_max =  D['M1']['toa_max'][0][0]
#    toa_th =  D['M1']['toa_th'][0][0]
#    toa_cum =  D['M1']['toa_cum'][0][0]
#    #toa_new2 = D['M1'][0][0].toa_new2
#    toa_win = D['M1']['toa_win'][0][0]
#    tau0    =  D['M1']['distance'][0][0]/0.3
#    LQI    = D['M1']['LQI1'][0][0]
#    index_conv = D['M1']['index_conv'][0][0]
#    los_cond = D['M1']['los_cond'][0][0][0][0]  
#  
#    u   = Txs[0:2,:].T
#    nd1 = find(LQI[0]<=10)
#    nd2 = find(LQI[1]<=10)
#    nd3 = find(LQI[2]<=10)
#    nd4 = find(LQI[3]<=10)
#
#    nu1 = find(LQI[0]>10)
#    nu2 = find(LQI[1]>10)
#    nu3 = find(LQI[2]>10)
#    nu4 = find(LQI[3]>10)
#  
#    ###
#    ### find the number of the point whose LQI > 10 for all the 4 receptions
#    ###
#    nu12   = np.array([val for val in nu1 if val in nu2])
#    nu34   = np.array([val for val in nu3 if val in nu4])
#    nu1234 = np.array([val for val in nu12 if val in nu34])
#  
#    etoa_max = abs(toa_max - tau0)
#    etoa_cum = abs(toa_cum - tau0)
#    etoa_th  = abs(toa_th - tau0)
#    #etoa_new2 = abs(toa_new2 - tau0)
#    etoa_win = abs(toa_win - tau0)  
#
#    rx1los = los_cond['Rx1'][0][0]['los']-1
#    rx2los = los_cond['Rx2'][0][0]['los']-1
#    rx3los = los_cond['Rx3'][0][0]['los']-1
#    rx4los = los_cond['Rx4'][0][0]['los']-1
#
#    rx1nlos = los_cond['Rx1'][0][0]['nlos']-1
#    rx2nlos = los_cond['Rx2'][0][0]['nlos']-1
#    rx3nlos = los_cond['Rx3'][0][0]['nlos']-1
#    rx4nlos = los_cond['Rx4'][0][0]['nlos']-1
#  
#    rx1nlos2 = los_cond['Rx1'][0][0]['nlos2']-1
#    rx2nlos2 = los_cond['Rx2'][0][0]['nlos2']-1
#    rx3nlos2 = los_cond['Rx3'][0][0]['nlos2']-1
#    rx4nlos2 = los_cond['Rx4'][0][0]['nlos2']-1
#
#    Urx1los   = np.array([val for val in nu1234 if val in rx1los])
#    Urx1nlos  = np.array([val for val in nu1234 if val in rx1nlos])
#    Urx1nlos2 = np.array([val for val in nu1234 if val in rx1nlos2])
#  
#    Urx2los   = np.array([val for val in nu1234 if val in rx2los])
#    Urx2nlos  = np.array([val for val in nu1234 if val in rx2nlos])
#    Urx2nlos2 = np.array([val for val in nu1234 if val in rx2nlos2])
#
#    Urx3los   = np.array([val for val in nu1234 if val in rx3los])
#    Urx3nlos  = np.array([val for val in nu1234 if val in rx3nlos])
#    Urx3nlos2 = np.array([val for val in nu1234 if val in rx3nlos2])
#
#    Urx4los   = np.array([val for val in nu1234 if val in rx4los])
#    Urx4nlos  = np.array([val for val in nu1234 if val in rx4nlos])
#    Urx4nlos2 = np.array([val for val in nu1234 if val in rx4nlos2])
####
#### Modèle MultiWALL Dans Le Batiment Siradel
####
#
#    Emax1 = Emax[0,:]
#    Emax2 = Emax[1,:]
#    Emax3 = Emax[2,:]
#    Emax4 = Emax[3,:]
#
#    Etot1 = Etot[0,:]
#    Etot2 = Etot[1,:]
#    Etot3 = Etot[2,:]
#    Etot4 = Etot[3,:]
#
#    Etau01 = Etau0[0,:]
#    Etau02 = Etau0[1,:]
#    Etau03 = Etau0[2,:]
#    Etau04 = Etau0[3,:]
#
#    Efirst1 = Efirst[0,:]
#    Efirst2 = Efirst[1,:]
#    Efirst3 = Efirst[2,:]
#    Efirst4 = Efirst[3,:]
#
#    Rx1 = Rx[1,0:2]
#    Rx2 = Rx[2,0:2]
#    Rx3 = Rx[3,0:2]
#    Rx4 = Rx[4,0:2]
#
#    D21  = Dgrid_points(u,Rx1)
#    D22  = Dgrid_points(u,Rx2)
#    D23  = Dgrid_points(u,Rx3)
#    D24  = Dgrid_points(u,Rx4)
#
##
## set frequency 
##
#
#    f      = 5.07
#
#    PL1    = OneSlopeMdl(D21,2,f)
#    PL2    = OneSlopeMdl(D22,2,f)
#    PL3    = OneSlopeMdl(D23,2,f)
#    PL4    = OneSlopeMdl(D24,2,f)
#
#    Lw1o,Lw1p  = Loss0_v2(S,u,f,Rx1)
#    Lw2o,Lw2p  = Loss0_v2(S,u,f,Rx2)
#    Lw3o,Lw3p  = Loss0_v2(S,u,f,Rx3)
#    Lw4o,Lw4p  = Loss0_v2(S,u,f,Rx4)
#
#
#
#
####
####  Calculate the diffraction loss due to the obstacle metal
####
#
#    LM = np.array([])
#
#    for i in range(302):
#        for j in range(4):
#            xt = u[i][0]
#            yt = u[i][1]
#            xr = Rx[:,0:2][j+1][0]
#            yr = Rx[:,0:2][j+1][1]
#            for k in range(17):
#                Slist = Interline(xt,yt,xr,yr,FL[k])
#                Lm = Loss_obstacle(Slist,xt,yt,xr,yr,FL[k])
#                LM = hstack((LM,Lm))
#
#    LI = np.array([])
#    LJ = np.array([])
#    NU = np.array([])
#    MM = LM.reshape(302,68)
#    for i in range (302):
#        G = MM[i].reshape(4,17)
#        for j in range(4):
#            Nu = find(G[j]!=0)
#            if len(Nu) >1:
#                print i,j,Nu
#                LI = hstack((LI,i))
#                LJ = hstack((LJ,j))
#                NU = hstack((NU,Nu))
#  
#    LT = np.array([])
#    for i in range(len(LI)):
#      
#      
#        xt = u[LI[i]][0]
#        yt = u[LI[i]][1]
#        xr = Rx[:,0:2][LJ[i]+1][0]
#        yr = Rx[:,0:2][LJ[i]+1][1]
#      
#
#        Lt = Loss_2obstacle(xt,yt,xr,yr,FL[int(NU[(2*i)])],FL[int(NU[(2*i+1)])])  
#        LT = np.hstack((LT,Lt))
#
#    Loss_lien = np.array([])
#    for i in range(302):
#        LSR = MM[i].reshape(4,17)
#        Ls = sum(LSR,axis=1)
#        Loss_lien = hstack((Loss_lien,Ls))
#      
#    Loss_lien = Loss_lien.reshape(302,4)
#  
#    for i in range(56):
#      
#        Loss_lien[LI[i],LJ[i]]= LT[i]
#
#    traj = np.array([1,2,3,5,6,7,8,9,10,11,14,19,20,37,38,43,44,49,50,55,56,61,62,67,68,83,
#        82,81,80,98,99,102,103,107,109,108,106,104,101,100,97,98,80,79,114,115,
#        119,121,123,131,135,140,143,148,147,146,153,154,159,158,157,168,174,178,
#        181,186,188,191,192,195,196,199,200,203,205,206,212,213,214,217,218,221,
#        222,225,226,227,224,223,220,219,216,215,213,212,209,210,228,229,232,234,
#        236,238,240,242,244,243,241,239,237,235,233,231,230,211,208,207,204,202,
#        201,198,197,194,193,190,189,185,182,177,174,169,166,161,151,147,144,139,
#        136,130,124,120,119,116,112,77,76,71,70,65,64,59,58,53,52,47,46,41,40,35,
#        34,33,32,31,30,29,26,25,24,23,22,21,36,20,19,18,17,16])
#
#
#
####
#### Pass loss calculate
####      
#
##        Pt    =  6.368
#    Pt = 0
#
#    Pregle = -1.7227573189526098
#    """
#    Power1o = Pt - PL1 - Lw1o
#    Power2o = Pt - PL2 - Lw2o
#    Power3o = Pt - PL3 - Lw3o
#    Power4o = Pt - PL4 - Lw4o
#    """
#    Power1o = Pt - PL1 - Lw1o - Loss_lien[:,0]
#    Power2o = Pt - PL2 - Lw2o - Loss_lien[:,1]
#    Power3o = Pt - PL3 - Lw3o - Loss_lien[:,2]
#    Power4o = Pt - PL4 - Lw4o - Loss_lien[:,3]
#    """  
#    Power1p = Pt - PL1 - Lw1p
#    Power2p = Pt - PL2 - Lw2p
#    Power3p = Pt - PL3 - Lw3p
#    Power4p = Pt - PL4 - Lw4p
#    """
#    Power1p = Pt - PL1 - Lw1p - Loss_lien[:,0]+Pregle
#    Power2p = Pt - PL2 - Lw2p - Loss_lien[:,1]+Pregle
#    Power3p = Pt - PL3 - Lw3p - Loss_lien[:,2]+Pregle
#    Power4p = Pt - PL4 - Lw4p - Loss_lien[:,3]+Pregle
#  
#    plot(Power1p[traj-1])
#    plot(Power2p[traj-1])
#    plot(Power3p[traj-1])
#    plot(Power4p[traj-1])
#
#    eMW_Ef1 = Power1p - Efirst1
#    eMW_Ef2 = Power2p - Efirst2
#    eMW_Ef3 = Power3p - Efirst3
#    eMW_Ef4 = Power4p - Efirst4
#
#    ### reglementation -1.7228dB
#
#    REG = -1.7228
#
#    OSM1 =  - PL1 + REG
#    OSM2 =  - PL2 + REG
#    OSM3 =  - PL3 + REG
#    OSM4 =  - PL4 + REG
#
#    MWM1 = -PL1-Lw1p + REG
#    MWM2 = -PL2-Lw2p + REG
#    MWM3 = -PL3-Lw3p + REG
#    MWM4 = -PL4-Lw4p + REG
#
#    MWM21 = -PL1-Lw1p-Loss_lien[:,0] + REG
#    MWM22 = -PL2-Lw2p-Loss_lien[:,1] + REG
#    MWM23 = -PL3-Lw3p-Loss_lien[:,2] + REG
#    MWM24 = -PL4-Lw4p-Loss_lien[:,3] + REG
#
#  
#    OSM_moy = mean(hstack((abs((OSM1 - Efirst1)[nu1234]),abs((OSM2 - Efirst2)[nu1234]),abs((OSM3 - Efirst3)[nu1234]),abs((OSM4 - Efirst4)[nu1234]))))
#    OSM_std = std(hstack((abs((OSM1 - Efirst1)[nu1234]),abs((OSM2 - Efirst2)[nu1234]),abs((OSM3 - Efirst3)[nu1234]),abs((OSM4 - Efirst4)[nu1234]))))  
#    MWM_moy = mean(hstack((abs((MWM1 - Efirst1)[nu1234]),abs((MWM2 - Efirst2)[nu1234]),abs((MWM3 - Efirst3)[nu1234]),abs((MWM4 - Efirst4)[nu1234]))))
#    MWM_std = std(hstack((abs((MWM1 - Efirst1)[nu1234]),abs((MWM2 - Efirst2)[nu1234]),abs((MWM3 - Efirst3)[nu1234]),abs((MWM4 - Efirst4)[nu1234]))))  
#    MWM2_moy = mean(hstack((abs((MWM21 - Efirst1)[nu1234]),abs((MWM22 - Efirst2)[nu1234]),abs((MWM23 - Efirst3)[nu1234]),abs((MWM24 - Efirst4)[nu1234]))))
#    MWM2_std = std(hstack((abs((MWM21 - Efirst1)[nu1234]),abs((MWM22 - Efirst2)[nu1234]),abs((MWM23 - Efirst3)[nu1234]),abs((MWM24 - Efirst4)[nu1234]))))  
#    """
#  
#    OSM_moy = mean(hstack(((OSM1 - Efirst1)[nu1234],(OSM2 - Efirst2)[nu1234],(OSM3 - Efirst3)[nu1234],(OSM4 - Efirst4)[nu1234])))
#    :x
#    MWM_moy = mean(hstack(((MWM1 - Efirst1)[nu1234],(MWM2 - Efirst2)[nu1234],(MWM3 - Efirst3)[nu1234],(MWM4 - Efirst4)[nu1234])))
#    MWM_std = std(hstack(((MWM1 - Efirst1)[nu1234],(MWM2 - Efirst2)[nu1234],(MWM3 - Efirst3)[nu1234],(MWM4 - Efirst4)[nu1234])))  
#    MWM2_moy = mean(hstack(((MWM21 - Efirst1)[nu1234],(MWM22 - Efirst2)[nu1234],(MWM23 - Efirst3)[nu1234],(MWM24 - Efirst4)[nu1234])))
#    MWM2_std = std(hstack(((MWM21 - Efirst1)[nu1234],(MWM22 - Efirst2)[nu1234],(MWM23 - Efirst3)[nu1234],(MWM24 - Efirst4)[nu1234])))  
#    """
#    """
##
## export the results of OSM, MWM and MWM2  
##
#
#    F = {}
#    F['OSM']  = np.array([])
#    F['MWM']  = np.array([])
#    F['MWM2'] = np.array([])
#    F['x_Tx'] = np.array([])
#    F['y_Tx'] = np.array([])
#    F['x_Rx'] = np.array([])
#    F['y_Rx'] = np.array([])
#    F['lqi'] = np.array([])  
#    F['id'] = np.array([])  
#    F['Efirst'] = np.array([])
#    F['distance'] = np.array([])  
#    for i in range(302):
#        F['OSM'] = np.hstack((F['OSM'],np.array([OSM1[i],OSM2[i],OSM3[i],OSM4[i]])))
#        F['MWM'] = np.hstack((F['MWM'],np.array([MWM1[i],MWM2[i],MWM3[i],MWM4[i]])))
#        F['MWM2'] = np.hstack((F['MWM2'],np.array([MWM21[i],MWM22[i],MWM23[i],MWM24[i]])))
#        F['x_Tx'] = np.hstack((F['x_Tx'],np.array([u[i][0],u[i][0],u[i][0],u[i][0]])))
#        F['y_Tx'] = np.hstack((F['y_Tx'],np.array([u[i][1],u[i][1],u[i][1],u[i][1]])))
#        F['x_Rx'] = np.hstack((F['x_Rx'],np.array([Rx[1][0],Rx[2][0],Rx[3][0],Rx[4][0]])))
#        F['y_Rx'] = np.hstack((F['y_Rx'],np.array([Rx[1][1],Rx[2][1],Rx[3][1],Rx[4][1]])))
#        F['lqi'] = np.hstack((F['lqi'],np.array([LQI[0][i],LQI[1][i],LQI[2][i],LQI[3][i]])))
#        F['id'] = np.hstack((F['id'],np.array([1,2,3,4])))
#        F['Efirst'] = hstack((F['Efirst'],np.array([Efirst1[i],Efirst2[i],Efirst3[i],Efirst4[i]])))
#        F['distance'] = hstack((F['distance'],np.array([tau0[0][i]*0.3,tau0[1][i]*0.3,tau0[2][i]*0.3,tau0[3][i]*0.3])))
#      
#    io.savemat('M1_h1_Multiwall',F)
#
#"""
#  
#  
