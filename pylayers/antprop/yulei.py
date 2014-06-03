# Deprecated
def Dis(x1,y1,x2,y2,x3,y3):
    """ Distance between a point and a line

    Parameters
    ----------

    point (x1,y1)
    line  {(x2,y2),(x3,y3)}

    Returns
    -------

    h


    See Also
    --------

    Loss_Obstacle

    OBSOLETE

    """
    A = (y3-y2)/(x3-x2)
    B = (x3*y2-x2*y3)/(x3-x2)

    # h: distance
    # ax+by+c=0   h = |a*xo+b*yo+c|/np.sqrt(a*a+b*b)
    #

    # ????
    h = abs(A*x1-y1+B)/np.sqrt(A*A+1)

    return h
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

    OBSOLETE

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

def Interline(x1,y1,x2,y2,Obstacle):
    """ Intersection furniture with line
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
def Loss_obstacle(SS,x1,y1,x2,y2,Obstacle):
    """ calculate loss due to furniture

    Parameters
    ----------

    SS
    x1
    y1
    x2
    y2
    Obstacle


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
