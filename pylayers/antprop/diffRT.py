import numpy as np
import scipy.special as sps
import matplotlib.pyplot as plt
import pdb
def diff(fGHz,phi0,phi,si,sd,N,mat0,matN,beta=np.pi/2,debug=False):
    """ Luebbers Diffration coefficient
    for Ray tracing 



    Parameters
    ----------

    Nf : number of frequencies
    Nr : number of rays

    fGHz : np.array (Nf)
    phi0 : np.array (Nr)
    phi : np.array (Nr)
    si : np.array (Nr)
    sd : np.array (Nr)
    N: np.array (Nb)
    mat0 : Mat
    matN : Mat
    beta : np.array (Nb)
        skew incidence angle (rad)

    Return
    ------

    Ds : numpy array 
        Diffraction soft
    Dh : numpy array 
        Diffraction hard

    Examples
    --------


    >>> import numpy as np
    >>> from pylayers.antprop.slab import *
    >>> Nf=3
    >>> Nr=10
    >>> Nb=5
    >>> fGHz = np.linspace(0,10,Nf)
    >>> N = np.linspace(1,10,Nb)#320/180.
    >>> phi0 = np.linspace(0.01,2*np.pi-0.01,Nr)#40*np.pi/180.
    >>> phi = np.linspace(0.01,2*np.pi-0.01,Nr)
    >>> dm = MatDB()
    >>> mat0 = dm['METAL']
    >>> matN = dm['METAL']
    >>> si = 10000.*np.ones(Nr)
    >>> sd = 1.*np.ones(Nr)
    >>> plt.ion()
    >>> Ds,Dh,D1,D2,D3,D4 = diff(fGHz,phi0,phi,si,sd,N,mat0,matN)

    """

    if not isinstance(fGHz,np.ndarray):
        fGHz = np.array([fGHz])
    if not isinstance(phi0,np.ndarray):
        phi0 = np.array([phi0])
    if not isinstance(phi,np.ndarray):
        phi = np.array([phi])
    if not isinstance(si,np.ndarray):
        si = np.array([si])
    if not isinstance(sd,np.ndarray):
        sd = np.array([sd])
    if not isinstance(N,np.ndarray):
        N = np.array([N])
    if not isinstance(beta,np.ndarray):
        beta = np.array([beta])

    fGHz  = fGHz[:,None]
    phi0  = phi0[None,:]
    phi   = phi[None,:]
    si    = si[None,:]
    sd    = sd[None,:]
    N     = N[None,:]
    beta  = beta[None,:]

    L     = si*sd/(si+sd)
    k     = 2*np.pi*fGHz/0.3

#--------------------------------------------------
# R on faces 'o' and 'n'
#--------------------------------------------------


    tho  = np.empty((fGHz.shape[0],phi.shape[1]))
    thn  = np.empty((fGHz.shape[0],phi.shape[1]))
    # PHI0 = phi0 * np.ones(phi.shape)
    # PHI  = np.ones(phi0.shape)*phi
    # BN   = np.ones(phi0.shape)*N



    c1 = phi>phi0
    c2 = ~c1
    tho[:,c1[0,:]] = phi0[:,c1[0,:]]
    thn[:,c1[0,:]] = N[:,c1[0,:]]*np.pi-phi[:,c1[0,:]]
    tho[:,c2[0,:]] = phi[:,c2[0,:]]
    thn[:,c2[0,:]] = N[:,c2[0,:]]*np.pi-phi0[:,c2[0,:]]



    er0  = np.real(mat0['epr'])
    err0 = np.imag(mat0['epr'])
    ur0  = np.real(mat0['mur'])
    urr0 = np.imag(mat0['mur'])
    sigma0 = mat0['sigma']
    deltah0 = mat0['roughness']

    erN  = np.real(matN['epr'])
    errN = np.imag(matN['epr'])
    urN  = np.real(mat0['mur'])
    urrN = np.imag(mat0['mur'])
    sigmaN = matN['sigma']
    deltahN = matN['roughness']


    Rsofto,Rhardo = R(tho,k,er0,err0,sigma0,ur0,urr0,deltah0)
    Rsoftn,Rhardn = R(thn,k,erN,errN,sigmaN,urN,urrN,deltahN)

#--------------------------------------------------
# grazing angle Go et Gn
#--------------------------------------------------

    Gsofto,Gsoftn = G(N,phi0,Rsofto,Rsoftn)

    Ghardo,Ghardn = G(N,phi0,Rhardo,Rhardn)

#--------------------------------------------------
#calcul des 4 termes du coeff diff
#--------------------------------------------------
    sign  =  1.0
    D1     = Dfunc(sign,k,N,phi-phi0,si,sd,beta)

    sign  =  -1.0
    D2     = Dfunc(sign,k,N,phi-phi0,si,sd,beta)

    sign  =  +1.0
    D3     = Dfunc(sign,k,N,phi+phi0,si,sd,beta)

    sign  =  -1.0
    D4     = Dfunc(sign,k,N,phi+phi0,si,sd,beta)

#--------------------------------------
#n>=1 : exterior wedge
#--------------------------------------
    Dsoft =np.empty(np.shape(D1),dtype=complex)
    Dhard =np.empty(np.shape(D1),dtype=complex)

    #c1 = BN>=1.0

    Dsoft = D1+D2+Rsoftn*D3+Rsofto*D4
    Dhard = D1+D2+Rhardn*D3+Rhardo*D4
#    Dsoft = D2-D4
#    Dhard = D2+D4
    #Dsoft = D1+D2-D3-D4
    #Dhard = D1+D2+D3+D4
#    Dsoft = Gsoftn*(D1+Rsoftn*D3)+Gsofto*(D2+Rsofto*D4)
#    Dhard = Ghardn*(D1+Rhardn*D3)+Ghardo*(D2+Rhardo*D4)
#    c1 = abs(Gsoftn+1.0) < 1e-6
#    c2 = abs(Gsofto+1.0) < 1e-6
#    c3 = abs(Ghardn+1.0) < 1e-6
#    c4 = abs(Ghardo+1.0) < 1e-6
#
#    Dsoft[c1]= 0.5*(D1[c1]+D3[c1])+Gsofto[c1]*(D2[c1]+Rsofto[c1]*D4[c1])
#    Dsoft[c2]= Gsoftn[c2]*(D1[c2]+Rsoftn[c2]*D3[c2])+0.5*(D2[c2]+D4[c2])
#    Dhard[c3]= 0.5*(D1[c3]+D3[c3])+Ghardo[c3]*(D2[c3]+Rhardo[c3]*D4[c3])
#    Dhard[c4]= Ghardn[c4]*(D1[c4]+Rhardn[c4]*D3[c4])+0.5*(D2[c4]+D4[c4])
#--------------------------------------
#traitement des cas ou Go (ou Gn) = -1
#--------------------------------------

#        if (abs(Gsoftn+1.0) < 1e-6):
#            DTsoft = 0.5*(D1+D3)+Gsofto*(D2+Rsofto*D4)
#
#        if (abs(Gsofto+1.0)<1e-6):
#            DTsoft = Gsoftn*(D1+Rsoftn*D3)+0.5*(D2+D4)
#
#        if (abs(Ghardn+1.0) < 1.0e-6):
#            DThard = 0.5*(D1+D3)+Ghardo*(D2+Rhardo*D4)
#
#        if (abs(Ghardo+1.0)<1e-6):
#            DThard = Ghardn*(D1+Rhardn*D3)+0.5*(D2+D4)
#
##--------------------------------------
##cas ou n<1 : interior wedge
##--------------------------------------
#    else:
#
#        thoz  = N*np.pi-tho
#        thnz  = N*np.pi-thn
#
#
#        [Rsoftnz,Rhardnz] = R(thnz,k,ero,erro,condo,uro,deltaho)
#        [Rsoftoz,Rhardoz] = R(thoz,k,ern,errn,condn,urn,deltahn)
#
#        DTsoft = Rsoftoz*Rsoftnz*D1+Rsoftn*D3+(Rsofto*Rsoftn*D2+Rsofto*D4)
#
#        DThard = Rhardoz*Rhardnz*D1+Rhardn*D3+(Rhardo*Rhardn*D2+Rhardo*D4)
    if debug:
        return Dsoft,Dhard,D1,D2,D3,D4
    else :
        return Dsoft,Dhard#,D1,D2,D3,D4


def G(N,phi0,Ro,Rn):
    """ grazing angle correction

    Parameters
    ----------

    N : wedge parameter
    phi0 : incidence angle (rad)
    Ro : R coefficient on face o
    Rn : R coefficient on face n

    Luebbers 89 "a heuristique UTD slope diffraction coefficient for
                rough lossy wedges"
    """


    if not isinstance(phi0,np.ndarray):
        phi0 = np.array([phi0])
    if not isinstance(N,np.ndarray):
        N = np.array([N])

    PHI0 = phi0 * np.ones(Ro.shape)
    BN   = N * np.ones(Ro.shape)

# face o

    Go = np.ones(np.shape(Ro))

    c1 = (abs(PHI0) < 1.0e-6) * (abs(Ro+1.0)>1.0e-6)
    c2 = (abs(PHI0) < 1.0e-6) * (abs(Ro+1.0)<1.0e-6)
    c3 = abs(PHI0-BN*np.pi) < 1.0e-6

    Go[c1] = 1.0/(1.0+Ro[c1])
    Go[c2] = -1.
    Go[c3] = 0.5

# face n
    Gn = np.ones(np.shape(Rn))

    c1 = (abs(PHI0-BN*np.pi) < 1.0e-6) * (abs(Rn+1.0)>1.0e-6)
    c2 = (abs(PHI0-BN*np.pi) < 1.0e-6) * (abs(Rn+1.0)<1.0e-6)
    c3 = abs(PHI0) < 1.0e-6

    Gn[c1] = 1.0/(1.0+Rn[c1])
    Gn[c2] = -1.
    Gn[c3] = 0.5

    return Go,Gn

def Dfunc(sign,k,N,dphi,si,sd,beta=np.pi/2):
    """

    Parameters
    ----------

    sign : int
        +1 | -1
    k : wave number
    N : wedge parameter
    dphi : phi-phi0 or phi+phi0
    si : distance source-D
    sd : distance D-observation
    beta : skew incidence angle

    Reference
    ---------

    [1] KOUYOUMJIAN-PATHAK a uniform geometrical theory of diffraction for an edge
    in a perfectly conducting surface" IEEE AP nov 74 vol 62 N11

    Notes
    -----

            e-jnp.pi/4           1
    Di= ------------------  *  ----------- * F(kla)    ([1] eq 25)
        2n*racine(2*np.pi*k)    np.tan(dphi/n)sin(beta)

    """

    cste = (1.0-1.0*1j)*(1.0/(4.0*N*np.sqrt(k*np.pi)*np.sin(beta)))
    rnn = (dphi+np.pi*sign)/(2.0*N*np.pi)
    nn  =  np.zeros(np.shape(rnn))

    nn[rnn>0.5] = 1
    nn[rnn>1.5] = 2
    nn[rnn<-0.5] = -1
    nn[rnn<-1.5] = -2

# KLA  ref[1] eq 27
    L   = ((si*sd)*np.sin(beta)**2)/(1.*(si+sd))
    AC  = np.cos( (2.0*N*nn*np.pi-dphi) / 2.0 )
    A   = 2*AC**2
    KLA = k * L * A

    epsi  = AC*2.0
    angle = (np.pi+sign*dphi)/(2.0*N)
    tan   = np.tan(angle)

    Di = np.empty(KLA.shape)
    Fkla,ys,yL = FreF(KLA)
    # 4.56 Mac Namara
    Di = -cste*Fkla/tan

    c5 = np.where(np.abs(tan)<1e-9)
    BL = np.ones(Di.shape)*L
    Di[c5] = 0.5*np.sqrt(BL[c5])

    return(Di)

def  FresnelI(x) :
    """ calculates Fresnel integral

    Parameters
    ----------

    x : array
        real argument

    """
    print x.shape

    v  = np.empty(x.shape,dtype=complex)
    y  = np.abs(x)
    z  = .25*y

    u1 = np.where(z>1)
    u2 = np.where(z<=1)

    y1 = y[u1]
    y2 = y[u2]

    d1  = np.cos(y1)
    d2  = np.cos(y2)

    e1  = np.sin(y1)
    e2  = np.sin(y2)

    z1  = z[u1]
    z2  = z[u2]

    c1  = np.sqrt(z1)
    c2  = np.sqrt(z2)

# ----------------------------------------
#  x>4, z>1
# ----------------------------------------

    v1 = 0.5 - 0.5*1j

    c1 = (1.0)/c1
    z1  = c1*c1
    # import ipdb
    # ipdb.set_trace()
    a1=((((((((((
      .23393900e-3*z1 -.12179300e-2)*z1   +.21029670e-2)*z1
      +.2464200e-3)*z1 -.67488730e-2)*z1   +.11948809e-1)*z1
      -.9497136e-2)*z1 +.68989200e-3)*z1   +.57709560e-2)*z1
      +.3936000e-5)*z1 -.24933975e-1)*z1*c1
    # a1=((((((((((.23393900e-3*z1 -.12179300e-2)*z1+.21029670e-2)*z1+.2464200e-3)*z1 -.67488730e-2)*z1   +.11948809e-1)*z1-.9497136e-2)*z1 +.68989200e-3)*z1   +.57709560e-2)*z1+.3936000e-5)*z1 -.24933975e-1)*z1*c1

    b1=(((((((((((
       .838386000e-3*z1  -.55985150e-2)*z1  +.16497308e-1)*z1
      -.27928955e-1)*z1  +.29064067e-1)*z1  -.17122914e-1)*z1
      +.19032180e-2)*z1  +.48514660e-2)*z1  +.23006000e-4)*z1
      -.93513410e-2)*z1  +.23000000e-7)*z1  +.19947114000)*c1

# ----------------------------------------
# x<4, z<1
# ----------------------------------------


    a2=(((((((((((
       0.34404779e-1  *z2 - 0.15023096)*z2 - 0.25639041e-1)*z2
      +0.850663781 )*z2 - 0.75752419e-1 )*z2 - 0.305048566e1)*z2
      -0.16898657e-1 )*z2 + 0.6920691902e1)*z2 - 0.576361e-3 )*z2
      -0.6808568854e1)*z2 - 0.1702e-5)*z2 + 0.159576914e1)*c2

    b2=(((((((((((
       .19547031e-1  *z2 -.216195929e0 )*z2 +.702222016e0)*z2
      -.4033492760e0)*z2 -.1363729124e1)*z2 -.138341947e0)*z2
      +.5075161298e1)*z2 -.952089500e-2)*z2 -.778002040e1)*z2
      -.928100000e-4)*z2 +.4255387524e1)*z2 -.33000000e-7)*c2


    w1    = a1*d1+b1*e1+ 1j*(b1*d1-a1*e1) + v1
    w2    = a2*d2+b2*e2+ 1j*(b2*d2-a2*e2)

    v[u1] = w1
    v[u2] = w2

    y = v*(np.sqrt(np.pi/2.0))

    return y


def FreF(x) :
    """ F function from Pathack

    Parameters
    ----------

    x : array
        real argument

    Examples
    --------

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> x = np.logspace(-4,2,400);
    >>> F = FreF(x)
    >>> plt.semilogx(x,,np.abs(F))
    >>> plt.grid()

    """
    ejp4  = np.exp(1j*np.pi/4)
    emjp4 = np.exp(-1j*np.pi/4)
    y     = np.empty(x.shape,dtype=complex)

    u1    = np.where(x>10)[0]
    u2    = np.where(x<=10)[0]
    xu1   = x[u1]
    xu2   = x[u2]


    x2    = xu1*xu1
    x3    = x2*xu1
    x4    = x3*xu1
    w1    = 1-0.75/x2+4.6875/x4 + 1j*( 0.5/xu1 -1.875/x3)

    cst   = (1.0 - 1j )*0.5*np.sqrt(np.pi/2)
    carx  = abs(xu2)
    racx  = np.sqrt(carx)
    modx  = np.mod(xu2,2*np.pi)
    expjx = np.exp(1j*modx)
    fr    = FresnelI(carx)
    into  = cst - fr
    w2    = 2.0*racx*1j*expjx*into

    y[u1] = w1
    y[u2] = w2

    # [1] eq 30
    ys = (np.sqrt(np.pi*x)-2*x*ejp4-(2/3.)*x**2*emjp4)*np.exp(1j*(np.pi/4+x))
    yl = 1-0.75/(x*x)+4.6875/(x*x*x*x) + 1j*( 0.5/x -1.875/(x*x*x))

    return y,ys,yl

def FreF2(x):
    """ F function using numpy fresnel function

    Parameters
    ----------
    Not working for large argument

    """
    y     = np.empty(x.shape,dtype=complex)
    u1    = np.where(x>5)[0]
    u2    = np.where(x<=5)[0]
    xu1   = x[u1]
    xu2   = x[u2]
    x2    = xu1*xu1
    x3    = x2*xu1
    x4    = x3*xu1
    w1    = 1-0.75/x2+4.6875/x4 + 1j*( 0.5/xu1 -1.875/x3)
    cst   = np.sqrt(np.pi/2.)
    sF,cF = sps.fresnel(np.sqrt(xu2/cst))
    Fc    = (0.5-cF)*cst
    Fs    = (0.5-sF)*cst
    modx  = np.mod(xu2,2*np.pi)
    expjx = np.exp(1j*modx)
    w2    = 2*1j*np.sqrt(xu2)*expjx*(Fc-1j*Fs)
    y[u1] = w1
    y[u2] = w2
    return(y)


def R(th,k,er,err,sigma,ur,urr,deltah):
    """ R coeff

    Parameters
    ----------

    th : np.array
        incidence angle  (axe 0)
    k    : np.array
        wave number      (axe 1)
    er   : real part of permittivity
    err  : imaginary part of permittivity
    sigma : conductivity
    ur   : real part of permeability
    urr  : imaginary part of permeability
    deltah : height standard deviation

    Examples
    --------

    >>> import numpy as np
    >>> th = np.linspace(0,np.pi/2,180)[None,:]
    >>> fGHz = 0.3
    >>> lamda = 0.3/fGHz
    >>> k = np.array([2*np.pi/2])[:,None]
    >>> Rs,Rh = R(th,k,9,0,0.01,1,0,0)

    """

    cel = 299792458
    #--------------------------------------------
    #cas des surfaces dielectriques (sinon er=-1)
    #--------------------------------------------

    if (er >= 0.0 ):
        if ( (( ur-1.0)<1e-16) & ((er-1.0)<1e-16) ):
            Rs = np.zeros(len(th))
            Rh = np.zeros(len(th))

        u1 = np.where(th >= 1.5*np.pi)
        u2 = np.where(th >= np.pi )
        u3 = np.where(th >= 0.5*np.pi)

        th[u1] = 2.0*np.pi - th[u1]
        th[u2] = th[u2] - np.pi
        th[u3] = np.pi - th[u3]

        #if (th >= 1.5*np.pi ):
        #    th = 2.0*np.pi - th
        #elif (th >= np.pi ):
        #    th = th - np.pi
        #elif (th >= 0.5*np.pi):
        #    th = np.pi - th

        uo   = 4.0*np.pi*1e-7
        eo   = 1.0/(uo*cel*cel)

        pulse   = k*cel
        permi   = (er-1j*err)-(1j*sigma)/(pulse*eo)

        perme   = ur - 1j*urr

        yy      = (permi/perme)

        st      = np.sin(th)
        ct      = np.cos(th)

        bb      = np.sqrt(yy-ct**2)

        Rs  = (st - bb) / (st + bb )
        Rh  = (yy*st-bb)/(yy*st+bb)

    else: # metalic case
        Rs = -np.ones(th.shape)
        Rh =  np.ones(th.shape)

    roughness = 1.0

    Rs  = Rs* roughness
    Rh  = Rh* roughness
    return Rs,Rh
