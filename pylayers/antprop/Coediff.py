import numpy as np
import scipy.special as sps
import matplotlib.pyplot as plt
import pdb

def Coediff(k,N,phi0,phi,si,sd,sf,ero,erro,condo,uro,urro,deltaho,
            ern,errn,condn,urn,urrn,deltahn,beta0=np.pi/2):
    """ Luebbers Diffration coefficient


    Parameters
    ----------

    k      -1
    N
    phi0   0
    phi    1
    si
    sd
    sf
    ero
    erro
    sigmao
    uro
    urro
    deltaho
    ern
    errn
    sigman
    urn
    urnn
    deltahn

    """

#--------------------------------------------------
# reflection on faces 'o' and 'n'
#--------------------------------------------------

    if (phi > phi0):
        tho = phi0
        thn = N*np.pi-phi
    else:
        tho = phi
        thn = N*np.pi-phi0

    Rsofto,Rhardo = reflection(tho,k,ero,erro,sigmao,uro,deltaho)
    Rsoftn,Rhardn = reflection(thn,k,ern,errn,sigman,urn,deltahn)


#--------------------------------------------------
# grazing angle Go et Gn
#--------------------------------------------------

    Gsofto,Gsoftn = paramG(N,phi0,Rsofto,Rsoftn)

    Ghardo,Ghardn = paramG(N,phi0,Rhardo,Rhardn)

#--------------------------------------------------
#calcul des 4 termes du coeff diff
#--------------------------------------------------

    sign  =  1.0
    D1     = Dfunc(sign,k,N,phi-phi0,si,sd,beta0)

    sign  =  -1.0
    D2     = Dfunc(sign,k,N,phi-phi0,si,sd,beta0)

    sign  =  +1.0
    D3     = Dfunc(sign,k,N,phi+phi0,si,sd,beta0)

    sign  =  -1.0
    D4     = Dfunc(sign,k,N,phi+phi0,si,sd,beta0)

#--------------------------------------
#n>=1 : exterior wedge
#--------------------------------------
    if (N >= 1.0):
        DTsoft = Gsoftn*(D1+Rsoftn*D3)+Gsofto*(D2+Rsofto*D4)
        DThard = Ghardn*(D1+Rhardn*D3)+Ghardo*(D2+Rhardo*D4)

#--------------------------------------
#traitement des cas ou Go (ou Gn) = -1
#--------------------------------------

        if (abs(Gsoftn+1.0) < 1e-6):
            DTsoft = 0.5*(D1+D3)+Gsofto*(D2+Rsofto*D4)

        if (abs(Gsofto+1.0)<1e-6):
            DTsoft = Gsoftn*(D1+Rsoftn*D3)+0.5*(D2+D4)

        if (abs(Ghardn+1.0) < 1.0e-6):
            DThard = 0.5*(D1+D3)+Ghardo*(D2+Rhardo*D4)

        if (abs(Ghardo+1.0)<1e-6):
            DThard = Ghardn*(D1+Rhardn*D3)+0.5*(D2+D4)

#--------------------------------------
#cas ou n<1 : interior wedge
#--------------------------------------
    else:

        thoz  = N*np.pi-tho
        thnz  = N*np.pi-thn


        [Rsoftnz,Rhardnz] = reflection(thnz,k,ero,erro,condo,uro,deltaho)
        [Rsoftoz,Rhardoz] = reflection(thoz,k,ern,errn,condn,urn,deltahn)

        DTsoft = Rsoftoz*Rsoftnz*D1+Rsoftn*D3+(Rsofto*Rsoftn*D2+Rsofto*D4)

        DThard = Rhardoz*Rhardnz*D1+Rhardn*D3+(Rhardo*Rhardn*D2+Rhardo*D4)

    return DTsoft,DThard


def G(N,phio,Ro,Rn):
    """ grazing angle correction

    Parameters
    ----------

    N : wedge parameter
    phio : incidence angle (rad)
    Ro : reflection coefficient on face o
    Rn : reflection coefficient on face n

    Luebbers 89 "a heuristique UTD slope diffraction coefficient for
                rough lossy wedges"
    """

#----------------------------------------------------
# phio=0
#----------------------------------------------------

    Go = 1.0;
    if (abs(phio) < 1.0e-6 ):
        if (abs(Ro+1.0)>1.0e-6):
            Go = 1.0/(1.0+Ro)
        else:
            Go = -1.0
    elif (abs(phio-N*np.pi)<1.0e-6):
            Go = 0.5

#----------------------------------------------------
# phio=N.PI
#----------------------------------------------------
    Gn = 1.0
    if (abs(phio-N*np.pi) < 1.0e-6 ):
        if (abs(Ro+1.0)>1.0e-6):
            Gn = 1.0/(1.0+Rn)
        else:
            Gn = -1.0
    elif (abs(phio)<1.0e-6):
            Gn = 0.5


def Dfunc(sign,k,N,dphi,si,sd,beta0=np.pi/2):
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

    Reference
    ---------

    [1] KOUYOUMJIAN-PATHAK a uniform geometrical theory of diffraction for an edge
    in a perfectly conducting surface" IEEE AP nov 74 vol 62 N11

    Notes
    -----

            e-jnp.pi/4                1
    Di= ------------------ *    ----------- * F(kla)    ([1] eq 25)
        2n*racine(2*np.pi*k)     np.tan(dphi/n)

    """

    cste = (1.0-1.0*1j)*(1.0/(4.0*N*np.sqrt(k*np.pi)*np.sin(beta0)))
    rnn = (dphi+np.pi*sign)/(2.0*N*np.pi)
    nn  = 0.0

    if (rnn > 0.5):
        nn = 1.0
    if (rnn > 1.5):
        nn = 2.0
    if (rnn < -0.5):
        nn = -1.0
    if (rnn < -1.5 ):
        nn  = -2.0

# KLA

    # [1] eq 27
    L   = (si*sd)/(si+sd)
    AC  = np.cos( (2.0*N*nn*np.pi-dphi) / 2.0 )
    A   = 2*AC*AC
    KLA = k * L * A

    epsi    = AC*2.0
    angle   = (np.pi+sign*dphi)/(2.0*N)
    cot     = np.tan(angle)

    if (abs(cot) > 1e-9 ):
        Fkla = FreF(KLA)
        Di   = -cste*Fkla/cot
    else:
        Di   = 0.5*np.sqrt(L)
    return(Di)

def  FresnelI(x) :
    """ calculates Fresnel integral

    Parameters
    ----------

    x : array
        real argument

    """


    v  = np.zeros(len(x),dtype=complex)
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

    a1=((((((((((
      .23393900e-3*z1 -.12179300e-2)*z1   +.21029670e-2)*z1
      +.2464200e-3)*z1 -.67488730e-2)*z1   +.11948809e-1)*z1
      -.9497136e-2)*z1 +.68989200e-3)*z1   +.57709560e-2)*z1
      +.3936000e-5)*z1 -.24933975e-1)*z1*c1

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
    >>> plt.semilogx(x,np.abs(F))
    >>> plt.grid()

    """

    y     = np.zeros(len(x),dtype=complex)

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

    return y

def FreF2(x):
    """ F function using numpy fresnel function

    Parameters
    ----------
    Not working for large argument

    """
    y     = np.zeros(len(x),dtype=complex)
    u1    = np.where(x>1000)[0]
    u2    = np.where(x<=1000)[0]
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
    return(y,Fc,Fs)


def reflection(th,k,er,err,sigma,ur,urr,deltah):
    """ reflection coeff

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
    >>> th = np.linspace(0,np.pi/2,180)
    >>> fGHz = 0.3
    >>> lamda = 0.3/fGHz
    >>> k = 2*np.pi/2
    >>> Rs,Rh = reflection(th,k,9,0,0.01,1,0,0)

    """

    cel = 2.997925e8
    if not isinstance(k,np.ndarray):
        k = np.array([k])
    if not isinstance(th,np.ndarray):
        th = np.array([th])
    #pdb.set_trace()
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

        yy      = (permi/perme)[:,None]

        st      = np.sin(th)[None,:]
        ct      = np.cos(th)[None,:]

        bb      = np.sqrt(yy-ct**2)

        Rs  = (st - bb) / (st + bb )
        Rh  = (yy*st-bb)/(yy*st+bb)

    else:
        Rs = -ones(len(th))[None,:]
        Rh = ones(len(th))[None,:]

    roughness = 1.0

    Rs  = Rs* roughness
    Rh  = Rh* roughness
    return Rs,Rh
