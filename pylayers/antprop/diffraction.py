# coding: utf-8

# This notebook is an implementation from the paper ["Fast UTD Diffraction
# Coefficient Using Only One Suare
# Root"](http://onlinelibrary.wiley.com/doi/10.1002/mop.28298/pdf) from Jean
# Francois Legendre and Thierry Marsault.
# $N$ est compris entre 0 et 1. 0 est la lame de coupeau et 1 est le demi-plan

import numpy as np
import matplotlib.pyplot as plt
import pdb

def diff(**kwargs):
    """this methode handles coefficient diffraction using one square root
    according to work of J.F.Legendre and T.Marsault

    Parameters
    ----------

    beta : exterior angle of wedge (axis 0)
    phi0 : incident angle on face 0 (axis 1)
    phi  : angle of diffraction from face 0 (axis 2)
    si   : distance from the source (axis 3)
    sd   : distance from the diffraction (axis 4)
    N    : aperture of wedge (axis 5)
    fGHz : frequency (axis 6)

    """

    # MDA reshaping
    # nf,9
    defaults = { 'beta' : np.array([np.pi/2.]),
                 'phi0' : np.array([np.pi/2.]),
                 'phi'  : np.linspace(0.01,np.array([150])[0]*np.pi,150),
                 'si'   : np.array([100]),
                 'sd'   : np.array([100]),
                 'N'    : np.array([1.5]),
                 'fGHz' :np.linspace(2.4,3.6,10)
                   }

    for k in defaults:
        if k not in kwargs:
            kwargs[k]=defaults[k]

    beta = kwargs.pop('beta')
    phi0 = kwargs.pop('phi0')
    phi  = kwargs.pop('phi')
    si   = kwargs.pop('si')
    sd   = kwargs.pop('sd')
    N    = kwargs.pop('N')
    fGHz = kwargs.pop('fGHz')


    beta  = beta[:,None,None,None,None,None,None]
    phi0  = phi0[None,:,None,None,None,None,None]
    phi   = phi[None,None,:,None,None,None,None]
    si    = si[None,None,None,:,None,None,None]
    sd    = sd[None,None,None,None,:,None,None]
    N     = N[None,None,None,None,None,:,None]
    fGHz  = fGHz[None,None,None,None,None,None,:]

    L     = si*sd/(si+sd)
    k     = 2*np.pi*fGHz/0.3
    alpha = N*np.sqrt(L*2*k)
    ps4   = np.pi/4
    #ps4   = np.pi


    zeta11 = (np.pi+phi-phi0)/(2*N)
    zeta12 = (np.pi-phi+phi0)/(2*N)
    zeta21 = (np.pi-phi-phi0)/(2*N)
    zeta22 = (phi+phi0-(2*N-1)*np.pi)/(2*N)

    #zeta11 = (ps4+phi-phi0)/(2*N)
    #zeta12 = (ps4-phi+phi0)/(2*N)
    #zeta21 = (ps4-phi-phi0)/(2*N)
    #zeta22 = (phi+phi0-(2*N-1)*ps4)/(2*N)

    c11m = np.where(zeta11<0)
    c11p = np.where(zeta11>(ps4))
    zeta11[c11m] = zeta11[c11m]*(-1)
    zeta11[c11p] = zeta11[c11p]-ps4

    c12m = np.where(zeta12<0)
    c12p = np.where(zeta12>(ps4))
    zeta12[c12m]=zeta12[c12m]*(-1)
    zeta12[c12p]=zeta12[c12p]-ps4


    c21m = np.where(zeta21<0)
    c21p = np.where(zeta21>(ps4))
    zeta21[c21m]=zeta21[c21m]*(-1)
    zeta21[c21p]=zeta21[c21p]-ps4


    c22m = np.where(zeta22<0)
    c22p = np.where(zeta22>(ps4))
    zeta22[c22m]=zeta22[c22m]*(-1)
    zeta22[c22p]=zeta22[c22p]-ps4


    #sg11 = np.ones(np.shape(zeta11)[1:3])[None,:,:,None,None,None,None]
    sg11 = np.ones(np.shape(zeta11)[1:4])[None,:,:,None,None,:,None]
    sg11[c11m]=-sg11[c11m]
    sg11[c11p]=-sg11[c11p]


    #y11 = alpha * np.sin(zeta11)
    sinq11 = zeta11*(1-(1/6.*zeta11*zeta11)*(1-(1/20.*zeta11*zeta11)*(1-(1/7.*zeta11*zeta11))))
    #sinq11 = zeta11*((1-1/6)*zeta11*zeta11*(1-0.05*zeta11*zeta11*(1-0.1428571429*zeta11*zeta11)))
    y11 = alpha * sinq11

    cy11m = np.where(y11<0)
    y11[cy11m] = y11[cy11m]*(-1)
    sg11[cy11m] = sg11[cy11m]*(-1)

    #sg12 = np.ones(np.shape(zeta12)[1:3])[None,:,:,None,None,None,None]
    sg12 = np.ones(np.shape(zeta12)[1:4])[None,:,:,None,None,:,None]
    sg12[c12m]=-sg12[c12m]
    sg12[c12p]=-sg12[c12p]


    #sinq12 = zeta12*((1-1/6)*zeta12*zeta12*(1-0.05*zeta12*zeta12*(1-0.1428571429*zeta12*zeta12)))
    sinq12 = zeta12*(1-(1/6.*zeta12*zeta12)*(1-(1/20.*zeta12*zeta12)*(1-(1/7.*zeta12*zeta12))))
    #y12 = alpha * np.sin(zeta12)
    y12 = alpha * sinq12

    cy12m = np.where(y12<0)
    y12[cy12m] = y11[cy12m]*(-1)
    sg12[cy12m] = sg11[cy12m]*(-1)

    #sg21 = np.ones(np.shape(zeta21)[1:3])[None,:,:,None,None,None,None]
    sg21 = np.ones(np.shape(zeta21)[1:4])[None,:,:,None,None,:,None]
    sg21[c21m]=-sg21[c21m]
    sg21[c21p]=-sg21[c21p]


    #y12 = alpha * np.sin(zeta12)
    sinq21 = zeta21*(1-(1/6.*zeta21*zeta21)*(1-(1/20.*zeta21*zeta21)*(1-(1/7.*zeta21*zeta21))))
    y21 = alpha * sinq21

    cy21m = np.where(y21<0)
    y21[cy21m] = y11[cy21m]*(-1)
    sg21[cy21m] = sg11[cy21m]*(-1)

    #sg22 = np.ones(np.shape(zeta22)[1:3])[None,:,:,None,None,None,None]
    sg22 = np.ones(np.shape(zeta22)[1:4])[None,:,:,None,None,:,None]
    sg22[c22m]=-sg22[c22m]
    sg22[c22p]=-sg22[c22p]


    sinq22 = zeta22*(1-(1/6.*zeta22*zeta22)*(1-(1/20.*zeta22*zeta22)*(1-(1/7.*zeta22*zeta22))))
    #y22 = alpha * np.sin(zeta22)
    y22 = alpha * sinq22

    cy22m = np.where(y22<0)
    y22[cy22m] = y11[cy22m]*(-1)
    sg22[cy22m] = sg11[cy22m]*(-1)


    y11c = y11*y11
    cinf11 = np.where(y11c<0.5)
    csup11 = np.where(y11c>=0.5)
    #cinf11 = np.where(y11c<0.005)
    #csup11 = np.where(y11c>=0.005)
    yinf11 = y11[cinf11]
    ysup11 = y11[csup11]
    ginf11 = -0.5+0.39384228*yinf11*(1+1j)
    tmp11 = 0.5/y11c[csup11]
    gsup11 = -(0.1994711402/ysup11)*((tmp11+1)+1j*(tmp11-1))
    g11 = np.zeros(np.shape(y11),dtype=complex)
    g11[cinf11]=ginf11
    g11[csup11]=gsup11


    cosq11 = (1-(0.5*zeta11*zeta11)*(1-(1/12.*zeta11*zeta11)*(1-(1/30.*zeta11*zeta11))*(1-(1/56.*zeta11*zeta11))))
    #g11 = g11*np.cos(zeta11)*sg11*(1+1j)
    g11 = g11*cosq11*sg11*(1+1j)

    y12c = y12*y12
    cinf12 = np.where(y12c<0.5)
    csup12 = np.where(y12c>=0.5)
    yinf12 = y12[cinf12]
    ysup12 = y12[csup12]
    ginf12 = -0.5+0.39384228*yinf12*(1+1j)
    tmp12 = 0.5/y12c[csup12]
    gsup12 = -(0.1994711402/ysup12)*((tmp12+1)+1j*(tmp12-1))
    g12 = np.zeros(np.shape(y12),dtype=complex)
    g12[cinf12]=ginf12
    g12[csup12]=gsup12

    cosq12 = (1-(0.5*zeta12*zeta12)*(1-(1/12.*zeta12*zeta12)*(1-(1/30.*zeta12*zeta12))*(1-(1/56.*zeta12*zeta12))))
    g12 = g12*cosq12*sg12*(1+1j)
    #g12 = g12*np.cos(zeta12)*sg12*(1+1j)

    y21c = y21*y21
    cinf21 = np.where(y21c<0.5)
    csup21 = np.where(y21c>=0.5)
    yinf21 = y21[cinf21]
    ysup21 = y21[csup21]
    ginf21 = -0.5+0.39384228*yinf21*(1+1j)
    tmp21 = 0.5/y21c[csup21]
    gsup21 = -(0.1994711402/ysup21)*((tmp21+1)+1j*(tmp21-1))
    g21 = np.zeros(np.shape(y21),dtype=complex)
    g21[cinf21]=ginf21
    g21[csup21]=gsup21

    cosq21 = (1-(0.5*zeta21*zeta21)*(1-(1/12.*zeta21*zeta21)*(1-(1/30.*zeta21*zeta21))*(1-(1/56.*zeta21*zeta21))))
    g21 = g21*cosq21*sg21*(1+1j)
    #g21 = g21*np.cos(zeta21)*sg21*(1+1j)

    y22c = y22*y22
    cinf22 = np.where(y22c<0.5)
    csup22 = np.where(y22c>=0.5)
    yinf22 = y22[cinf22]
    ysup22 = y22[csup22]
    ginf22 = -0.5+0.39384228*yinf22*(1+1j)
    tmp22 = 0.5/y22c[csup22]
    gsup22 = -(0.1994711402/ysup22)*((tmp22+1)+1j*(tmp22-1))
    g22 = np.zeros(np.shape(y22),dtype=complex)
    g22[cinf22]=ginf22
    g22[csup22]=gsup22

    cosq22 = (1-(0.5*zeta22*zeta22)*(1-(1/12.*zeta22*zeta22)*(1-(1/30.*zeta22*zeta22))*(1-(1/56.*zeta22*zeta22))))
    g22 = g22*cosq22*sg22*(1+1j)
    #g22 = g22*np.cos(zeta22)*sg22*(1+1j)

    #return(g11,g12,g21,g22)

    #computation of field
    #pdb.set_trace()
    s0 = np.array([100])
    E = np.exp(-1j*k*(sd+si-s0))
    cste = s0/sd+si

    #Amplitude field

    A11 = cste*g11*E
    A12 = cste*g12*E
    A21 = cste*g21*E
    A22 = cste*g22*E

    return(A11,A12,A21,A22,g11,g12,g21,g22)


if __name__ == "__main__":

    fGHz  = np.linspace(2.4,3.6,10)
    Nphi0 = 100
    Nphi  = 150
    #Nsi   = 150
    #Nsd   = 60

    #fGHz = np.linspace(1,10,100)
    #phi0 = np.array([np.pi/2.])
    phi0 = np.array([np.pi/2.])
    N     = np.array([1.5])
    phi  = np.linspace(0.01,N[0]*np.pi,Nphi)
    si   = np.array([100])
    sd   = np.array([100])
    #sd   = np.linspace(2,100,100)
    beta = np.array([np.pi/2.])
    #si   = np.linspace(0.1,4,Nsi)
    #sd   = np.linspace(0.1,10,Nsd)


    #g11,g12,g21,g22 = diff(fGHz,phi0,phi,si,sd,N,beta)
    #g11,g12,g21,g22 = diff(beta=beta,phi0=phi0,phi=phi,si=si,sd=sd,N=N,fGHz=fGHz)
    A11,A12,A21,A22,g11,g12,g21,g22 = diff(beta=beta,phi0=phi0,phi=phi,si=si,sd=sd,N=N,fGHz=fGHz)


    D = g11+g12+g21+g22
    A = A11+A12+A21+A22

    a = 20*np.log10(np.abs(A[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,a)
    #plt.xlabel(r'$\phi(degrees)')
    #plt.ylabel('Amplitude')
    #si    = si[None,None,None,:,None,None,None]
    #sd    = sd[None,None,None,None,:,None,None]
    #k     = 2*np.pi*fGHz/0.3
    #s0 = np.array([100])
    #E = np.exp(-1j*k*(sd+si-s0))
    #cste = s0/sd+si

    #plt.plot(phi*180/np.pi,np.abs(A[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,np.abs(D[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g11[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g12[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g21[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g22[0,0,:,0,0,0,0]))
    #plt.show()
