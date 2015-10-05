# coding: utf-8

# This notebook is an implementation from the paper ["Fast UTD Diffraction
# Coefficient Using Only One Suare
# Root"](http://onlinelibrary.wiley.com/doi/10.1002/mop.28298/pdf) from Jean
# Francois Legendre and Thierry Marsault.
# $N$ est compris entre 0 et 1. 0 est la lame de coupeau et 1 est le demi-plan

import numpy as np
import matplotlib.pyplot as plt
import pdb

def diff(fGHz,phi0,phi,si,sd,N,beta):
    """

    axis 0 : fGHz
    axis 1 : phi0 (incident angle on face 0)
    axis 2 : phi  (diffraction from face 0 )
    axis 3 : distance si
    axis 4 : distance sd
    axis 5 : N
    axis 6 : beta

    """

    # MDA reshaping
    # nf,9
    fGHz  = fGHz[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    phi0  = phi0[np.newaxis,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    phi   = phi[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    si    = si[np.newaxis,np.newaxis,np.newaxis,:,np.newaxis,np.newaxis,np.newaxis]
    sd    = sd[np.newaxis,np.newaxis,np.newaxis,np.newaxis,:,np.newaxis,np.newaxis]
    N     = N[np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis,:,np.newaxis]
    beta  = beta[np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis,:]

    L     = si*sd/(si+sd)
    k     = 2*np.pi*fGHz/0.3
    alpha = N*np.sqrt(L*2*k)
    ps4   = np.pi/4
    ps4   = np.pi


    xsi11 = (ps4+phi-phi0)/(2*N)
    xsi12 = (ps4-phi+phi0)/(2*N)
    xsi21 = (ps4-phi-phi0)/(2*N)
    #xsi22 = (phi+phi0-(2*N-1)*ps4)/(2*N)
    xsi22 = (ps4+phi+phi0)/(2*N)

    c11m = np.where(xsi11<0)
    c11p = np.where(xsi11>ps4)
    xsi11[c11m] = xsi11[c11m]*(-1)
    xsi11[c11p] = xsi11[c11p]-ps4

    c12m = np.where(xsi12<0)
    c12p = np.where(xsi12>ps4)
    xsi12[c12m]=xsi12[c12m]*(-1)
    xsi12[c12p]=xsi12[c12p]-ps4


    c21m = np.where(xsi21<0)
    c21p = np.where(xsi21>ps4)
    xsi21[c21m]=xsi21[c21m]*(-1)
    xsi21[c21p]=xsi21[c21p]-ps4


    c22m = np.where(xsi22<0)
    c22p = np.where(xsi22>ps4)
    xsi22[c22m]=xsi22[c22m]*(-1)
    xsi22[c22p]=xsi22[c22p]-ps4


    sg11 = np.ones(np.shape(xsi11)[1:3])[np.newaxis,:,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    sg11[c11m]=-sg11[c11m]
    sg11[c11p]=-sg11[c11p]

    y11 = alpha * np.sin(xsi11)
    cy11m = np.where(y11<0)
    y11[cy11m] = y11[cy11m]*(-1)
    sg11[cy11m] = sg11[cy11m]*(-1)

    sg12 = np.ones(np.shape(xsi12)[1:3])[np.newaxis,:,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    sg12[c12m]=-sg12[c12m]
    sg12[c12p]=-sg12[c12p]

    y12 = alpha * np.sin(xsi12)
    cy12m = np.where(y12<0)
    y12[cy12m] = y11[cy12m]*(-1)
    sg12[cy12m] = sg11[cy12m]*(-1)

    sg21 = np.ones(np.shape(xsi21)[1:3])[np.newaxis,:,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    sg21[c21m]=-sg21[c21m]
    sg21[c21p]=-sg21[c21p]

    y21 = alpha * np.sin(xsi21)
    cy21m = np.where(y21<0)
    y21[cy21m] = y11[cy21m]*(-1)
    sg21[cy21m] = sg11[cy21m]*(-1)

    sg22 = np.ones(np.shape(xsi22)[1:3])[np.newaxis,:,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    sg22[c22m]=-sg22[c22m]
    sg22[c22p]=-sg22[c22p]

    y22 = alpha * np.sin(xsi22)
    cy22m = np.where(y22<0)
    y22[cy22m] = y11[cy22m]*(-1)
    sg22[cy22m] = sg11[cy22m]*(-1)


    y11c = y11*y11
    cinf11 = np.where(y11c<0.5)
    csup11 = np.where(y11c>=0.5)
    yinf11 = y11[cinf11]
    ysup11 = y11[csup11]
    ginf11 = -0.5+0.39384228*yinf11*(1+1j)
    tmp11 = 0.5/y11c[csup11]
    gsup11 = -(0.1994711402/ysup11)*((tmp11+1)+1j*(tmp11-1))
    g11 = np.zeros(np.shape(y11),dtype=complex)
    g11[cinf11]=ginf11
    g11[csup11]=gsup11

    g11 = g11*np.cos(xsi11)*sg11*(1+1j)

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

    g12 = g12*np.cos(xsi12)*sg12*(1+1j)

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

    g21 = g21*np.cos(xsi21)*sg21*(1+1j)

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
    g22 = g22*np.cos(xsi22)*sg22*(1+1j)

    return(g11,g12,g21,g22)


if __name__ == "__main__":

    fGHz = np.linspace(2.4,3.6,10)
    Nphi0 = 100
    Nphi = 550
    Nsi = 150
    Nsd = 60

    fGHz = np.linspace(1,10,100)
    phi0 = np.array([np.pi/4.])
    #phi0 = np.linspace(0,N*np.pi,Nphi0)
    N     = np.array([2])
    phi  = np.linspace(0.01,N[0]*np.pi,Nphi)
    si   = np.array([100])
    sd   = np.array([100])
    #sd   = np.linspace(2,100,100)
    beta = np.array([np.pi/2])
    #si   = np.linspace(0.1,4,Nsi)
    #sd   = np.linspace(0.1,10,Nsd)


    g11,g12,g21,g22 = diff(fGHz,phi0,phi,si,sd,N,beta)
    D = g11+g12+g21+g22
    plt.plot(phi*180/np.pi,abs(D[0,0,:,0,0,0,0]))
    plt.plot(phi*180/np.pi,abs(D[10,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g11[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g12[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g21[0,0,:,0,0,0,0]))
    #plt.plot(phi*180/np.pi,abs(g22[0,0,:,0,0,0,0]))
    #plt.show()
