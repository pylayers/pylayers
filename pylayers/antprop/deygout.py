import numpy as np
import time
import pdb
from numba import jit
@jit
def cover(X,Y,Z,Ha,Hb,fGHz):
    Nphi,Nr = Z.shape
    #pdb.set_trace()
    #print Nphi,Nl
    L = np.zeros((Nphi,Nr-2))
    for ip in xrange(Nphi):
        for il in xrange(Nr-1):
            uk = np.arange(0,il+1)
            z = np.empty(len(uk))
            x = X[ip,uk]
            y = Y[ip,uk]
            z[uk] = Z[ip,uk]
            d = np.sqrt((x-x[0])**2+(y-y[0])**2)
            z[0]  = z[0] + Ha
            z[-1] = z[-1] + Hb
            #plt.plot(d,z)
            v = Deygout(z,d,fGHz,0,0)
            #print v
            L[ip,il-1]=v
    return(L)

@jit
def Deygout(z,d,fGHz,L,depth):
    lmbda = 0.3/fGHz
    depth = depth+1
    if depth <3:
        if len(z)>3:
            u = np.arange(len(z))/(len(z)-1.0)
            l = (z[0])*(1-u)+(z[-1])*u
            h = z[1:-1]-l[1:-1]
            nu = h*np.sqrt((2/lmbda)*(1/d[1:-1]+1/(d[-1]-d[1:-1])))
            imax = np.nanargmax(nu)
            numax = nu[imax]
        else:
            numax = -10
        if numax>-0.78:
            w  = numax -0.1
            L  = L + 6.9 + 20*np.log10(np.sqrt(w**2+1)+w)
            z1 = z[0:imax]
            d1 = d[0:imax]
            Ll = Deygout(z1,d1,fGHz,0,depth)
            z2 = z[imax:]
            d2 = d[imax:]
            Lr = Deygout(z2,d2,fGHz,0,depth)
            L  = L+Lr+Ll
    return(L)
