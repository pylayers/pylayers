import scipy.signal as si
import numpy as np
import scipy as sp


def weight(N,**kwags):
    default = {'sll':20,
               'typ','uniform'
              }

    for k in defaults:
        fi k not in kwargs:
            kwargs[k]=defaluts[k]

    if typ=='uniform':
        pass

def rejection(w,theta,fGHz=10,d=[],thresh=0.1):
    """
    Calculate rejection for a given arbitrary array

    Parameters
    ----------

    w     : weighting coefficient (can be complex)
    theta : theta interval
    fGHz  : frequency in GHz
    d     : interelement distance (default lambda/2)
    thresh : threshold for null first derivative evaluation

    Notes
    -----


    In order to expressed criteria on rejection it is important to give a precise mathematical definition
    of what exactly rejection is.
    Let :math:`\mathbf{w}^T` be the 1xN vector of antenna array weights.
    The complex array factor is given by :

. math::

\mathbf{F}(\theta)= \mathbf{w}^{\dagger} . \mathbf{S}(\theta)

Let defines the  :math:`N\times N` matrix

.. math::

\mathbf{W}=\mathbf{w}^{\dagger}\mathbf{w}

..math::

|F(\theta)|^2 = \textrm{diag}(\mathbf{F}^{\dagger}(\theta)\mathbf{F})
|F(\theta)|^2 = \textrm{diag}( \mathbf{S}^{\dagger} \mathbf{W} \mathbf{S} )
\frac{d}{d\theta}|F(\theta)|^2 = \textrm{diag}( \mathbf{U}^{\dagger} \mathbf{W}  \mathbf{S} )+\textrm{diag}( \mathbf{S}^{\dagger} \mathbf{W}  \mathbf{U} )$$
\frac{d^2}{d\theta^2}|F(\theta)|^2 =\textrm{diag}( \mathbf{R}^{\dagger} \mathbf{W}  \mathbf{S} )+ 
                                       \textrm{diag}( \mathbf{U}^{\dagger} \mathbf{W}  \mathbf{U} )+
                                       \textrm{diag}( \mathbf{U}^{\dagger} \mathbf{W}  \mathbf{U})+
                                       \textrm{diag}( \mathbf{S}^{\dagger} \mathbf{W}  \mathbf{R} )

 \frac{d^2}{d\theta^2}|F(\theta)|^2 =\textrm{diag}( \mathbf{R}^{\dagger} \mathbf{W}  \mathbf{S} )+ 
                                       2\textrm{diag}( \mathbf{U}^{\dagger} \mathbf{W}  \mathbf{U} )+
                                       \textrm{diag}( \mathbf{S}^{\dagger} \mathbf{W}  \mathbf{R} )
    """
    N  = len(w)
    n  = arange(N)
    lam   = 0.3/fGHz
    k     = 2*pi/lam
    if d ==[]:
        d     = lam/2

    W  = outer(conj(w.T),w)
    u  = 1j*k*d*outer(n,sin(theta))
    v  = 1j*k*d*outer(n,cos(theta))
    S  = exp(u)
    U  = v*exp(u)
    R  = -u*exp(v)+v*v*exp(u)
    T  = S+U
    F  = real(dot(conj(S.T),dot(W,S)))
    G  = real(dot(conj(U.T),dot(W,S))+dot(conj(S.T),dot(W,U)))
    H  = real(dot(conj(R.T),dot(W,S))+2*dot(conj(U.T),dot(W,U))+dot(conj(S.T),dot(W,R)))
    f  = diag(F)/max(diag(F))
    g  = diag(G)/max(diag(G))
    h  = diag(H)/max(diag(H))
    # max condition (first derivative absolute value below threshold and second derivative <0)
    z1     = nonzero((abs(g)<thresh) & (h < 0))[0]
    # find mainlobe
    ml = nonzero((abs(f)>=0.5) & (h <0))[0]
    # exclude mainlobe from maxima
    z    = setdiff1d(z1,ml)
    rejdB = log10(max(abs(f[z]))/max(abs(f)))
    bw = theta[ml[-1]]-theta[ml[0]]
    bwdeg = bw*180/pi
    #print max(abs(f[z])),max(abs(f)),rejdB
    return(f,g,h,z,ml,rejdB,bwdeg)


def visurej(f,g,h,z,main,rejdB,bwdeg,titre1):
    plot(thetadeg,log10(abs(f)),'k')
    plot(thetadeg,g,'r')
    plot(thetadeg,sign(h),'b')
    title('Attempt to identify the sidelobes maxima : imprecise approach')
    z     = nonzero((abs(g)<0.01) & (h < 0))[0]
    plot(thetadeg[z],log10(abs(f[z])),'ro')
    plot(thetadeg[main],log10(abs(f[main])),'go')
    #rej=log10(max(abs(f[z]))/max(abs(f)))
    axis((-90,90,-5,2))
    legend(('F(dB)/10 ','normalized derivative (linear)','sign of 2nd der'))
    titre = titre1+'beamwidth (deg) : '+ str(round(bwdeg*100)/100)+' achieved rejection (dB) : ' + str(round(rejdB*100)/10)
    title(titre)
