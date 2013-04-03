import pdb
import doctest
import pylayers.antprop.antenna as ant
from pylayers.antprop.spharm import *
from sphere import spherepack, Wrapec, mathtogeo
import numpy as np

def vsh(A, dsf=1):
    """

    Parameters
    ----------
    A   :  antenna
    dsf :  int
        down sampling factor  'default 1'

    Summary
    -------

    This function calculates the Vector Spherical Harmonics coefficients
    It makes use of the spherepack function vha

        m : phi    longitude
        n : theta  latitude

    Antenna pattern are stored       (f theta phi)
    Coeff are stored with this order (f , n , m )

    The vsh coefficient are organized differently
    should be better for compression along frequency axis


    """

    th = A.theta[::dsf]
    ph = A.phi[::dsf]

    nth = len(th)
    nph = len(ph)
    nf = A.Nf

    if (nph % 2) == 1:
        mdab = min(nth, (nph + 1) / 2)
    else:
        mdab = min(nth, nph / 2)

    ndab = nth

    Br = 1j * np.zeros((nf, ndab, mdab))
    Bi = 1j * np.zeros((nf, ndab, mdab))
    Cr = 1j * np.zeros((nf, ndab, mdab))
    Ci = 1j * np.zeros((nf, ndab, mdab))

    gridComp = Wrapec()
    wvha, lvha = gridComp.vhai(nth, nph)

    for k in range(nf):
        #
        # Real part
        #
        Fpr = A.Fphi[k][::dsf, ::dsf].real
        Ftr = A.Ftheta[k][::dsf, ::dsf].real
        #
        # Fpr     Ntheta,Nphi
        #
        brr, bir, crr, cir = gridComp.vha(nth, nph, 1,
                                          lvha, wvha,
                                          np.transpose(Fpr),
                                          np.transpose(Ftr))
        #
        # Imaginary part
        #
        Fpi = A.Fphi[k][::dsf, ::dsf].imag
        Fti = A.Ftheta[k][::dsf, ::dsf].imag
        bri, bii, cri, cii = gridComp.vha(nth, nph, 1,
                                          lvha, wvha,
                                          np.transpose(Fpi),
                                          np.transpose(Fti))

        Br[k, :, :] = brr + 1j * bri
        Bi[k, :, :] = bir + 1j * bii
        Cr[k, :, :] = crr + 1j * cri
        Ci[k, :, :] = cir + 1j * cii

    #
    # m=0 row is multiplied by 0.5
    #

    Br[:, :, 0] = 0.5 * Br[:, :, 0]
    Bi[:, :, 0] = 0.5 * Bi[:, :, 0]
    Cr[:, :, 0] = 0.5 * Cr[:, :, 0]
    Ci[:, :, 0] = 0.5 * Ci[:, :, 0]

    #print "A.fa[0] = ",A.fa[0]
    #print "A.fa[-1] = ",A.fa[-1]

    Br = ant.SHCoeff(typ='s1', fmin=A.fa[0], fmax=A.fa[-1], data=Br)
    Bi = ant.SHCoeff(typ='s1', fmin=A.fa[0], fmax=A.fa[-1], data=Bi)
    Cr = ant.SHCoeff(typ='s1', fmin=A.fa[0], fmax=A.fa[-1], data=Cr)
    Ci = ant.SHCoeff(typ='s1', fmin=A.fa[0], fmax=A.fa[-1], data=Ci)

    A.C = ant.VSHCoeff(Br, Bi, Cr, Ci)
    return(A)

if (__name__ == "__main__"):
    doctest.testmod()

