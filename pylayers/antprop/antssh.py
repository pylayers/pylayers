from pylayers.antprop.antenna import *
from pylayers.antprop.spharm import *
from pylayers.antprop.antenna import *
from pylayers.antprop.antenna import *
import numpy as np
import scipy as sp
import pdb
import sys
import matplotlib.pyplot as plt
import doctest


def SSHFunc(L, theta,phi):
	
	"""
	L : integer, spherical harmonics order
	theta: numpy array(1, nth)
	phi: numpy array(1,nph)
	Compute the spherical harmonic functions for the order L
	return a spherical matrix ((1+L)*(2+L)/2,nth*nph) and the index (l,m) of the shperical harmonics 	
	"""
		
	l = np.arange(0,1+L).reshape(1,(1+L))
	m = np.arange(0,1+L).reshape(((1+L),1))
	# normalize the associated legendre polynoms 
	NRM = np.sqrt((2*l+1)*factorial(l-abs(m))/(4*np.pi*factorial(l+abs(m)))) 
	NRM = NRM.reshape((1+L,1+L,1))
	# compute the associated legendre polynoms part Plm(cos(theta))
	ll = l.reshape(1,(1+L),1)
	mm = m.reshape(((1+L),1,1))
	x  = np.cos(theta).reshape((1,1, len(theta)))
	PLM = sp.special.lpmv(mm,ll,x)
	# Normalize
	NPLM = NRM*PLM
	NPLM = NPLM.reshape((1+L,1+L,len(theta),1))
	# compute the exp(j*m*phi) part
	PHI = phi.reshape((1,len(phi)))
	mm = m.reshape(((1+L),1))
	EXP = np.exp(1j*mm*PHI)
	EXP = EXP.reshape((1+L,1,1,len(phi)))
	# Compute Y : the spherical harmonics matrix and reshaping it
	Yi= NPLM*EXP
	Yi = Yi.reshape(((1+L)**2,len(theta)*len(phi)))
	#~ Y = Yi
	nzero_rows = Yi.any(axis = 1)
	Y = Yi[nzero_rows] # eliminating the non defined functions (Y01) 
	ll = (l*np.ones((1+L,1))).reshape(((1+L)**2))
	mm = (m*np.ones((1,1+L))).reshape(((1+L)**2))
	# spherical harmonics index 
	#~ indx = np.array([ll[nzero_rows],mm[nzero_rows]]).T
	indx = np.array([ll[nzero_rows],mm[nzero_rows]]).T
	Y2 = ((-1)**indx[1+L:,1]).reshape((len(indx[1+L:,1]),1))*np.conj(Y[1+L:,:])
	Y = np.append(Y,Y2, axis  = 0)
	indx2 =  np.array([indx[1+L:,0],(-1)*indx[1+L:,1]]).T
	indx = np.append(indx,indx2, axis = 0)
	return Y,  indx



def SphereToCart (theta, phi, eth, eph, bfreq ):
    """
    Convert from spherical to cartesian coordinates 
    bfreq : boolean parameter to indicate if the conversion is done for all frequencies of only one.
    """

    if bfreq == False:
        PHI = phi.reshape((1,len(phi)))
        THETA = theta.reshape((len(theta),1))
        ec = ndarray(shape = (3, len(theta),len(phi)) , dtype  = complex  )

    else:
        PHI = phi.reshape((1,1,len(phi)))
        THETA = theta.reshape((1,len(theta),1))
        ec = np.ndarray(shape = (3, eth.shape[0], len(theta),len(phi)) , dtype  = complex  )


    ec[0] = np.cos(THETA)*np.cos(PHI)*eth -np.sin(PHI)*eph
    ec[1] = np.cos(THETA)*np.sin(PHI)*eth +np.cos(PHI)*eph
    ec[2] = -np.sin(THETA)*eth

    return ec

def CartToSphere (theta, phi, ex, ey,ez, bfreq = True ):

    """
    Convert from cartesian to spherical  coordinates 
    bfreq : boolean parameter to indicate if the conversion is done for all frequencies of only one.
    """


    if bfreq == False:
        PHI = phi.reshape((1,len(phi)))
        THETA = theta.reshape((len(theta),1))
        es = np.ndarray(shape = (2, len(theta),len(phi)) , dtype  = complex  )

    else:
        PHI = phi.reshape((1,1,len(phi)))
        THETA = theta.reshape((1,len(theta),1))
        es = np.ndarray(shape = (2, ex.shape[0], len(theta),len(phi)) , dtype  = complex  )

    es[0] = np.cos(THETA)*np.cos(PHI)*ex + np.cos(THETA)*np.sin(PHI)*ey -np.sin(THETA)*ez
    es[1] = -np.sin(PHI)*ex +np.cos(PHI)*ey

    return es[0],es[1]

     


def ssh(A, dsf=1, L= 20):
	"""

	Parameters
	----------

	A   :  antenna
	dsf :  int
		down sampling factor  'default 1'

	Summary
	-------

	This function calculates the Scalar Spherical Harmonics coefficients  

		m : phi    longitude
		l : theta  latitude

	Antenna pattern are stored       (f theta phi)
	Coeff are stored with this order (f , l , m )

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


	Etheta  =  A.Ftheta[:,::dsf,::dsf]
	Ephi    =  A.Fphi [:,::dsf,::dsf] 

	# compute the spherical harmonics fucntions at the order L 
	Y,ssh_index  = SSHFunc(L,th,ph)
	# Compute the pseudo inverse of Y
	Ypinv = np.linalg.pinv(Y) 

	# convert the field from spherical to cartesian coordinates system
	Ex, Ey,Ez =SphereToCart (th, ph, Etheta, Ephi, True) 
	#
	Ex = Ex.reshape((nf,nth*nph))
	Ey = Ey.reshape((nf,nth*nph))
	Ez = Ez.reshape((nf,nth*nph))
	#pdb.set_trace()
	cx  =  np.dot(Ex,Ypinv)
	cy  =  np.dot(Ey,Ypinv)
	cz  =  np.dot(Ez,Ypinv)
	lmax = L

	Cx = SCoeff(typ='s2', fmin=A.fa[0], fmax=A.fa[-1],lmax = lmax, data=cx,ind =ssh_index)
	Cy = SCoeff(typ='s2', fmin=A.fa[0], fmax=A.fa[-1],lmax = lmax, data=cy,ind =ssh_index)
	Cz = SCoeff(typ='s2', fmin=A.fa[0], fmax=A.fa[-1],lmax = lmax, data=cz,ind =ssh_index)

	A.S = SSHCoeff(Cx,Cy,Cz)

	return(A)


if (__name__=="__main__"):
	doctest.testmod()
   









