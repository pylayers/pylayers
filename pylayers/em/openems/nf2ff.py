import h5py
import numpy as np
import matplotlib.pyplot as plt

f = h5py.File('nf2ff.h5','r')
hdf_mesh = f['Mesh']
r     = np.array(hdf_mesh['r'])
theta = np.array(hdf_mesh['theta'])
phi   = np.array(hdf_mesh['phi'])
nf2ff = f['nf2ff']

attrs = nf2ff.attrs.items()
freq = attrs[0][1]
Prad = attrs[1][1]
Dmax = attrs[2][1]

try:
    Eps_r = nf2ff['Eps_r']
except:
    Eps_r = np.ones(freq.shape)
try:
    Mue_r = nf2ff['Mue_r']
except:
    Mue_r = np.ones(freq.shape)

for k,f in enumerate(freq):
    Er_theta = np.array(nf2ff['E_theta']['FD']['f'+str(k)+'_real'])
    Ei_theta = np.array(nf2ff['E_theta']['FD']['f'+str(k)+'_imag'])
    Er_phi = np.array(nf2ff['E_phi']['FD']['f'+str(k)+'_real'])
    Ei_phi = np.array(nf2ff['E_phi']['FD']['f'+str(k)+'_imag'])
    E_theta = Er_theta+1j*Ei_theta
    E_phi   = Er_phi+1j*Ei_phi
    E_norm  = np.sqrt(E_theta*np.conj(E_theta)+E_phi*np.conj(E_phi))
    P_rad   = np.array(nf2ff['P_rad']['FD']['f'+str(k)])

plt.ion()
plt.polar(theta,Er_theta[0,:],'r')
plt.polar(theta,Ei_theta[0,:],'r')
plt.polar(theta,Er_theta[1,:],'b')
plt.polar(theta,Ei_theta[1,:],'b')
plt.figure()
plt.ion()
plt.polar(theta,E_norm[0,:],'b')
plt.polar(theta,E_norm[1,:],'r')

plt.show()

## Calculation of right- and left-handed circular polarization
## adopted from
## 2012, Tim Pegg <teepegg@gmail.com>
##
##% Setup vectors for converting to LHCP and RHCP polarization senses
#for k,f in enumerate(freq):
#    E_cprh[k] = (cos(phi)+1j*sin(phi))*(E_theta[k]+1j*E_phi[k])/np.sqrt(2);
#    E_cplh[k] = (cos(phi)-1j*sin(phi))*(E_theta[k]-1j*E_phi[k])/np.sqrt(2);
