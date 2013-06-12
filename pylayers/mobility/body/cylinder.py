import numpy as np 
import scipy as sp 
import matplotlib.pylab as plt 

# Cylinder radius
R = 2.8
# Frequency 
fMHz = 1000.
# wavelength 
lamda = 300/fMHz
#  Impedance
Z0 = 50 
# radial distance 
dst = 0 

#
# Transmitter characteristics
#
Feeder_Loss_dB_Tx = 0 
#
# length 
# 
l_em = 0.063
r_em = 0.004
g_em_dipole_dBi = 1.852
g_em_monopole_dBi = 4.857 
Gmax_dBi = 3.9
#
# 
#
freq_accord = 300/(4*l_em)
kl_em = 2*np.pi*l_em/lamda
#
#
Ze_r = 20 * (kl_em*kl_em)
Ze_i = 120/np.tan(kl_em)*(np.log(l_em/r_em)-1)
#
#
S11e = np.sqrt(((Ze_r-Z0)**2 + Ze_i**2)/((Ze_r+Z0)**2 + Ze_i**2))
S21e = np.sqrt(1-S11e**2)
#
#
Eff_em = 20*np.log10(S21e)- Feeder_Loss_dB_Tr

# Receiver Characteristics
#
Feeder_Loss_dB_Rx = 0 
#
# Length
#
l_rc = 0.063
r_rc = 0.004
#
#
#
g_rc_dip  = 1.854
g_rc_mono = 4.854

g_rc_max  = 3.9
Freq_accord = 300/(4*l_rc)

kl_rc = (2.*pi/lamda)*l_rc

Zr_r = 20 * (kl_rc*kl_rc)
Zr_i = 120/np.tan(kl_rc)*(log(l_rc/r_rc)-1)

S11r = np.sqrt(((Ze_r-Z0)**2 + Zr_i**2)/((Zr_r+Z0)**2 + Zr_i**2))
S21r = np.sqrt(1-S11**2)


#
# Boithias formule 7.17 page 166 
#
delta2 = delta*conj(delta)
num = delta2  + 1.6*np.sqrt(delta2)+0.75
den = delta2  + 4.5*np.sqrt(delta2)+1.35
belta = num/den
