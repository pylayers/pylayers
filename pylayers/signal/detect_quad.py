from scipy.signal import *
from pylayers.signal.DF import DF
import numpy as np
import scipy.fftpack as fft
import matplotlib.pyplot as plt


#def qdetect(w,hz,beta,N):
beta = 1
#"""
#n,x,hmoy
#"""
#
# filtered noise H(z)
#
N  = 10000
w  = np.random.randn(1,N)
wt  = 0.01
fe  = 10
fN  = fe/2.0
ws = np.array([3.168,3.696])/fN
wp = [ws[0]+wt,ws[1]-wt]
lowpass = DF()
gpass = 0.5
gstop = 40
lowpass.ellip_bp(wp,ws,gpass,gstop)
#
# Phi_n(f) = N0/2 | H(f) | **2
#
noise = lowpass.filter(w)
Noise = fft.fft(noise)
PhiN = Noise*np.conj(Noise)/N
varN = sum(PhiN)/N
#
# Rn
#
Rnoise = np.real(fft.ifft(PhiN))
Phi_h2 = 2*Rnoise*Rnoise
#Phi_x  = np.real(ffr(Phi_h2))
Rn0 = max(Rnoise)
#
#  quadratic operator
#
#
x     = (beta*noise)**2
#
# Integrator
#
N = 100
a  = (1./N)*np.ones(N);
integrator = DF(a=a,b=np.array([1]))
y  = integrator.filter(x)

my = np.mean(y)
vy = np.var(y)
