from pylayers.signal.invfreqs import *
from scipy.signal import *
import matplotlib.pyplot as plt 

b,a = butter(5,0.25)

w,h = freqz(b,a)
plt.subplot(221)
plt.plot(w,abs(h))
plt.subplot(222)
plt.plot(np.real(h),np.imag(h))

bb,aa = invfreqz(h,w,5,5)
w,h2 = freqz(bb,aa)
plt.subplot(223)
plt.plot(w,abs(h2))
plt.subplot(224)
plt.plot(np.real(h2),np.imag(h2))
plt.show()
