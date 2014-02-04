import numpy as np
import matplotlib.pyplot as plt
fGHz=np.arange(2,11,0.1)
fcenter = fGHz[len(fGHz)/2]
dist = np.arange(1,9,0.4)
tEtt = np.load('tEtt.npy')
tEtp = np.load('tEtp.npy')
tEpt = np.load('tEpt.npy')
tEpp = np.load('tEpp.npy')
plt.plot(dist,10*np.log10(tEtt))
plt.plot(dist,10*np.log10(tEtp))
plt.plot(dist,10*np.log10(tEpt))
plt.plot(dist,10*np.log10(tEpp))
plt.xlabel('Tx Rx distance (meters)')
plt.ylabel('Attenuation (dB)')
plt.legend((u'$\\theta\\theta$',u'$\\theta\phi$',u'$\phi\\theta$',u'$\phi\phi$'))
plt.show()
