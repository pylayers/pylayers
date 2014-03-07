import numpy as np
import pylayers.util.plotutil as plu
x = np.arange(0,10,.01)
z1 = np.cos(2*x)*np.sin(10*x) + 1j * np.cos(3*x)*np.sin(11*x)
z2 = np.cos(3*x)*np.sin(11*x) + 1j * np.cos(4*x)*np.sin(10*x)
z3 = np.cos(4*x)*np.sin(12*x) + 1j * np.cos(5*x)*np.sin(12*x)
z4 = np.cos(5*x)*np.sin(13*x) + 1j * np.cos(6*x)*np.sin(13*x)
y = np.vstack((z1,z2,z3,z4))
plu.mulcplot(x,y)
plt.show()
plu.mulcplot(x,y,typ='v',)
plt.show()
plu.mulcplot(x,y,typ='r',ncol=2,nlin=2,color='k',linewidth=2)
plt.show()
