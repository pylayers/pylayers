from pylayers.antprop.aarray import *
from pylayers.antprop.antenna import *

# A   = AntArray(typant='Gauss') ok
# A   = AntArray(typant='Hertz') ok 
# A   = AntArray(typant='Huygens') ok
# A   = AntArray(typant='Omni') ok
# A   = AntArray(typant='Array') doesn't work!
# A   = AntArray(typant='wireplate') doesn't work!
# A   = AntArray(typant='sh3') doesn't work!
# A   = AntArray(typant='vsh3') doesn't work
# A   = AntArray(typant='3gpp') doesn't work
f = 2.4
ld = 0.3/f
da = ld/2
A   = AntArray(typant='Huygens',dm=[da,0,0],N=[10,1,1])
# A1 = AntArray(dm=[0.0025,0.0025,0],N=[8,2,1],typant='Huygens')
npts = 10
#w   = 0.5*np.random.randn(npts)+1j*0.5*np.random.randn(npts) # (npts,)
w   = np.ones(npts)
wn  = w[:,None] # (npts,1)
k = np.arange(npts)
#w = np.exp(1j*k[:,None]*0.08*A.fGHz[None,:])
A.w = wn # (npts,1)
A.eval()
A.plotG()
A.show3()
# A.plot3d() doesn't work
