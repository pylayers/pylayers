from pylayers.antprop.aarray import *
import matplotlib.pyplot as plt 
print '--------------'
print 'antprop/test_aarray.py'
print '--------------'
#
# This is a uniform array (10 x 10) 
#
A = AntArray(tarr='UA',N=[10,10,1],dm=[0.5,0.5,0],fGHz=0.3)
A.eval()
f1,a1 = A.plotG()
f1.savefig('ArrayDiag.png')
f2,a2 = A.show()
f2.savefig('ArrayConfig.png')

Nrays = 100
theta = np.pi*np.random.rand(Nrays)
phi = 2*np.pi*np.random.rand(Nrays)
A.eval(grid=False,th=theta,ph=phi)
plt.figure()
plt.imshow(np.angle(A.Ft[:,:,0]),cmap='jet',interpolation='nearest')
plt.colorbar()
plt.figure()
plt.imshow(np.angle(A.Fp[:,:,0]),cmap='jet',interpolation='nearest')
plt.colorbar()
#plt.axis('auto')
