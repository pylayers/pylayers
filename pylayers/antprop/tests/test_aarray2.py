from pylayers.antprop.aarray import *
import matplotlib.pyplot as plt 
print('--------------')
print('antprop/test_aarray2.py')
print('--------------')
#
# This is an array obtained by juxtaposition of 2 uniform array
#   - the first along the z axis
#   - the second along the x axis
#  Antenna spacing is 4cm 0.04m
#
fGHz = 6
lamda = 0.3/6
UA1 = UArray(N=[10,1,1],dm=[0.04,0.0,0])
UA2 = UArray(N=[1,1,10],dm=[0,0.0,0.04])
#A1.eval()
#f1,a1 = UA1.plotG()
#f1.savefig('ArrayDiag.png')
#f2,a2 = UA1.show()
#f2.savefig('ArrayConfig.png')
#
#Nrays = 100
#theta = np.pi*np.random.rand(Nrays)
#phi = 2*np.pi*np.random.rand(Nrays)
#A1.eval(grid=False,th=theta,ph=phi)
#plt.figure()
#plt.imshow(np.angle(A1.Ft[:,:,0]),cmap='jet',interpolation='nearest')
#plt.colorbar()
#plt.figure()
#plt.imshow(np.angle(A1.Fp[:,:,0]),cmap='jet',interpolation='nearest')
#plt.colorbar()
##plt.axis('auto')
