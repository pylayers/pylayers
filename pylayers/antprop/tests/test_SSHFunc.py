import numpy as np 
import matplotlib.pyplot as plt
from pylayers.antprop.antssh import *
from pylayers.antprop.antenna import *
from sklearn.linear_model import Lasso
from scipy.optimize import fmin,minimize
plt.ion()
# reference antenna gain
A = Antenna('hplanesectoralhorn',fGHz=26)
A.eval()
# G est reel
G = A.G[:,:,0]
nth = G.shape[0]
nph = G.shape[1]
L  = 45
theta = np.linspace(0,np.pi,nth)
phi = np.linspace(0,2*np.pi,nph)
# Spherical Harmonics matrix
Y,idx = SSHFunc(L,theta,phi)
Ypinv = np.linalg.pinv(Y)
cg = np.dot(G.ravel(),Ypinv)
Gr = np.real(np.dot(cg,Y).reshape(nth,nph))
#plt.subplot(121)
#plt.imshow(np.abs(Y.T),cmap='jet')
#plt.axis('auto')
#plt.subplot(122)
#plt.imshow(np.angle(Y.T),cmap='jet')
#plt.axis('auto')
#
Gt = G[:,[0,45,90,135]]
Y2,idx = SSHFunc(L,theta,phi[[0,45,90,135]])

def fun(x,Y2=Y2,gt=Gt.ravel(),alpha=0.1):
    c  = x[::2]+1j*x[1::2]
    gr = np.real(np.dot(c,Y2))
    L2 = np.sqrt(np.dot(gr-gt,gr-gt))
    L1 = np.sum(np.abs(x))
    M  = L2+alpha*L1
    return(M)

Y2pinv = np.linalg.pinv(Y2)
cg2 = np.dot(Gt.ravel(),Y2pinv)
xg2 = np.array(zip(np.real(cg2),np.imag(cg2))).ravel()
xg = np.array(zip(np.real(cg),np.imag(cg))).ravel()
print "Start optimization"
ropt = minimize(fun,xg2,method='CG')
xopt = fmin(fun,xg2)
xopt = ropt.x
copt = xopt[::2]+1j*xopt[1::2]
Gr2 = np.real(np.dot(copt,Y).reshape(90,181))
plt.subplot()
plt.subplot(311)
plt.imshow(10*np.log10(np.abs(G)),cmap='jet',vmin=-40)
plt.colorbar()
plt.subplot(312)
plt.imshow(10*np.log10(np.abs(Gr)),cmap='jet',vmin=-40)
plt.colorbar()
plt.subplot(313)
plt.imshow(10*np.log10(np.abs(Gr2)),cmap='jet',vmin=-40,vmax=18)
plt.colorbar()
plt.show()
#alpha = 0.1
#lasso = Lasso(alpha=alpha)
#h=lasso.fit(Y2.T,Gt.ravel).predict(Y2.T)
