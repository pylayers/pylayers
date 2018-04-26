from pylayers.antprop.loss import *
import matplotlib.pyplot as plt
fGHz=2.4
p0=np.array(([0,0,20]))
p1=np.array(([0,1,20]))
p0=p0.reshape(3,1)
p1=p1.reshape(3,1)
TRF = [] #Two Ray model on flat earth
TRC = [] #Two Ray model on curved earth
PLoss=[]
for d in np.arange(1,10000,1):
    p1[1,:]=d
    TRF.append(two_rays_flatearth(p0[:,0],p1[:,0],fGHz,GtdB=0.,GrdB=0.,))
    TRC.append(two_rays_curvedearth(d,p0[2,:],p1[2,:],fGHz))
    PLoss.append(PL(fGHz, p0[:,0],p1[:,0], n=2.0, dB=True, d0=np.array([1])))
PLoss=np.array(PLoss)[:,0,0]
plt.semilogx(TRF,label='two-rays model flat earth')
plt.semilogx(TRC,label='two-rays model curved earth')
plt.semilogx(-PLoss,label='Path Loss')
plt.legend()
plt.show()
