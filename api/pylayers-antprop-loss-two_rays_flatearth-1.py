from pylayers.antprop.loss import *
NPT=10000
x=np.array([0,0,8])
x=x.reshape(3,1)
y = np.ones((3,NPT))
y[0,:]=0
y[1,:]=np.arange(NPT)
y[2,:]=2
g0=1
g1=1
fGHz=2.4
PL2R=two_rays_flatearth(p0=x,p1=y,fGHz=fGHz,GtdB=g0,GrdB=g1)
PL1R = PL(fGHz,x,y,2)
plt.semilogx(PL2R,label='two-ray model')
plt.semilogx(-PL1R[0,:],label='one slope model')
plt.axis([10,NPT,-150,-50])
plt.title('Loss 2-rays model vs one slope model')
plt.xlabel('distance (m)')
plt.ylabel('Loss Pr/Pt (dB)')
plt.legend()
plt.show()

d=np.arange(1,1000)
PL2Rd = two_rays_flatearth(d=d,ht=np.array([5]),hr=np.array([10]),fGHz=fGHz,GtdB=g0,GrdB=g1)
plt.semilogx(PL2Rd,label='two-ray model')
plt.semilogx(-PL1R[0,:],label='one slope model')
plt.axis([10,NPT,-150,-50])
plt.title('Loss 2-rays model vs one slope model')
plt.xlabel('distance (m)')
plt.ylabel('Loss Pr/Pt (dB)')
plt.legend()
plt.show()
