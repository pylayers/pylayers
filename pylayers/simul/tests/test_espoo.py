from pylayers.simul.link import *
from pylayers.antprop.channel import *
from pylayers.antprop.antenna import *
fGHz = np.linspace(32.6-1.6/2,32.6+1.6/2,2000)
DL=DLink()
DL.fGHz=fGHz
L=Layout('espoo.ini')
DL.L=L
DL.a=np.array([150,150,1.75])
DL.b=np.array([220,185,6])
DL.Aa = Antenna('Omni')
DL.Ab = Antenna('aperture')
#DL.eval(force=1,ra_ceil_H=0)
DL.eval(force=True,bt=False,cutoff=5,threshold=0.7,ra_vectorized=True,debug=True,diffraction=True)
#DL.C.cut(threshold_dB=80)
# angular range 
phi = np.arange(np.pi/2,3*np.pi/2,5*np.pi/180)-np.pi/3.
# angular frequency profile 
#afp=DL.afp(phi)
# angular delay profile 
#adp=afp.toadp()
#adp.imshow(cmap=cm.jet)
#plt.figure()
#adp.imshow(cmap=cm.jet)
#plt.figure()
#adp.polarplot(vmin=-160)
#DL.R.show(L=DL.L,rlist=DL.C.selected)
