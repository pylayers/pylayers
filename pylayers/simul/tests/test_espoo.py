from pylayers.simul.link import *
from pylayers.antprop.channel import *
from pylayers.antprop.antenna import *
from pylayers.util.geomutil import *

# parameters
fmin  = 31.8 # GHz
fmax  = 33.4         # GHz
bw    = fmax - fmin  # bandwidth in GHz
nf    = 10001           # sweep frequency points
fc    = (fmin + fmax)/2.
hms   = 1.65 # height of the MS
hbs   = 6    # height of the BS
fGHz  = np.linspace(fmin,fmax,nf)
fonts = 20 # fontsize

# Initialization of the Layout and the link.
DL      = DLink()
DL.fGHz = fGHz

L       = Layout('espoo.lay',bbuild=1)
DL.L    = L

# coordinates of the MS and the BS in the scene.
ms10 = np.array([188.2,199.5,hms])
ms11 = np.array([170.6,192,hms])
ms12 = np.array([208,208.3,hms])
ms13 = np.array([224,170.4,hms])
ms1  = np.array([197,189.5,hms])
ms2  = np.array([201.7,179.7,hms])
ms3  = np.array([185.4,171.1,hms])
ms4  = np.array([195.1,159.9,hms]) # np.array([198,161,hms])
ms5  = np.array([232.5,148.5,hms])
ms6  = np.array([176.5,179,hms])
ms7  = np.array([163.2,165.8,hms])
ms8  = np.array([148.4,154.9,hms])
ms9  = np.array([191.9,209.6,hms])
ms   = np.vstack((ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8,ms9,ms10,ms11,ms12,ms13))

bs   = np.array([220,185,hbs])

# coordinates of the MS and the BS antennas
DL.a  = ms12
DL.b  = bs

# MS antenna
DL.Aa = Antenna('Omni')
# BS antenna
DL.Ab = Antenna('aperture',fGHz=fGHz)

# orientation of the Horn antenna
# alpha = 0.37 rads (z axis) 
# beta  = 0 (y axis)
# gamma = - (40deg +  pi/2 (x axis) ;  40 deg  = 0.698 rad
DL.Tb = MEulerAngle(0.37,0,-(0.698+np.pi/2.))

plot_pos = False

if plot_pos:
	DL.L.show(L=DL.L)
	plt.axis('on')
	plt.plot(bs[0],bs[1],'or')
	for k in range(13):
	    plt.plot(ms[k][0],ms[k][1],'ob')
	    plt.annotate(xy=(ms[k][0],ms[k][1]),s='ms'+str(k+1),fontsize=fonts)
	plt.title('Layout of the Umi environment',fontsize=fonts)
	plt.show()


# DL.L._show3() # to see the scene of the layout
# DL._show3()   # to see the scene + antennas in the layout
# L.showG('s',labels=1) # enable to see the no of segments.

#cutoff : jusqu'ou on peut aller en profondeur dans le graph
#threshold  : 
#	proche de 1 = il va chercher les trajets equivalents du LOS 
# 	proche de 0 = il va chercher les trajets equivalents au NLOS  (ca va prendre plus de temps)

# ra_ceil_H only the ground reflection
DL.eval(force=1,cutoff=3,threshold=0.4,nD=1,ra_ceil_H=0)
#DL.C.cut(threshold_dB=90)

# angular range
phimin  = np.pi/4.     # 45 deg
phimax  = 3*np.pi/2.   # 270 deg
phistep = 5*np.pi/180. # 5 deg
phi = np.arange(phimin,phimax,phistep)

# angular frequency profile 
afp = DL.afp(phi)

# angular delay profile 
adp = afp.toadp()
# adp.imshow(cmap=cm.jet)


# plt.figure()
# adp.imshow(cmap=cm.jet)

# polarplot
# plt.figure()
adp.polarplot(vmin=-130,title='PADP of BS-MS1')#+str())
plt.show()
# DL.R.show(L=DL.L,rlist=DL.C.selected)

