from pylayers.simul.link import *
from pylayers.antprop.channel import *
from pylayers.antprop.antenna import *
from pylayers.util.geomutil import *

# parameters settings RT
fmin  = 31.8 # GHz
fmax  = 33.4         # GHz
bw    = fmax - fmin  # bandwidth in GHz
nfrt  = 10001           # sweep frequency points
fc    = (fmin + fmax)/2.
hms   = 1.65 # height of the MS
hbs   = 6    # height of the BS
fGHzrt  = np.linspace(fmin,fmax,nfrt)
fonts = 20 # fontsize

# Parameters settings meas.
fcGHz = 28.5
BGHz  = 1.6
Nf    = 10001
fGHzm  = np.linspace(fcGHz-BGHz*0.5,fcGHz+BGHz*0.5,Nf)
tauns_meas = np.linspace(0,1/(fGHzm[1]-fGHzm[0]),Nf)
fGHz  = np.array([32.6])
thres  = 0.1 # for the RT
cutoff = 3   # for the RT
numsnap = 10
d = {}


print "---------------------------"
print "Beginning of the RT & MEAS."
print "---------------------------"
print "Freq (GHz)       :",fcGHz
print "BW   (GHz)       :",BGHz
print "Freq. points     :",Nf
print "threshold        :",thres
print "cutoff           :",cutoff
print "BS-MS            ",numsnap
print "------------------------"

# Initialization of the Layout and the link.
L       = Layout('espoo_journal.ini',bbuild=1,bdiffraction=True)
DL = DLink(L=L,fGHz=fGHzrt,outdoor=True)
DL.fGHz = fGHzrt

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
DL.a  = ms10
DL.b  = bs

# MS antenna
DL.Aa = Antenna(typ='Omni',param = { 'pol' : 't', 'GmaxdB': 2 })
# BS antenna
DL.Ab = Antenna('aperture')

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


DL.eval(force=1,bt=1,cutoff=3,threshold=thres,nD=1,diffraction=False)
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

#
# Measurement
# 

#directory   = '/home/dialounke/code/Aalto/32GHz'
directory   = '/home/dialounke/code/pylayers/pylayers/measures/Aalto'
pathdata    = '/all_dams_data/'
calibration = '/B2B/B2B.mat'
C           = io.loadmat(directory+calibration)
cal_trf     = C['cal_trf'][:,0]
filename    = directory+pathdata+'Tx'+str(numsnap)+'Rx1_10.txt'
# import ipdb
# ipdb.set_trace()
D = np.loadtxt(filename,skiprows=2)
amp = D[:,2::2]
ang = D[:,3::2]
rot = (D[:,0]-360)*np.pi/180.

rotdeg =  -(rot * 180.)/ np.pi - 45
rot_step = rot[1]-rot[0]
window = np.ones(Nf)

trf = amp*np.exp(1j*ang*np.pi/180)*cal_trf[None,:]*window
cir = np.fft.ifft(trf,axis=1)
d[numsnap] = {'cal':trf,'cir':cir}

padp_meas   = np.abs(d[numsnap]['cir'])       # (ang,Nf)
padpdB_meas = 20*np.log10(padp_meas)

pdp_meas   = np.sum(padp_meas,axis=0)/len(padp_meas)         # (Nf,)
pdpdB_meas = 20*np.log10(np.abs(pdp_meas))

padp_meas   = np.abs(d[numsnap]['cir'])   # (ang,Nf)
padpdB_meas = 20*np.log10(padp_meas)
pdp_meas   = np.sum(padp_meas,axis=0)/len(padp_meas)     # (Nf,)
pdpdB_meas = 20*np.log10(np.abs(pdp_meas))


polar = 0
plotpdp = 1
plotpadp = 1

plt.ion()

if polar:
	adp.polarplot(vmin=-130,title='PADP of BS-MS1')#+str())

if plotpdp:

	# RT
	cir = DL.H.getcir(BWGHz=BGHz,Nf=Nf,fftshift=1)
	cir.plot()

	# Meas
	plt.plot(tauns_meas,pdpdB_meas,label='MEAS: BS-MS'+str(numsnap),linewidth=2)
	
	plt.title('Meas vs RT: link BS-MS'+ ' '+str(numsnap))
	plt.xlim(0,500)
	plt.ylim(-180,-70)
	plt.xlabel('Propagation delay [ns]',fontsize=fonts)
	plt.ylabel('CIR [dB]',fontsize=fonts)
	plt.xticks(fontsize=fonts)
	plt.yticks(fontsize=fonts)
	plt.grid()
	plt.legend(fontsize=fonts)	

	# DL.R.show(L=DL.L,rlist=DL.C.selected)


if plotpadp:

	# Meas.
	plt.figure()
	imin = 0
	imax = 1000
	plt.imshow(padpdB_meas[:,imin:imax].T,origin='lower',cmap='jet',extent=(rotdeg[-1],rotdeg[0],tauns_meas[imin],tauns_meas[imax]),vmax=-65,vmin=-120,interpolation='nearest')
	plt.axis('auto')
	plt.colorbar()
	plt.ylabel('Propagation delay [ns]',fontsize=fonts)
	plt.xlabel('Angle[deg]',fontsize=fonts)
	plt.title(' MEAS: BS-MS'+str(numsnap),fontsize=fonts)
	plt.xticks(fontsize=fonts)
	plt.yticks(fontsize=fonts)
	plt.show()

	# RT
	plt.figure()
	imin = 0
	imax = 1000
	phideg = (phi*180)/np.pi
	#plt.imshow(20*np.log10(np.abs(adp.y[:,imin:imax].T)),origin='lower',cmap='jet',extent=(phideg[0],phideg[-1],tauns_meas[imin],tauns_meas[imax]),vmax=-65,vmin=-120,interpolation='nearest')
	plt.imshow(20*np.log10(np.abs(adp.y[:,imin:imax].T)),origin='lower',cmap='jet',extent=(phideg[0],phideg[-1],tauns_meas[imin],tauns_meas[imax]),interpolation='nearest')
	plt.axis('auto')
	plt.xlabel('Angle (deg)',fontsize=fonts)
	plt.ylabel('Propagation delay (ns)',fontsize=fonts)
	plt.xticks(fontsize=fonts)
	plt.yticks(fontsize=fonts)
	plt.title('RT: link BS-MS'+str(numsnap),fontsize=fonts)
	plt.colorbar()



