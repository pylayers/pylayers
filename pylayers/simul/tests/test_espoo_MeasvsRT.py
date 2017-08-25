from pylayers.simul.link import *
from pylayers.antprop.channel import *
from pylayers.antprop.antenna import *
from pylayers.util.geomutil import *


plt.close('all')


# parameters settings RT
fmin       = 31.8 # GHz
fmax       = 33.4         # GHz
BGHz       = fmax - fmin  # bandwidth in GHz
Nf         = 10001           # sweep frequency points
fcGHz      = (fmin + fmax)/2.
hms        = 1.65 # height of the MS
hbs        = 6    # height of the BS
fGHz       = np.linspace(fmin,fmax,Nf)
tauns_meas = np.linspace(0,1/(fGHz[1]-fGHz[0]),Nf)
fonts      = 20 # fontsize
thres      = 0.1
cutoff     = 1

d = {}

# Initialization of the Layout and the link.
L  = Layout('espoo.lay',bbuild=1)
DL = DLink(L=L,fGHz=fGHz,outdoor=True)
DL.fGHz = fGHz

# coordinates of the MS and the BS in the scene.
ms1  = np.array([197,189.8,hms]) # GOOD MATCH!!! peak search aki losdelay = 79.3750 ns losang = 215 deg - gain omni = 11 dB
ms2  = np.array([201.8,179.4,hms]) # GOOD MATCH!!! peak search aki losdelay = 65 ns losang = 180 deg - gain omni = 3.51 dB
ms3  = np.array([186,169.9,hms]) # GOOD MATCH!!! peak search aki losdelay = 125 ns losang = 180 deg - gain omni = 8.79 dB
ms4  = np.array([197.8,157.6,hms]) # QUITE GOOD MATCH!!! peak search aki losdelay = 118.75 ns losang = 150 deg - gain omni = 12.857 dB
ms5  = np.array([232.5,148.5,hms]) # GOOD MATCH!!! peak search aki losdelay = 129.375 ns losang = 95 deg - gain omni = 4.55 dB
ms6  = np.array([177.2,177.9,hms]) # GOOD MATCH!!! peak search aki losdelay = 146.875 ns losang = 190 deg - gain omni =  7.19 dB
ms7  = np.array([162.8,167,hms]) # GOOD MATCH!!! peak search aki losdelay = 200.625 ns losang = 185 deg - gain omni = 8.46 dB
ms8  = np.array([148.3,155,hms]) # GOOD MATCH!!! peak search aki losdelay = 259.375 ns losang = 180 deg - gain omni = 8.7 dB
ms9  = np.array([191.9,209.5,hms]) #GOOD MATCH!!! peak search aki losdelay = 125 ns losang = 240 deg - gain omni = 1.76 dB
ms10 = np.array([188.1,199.5,hms]) #GOOD MATCH!!! peak search aki losdelay = 117.5 ns losang = 230 deg - gain omni = 2.9 dB
ms11 = np.array([170.5,191.5,hms]) #QUITE GOOD MATCH!!! peak search aki losdelay = 166.875 ns losang = 215 deg - gain omni = 9.87 dB
ms12 = np.array([206,205.9,hms]) #QUITE GOOD MATCH!!! peak search aki losdelay = 88.75 ns losang = 255 deg - gain omni = 1.75 dB
ms13 = np.array([224,170.4,hms]) #GOOD MATCH!!! peak search aki losdelay = 52.5 ns losang = 95 deg - gain omni = 5.76 dB
ms   = np.vstack((ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8,ms9,ms10,ms11,ms12,ms13))
bs   = np.array([220,185,hbs])

# coordinates of the MS and the BS antennas
DL.a  = ms6
DL.b  = bs
numsnap = 6

# MS antenna
DL.Aa = Antenna(typ='Omni',param = { 'pol' : 't', 'GmaxdB': 0})
DL.Ab = Antenna(typ='Omni',param = { 'pol' : 't', 'GmaxdB': 0})
# BS antenna

# DL.Ab = Antenna(typ='aperture', fGHz = np.array([32.6]) ,param = { 'HPBW_x_deg': 40, 'HPBW_y_deg':10,
#     'Gfactor':27000, 'fcGHz': fcGHz, 'polar':'x'})

# orientation of the Horn antenna
# alpha = 0.37 rads (z axis) 
# beta  = 0 (y axis)
# gamma = - (40deg +  pi/2 (x axis) ;  40 deg  = 0.698 rad
# DL.Tb = MEulerAngle(0.37,0,-(0.698+np.pi/2.))

plot_pos = 0

plt.ion()

if plot_pos:
	DL.L.show(L=DL.L)
	plt.axis('on')
	plt.plot(bs[0],bs[1],'or')
	for k in range(13):
	    plt.plot(ms[k][0],ms[k][1],'ob')
	    plt.annotate(xy=(ms[k][0],ms[k][1]),s='ms'+str(k+1),fontsize=fonts)
	plt.title('Layout of the Umi environment',fontsize=fonts)
	plt.show()


print "---------------------------"
print "Beginning of the RT & MEAS."
print "---------------------------"
print "Freq (GHz)       :",fcGHz
print "BW   (GHz)       :",BGHz
print "Freq. points     :",Nf
print "threshold        :",thres
print "cutoff           :",cutoff
print "BS-MS            :",numsnap
print "------------------------"


DL.eval(force = 1,bt=1,cutoff = cutoff,threshold=thres,nD=1,diffraction=False)
#DL.C.cut(threshold_dB=90)

# angular range
phimin  = np.pi/4.     # 45 deg
phimax  = 3*np.pi/2.   # 270 deg
phistep = 5*np.pi/180. # 5 deg
phi = np.arange(phimin,phimax,phistep)
phideg = (phi*180)/np.pi

# angular frequency profile 
tilt = 0.698132 # 40 deg #0.698132 # 0.51; 0.53 ; augmenter = tilt bas et diminuer = tilt haut
phis = 0.37 - phi
beta = 0
gamma = -( tilt + np.pi/2.) 
#DL.Tb = MEulerAngle(phi[26],beta=beta,gamma=gamma)
afp = DL.afp(phi=phis,beta=beta,gamma=gamma)
#afp = DL.afp(phi=[0,0],beta=beta,gamma=0)

# # angular delay profile 
adp = afp.toadp() # (angle x Nf)

padp_rt   = adp.y # (angle x Nf)
padpdB_rt = 20*np.log10(np.abs(padp_rt)) # (angle x Nf)

# deriving pdp from padp: pdp = marginal integral of padp
pdp_rt = np.mean(adp.y*np.conj(adp.y),axis=0) #mean over angles.
#pdp_rt = np.mean(np.abs(adp.y)**2,axis=0) # mean over angles.
pdpdB_rt = 10*np.log10(np.abs(pdp_rt))


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

trf = amp*np.exp(1j*ang*np.pi/180)*cal_trf[None,:]#*window
cir = np.fft.ifft(trf,axis=1)
d[numsnap] = {'cal':trf,'cir':cir}

padp_meas   = np.abs(d[numsnap]['cir'])       # (ang,Nf)
padpdB_meas = 20*np.log10(padp_meas)
pdp_meas   = np.mean(padp_meas*np.conj(padp_meas),axis=0) # mean over angles
pdpdB_meas = 10*np.log10(np.abs(pdp_meas))


polar    = 0
plotpdp  = 1
plotpadp = 1

if polar:
	adp.polarplot(vmin=-130,title='PADP of BS-MS1')#+str())

if plotpdp:

	figpdp = plt.figure(figsize=(10,8))

	tauns_meas = np.linspace(0,1/(fGHz[1]-fGHz[0]),Nf)

	d_m = tauns_meas * 0.3
	GtdB = 2
	#GrdB = 5.38
	
	# ouverture mi puissance et pas angulaire
	# GmdB : mean gain depending on the aperture at mid power (19 - 3 dB x 2) nd the angular power step
	GrdB = 10*np.log10((10**1.9 + 10**1.6 + 10**1.6)/46.)
	GmdB = GrdB + GtdB
	PLdB = -32.4 - 20*np.log10(32.6) - 20*np.log10(d_m) + GmdB

	#pdprt1        = DL.C.Ctt.energy(Friis=1)
	pdprt_dB     = 10*np.log10(DL.C.Ctp.energy(Friis=1))
	uu           = np.where(pdprt_dB<-140)
	pdprt_dB[uu] = 0

	# pdprt        = DL.H.energy()[:,0,0]
	# pdprt_dB     = 10*np.log10(np.abs(pdprt))
	# uu           = np.where(pdprt_dB<-140)
	# pdprt_dB[uu] = 0
	
	# plotting
	# plt.plot(DL.H.taud,pdprt_dB,'or',label='from RT-NRJ')
	plt.plot(tauns_meas,pdpdB_rt,'k',label='RT: BS-MS'+str(numsnap),linewidth=2)
	plt.plot(tauns_meas,pdpdB_meas,'b',label='MEAS: BS-MS'+str(numsnap),linewidth=2)
	#plt.semilogx(tauns_meas,PLdB,'r',label='FRIIS',linewidth=2)
	plt.title('Meas vs RT: link BS-MS'+ ' '+str(numsnap),fontsize=fonts)
	plt.xlim(0,700) #plt.xlim(0,700)
	plt.ylim(-170,-60)#plt.ylim(-170,-80)
	plt.xlabel('Propagation delay [ns]',fontsize=fonts)
	plt.ylabel('CIR [dB]',fontsize=fonts)
	plt.xticks(fontsize=fonts)
	plt.yticks(fontsize=fonts)
	plt.grid()
	plt.legend(fontsize=fonts)
	# DL.R.show(L=DL.L,rlist=DL.C.selected)
	savefig = '/home/dialounke/slides_presentation/Espoo/RT_figs/'
	pdpfign = 'BSMS'+ str(numsnap)+'_'+'pdp.png'
	#figpdp.savefig(pdpfign)



if plotpadp:

	# Meas.
	figpadpmeas = plt.figure(figsize=(10,8))
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
	figpadprt = plt.figure(figsize=(10,8))
	imin = 0
	imax = 1000
	phideg = (phi*180)/np.pi
	plt.imshow(padpdB_rt[:,imin:imax].T,origin='lower',cmap='jet',extent=(phideg[0],phideg[-1],tauns_meas[imin],tauns_meas[imax]),vmax=-65,vmin=-120,interpolation='nearest')
	plt.axis('auto')
	plt.xlabel('Angle (deg)',fontsize=fonts)
	plt.ylabel('Propagation delay (ns)',fontsize=fonts)
	plt.xticks(fontsize=fonts)
	plt.yticks(fontsize=fonts)
	plt.title('RT: link BS-MS'+str(numsnap),fontsize=fonts)
	plt.colorbar()

	padpmfign = 'BSMS'+ str(numsnap)+'_'+'padp_meas.png'
	#figpadpmeas.savefig(padpmfign)

	padprfign = 'BSMS'+ str(numsnap)+'_'+'padp_rt.png'
	#figpadprt.savefig(padprfign)


umpdp     = np.where(pdpdB_meas==np.max(pdpdB_meas))[0]
urtpdp    = np.where(pdpdB_rt==np.max(pdpdB_rt))[0]
umpadp    = np.where(padpdB_meas==np.max(padpdB_meas))[0]
urtpadp   = np.where(padp_rt==np.max(padp_rt))[0]

print "-------------------------------------------"
print "PADP Meas. [dB]     :",padpdB_meas.max()
print "PADP RT [dB]        :",padpdB_rt.max()
print "PADP Offset [dB]    :",padpdB_meas.max() - padpdB_rt.max() 
print "-------------------------------------------"
print "PDP Meas. [dB]      :",pdpdB_meas.max()
print "PDP RT [dB]         :",pdpdB_rt.max()
print "PDP Offset [dB]     :",pdpdB_meas.max() - pdpdB_rt.max()
print "-------------------------------------------"
print "peak los Meas. [ns] :",tauns_meas[umpdp]
print "peak los RT    [ns] :",tauns_meas[urtpdp]
print "Offset         [ns] :",tauns_meas[urtpdp] - tauns_meas[umpdp]
print "-------------------------------------------"
print "peak ang Meas. [ns] :",phi[umpadp]*180/np.pi
print "peak ang RT    [ns] :",phi[urtpadp]*180/np.pi
print "Offset         [ns] :",((phi[urtpadp] - phi[umpadp])*180)/np.pi