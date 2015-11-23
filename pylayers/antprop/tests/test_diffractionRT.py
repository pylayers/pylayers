import numpy as np
import pdb
from pylayers.antprop.slab import *
from pylayers.antprop.diffRT import *

#
# Metalic case : MacNamara Page 202
#
Nf=1
Nr=800
fGHz = np.linspace(1,10,Nf)
N = 320/180.*np.ones(Nr)#320/180.
phi0 = np.ones(Nr)#np.linspace(0.01,2*np.pi-0.01,Nr)#40*np.pi/180.
phi = np.linspace(0.01,2*np.pi-0.01,Nr)
#phi = np.linspace(0,3*np.pi/2,10)
dm = MatDB()
mat0 = dm['METAL']
matN = dm['METAL']
si = 10000.*np.ones(Nr)
sd = 1.*np.ones(Nr)
plt.ion()
Ds,Dh,D1,D2,D3,D4 = diff(fGHz,phi0,phi,si,sd,N,mat0,matN,debug=True)
plt.legend()
plt.grid()
UI  = np.zeros(np.shape(Ds))
UII = np.zeros(np.shape(Ds))
UI[0,phi<(np.pi+phi0)]=1
UII[0,phi<(np.pi-phi0)]=1.
di  = np.sqrt(sd**2+si**2-2*si*sd*np.cos(phi-phi0))
dr  = np.sqrt(sd**2+si**2-2*si*sd*np.cos(phi+phi0))
ds  = si + sd
#di  = np.cos(phi[None,None,:,None,None,None,None]-phi0)
#dr  = np.cos(phi[None,None,:,None,None,None,None]+phi0)
#ds  = sd
Ei  = np.exp(-1j*2*np.pi*fGHz[:,None]*di/0.3)*UI
Ers = -np.exp(-1j*2*np.pi*fGHz[:,None]*dr/0.3)*UII
Erh = np.exp(-1j*2*np.pi*fGHz[:,None]*dr/0.3)*UII
Ed2 = D2*np.exp(-1j*2*np.pi*fGHz[:,None]*ds/0.3)
Ed4 = D4*np.exp(-1j*2*np.pi*fGHz[:,None]*ds/0.3)
#Ets = Ei+Ers+Ed1
Ets = Ei + Ers + (Ed2-Ed4)
Eth = Ei + Erh + (Ed2+Ed4)
Ethf = Ei + Erh + Dh*np.exp(-1j*2*np.pi*fGHz[:,None]*ds/0.3)
Etsf = Ei + Ers + Ds*np.exp(-1j*2*np.pi*fGHz[:,None]*ds/0.3)
Ec2 = Ei + Ed2
Ec4 = Ers - Ed4
#Ec24 = Ed2+Ei+Ed4+Ers
#Eth = Ei+Erh
#plt.figure()
plt.plot(phi*180/np.pi,np.abs(Ed2[0,:]),'g',label='ISBinf')
plt.plot(phi*180/np.pi,np.abs(Ei[0,:]),'g')
plt.plot(phi*180/np.pi,np.abs(Ec2[0,:]),'g')
plt.plot(phi*180/np.pi,np.abs(Ed4[0,:]),'r',label='ISBinf')
plt.plot(phi*180/np.pi,np.abs(Ers[0,:]),'r')
plt.plot(phi*180/np.pi,np.abs(Ec4[0,:]),'r')
#plt.plot(phi*180/np.pi,np.abs(Ec24[0,0,:,0,0,0,0]),'r')
plt.figure()
#plt.plot(phi*180/np.pi,20*np.log10(np.abs(Ets[0,0,:,0,0,0,0])),'b',linewidth=2.5)
#plt.plot(phi*180/np.pi,20*np.log10(np.abs(Eth[0,0,:,0,0,0,0])),'r',linewidth=2.5)
plt.plot(phi*180/np.pi,20*np.log10(np.abs(Etsf[0,:])),'b',linewidth=1.5,label='soft')
plt.plot(phi*180/np.pi,20*np.log10(np.abs(Ethf[0,:])),'r',linewidth=1.5,label='hard')
plt.plot(phi*180/np.pi,20*np.log10(np.abs(Ds[0,:])),'b.')
plt.plot(phi*180/np.pi,20*np.log10(np.abs(Dh[0,:])),'r.')
plt.ylim([-40,20])
plt.xlim([0,320])
plt.xlabel(u'Angle $\phi$')
plt.ylabel(u'Magnitude (dB)')
plt.title(u'$\\alpha=4^{\circ},\phi_0=55^{\circ},f= 3 GHz, sd=1m$',fontsize=14)
plt.grid()