import numpy as np
import pdb
from pylayers.antprop.slab import *
from pylayers.antprop.diffRT import *

#
# Metalic case : MacNamara Page 202
#
Nf=3
Nr=10
fGHz = np.linspace(1,10,Nf)
N = np.linspace(1,10,Nr)#320/180.
phi0 = np.linspace(0.01,2*np.pi-0.01,Nr)#40*np.pi/180.
phi = np.linspace(0.01,2*np.pi-0.01,Nr)
#phi = np.linspace(0,3*np.pi/2,10)
dm = MatDB()
mat0 = dm['METAL']
matN = dm['METAL']
si = 10000.*np.ones(Nr)
sd = 1.*np.ones(Nr)
plt.ion()
Ds,Dh = diff(fGHz,phi0,phi,si,sd,N,mat0,matN)
#plt.plot(phi*180/np.pi,np.real(Ds[0,0,:,0,0,0,0]),'b',label='Soft -')
#plt.plot(phi*180/np.pi,np.imag(Ds[0,0,:,0,0,0,0]),'b.',label='Soft -')
#plt.plot(phi*180/np.pi,np.real(Dh[0,0,:,0,0,0,0]),'r',label='Hard +')
#plt.plot(phi*180/np.pi,np.imag(Dh[0,0,:,0,0,0,0]),'r.',label='Hard +')
#plt.plot(phi*180/np.pi,np.real(D1[0,0,:,0,0,0,0]),'k',label='ISBinf')
#plt.plot(phi*180/np.pi,np.imag(D1[0,0,:,0,0,0,0]),'k',label='ISBinf')
#plt.plot(phi*180/np.pi,np.real(D2[0,0,:,0,0,0,0]),'g',label='ISBsup')
#plt.plot(phi*180/np.pi,np.imag(D2[0,0,:,0,0,0,0]),'g.',label='ISBsup')
#plt.plot(phi*180/np.pi,20*np.log10(np.abs(D3[0,0,:,0,0,0,0])),'c',label='RSBn')
#plt.plot(phi*180/np.pi,np.real(D4[0,0,:,0,0,0,0]),'m',label='RSBo')
#plt.plot(phi*180/np.pi,np.imag(D4[0,0,:,0,0,0,0]),'m.',label='RSBo')
# plt.legend()
# plt.grid()
# UI  = np.zeros(np.shape(Ds))
# UII = np.zeros(np.shape(Ds))
# UI[0,0,phi<(np.pi+phi0),0,0,0,0]=1
# UII[0,0,phi<(np.pi-phi0),0,0,0,0]=1.
# di  = np.sqrt(sd**2+si**2-2*si*sd*np.cos(phi[None,None,:,None,None,None,None]-phi0))
# dr  = np.sqrt(sd**2+si**2-2*si*sd*np.cos(phi[None,None,:,None,None,None,None]+phi0))
# ds  = si + sd
# #di  = np.cos(phi[None,None,:,None,None,None,None]-phi0)
# #dr  = np.cos(phi[None,None,:,None,None,None,None]+phi0)
# #ds  = sd
# Ei  = np.exp(-1j*2*np.pi*fGHz*di/0.3)*UI
# Ers = -np.exp(-1j*2*np.pi*fGHz*dr/0.3)*UII
# Erh = np.exp(-1j*2*np.pi*fGHz*dr/0.3)*UII
# Ed2 = D2*np.exp(-1j*2*np.pi*fGHz*ds/0.3)
# Ed4 = D4*np.exp(-1j*2*np.pi*fGHz*ds/0.3)
# #Ets = Ei+Ers+Ed1
# Ets = Ei + Ers + (Ed2-Ed4)
# Eth = Ei + Erh + (Ed2+Ed4)
# Ethf = Ei + Erh + Dh*np.exp(-1j*2*np.pi*fGHz*ds/0.3)
# Etsf = Ei + Ers + Ds*np.exp(-1j*2*np.pi*fGHz*ds/0.3)
# Ec2 = Ei + Ed2
# Ec4 = Ers - Ed4
# #Ec24 = Ed2+Ei+Ed4+Ers
# #Eth = Ei+Erh
# #plt.figure()
# plt.plot(phi*180/np.pi,np.abs(Ed2[0,0,:,0,0,0,0]),'g',label='ISBinf')
# plt.plot(phi*180/np.pi,np.abs(Ei[0,0,:,0,0,0,0]),'g')
# plt.plot(phi*180/np.pi,np.abs(Ec2[0,0,:,0,0,0,0]),'g')
# plt.plot(phi*180/np.pi,np.abs(Ed4[0,0,:,0,0,0,0]),'r',label='ISBinf')
# plt.plot(phi*180/np.pi,np.abs(Ers[0,0,:,0,0,0,0]),'r')
# plt.plot(phi*180/np.pi,np.abs(Ec4[0,0,:,0,0,0,0]),'r')
# #plt.plot(phi*180/np.pi,np.abs(Ec24[0,0,:,0,0,0,0]),'r')
# plt.figure()
# #plt.plot(phi*180/np.pi,20*np.log10(np.abs(Ets[0,0,:,0,0,0,0])),'b',linewidth=2.5)
# #plt.plot(phi*180/np.pi,20*np.log10(np.abs(Eth[0,0,:,0,0,0,0])),'r',linewidth=2.5)
# plt.plot(phi*180/np.pi,20*np.log10(np.abs(Etsf[0,0,:,0,0,0,0])),'b',linewidth=1.5,label='soft')
# plt.plot(phi*180/np.pi,20*np.log10(np.abs(Ethf[0,0,:,0,0,0,0])),'r',linewidth=1.5,label='hard')
# plt.plot(phi*180/np.pi,20*np.log10(np.abs(Ds[0,0,:,0,0,0,0])),'b.')
# plt.plot(phi*180/np.pi,20*np.log10(np.abs(Dh[0,0,:,0,0,0,0])),'r.')
# plt.ylim([-40,20])
# plt.xlim([0,320])
# plt.xlabel(u'Angle $\phi$')
# plt.ylabel(u'Magnitude (dB)')
# plt.title(u'$\\alpha=4^{\circ},\phi_0=55^{\circ},f= 3 GHz, sd=1m$',fontsize=14)
# plt.grid()
# #plt.plot(phi*180/np.pi,np.real(Ei[0,0,:,0,0,0,0]),'g')
# #plt.plot(phi*180/np.pi,20*np.log10(np.abs(Ets[0,0,:,0,0,0,0])),'g.')
# #plt.plot(phi*180/np.pi,20*np.log10(np.abs(Eth[0,0,:,0,0,0,0])),'r')
# #plt.plot(phi*180/np.pi,np.real(Eth[0,0,:,0,0,0,0]),'r.')
