# -*- coding: latin1 -*-

from pylayers.simul.link import *
import pylayers.signal.waveform as wvf
import pdb
import warnings 

fcGHz=5
WMHz = 3000
Nf = 400 
fGHz = np.linspace(fcGHz-WMHz*0.5e-3,fcGHz+WMHz*0.5e-3,Nf)




L=Layout('defstr.lay',bgraphs=1)
L=Layout('espoo.lay',bgraphs=1)
print('DL has removed airwalls') 
DL=DLink(L=L,fGHz=fGHz,outdoor=False,applywav=True)
DL.ca= 2
DL.cb= 5
DL.ca= 106 # uncomment to test with espoo.lay
DL.cb= 155 # uncomment to test with espoo.lay
DL.a = DL.a+0.8
DL.b = DL.b+0.8
DL.Aa=Antenna(typ='Omni')
DL.Ab=Antenna(typ='Omni')
DL.eval(verbose=True,force=True,bt=False,cutoff=4,threshold=0.8,ra_vectorized=True,nD=1)

# print('DL2 still has airwalls') 
# DL2=DLink(L=L,fGHz=fGHz,outdoor=False,applywav=True)
# # DL2.ca= 106 # uncomment to test with espoo.lay
# # DL2.cb= 155 # uncomment to test with espoo.lay
# DL2.ca= 2
# DL2.cb= 5
# DL2.a = DL2.a+0.8
# DL2.b = DL2.b+0.8
# DL2.Aa=Antenna(typ='Omni')
# DL2.Ab=Antenna(typ='Omni')
# DL2.eval(verbose=True,force=True,bt=False,cutoff=4,threshold=0.8,ra_vectorized=True,nD=1,rm_aw=False)

# # ray index have been modified due to airwall removal
# # DL.R._rayidx_aw keep trance of this change
# u_aw = DL.R._rayidx_aw

# print (np.where(DL.C.Cpp.y-DL2.C.Cpp.y[u_aw]>1e-5))
# print (np.where(DL.C.Ctt.y-DL2.C.Ctt.y[u_aw]>1e-5))
# print (np.where(DL.C.Ctp.y-DL2.C.Ctp.y[u_aw]>1e-5))
# print (np.where(DL.C.Cpt.y-DL2.C.Cpt.y[u_aw]>1e-5))
