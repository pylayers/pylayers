# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>



from pylayers.antprop.slab import *
from numpy import *


sl = SlabDB()



sl.mat.add(name='TESS-p50',cval=3+0j,sigma=0.06,typ='epsr')



sl.add(name='TESS-p50-5cm',lname=['TESS-p50'],lthick=[0.05])
sl.add(name='TESS-p50-10cm',lname=['TESS-p50'],lthick=[0.10])
sl.add(name='TESS-p50-15cm',lname=['TESS-p50'],lthick=[0.15])

# <markdowncell>

# These Tessereau page 50 



fGHz=4
theta = np.arange(0,pi/2,0.01)



sl['TESS-p50-5cm'].ev(fGHz,theta,compensate=True)
sl['TESS-p50-10cm'].ev(fGHz,theta,compensate=True)
sl['TESS-p50-15cm'].ev(fGHz,theta,compensate=True)


# by default var='a' and kv = 0 
fig,ax = sl['TESS-p50-5cm'].plotwrt(color='k',ncol=2,nlin=2,labels=[''])
fig,ax = sl['TESS-p50-10cm'].plotwrt(color='k',ncol=2,nlin=2,labels=[''],linestyle='dashed',fig=fig,ax=ax)
fig,ax = sl['TESS-p50-15cm'].plotwrt(color='k',ncol=2,nlin=2,labels=[''],linestyle='dashdot',fig=fig,ax=ax)



fGHz = np.arange(2,16,0.1)
theta = 0 

sl['TESS-p50-5cm'].ev(fGHz,theta,compensate=False)
sl['TESS-p50-10cm'].ev(fGHz,theta,compensate=False)
sl['TESS-p50-15cm'].ev(fGHz,theta,compensate=False)
    
fig,ax = sl['TESS-p50-5cm'].plotwrt('f',coeff='T',types=['ru'],labels=[''],color='k')
fig,ax = sl['TESS-p50-10cm'].plotwrt('f',coeff='T',types=['ru'],labels=[''],color='k',linestyle='dashed',fig=fig,ax=ax)
fig,ax = sl['TESS-p50-15cm'].plotwrt('f',coeff='T',types=['ru'],labels=[''],color='k',linestyle='dashdot',fig=fig,ax=ax)

sl['TESS-p50-5cm'].ev(fGHz,theta,compensate=True)
sl['TESS-p50-10cm'].ev(fGHz,theta,compensate=True)
sl['TESS-p50-15cm'].ev(fGHz,theta,compensate=True)

fig,ax = sl['TESS-p50-5cm'].plotwrt('f',coeff='T',types=['ru'],labels=['5cm compensated',''],color='r',fig=fig,ax=ax)
fig,ax = sl['TESS-p50-10cm'].plotwrt('f',coeff='T',types=['ru'],labels=['10cm compensated',''],color='r',linestyle='dashed',fig=fig,ax=ax)
fig,ax = sl['TESS-p50-15cm'].plotwrt('f',coeff='T',types=['ru'],labels=['15cm not compensated',''],color='r',linestyle='dashdot',fig=fig,ax=ax) 

fig,ax = sl['TESS-p50-5cm'].plotwrt('f',coeff='T',labels=[''],color='k')
fig,ax = sl['TESS-p50-10cm'].plotwrt('f',coeff='T',labels=[''],color='k',linestyle='dashed',fig=fig,ax=ax)
fig,ax = sl['TESS-p50-15cm'].plotwrt('f',coeff='T',labels=[''],color='k',linestyle='dashdot',fig=fig,ax=ax)

