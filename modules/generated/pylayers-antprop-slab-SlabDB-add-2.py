from pylayers.antprop.slab import *
import numpy as np
import matplotlib.pylab as plt
sl = SlabDB('matDB.ini','slabDB.ini')
sl.mat.add(name='CoatingPilkington',cval=1,sigma=2.5e6,typ='epsr')
sl.mat.add(name='GlassPilkington',cval = 6.9,sigma = 5e-4,typ='epsr')
sl.add('Optitherm382',['CoatingPilkington',
# 'GlassPilkington'],[100e-9,0.00382])
fGHz  = np.linspace(0.9,2.2,50)
theta = np.linspace(0,np.pi/2,100)
sl['Optitherm382'].ev(fGHz,theta)
sl['Optitherm382'].pcolor(dB=True)
