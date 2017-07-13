import mayavi.mlab as mlab
from pylayers.simul.simulem import * #Obsolete
from pylayers.simul.link import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
from pylayers.antprop.signature import *
import pylayers.util.pyutil as pyu
from pylayers.gis.layout import *
from pylayers.util.project import *
import pylayers.signal.bsignal as bs
from datetime import datetime
import time
import pdb
import pickle
import numpy as np
import pylayers.signal.waveform as wvf
#import matplotlib.pyplot as plt

# Set the frequency range
fGHz = np.arange(2,6,0.01)
# Define the Layout
L = Layout('DLR2.lay')
# Create a link for the Layout
S = DLink(L=L,fGHz=fGHz)
ak,tauk=S.eval(force=['Ct','H'])
#S._show3()
S.H.plot()
wav = wvf.Waveform(fcGHz=4,bandGHz=4)
h  = S.H.totime()
## applyu waveform
#hw = S.H.applywav(wav.sfg)
#
#Nray = h.y.shape[0]
#ir = 0
#it = 0
#plt.figure()
#tEk1 = []
#tEk2 = []
#for k in range(Nray):
#    Ek = np.sum(h.y[k,ir,it,:]*h.y[k,ir,it,:])
#    tEk1.append(Ek)
#alpha=np.sqrt(tEk1)/max(np.sqrt(tEk1))
#for k in range(Nray):
#    plt.plot(h.x,h.y[k,ir,it,:]+0.001*k,'k',alpha=alpha[k])
#plt.figure()
#for k in range(Nray):
#    Ek = np.sum(hw.y[k,:]*hw.y[k,:])
#    tEk2.append(Ek)
#    plt.plot(hw.x,hw.y[k,:],'k',alpha=alpha[k])
