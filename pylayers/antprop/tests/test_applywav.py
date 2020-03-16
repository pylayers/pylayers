import numpy as np
from pylayers.simul.link import *
import matplotlib.pyplot as plt
import time
print("=======================")
print(" start test_applywav.py ")
print("=======================")

fGHz=np.arange(2,11,0.1)
DL=DLink(L='defstr.ini',fGHz=fGHz)
DL.a=np.array([759,1114,1.0])
DL.b=np.array([767,1114,1.5])
DL.wav = wvf.Waveform(fcGHz=5,bandGHz=3)
DL.eval(diffraction=True,force=True,cutoff=6,ra_vectorized=True,applywav=False)
cir = DL.H.get_cir(DL.wav.sfg)



cir.plot(typ='v')
