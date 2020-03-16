#-*- coding:Utf-8 -*-
from pylayers.simul.link import *
import numpy as np
import matplotlib.pyplot as plt 
import time
print("=======================")
print(" start test_reciprocity.py ")
print("=======================")

DL=DLink(L='defstr.ini')

DL.a=np.array([759,1114,1.0])
DL.b=np.array([767,1114,1.5])
DL.fGHz=np.arange(2,11,0.1)
DL.wav = wvf.Waveform(fcGHz=5,bandGHz=3)
# DL.eval(diffraction=True,force=True,ra_vectorized=True,alg=20152,si_reverb=2)
ra_vectorized = False

DL.eval(diffraction=True,force=True,applywav=True)
DL2=DLink(L='defstr.ini')

DL2.a=np.array([759,1114,1.0])
DL2.b=np.array([767,1114,1.5])
DL2.fGHz=np.arange(2,11,0.1)
DL2.wav = wvf.Waveform(fcGHz=5,bandGHz=3)
# ##########################################################################
# #!!never use the following part except for this particular test file !!!!
# ##########################################################################
DL2.load(DL2.Si,DL.dexist['sig']['grpname'])
if ra_vectorized:
    r = DL2.Si.raysv(DL2.a,DL2.b)
else:
    r = DL2.Si.rays(DL2.a,DL2.b)
rr=r.reciprocal()
DL2.R = rr.to3D(DL2.L)
DL2.R.locbas(DL2.L)
DL2.R.fillinter(DL2.L)
DL2.C=DL2.R.eval(DL2.fGHz)
DL2.a,DL2.b = DL2.b,DL2.a

# ##########################################################################

DL.R.check_reciprocity(DL2.R)
DL.C.check_reciprocity(DL2.C)

