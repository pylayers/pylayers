# -*- coding: utf-8 -*-
r"""
================================
UWB Ray tracing simulation  in outdoor scenario
================================

Example taken from 

Joseph W. Schuster, R.J Luebbers
[Techniques for Evaluating the Accuracy of Ray-Tracing Propagation Models for Microcells](http://tinyurl.com/z2plo2v) 
"""
from pylayers.simul.link import *
import matplotlib.cm as cm
L = Layout('Luebbers.ini')
L._visual_check()

# Creating a Link 

DL=DLink(L=L)
DL.fGHz= np.arange(1,6,0.01)
#DL.a = np.array(([25,21.67,2.]))
DL.a = np.array(([37.5,5,2.]))
DL.b = np.array(([12.5,30.,2.]))
#
# Set the antennas
#

DL.Aa=Antenna(typ='Omni')
DL.Ab=Antenna(typ='Omni')

DL.show()

#
# Link evaluation 
#

DL.eval(force=1,diffraction=1)

DL.R.show(L=DL.L)

DL.ir.plot(typ=['l20'],figsize=(8,4))
plt.xlim(100,300)
plt.ylim(-200,-50)
plt.show()
