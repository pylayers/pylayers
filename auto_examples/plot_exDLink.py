# -*- coding: utf-8 -*-
r"""
=================================================
Evaluation of a radio link DLink 
=================================================

The example load a Layout, set the transmitter iand receiver 
position and evaluates the link
"""
from pylayers.simul.link import *
import pylayers.signal.waveform as wvf
import pdb

# set the frequency range
fGHz=np.arange(2.41,10.,0.05)
# set the layout
L=Layout('defstr.ini')
# set the link
DL=DLink(L=L,fGHz=fGHz)
DL.b=np.array([755,1110,1.2])
#DL.Aa=Antenna(typ='Omni')
#DL.Ab=Antenna(typ='Omn')
# Point outside
#DL.eval(force=['sig','ray','Ct','H'],ra_vectorized=True,diffraction=True,ra_ceil_H=0)
ak,tauk = DL.eval(force=['sig','ray','Ct','H'],ra_vectorized=True,diffraction=True)
DL.H.show()
plt.show()
