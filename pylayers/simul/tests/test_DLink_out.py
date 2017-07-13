from pylayers.simul.link import *
import pylayers.signal.waveform as wvf
import pdb
fGHz=np.arange(2.41,10.,0.05)
L=Layout('defstr.lay')
DL=DLink(L=L,fGHz=fGHz)
#DL.Aa=Antenna(typ='Omni')
#DL.Ab=Antenna(typ='Omn')

DL.b=np.array([752,1118,0.5])
DL.eval(force=['sig','ray','Ct','H'],ra_vectorized=True,diffraction=True,ra_number_mirror_cf=2)
dist_a_b = np.sqrt(np.sum((DL.a-DL.b)**2))
#
if DL.R.los:
    ak0 = DL.H.ak[0]
    tk0 = DL.H.tk[0]
    assert tk0*0.3 == dist_a_b, 'invalid distance'
    lak0 = 20* np.log10(ak0)
    Friss= 20*np.log10(2.4)+20*np.log10(dist_a_b) + 32.4
    assert np.allclose(-lak0,Friss,0.1), 'issue in Friss'

