from pylayers.simul.link import *

DL=DLink(L=Layout('defstr.ini'))
DL.fGHz=np.arange(2.41,2.43,0.01)
DL.b=DL.b+1
DL.Aa=Antenna(typ='Omni')
DL.Ab=Antenna(typ='Omni')
DL.eval() 
dist_a_b = np.sqrt(np.sum((DL.a-DL.b)**2))

if DL.R.los:
    ak0 = DL.H.ak[0]
    tk0 = DL.H.tk[0]
    assert tk0*0.3 == dist_a_b, 'invalid distance'
    lak0 = 20* np.log10(ak0)
    Friss= 20*np.log10(2.4)+20*np.log10(dist_a_b) + 32.4
    assert np.allclose(-lak0,Friss,0.1), 'issue in Friss'