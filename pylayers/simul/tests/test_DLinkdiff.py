from pylayers.simul.link import *
DL=DLink(L=Layout('defdiff.ini'))
DL.a[0]=DL.a[0]+3
DL.eval(force =True,diffraction=True)
uy = np.where(DL.H.y!=0)

H = Tchannel()