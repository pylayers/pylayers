from pylayers.antprop.antenna import *
A = Antenna('S2R2.sh3')
A.eval()
tau = A.getdelay()
A.elec_delay(tau)
