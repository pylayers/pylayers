from pylayers.antprop.antenna import *
A = Antenna('S2R2.sh3')
A.Fsynth()
tau = A.getdelay()
A.elec_delay(tau)
