from pylayers.signal.bsignal import *
import numpy as np
alphak = 10.np.random.rand(7)
tauk = 100.np.random.rand(7)
tau = np.arange(0,150,0.1)
y = np.zeros(len(tau))
CIR = TUsignal(alphak,tauk)
CIR.aggcir(alphak,tauk)
CIR.plot(typ='v')
