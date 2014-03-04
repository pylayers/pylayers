from pylayers.util.project import *
from pylayers.measures.mesuwb import *
import matplotlib.pylab as plt
M  = UWBMesure(1)
T  = M.tdd
freq,pl = T.PL(3,7,10)
