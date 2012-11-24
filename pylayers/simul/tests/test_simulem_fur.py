from pylayers.simul.simulem import *
from pylayers.signal.bsignal import *
from pylayers.measures.mesuwb import *
import matplotlib.pyplot as plt



M = UWBMesure(13)
cir = TUsignal()
cirf = TUsignal()
cir.readcir("where2cir-tx002-rx012.mat","Tx002")
cirf.readcir("where2-furcir-tx002-rx012.mat","Tx002")

plt.ion()
M.tdd.ch2.plot()
cir.plot(col='blue')
cirf.plot(col='red')
plt.xlim(0,200)
plt.show()
