from pylayers.simul.simulem import *
from pylayers.signal.bsignal import *
from pylayers.measures.mesuwb import *
import matplotlib.pyplot as plt
from pylayers.gis.layout import * 



M=UWBMesure(173)

cir=TUsignal()
cirf=TUsignal()
cir.readcir("where2cir-tx001-rx145.mat","Tx001")
cirf.readcir("where2-furcir-tx001-rx145.mat","Tx001")

#cir.readcir("where2cir-tx002-rx012.mat","Tx002")
#cirf.readcir("where2-furcir-tx002-rx012.mat","Tx002")

plt.ion()
fig = plt.figure()
ax1 = fig.add_subplot(411)
cir.plot(col='blue')
ax2 = fig.add_subplot(412,sharex=ax1,sharey=ax1)
cirf.plot(col='red')
ax3 = fig.add_subplot(413,sharex=ax1,sharey=ax1)
M.tdd.ch2.plot()
ax4 = fig.add_subplot(414)
L=Layout()
L.load('siradel-cut-fur.str')
L.showGs(fig=fig,ax=ax4)
ax4.plot(M.tx[0],M.tx[1],'or')
ax4.plot(M.rx[1][0],M.rx[1][1],'ob')
plt.show()

