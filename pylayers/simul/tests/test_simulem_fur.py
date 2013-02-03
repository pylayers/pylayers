from pylayers.simul.simulem import *
from pylayers.signal.bsignal import *
from pylayers.measures.mesuwb import *
import matplotlib.pyplot as plt
from pylayers.gis.layout import * 



#M=UWBMesure(173)
M=UWBMesure(13)
#M=UWBMesure(1)

cir=TUsignal()
cirf=TUsignal()
#cir.readcir("where2cir-tx001-rx145.mat","Tx001")
#cirf.readcir("where2-furcir-tx001-rx145.mat","Tx001")

cir.readcir("where2cir-tx002-rx012.mat","Tx002")
#cirf.readcir("where2-furcir-tx002-rx012.mat","Tx002")

#cir.readcir("where2cir-tx001-rx001.mat","Tx001")
#cirf.readcir("where2-furcir-tx001-rx001.mat","Tx001")

plt.ion()
fig = plt.figure()
fig.subplots_adjust(hspace=0.5)

ax1 = fig.add_subplot(411,title="points and layout")
L=Layout()
L.load('siradel-cut-fur.str')
L.showGs(fig=fig,ax=ax1)
ax1.plot(M.tx[0],M.tx[1],'or')
#ax1.plot(M.rx[1][0],M.rx[1][1],'ob')
ax1.plot(M.rx[2][0],M.rx[2][1],'ob')


ax2 = fig.add_subplot(412,title="Measurement")
M.tdd.ch2.plot()

#ax3 = fig.add_subplot(413,title="Simulation with furniture",sharex=ax2,sharey=ax2)
#cirf.plot(col='red')

ax4 = fig.add_subplot(414,title="Simulation",sharex=ax2,sharey=ax2)
cir.plot(col='blue')

plt.show()

