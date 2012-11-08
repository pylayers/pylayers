from pylayers.util import project
from pylayers.signal.bsignal import *
import matplotlib.pylab as plt

cir1 = TUsignal()
cir2 = TUsignal()
cir3 = TUsignal()
cir1.readcir('defaultcir-tx001-rx001.mat','Tx001')
cir2.readcir('defaultcir-tx001-rx002.mat','Tx001')
plt.ion()
fig = plt.figure()
cir1.show(fig)
fig = plt.figure()
cir2.show(fig)
