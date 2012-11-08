from pylayers.util import project
from pylayers.signal.bsignal import *
import matplotlib.pylab as plt

cir1 = TUsignal()
cir2 = TUsignal()
cir3 = TUsignal()
cir1.readcir('defaultcir-1-2-p001.mat','1')
cir2.readcir('defaultcir-1-2-p002.mat','1')
cir3.readcir('defaultcir-1-2-p003.mat','1')
plt.ion()
fig = plt.figure()
cir1.show(fig)
fig = plt.figure()
cir2.show(fig)
fig = plt.figure()
cir3.show(fig)
