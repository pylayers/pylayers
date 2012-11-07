from pylayers.util import project
from pylayers.signal.bsignal import *
import matplotlib.pylab as plt

cir1 = TUsignal()
cir2 = TUsignal()
cir3 = TUsignal()
cir1.readcir('defaultcir-ap6-ag1-p001.mat','6')
cir2.readcir('defaultcir-ap6-ag1-p002.mat','6')
cir3.readcir('defaultcir-ap6-ag1-p003.mat','6')

fig = plt.figure()
cir1.show(fig)
cir2.show(fig)
cir3.show(fig)
