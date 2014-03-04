from pylayers.util.project import *
from pylayers.measures.mesuwb import *
import matplotlib.pylab as plt
M  = UWBMesure(1)
F  = M.fdd
fig = plt.figure()
F.plot('moddB')
plt.tight_layout()
fig = plt.figure()
F.plot('mod')
plt.tight_layout()
fig = plt.figure()
F.plot('ang')
plt.tight_layout()
plt.show()
