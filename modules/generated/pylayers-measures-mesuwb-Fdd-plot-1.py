from pylayers.util.project import *
from pylayers.measures.mesuwb import *
import matplotlib.pylab as plt
M  = UWBMeasure(1)
F  = M.fdd
fig = plt.figure()
F.plot('moddB')
fig = plt.figure()
F.plot('mod')
F.plot('ang')
