from pylayers.measures.mesuwb import *
import matplotlib.pylab as plt
ntx = 2
M  = UWBMeasure(ntx)
T  = M.tdd
fig = plt.figure()
t = plt.title('test Tdd.show  Tx='+str(ntx))
T.show()
