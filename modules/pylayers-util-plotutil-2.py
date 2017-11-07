import scipy as sp
import matplotlib.pyplot as plt
from pylayers.util.plotutil import *
N   = 10
pt  = sp.rand(2,N)
ph  = sp.rand(2,N)
f,a = displot(pt,ph)
txt = plt.title('pylayers.util.geomutil.displot(pt,ph) : plot 10 random segments')
plt.show()
