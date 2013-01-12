from pylayers.simul.simulem import *
from pylayers.gis.layout import *
from numpy import *
import matplotlib.pylab as plt

S = Simul('example.ini')
S.L.showGs()
plt.show()