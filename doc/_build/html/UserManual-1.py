from Simul import *
from Layout import *
from numpy import *
import matplotlib.pylab as plt

S = Simul('where2.ini')
S.L.loadfur('FurSiradel.ini')
S.show(furniture=True)
plt.show()