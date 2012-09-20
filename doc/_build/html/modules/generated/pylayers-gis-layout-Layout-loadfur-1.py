import matplotlib.pyplot as plt
from pylayers.gis.layout import *
L = Layout()
L.load('sircut.str')
L.loadfur('Furw1.ini')
f,a = L.showGs()
plt.show()
