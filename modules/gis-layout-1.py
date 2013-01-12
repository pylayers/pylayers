import matplotlib.pyplot as plt
from pylayers.gis.layout import *
L = Layout()
L.load('Lstruc.str')
L.loadfur('Furw1.ini')
ax = L.showGs()
plt.show()
