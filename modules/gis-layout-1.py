import matplotlib.pyplot as plt
from pylayers.gis.layout import *
L = Layout('WHERE1.ini')
L.loadfur('Furw1.ini')
fig = plt.figure()
ax = L.showGs(fig=fig,furniture=True)
ti = plt.title('loadfur')
plt.show()
