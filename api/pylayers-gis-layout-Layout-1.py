import matplotlib.pyplot as plt
from pylayers.gis.layout import *
L = Layout('WHERE1.ini')
L.loadfur('Furw1.ini')
fig = plt.figure()
ax = fig.gca()
fig,ax = L.showGs(fig=fig,ax=ax,furniture=True)
ti = plt.title('loadfur')
plt.show()
