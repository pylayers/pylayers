import matplotlib.pyplot as plt
from pylayers.gis.layout import *
L = Layout('WHERE1.ini')
loadfur('Furw1.ini')
fig = plt.figure()
ax = showGs(fig=fig,furniture=True)
ti = plt.title('loadfur')
plt.show()
