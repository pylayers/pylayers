from pylayers.gis.layout import *
L = Layout('example.ini')
L.build()
fig = plt.figure()
fig,ax = L.showGs(fig=fig)
ax = L.showGv(ax=ax)
ti = plt.title('Show Gv')
t = plt.axis('off')
plt.show()
