from pylayers.gis.layout import *
L = Layout('exemple.str')
build()
fig = plt.figure()
fig,ax = showGs(fig=fig)
ax = showGv(ax=ax)
ti = plt.title('Show Gv')
t = plt.axis('off')
plt.show()
