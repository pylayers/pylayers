from pylayers.gis.layout import Layout
import matplotlib.pyplot as plt 
import doctest 

#doctest.testmod(layout)


L = Layout()
L.load('TA-Office.str')
#L.editor()
fig = plt.gcf()
ax1  = fig.add_subplot(221)
fig,ax1  = L.showGs(fig=fig,ax=ax1)
#plt.savefig('graphGs.png')
# build topological graph 
L.buildGt()
ax2 = fig.add_subplot(222)
L.showG(fig=fig,ax=ax2,graph='t')
plt.title('Topological graph')
#plt.savefig('graphGt.png')
# build graph of rooms
L.buildGr()
ax3 = fig.add_subplot(223)
L.showG(fig=fig,ax=ax3,graph='r')
#plt.savefig('graphGr.png')
plt.title('Graph of rooms')
L.buildGv()
ax4 = fig.add_subplot(224)
L.showG(fig=fig,ax=ax4,graph='v')
#L.showG('v',figsize=(20,20))
plt.title('Visibility graph')
plt.savefig('graphs.png')
