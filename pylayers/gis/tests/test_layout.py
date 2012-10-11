from pylayers.gis.layout import Layout
import matplotlib.pyplot as plt 
import doctest 

#doctest.testmod(layout)


L = Layout()
L.load('TA-Office.str')
#L.editor()
fig = plt.gcf()
ax  = plt.gca()
fig,ax  = L.showGs(fig=fig,ax=ax)
plt.savefig('graphGs.png')
# build topological graph 
L.buildGt()
L.showG('t')
plt.title('topological graph')
plt.savefig('graphGt.png')
# build graph of rooms
L.buildGr()
L.showG('r')
plt.savefig('graphGr.png')
#plt.title('Graph of rooms')
L.buildGv()
L.showG('v',figsize=(20,20))
plt.title('Visibility graph')
plt.savefig('graphGv.png')
