from pylayers.gis.layout import Layout
import matplotlib.pyplot as plt 
import doctest 

#doctest.testmod(layout)


L = Layout()
L.load('Lstruc.str2')
L.editor()
#ax  = plt.gca()
#ax  = L.showGs(ax=ax)
#plt.show()
# build topological graph 
#L.buildGt()
#L.showG('t')
#plt.title('topological graph')
# build graph of rooms
#L.buildGr()
#L.showG('r')
#plt.title('Graph of rooms')
#plt.show()
#L.buildGv()
#L.showG('v')
#plt.title('Visibility graph')
#plt.savefig('graphGv.png')
