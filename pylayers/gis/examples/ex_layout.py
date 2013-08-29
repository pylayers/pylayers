from pylayers.gis.layout import Layout
import matplotlib.pyplot as plt 
import doctest 

#doctest.testmod(layout)


#L = Layout('TA-Office.ini')
L = Layout('WHERE1.ini')
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()
#L.editor()
fig = plt.gcf()
#ax1  = fig.add_subplot(221)
ax1  = fig.add_subplot(321)
L.display['thin']=True
fig,ax1  = L.showG(graph='s',fig=fig,ax=ax1)
#L.display['edlabel']=True
#L.display['edlblsize']=50
# display selected segments
L.display['thin']=True
L.showG(fig=fig,ax=ax1,graph='t')
fig = plt.gcf()
ax1 = plt.gca()
fig,ax1 =  L.showGs(fig=fig,ax=ax1,edlist=[125],width=4)
ax11 = fig.add_subplot(322)
L.showG(fig=fig,ax=ax11,graph='s')
#plt.savefig('graphGs.png')
#build topological graph 
ax2 = fig.add_subplot(323)
L.showG(fig=fig,ax=ax2,graph='t')
plt.title('Topological graph')
#plt.savefig('graphGt.png')
# build graph of rooms
ax3 = fig.add_subplot(324)
L.showG(fig=fig,ax=ax3,graph='r')
#plt.savefig('graphGr.png')
plt.title('Graph of rooms')
ax4 = fig.add_subplot(325)
L.showG(fig=fig,ax=ax4,graph='v')
plt.title('Visibility graph')
ax5 = fig.add_subplot(326)
L.showG(fig=fig,ax=ax5,graph='i')
plt.title('Interaction graph')
plt.show()
#plt.savefig('graphs.png')
