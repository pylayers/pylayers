from pylayers.gis.layout import *
plt.ion()
L = Layout('defstr.ini')
# single point 
#p1 = np.array([3,3])
#p1 =np.append(L.Gt.pos[5],1) 
p1 =L.Gt.pos[5]
# list of points 
#p2 = np.array([18,3]) 
#p2 = np.append(L.Gt.pos[2],2.7)
p2 = L.Gt.pos[2]
#tp2 = np.vstack((p2,p2+np.array([1,1,0]))).T
tp2 = np.vstack((p2,p2+np.array([1,1]))).T
#sgl = L.angleonlink3(p1,tp2)
sgl = L.angleonlink(p1,tp2)
L.showG('s',labels=1)
plt.plot(np.array([p1[0],p2[0]]),np.array([p1[1],p2[1]]))
plt.show()
#L.sla[sgl['s']]
