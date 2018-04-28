from pylayers.gis.layout import *
L = Layout('TA-Office.ini')
p1 = np.array([[0,0,0],[0,0,0]])
p2 = np.array([[10,10,10],[10,10,10]])
seglist = L.seginframe2(p1,p2)
edlist  = map(lambda x: L.tsg[x],seglist)
fig,ax = L.showG('s',edlist=edlist)
