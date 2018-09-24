from pylayers.gis.layout import *
L = Layout('TA-Office.lay')
p1 = np.array([[0,0,0],[0,0,0]])
p2 = np.array([[10,10,10],[10,10,10]])
seglist = L.seginframe2(p1,p2)
edlist  = [ L.tsg[x] for x in  seglist ]
fig,ax = L.showG('s',edlist=edlist)
