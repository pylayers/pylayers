from pylayers.gis.layout import *
import pylayers.util.geomutil as geu
L=Layout('WHERE1.ini')
L.build()

lc = [(n,L.Gt.node[n]['polyg'].isconvex()) for n in L.Gt.nodes()]
cnc = [n for n in L.Gt.nodes() if not L.Gt.node[n]['polyg'].isconvex()]
fig,ax=L.showG('st',labels=True)
for cy,c in lc:
    if c:
        print cy,  
        fig,ax= L.Gt.node[cy]['polyg'].plot(color='blue',alpha=0.5,fig=fig,ax=ax)
    else:
        fig, ax = L.Gt.node[cy]['polyg'].plot(color='red', alpha=0.5,fig=fig,ax=ax)