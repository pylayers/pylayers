from pylayers.util.cone import *
from pylayers.gis.layout import *

L=Layout('defstr.ini')

def get_tahe(nseg):
    return L.pt[:,L.tahe[:,L.tgs[nseg]]]


th4=get_tahe(4)
th9=get_tahe(9)
th8=get_tahe(8)
th8[:,1]=th8[:,1]+0.5
th7=get_tahe(7)
th2=get_tahe(2)

cn1 = cone.Cone()
cn1.from2csegs(th4,th9)
cn2 = cone.Cone()
cn2.from2csegs(th7,th8)


Mpta9 = geu.mirror(th2[:,0],th9[:,0],th9[:,1])
Mphe9 = geu.mirror(th2[:,1],th9[:,0],th9[:,1])
typ9,prob9 = cn1.belong_seg(Mpta9,Mphe9)

Mpta8 = geu.mirror(th2[:,0],th8[:,0],th8[:,1])
Mphe8 = geu.mirror(th2[:,1],th8[:,0],th8[:,1])
typ8,prob8 = cn2.belong_seg(Mpta8,Mphe8)




