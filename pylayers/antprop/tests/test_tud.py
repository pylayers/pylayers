from pylayers.simul.simulem import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
import pylayers.util.pyutil as pyu
# create a Simul object

S = Simul()
# loading a layout 
filestr = 'defstr'
S.layout(filestr+'.str','matDB.ini','slabDB.ini')
try:
    S.L.dumpr()
except:
    S.L.build()
    S.L.dumpw()


S.tx = RadioNode(typ='tx')
S.tx.point([1.2,1,1.4])


# TX / RX
itx=1
irx=1
S.rx = RadioNode(typ='rx')
S.rx.point([8,-1.2,1.5])
S.save()

S.run(itx,irx)


C=Ctilde()
C.load(pyu.getlong(S.dfield[itx][irx],pstruc['DIRTRA']))

#G3=GrRay3D()
#G3.load(S.dtra[itx][irx],S.L)
##fig,ax=S.L.showG('')
##Gr3.show(fig=fig,ax=ax,rayset=[1,2])
#r3=G3.ray3d[1]
#r3.locbas(S.L)

Gt=GrRayTud()
Gt.load(S.dtud[itx][irx],S.dtang[itx][irx],S.drang[itx][irx],S.sl)

rt=Gt.rayTud[1]

#thetas=Gt.get_thetas()
#mat=Gt.get_mat()




#B=rt.inter[0]
#Be=B.eval()

#T=rt.inter[3]
#Te=T.eval()


