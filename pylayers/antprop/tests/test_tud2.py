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

Gt=GrRayTud()
Gt.load(S.dtud[itx][irx],S.dtang[itx][irx],S.drang[itx][irx],S.sl)

# dictionnary of length of interactions
# keys are the number of interactions.
k=Gt.dli.keys()
# Gt.dli is a dictionnary of dictionnary
Gt.dli[k[0]].keys()

print 'evaluation of all rays & interactions'
Gt.eval()

print 'memory size occupied by Interaction matrix = ',Gt.I.I.nbytes/1e6,'MB'
print 'memory size occupied by Ctilde matrix = ',Gt.Ctilde.nbytes/1e6,'MB'

CC = Gt.Ctilde

C=Ctilde()
C.load(pyu.getlong(S.dfield[itx][irx],pstruc['DIRTRA']))




#Gt.I.eval()


#rt=Gt.rayTud[1]
    
#thetas=Gt.get_thetas()
#mat=Gt.get_mat()




#B=rt.inter[0]
#Be=B.eval()

#T=rt.inter[3]
#Te=T.eval()


