from pylayers.simul.simulem import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
import pylayers.util.pyutil as pyu
import pylayers.signal.bsignal as bs
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
Gt.load(S.dtud[itx][irx],S.dtang[itx][irx],S.drang[itx][irx],S.slab)


# dictionnary of length of interactions
# keys are the number of interactions.
k=Gt.dli.keys()
# Gt.dli is a dictionnary of dictionnary
Gt.dli[k[0]].keys()

print 'evaluation of all rays & interactions'
Gt.eval()

print 'memory size occupied by Interaction matrix = ',Gt.I.I.nbytes/1e6,'MB'
print 'memory size occupied by Ctilde matrix = ',Gt.Ctilde.nbytes/1e6,'MB'


C=Ctilde()
C.load(pyu.getlong(S.dfield[itx][irx],pstruc['DIRTRA']))

freq=Gt.I.f
nfreq=Gt.I.nf
nray=np.array(([Gt.nray]))

Cr=Gt.Ctilde.reshape(nray, nfreq,2,2)
Cr=np.transpose(Gt.Ctilde,(1,0,2,3))

c11 = Cr[:,:,0,0]
c12 = Cr[:,:,0,1]
c21 = Cr[:,:,1,0]
c22 = Cr[:,:,1,1]


C2=Ctilde()
C2.Ctt = bs.FUsignal(freq, c11)
C2.Ctp = bs.FUsignal(freq, c12)
C2.Cpt = bs.FUsignal(freq, c21)
C2.Cpp = bs.FUsignal(freq, c22)
C2.nfreq = Gt.I.nf
C2.nray = Gt.nray
C2.tauk=Gt.delays

plt.ion()

C.show()
C2.show()
























#S2 = Simul()
## loading a layout 
#filestr = 'defstr'
#S2.layout(filestr+'.str','matDB.ini','slabDB.ini')
#try:
#    S2.L.dumpr()
#except:
#    S2.L.build()
#    S2.L.dumpw()


#S2.rx = RadioNode(typ='rx')
#S2.rx.point([1.2,1,1.4])


## TX / RX
#itx=1
#irx=1
#S2.tx = RadioNode(typ='tx')
#S2.tx.point([8,-1.2,1.5])
#S2.save()

#S2.run(itx,irx)

#Gt2=GrRayTud()
#Gt2.load(S2.dtud[itx][irx],S2.dtang[itx][irx],S2.drang[itx][irx],S2.sl)


## dictionnary of length of interactions
## keys are the number of interactions.
#k=Gt2.dli.keys()
## Gt.dli is a dictionnary of dictionnary
#Gt2.dli[k[0]].keys()

#print 'evaluation of all rays & interactions'
#Gt2.eval()


#C2=Ctilde()
#C2.load(pyu.getlong(S2.dfield[itx][irx],pstruc['DIRTRA']))



#CC = Gt.Ctilde

#C=Ctilde()
#C.load(pyu.getlong(S.dfield[itx][irx],pstruc['DIRTRA']))




#Gt.I.eval()


#rt=Gt.rayTud[1]
    
#thetas=Gt.get_thetas()
#mat=Gt.get_mat()




#B=rt.inter[0]
#Be=B.eval()

#T=rt.inter[3]
#Te=T.eval()


