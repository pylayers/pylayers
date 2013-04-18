from pylayers.simul.simulem import *
from pylayers.antprop.raysc import *
from pylayers.antprop.channel import *
import pylayers.util.pyutil as pyu
import pylayers.signal.bsignal as bs
import time
# create a Simul object

def plotray(r):

    plt.ion()
    plt.close('all')
    fig=plt.figure('Cpp')
    f,a=C.Cpp.plot(fig=fig,iy=np.array(([r])))
    f,a,Cn.Cpp.plot(fig=fig,iy=np.array(([r])))
    a[0].legend(('Fried','new'))

    fig2=plt.figure('Ctt')
    f,a=C.Ctt.plot(fig=fig2,iy=np.array(([r])))
    f,a,Cn.Ctt.plot(fig=fig2,iy=np.array(([r])))
    a[0].legend(('Fried','new'))

    Gt.info(r)

    plt.show()


###################################
#   Simulation creation 
#
#################################
S = Simul()
# loading a layout 
filestr = 'TA-Office'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
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


###################################
#   New load function 
#   load pulray tud file
#
#################################


Gt=GrRayTud()
Gt.load(S.dtud[itx][irx],S.dtang[itx][irx],S.drang[itx][irx],S.slab)



a=time.time()
print 'evaluation of all rays & interactions'
Gt.eval()
b=time.time()

#print 'memory size occupied by Interaction matrix = ',Gt.I.I.nbytes/1e6,'MB'
#print 'memory size occupied by Ctilde matrix = ',Gt.Ctilde.nbytes/1e6,'MB'
print 'evaluation in ',(b-a) ,'seconds'

C=Ctilde()
C.load(pyu.getlong(S.dfield[itx][irx],pstruc['DIRTRA']))

freq=Gt.I.f
nfreq=Gt.I.nf
nray=np.array(([Gt.nray]))


Cr=np.swapaxes(Gt.Ctilde,1,0)
#Cr=Gt.Ctilde.reshape(nray, nfreq,2,2)
#Cr=np.transpose(Gt.Ctilde,(1,0,2,3))



c11 = Cr[:,:,0,0]
c12 = Cr[:,:,0,1]
c21 = Cr[:,:,1,0]
c22 = Cr[:,:,1,1]

Cn=Ctilde()
Cn.Cpp = bs.FUsignal(freq, c11)
Cn.Ctp = bs.FUsignal(freq, c12)
Cn.Cpt = bs.FUsignal(freq, c21)
Cn.Ctt = bs.FUsignal(freq, c22)
Cn.nfreq = Gt.I.nf
Cn.nray = Gt.nray
Cn.tauk=Gt.delays


r=0
plotray(r)



















### RECIPROCITY TEST

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
#Gt2.load(S2.dtud[itx][irx],S2.dtang[itx][irx],S2.drang[itx][irx],S.slab)



#a=time.time()
#print 'evaluation of all rays & interactions'
#Gt2.eval()
#b=time.time()

#print 'memory size occupied by Interaction matrix = ',Gt.I.I.nbytes/1e6,'MB'
#print 'memory size occupied by Ctilde matrix = ',Gt.Ctilde.nbytes/1e6,'MB'
#print 'evaluation in ',(b-a) ,'seconds'

#C2=Ctilde()
#C2.load(pyu.getlong(S.dfield[itx][irx],pstruc['DIRTRA']))

#freq=Gt2.I.f
#nfreq=Gt2.I.nf
#nray=np.array(([Gt2.nray]))

#Cr2=Gt2.Ctilde.reshape(nray, nfreq,2,2)
#Cr2=np.transpose(Gt2.Ctilde,(1,0,2,3))

#c11_ = Cr2[:,:,0,0]
#c12_ = Cr2[:,:,0,1]
#c21_ = Cr2[:,:,1,0]
#c22_ = Cr2[:,:,1,1]


#Cn2=Ctilde()
#Cn2.Ctt = bs.FUsignal(freq, c11_)
#Cn2.Ctp = bs.FUsignal(freq, c12_)
#Cn2.Cpt = bs.FUsignal(freq, c21_)
#Cn2.Cpp = bs.FUsignal(freq, c22_)
#Cn2.nfreq = Gt2.I.nf
#Cn2.nray = Gt2.nray
#Cn2.tauk=Gt2.delays

##plt.ion()



#plt.figure('2Cpp - Fried')
#C.Cpp.plot()
#plt.figure('2Cpp - new')
#Cn2.Cpp.plot()

#plt.figure('2Ctt - Fried')
#C.Ctt.plot()
#plt.figure('2Ctt - new')
#Cn2.Ctt.plot()






############################################################################












#plt.figure('Ctp - Fried')
#C.Ctp.plot()
#plt.figure('Ctp - new')
#C2.Ctp.plot()

#plt.figure('Cpt - Fried')
#C.Cpt.plot()
#plt.figure('Cpt - new')
#C2.Cpt.plot()


#raw_input('press any key to close figure')
#plt.close('all')














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


