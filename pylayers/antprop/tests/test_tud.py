from pylayers.simul.simulem import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
import pylayers.util.pyutil as pyu
#import pylayers.antprop. as pyu
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
Gt.load(S.dtud[itx][irx],S.dtang[itx][irx],S.drang[itx][irx],S.L.sl)

rt=Gt.rayTud[1]

#C = Ctilde()
#C.load(pyu.getlong(S.dfield[itx][irx],'output'))

#rt.eval(freq)

ir = 1
freq = C.Ctt.x
tang = C.tang[ir]
rang = C.rang[ir]
tauk = C.tauk[ir]
Ctt  = C.Ctt.y[ir,:]
Ctp  = C.Ctp.y[ir,:]
Cpt  = C.Cpt.y[ir,:]
Cpp  = C.Cpp.y[ir,:]
plt.ion()
ax1 = plt.subplot(221)
plt.plot(freq,abs(Ctt))
plt.title('Ctt')
plt.subplot(222,sharex=ax1,sharey=ax1)
plt.plot(freq,abs(Ctp))
plt.title('Ctp')
plt.subplot(223,sharex=ax1,sharey=ax1)
plt.plot(freq,abs(Cpt))
plt.title('Cpt')
plt.subplot(224,sharex=ax1,sharey=ax1)
plt.plot(freq,abs(Cpp))
plt.title('Cpp')
## print launching parameters
#S.palch.info()

#S2 = Simul()
## loading a layout 
#filestr = 'defstr'
#S2.layout(filestr+'.str','matDB.ini','slabDB.ini')
#S2.tx = RadioNode(typ='tx')
#S2.tx.point([8,-1.2,1.5])
## TX / RX
#itx=1
#irx=1
#S2.rx = RadioNode(typ='rx')
#S2.rx.point([1.2,1,1.4])
#S2.save()
#S2.run(itx,irx)
#
#
#C=Ctilde()
#C.load(pyu.getlong(S2.dfield[itx][irx],pstruc['DIRTRA']))
#
#G3=GrRay3D()
#G3.load(S2.dtra[itx][irx],S2.L)
##fig,ax=S2.L.showG('')
##Gr3.show(fig=fig,ax=ax,rayset=[1,2])
#
#
#Gt=GrRayTud()
#Gt.load(S2.dtud[itx][irx],S2.dtang[itx][irx],S2.drang[itx][irx],S2.L.sl)
#
#r3=G3.ray3d[1]
#r3.locbas(S2.L)
##r3.show3()
#rt=Gt.rayTud[1]
#
#C = Ctilde()
#C.load(pyu.getlong(S2.dfield[itx][irx],'output'))
#
##rt.eval(freq)
#
#ir = 1
#freq = C.Ctt.x
#tang = C.tang[ir]
#rang = C.rang[ir]
#tauk = C.tauk[ir]
#Ctt  = C.Ctt.y[ir,:]
#Ctp  = C.Ctp.y[ir,:]
#Cpt  = C.Cpt.y[ir,:]
#Cpp  = C.Cpp.y[ir,:]
#plt.ion()
#plt.figure()
#ax1 = plt.subplot(221)
#plt.plot(freq,abs(Ctt))
#plt.title('Ctt')
#plt.subplot(222,sharex=ax1,sharey=ax1)
#plt.plot(freq,abs(Ctp))
#plt.title('Ctp')
#plt.subplot(223,sharex=ax1,sharey=ax1)
#plt.plot(freq,abs(Cpt))
#plt.title('Cpt')
#plt.subplot(224,sharex=ax1,sharey=ax1)
#plt.plot(freq,abs(Cpp))
#plt.title('Cpp')
### print launching parameters
##S2.palch.info()
#
### ang Tx : angular step from Tx
##S2.palch.angTx  = 1
#
### ISB ang Incident Shadow Boundary angle (degree) 
##S2.palch.ISBang = 90  
#
### ray elimination Threshold 
##S2.palch.ethreshold = 0.001
#
### maximum depth
##S2.palch.maxdeep  = 10
#
### typealgo = 0 (include diffraction) 1 (no diffraction)
##S2.palch.typalgo = 1
##title = str(S2.palch.angTx) + '-' +\
##        str(S2.palch.ISBang) + '-' +\
##        str(S2.palch.ethreshold) + '-' + \
##        str(S2.palch.maxdeep) + '-' + \
##        str(S2.palch.typalgo)
#
##S2.palch.save()
##S2.pafreq.fghzmin=2
##S2.pafreq.fghzmax=11
##S2.pafreq.nf=181
##S2.pafreq.save()
### showing the simulation 
##print "Launching "
##print "-----------------"
##S2.launching(1)
#
### retrieve the launching tree
#
##L1 = S2.getlaunch(1)
#
### display the launching tree for different depths
#
##fig = plt.figure(figsize=(10,10))
##plt.title('launching parameters '+title+' '+filestr )
##plt.axis('off')
##N = S2.palch.maxdeep
##M = N/2
###
##for k in range(N):
##    ax = fig.add_subplot(M,2,k+1)
##    fig,ax = L1.show(S2.L,k+1,f=fig)
#
##fig.savefig(pylayersdir+'/doc/auto_examples/simul/'+filestr+'-launching.png')    
##print "Tracing "
##print "-----------------"
##print "purc :",S2.config.get('tud','purc')
##fig = plt.figure()
##S2.tracing(1,1)
##gr = GrRay3D()
##gr.load(S2.dtra[1][1],S2.L)
##f,a = S2.L.showGs(fig=fig)
###plt.axis('on')
##gr.show(fig=f,ax=a,rayset=np.arange(100))
##print "Tratotud "
##print "-----------------"
##print "purc :",S2.config.get('tud','purc')
##S2.tratotud(1,1)
##gt = GrRayTud()
### loading rays in tud format 
##gt.load(S2.dtud[1][1],S2.dtang[1][1],S2.drang[1][1],S2.L.sl)
##print "Evalfield "
##print "-----------------"
##S2.field(1,1)
##S2.cir(1,1)
##f = plt.figure()
##S2.pltcir(1,1,fig=f)
### ang Tx : angular step from Tx
##S2.palch.angTx  = 1
#
### ISB ang Incident Shadow Boundary angle (degree) 
##S2.palch.ISBang = 90  
#
### ray elimination Threshold 
##S2.palch.ethreshold = 0.001
#
### maximum depth
##S2.palch.maxdeep  = 10
#
### typealgo = 0 (include diffraction) 1 (no diffraction)
##S2.palch.typalgo = 1
##title = str(S2.palch.angTx) + '-' +\
##        str(S2.palch.ISBang) + '-' +\
##        str(S2.palch.ethreshold) + '-' + \
##        str(S2.palch.maxdeep) + '-' + \
##        str(S2.palch.typalgo)
#
##S2.palch.save()
##S2.pafreq.fghzmin=2
##S2.pafreq.fghzmax=11
##S2.pafreq.nf=181
##S2.pafreq.save()
### showing the simulation 
##print "Launching "
##print "-----------------"
##S2.launching(1)
#
### retrieve the launching tree
#
##L1 = S2.getlaunch(1)
#
### display the launching tree for different depths
#
##fig = plt.figure(figsize=(10,10))
##plt.title('launching parameters '+title+' '+filestr )
##plt.axis('off')
##N = S2.palch.maxdeep
##M = N/2
###
##for k in range(N):
##    ax = fig.add_subplot(M,2,k+1)
##    fig,ax = L1.show(S2.L,k+1,f=fig)
#
##fig.savefig(pylayersdir+'/doc/auto_examples/simul/'+filestr+'-launching.png')    
##print "Tracing "
##print "-----------------"
##print "purc :",S2.config.get('tud','purc')
##fig = plt.figure()
##S2.tracing(1,1)
##gr = GrRay3D()
##gr.load(S2.dtra[1][1],S2.L)
##f,a = S2.L.showGs(fig=fig)
###plt.axis('on')
##gr.show(fig=f,ax=a,rayset=np.arange(100))
##print "Tratotud "
##print "-----------------"
##print "purc :",S2.config.get('tud','purc')
##S2.tratotud(1,1)
##gt = GrRayTud()
### loading rays in tud format 
##gt.load(S2.dtud[1][1],S2.dtang[1][1],S2.drang[1][1],S2.L.sl)
##print "Evalfield "
##print "-----------------"
##S2.field(1,1)
##S2.cir(1,1)
##f = plt.figure()
##S2.pltcir(1,1,fig=f)
#thetas=Gt.get_thetas()
#mat=Gt.get_mat()




#B=rt.inter[0]
#Be=B.eval()

#T=rt.inter[3]
#Te=T.eval()


