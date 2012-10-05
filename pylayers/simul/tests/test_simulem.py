from pylayers.simul.simulem import *

# create a Simul object

S = Simul()

# loading a layout 
filestr = 'defstr'
S.layout(filestr+'.str','matDB.ini','slabDB.ini')
# setting transmitter

S.tx = RadioNode(typ='tx')
S.tx.point([1.2,1,1.4])

# setting receiver

S.rx = RadioNode(typ='rx')
S.rx.point([8,-1.2,1.5])

S.save()

# print launching parameters
S.palch.info()

# ang Tx : angular step from Tx
S.palch.angTx  = 1

# ISB ang Incident Shadow Boundary angle (degree) 
S.palch.ISBang = 90  

# ray elimination Threshold 
S.palch.ethreshold = 0.001

# maximum depth
S.palch.maxdeep  = 10

# typealgo = 0 (include diffraction) 1 (no diffraction)
S.palch.typalgo = 1
title = str(S.palch.angTx) + '-' +\
        str(S.palch.ISBang) + '-' +\
        str(S.palch.ethreshold) + '-' + \
        str(S.palch.maxdeep) + '-' + \
        str(S.palch.typalgo)

S.palch.save()
S.pafreq.fghzmin=2
S.pafreq.fghzmax=11
S.pafreq.nf=181
S.pafreq.save()
# showing the simulation 
print "Launching "
print "-----------------"
S.launching(1)

# retrieve the launching tree

L1 = S.getlaunch(1)

# display the launching tree for different depths

fig = plt.figure(figsize=(10,10))
plt.title('launching parameters '+title+' '+filestr )
plt.axis('off')
N = S.palch.maxdeep
M = N/2
#
for k in range(N):
    ax = fig.add_subplot(M,2,k+1)
    fig,ax = L1.show(S.L,k+1,f=fig)

fig.savefig(pylayersdir+'/doc/auto_examples/simul/'+filestr+'-launching.png')    
print "Tracing "
print "-----------------"
print "purc :",S.config.get('tud','purc')
fig = plt.figure()
S.tracing(1,1)
gr = GrRay3D()
gr.load(S.dtra[1][1],S.L)
f,a = S.L.showGs(fig=fig)
#plt.axis('on')
gr.show(fig=f,ax=a,rayset=np.arange(100))
print "Tratotud "
print "-----------------"
print "purc :",S.config.get('tud','purc')
S.tratotud(1,1)
gt = GrRayTud()
# loading rays in tud format 
gt.load(S.dtud[1][1],S.dtang[1][1],S.drang[1][1],S.L.sl)
print "Evalfield "
print "-----------------"
S.field(1,1)
S.cir(1,1)
f=plt.figure()
S.pltcir(1,1,fig=f)
