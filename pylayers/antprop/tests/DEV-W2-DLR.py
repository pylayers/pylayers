# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Simulation on WHERE2 DLR environment

# <codecell>

from pylayers.simul.simulem import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
from pylayers.antprop.signature import *
import pylayers.util.pyutil as pyu
from pylayers.gis.layout import *
import pylayers.signal.bsignal as bs
from datetime import datetime
import time

# <markdowncell>

# Create the simulation object with `DLR.ini` layout  

# <codecell>

S = Simul()
# loading a layout 
filestr = 'DLR'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
try:
    S.L.dumpr()
except:
    S.L.build()
    S.L.dumpw()
S.L.build()
S.L.dumpw()
S.L.display['ednodes']=False
S.L.display['nodes']=False
S.L.display['title']='DLR WP4 WHERE2 measurement site'
S.L.display['overlay']=False
fig,ax = S.L.showGs()    

# <markdowncell>

# Select Tx and Rx positions

# <headingcell level=2>

# 
# Measures from WHERE2 WP4 measurement campaign @ DLR

# <codecell>

MT_DLR_RDTMaster={ 2:[2.12,0,1.275], 3:[3.12,0,1.275], 4:[4.12,0,1.275], 5:[5.12,0,1.275], 6:[6.12,0,1.275], 7:[7.12,0,1.275], 8:[8.12,0,1.275], 9:[9.12,0,1.275], 10:[10.12,0,1.275], 11:[11.12,0,1.275], 12:[12.12,0,1.275], 13:[13.12,0,1.275], 14:[14.12,0,1.275], 15:[15.12,0,1.275], 16:[16.12,0,1.275], 17:[17.12,0,1.275], 18:[18.12,0,1.275], 19:[19.12,0,1.275], 20:[20.12,0,1.275], 21:[21.12,0,1.275], 22:[22.12,0,1.275], 23:[23.12,0,1.275], 24:[24.12,0,1.275], 25:[25.12,0,1.275], 26:[26.12,0,1.275], 27:[27.12,0,1.275], 28:[28.12,0,1.275], 29:[29.12,0,1.275], 30:[30.12,0,1.275], 31:[31.12,0,1.275], 33:[30.12,0,1.275], 62:[24.12,0,1.275], 63:[24.62,0,1.275], 64:[25.12,0,1.275], 65:[25.62,0,1.275], 66:[26.12,0,1.275], 67:[26.62,0,1.275], 68:[27.12,0,1.275], 69:[27.62,0,1.275], 70:[28.12,0,1.275], 71:[28.62,0,1.275], 72:[29.12,0,1.275], 73:[29.62,0,1.275], 34:[30.12,0,1.275], 35:[30,1.38,1.275], 36:[30,1.88,1.275], 37:[30,2.38,1.275], 38:[30,2.88,1.275], 39:[30,3.88,1.275], 40:[30,4.88,1.275], 41:[30,5.88,1.275]}

# <codecell>

TrolleyMT_ACO_04={ 2:60, 3:59, 4:58, 5:58, 6:58, 7:58, 8:58, 9:58, 10:58, 11:58, 12:58,
13:57, 14:52, 15:12, 16:51, 17:11, 18:50, 19:10, 20:53, 21:54, 22:55, 23:56, 24:56, 25:56,
26:56, 27:56, 28:56, 29:56, 30:56, 31:56, 33:56, 62:56, 63:56, 64:56, 65:56, 66:56, 67:56,
68:56, 69:56, 70:56, 71:56, 72:56, 73:56, 34:56, 35:56, 36:56, 37:56, 38:56, 39:56, 40:56, 41:56} 

# <codecell>

MT_ACO_04 = { 
    60: [12.5, 2, 1.28],
    59 :[12.5, 1.5, 1.28],
    58 :[12.5, 1, 1.28 ],
    57 :[12.5, 0.5, 1.28],
    52 :[12.5, 0, 1.28],
    12 :[12,  0,  1.28],
    51 :[ 11.5, 0, 1.28],
    11 :[ 11, 0, 1.28],
    50 : [10.5, 0, 1.28],
    10 :[10, 0, 1.28],
    53 : [10, -0.5, 1.28],
    54 : [10, -1, 1.28],
    55 : [10, -1.5, 1.28],
    56 : [10, -2, 1.28],
}

# <codecell>

Dongle = 389
AnchorNodes = {390:{'name':'MT_ACO_05','coord':[6,0.81,1.64]},
               386:{'name':'MT_ACO_08','coord':[30.574,2.8,1.291]},
               391:{'name':'MT_ACO_07','coord':[11.78,-5.553,1.5]},
               385:{'name': 'MT_ACO_01','coord':[19.52,-0.69,1.446]},
               387:{'name':'MT_ACO_03','coord':[28.606,-0.74,1.467]},
               400:{'name':'MT_ACO_02','coord':[30.574,2.8,1.291]},
               1:{'name':'MT_DLR_RTDSlave','coord':[0.85,0,1.18]}
              }

# <codecell>

#
#  Define the link to selected
#
A = 390
B = 386

# <headingcell level=2>

# Loading the Layout

# <codecell>

L = Layout('DLR.ini')
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()
#L.showGs()
fig,ax = L.showG('r')
plt.axis('off')

# <codecell>

path = '/home/uguen/Bureau/WHERE2/WHERE2-DLR-Simulations/'
ZigbeeNode = np.loadtxt(path+'Coord_Zigbee_Nodes.csv',delimiter=',')
# Dongle 
Dongle = np.loadtxt(path+'Coord_dongle.csv',delimiter=',')

# <codecell>

fig,ax=L.showG('s',nodes=False)
axis('off')
S.tx.clear()
S.rx.clear()
da ={}
dm ={}
for c,k in enumerate(AnchorNodes):
    pta = array([AnchorNodes[k]['coord'][0],AnchorNodes[k]['coord'][1],AnchorNodes[k]['coord'][2]]).reshape(3,1)
    S.tx.point(pta,mode="add")
    da[c]=k
    plt.plot(pta[0,:],pta[1,:],'or')
for c,k in enumerate(MT_DLR_RDTMaster):
    ptm = array([MT_DLR_RDTMaster[k][0],MT_DLR_RDTMaster[k][1],MT_DLR_RDTMaster[k][2]]).reshape(3,1)
    dm[c]=k
    S.rx.point(ptm,mode="add")
    plt.plot(ptm[0,:],ptm[1,:],'ob')
#plt.scatter(Dongle[:,1],Dongle[:,2])
#plt.scatter(ZigbeeNode[:,1],ZigbeeNode[:,2],color='r')

# <codecell>

print S.tx.position

# <codecell>

print S.rx.position

# <codecell>

S.show(s=20)

# <codecell>

p=S.L.Gr.node[0]['polyg']

# <codecell>

tx = S.tx.position[:,4]
Rtx = S.L.pt2ro(tx)
print "transmitter :",tx," is in room ",Rtx

# <codecell>

rx = S.rx.position[:,28]
Rrx = S.L.pt2ro(rx)
print "mobile node :",rx," is in room ",Rrx

# <codecell>

ctx = L.pt2cy(tx)
crx = L.pt2cy(rx)
Si = Signatures(S.L,ctx,crx)

# <codecell>

Si.run1(cutoff=4)
Si

# <codecell>

r2d=Si.rays(tx,rx)

# <codecell>

r2d.show(L)

# <codecell>

def showr2d(L,r2d,tx,rx):
    """
    r2d['pt'] : nd,ni,nr
    """
    L.display['thin']=True
    col = ['r','b','g','c','m','k','y']
    fig,ax = L.showGs()
    for k in r2d:
        r = r2d[k]
        pts = r['pt']
        sh = np.shape(pts)
        for r in range(sh[2]):
            x = np.hstack((tx[0],pts[0,:,r],rx[0]))
            y = np.hstack((tx[1],pts[1,:,r],rx[1]))
            plt.plot(x,y,col[k])
    

# <codecell>

showr2d(L,r2d,tx,rx)

# <codecell>

S.L.buildGi()
S.L.Gi.edge

# <codecell>

# ACOLinkMeasurementID;NodeID;ACO_PeerID;ACO_RSSI;Timestamp
fd = open(path+'mes/ACO_link_meas.txt')
lis = fd.readlines()
#Nodes = np.loadtxt('mes/Nodes.txt',delimiter=';')
tl  = {}
for li in lis:
    li = li.replace('\r\n','')
    s = li.split(';')
    try:
        date = datetime.strptime(s[4],'%Y-%m-%d %H:%M:%S')
        timestamp = time.mktime(date.timetuple())
    except:
        timestamp = 0
    cle= s[1]+'-'+s[2]
    c = np.array([timestamp,eval(s[3])])
    #c = np.array([timestamp,eval(s[0]),eval(s[1]),eval(s[2]),eval(s[3])])
    try:
        tl[cle] = np.vstack((tl[cle],c))
    except:
        tl[cle] = c

# <codecell>

print "first data",lis[0]
print "last data", lis[-1]
print "available links", tl.keys()
print "raw integer timestamp",tl['387-389'][:,0]
print "time in second",tl['387-389'][:,0]-tl['387-389'][0,0]

N=20243/3600.
Nh=floor(N)
Nm=(20243-Nh*3600)/60 
print N
print Nh,"hours and",Nm,"minutes"

# <codecell>

tstart = 10+48/60.
for k in range(36):
    cle = tl.keys()[k]
    L1 = cle.split('-')[0]
    L2 = cle.split('-')[1]
    plt.subplot(6,6,k+1)
    ts = (tl[cle][:,0]-tl[cle][0,0])
    th = tstart+ ts/3600
    plt.plot(th,tl[cle][:,1],'.')
    plt.axis([10,17,-84,-38])
    plt.title(cle)
#plt.savefig('data.png')

# <codecell>

#plt.figure()
#plt.plot(tc[:,0],tc[:,1],'.')
#plt.xlabel('Integer timestamp = time.mktime(date.timetuple())')
#plt.ylabel('Received Power (dBm)')
#plt.show()
#plt.figure()
#plt.hist(tc[:,1],150)
#plt.title('Empirical distribution of Received Power - All data')
#plt.show()
plt.figure(figsize=(10,10))
v1 = tl['391-390']
v2 = tl['391-386']
plt.subplot(211)
plt.plot((v1[:,0]-v1[0,0])/100,v1[:,1],'.')
plt.plot((v2[:,0]-v2[0,0])/100,v2[:,1]+0.5,'.r')
plt.xlabel('Normalized timestamp')
plt.ylabel('Received Power (dBm)')
plt.title('comparison of recorded activity on 2 correlated multi stable links')
plt.legend(('Link 391-390', 'Linkl 391-386 + 0.5dBm(offset) ' ),
                      'upper right', shadow=True, fancybox=True)
plt.subplot(212)
plt.hist(v1[:,1],50,alpha=0.5)
plt.hist(v2[:,1],50,alpha=0.5,color='r')
plt.xlabel('Power (dBm)')
plt.savefig('391-390vs391-386.png')

# <codecell>

itx=1
irx=1
S.tx = RadioNode(typ='tx')
S.tx.point([1.2,1,1.4])
S.rx = RadioNode(typ='rx')
S.rx.point([8,-1.2,1.5])
S.save()

# <markdowncell>

# Run raytracing determination between itx and irx  

# <codecell>

S.run(itx,irx)

# <markdowncell>

# Load Tud object from generated file. A Tud object is gathering all the necessary information for avaluating the field along the rays of the link (itx,irx)

# <codecell>

Gt=GrRayTud()
Gt.load(S.dtud[itx][irx],S.dtang[itx][irx],S.drang[itx][irx],S.slab)

# <markdowncell>

# # How to use the GrRayTud object? 

# <markdowncell>

# ## Gt is an object from class GrRayTud (Group of rays in TUD format) 

# <markdowncell>

# ### Gt.eval
# Thhis method evaluates the field along the rays. 

# <codecell>

Gt.eval()

# <markdowncell>

# ### Gt.info(r)
# 
# returns the information associated to a given ray number r

# <codecell>

Gt.info(7)

# <markdowncell>

# ### Gt.ray(r)
# 
# return the index of the interactions of ray r into the Gt.I.I matrix 

# <codecell>

Gt.ray(0)

# <markdowncell>

# ### Gt.typ(r)
# 
# return the type of the interactions of ray r 

# <codecell>

Gt.typ(0)

# <markdowncell>

# ## Gt Attributes

# <markdowncell>

# ### Gt.I

# <markdowncell>

# Gt has an attribute I, which is a class Interactions.
# 
# This class Interactions gather all basis interactions : B / T / R / L 

# <codecell>

print Gt.I.B
print Gt.I.R
print Gt.I.T
print Gt.I.L

# <markdowncell>

# All basis interactions B / T / R / L have the same attributes:
# 
# * idx : the index of the interaction 
# * data : which is a np.shape(self.data) = len(self.idx) , 3 
# * data[:,0] = theta
# * data[:,1] = si
# * data[:,2] = sr or st ( named sout in the following)

# <codecell>

print Gt.I.R.idx,'\n'
print Gt.I.R.data

# <markdowncell>

# T and R basis interractions have also an extra attribute:
# 
# * dusl : a dictonnary of used slab ( the number is the **position** into the self.idx, not direcly the index)

# <codecell>

print Gt.I.R.dusl

# <markdowncell>

# ### Gt.dli

# <codecell>

Gt.dli

# <markdowncell>

# the Gt object has an attribute : **dli => dictionnary of length of interaction**
# 
# The key of this dictionnary is the number of interaction

# <codecell>

Gt.dli.keys()

# <markdowncell>

# Thus, Lets see what contains Gt.dli for 3 interractions

# <codecell>

Gt.dli[3]

# <markdowncell>

# This is a new dictionnary which gives:
# 
# * 'nbrays' : The number of ray which have 3 interractions ( here only 1)
# 
# * 'rayidx' : The index of the ray ( here only the ray number 0)
# 
# * 'rays' : a np array which contains the index of the interraction maxtrix

