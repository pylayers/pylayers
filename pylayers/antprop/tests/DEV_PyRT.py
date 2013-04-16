# -*- coding: utf-8 -*-

from pylayers.simul.simulem import *
from pylayers.antprop.rays import *
from pylayers.antprop.channel import *
from pylayers.antprop.signature import *
import pylayers.util.pyutil as pyu
from pylayers.gis.layout import *
import pylayers.signal.bsignal as bs
from datetime import datetime
import time
import ConfigParser
import pickle
def showr2d(L,r2d,tx,rx):
    """
    r2d['pt'] : nd,ni,nr
    """
    L.display['thin']=False
    col = ['r','b','g','c','m','k','y']
    fig,ax = L.showG('')
    for k in r2d:
        r = r2d[k]
        pts = r['pt']
        sh = np.shape(pts)
        for r in range(sh[2]):
            x = np.hstack((tx[0],pts[0,:,r],rx[0]))
            y = np.hstack((tx[1],pts[1,:,r],rx[1]))
            plt.plot(x,y,col[eval(k)])

def showr(L,r2d,tx,rx,k,l,color='b'):
    """
    r2d['pt'] : nd,ni,nr
    """
    L.display['thin']=False
    fig,ax = L.showG('')
    r = r2d[str(k)]
    pts = r['pt']
    sh = np.shape(pts)
    x = np.hstack((tx[0],pts[0,:,l],rx[0]))
    y = np.hstack((tx[1],pts[1,:,l],rx[1]))
    plt.plot(x,y,color=color)




S = Simul()
filestr = 'DLR'
S.layout(filestr+'.ini','matDB.ini','slabDB.ini')
try:
    S.L.dumpr()
except:
    S.L.build()
    S.L.dumpw()
S.L.display['ednodes']=False
S.L.display['nodes']=False
S.L.display['title']='DLR WP4 WHERE2 measurement site'
S.L.display['overlay']=False
S.L.showGs() 
S.L.showG('i')


MT_DLR_RDTMaster={ 2:[2.12,0,1.275], 3:[3.12,0,1.275], 4:[4.12,0,1.275], 5:[5.12,0,1.275], 6:[6.12,0,1.275], 7:[7.12,0,1.275], 8:[8.12,0,1.275], 9:[9.12,0,1.275], 10:[10.12,0,1.275], 11:[11.12,0,1.275], 12:[12.12,0,1.275], 13:[13.12,0,1.275], 14:[14.12,0,1.275], 15:[15.12,0,1.275], 16:[16.12,0,1.275], 17:[17.12,0,1.275], 18:[18.12,0,1.275], 19:[19.12,0,1.275], 20:[20.12,0,1.275], 21:[21.12,0,1.275], 22:[22.12,0,1.275], 23:[23.12,0,1.275], 24:[24.12,0,1.275], 25:[25.12,0,1.275], 26:[26.12,0,1.275], 27:[27.12,0,1.275], 28:[28.12,0,1.275], 29:[29.12,0,1.275], 30:[30.12,0,1.275], 31:[31.12,0,1.275], 33:[30.12,0,1.275], 62:[24.12,0,1.275], 63:[24.62,0,1.275], 64:[25.12,0,1.275], 65:[25.62,0,1.275], 66:[26.12,0,1.275], 67:[26.62,0,1.275], 68:[27.12,0,1.275], 69:[27.62,0,1.275], 70:[28.12,0,1.275], 71:[28.62,0,1.275], 72:[29.12,0,1.275], 73:[29.62,0,1.275], 34:[30.12,0,1.275], 35:[30,1.38,1.275], 36:[30,1.88,1.275], 37:[30,2.38,1.275], 38:[30,2.88,1.275], 39:[30,3.88,1.275], 40:[30,4.88,1.275], 41:[30,5.88,1.275]}


TrolleyMT_ACO_04={ 2:60, 3:59, 4:58, 5:58, 6:58, 7:58, 8:58, 9:58, 10:58, 11:58, 12:58,
13:57, 14:52, 15:12, 16:51, 17:11, 18:50, 19:10, 20:53, 21:54, 22:55, 23:56, 24:56, 25:56,
26:56, 27:56, 28:56, 29:56, 30:56, 31:56, 33:56, 62:56, 63:56, 64:56, 65:56, 66:56, 67:56,
68:56, 69:56, 70:56, 71:56, 72:56, 73:56, 34:56, 35:56, 36:56, 37:56, 38:56, 39:56, 40:56, 41:56} 


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

Dongle = 389
AnchorNodes = {390:{'name':'MT_ACO_05','coord':[6,0.81,1.64]},
               386:{'name':'MT_ACO_08','coord':[30.574,2.8,1.291]},
               391:{'name':'MT_ACO_07','coord':[11.78,-5.553,1.5]},
               385:{'name': 'MT_ACO_01','coord':[19.52,-0.69,1.446]},
               387:{'name':'MT_ACO_03','coord':[28.606,-0.74,1.467]},
               400:{'name':'MT_ACO_02','coord':[30.574,2.8,1.291]},
               1:{'name':'MT_DLR_RTDSlave','coord':[0.85,0,1.18]}
              }

#
#  Define the link to selected
#
A = 390
B = 386


L = Layout('DLR.ini')
try:
    L.dumpr()
except:
    L.build()
    L.dumpw()
L.showG('r')
#path = '/private/staff/n/en/buguen/Bureau/WHERE2/WHERE2-DLR-Simulations/'
#ZigbeeNode = np.loadtxt(path+'Coord_Zigbee_Nodes.csv',delimiter=',')
#Dongle = np.loadtxt(path+'Coord_dongle.csv',delimiter=',')

#
# Reading Anchor Nodes and Mobile Nodes
#
S.tx.clear()
S.rx.clear()
da ={}
dm ={}
for c,k in enumerate(AnchorNodes):
    pta = array([AnchorNodes[k]['coord'][0],AnchorNodes[k]['coord'][1],AnchorNodes[k]['coord'][2]]).reshape(3,1)
    S.tx.point(pta,mode="add")
    da[c]=k
    # plt.plot(pta[0,:],pta[1,:],'or')
for c,k in enumerate(MT_DLR_RDTMaster):
    ptm = array([MT_DLR_RDTMaster[k][0],MT_DLR_RDTMaster[k][1],MT_DLR_RDTMaster[k][2]]).reshape(3,1)
    dm[c]=k
    S.rx.point(ptm,mode="add")
    # plt.plot(ptm[0,:],ptm[1,:],'ob')

#
# Select a link 
#
tx = S.tx.position[:,4]
Rtx = S.L.pt2ro(tx)
print "transmitter :",tx," is in room ",Rtx

rx = array([15,3,2.5])
#rx = S.rx.position[:,28]
Rrx = S.L.pt2ro(rx)
print "mobile node :",rx," is in room ",Rrx

print tx
print rx
a=time.time()

if not os.path.exists('r2d.pickle'):
    Si = Signatures(S.L,tx,rx)
    Si.run(tx,rx,4)
    r2d = Si.rays()
    file=open("r2d.pickle","w")
    pickle.dump(r2d,file)
    file.close()
else:
    file = open("r2d.pickle","r")
    r2d = pickle.load(file)
r3d = r2d.to3D()
r3d.locbas(L)
r3d.fillinter(L)
#<<<<<<< HEAD
#r3d.eval()

#c11 = r3d.Ctilde[:,:,0,0]
#c12 = r3d.Ctilde[:,:,0,1]
#c21 = r3d.Ctilde[:,:,1,0]
#c22 = r3d.Ctilde[:,:,1,1]



#Cn=Ctilde()
#Cn.Cpp = bs.FUsignal(r3d.I.f, c11)
#Cn.Ctp = bs.FUsignal(r3d.I.f, c12)
#Cn.Cpt = bs.FUsignal(r3d.I.f, c21)
#Cn.Ctt = bs.FUsignal(r3d.I.f, c22)
#Cn.nfreq = r3d.I.nf
#Cn.nray = r3d.nray
#Cn.tauk=r3d.delays

#raynumber = 4

#fig=plt.figure('Cpp')
#f,ax=Cn.Cpp.plot(fig=fig,iy=np.array(([raynumber])))

#r3d.info(raynumber)
## plt.show()

#=======

config = ConfigParser.ConfigParser()
_filesimul = 'default.ini'
filesimul = pyu.getlong(_filesimul, "ini")
config.read(filesimul)
fGHz = np.linspace(eval(config.get("frequency", "fghzmin")), 
                     eval(config.get("frequency", "fghzmax")), 
                     eval(config.get("frequency", "nf")))

r3d.eval(fGHz)
#
#c11 = r3d.Ctilde[:,:,0,0]
#c12 = r3d.Ctilde[:,:,0,1]
#c21 = r3d.Ctilde[:,:,1,0]
#c22 = r3d.Ctilde[:,:,1,1]
#
#
#
#Cn=Ctilde()
#Cn.Cpp = bs.FUsignal(r3d.I.f, c11)
#Cn.Ctp = bs.FUsignal(r3d.I.f, c12)
#Cn.Cpt = bs.FUsignal(r3d.I.f, c21)
#Cn.Ctt = bs.FUsignal(r3d.I.f, c22)
#Cn.nfreq = r3d.I.nf
#Cn.nray = r3d.nray
#Cn.tauk=r3d.delays
#
#raynumber = 4
#
#fig=plt.figure('Cpp')
#f,ax=Cn.Cpp.plot(fig=fig,iy=np.array(([raynumber])))
#

