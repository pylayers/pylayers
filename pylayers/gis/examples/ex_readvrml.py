from pylayers.gis.readvrml import *
import pylayers.util.pyutil as pyu
import matplotlib.pyplot as plt
import numpy as np

filename = pyu.getlong('B11D-E1.wrl','struc')
VL = VLayout()
VL.load(filename)
plt.figure(figsize=(15,5))

# 0 1 WALL 2 3 4
VL.show(1)
#
# Analysis of walls
#
dwall = VL.wallanalysis()
#
# Visualization
#
dpt = {}
k=0
for iw in dwall:
    seg = dwall[iw]['seg']
    thick = dwall[iw]['thickness']
    bdoor = dwall[iw]['door']
    x,y = seg.xy
    pt = np.array((x[0],y[0]))
    ph = np.array((x[1],y[1]))
    dpt[2*k]=pt
    dpt[2*k+1]=ph+1
    k = k+2
    dwall[iw]['tail']=2*k
    dwall[iw]['head']=2*k+1
    try:
        tpt = np.vstack((tpt,pt,ph))
    except:
        tpt = np.vstack((pt,ph))

    if bdoor:
        plt.plot(x,y,color='r',linewidth=thick*10,alpha=1)
    else:
        plt.plot(x,y,color='k',linewidth=thick*10,alpha=1)

plt.axis('scaled')
savestr2(dpt,dwall,_filename='struc.str2')
#vrml2geom(tg,'11DE1')
#lfi = glob.glob('/private/staff/n/en/buguen/Pyproject/struc/bat11/*BAT11D*E1.wrl')
#for fi in lfi:
#    if fi.find('furniture')==-1:
#        tmp  = fi.split('/')
#        name = tmp[-1].split('.')[0]
#        print name
#        tg   = parsevrml(fi)
#        vrml2geom(tg,name)
#os.system('ls /private/staff/n/en/Pyproject/struc/bat11')
#
##self.entity = {}
#for t in dg:     # WALL, COLUMN, DOOR , STAIR , SPACE
#    self.entity[t] = {}
#    for ID in dg[t]:
#        c = dg[t][ID]['coord']
#        self.entity[t][ID] = {}
#        self.entity[t][ID]['coord']=c
#        l = dg[t][ID]['index']
#        dp = {}
#        p  = []
#        k  = 0
#        for il in l:
#            if il == -1:
#                dp[k] = p
#                p     = []
#                k     = k + 1
#            else:
#                p.append(il)
#
#        self.entity[t][ID]['index'] = dp
#
#
#plt.axis('equal')
##doctest.testmod()

a = 0.07
rx  = np.array([[0,0,0,0,0,0,0,0],[0,a,2*a,3*a,4*a,5*a,6*a,7*a],[1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]])
Orx = np.array([11.5,6.13,0])
Rx  = rx + np.outer(Orx,np.ones(8))
plt.plot(Rx[0,:],Rx[1,:],'or')





