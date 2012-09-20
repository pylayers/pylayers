from pylayers.gis.lirevrml import *
import pylayers.util.pyutil as pyu
import matplotlib.pyplot as plt 
import numpy as np 

filename = pyu.getlong('B11D-R0.wrl','struc')
VL = VLayout()
VL.load(filename)
plt.figure(figsize=(15,5))
# 0 1 WALL 2 3 4 
VL.show(1)
dwall = VL.wallanalysis()
for iw in dwall:
    seg = dwall[iw]['seg']
    thick = dwall[iw]['thickness']
    bdoor = dwall[iw]['door']
    x,y = seg.xy
    if bdoor:
        plt.plot(x,y,color='r',linewidth=thick*10,alpha=1)
    else:    
        plt.plot(x,y,color='k',linewidth=thick*10,alpha=1)
plt.axis('scaled')    
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
rx  = np.array([[0,0,0,0,0,0,0,0],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7],[1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]])
Orx = np.array([10.64-0.78,6.13,0])
Rx  = rx + np.outer(Orx,np.ones(8))
plt.plot(Rx[0,:],Rx[1,:],'or')
