from pylayers.util.geomutil import *
import numpy as np 
Nplane = 1
r0 = np.array([[1,2/3.,1]]).T
r1 = np.array([[1,1/3.,2]]).T
p  = np.concatenate((r0[:,:,None],r0[:,:,None]),axis=2)
p1 = np.array([[0,0] ,[1,0],[0,1]])  #yz
p2 = np.array([[1,0] ,[0,0],[0,1]])  #xz 
#p3 = np.array([[1,0] ,[0,1],[0,0]])  #xy 
aplane = np.concatenate((p1[:,None,:],p2[:,None,:]),axis=1) 
#aplane = p1
q1 = np.array([[0,0,0]]).T
q2 = np.array([[0,1,0]]).T
pplane = np.hstack((q1,q2)) 
#pplane = q1 
Ip1 = mirror3(r0,aplane,pplane)
Ip2 = mirror3b(p,aplane,pplane)
