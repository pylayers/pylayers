from pylayers.util.geomutil import *
import numpy as np 
Nplane = 1
r0 = np.array([[1,2/3.,1]]).T
r1 = np.array([[1,1/3.,2]]).T
p  = np.concatenate((r0[:,:,None],r1[:,:,None]),axis=2)
pp = p[...,None] 
p1 = np.array([[0,0] ,[1,0],[0,1]])  #yz
p2 = np.array([[1,0] ,[0,0],[0,1]])  #xz 
p3 = np.array([[1,0] ,[0,0],[0,1]])  #xz 
p4 = np.array([[0,0] ,[1,0],[0,1]])  #yz
#p3 = np.array([[1,0] ,[0,1],[0,0]])  #xy 
aplane1 = np.concatenate((p1[:,None,:],p2[:,None,:]),axis=1) 
aplane2 = np.concatenate((p3[:,None,:],p4[:,None,:]),axis=1) 
#aplane = p1
q1 = np.array([[0,0,0]]).T
q2 = np.array([[0,1,0]]).T
q3 = np.array([[1,1,0]]).T
q4 = np.array([[2,0,0]]).T
pplane1 = np.hstack((q1,q2)) 
pplane2 = np.hstack((q3,q4)) 
#pplane = q1 
Ip1 = mirror3(r0,aplane1,pplane1)
Ip2 = mirror3b(p,aplane1,pplane1)
# s x v  x p x n 
# space x vector x plane x nsig 
aplane = np.concatenate((aplane1[...,None],aplane2[...,None]),axis=3)
pplane = np.concatenate((pplane1[...,None],pplane2[...,None]),axis=2)
Ip3 = mirror3c(pp,aplane,pplane)
