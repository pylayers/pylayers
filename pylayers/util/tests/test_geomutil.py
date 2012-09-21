import os
import scipy as sp 
import numpy as np 
from pylayers.util.geomutil import *
import matplotlib.pylab as plt 

gv = GeomVect()

#
# Displaying 2D or 3D segments 
#
ds2 = {1:([0,0],[1,1]),
       2:([1,1],[1,0])}

ds3 = {1:([0,0,0],[1,1,4]),
       2:([1,1,4],[1,0,8])}
gv.segments(ds2)
gv.show3()

gv.segments(ds3,False)
gv.show3()
#
#
#
gv.points([0,0,0])
#gv.show3()
#M  = sp.rand(3,3)
#pt = np.array([0.,0.,0.])
#filebase = "geomBase_demo"
#fileell  = "ellipse"
#col = np.array([[1,0.25,0.25],
#            [1,0.5,0.5],
#            [1,0.75,0.75]]) 
#
#g1 = GeomVect(filebase)
#g1.geomBase(M,pt,col)
##
#g2 = GeomVect(fileell) 
#p  = np.array((0,0,1))
#N  = 20 
#Eth = 1+1j
#Eph = 3-1j
##       
#th  = np.pi/4
#ph  = np.pi/4
##
#B   = vec_sph(th,ph)
#os.system("geomview geomBase_demo.vect 2>/dev/null &")
#ellipse(B[2,:],B[0,:],B[1,:],Eth,Eph,N)
##
#os.system("geomview -b 1 1 1 geomBase_demo.vect 2>/dev/null &")
###
### Angles of Arrival  (Global Frame)
###
#ag = array(((0,0.),(0,pi/3),(pi/4,pi/4),(pi/5,pi/5)))
#ag = pi*sp.rand(100,2)
#N  = shape(ag)[0]
#M  = SphericalBasis(ag)
##
#tha = M[0,:,:]
#pha = M[1,:,:]
#sa  = M[2,:,:]
#
#Ba_Rg  = dstack((tha,pha)).transpose((0,2,1))
##
##   Construction d'une matrice de rotation de l'antenne Ta 
##
#alpha = np.pi/2 
#beta  = np.pi/3 
#gamma = np.pi/4 
#
#Ta = MEulerAngle(alpha,beta,gamma)
##
#tha_Ra = np.dot(Ta.T,M[0,:,:]).T
#pha_Ra = np.dot(Ta.T,M[1,:,:]).T
#sa_Ra  = np.dot(Ta.T,M[2,:,:]).T
##
##
## On calcule les angles a partir de sa_Ra
##
#a_new = AngleDir(sa_Ra)
#K     = SphericalBasis(a_new)
#tha2  = K[0,:,:]
#pha2  = K[1,:,:]
##
#Ba_RaT =  dstack((tha2,pha2)).transpose((2,0,1))
#R2 = mul3(Ba_RaT,mul3(Ta.T,Ba_Rg))
##
##
#p  = np.array([[0.,0.],[5.,0.],[10.,0.],[10.,-2.],[0.,-2.]])
#lr = shg.asLineString(p)
#pol1 = shg.Polygon(lr)
#pol2 = shrinkPolygon(pol1)
#pol2 = simplifyPolygon(pol1)
