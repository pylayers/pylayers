from pylayers.util.geomutil import *

line0 = (np.array([0,0]),np.array([1,1]))
line1 = (np.array([0,0]),np.array([1,-1]))
#
#
pta = np.array([0,0])
phe = np.array([3,1])
seg  = (pta,phe)
pdb.set_trace()
tahe,r =  intersect_cone_seg(line0,line1,seg,bvis=True,bbool=True)
#
#
#
pta = np.array([10,3])
phe = np.array([-10,2.9])
seg  = (pta,phe)
tahe,r =  intersect_cone_seg(line0,line1,seg,bvis=True,bbool=True)
#print tahe,r
print seg[0],seg[1],tahe,r
#
# case baa 
#
pta = np.array([3,1001])
phe = np.array([3,1000])
seg  = (pta,phe)
tahe,r =  intersect_cone_seg(line0,line1,seg,bvis=False)
print seg[0],seg[1],tahe,r
#assert tahe==[]
#assert r==0
seg  = (phe,pta)
tahe,r =  intersect_cone_seg(line1,line0,seg,bvis=False)
print seg[0],seg[1],tahe,r
#assert tahe==[]
#assert r==0
#
# case bai 
#
pta = np.array([3,10])
phe = np.array([3,0])
seg = (phe,pta)
tahe,r = intersect_cone_seg(line0,line1,seg,bvis=False)
print seg[0],seg[1],tahe,r
#print tahe,r
#assert np.isclose((tahe[0]-np.array([3,3])).all(),0)
#assert np.isclose((tahe[1]-np.array([3,0])).all(),0)
#assert np.isclose(r,0.5)
#print tahe,r
seg  = (pta,phe)
tahe,r = intersect_cone_seg(line0,line1,seg,bvis=False)
print seg[0],seg[1],tahe,r
#assert np.isclose((tahe[0]-np.array([3,3])).all(),0)
#assert np.isclose((tahe[1]-np.array([3,0])).all(),0)
#assert np.isclose(r,0.5)
#print tahe,r
#
#
#
pta = np.array([3,3])
phe = np.array([3,-3])
seg  = (phe,pta)
tahe,r = intersect_cone_seg(line0,line1,seg,bvis=False)
print seg[0],seg[1],tahe,r
#assert np.isclose((tahe[0]-np.array([3,3])).all(),0)
#assert np.isclose((tahe[1]-np.array([3,-3])).all(),0)
#print tahe,r
#assert np.isclose(r,1)
pta = np.array([3,2])
phe = np.array([3,-2])
seg  = (phe,pta)
tahe,r =  intersect_cone_seg(line0,line1,seg,bvis=False)
print seg[0],seg[1],tahe,r
#assert np.isclose((tahe[0]-np.array([3,2])).all(),0)
#assert np.isclose((tahe[1]-np.array([3,-2])).all(),0)
#assert np.isclose(r,2/3.)
pta = np.array([3,3])
phe = np.array([3,4])
seg  = (phe,pta)
tahe,r =  intersect_cone_seg(line0,line1,seg,bvis=False)
print seg[0],seg[1],tahe,r
#assert np.isclose((tahe[0]-np.array([3,2])).all(),0)
#assert np.isclose((tahe[1]-np.array([3,-2])).all(),0)
#assert np.isclose(r,0)
#print tahe,r
#seg2  = (phe2,pta2)
#tahe,r2 =  intersect_cone_seg(line0,line1,seg2)
#print tahe , r
#seg3  = (phe3,pta3)
#tahe,r3 =  intersect_cone_seg(line0,line1,seg3)
#print tahe ,r3
#seg4  = (phe4,pta4)
#tahe,r4 =  intersect_cone_seg(line0,line1,seg4)
#print tahe ,r4
#seg6  = (phe6,pta6)
#tahe,r6 =  intersect_cone_seg(line0,line1,seg6)
#print tahe,r6
##
###[array([ 3.,  3.]), array([3, 0])] 0.5             
###[array([3, 0]), array([3, 2])] 0.333333333333 
###[array([ 3., -3.]), array([3, 0])] 0.5      
###[array([ 3.,  3.]), array([ 3., -3.])] 1.0         
###[array([ 3.,  3.]), array([ 3., -3.])] 1.0 
line0 = (np.array([0,0]),np.array([10,0]))
line1 = (np.array([0,0]),np.array([10,-10]))
pta = np.array([3,-3])
phe = np.array([30,-3])
seg  = (phe,pta)
tahe,r = intersect_cone_seg(line0,line1,seg,bvis=False)
print seg[0],seg[1],tahe,r
#assert np.isclose((tahe[0]-np.array([3,-3])).all(),0)
#assert np.isclose((tahe[1]-np.array([30,-3])).all(),0)
#assert np.isclose(r,1)
phe = np.array([30,-2.999999])
seg = (pta,phe)
tahe,r =  intersect_cone_seg(line0,line1,seg,bvis=False,bbool=False)
print seg[0],seg[1],tahe,r
#assert np.isclose((tahe[0]-np.array([3,-3])).all(),0)
#assert np.isclose((tahe[1]-np.array([30,-2.999999])).all(),0)

line0 = (np.array([ 256.85347058,  192.31004176]), np.array([ 0.97752571, -0.21081625]))
line1 = (np.array([ 256.85347058,  192.31004176]), np.array([ 0.89570579, -0.44464721]))
seg   = (np.array([ 264.925,  189.851]), np.array([ 261.676,  191.27 ]))
print seg[0],seg[1],tahe,r
#pdb.set_trace()
tahe,r = intersect_cone_seg(line0,line1,seg,bvis=False,bbool=False)
#print tahe,r
line0 = (np.array([ 190.528,  206.614]), np.array([ 0.3597124 , -0.93306323])) 
line1 = (np.array([ 190.528,  206.614]), np.array([-0.91274561, -0.40852839]))
seg =(np.array([ 173.585,  196.352]), np.array([ 174.041,  195.169]))  
tahe,r = intersect_cone_seg(line0,line1,seg,bvis=False,bbool=False)
print seg[0],seg[1],tahe,r
#print tahe,r


line0 = (np.array([ 210.254,  215.443]), np.array([-0.99821912,  0.05965383]))
line1 = (np.array([ 210.254,  215.443]), np.array([ 0.91274561,  0.40852839]))
seg = (np.array([ 199.293,  243.879]), np.array([ 205.758,  227.107]))
tahe,r = intersect_cone_seg(line0,line1,seg,bvis=False,bbool=False)
print seg[0],seg[1],tahe,r
