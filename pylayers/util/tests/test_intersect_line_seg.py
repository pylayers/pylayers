from pylayers.util.geomutil import *

line = (np.array([0,0]),np.array([3,0]))
pta = np.array([3,1])
phe = np.array([4,1])
seg  = (pta,phe)

x, P   =  intersect_line_seg(line,seg)

print "pta : ",pta
print "phe : ",phe
print x,P
