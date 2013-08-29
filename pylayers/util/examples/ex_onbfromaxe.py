from pylayers.util.geomutil import *

A = np.array([[0,0,0,0],[1,2,3,4],[0,0,0,0]])
B = np.array([[0,0,0,0],[1,2,3,4],[10,10,10,10]])
T = onbfromaxe(A,B)

for k in range(4):
    filename = 'onbasis'+str(k) 
    gv = GeomVect(filename)
    gv.geomBase(T[k,:,:])
   
