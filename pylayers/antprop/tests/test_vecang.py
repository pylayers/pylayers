#!/usr/bin/python
# -*- coding: latin1 -*-
import pylayers.util.geomutil as geu
import numpy as np


A= np.array(([-1,-1]))
B= np.array(([1,-1]))
C= np.array(([1,1]))
D= np.array(([-1,1]))


vAB = (B-A)/ np.sqrt(np.sum((B-A)*(B-A)))
vBA = (A-B)/ np.sqrt(np.sum((A-B)*(A-B)))

vAD = (D-A)/ np.sqrt(np.sum((D-A)*(D-A)))
vDA = (A-D)/ np.sqrt(np.sum((A-D)*(A-D)))
    
vBC = (C-B)/ np.sqrt(np.sum((C-B)*(C-B)))
vCB = (B-C)/ np.sqrt(np.sum((B-C)*(B-C)))

vCD = (D-C)/ np.sqrt(np.sum((D-C)*(D-C)))
vDC = (C-D)/ np.sqrt(np.sum((C-D)*(C-D)))

#input
vBDu = np.array([-1,1])
vBDu = (vBDu)/ np.sqrt(np.sum((vBDu)*(vBDu)))

vCAu = np.array([-1,-1])
vCAu = (vCAu)/ np.sqrt(np.sum((vCAu)*(vCAu)))

vACu = np.array([1,1])
vACu = (vACu)/ np.sqrt(np.sum((vACu)*(vACu)))

vDBu = np.array([1,-1])
vDBu = (vDBu)/ np.sqrt(np.sum((vDBu)*(vDBu)))


#point A
#diff dict : keys are diffraction point
#values are the 2 possible 0 or n associated face
diff={'A':['vBA','vDA'],
      'B':['vAB','vCB'],
      'C':['vDC','vBC'],
      'D':['vCD','vAD'] }   
#tri is the associated edge pa--pt--pb
tri = {'A':['BAD','DAB'],
       'B':['ABC','CBA'],
       'C':['DCB','BCD'],
       'D':['CDA','ADC']}   

print ' D          C'
print '  x--------x'
print '  |        |'
print '  |        |'
print '  |        |'
print '  x--------x'
print ' A          B'

for d in diff.keys():
    for it,t in enumerate(tri[d]):
        print '\ndiff point :',d,' tri:',t
        x=eval(t[0])
        y=eval(t[1])
        z=eval(t[2])
        #is left(a,b,c) (Test point c is at left of the vector a-->b)
        uleft = geu.isleft(x[:,None],y[:,None],z[:,None])[0]
        if uleft : 
            print t[2],' is  at left of vect ',t[0],t[1]
            print '-------------------------------------'
            vin = eval(diff[d][it])
            print diff[d][it],u'\u2196',geu.vecang(vin,vBDu)*180/np.pi
            print diff[d][it],u'\u2199',geu.vecang(vin,vCAu)*180/np.pi
            print diff[d][it],u'\u2197',geu.vecang(vin,vACu)*180/np.pi
            print diff[d][it],u'\u2198',geu.vecang(vin,vDBu)*180/np.pi
        else: 
            print t[2],' is  at right of vect ',t[0],t[1]
            print '-------------------------------------'
            vin = eval(diff[d][it])
            print diff[d][it],u'\u2196',geu.vecang(vBDu,vin)*180/np.pi
            print diff[d][it],u'\u2199',geu.vecang(vCAu,vin)*180/np.pi
            print diff[d][it],u'\u2197',geu.vecang(vACu,vin)*180/np.pi
            print diff[d][it],u'\u2198',geu.vecang(vDBu,vin)*180/np.pi

        
#     if uleft : 
#         print geu.vecang(vBA,vBDu)*180/np.pi
#         print geu.vecang(vBA,vCAu)*180/np.pi
#         print geu.vecang(vBA,vACu)*180/np.pi
#         print geu.vecang(vBA,vDBu)*180/np.pi
#     else : 
#         print geu.vecang(vBDu,vBA)*180/np.pi
#         print geu.vecang(vCAu,vBA)*180/np.pi
#         print geu.vecang(vACu,vBA)*180/np.pi
#         print geu.vecang(vDBu,vBA)*180/np.pi



# uleft = geu.isleft(B[:,None],A[:,None],D[:,None])[0]
# if uleft : 
#     print geu.vecang(vBA,vBDu)*180/np.pi
#     print geu.vecang(vBA,vCAu)*180/np.pi
#     print geu.vecang(vBA,vACu)*180/np.pi
#     print geu.vecang(vBA,vDBu)*180/np.pi
# else : 
#     print geu.vecang(vBDu,vBA)*180/np.pi
#     print geu.vecang(vCAu,vBA)*180/np.pi
#     print geu.vecang(vACu,vBA)*180/np.pi
#     print geu.vecang(vDBu,vBA)*180/np.pi