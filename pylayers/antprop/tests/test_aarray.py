from pylayers.antprop.aarray import *
import matplotlib.pyplot as plt 
print '--------------'
print 'antprop/test_aarray.py'
print '--------------'
A = AntArray(tarr='UA',N=[10,10,1],dm=[0.5,0.5,0],fGHz=0.3)
A.eval()
f1,a1 = A.plotG()
f1.savefig('ArrayDiag.png')
f2,a2 = A.show()
f2.savefig('ArrayConfig.png')
