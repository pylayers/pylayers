import matplotlib
from pylayers.antprop.aarray import *
print '--------------'
print 'antprop/test_aarray.py'
print '--------------'
A = AntArray()
A.eval()
f1,a1 = A.plotG()
f1.savefig('ArrayDiag.png')
f2,a2 = A.show()
f2.savefig('ArrayConfig.png')
