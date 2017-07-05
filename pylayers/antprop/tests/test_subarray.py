from pylayers.antprop.aarray import *
import matplotlib.pyplot as plt 
import pdb 
print '--------------'
print 'antprop/test_subarray.py'
print '--------------'
fcGHz = 60
lamda = 0.3/fcGHz
N1   = [ 4,4,1]
N2   = [ 2,2,1]
dm1 = [lamda/2.,lamda/2.,0]
dm2 = [3*lamda,3*lamda,0]
A1 = AntArray(fGHz=np.array([fcGHz]),N=N1,dm=dm1,typant='Omni')
A2 = AntArray(fGHz=np.array([fcGHz]),N=N2,dm=dm2,array=A1)
#A1.eval()
