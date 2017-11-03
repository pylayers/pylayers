from pylayers.util.geomutil import *
from pylayers.util.plotutil import *
import shapely.geometry as shg
import numpy as np
import scipy as sp
import pdb
from numpy.testing import ( TestCase, assert_almost_equal, assert_raises, assert_equal, assert_, run_module_suite)

print "test_intersect3"
a = np.array([[1,0,1]]).T
b = np.array([[10,0,1]]).T
pg = np.array([[5,0,0]]).T
u1 = np.array([[0,1,0]]).T
u2 = np.array([[0,0,1]]).T
l1 = np.array([3])
l2 = np.array([3])
bo = intersect3(a,b,pg,u1,u2,l1,l2)
print "Occultation is True",bo 
pg = np.array([[5,-4,0]]).T
bo,pinter1 = intersect3(a,b,pg,u1,u2,l1,l2)
print "Occultation is False",bo 
a = np.array([[1,0,1]]).T
b = np.array([[2,-1,1]]).T
u1 = np.array([[np.sqrt(2)/2,-np.sqrt(2)/2,0]]).T
u2 = np.array([[0,0,1]]).T
print "Matrix is Singular - should return False",bo 
bo,pinter2 = intersect3(a,b,pg,u1,u2,l1,l2)
print bo 
