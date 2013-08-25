from pylayers.util.geomutil import *
from pylayers.util.plotutil import *
import shapely.geometry as shg
import numpy as np
import scipy as sp
from numpy.testing import ( TestCase, assert_almost_equal, assert_raises, assert_equal, assert_, run_module_suite)

class Tesonb(TestCase):
    def test_onbfrmaxe(self):
        A = np.array([[0,0,0,0],[1,2,3,4],[0,0,0,0]])
        B = np.array([[0,0,0,0],[1,2,3,4],[10,10,10,10]])
        T = onbfromaxe(A,B)
        print np.shape(T)
        print T[:,0,:]
        print T[:,1,:]
        print T[:,2,:]
        assert_equal(np.shape(T),(4,3,3))

if __name__ == "__main__":
    run_module_suite()
