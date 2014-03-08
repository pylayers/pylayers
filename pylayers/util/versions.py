import sys
import IPython
import numpy
import scipy
import matplotlib
import networkx
import shapely
import doctest
import sphere
import SimPy
import pandas


b = eval(IPython.sys_info())
for k in b:
    print k, b[k]
print 'Python   : '   + str(sys.version)
print 'numpy    : '   + str(numpy.__version__)
print 'scipy    : '   + str(scipy.__version__)
print 'matplotlib : ' + str(matplotlib.__version__)
print 'networkx : '   + str(networkx.__version__)
print 'SImPy : '   + str(SimPy.__version__)
print 'Pandas : '   + str(pandas.__version__)
