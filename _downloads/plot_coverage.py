# -*- coding: utf-8 -*-
r"""
==================================
Indoor Radio Coverage with Motley Keenan
==================================

"""
from pylayers.antprop.coverage import *
import matplotlib.pyplot as plt
import time
C = Coverage()
C.L._filename
C.tx = np.array((39,1))
start = time.time()
C.cover()
finish = time.time()
C.show()
#print 'Tx position: ',C.tx 
#print 'Coverage in %1.2f seconds' % (finish-start)
plt.show()
