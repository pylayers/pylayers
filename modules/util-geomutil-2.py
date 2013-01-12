from pylayers.util.geomutil import * 
import matplotlib.pyplot as plt
import numpy as np
p1 = np.array([[0,1,1,0],[0,0,1,1]])
P1 = Polygon(p1)
p2 = [[3,4,4,3],[1,1,2,2]]
P2 = Polygon(p2)
p3 = [np.array([10,10]),np.array([11,10]),np.array([11,11]),np.array([10,11])]
P3 = Polygon(p3)
P1.plot('red',alpha=0.3)
P2.plot('blue',alpha=0.7)
P3.plot('green',alpha=1)
title = plt.title('test plotting polygons')
