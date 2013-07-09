from pylayers.util.geomutil import * 
import matplotlib.pyplot as plt
import numpy as np
l1 = np.array([[0,1,1,0],[0,0,1,1]])
L1 = LineString(l1)
l2 = [[3,4,4,3],[1,1,2,2]]
L2 = LineString(l2)
fig,ax = L1.plot(color='red',alpha=0.3,linewidth=3)
fig,ax = L2.plot(fig=fig,ax=ax,color='blue',alpha=0.7,linewidth=2)
title = plt.title('test plotting LineString')
