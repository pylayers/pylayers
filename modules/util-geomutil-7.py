from pylayers.util.geomutil import *
import matplotlib.pyplot as plt
fig = plt.figure()
ax  = fig.gca()
p1 = np.array([0,0])
p2 = np.array([1,0])
p3 = np.array([0,1])
p4 = np.array([1,1])
ax = linet(ax,p1,p2,al=0.7,color='red',linewidth=3)
ax = linet(ax,p2,p3,al=0.8,color='blue',linewidth=2)
ax = linet(ax,p3,p4,al=0.9,color='green',linewidth=1)
ax = linet(ax,p4,p1,al=1,color='cyan',linewidth=0.2)
