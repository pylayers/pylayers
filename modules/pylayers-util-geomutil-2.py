from pylayers.util.geomutil import *
import shapely.geometry as shg
import matplotlib.pyplot as plt
points  = shg.MultiPoint([(0, 0), (0, 1), (2.5,1), (2.5, 2),                                           (2.8,2), (2.8, 1.1), (3.2, 1.1),                                           (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
polyg   = Polygon(points)
Gv      = polyg.buildGv(show=True)
plt.axis('off')
# (-0.5, 4.0, -0.5, 2.5)
title = plt.title('Testing buildGv')
