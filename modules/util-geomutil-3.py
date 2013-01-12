from pylayers.util.geomutil import *
import shapely.geometry as shg
import matplotlib.pyplot as plt
points  = shg.MultiPoint([(0, 0), (0, 1), (3.2, 1), (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
N = len(points)
polyg   = Polygon(points)
tcc,n   = polyg.ptconvex()
k = 0
for p in points:
  if tcc[k] == 1 :
      plt.plot(p.x, p.y, 'o', color='red',alpha=1)
  else:
      plt.plot(p.x, p.y, 'o', color='blue',alpha=0.3)
  k = k+1
polyg.plot()
plt.figure()
points  = shg.MultiPoint([(0, 0), (1, 1), (2, 0), (1, 0)])
poly    = Polygon(points)
tcc,n   = polyg.ptconvex()
poly.plot()
