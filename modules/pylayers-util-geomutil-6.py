from pylayers.util.geomutil import *
import shapely.geometry as shg
import matplotlib.pyplot as plt
points  = shg.MultiPoint([(0, 0), (0, 1), (3.2, 1), (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
polyg1   = Polygon(points)
cvex,ccave   = polyg.ptconvex2() 
points  = shg.MultiPoint([(0, 0), (0, 1), (-3.2, 1), (-3.2, 0.7), (-0.4, 0.7), (-0.4, 0)])
polyg1   = Polygon(points)
cvex,ccave   = polyg.ptconvex2() 
