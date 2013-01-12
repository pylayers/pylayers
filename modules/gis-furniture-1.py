from pylayers.gis.furniture import *
import matplotlib.pylab as plt 
F = Furniture()
F.load('Furw1.ini','R1_A')
F.show()
axis = plt.axis('scaled')
plt.show() 
