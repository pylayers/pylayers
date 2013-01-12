import matplotlib.pyplot as plt
from pylayers.gis.layout import *
L = Layout()
for _filename in L.ls():
   plt.figure()
   L.load(_filename)
   fig,ax = L.showGs()
   plt.title(_filename)
plt.show()
