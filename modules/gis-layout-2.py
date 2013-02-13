import matplotlib.pyplot as plt
from pylayers.gis.layout import *
L = Layout()
fillist = L.ls()
for _filename in filelist:
   plt.figure()
   L.load(_filename)
   fig,ax = L.showGs()
   plt.title(_filename)
plt.show()
