from pylayers.gis.layout import *
import matplotlib.pylab as plt

L = Layout()
#L.load('TA-Office.ini')
L.load('defstr.ini')
L.buildGt()
L.buildGr()
L.buildGv()
L.buildGi()
L.showG('v')
sf1, sf2 = L.signature(0, 1)
plt.show()
