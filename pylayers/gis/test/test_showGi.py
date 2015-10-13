from pylayers.gis.layout import Layout
import networkx as nx
import matplotlib.pyplot as plt
import doctest

plt.ion()
#doctest.testmod(layout)
#L = Layout('TA-Office.ini')
L = Layout('WHERE1.ini')
#L= Layout('11Dbibli.ini')
#L.show()
#L = Layout('PTIN.ini')
#L = Layout('DLR.ini')
#L.build()
L.dumpr()
L.showG('i',en=11)
#Ga = L.buildGr()
#L.showGs()
#nx.draw(Ga,Ga.pos)
