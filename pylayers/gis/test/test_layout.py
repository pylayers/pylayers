from pylayers.gis.layout import Layout
import networkx as nx
import matplotlib.pyplot as plt

#doctest.testmod(layout)
#L = Layout('TA-Office.ini')
L = Layout()
lL = L.ls()
for tL in lL:
    print tL 
#f = plt.figure(figsize=(20,10))
#plt.axis('off')
#f,a = L.showG('s',nodes=False,fig=f)
#plt.show()
#f,a = L.showG('r',edge_color='b',linewidth=4,fig=f)
#L= Layout('11Dbibli.ini')
#L.show()
#L = Layout('PTIN.ini')
#L = Layout('DLR.ini')
#L.buildGt()
#Ga = L.buildGr()
#L.showGs()
#nx.draw(Ga,Ga.pos)
