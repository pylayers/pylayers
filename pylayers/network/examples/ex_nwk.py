import numpy as np
import matplotlib.pylab as plt
import networkx as nx
# Berkeley WSN-SIM module
import pylayers.network.net as net


numMotes = 500
numGateways =  10
loadBalance  = 60
dx = 316
dy = 316
numOffsets = 15
numTimeSlots = 334
# create a network
network = net.Network(numMotes, numGateways, loadBalance)

network.connect(fGHz=2.40,PtdBm=0,n=2)

network.routing()

network.show(connectivity=False,routes=True)
plt.figure()
network.show(connectivity=True,routes=False)
plt.show()
#network.scheduling(numTimeSlots,numOffsets,export=True)

#G = nx.Graph()
#G.add_nodes_from(network)
#for n in network.P:
#    G.add_edge(n,network.P[n])
#
#nx.draw(G,network.pos)
#plt.show()
