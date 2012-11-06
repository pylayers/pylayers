import networkx as nx
import numpy as np 
import math
import random
import doctest
import matplotlib.pylab as plt 
import pdb
from prioritydict import priorityDictionary
#
# This module is adapted from WSN-SIM from Berkeley
#
#   http://wsn.eecs.berkeley.edu/trac/simulation
#
col = ['b','g','r','c','m','y','k','w','b','g','r','c','m','y','k','w']
class Network(nx.Graph):
    def __init__(self, numMotes=100, numGateways=10, loadBalance=10,dx=316,dy=316):
        nx.Graph.__init__(self)
        self.numMotes = numMotes
        self.numGateways = numGateways
        self.loadFactor = loadBalance

        self.packetCount = 0

        #Set up Variables that will be filled in.
        self.schedule = {}
        self.conModel = {}
        self.cxns = {}
        self.gtwys = {}
        self.routes = {}
        self.schedule = {}

        #Packet Simulator Variables
        self.linkPDR = {}
        self.activity = []
        self.queues = []
        self.droppedPkts = []

        self.pos = {}
        for i in range(0, numGateways + numMotes):
            pos = (math.ceil(random.random() * dx), math.ceil(random.random() * dy))
            self.pos[i] = pos
            if i < numGateways:
                self.add_node(i,gateway=True)
            else:
                self.add_node(i,gateway=False)

    def connect(self,dmax=175,pdrValue = 0.80,n=2,sensitivitydBm=-85,fGHz=2.4,PtdBm=0):
        """  connect nodes based on a connectivity model
        Parameters
        ----------

        dmax
        pdrValue :
            Packet Delivery Rate Value  (default 0.8)   (= 1 - PER)

        """
        for i in self.nodes():
            pA = self.pos[i]
            for j in self.nodes():
                if i != j:
                    pB = self.pos[j]
                    dX = pA[0] - pB[0]
                    dY = pA[1] - pB[1]
                    d = np.sqrt(dX**2 + dY**2)
                    d = np.ceil(d + 0.1)

                    if (i in self.cxns) == False:
                        self.cxns[i] = {}

                    if (j in self.cxns) == False:
                         self.cxns[j] = {}

                    if d < dmax:
                        prob,rssi = (np.ceil(connectModel(d,n=n,
                                         sensitivitydBm=sensitivitydBm,
                                         fGHz=fGHz,
                                         PtdBm=PtdBm)))

                        flip = random.random() * 100

                        if flip <= prob:
                            self.add_edge(i,j,rssi=rssi,etx=1./pdrValue)
                            self.cxns[i][j] = 1./pdrValue
                            self.cxns[j][i] = 1./pdrValue

    def disconnect(self):
        self.remove_edges_from(self.edges())


    def routing(self):

        self.D = {}  # dictionary of node -> final Costs
        self.G = {}  # dictionary of node -> Gateway
        self.L = {}  # dictionary of gateway -> Load

        self.Q = priorityDictionary()   # est.dist. of non-final vert.

        self.H = {}  # dictionary of node -> true Hop Counts
        self.P = {}  # dictionary of node -> Predecessor

        for gateway in range(0,self.numGateways):
            self.Q[gateway] = 0
            self.G[gateway] = gateway
            self.L[self.G[gateway]] = 0
            self.H[gateway] = 0

        for v in self.Q:
            self.D[v] = self.Q[v]
            if len(self.D) == self.number_of_nodes():
                break
            #if v in self.cxns:
            neigh_v = self.neighbors(v)
            if len(neigh_v) >0:
                #for w in self.cxns[v]:
                for w in neigh_v:
                    #Check that it is connected
                    if self.cxns[v][w]:
                        gatewayLoad = self.L[self.G[v]]

                        cost = self.D[v] + self.cxns[v][w] # old method to calculate using D table.
                        cost = cost + self.loadFactor * (gatewayLoad) / (200.)
                        if w in self.D:
                            if cost < self.D[w]:
                                raise ValueError("Error Dijkstra: Better path to already-final vertex")
                        elif w not in self.Q or cost < self.Q[w]:
                            self.Q[w] = cost  # assign new cost to Q dictionary
                            self.P[w] = v  # assign predecessor for Parent

                            if self.G.get(w, -1) != -1:
                                self.L[self.G[w]] = self.L[self.G[w]] - 1

                            self.G[w] = self.G[v]  # assign same Gateway as Parent
                            self.L[self.G[v]] = self.L[self.G[v]] + 1  # increment Gateway load
                            self.H[w] = self.H[v] + 1  # assign true Hop Count as increment from Parent
#
#        print "Hops[node]", H
#        print "D: ", D
#        HopsCounter = {}
#
#        for node in H.keys():
#    #               print "Testing node: ", n
#    #               print "Adding 1 to: ", H[n]
#    #       for node in H
#            if H[node] != 0:
#                HopsCounter[H[node]] = HopsCounter.get(H[node], 0) + 1
#    #       print "Added, new total = ", HopsCounter[H[n]]
#
#        return (P, G)

    def scheduling(self, numTimeSlots, numOffsets,export=False):
        """
        Parameters
        ----------
        numTimeSlots
        numOffsets

        """
        paths = []
        pathsByHop = {}
        numSlots = 0
        routes = self.P
        #print "routes"
        #print routes
        #Step 1: For each node in network. Establish a list of links for packet delivery
        for node in routes:
            path = []
            hopNode = node
            while hopNode >= self.numGateways:
                link = (hopNode, routes[hopNode])
                path.append(link)
                numSlots += 1
                hopNode = routes[hopNode]
                #print "link added", link
            pathLen = len(path)
            if pathsByHop.get(pathLen, -1) == -1:
                pathsByHop[pathLen] = []
            pathsByHop[pathLen].append(path)
        #print "PATHS (DICT): ", pathsByHop
        pathKeys = pathsByHop.keys()
        #print "PATH KEYS: ", pathKeys
        for i in range(0, len(pathKeys)):
            pathList = pathsByHop.get(pathKeys[i])
            #print "PATH LIST: ", pathList
            for link in pathList:
                paths.insert(0, link)
        #print "PATH LIST: ", paths
        #print links
        #print "Links to Schedule:", numSlots
        #TextBox.insert(END, "Links to Schedule: %d \n" % numSlots)
        #root.update()
        #interface.printText("Links to Schedule: %d \n", numSlots, 0)
        #global schedule
        self.schedule = {}
        numPass = 0
        maxPass = 3000
        #Step 2: Running through slots, assign schedules	
        pathSelect = 0
        prevPacketTime = -1
        #OPTION TO PACK OR SPREAD:
        #global pack
        #pack1 = pack.get()
        prevPacketTime = -1
        maxLinks = 1
        curLinks = 0
        for t in range(0, numTimeSlots):
            self.schedule[t] = {}
        #2a: Loop through all paths.
        while len(paths) != 0:
            #interface.root.update()
            #print "l:", len(paths)
            #print "pathSelect: ", pathSelect
            path = paths[pathSelect]
            #2b: Loop through all links in the path.
            for link in path:
                linkSet = False
                passes = 0
                slotsChecked = 0
                while (passes < 2):
                    if passes == 1:
                        prevPacketTime = -1
                    #2c: Visit time slots
                    for time in range(prevPacketTime+1, numTimeSlots):
                        #2d: Check if nodes are available during this time slot. Set used = False, change to True if not possible.
                        used = False
                        for key in self.schedule[time].keys():
                            for i in range(0, len(self.schedule[time][key])):
                                curLink = self.schedule[time][key][i]
                                #print "curlink", curLink
                                a = curLink[0]
                                b = curLink[1]
                                if (a == link[0] or b == link[0] or a == link[1] or b == link[1]):
                                    used = True
                                    break
                            if used == True:
                                break
                        #2e: Nodes are available, locate first available offset.
                        if used == False:
                            #Offset positions are available. Fill directly.
                            offset = len(self.schedule[time].keys())
                            if (offset < numOffsets):
                                self.schedule[time][offset] = []
                                linkSet = True
                            #Offset positions are not available. 
                            #Pick minimum slotted and check for interference
                            #Have to check either 0 or 1 for interference/connected links
                            else:
                                offsetCount = []
                                #Loop through all offsets and count current links in the offset.
                                for i in range(0, numOffsets):
                                    count = len(self.schedule[time][i])
                                    offsetCount.append(count)
                                while 1:
                                    offsetC = min(min(offsetCount), 999)
                                    #Not ready to schedule extra link, keep searching...
                                    if offsetC > maxLinks:
                                        if (slotsChecked < numTimeSlots * 15):
                                            slotsChecked += 1
                                            break
                                        else:
                                            slotsChecked = 0
                                            maxLinks += 1
                                    if offsetC == 999:
                                        break
                                    offset = offsetCount.index(offsetC)
                                    #2f For every link already in this time-offset slot, check for interference.
                                    interference = False
                                    for setLink in self.schedule[time][offset]:
                                        for a in range(0, 2):
                                            for b in range(0, 2):
                                                r1 = cxns.get(setLink[a], -1)
                                                #print "setlink/a/setlink[a]", setLink, " ", a," ", setLink[a], " r1: ", r1
                                                if r1 != -1:
                                                    if r1.get(link[b], -1) != -1:
                                                        interference = True
                                    if interference == False:
                                        linkSet = True
                                        break
                                    offsetCount[offset] = 1000 	#Some large number to never try again.
                            if linkSet == True:
                                self.schedule[time][offset].append(link)
                                numScheduled = len(self.schedule[time][offset])
                                #interface.colorCells(time, offset, numScheduled)
                                linkSet = True
                                prevPacketTime = time
                                slotsChecked = 0
                                break
                    passes = 1
                    if linkSet == True:
                        break
            if linkSet == True:
                paths.remove(path)
                pathSelect = 0
            if pathSelect == len(paths):
                #print "Pass ", numPass, " complete."
                numPass += 1
                pathSelect = 0
                prevPacketTime = -1
                if numPass == maxPass:
                    break
            #print "Remaining Links: ", len(links)
        slots = 0
        links = 0
        for time in range(0, 332):
            s = self.schedule.get(time, -1)
            if s == -1:
                continue
            else:
                #print "t: ", time
                for o in range(0, len(s)):
                    l = self.schedule[time][o]
                    #print "offset: ", o, " links: ", len(l)
                    #print len(l)
                    links += len(l)
            slots += len(s)
        maxTime = len(self.schedule)
        avgLpS = links * 1.0 / slots
        if export:
            f = open('schedule.txt', 'w')
            f.write(str(self.numMotes))
            f.write('\n')
            f.write(str(self.numGateways))
            f.write('\n')

            for time in range(0, numTimeSlots):
                first = 1
                timeSchedule = self.schedule.get(time, -1)
                if timeSchedule != -1:
                    for offset in timeSchedule.keys():
                        for link in timeSchedule[offset]:
                            link = str(link)
                            link = link.replace(" ", "")
                            #print "Time: ",time, " Offset: ",offset," Schedule: ", link

                            if first == 0:
                                f.write(';')
                            f.write(link)
                            first = 0

                    f.write('\n')

            f.close()


    def show(self,connectivity=True,routes=False):
        """
        Parameters
        ----------
        conectivity
        routes

        """
        lm = list(np.arange(self.numMotes) + self.numGateways)
        lg = list(np.arange(self.numGateways))
        if connectivity:
            nx.draw_networkx_nodes(self,self.pos,nodelist=lg,node_size=50,node_color='r',alpha=0.5)
            nx.draw_networkx_nodes(self,self.pos,nodelist=lm,node_size=30,node_color='b',alpha=0.2)
            nx.draw_networkx_edges(self,self.pos,node_size=30,alpha=0.1)
            nx.draw_networkx_labels(self,self.pos)
        if routes:
            G = nx.Graph()
            G.add_nodes_from(self)
            for n in self.P:
                G.add_edge(n,self.P[n])
            tG = nx.connected_component_subgraphs(G)
            N  = np.int(np.floor(np.sqrt(len(tG))))
            M  = np.int(np.ceil(len(tG)/(1.*N)))
            for k,g in enumerate(tG):
                #nx.draw(g,self.pos,node_color=col[k])
                plt.subplot(M,N,k+1)
                nx.draw(g,node_color=col[k])
            #for ig in range(self.numGateways):
            #    succcig = self.neighbors(ig)
            #    for 
        #nx.draw_networkx_labels(self,self.pos)
        #plt.show()

#For all potentially interfering nodes (under 175 m) cxns of link = 2. For all connected, this value is 1.


def connectModel(d,fGHz=2.4,sensitivitydBm=-85,n=2,PtdBm=0):
    """ connectivity model

    Parameters
    ----------

    d  : distance in meters
    fGHz : 
    sensitivitydBm : 
    ptdBm : 

    """
    lda = 0.3/fGHz
    PL0 = -20*np.log10(lda/(4*np.pi))
    Pr  = PtdBm - (PL0 + 10*n*np.log10(d))

    rssi = Pr - sensitivitydBm
    if ( rssi > 40):
        cdf = 100

    elif ( rssi  < 0):
        cdf = 0

    else:
        cdf = (rssi / 40.) * 100

    return cdf,rssi


if __name__=="__main__":
    doctest.testmod()
