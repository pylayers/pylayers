import numpy as np
import matplotlib.pyplot as plt 
import networkx as nx
import doctest

def edgetype(G):
    """

    .. plot::
        :include-source: 

        >>> from pylayers.util.geomutil import *
        >>> import shapely.geometry as shg 
        >>> import matplotlib.pyplot as plt 
        >>> points = shg.MultiPoint([(0, 0),(0, 1),(1,1),(1.5,1),(2.5,1),(2.5,2),(2.8,2),(2.8,1.1),(3.2, 1.1), (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
        >>> polyg  = Polygon(points)
        >>> Gv     = polyg.buildGv(show=True) 
        >>> plt.show()

    """
    edges = np.array(G.edges())
    tedg  = np.array(G.edges())
    eprod = np.prod(tedg,axis=1)
    esum  = np.sum(tedg,axis=1)
    inded = np.nonzero(eprod<0)[0] 
    ekeep = np.nonzero(eprod>0)[0]
    ieded = np.nonzero(esum>0)[0]
    indnd = np.nonzero(esum<0)[0]
    u1 = np.in1d(ieded,ekeep)
    u2 = np.in1d(indnd,ekeep)

    nded = list(edges[inded])
    eded = list(edges[ieded[u1]])
    ndnd = list(edges[indnd[u2]])
    return(ndnd,nded,eded)


def find_all_paths(graph, start, end):
    path  = []
    paths = []
    queue = [(start, end, path)]
    while queue:
        start, end, path = queue.pop()
        #print 'PATH', path

        path = path + [start]
        if start == end:
            paths.append(path)
        for node in set(graph[start]).difference(path):
            queue.append((node, end, path))
    return paths


if __name__=="__main__":

    doctest.testmod()
#    points  = shg.MultiPoint([(0, 0),(0, 1),(2.5,1),(2.5,2),(2.8,2),(2.8,1.1),(3.2, 1.1), (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
#    polyg   = Polygon(points)
#    Gv      = polyg.buildGv() 
#    plt.title('Testing buildGv : wrong design ')
#    #plt.show()
#    A=edgetype(Gv)


