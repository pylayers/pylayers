import numpy as np
import matplotlib.pyplot as plt 
import networkx as nx
import doctest


def draw(G,**kwargs):
    """ draw a networkx graph
    G : Graph with pos (geometric graph)


    """

    defaults = {'show': False,
                'fig': [],
                'ax': [],
                'nodes':True,
                'edges':True,
                'labels':True,
                'linewidth': 2,
                'node_color':'w',
                'edge_color':'k',
                'node_size': 200,
                'linewidth': 2,
                'font_size': 30,
                'alphan': 0.8,
                'alphae': 1.0,
                'nodelist': [],
                'edgelist': [],
                'figsize': (8,8)
                 }
    #
    # update default values
    #
    for key, value in defaults.items():
        if key not in kwargs:
            kwargs[key] = value
    #
    # getting fig and ax
    #
    if kwargs['fig'] == []:
        fig = plt.figure(figsize=kwargs['figsize'])
        fig.set_frameon(True)
    else:
        fig = kwargs['fig']

    if kwargs['ax'] == []:
        ax = fig.gca()
    else:
        ax = kwargs['ax']

    #
    # edges list and nodes list
    #

    if kwargs['nodelist']==[]:
        nodelist =  G.nodes()
    if kwargs['edgelist']==[]:
        edgelist =  G.edges()

    if kwargs['nodes']:
        nx.draw_networkx_nodes(G, G.pos,
                               nodelist = nodelist,
                               node_color = kwargs['node_color'],
                               node_size  = kwargs['node_size'],
                               alpha = kwargs['alphan'])
    if kwargs['labels']:
        nx.draw_networkx_labels(G, G.pos,
                                font_size=kwargs['font_size'])

    if kwargs['edges']:
        nx.draw_networkx_edges(G, G.pos,
                               edgelist = edgelist,
                               edge_color = kwargs['edge_color'],
                               linewidth = kwargs['linewidth'],
                               alpha = kwargs['alphae'])
    if kwargs['show']:
        plt.show()

    return fig,ax

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
    plt.ion()
    doctest.testmod()
#    points  = shg.MultiPoint([(0, 0),(0, 1),(2.5,1),(2.5,2),(2.8,2),(2.8,1.1),(3.2, 1.1), (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
#    polyg   = Polygon(points)
#    Gv      = polyg.buildGv() 
#    plt.title('Testing buildGv : wrong design ')
#    #plt.show()
#    A=edgetype(Gv)


