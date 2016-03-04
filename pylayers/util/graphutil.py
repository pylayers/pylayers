#-*- coding:Utf-8 -*-
"""
.. currentmodule:: pylayers.util.graphutil

This class handle the description of an Indoor layout

Utility functions
-----------------

.. autosummary::
    :toctree: generated/

    draw
    edgetype
    find_all_paths

"""
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import doctest
import pdb

def draw(G,**kwargs):
    """ draw a networkx graph


    Parameters
    ----------

    G : Graph with pos (geometric graph)
    show : False
    fig : []
    ax : []
    'nodes':True,
    'edges':True,
    'airwalls':False,
    'labels':True,
    'width': 2,
    'node_color':'w',
    'edge_color':'k',
    'node_size': 200,
    'font_size': 30,
    'alphan': 0.8,
    'alphae': 1.0,
    'nodelist': [],
    'edgelist': [],
    'figsize': (8,8)

    See Also
    --------

    pylayers.gis.layout.showG

    """

    defaults = {'show': False,
                'fig': [],
                'ax': [],
                'arrows':False,
                'nodes':True,
                'edges':True,
                'airwalls':False,
                'labels':True,
                'width': 2,
                'node_color':'w',
                'edge_color':'k',
                'posnode_color':'k',
                'negnode_color':'b',
                'node_size': 200,
                'width': 2,
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
        ax = fig.add_subplot(111,axisbg='white')
    else:
        ax = kwargs['ax']

    #
    # edges list and nodes list
    #

    if kwargs['nodelist']==[]:
        nodelist =  G.nodes()
    else:
        nodelist=kwargs['nodelist']

    if kwargs['edgelist']==[]:
        edgelist =  G.edges()
    else:
        edgelist = map(lambda x : G.edges()[x],kwargs['edgelist']) # for nx an edge list is a list of tuple
    # remove airwalls
    try:
        #pno = filter(lambda x : G.nodes()[x]>0,nodelist)
        pno = filter(lambda x : x>0,nodelist)
        # node == air
        na = filter(lambda x : G.node[x]['name']=='AIR',pno)
        # edge == air
        ea=[]
        [[ea.append((n1,n2)) for n2 in G.edge[n1].keys()] for n1 in na]
        [[ea.append((n2,n1)) for n2 in G.edge[n1].keys()] for n1 in na]

        nodelista = filter(lambda x: x in na, nodelist)
        edgelista = filter(lambda x: x in ea, edgelist)
        nodelist = filter(lambda x: x not in na, nodelist)
        edgelist = filter(lambda x: x not in ea, edgelist)
    except:
        pass
    if kwargs['nodes']:
        nx.draw_networkx_nodes(G, G.pos,
                               nodelist = nodelist,
                               node_color = kwargs['node_color'],
                               node_size  = kwargs['node_size'],
                               alpha = kwargs['alphan'],ax=ax)
    if kwargs['labels']:
        nlp = filter(lambda x: x>0, nodelist)
        nln = filter(lambda x: x<0, nodelist)
        nx.draw_networkx_labels(G, G.pos,
                                labels={n:n for n in nln},
                                font_color=kwargs['negnode_color'],
                                font_size=kwargs['font_size'],ax=ax)
        nx.draw_networkx_labels(G, G.pos,
                                labels={n:n for n in nlp},
                                font_color=kwargs['posnode_color'],
                                font_size=kwargs['font_size'],ax=ax)

    if kwargs['edges']:
        if kwargs['edge_color'] == '':
            kwargs['edge_color']='k'

        nx.draw_networkx_edges(G, G.pos,
                               edgelist = edgelist,
                               edge_color = kwargs['edge_color'],
                               width = kwargs['width'],
                               arrows= kwargs['arrows'],
                               alpha = kwargs['alphae'],ax=ax)
        if kwargs['airwalls']:
            try:
                nx.draw_networkx_edges(G, G.pos,
                                       edgelist = edgelista,
                                       edge_color = kwargs['edge_color'],
                                       width = kwargs['width'],
                                       arrows= kwargs['arrows'],
                                       alpha = kwargs['alphae'],
                                       style='dotted',
                                       ax=ax)
            except:
                pass
    if kwargs['show']:
        plt.show()

    return fig,ax

def edgetype(G):
    """  edge type

    Examples
    --------

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
    """ find all paths

    Parameters
    ----------

    graph :
    start:
    end :

    """
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


