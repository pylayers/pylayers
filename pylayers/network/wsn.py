import networkx as nx
import matplotlib.pylab as plt 

def shortest_path_tree(G,nroot):
    SPT       = nx.Graph()
    SPT.pos   = G.pos
    SPT.node  = G.node
    for n in G.nodes():
        if n!=nroot:
            p = nx.dijkstra_path(G,nroot,n,True)
            SPT.add_path(p)
    return(SPT)

def wsngraph():
    """
    return a weighted graph
    """
    G = nx.Graph()
    G.add_node(1)
    G.add_node(2)
    G.add_node(3)
    G.add_node(4)
    G.add_node(5)
    G.add_node(6)
    G.add_node(7)
    G.add_node(8)
    G.add_node(9)
    G.add_node(10)
    G.add_node(11)
    G.add_node(12)
    G.add_edge(1,3,weight=1)
    G.add_edge(1,2,weight=6)
    G.add_edge(1,12,weight=16)
    G.add_edge(2,11,weight=12)
    G.add_edge(2,6,weight=10)
    G.add_edge(2,5,weight=11)
    G.add_edge(3,4,weight=10)
    G.add_edge(3,7,weight=11)
    G.add_edge(3,8,weight=14)
    G.add_edge(3,9,weight=11)
    G.add_edge(4,7,weight=9)
    G.add_edge(5,6,weight=7)
    G.add_edge(5,9,weight=12)
    G.add_edge(6,9,weight=9)
    G.add_edge(7,10,weight=10)
    G.add_edge(8,10,weight=2)
    G.add_edge(8,11,weight=11)
    G.add_edge(8,9,weight=12)
    G.add_edge(9,11,weight=8)
    G.add_edge(10,12,weight=3)
    G.pos={}
    G.pos[1]=(6,4)
    G.pos[2]=(-1,3.7)
    G.pos[3]=(4.7,3.5)
    G.pos[4]=(5.3,3.2)
    G.pos[5]=(0,3)
    G.pos[6]=(1.4,3.4)
    G.pos[7]=(5,2.6)
    G.pos[8]=(4.7,0)
    G.pos[9]=(1.4,2.4)
    G.pos[10]=(5.2,0.5)
    G.pos[11]=(1.3,0)
    G.pos[12]=(6,2.4)
    elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] > 8]
    esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] <= 8]
    nx.draw_networkx_nodes(G,G.pos,node_color='w')
    nx.draw_networkx_edges(G,G.pos,elarge,width=3,edge_color='r')
    nx.draw_networkx_edges(G,G.pos,esmall,width=1,edge_color='b')
    nx.draw_networkx_labels(G,G.pos)
    return(G) 

#if ( __name__=="__main__"):

G  = wsngraph()
T  = nx.minimum_spanning_tree(G)
nx.draw(T,G.pos)
T.edge.values()

SPT = shortest_path_tree(G,1)




