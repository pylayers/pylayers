import pylayers.util.geomutil as geu
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from imposm.parser import OSMParser
import networkx as nx
import numpy as np
import pdb

# Simple class that handles the parsed OSM data.
class Way(object):
    def __init__(self,refs,tags,coords):
        self.refs  = refs
        self.tags = tags
        N = len(refs)
        p = np.zeros((2, N))
        for k, nid in enumerate(refs):
            p[0, k] = coords.xy[nid][0]
            p[1, k] = coords.xy[nid][1]
        if (N>=4) & (refs[0]==refs[-1]):
            self.shp = geu.Polygon(p)
        else:
            self.shp = geu.LineString(p)

    def show(self,fig=[],ax=[]):
        """
        """
        fig,ax = self.shp.plot(fig=fig,ax=ax)
        return(fig,ax)

class Coords(object):
    """
    """
    latlon = {}
    cpt = 0
    xy ={}
    minlon = 1000
    maxlon = -1000
    minlat = 1000
    maxlat = -1000

    def coords(self, coords):
        # callback method for coords
        for osmid, lon, lat in coords:
            self.latlon[osmid] = np.array([lon, lat])
            # find extrema
            self.minlon = min(lon,self.minlon)
            self.maxlon = max(lon,self.maxlon)
            self.minlat = min(lat,self.minlat)
            self.maxlat = max(lat,self.maxlat)

            self.cpt += 1
        self.boundary=np.array([self.minlat,self.minlon,self.maxlat,self.maxlon])

    def cartesian(self):
        """
        """
        bd = self.boundary
        lon_0 = (bd[1]+bd[3])/2.
        lat_0 = (bd[0]+bd[2])/2.

        m = Basemap(llcrnrlon=bd[1]-0.01, llcrnrlat=bd[0]-0.01,
                    urcrnrlon=bd[3]+0.01, urcrnrlat=bd[2]+0.01,
                resolution='i', projection='cass', lon_0=lon_0, lat_0=lat_0)

        for id in self.latlon:
            x, y = m(self.latlon[id][0], self.latlon[id][1])
            self.xy[id]  = np.array([x,y])

class Nodes(object):
    """
    """

    node = {}
    cpt = 0

    def nodes(self,nodes):
        """  parse tagged nodes
        """
        for osmid,tags,coords in nodes:
            self.node[osmid] = {}
            self.node[osmid]['tags'] = tags
            self.node[osmid]['lonlat'] = coords
            lon = coords[0]
            lat = coords[1]
            self.cpt += 1

class Ways(object):
    w = {}
    cpt = 0

    def ways(self, ways):
        for osmid, tags, refs in ways:
                self.w[osmid] = (refs,tags)
                self.cpt += 1
    
    def eval(self,coords):
        for osmid in self.w:
            refs = self.w[osmid][0]
            tags = self.w[osmid][1]
            self.way = Way(refs,tags,coords)

    def show(self,fig=[],ax=[]):
        """ show ways
        """
        if fig==[]:
            fig = plt.figure()
        elif ax==[]:
            ax = fig.gca()

        for b in self.way:
            tags  = self.way[b].tags
            if ('building' in tags) | ('buildingpart' in tags):
                try:
                    poly  = self.way[b].poly
                    if 'building:roof:colour' in tags:
                        col = '#'+tags['building:roof:colour']
                    else:
                        col = '#abcdef'
                    fig, ax = poly.plot(fig=fig, ax=ax,color=col)
                except:
                    print "building: ",b," is not a polygon"
        plt.axis('scaled')
        return(fig,ax)

class Relations(object):
    relation = {}
    cpt = 0
    def relations(self,rels):
        for osmid,tags,member in rels:
                self.relation[osmid]={}
                self.relation[osmid]['tags']=tags
                self.relation[osmid]['members']=member
                self.cpt = self.cpt+1

class Building(nx.DiGraph):
    """
    """

    def __init__(self,rootid):
        nx.DiGraph.__init__(self)
        self.rootid=rootid
    def __repr__(self):
        st = str(self.rootid)+'\n'
        levels = nx.neighbors(self,self.rootid)
        st = st + '---'+'\n' 
        for l in levels:
            nw = len(nx.neighbors(self,l))
            st = st + str(l)+' : '+ str(nw) + '\n'
        return(st)



    def build(self,typ,eid):
        """
        """
        if typ=='relation':
            tags = relations.relation[eid]['tags']
            members  =  relations.relation[eid]['members']
        if typ=='way':
            tags = ways.way[eid].tags
            members  =  ways.way[eid].refs
        if typ=='node':
            try:
                tags = nodes.node[eid]['tags']
            except:
                tags = None
            members = None

        self.add_node(eid,tags=tags,type=typ)
        if members is not None:
            for m in members:
                if typ=='relation':
                    eidn = m[0]
                    typn = m[1]
                if typ=='way':
                    eidn = m
                    typn = 'node'

                self.add_edge(eid,eidn)
                self = self.build(typ=typn,eid=eidn)

        return self
    def show(self,nid=None,fig=[],ax=[]):
        """
        """
        if fig==[]:
            fig = plt.figure()
        elif ax==[]:
            ax = fig.gca()

        if nid ==None:
            nid = self.rootid
        nb = nx.neighbors(self,nid)
        for k in nb:
            #print k,self.node[k]['type'],self.node[k]['tags']
            if self.node[k]['type']=='way':
                fig,ax = ways.way[k].show(fig=fig,ax=ax)
            else:
                fig,ax = self.show(k,fig=fig,ax=ax)
        return fig,ax

def osmparse(filename):
    """
    """
    
    coords = Coords()
    nodes = Nodes()
    ways = Ways()
    relations = Relations()

    coords_parser = OSMParser(concurrency=4, coords_callback=coords.coords)
    nodes_parser = OSMParser(concurrency=4, nodes_callback=nodes.nodes)
    ways_parser = OSMParser(concurrency=4, ways_callback=ways.ways)
    relations_parser = OSMParser(concurrency=4,relations_callback=relations.relations)

    print "parsing coords"
    coords_parser.parse(filename)
    coords.cartesian()
    print str(coords.cpt)

    print "parsing nodes"
    nodes_parser.parse(filename)
    print str(nodes.cpt)

    print "parsing ways"
    ways_parser.parse(filename)
    print str(ways.cpt)
    ways.eval(coords)

    print "parsing relations"
    relations_parser.parse(filename)
    print str(relations.cpt)

    return coords,nodes,ways,relations

def buildingsparse(filename):
    """
    """
    coords,nodes,ways,relations = osmparse(filename)
    for bid in relations.relation:
        tags = relations.relation[bid]['tags']
        if tags['type']=='building':
            print "Constructing Indoor building ", bid
            bdg = Building(bid)
            bdg = bdg.build(typ='relation',eid=bid)
    return bdg

