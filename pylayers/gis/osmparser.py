from pylayers.util.geomutil import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from imposm.parser import OSMParser
import numpy as np
import pdb

# simple class that handles the parsed OSM data.
# class Nodes(object):
#    nodes = {}
#    def p1(self,nodes):
#        for n in nodes:
#             print n[0],n[2]

# Simple class that handles the parsed OSM data.


class Nodes(object):
    """
    """
    coordsNum = 0
    latlon = {}
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
            #if 'door' in tags:
            #    self.typ[osmid] = "door"
            #    self.height[osmid] = "door"
            #    self.width[osmid] = "door"

            self.coordsNum += 1
        self.boundary=np.array([self.minlat,self.minlon,self.maxlat,self.maxlon])

class Relation(object):

    def relation(self,rels):
        for osmid,tags,member in rels:
            print "osmid",osmid
            print "tags",tags
            print "member",member

class Buildings(object):
     building = {}
    cpt = 0

    def ways(self, ways):
        # callback method for ways
        for osmid, tags, refs in ways:
            if ('building' in tags) or ('buildingpart' in tags):
                self.building[osmid] = {}
                self.building[osmid]['ref'] = refs
                self.building[osmid]['tags'] = tags
                self.cpt += 1

    def cartesian(self,nodes):
        """
        """
        for bid in self.building:
            bref  = self.building[bid]['ref']
            N = len(bref)
            p = np.zeros((2, N))
            for k, id in enumerate(bref):
                x, y = m(nodes.latlon[id][0], nodes.latlon[id][1])
                p[0, k] = x
                p[1, k] = y
            self.building[bid]['poly'] = Polygon(p)


    def show(self,fig=[],ax=[]):
        """
        """
        if fig==[]:
            fig = plt.figure()
        elif ax==[]:
            ax = fig.gca()

        for b in self.building:
            tags  = self.building[b]['tags']
            poly  = self.building[b]['poly']
            if 'building:roof:colour' in tags:
                col = '#'+tags['building:roof:colour']
            else:
                col = '#abcdef'
            fig, ax = poly.plot(fig=fig, ax=ax,color=col)
        plt.axis('scaled')
        return(fig,ax)

#
# instantiate counter and parser and start parsing

# Instantiate counter and  parser and start parsing
# counter = CoordsCounter()
# p = OSMParser(concurrency=4, coords_callback=counter.coords)
#
#

bats = Buildings()
nodes = Nodes()
rels = Relations()

building_parser = OSMParser(concurrency=4, ways_callback=bats.ways)
node_parser = OSMParser(concurrency=4, coords_callback=nodes.coords)
relation_parser = OSMParser(concurrency=4, member_callback=rels.relation)

#m = Basemap(
#    llcrnrlon=-1.65263, llcrnrlat=48.1127, urcrnrlon=-1.62759, urcrnrlat=48.12547,
#    resolution='i', projection='cass', lon_0=-1.63, lat_0=48.115)

#_filename = 'map.osm'
_filename = 'indoor.osm'
# q.parse('map.osm')
# p.parse('map.osm')

#
# 1 - Parse OSM nodes
# 2 - Determine lat-lon bounding box
#
node_parser.parse(_filename)

bd = nodes.boundary
lon_0 = (bd[1]+bd[3])/2.
lat_0 = (bd[0]+bd[2])/2.

m = Basemap(llcrnrlon=bd[1], llcrnrlat=bd[0], urcrnrlon=bd[3], urcrnrlat=bd[2],
    resolution='i', projection='cass', lon_0=lon_0, lat_0=lat_0)


building_parser.parse(_filename)
bats.cartesian(nodes)


relation_parser()

fig = plt.figure()
ax = fig.gca()
#fig,ax = bats.show(fig=fig,ax=ax)
#plt.show()

indoor ={}
nbat = 0
#
# Get the shell of the building
#
# Ajouter les nodes id impliques
#
for b in bats.building:
    if 'buildingpart' in bats.building[b]['tags']:
        if  bats.building[b]['tags']['buildingpart']=='shell':
            pshell = bats.building[b]['poly']
            indoor[nbat]={'id':b,'shell':pshell}
            nbat +=1
            fig,ax =pshell.plot(fig=fig,ax=ax)

#
# Get the room included within the shell
#
for bid in indoor:
    indoor[bid]['level']={}
    pshell = indoor[bid]['shell']
    for b in bats.building:
        tags = bats.building[b]['tags']
        if b != indoor[bid]['id']:
            if 'buildingpart' in tags:
                print tags
                try :
                    level = tags['level']
                except:
                    level = 0
                if (tags['buildingpart']=='room') | \
                   (tags['buildingpart']=='corridor') | \
                   (tags['buildingpart']=='hall') | \
                   (tags['buildingpart']=='verticalpassage'):
                    proom  = bats.building[b]['poly']
                    if proom.within(pshell):
                        try:
                            indoor[bid]['level'][level].append(proom)
                        except:
                            indoor[bid]['level'][level]=[proom]

for bid in indoor:
    for level in indoor[bid]['level']:
        for r in indoor[bid]['level'][level]:
            fig,ax = r.plot(fig=fig,ax=ax)
