from pylayers.util.geomutil import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from imposm.parser import OSMParser
import numpy as np 
import pdb

# simple class that handles the parsed OSM data.
#class Nodes(object):
#    nodes = {}
#    def p1(self,nodes):
#        for n in nodes:
#             print n[0],n[2]

# Simple class that handles the parsed OSM data. 
class Nodes(object): 
    coordsNum = 0 
    latlon ={}
    def coords(self, coords): 
        # callback method for coords 
        for osm_id, lon, lat in coords: 
            self.latlon[osm_id]= np.array([lon,lat])
            print osm_id,lon,lat
            self.coordsNum += 1 

class Building(object):
    building = {} 
    cpt = 0
    def ways(self, ways):
        # callback method for ways
        for osmid, tags, refs in ways:
            if 'building' in tags:
                self.building[osmid]={}
                self.building[osmid]['ref']=refs
                self.building[osmid]['tags']=tags
                self.cpt += 1

# instantiate counter and parser and start parsing

# Instantiate counter and  parser and start parsing 
#counter = CoordsCounter() 
#p = OSMParser(concurrency=4, coords_callback=counter.coords) 
a    = Building()
node = Nodes()
p    = OSMParser(concurrency=4, ways_callback=a.ways)
q    = OSMParser(concurrency=4, coords_callback=node.coords)

m = Basemap(llcrnrlon=-1.65263,llcrnrlat=48.1127,urcrnrlon=-1.62759,urcrnrlat=48.12547,
            resolution='i',projection='cass',lon_0=-1.63,lat_0=48.115)
q.parse('map.osm')
p.parse('map.osm')

# done

for b in a.building:
    print '-------'
    f =  a.building[b]['ref']
    N = len(f)
    print N
    p = np.zeros((2,N))
    for k,id in enumerate(f):
        x,y = m(node.latlon[id][0],node.latlon[id][1])
        p[0,k] = x
        p[1,k] = y
    #print p        
    Pol=Polygon(p)
    Pol.plot()
    del p
plt.show()    
