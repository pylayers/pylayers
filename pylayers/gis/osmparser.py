from imposm.parser import OSMParser
import pdb

# simple class that handles the parsed OSM data.
class Nodes(object):
    def p1(self,nodes):
        for n in nodes:
            print n 

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

a    = Building()
node = Nodes()
p    = OSMParser(concurrency=4, ways_callback=a.ways)
#p   = OSMParser(concurrency=4, nodes_callback=node.p1)

p.parse('map.osm')

# done
for b in a.building:
    print '-------'
    print a.building[b]['ref']
