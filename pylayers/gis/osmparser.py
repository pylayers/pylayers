#-*- coding:Utf-8 -*-
"""

Way Class
==========

.. autosummary::
    :toctree: generated/

     Way.__init__
     Way.__repr__
     Way.show

Coords Class
============

.. autosummary::
    :toctree: generated/

    Coords.__repr__
    Coords.clean
    Coords.coords
    Coords.cartesian

Nodes Class
============

.. autosummary::
    :toctree: generated/

    Nodes.nodes
    Nodes.clean

Ways Class
============

.. autosummary::
    :toctree: generated/

    Ways.ways
    Ways.clean
    Ways.building
    Ways.eval
    Ways.show
    Ways.tomaska
    Ways.showold

Relations Class
===============

.. autosummary::
    :toctree: generated/

     Relations.relations
     Relations.clean

FloorPlan Class
===============

.. autosummary::
    :toctree: generated/

     FloorPlan.__init__
     FloorPlan.__repr__
     FloorPlan.build
     FloorPlan.show

Utility Functions
=================

.. autosummary::
    :toctree: generated/

     osmparse
     extract
     getbdg
     buildingsparse

"""
#
# Module OSMParser
#
# This module provides classes to handle open stree map objects
#
#


from osmapi import OsmApi
import urllib
import urllib2 as url
from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
# imposm is required for handling osm files
# the installation of imposm is not straightforward
try:
    from imposm.parser import OSMParser
except:
    print "Warning : OSM Parser seems to be not installed"
import networkx as nx
import numpy as np
import pdb

# classes that handle the OSM data file format.
class Way(object):
    """

    A Way is a polyline or a Polygon (if closed)

    typ : 0 Polygon
          1 LineString

    """
    def __init__(self,refs,tags,coords):
        """ object constructor

        Parameters
        ----------

        refs  :
        tags  :
        coords :

        """
        self.refs  = refs
        self.tags = tags
        N = len(refs)
        p = np.zeros((2, N))
        self.valid = True
        for k, nid in enumerate(refs):
            try:
                p[0, k] = coords.xy[nid][0]
                p[1, k] = coords.xy[nid][1]
            except:
                self.valid=False
                break
        # closed way or open way
        if self.valid:
            if (N>=4) & (refs[0]==refs[-1]):
                self.shp = geu.Polygon(p)
                self.typ = 0
            else:
                self.shp = geu.LineString(p)
                self.typ = 1

    def __repr__(self):
        st = ''
        st = st + str(self.tags) + ':' + str(self.refs)
        return(st)

    def show(self,fig=[],ax=[]):
        """ show way

        Parameters
        ----------

        fig : matplotlib figure
        ax  : axes

        """
        fig,ax = self.shp.plot(fig=fig,ax=ax)
        return(fig,ax)

class Coords(object):
    """
    Coords is a point in OSM

    Attributes
    ----------

    xy :
    latlon :
    cpt :
    minlon :
    maxlon :
    minlat :
    maxlat :
    boundary : np.array
        (minlon,minlat,maxlon,maxlat) 

    Notes
    -----


    """
    cpt = 0
    latlon = {}
    xy = {}
    minlon = 1000
    maxlon = -1000
    minlat = 1000
    maxlat = -1000

    def __repr__(self):
        st = ''
        for k in self.xy:
            st = st + str(k)+ ':' + str(self.xy[k])+'\n'
        st = st+ 'Ncoords = '+ str(len(self.xy))+'\n'
        return(st)

    def clean(self):
        self.cpt = 0
        self.latlon={}
        self.xy = {}
        self.minlon =1000
        self.maxlon =-1000
        self.minlat =1000
        self.maxlat =-1000

    def coords(self, coords):
        """ calculates extrema of coords

        Parameters
        ----------

        coords

        """
        # callback method for coords
        for osmid, lon, lat in coords:
            self.latlon[osmid] = np.array([lon, lat])
            # find extrema
            self.minlon = min(lon,self.minlon)
            self.maxlon = max(lon,self.maxlon)
            self.minlat = min(lat,self.minlat)
            self.maxlat = max(lat,self.maxlat)

            self.cpt += 1
        self.boundary=np.array([self.minlon,self.minlat,self.maxlon,self.maxlat])

    def cartesian(self,cart=False,delta=0):
        """ convert Latitude/Longitude in cartesian

        Parameters
        ----------
        cart : Boolean 
            conversion to cartesian 
        delta : offset 
 +            default 0 : in this case the origin corresponds to the lower left point
 
        Notes
        -----

        This method converts latlon coordinates into cartesian x,y coordinates in
        Cassini projection relatively to specified latlon boundary 
        The basemap objet for back and forth coordinates.
        conversion is returned.

        The transformation is centered on the mean of latitude and longitude. 
        The cartesian origin (0,0) correspond to the lower left corner (lonmin,latmin) 

        Returns
        -------

        m : Basemap converter



        Warning
        -------

        If boundary is modified cartesian coordinates change.


        """
        bd = self.boundary
        lon_0 = (bd[0]+bd[2])/2.
        lat_0 = (bd[1]+bd[3])/2.

        m = Basemap(llcrnrlon=bd[0]-delta, llcrnrlat=bd[1]-delta,
                    urcrnrlon=bd[2]+delta, urcrnrlat=bd[3]+delta,
                resolution='i', projection='cass', lon_0=lon_0, lat_0=lat_0)

        for id in self.latlon:
            if cart:
                x, y = m(self.latlon[id][0], self.latlon[id][1])
            else:
                x, y = (self.latlon[id][0], self.latlon[id][1])

            self.xy[id]  = np.array([x,y])

        return(m)

class Nodes(object):
    """

    osm Nodes container

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

    def clean(self):
        self.node= {}
        self.cpt = 0

class Ways(object):
    """

    Attributes
    ----------

    w : dict
    way : dict
    cpt : int

    Methods
    -------

    ways
    eval
    show

    """
    w = {}
    way = {}
    cpt = 0

    def ways(self, ways):
        """ general callback function
        """
        for osmid, tags, refs in ways:
            self.w[osmid] = [refs,tags]
            self.cpt += 1

    def clean(self):
        self.w = {}
        self.way = {}
        self.cpt = 0

    def building(self, ways):
        """ building callback function
        """
        for osmid, tags, refs in ways:
            if 'building' in tags:
                if 'height' in tags:
                    tags = {'height':tags['height']}
                else:
                    tags = {'height':8.5}

                self.w[osmid] = [refs,tags]
                self.cpt += 1

    def eval(self,coords):
        """ convert into a Way object

        Parameters
        ----------

        coords : osm coordinates

        """

        for osmid in self.w:
            refs = self.w[osmid][0]
            tags = self.w[osmid][1]
            away =  Way(refs,tags,coords)
            if away.valid:
                self.way[osmid] = away

    def show(self,typ=2,**kwargs):
        """ show all way

        Parameters
        ----------

        fig : figure
        ax  : axe
        typ : 0|1|2 (default)
                0 : display only way of typ 0  (Polygon)
                1 : display only way of typ 1  (Linestring)
                2 : display all way (default)
        """
        if 'fig' not in kwargs:
            fig = plt.figure(**kwargs)
        else:
            fig = kwargs['fig']

        if 'ax' not in kwargs:
            ax = fig.gca()
        else:
            ax = kwargs['ax']

        lonmin = 360
        lonmax = -360
        latmin = 360
        latmax = -360
        for b in self.way:
            if typ==0:
                if self.way.typ==0:
                    p =ways.way[b].shp
            if typ==1:
                if self.way.typ==1:
                    p = self.way[b].shp
            if typ==2:
                p = self.way[b].shp

            lpoly.append(p)
        city = PolyCollection(lpoly,closed=False)

        ax.axis((lonmin,lonmax,latmin,latmax))
        ax.add_collection(city)
        ax.autoscale_view()
        plt.axis('scaled')
        return(fig,ax)

    def tomaska(self,m,lonlat=True):
        """ convert to masked array

        Parameters
        ----------

        m : Basemap object
            for converting to and from map projection coordinates
        lonlat : boolean
            returns in WGS84 format if True

        Returns
        -------

        ptma : masked array

        """

        tpt  = np.empty((2,))
        mask = np.ones((2,))
        N = len(self.way.keys())
        for k,b in enumerate(self.way):
            # retrieve PolyGon or LineString
            # progress bar
            if k%1000==0:
                print k,N
            shp = self.way[b].shp
            if type(shp)==geu.Polygon:
                pa = self.way[b].shp.ndarray()
                Np = np.shape(pa)[1]
                for ip in range(Np+1):
                    tpt  = np.vstack((tpt,pa[:,ip%Np]))
                    mask = np.vstack((mask,np.array([[0,0]])))
                tpt = np.vstack((tpt,np.array([[0,0]])))
                mask = np.vstack((mask,np.array([[1,1]])))

        if lonlat:
            (lon,lat) = m(tpt[:,0],tpt[:,1],inverse=True)
            tpt = np.vstack([lon,lat]).T

        #vertices = np.ma.masked_array(tpt, mask)
        #return(vertices)
        return(tpt,mask)

    def showold(self,fig=[],ax=[]):
        """ show ways

        Parameters
        ----------

        fig
        ax

        """
        if fig==[]:
            fig = plt.figure()
        elif ax==[]:
            ax = fig.gca()

        for b in self.way:
            tags  = self.way[b].tags
            if ('building' in tags) | ('buildingpart' in tags):
                print "buildingpart found"
                try:
                    poly  = self.way[b].shp
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
    def clean(self):
        self.relation= {}
        self.cpt = 0

class FloorPlan(nx.DiGraph):
    """ FloorPlan class derived from nx.DigGraph

    """

    def __init__(self,rootid,coords,nodes,ways,relations):
        """ object constructor

        Parameters
        ----------

        rootid
        coords
        nodes
        ways
        relations

        """
        nx.DiGraph.__init__(self)
        self.rootid=rootid
        self.coords = coords
        self.nodes = nodes
        self.ways = ways
        self.relations = relations

    def __repr__(self):

        st = str(self.rootid)+'\n'
        levels = nx.neighbors(self,self.rootid)
        st = st + '---'+'\n'
        for l in levels:
            nw = len(nx.neighbors(self,l))
            st = st + str(l)+' : '+ str(nw) + '\n'
        return(st)



    def build(self,typ,eid):
        """ Notes : recursive construction

        Parameters
        ----------

        typ : string
            'relation' | 'way' | 'node'

        """
        if typ=='relation':
            tags = self.relations.relation[eid]['tags']
            members  =  self.relations.relation[eid]['members']
        if typ=='way':
            tags = self.ways.way[eid].tags
            members  =  self.ways.way[eid].refs
        if typ=='node':
            try:
                tags = self.nodes.node[eid]['tags']
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
                #self = self.build(typ=typn,eid=eidn)
                self.build(typ=typn,eid=eidn)


    def show(self,nid=None,fig=[],ax=[]):
        """ show the floorplan

        Parameters
        ----------
        nid :
        fig :
        ax  :

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
                fig,ax = self.ways.way[k].show(fig=fig,ax=ax)
            else:
                fig,ax = self.show(k,fig=fig,ax=ax)
        return fig,ax
#
#  Functions
#     osmparse
#     getbdg
#
#
def osmparse(_filename,typ='floorplan',verbose=False,c=True,n=True,w=True,r=True,cart=False):
    """ parse osm files

    Parameters
    ----------

    typ : string
        floorplan | building
    verbose : boolean
        default : False
    c : boolean
        read coords
    n : boolean
        read nodes
    w : boolean
        read  ways
    r : boolean
        read relations

    Returns
    -------

    coords :
    nodes  :
    ways   :
    relations :

    """

    if (('/' in _filename) or ('//' in _filename)):
        filename = _filename
    else:
        filename = pyu.getlong(_filename,os.path.join('gis','osm'))
    #
    # Read coords and create basemap converter
    #
    if c:
        coords = Coords()
        coords.clean()
        coords_parser = OSMParser(concurrency=4, coords_callback=coords.coords)
        if verbose:
            print "parsing coords"

        coords_parser.parse(filename)

        m = coords.cartesian(cart=cart)

        if verbose:
            print str(coords.cpt)
    else:
        coords = None
    #
    # Read nodes
    #
    if n:
        nodes = Nodes()
        nodes.clean()
        nodes_parser = OSMParser(concurrency=4, nodes_callback=nodes.nodes)

        if verbose:
            print "parsing nodes"

        nodes_parser.parse(filename)

        if verbose:
            print str(nodes.cpt)
    else:
        nodes = None

    #
    # Read ways
    #
    if w:
        ways = Ways()
        ways.clean()
        if typ=='building':
            ways_parser = OSMParser(concurrency=4, ways_callback=ways.building)
        if typ=='floorplan':
            ways_parser = OSMParser(concurrency=4, ways_callback=ways.ways)
        if verbose:
            print "parsing ways"
        ways_parser.parse(filename)

        if verbose:
            print str(ways.cpt)

        # convert lat,lon in cartesian
        ways.eval(coords)
    else:
        ways = None

    #
    # Read realtions
    #
    if r:
        relations = Relations()
        relations.clean()
        relations_parser = OSMParser(concurrency=4,relations_callback=relations.relations)

        if verbose:
            print "parsing relations"
        relations_parser.parse(filename)

        if verbose:
            print str(relations.cpt)
    else:
        relations = None

    return coords,nodes,ways,relations,m


def extract(alat,alon,fileosm,fileout):
    """ extraction of an osm sub region using osmconvert

    Parameters
    ----------

    alat : array of latitude (1xn)
    alon : array of longitude (1xn)
    fileosm : source osm file
    filout  : output osm file

    Returns
    -------

    m : Basemap oject for coordinates conversion



    Notes
    -----

    This function takes two (1xn) arrays of latitude an longitude values
    Calculates extrema of those values.
    Invokes osmconvert script on a source fileosm and extract the
    corresponding zone in the fileout file.

    The functions returns a basemap object for coordinates conversion on this
    file.



    """
    latmax = alat.max()
    latmin = alat.min()
    lonmax = alon.max()
    lonmin = alon.min()
    lon_0=(lonmax+lonmin)/2
    lat_0=(latmax+latmin)/2

    command = 'osmconvert -b='+str(lonmin)+','\
            + str(latmin)+','+str(lonmax)+','\
            + str(latmax)+' '+fileosm +' > '+ fileout+'.osm'
    print command
    os.system(command)

    m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,
            resolution='i',projection='cass',lon_0=lon_0,lat_0=lat_0)

    return(m)

def getbdg(fileosm,verbose=False):
    """ get building from osm file

    Parameters
    ----------

    fileosm : string

    Returns
    -------

    zone : list of Polygon

    """

    coords,nodes,ways,relation,m = osmparse(fileosm,typ='building',verbose=verbose)
    zone = []
    for w in ways.way:
        zone.append(ways.way[w].shp)
    return(zone)

def buildingsparse(filename):
    """ parse buildings

    Parameters
    ----------

    filename : string

    """
    coords,nodes,ways,relations,m = osmparse(filename)
    for bid in relations.relation:
        tags = relations.relation[bid]['tags']
        if tags['type']=='building':
            print "Constructing Indoor building ", bid
            bdg = FloorPlan(bid,coords,nodes,ways,relations)
            bdg.build(typ='relation',eid=bid)
    return bdg,m

