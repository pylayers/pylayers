# -*- moding: utf-8 -*-
"""
Module OSMParser

.. currentmodule:: pylayers.gis.osmparser

.. autosummary::

"""
from osmapi import OsmApi
import geocoder as geo
import sys
import urllib

if sys.version_info.major == 2:
    from  urllib2 import urlopen
else:
    from  urllib.request import urlopen

from pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import xml.etree.ElementTree as xml

import pdb

# classes that handle the OSM data file format.
class Way(object):
    """

    A Way is a polyline or a Polygon (if closed)

    typ : 0 Polygon
          1 LineString

    """
    def __init__(self, refs, tags, coords):
        """ object constructor

        Parameters
        ----------

        refs  :
        tags  :
        coords :
        nodes_sign : int
            if data comes from osm nodes are >0 in ways sequence
            if data comes from josm editor nodes are <0 ways sequence

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
        fig, ax = self.shp.plot(fig=fig, ax=ax)
        return fig, ax

class Coords(object):
    """ Coords describes a set of points in OSM

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

    def __init__(self,idx=[],latlon=[]):
        """
        """
        if latlon!=[]:
            for k,ix in enumerate(idx):
                self.latlon[k]=np.array(latlon[k][0],latlon[k][1])

    def __repr__(self):
        st = ''
        for k in self.xy:
            st = st + str(k)+ ':' + str(self.xy[k])+'\n'
        st = st+ 'Ncoords = '+ str(len(self.xy))+'\n'

        return(st)

    def filter(self,lexcluded):
        for ix in lexcluded:
            self.latlon.pop(-np.abs(ix))
            self.xy.pop(-np.abs(ix))
            self.cpt-=1

    def clean(self):
        """ reset coordinates

        """
        self.cpt = 0
        self.latlon={}
        self.xy = {}
        self.minlon = 1000
        self.maxlon = -1000
        self.minlat = 1000
        self.maxlat = -1000

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

        self.boundary = np.array([self.minlon,self.minlat,self.maxlon,self.maxlat])

    def cartesian(self, cart=False, delta=0, projection='cass'):
        """ convert Latitude/Longitude in cartesian

        Parameters
        ----------

        cart : Boolean
            conversion to cartesian
        delta : offset
            default 0 : in this case the origin corresponds to the lower left point
        projection = ['aeqd','gnom','ortho','geos','nsper','moll','lcc','laea','cass']

        Notes
        -----

        This method converts latlon coordinates into cartesian x,y coordinates in
        a given projection (Default Cassini) relatively to specified latlon boundary
        The basemap objet for back and forth coordinates conversion is returned.

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
        print(bd)
        lon_0 = (bd[0]+bd[2])/2.
        lat_0 = (bd[1]+bd[3])/2.

        m = Basemap(llcrnrlon=bd[0]-delta, llcrnrlat=bd[1]-delta,
                    urcrnrlon=bd[2]+delta, urcrnrlat=bd[3]+delta,
                resolution='i', projection = projection, lon_0=lon_0, lat_0=lat_0)

        for kid in self.latlon:
            if cart:
                x, y = m(self.latlon[kid][0], self.latlon[kid][1])
            else:
                x, y = (self.latlon[kid][0], self.latlon[kid][1])

            self.xy[kid]  = np.array([x,y])

        return(m)

    def from_nodes(self,nodes):
        """ read coordinates from nodes

        Parameters
        ----------

        nodes : Nodes


        """
        for osmid in nodes.node:
            lon = nodes.node[osmid]['lonlat'][0]
            lat = nodes.node[osmid]['lonlat'][1]
            self.latlon[osmid] = np.array([lon,lat])
            # find extrema
            self.minlon = min(lon,self.minlon)
            self.maxlon = max(lon,self.maxlon)
            self.minlat = min(lat,self.minlat)
            self.maxlat = max(lat,self.maxlat)

            self.cpt += 1

        self.boundary = np.array([self.minlon,
                                  self.minlat,
                                  self.maxlon,
                                  self.maxlat])



class Nodes(object):
    """

    osm Nodes container

    """

    node = {}
    cpt = 0

    def nodes(self,nodes):
        """  parse tagged nodes
        """
        for osmid, tags, coords in nodes:
            self.node[osmid] = {}
            self.node[osmid]['tags'] = tags
            self.node[osmid]['lonlat'] = coords
            lon = coords[0]
            lat = coords[1]
            self.cpt += 1

    def __repr__(self):
        st = ''
        for kid in self.node:
            st = st+str(kid) + ' : ' + str(self.node[kid]['lonlat'])+'\n'
        return(st)

    def clean(self):
        self.node= {}
        self.cpt = 0

    def readmap(self,osmmap):
        """ read nodes from a map

        """
        for item in osmmap:
            if item['type']=='node':
                osmid = -item['data']['id']
                if osmid>0:
                    osmid = -osmid
                lon = item['data']['lon']
                lat = item['data']['lat']
                self.node[osmid]={}
                self.node[osmid]['tags'] = item['data']['tag']
                self.node[osmid]['lonlat'] = (lon,lat)
                self.cpt += 1

class Ways(object):
    """

    Attributes
    ----------

    w : dict
    way : dict
    cpt : int
        way counter


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

    def __repr__(self):
        st = ''
        for kid in self.w:
            st = st + str(self.w[kid])+'\n'
        return st

    def clean(self):
        """ clean ways

        """
        self.w = {}
        self.way = {}
        self.cpt = 0

    def building(self, ways, height=8.5, min_heigh=0):
        """ building callback function

        Parameters
        ----------

        ways : Ways
        height : float

        """
        for osmid, tags, refs in ways:
            if 'building' in tags:
                ntags ={}
                # height : from ground to roof top
                if 'height' in tags:
                    ntags['height'] = tags['height']
                elif 'building:height' in tags:
                    ntags['height']  = tags['building:height']
                else:
                    ntags['height'] = height

                # min_height : from ground to roof top
                if 'building:min_height' in tags:
                    ntags['min_height'] = tags['building:min_height']
                elif 'min_height' in tags:
                    ntags['min_height'] = tags['min_height']
                else:
                    ntags['min_height'] = min_height

                # material : from ground to roof top
                if 'building:material' in tags:
                    ntags['material'] = tags['building:material']
                elif 'material' in tags:
                    ntags['material'] = tags['material']
                else:
                    ntags['material'] = 'WALL'

                self.w[osmid] = [refs, ntags]
                self.cpt += 1

    def toway(self,coords):
        """ convert into a Way object

        Parameters
        ----------

        coords : osm coordinates

        """
        for osmid in self.w:
            refs = self.w[osmid][0]
            tags = self.w[osmid][1]
            away =  Way(refs, tags, coords)
            if away.valid:
                self.way[osmid] = away

    def show(self, typ=2, **kwargs):
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

        lpoly  = []
        lonmin = 360
        lonmax = -360
        latmin = 360
        latmax = -360
        for b in self.way:
            if typ==0:
                if self.way.typ == 0:
                    p =ways.way[b].shp
            if typ==1:
                if self.way.typ == 1:
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
                print(k,N)
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
                print("buildingpart found")
                try:
                    poly  = self.way[b].shp
                    if 'building:roof:colour' in tags:
                        col = '#'+tags['building:roof:colour']
                    else:
                        col = '#abcdef'
                    fig, ax = poly.plot(fig=fig, ax=ax,color=col)
                except:
                    print("building: ",b," is not a polygon")
        plt.axis('scaled')
        return(fig,ax)

    def readmap1(self, osmmap, coords):
        """ read ways from a map

        osmmap : OSM Map from josm file
        coords : coords object previously parsed

        """
        for item in osmmap:
            if item['type']=='way':
                way = item['data']
                osmid = way['id']
                refs_neg = way['nd']
                # nodes should have negative index (PyLayers convention)
                tags  = way['tag']
                if 'z' in tags:
                   z = tags['z']
                   if type(z) == str:
                       z = eval(z)
                   if type(z[0])==str:
                       z = (eval(z[0]),eval(z[1]))
                   tags['z'] = z

                self.w[osmid] = [refs_neg,tags]
                self.cpt += 1

        self.toway(coords)


    def readmap2(self, osmmap, coords, typ='building'):
        """ read ways from a map

        osmmap : OSM Map in json from OsmAPI
        coords : coords object previously parsed
        typ  : string
            type of ways to capture

        """
        ltree =[]
        for item in osmmap:
            # parse trees
            if item['type']=='node':
                node = item['data']
                tags = node['tag']
                if 'natural' in tags:
                    if tags['natural']=='tree':
                        ltree.append((node['lon'],node['lat']))

            if item['type']=='way':
                way = item['data']
                tags = way['tag']
                if typ in tags:
                    osmid = way['id']
                    refs_neg = way['nd']
                    refs_neg = [ -x for x in refs_neg if x >0]
                    ntags = {}
                    # nodes should have negative index (PyLayers convention)
                    if 'height' in tags:
                        ntags['height'] = tags['height']
                    elif 'building:height' in tags:
                        ntags['height']  = tags['building:height']

                    # min_height : from ground to roof top
                    if 'building:min_height' in tags:
                        ntags['min_height'] = tags['building:min_height']
                    elif 'min_height' in tags:
                        ntags['min_height'] = tags['min_height']

                    # material : from ground to roof top
                    if 'building:material' in tags:
                        ntags['material'] = tags['building:material']
                    elif 'material' in tags:
                        ntags['material'] = tags['material']
                    if 'z' in tags:
                       z = tags['z']
                       if type(z) == str:
                           z = eval(z)
                       if type(z[0])==str:
                           z = (eval(z[0]),eval(z[1]))
                       ntags['z'] = z

                    #print('osmparser readmap2',ntags)
                    self.w[osmid] = [refs_neg, ntags]
                    self.cpt += 1

        self.toway(coords)
        return(ltree)

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

    Methods
    -------

    build : recursive construction of floor plan
    show : show the floor plan

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

        nid : int 
        fig : plt.figure 
        ax  : plt.ax

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
#     getbdg
#
#
def getosm(**kwargs):
    """ get osm region from osmapi

    Parameters
    ----------

    address : string
    latlon : tuple or 0
    dist_m : float
    bcart : boolean
    level_height :  float
        typical level height for deriving building height from # levels
    typical_height : float
        typical height for building when no information

    Returns
    -------

        coords
        nodes
        ways
        dpoly
        m
        latlon : tuple or 0

    Notes
    -----

    There are 3 ways to read an Open Street Map structure

    1 - From an osm file ex : filename = 'B11.osm'
    2 - From an osm (lat,lon) string or tuple of float
    3 - From an osm address string

    if latlon tuple is precised, it has priority over the address string

    """

    filename = kwargs.pop('filename','')
    bcart = kwargs.pop('cart', False)
    typ = kwargs.pop('typ','indoor')
    level_height = kwargs.pop('level_height', 3.45)
    typical_height = kwargs.pop('typical_height', 10)
    bexcluded = kwargs.pop('bexcluded', False)

    # from coordinates
    if filename == '':
        address = kwargs.pop('address','Rennes')
        latlon = kwargs.pop('latlon', 0)
        dist_m = kwargs.pop('dist_m', 400)
        rad_to_deg = (180/np.pi)

        if latlon == 0:
            resp = geo.arcgis(address)
            try:
                lat = resp.geojson['features'][0]['properties']['lat']
                lon = resp.geojson['features'][0]['properties']['lng']
            except:
                print(place)
        else:
            lat = latlon[0]
            lon = latlon[1]

        rad_to_deg = (180/np.pi)
        # if latlon == 0:
        #     place = geo.osm(address)
        #     try:
        #         lat, lon = place.latlng
        #     except:
        #         print(place)
        # else:
        #     lat = latlon[0]
        #     lon = latlon[1]

        r_earth = 6370e3

        alpha = (dist_m/r_earth)*rad_to_deg
    #
    # get map from OsmApi (Database query)
    # same extension in longitude and latitude
    #
        Osm = OsmApi()
        osmmap = Osm.Map(lon-alpha, lat-alpha, lon+alpha, lat+alpha)

    else:
    #
    # get map from osm (xml) file
    # type : 'node'
    #        'ways'
    #
        latlon = 0
        e = xml.parse(filename).getroot()

        osmmap = []

        lnode_key = ['id', 'lat', 'lon']

        lway_key = ['id']

        for i in e:
            d = {}
            d['type'] = i.tag
            d['data'] = i.attrib
            #print(i.tag)

            if d['type'] == 'node':
                for k in lnode_key:
                    try:
                        d['data'][k] = eval(d['data'][k])
                    except:
                        pass
                    if k == 'id':
                        if not 'action' in d['data']:
                            d['data'][k] = -d['data'][k]
                d['data']['tag'] = {}

            elif d['type'] == 'way':
                lk = i.getchildren()
                nd = []
                tag = {}
                for k in lk:
                    if k.tag == 'nd':
                        nd.append(eval(k.get('ref')))
                    if k.tag == 'tag':
                        tag[k.get('k')] = k.get('v')
                d['data']['nd'] = nd
                d['data']['tag'] = tag

                # for k in lway_key:
                #     lk = k.get_children()
                #     print(lk)

                    # d['data'][k]=eval(d['data'][k])
                # d['data']['visible']=eval(d['data']['visible'])
            osmmap.append(d)

    nodes = Nodes()
    nodes.clean()
    nodes.readmap(osmmap)

    coords = Coords()
    coords.clean()
    coords.from_nodes(nodes)

    m = coords.cartesian(cart=bcart)
    ways = Ways()
    ways.clean()

    lat = coords.latlon[list(coords.latlon.keys())[0]][0]
    lon = coords.latlon[list(coords.latlon.keys())[0]][1]

    if typ == 'indoor':
        ways.readmap1(osmmap, coords)
    else:
        ltree = ways.readmap2(osmmap, coords)

    # list of nodes involved in buildings
    lnodes_id=[]
    for iw in ways.w:
        lnodes_id += ways.w[iw][0]

    # list of all nodes involved in buildings
    lnodes_id  = np.unique(np.array(lnodes_id))
    # list of all nodes
    lnodes_full = np.unique(np.array(list(coords.latlon.keys())))
    # intersection True if nodes not involved in building
    mask = np.in1d(lnodes_full, lnodes_id, invert=True)

    # nodes not involved in buildings
    if typ != 'indoor' and bexcluded:
        lexcluded = lnodes_full[mask]
        coords.filter(lexcluded)

    dpoly = {}
    for iw in ways.w:
        ways.way[iw].tags = {}
        #
        # slab and materials extraction
        #
        if 'material' in ways.w[iw][1]:
             ways.way[iw].tags['slab']=ways.w[iw][1]['material']
        elif 'building:material' in ways.w[iw][1]:
             ways.way[iw].tags['slab']=ways.w[iw][1]['building:material']
        else:
             ways.way[iw].tags['slab']='WALL'

        # min_height
        if 'building:min_height' in ways.w[iw][1]:
            min_height = eval(ways.w[iw][1]['building:min_height'])
        else:
            min_height = 0
        # height
        if 'height' in ways.w[iw][1]:
            st_height = ways.w[iw][1]['height'].replace('m','')
            ways.way[iw].tags['z'] = (min_height, eval(st_height))
        elif 'building:height' in ways.w[iw][1]:
            ways.way[iw].tags['z'] = (min_height, eval(ways.w[iw][1]['building:height']))
        elif 'building:levels' in ways.w[iw][1]:
            nb_levels = eval(ways.w[iw][1]['building:levels'])
            if type(nb_levels)!=int:
                try:
                    nb_levels = max(nb_levels)
                except:
                    nb_levels=2
                ways.way[iw].tags['z'] = (min_height, nb_levels*level_height)
            elif 'levels' in ways.w[iw][1]:
                nb_levels = eval(ways.w[iw][1]['levels'])
                if type(nb_levels)!=int:
                    try:
                        nb_levels=max(nb_levels)
                    except:
                        nb_levels=2
                    ways.way[iw].tags['z'] = (min_height, nb_levels*level_height)
            else:
                ways.way[iw].tags['z'] = (0,typical_height)

        ptpoly = [coords.xy[x] for x in ways.w[iw][0]]
        dpoly[iw] = geu.Polygon(ptpoly, vnodes=ways.w[iw][0])
        dpoly[iw].coorddeter()

    return coords, nodes, ways, m, (lat,lon), dpoly
#
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
    print(command)
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

    coords,nodes,ways,relation,m = osmparse(fileosm,typ='outdoor',verbose=verbose)
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
        if tags['type']=='outdoor':
            print("Constructing Indoor building ", bid)
            bdg = FloorPlan(bid,coords,nodes,ways,relations)
            bdg.build(typ='relation',eid=bid)
    return bdg,m

