#-*- coding:utf-8 -*-
#
#
#   ,Layout Module
#
#   unittesting in tests/test_layout_u.py
#
"""
.. currentmodule:: pylayers.gis.layout

.. autosummary::

"""
from __future__ import print_function
try:
    from tvtk.api import tvtk
    from mayavi import mlab
except:
    print(',Layout:Mayavi is not installed')

import pdb
import sys
import os
import copy
import glob
import time
import pickle
import tqdm
import numpy as np
import numpy.random as rd
import scipy as sp
import scipy.sparse as sparse
import doctest
import triangle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
import networkx as nx
import pandas as pd
from itertools import combinations, product
import ast
from networkx.readwrite import write_gpickle, read_gpickle
from mpl_toolkits.basemap import Basemap
import shapely.geometry as sh
import shapefile as shp
from shapely.ops import cascaded_union
from matplotlib.collections import PolyCollection,LineCollection
from numpy import array
import PIL.Image as Image
import hashlib
import pylayers.gis.kml as gkml
#from pathos.multiprocessing import ,ProcessingPool as Pool
#from pathos.multiprocessing import cpu_count
from functools import partial

if sys.version_info.major==2:
    from  urllib2 import urlopen
    import ConfigParser
else:
    from  urllib.request import urlopen
    import configparser as ConfigParser
# from c,StringIO import StringIO
# from multiprocessing import ,Pool

def _pickle_method(method):
	func_name = method.im_func.__name__
	obj = method.im_self
	cls = method.im_class
	if func_name.startswith('__') and not func_name.endswith('__'): #deal with mangled names
		cls_name = cls.__name__.lstrip('_')
		func_name = '_' + cls_name + func_name
	return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
	for cls in cls.__mro__:
		try:
			func = cls.__dict__[func_name]
		except KeyError:
			pass
		else:
			break
	return func.__get__(obj, cls)

import types
if sys.version_info.major==2:
    import copy_reg
    copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
else:
    import copyreg
    copyreg.pickle(types.MethodType, _pickle_method, _unpickle_method)

import pylayers.antprop.slab as sb
from pylayers.util import geomutil as geu
from pylayers.util import plotutil as plu
from pylayers.util import pyutil as pyu
from pylayers.util import graphutil as gru
from pylayers.util import cone

# Handle furnitures

import pylayers.gis.furniture as fur
import pylayers.gis.osmparser as osm
from pylayers.gis.selectl import SelectL
import pylayers.util.graphutil as gph
import pylayers.util.project as pro
from pylayers.util.project import logger

def pbar(verbose,**kwargs):
    if verbose:
        pbar=tqdm.tqdm(**kwargs)
        return pbar


class Layout(pro.PyLayers):
    """ Handling Layout

    Attributes
    ----------

    Gs     : Graph of points and segment (structure)
    Gt     : Graph of convex cycles      (topology)
    Gv     : Graph of visibility         (visibility)
    Gi     : Graph of interactions       (interactions)
    Gr     : Graph of rooms              (rooms)
    Nnode  : Number of nodes of Gs
    Nedge  : Number of edges of Gs
    pt     : points sequence
    tahe   : tail head


    Notes
    -----

    This class uses `networkx` to store Layout information

    Gs : structure
    Gt : topology
    Gv : visibility
    Gi : interaction
    Gr : room
    Gm :
    Gw : ways

    Np
    Ns
    Nss

    ax  : (xmin,ymin,xmax,ymax)
    axn : (0,Dx,0,Dy)

    filefur
    filegeom
    filematini
    fileslabini
    hasboundary
    segboundary
    min_sx
    min_sy
    max_sx
    max_sy
    labels
    lbltg
    lboundary
    listtransition
    loadosm
    lsss
    name
    normal
    p2pc
    pg

    pt : points coordinates
    tahe : segment tail head
    tgs : graph to segment
    tsg : segment to graph
    upnt : array of point index

    s2pc : segment to point coordinates
    s2pu : segment to point index
    sgsg

    sl


    typ  : 'indoor' | 'outdoor'
    coordinates : 'cart','lonlat'
    version
    _filename
    _hash

    _shseg : keys / segment index
             values / shapely LineString
    dca    : keys / Gt node
             values / list of air wall
    degree : keys / point degree
             values / array of index
    display : dictionnary for controling various visualization
    dsseg :

    indoor : if True allow indoor penetration
    isbuilt
    diffraction

    maxheight
    zceil
    zfloor
    zmin

    """
    def __init__(self,arg='',**kwargs):
        """ object constructor

        Parameters
        ----------

        arg : string or tuple
            layout file name, address or (lat,lon) or '(lat,lon)'
        mat :
            material dB file name
        slab :
            slab dB file name
        fur :
            furniture file name
        force : boolean
        check : boolean
        build : boolean
        verbose : boolean
        bcartesian : boolean
        xlim : '(xmin,xmax,ymin,ymax) | () default'
        dist_m : int
        typ : string
            'indoor' | 'outdoor'

        """
        self.arg = arg

        self._filematini = kwargs.pop('mat','matDB.ini')
        self._fileslabini = kwargs.pop('slab','slabDB.ini')
        self._filefur = kwargs.pop('fur','')
        self.bcheck = kwargs.pop('bcheck',False)
        self.bbuild = kwargs.pop('bbuild',False)
        self.bgraphs = kwargs.pop('bgraphs',False)
        self.bverbose = kwargs.pop('bverbose',False)
        self.bcartesian = kwargs.pop('bcartesian',True)
        self.xlim = kwargs.pop('xlim',())
        self.dist_m = kwargs.pop('dist_m',400)
        self.typ = kwargs.pop('typ','outdoor')
        self.bexcluded = kwargs.pop('bexcluded',False)

        self.labels = {}

        self.Np = 0
        self.Ns = 0
        self.Nss = 0
        self.lsss = []

        #
        # Initializing graphs
        #  Gs : structure
        #  Gt : topology
        #  Gi : interaction
        #  Gr : rooms
        #  Gm : mobility

        self.Gs = nx.Graph(name='Gs')
        self.Gr = nx.Graph(name='Gr')
        self.Gt = nx.Graph(name='Gt')
        self.Gm = nx.Graph(name='Gm')

        self.Gs.pos = {}
        self.tahe = np.zeros(([2, 0]), dtype=int)
        self.lbltg = []

        self.Gt.pos = {}
        self._shseg = {}

        self.hasboundary = False
        self.format = 'cart'
        self.version = '1.3'

        assert(self.typ in ['indoor','outdoor','floorplan'])

        self.isbuilt = False
        self.loadosm = False

        #
        # setting display option

        self.display = {}
        self.display['title'] = ''
        self.display['ticksoff'] = True
        self.display['nodes'] = False
        self.display['ndsize'] = 10
        self.display['ndlabel'] = False
        self.display['ndlblsize'] = 20
        self.display['edlblsize'] = 20
        self.display['fontsize'] = 20
        self.display['edlabel'] = False
        self.display['edges'] = True
        self.display['ednodes'] = False
        self.display['subseg'] = True
        self.display['isonb'] = True
        self.display['transition'] = True
        self.display['visu'] = False
        self.display['thin'] = False
        self.display['scaled'] = True
        self.display['alpha'] = 0.5
        self.display['layer'] = []
        self.display['clear'] = False
        self.display['activelayer'] = 'AIR'
        self.display['layers'] = []
        self.display['overlay'] = False
        self.display['overlay_flip'] = ""
        # self.display['overlay_file']="/home/buguen/Pyproject/data/image/"
        self.display['overlay_file'] = ""
        self.display['overlay_axis'] = ""
        # self.display['layerset'] = self.sl.keys()

        if self.xlim!=():
            self.display['box']= self.xlim
        else:
            self.display['box'] = (-50, 50, -50, 50)

        self.name = {}
        self.ax = self.display['box']
        self.zmin = 0
        self.maxheight = 3.

        newfile = False
        loadlay = False
        loadosm = False
        loadres = False

        #
        # Layout main argument
        #
        if type(self.arg)==tuple:
            self.arg = str(self.arg)
        #
        # Layout Point of Interest DataFrame
        #
        # A Layout is equipped with a DataFrame of Points of Interest
        #
        # dfpoi   type =['tree','human','tx','rx','support']
        #
        self.dfpoi = pd.DataFrame(columns=['name','type','lon','lat','alt','x',',y','z','radius'])
        if type(self.arg) is bytes:
            self.arg = self.arg.decode('utf-8')

        logger.info('analyze input string '+ self.arg)
        arg, ext = os.path.splitext(self.arg)

        if arg != '':
            if ext == '.ini':
                self._filename = self.arg
                loadlay = True
            if ext == '.lay':
                self._filename = self.arg
                loadlay = True
            elif ext == '.osm':
                self._filename = arg + '.lay'
                loadosm = True
            elif ext == '.res':
                self._filename = arg + '.lay'
                loadres = True
            else:
                self.typ = 'outdoor'
                self._filename = arg.replace(' ','_')+'.lay'
        else:  # No argument new file
            self._filename = 'newfile.lay'
            newfile = True
            self.sl = sb.SlabDB(fileslab=self._fileslabini, filemat=self._filematini)
            self.zfloor = 0.
            self.zceil = self.maxheight

        #
        # build directory to store graphs
        #

        if os.path.splitext(self._filename)[1] == '.lay':
            dirname = self._filename.replace('.lay','')
            path = os.path.join(pro.basename, 'struc', 'gpickle', dirname)

        #
        # load layout from file
        #
        if not newfile:
            if loadlay:
                filename = pyu.getlong(self._filename, pro.pstruc['DIRLAY'])
                if os.path.exists(filename):  # which exists
                    # check the .lay hash
                    fd = open(filename,'rb')
                    self._hash = hashlib.md5(fd.read()).hexdigest()
                    fd.close()
                    # check existence of gpickle graphs
                    if os.path.exists(os.path.join(path,'Gs.gpickle')):
                        if os.path.exists(os.path.join(path,'.hash')):
                            fd = open(os.path.join(path,'.hash'),'r')
                            lines = fd.readlines()
                            fd.close()
                            hashGs = lines[0].replace('\n','')
                            self.segboundary = eval(lines[1])
                            self.typ = lines[2]
                            if hashGs==self._hash:
                                logger.info('load from Gs graph')
                                self.hasboundary = True
                                self.dumpr('s')
                                self.lboundary = []
                                for ns in self.segboundary:
                                    for nh in nx.neighbors(self.Gs,ns):
                                        if nh not in self.lboundary:
                                            self.lboundary.append(nh)
                            else:
                                logger.info('load from .lay file')
                                self.load_fast()
                            self.get_boundary()
                        else:
                            logger.info('load from .lay file')
                            self.load_fast()
                    else:
                        logger.info('load from .lay file')
                        self.load_fast()

                    #self.load()
                else:  # which do not exist
                    newfile = True
                    print("new file - creating a void Layout", self._filename)
            elif loadosm:  # load .osm file
                self.importosm(fileosm=self.arg, cart=self.bcartesian, typ=self.typ)
                self.loadosm = True
            elif loadres:
                self.importres(_fileres=self.arg)
                self.sl = sb.SlabDB()
            elif '(' in str(arg):  # load from osmapi latlon (string or tuple
                latlon = eval(self.arg)
                self.importosm(latlon=latlon,
                               dist_m=self.dist_m,
                               cart=self.bcartesian,
                               typ=self.typ)
                self.loadosm = True
            else:  # load from address geocoding
                self.importosm(address=self.arg,
                        dist_m=self.dist_m,
                        cart=self.bcartesian,
                        typ=self.typ)

                self.loadosm = True

            # add boundary if it not exist
            if (not self.hasboundary) or (self.xlim != ()):
                logger.info('call boundary')
                self.boundary(xlim = self.xlim)

            logger.info('call subseg ')
            self.subseg()

            logger.info('call updateshseg ')
            self.updateshseg()

            try:
                self.geomfile()
            except:
                print("problem to construct geomfile")

            #
            # check layout
            #
            self.bconsistent = True
            if self.bcheck:
                self.bconsistent, dseg = self.check()

            # if Layout is correctly described
            # check if the graph gpickle files have been built
            if self.bconsistent:
                #
                # build and save graphs
                #
                if self.bbuild:
                    # ans = raw_input('Do you want to build the layout (y/N) ? ')
                    # if ans.lower()=='y'
                    self.build()
                    self.lbltg.append('s')
                    self.dumpw()
                #
                # load graphs from file
                #
                elif self.bgraphs:
                    if os.path.splitext(self._filename)[1]=='.lay':
                        dirname = self._filename.replace('.lay','')
                    path = os.path.join(pro.basename,
                                        'struc',
                                        'gpickle',
                                        dirname)
                    if os.path.exists(path):
                        # load graph Gt
                        # and compare the self._hash from ini file
                        # with the hash store in node 0 of Gt at time of the last build
                        # If they are different a rebuild is needeed
                        # Otherwise all the stored graphs are loaded
                        #
                        self.dumpr('t')
                        # If node 0 exists : the layout has been built

                        # If .ini file has changed rebuild
                        if self._hash == self.Gt.node[0]['hash']:
                            self.dumpr('stvirw')
                            self.isbuilt = True
                            bbuild = False
                        else:
                            print(".lay file has changed you must rebuild the grahs")
                    else:
                        # if graph are requested and it not exists a pickle of a graph
                        # they are built
                        self.build()
                        self.lbltg.append('s')
                        self.dumpw()
        #self.get_boundary()
        logger.info('end of __init__')

    def __repr__(self):
        st = '\n'
        st = st + "----------------\n"
        home = os.path.expanduser('~')
        with open(os.path.join(home, '.pylayers'),'r') as f:
            paths = f.readlines()
        uporj = paths.index('project\n')
        project = paths[uporj+1]
        st = st + "Project : " + project+'\n'
        if hasattr(self, '_hash'):
            st = st + self._filename + ' : ' + self._hash + "\n"
        else:
            st = st + self._filename + "\n"


        if self.isbuilt:
            st = st + 'Built with : ' + self.Gt.node[0]['hash'] + "\n"
        st = st + 'Type : ' + self.typ+'\n'

        if self.display['overlay_file'] != '':
            filename = pyu.getlong(
                self.display['overlay_file'], os.path.join('struc', 'images'))
            st = st + "Image('" + filename + "')\n"
        st = st + "Coordinates : " + self.format + "\n"
        if hasattr(self,'extent'):
            st = st + "----------------\n"
            st = st+ str(self.extent)+'\n'
        if hasattr(self,'extent_c'):
            st = st + "----------------\n"
            st = st+ str(self.extent_c)+'\n'
        if hasattr(self, 'Gs'):
            st = st + "----------------\n"
            st = st + "Gs : "+str(len(self.Gs.node))+"("+str(self.Np)+'/'+str(self.Ns)+'/'+str(len(self.lsss))+') :'+str(len(self.Gs.edges()))+'\n'
        if hasattr(self,'Gt'):
            st = st + "Gt : "+str(len(self.Gt.node))+' : '+str(len(self.Gt.edges()))+'\n'
        if hasattr(self,'Gv'):
            st = st + "Gv : "+str(len(self.Gv.node))+' : '+str(len(self.Gv.edges()))+'\n'
        if hasattr(self,'Gi'):
            st = st + "Gi : "+str(len(self.Gi.node))+' : '+str(len(self.Gi.edges()))+'\n'
        if hasattr(self,'Gr'):
            st = st + "Gr : "+str(len(self.Gr.node))+' : '+str(len(self.Gr.edges()))+'\n'
        if hasattr(self,'Gw'):
            st = st + "Gw : "+str(len(self.Gw.node))+' : '+str(len(self.Gw.edges()))+'\n'
        st = st + "----------------\n\n"
        # if hasattr(self, 'degree'):
        #     for k in self.degree:
        #         if (k < 2) or (k > 3):
        #             st = st + 'degree ' + \
        #                 str(k) + ' : ' + str(self.degree[k]) + "\n"
        #         else:
        #             st = st + 'number of node points of degree ' + \
        #                 str(k) + ' : ' + str(len(self.degree[k])) + "\n"
        st = st + "\n"
        st = st + "xrange : " + str(self.ax[0:2]) + "\n"
        st = st + "yrange : " + str(self.ax[2:]) + "\n"
        if hasattr(self,'pg'):
            st = st + "center : " + "( %.2f,%.2f)" % (self.pg[0],self.pg[1]) + "\n"
        if hasattr(self,'radius'):
            st = st + "radius : %.2f " % self.radius + "\n"
        # st = st + "\nUseful dictionnaries" + "\n----------------\n"
        # if hasattr(self,'dca'):
        #     st = st + "dca {cycle : []} cycle with an airwall" +"\n"
        # if hasattr(self,'di'):
        #     st = st + "di {interaction : [nstr,typi]}" +"\n"
        # if hasattr(self,'sl'):
        #     st = st + "sl {slab name : slab dictionary}" +"\n"
        # if hasattr(self,'name'):
        #     st = st + "name :  {slab :seglist} " +"\n"
        # st = st + "\nUseful arrays"+"\n----------------\n"
        # if hasattr(self,'pt'):
        #     st = st + "pt : numpy array of points " +"\n"
        # if hasattr(self,'normal'):
        #     st = st + "normal : numpy array of normal " +"\n"
        # if hasattr(self,'offset'):
        #     st = st + "offset : numpy array of offset " +"\n"
        # if hasattr(self,'tsg'):
        #     st = st + "tsg : get segment index in Gs from tahe" +"\n"
        # if hasattr(self,'isss'):
        #     st = st + "isss :  sub-segment index above Nsmax"+"\n"
        # if hasattr(self,'tgs'):
        #     st = st + "tgs : get segment index in tahe from self.Gs" +"\n"
        # if hasattr(self,'upnt'):
        #     st = st + "upnt : get point id index from self.pt"+"\n"
        # #if hasattr(self,'iupnt'):
        # #    st = st + "iupnt : get point index in self.pt from point id  "+"\n"
        # if hasattr(self,'lsss'):
        #     st = st + "lsss : list of segments with sub-segment"+"\n"
        # if hasattr(self,'sridess'):
        #     st = st + "stridess : stride to calculate the index of a subsegment" +"\n"
        # if hasattr(self,'sla'):
        #     st = st + "sla : list of all slab names (Nsmax+Nss+1)" +"\n"
        # if hasattr(self,'degree'):
        #     st = st + "degree : degree of nodes " +"\n"
        # st = st + "\nUseful tip" + "\n----------------\n"
        # st = st + "Point p in Gs => p_coord:\n"
        # #st = st + "p -> u = self.iupnt[-p] -> p_coord = self.pt[:,u]\n\n"
        #st = st + "Segment s in Gs => s_ab coordinates \n"
        #st = st + "s2pc : segment to point coordinates (sparse) [p1,p2] = L.s2pc.toarray().reshape(2,2).T \n"
        #st = st + \
        #    "s -> u = self.tgs[s] -> v = self.tahe[:,u] -> s_ab = self.pt[:,v]\n\n"
        return(st)

    def __add__(self, other):
        """ addition

        One can add either a numpy array or an other layout

        """
        Ls = copy.deepcopy(self)
        if type(other) == np.ndarray:
            for k in Ls.Gs.pos:
                Ls.Gs.pos[k] = Ls.Gs.pos[k] + other[0:2]
        else:
            offp = -min(Ls.Gs.nodes())
            offs = max(Ls.Gs.nodes())
            other.offset_index(offp=offp, offs=offs)
            Ls.Gs.node.update(other.Gs.node)
            Ls.Gs.edge.update(other.Gs.edge)
            Ls.Gs.adj.update(other.Gs.adj)
            Ls.Gs.pos.update(other.Gs.pos)
            Ls.Np = Ls.Np + other.Np
            Ls.Ns = Ls.Ns + other.Ns
            Ls.Nss = Ls.Nss + other.Nss

        return(Ls)

    def __mul__(self, alpha):
        """ scale the layout

        other : scaling factor (np.array or int or float)

        Returns
        -------

        Ls : Layout
            scaled layout

        """
        Ls = copy.deepcopy(self)
        Gs = Ls.Gs
        if type(alpha) != np.ndarray:
            assert((type(alpha) == float) or (
                type(alpha) == int)), " not float"
            alpha = np.array([alpha, alpha, alpha])
        else:
            assert(len(alpha) == 3), " not 3D"
        #
        # scaling x & y
        #
        x = np.array(Gs.pos.values())[:, 0]
        x = x * alpha[0]

        y = np.array(Gs.pos.values())[:, 1]
        y = y * alpha[1]

        xy = np.vstack((x, y)).T
        Ls.Gs.pos = dict(zip(Gs.pos.keys(), tuple(xy)))

        #
        # scaling z
        #

        nseg = filter(lambda x: x > 0, Gs.nodes())
        for k in nseg:
            Ls.Gs.node[k]['z'] = tuple(
                (np.array(Ls.Gs.node[k]['z']) - self.zmin) * alpha[2] + self.zmin)
            if 'ss_z' in Ls.Gs.node[k]:
                Ls.Gs.node[k]['ss_z'] = list(
                    (np.array(Ls.Gs.node[k]['ss_z']) - self.zmin) * alpha[2] + self.zmin)

        #
        # updating numpy array from graph
        #

        Ls.g2npy()
        return Ls

    def mapconv(xorlon,yorlat,inverse=False):
        """ map converter spherical to cartesian (default)

        Parameters
        ----------

        lonorx : np.array of float
            lon or x
        latory : np.array of float
            lat or y

        inverse : boolean
           True  (cartesian -> spherical)
           False (spherical -> cartesian)

        Examples
        --------

        lon = np.array([])
        lat = np.array([])
        p = L.mapconv(lon,lat)

        Info
        ----

        Lambert 93 (France) : Proj(init='epsg:2154')

        """
        # spherical projection
        psph = Proj(init='epsg:4326')
        # use Basemap (Deprecated)
        if hasattr(self,'m'):
            p = self.m(lonorx,latory,inverse=inverse)
        # use pyproj
        if hasattr(self,'cartproj'):
            if inverse:
                p = projection(self.cartproj,psph,lonorx,latory)
            else:
                p = projection(psph,self.cartproj,lonorx,latory)

    def switch(self):
        """ switch coordinates

        Info
        ----

        Use either basemap objet or pyproj object

        """
        if hasattr(self,'m'):
            if self.format == 'cart':
                for k in self.Gs.pos.keys():
                    p = self.m( self.Gs.pos[k][0], self.Gs.pos[k][1], inverse=True)
                    self.Gs.pos[k] = p
                if hasattr(self,'dpoly'):
                    if hasattr(self.dpoly,'_xy'):
                        for k in self.dpoly:
                        #self.dpoly[k].ndarray() = np.vstack(self.m(self.dpoly[k].ndarray()[0,:],self.dpoly[k].ndarray()[1,:],inverse=True))
                            self.dpoly[k]._xy = np.vstack(self.m(self.dpoly[k]._xy[0,:],self.dpoly[k]._xy[1,:],inverse=True))
                self.format = 'latlon'
            elif self.format == 'latlon':
                for k in self.Gs.pos.keys():
                    p = self.m( self.Gs.pos[k][0], self.Gs.pos[k][1])
                    self.Gs.pos[k] = p
                if hasattr(self,'dpoly'):
                    if hasattr(self.dpoly,'_xy'):
                        for k in self.dpoly:
                        #self.dpoly[k].ndarray() = np.vstack(self.m(self.dpoly[k].ndarray()[0,:],self.dpoly[k].ndarray()[1,:],inverse=True))
                            self.dpoly[k]._xy = np.vstack(self.m(self.dpoly[k]._xy[0,:],self.dpoly[k]._xy[1,:]))
                self.format ='cart'

            nodes = self.Gs.nodes()
            upnt = [n for n in nodes if n < 0]
            self.pt[0, :] = np.array([self.Gs.pos[k][0] for k in upnt])
            self.pt[1, :] = np.array([self.Gs.pos[k][1] for k in upnt])


    def _help(self):
        st = ''
        st = st + "\nUseful dictionnaries" + "\n----------------\n"
        if hasattr(self,'dca'):
            st = st + "dca {cycle : []} cycle with an airwall" +"\n"
        if hasattr(self,'di'):
            st = st + "di {interaction : [nstr,typi]}" +"\n"
        if hasattr(self,'sl'):
            st = st + "sl {slab name : slab dictionary}" +"\n"
        if hasattr(self,'name'):
            st = st + "name :  {slab :seglist} " +"\n"
        st = st + "\nUseful arrays"+"\n----------------\n"
        if hasattr(self,'pt'):
            st = st + "pt : numpy array of points " +"\n"
        if hasattr(self,'normal'):
            st = st + "normal : numpy array of normal " +"\n"
        if hasattr(self,'offset'):
            st = st + "offset : numpy array of offset " +"\n"
        if hasattr(self,'tsg'):
            st = st + "tsg : get segment index in Gs from tahe" +"\n"
        if hasattr(self,'isss'):
            st = st + "isss :  sub-segment index above Nsmax"+"\n"
        if hasattr(self,'tgs'):
            st = st + "tgs : get segment index in tahe from self.Gs" +"\n"
        if hasattr(self,'upnt'):
            st = st + "upnt : get point id index from self.pt"+"\n"

        st = st + "\nUseful Sparse arrays"+"\n----------------\n"
        if hasattr(self,'sgsg'):
            st = st + "sgsg : "+"get common point of 2 segment (usage self.sgsg[seg1,seg2] => return common point \n"
        if hasattr(self,'s2pc'):
            st = st + "s2pc : "+"from a Gs segment node to its 2 extremal points (tahe) coordinates\n"
        if hasattr(self,'s2pu'):
            st = st + "s2pu : "+"from a Gs segment node to its 2 extremal points (tahe) index\n"
        if hasattr(self,'p2pu'):
            st = st + "p2pu : "+"from a Gs point node to its coordinates\n"
        st = st + "\nUseful lists"+"\n----------------\n"
        #if hasattr(self,'iupnt'):
        #    st = st + "iupnt : get point index in self.pt from point id  "+"\n"
        if hasattr(self,'lsss'):
            st = st + "lsss : list of segments with sub-segment"+"\n"
        if hasattr(self,'sridess'): 
            st = st + "stridess : stride to calculate the index of a subsegment" +"\n"
        if hasattr(self,'sla'):
            st = st + "sla : list of all slab names (Nsmax+Nss+1)" +"\n"
        if hasattr(self,'degree'):
            st = st + "degree : degree of nodes " +"\n"
        st = st + "\nUseful tip" + "\n----------------\n"
        st = st + "Point p in Gs => p_coord: Not implemented\n"
        # st = st + "p -> u = self.upnt[-p] -> p_coord = self.pt[:,-u]\n\n"
        st = st + "Segment s in Gs => s_ab coordinates \n"
        st = st + \
            "s -> u = self.tgs[s] -> v = self.tahe[:,u] -> s_ab = self.pt[:,v]\n\n"
        print(st)

    def ls(self, typ='lay'):
        """ list the available file in dirstruc

        Parameters
        ----------

        typ : string optional
            {'lay'|'osm'|'wrl'}

        Returns
        -------

        lfile_s : list
            sorted list of all the .str file of strdir

        Notes
        -----

        strdir is defined in the Project module

        Examples
        --------

        Display all available structures


        >>> from pylayers.gis.layout import *
        >>> L = Layout()
        >>> fillist = L.ls()

        """

        if typ == 'lay':
            pathname = os.path.join(pro.pstruc['DIRLAY'], '*.' + typ)
        if typ == 'osm':
            pathname = os.path.join(pro.pstruc['DIROSM'], '*.' + typ)
        if typ == 'wrl':
            pathname = os.path.join(pro.pstruc['DIRWRL'], '*.' + typ)

        lfile_l = glob.glob(os.path.join(pro.basename, pathname))
        lfile_s = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_s.append(fis)
        lfile_s.sort()

        return lfile_s

    def offset_index(self, offp=0, offs=0):
        """ offset points and segment index

        Parameters
        ----------

        offp : offset points
        offs : offset segments

        See Also
        --------
        Portage vers networkx 2. inacheve

        __add__

        """

        newpoint = dict((k - offp, v) for k, v in self.Gs.node.items() if k < 0)
        assert (np.array(list(newpoint.keys())) < 0).all()

        newseg = dict((k + offs, v) for k, v in self.Gs.node.items() if k > 0)
        assert (np.array(list(newseg.keys())) > 0).all()

        newpoint.update(newseg)
        nx.set_node_attributes(self.Gs,newpoint)
        #self.Gs.node = newpoint

        newppoint = dict((k - offp, v) for k, v in self.Gs.pos.items() if k < 0)
        newpseg = dict((k + offs, v) for k, v in self.Gs.pos.items() if k > 0)
        newppoint.update(newpseg)

        self.Gs.pos = newppoint

        # adjascence list of segments
        ladjs = [self.Gs.adj[k] for k in self.Gs.adj.keys() if k > 0]
        # adjascence list of points
        ladjp = [self.Gs.adj[k] for k in self.Gs.adj.keys() if k < 0]

        nladjs = map(lambda x: dict((k - offp, v)
                                    for k, v in x.items()), ladjs)
        nladjp = map(lambda x: dict((k + offs, v)
                                    for k, v in x.items()), ladjp)

        lpt = [k - offp for k in self.Gs.adj.keys() if k < 0]
        lseg = [k + offs for k in self.Gs.adj.keys() if k > 0]

        dpt = dict(zip(lpt, nladjp))
        dseg = dict(zip(lseg, nladjs))
        dseg.update(dpt)
        print("Todo create a dictionnary of edges")
        nx.set_edge_attributes(self.Gs,dseg)
        # self.Gs.adj = dseg
        # self.Gs.edge = dseg

    def check(self, level=0, epsilon = 0.64):
        """ Check Layout consistency


        Parameters
        ----------

        level : int

        Returns
        -------

        consistent : Boolean
              True if consistent
        dseg : dictionnary of segments

        See Also
        --------

        GeomUtil.isBetween

        Notes
        -----

        For all segments
            get the 2 vertices
                for all the other vertices
                    check if it belongs to segment

        If there are points which are not valid they are displayed

        In red point with degree == 1 , In black points with degree == 0

        """

        assert(len(self.Gs.pos)==len(self.Gs.nodes()))
        bconsistent = True

        nodes = self.Gs.nodes()

        if len(nodes) > 0:
            #
            # points
            # segments
            # degree of segments

            useg = [ x for x in nodes if x > 0 ]
            upnt = [ x for x in nodes if x < 0 ]

            degseg = [nx.degree(self.Gs, x) for x in  useg ]

            #
            # 1)   all segments have degree 2
            #
            assert(np.all(array(degseg) == 2))

            #
            # degree of points
            # maximum degree of points
            #

            degpnt = [ nx.degree(self.Gs, x) for x in  upnt ]  # points absolute degrees
            degmin = min(degpnt)
            degmax = max(degpnt)

            #
            #  No isolated points (degree 0)
            #  No points of degree 1
            #
            if (degmin <= 1):
                f, a = self.showG('s', aw=1)
                deg0 = [ x for x in upnt if nx.degree(self.Gs, x) == 0]
                deg1 = [ x for x in upnt if nx.degree(self.Gs, x) == 1]

                if len(deg0) > 0:
                    logger.critical( "It exists degree 0 points :  %r", deg0 )
                    f, a = self.pltvnodes(deg0, fig=f, ax=a)
                    bconsistent = False

                if len(deg1) > 0:
                    logger.critical( "It exists degree 1 points :  %r", deg1 )
                    f, a = self.pltvnodes(deg1, fig=f, ax=a)
                    bconsistent = False

            # self.deg = {}
            # for deg in range(degmax + 1):
            #     num = filter(lambda x: degpnt[x] == deg, range(
            #         len(degpnt)))  # position of degree 1 point
            #     npt = map(lambda x: upnt[x], num)  # number of degree 1 points
            #     self.deg[deg] = npt

            #
            # check if there are duplicate points or segments
            #
            # TODO argsort x coordinate
            #

            # get all the nodes
            ke = list(self.Gs.pos.keys())
            lpos = list(self.Gs.pos.values())

            x = np.array([ pp[0] for pp in  lpos ] )
            y = np.array([ pp[1] for pp in  lpos ] )
            p = np.vstack((x, y))

            d1 = p - np.roll(p, 1, axis=1)
            sd1 = np.sum(np.abs(d1), axis=0)
            if not sd1.all() != 0:
                lu = np.where(sd1 == 0)[0]

                for u in lu:
                    # if ke[u]>0:
                    #     self.del_segment(ke[u])
                    if ke[u] < 0:
                        self.del_points(ke[u])

                nodes = self.Gs.nodes()
                # useg  = filter(lambda x : x>0,nodes)
                upnt = filter(lambda x: x < 0, nodes)

            # iterate on useg : list of segments
            # s : n1 <--> n2
            #
            # Is there a point different from (n1-n2) in betweeen of an existing segment s ?
            #
            # Not scalable. Double for loop

            dseg = {}
            if (self.typ == 'indoor') or (self.typ=='outdoor'):
                for s in useg:
                    # n1, n2 = np.array(self.Gs.neighbors(s))  # node s neighbors
                    n1, n2 = np.array(self.Gs[s])  # node s neighbors
                    p1 = np.array(self.Gs.pos[n1])           # p1 --- p2
                    p2 = np.array(self.Gs.pos[n2])  # s
                    #
                    # iterate on upnt : list of points
                    for n in upnt:
                        if (n1 != n) & (n2 != n):
                            p = np.array(self.Gs.pos[n])
                            if geu.isBetween(p1, p2, p,epsilon=epsilon):
                                if s in dseg:
                                    dseg[s].append(n)
                                else:
                                    dseg[s]=[n]
                                logger.critical("segment %d contains point %d", s, n)
                                bconsistent = False
                    if level > 0:
                        cycle = self.Gs.node[s]['ncycles']
                        if len(cycle) == 0:
                            logger.critical("segment %d has no cycle", s)
                        if len(cycle) == 3:
                            logger.critical(
                                "segment %d has cycle %s", s, str(cycle))

        #
        # check if Gs points are unique
        # segments can be duplicated
        #

        P = np.array([self.Gs.pos[k] for k in upnt])
        similar = geu.check_point_unicity(P)

        if len(similar) != 0:
            logger.critical("points at index(es) %s in self.Gs.pos are similar", str(similar))
            bconsistent = False

        return bconsistent, dseg

    def clip(self, xmin, xmax, ymin, ymax):
        """ return the list of edges which cross or belong to the clipping zone

        Parameters
        ----------

        xmin : float
        xmax : float
        ymin : float
        ymax : float

        Returns
        -------

        seglist : list of segment number

        Notes
        -----

          1) Determine all segments outside the clipping zone
          2) Union of the 4 conditions
          3) setdiff1d between the whole array of segments and the segments outside

        """

        p0 = self.pt[:, self.tahe[0, :]]
        p1 = self.pt[:, self.tahe[1, :]]

        maxx = np.maximum(p0[0, :], p1[0, :])
        maxy = np.maximum(p0[1, :], p1[1, :])

        minx = np.minimum(p0[0, :], p1[0, :])
        miny = np.minimum(p0[1, :], p1[1, :])

        nxp = np.nonzero(maxx < xmin)[0]
        nxm = np.nonzero(minx > xmax)[0]
        nyp = np.nonzero(maxy < ymin)[0]
        nym = np.nonzero(miny > ymax)[0]

        u = np.union1d(nxp, nxm)
        u = np.union1d(u, nyp)
        u = np.union1d(u, nym)

        iseg = np.arange(self.Ns)

        return np.setdiff1d(iseg, u)


    def check_Gi(self):

        for nit1 in self.Gi.nodes():
            if len(nit1)>1:
                cy1 = nit1[-1]
                for nint2 in self.Gi[nit1].keys():
                    if len(nint2) > 1 :
                        assert nint2[1] == cy1


        # for e0,e1 in self.Gi.edges():


    def g2npy(self,verbose=False):
        """ conversion from graphs to numpy arrays

        Parameters
        ----------

        verbose : boolean

        Notes
        -----

        This function updates the following arrays in self:

        + pt   (2xNp)
        + pg   center of gravity
        + tahe (2xNs)
        + tgs  : graph to segment
        + tsg  : segment to graph
        + dca  : dictionnary of cycle with an airwall (_AIR)
        + s2pu : sparse_lil_matrix
        + s2pc : sparse_lil_matrix
        + lsss : list of iso segments
        + maxheight :
        + normal :

        assert self.pt[self.iupnt[-1]] == self.pt[:,self.iupnt[-1]]

        See Also
        --------

        extrseg

        """

        nodes = self.Gs.nodes()
        # nodes include points and segments

        # segment index
        useg = [n for n in nodes if n > 0]

        # points index
        upnt = [n for n in nodes if n < 0]


        # sparse matrix segment-segment
        # usage
        # self.sgsg[seg1,seg2] => return common point

        mno = max(nodes)

        #self.sgsg = sparse.lil_matrix((mno+1,mno+1),dtype='int')

        # loop over segments
        # a segment is always connected to 2 nodes
        logger.info('g2npy : loop over segments')
        for s in useg:
            # get point index of the segment
            # s > 0
            # v1.1 lpts = [ x for x in nx.neighbors(self.Gs,s)]
            lpts = [ x for x in self.Gs[s] ]
            assert(len(lpts)==2)
            assert(lpts[0]<0)
            assert(lpts[1]<0)
            # get point 0 neighbors
            a = [ x for x in self.Gs[lpts[0]]]
            # a = self.Gs.edge[lpts[0]].keys()
            # get point 1 neighbors
            #b = self.Gs.edge[lpts[1]].keys()
            b = [ x for x in self.Gs[lpts[1]]]

            nsa = np.setdiff1d(a,b)
            nsb = np.setdiff1d(b,a)
            u = np.hstack((nsa,nsb))

            npta = [lpts[0]]*len(nsa)
            nptb = [lpts[1]]*len(nsb)
            ns = np.hstack((npta,nptb))

            #self.sgsg[s,u]=ns


        # conversion in numpy array
        self.upnt = np.array((upnt))
        self.useg = np.array((useg))

        l1 =  [(x, self.Gs.nodes[x]['name']) for x in self.useg ]
        for l in l1:
            try:
                self.name[l[1]].append(l[0])
            except:
                self.name[l[1]] = [l[0]]

        # association

        # utmp = np.array(zip(-self.upnt,np.arange(len(self.upnt))))
        # mutmp = max(utmp[:,0])
        # self.iupnt = -np.ones((mutmp+1),dtype='int')
        # self.iupnt[utmp[:,0]]=utmp[:,1]

        # degree of segment nodes
        degseg = [ nx.degree(self.Gs,x) for x in useg ]

        assert(np.all(np.array(degseg) == 2))  # all segments must have degree 2

        #
        # self.degree : dictionnary (point degree : list of point index)
        #

        # points absolute degrees
        degpnt = np.array([nx.degree(self.Gs, x) for x in  upnt])

        # lairwall : list of air wall segments

        lairwall = []

        if 'AIR' in self.name:
            lairwall += self.name['AIR']
        else:
            self.name['AIR'] = []

        if '_AIR' in self.name:
            lairwall += self.name['_AIR']
        else:
            self.name['_AIR'] = []

        # as self.name['AIR'] and self.name['_AIR'] are tested
        # we define them as void list if not defined

        #
        #  function to count airwall connected to a point
        #  probably this is not the faster solution
        #

        def nairwall(nupt):
            #v1.1 lseg = nx.neighbors(self.Gs, nupt)
            lseg = self.Gs[nupt]
            n = 0
            for ns in lseg:
                if ns in lairwall:
                    n = n + 1
            return n

        nairwall = np.array([ nairwall(x) for x in  upnt])

        if verbose:
            print('buildging nairwall : Done')
        #
        # if a node is connected to N air wall ==> deg = deg - N
        #

        degpnt = degpnt - nairwall

        try:
            degmax = max(degpnt)
        except:
            degmax = 1

        self.degree = {}
        logger.info('g2npy : Start node degree determination')
        for deg in range(degmax + 1):
            #num = filter(lambda x: degpnt[x] == deg, range(
            #    len(degpnt)))  # position of degree 1 point
            num = [ x for x in range(len(degpnt)) if degpnt[x] == deg ]
            # number of degree 1 points
            #npt = np.array(map(lambda x: upnt[x], num))
            npt = np.array([upnt[x] for x in num])
            self.degree[deg] = npt

        logger.info('g2npy : convert geometric information in numpy array')

        self.pt = np.array(np.zeros([2, len(upnt)]), dtype=float)
        self.tahe = np.array(np.zeros([2, len(useg)]), dtype=int)

        self.Np = len(upnt)
        self.Ns = len(useg)

        self.pt[0, :] = np.array([self.Gs.pos[k][0] for k in upnt])
        self.pt[1, :] = np.array([self.Gs.pos[k][1] for k in upnt])

        self.pg = np.sum(self.pt, axis=1) / np.shape(self.pt)[1]
        ptc = self.pt-self.pg[:,None]
        dptc = np.sqrt(np.sum(ptc*ptc,axis=0))
        self.radius  = dptc.max()
        self.pg = np.hstack((self.pg, 0.))

        if self.Ns>0:
            ntahe = np.array([ [n for n in nx.neighbors(self.Gs,x) ]  for x in useg ])
            ntail = ntahe[:,0]
            nhead = ntahe[:,1]

            # create sparse matrix from a Gs segment node to its 2 extremal points (tahe) index
            mlgsn = max(self.Gs.nodes())+1
            self.s2pu = sparse.lil_matrix((mlgsn,2),dtype='int')
            self.s2pu[useg,:] = ntahe
            # convert to compressed row sparse matrix
            # to be more efficient on row slicing
            self.s2pu = self.s2pu.tocsr()


        if self.Ns > 0:
            aupnt = np.array(upnt)
            logger.info('g2npy : build tail head tahe')
            self.tahe[0, :] = np.array([np.where(aupnt==x)[0][0] for x in ntail ])
            self.tahe[1, :] = np.array([np.where(aupnt==x)[0][0] for x in nhead ])

        #
        # transcoding array between graph numbering (discontinuous) and numpy numbering (continuous)
        #
        Nsmax = 0
        self.tsg = np.array(useg)
        try:
            Nsmax = max(self.tsg)
        except:
            logger.warning("No segments in Layout yet")

        #
        # handling of segment related arrays
        #

        if Nsmax > 0:
            self.tgs = -np.ones(Nsmax + 1, dtype=int)
            rag = np.arange(len(useg))
            self.tgs[self.tsg] = rag

            #
            # calculate normal to segment ta-he
            #
            # This could becomes obsolete once the normal will be calculated at
            # creation of the segment
            #

            X = np.vstack((self.pt[0, self.tahe[0, :]],
                           self.pt[0, self.tahe[1, :]]))
            Y = np.vstack((self.pt[1, self.tahe[0, :]],
                           self.pt[1, self.tahe[1, :]]))

            normx = Y[0, :] - Y[1, :]
            normy = X[1, :] - X[0, :]

            scale = np.sqrt(normx * normx + normy * normy)
            assert (scale.all() > 0), pdb.set_trace()
            self.normal = np.vstack(
                (normx, normy, np.zeros(len(scale)))) / scale

            # for ks in ds:
            #
            # lsss : list of subsegment
            #

            # nsmax = max(self.Gs.node.keys())

            # Warning
            # -------
            # nsmax can be different from the total number of segments
            # This means that the numerotation of segments do not need to be
            # contiguous.
            # stridess : length is equal to nsmax+1
            # sla is an array of string, index 0 is not used because there is
            # no such segment number.
            #
            self.lsss = [x for x in useg if len(self.Gs.nodes[x]['iso']) > 0]

            # self.isss = []

            # self.stridess = np.array(np.zeros(nsmax+1),dtype=int)
            # self.stridess = np.empty(nsmax+1,dtype=int)
            # +1 is for discarding index 0 (unused here)
            # self.offset = np.empty(nsmax+1+self.Nss,dtype=int)

            # Storing segment normals
            # Handling of subsegments
            #
            # index is for indexing subsegment after the nsmax value
            #
            # index = nsmax+1
            # for ks in useg:
            #     k = self.tgs[ks]                        # index numpy
            #     self.offset[k] = self.Gs.node[ks]['offset']
            #     self.Gs.node[ks]['norm'] = self.normal[:,k]  # update normal
            #     nameslab  = self.Gs.node[ks]['name']   # update sla array
            #     assert nameslab!='', "segment "+str(ks)+ " is not defined"
            #     self.sla[ks] = nameslab
            #     # stridess is different from 0 only for subsegments
            #     self.stridess[ks] = 0                   # initialize stridess[ks]
            #     #if index==155:
            #     if self.Gs.node[ks].has_key('ss_name'): # if segment has sub segment
            #         nss = len(self.Gs.node[ks]['ss_name'])  # retrieve number of sseg
            #         self.stridess[ks]=index-1           # update stridess[ks] dict
            #         for uk,slabname in enumerate(self.Gs.node[ks]['ss_name']):
            #             self.lsss.append(ks)
            #             self.sla[index] = slabname
            #             self.isss.append(index)
            #             self.offset[index] = self.Gs.node[ks]['ss_offset'][uk]
            #             index = index+1

        # append sub segment normal to normal

        # create sparse matrix from a Gs segment node to its 2 extremal points (tahe) coordinates
        if self.Ns >0:
            self.s2pc = sparse.lil_matrix((mlgsn,4))

            ptail = self.pt[:,self.tahe[0,:]]
            phead = self.pt[:,self.tahe[1,:]]
            A = np.vstack((ptail,phead)).T
            self.s2pc[self.tsg,:]=A

            # convert to compressed row sparse matrix
            # to be more efficient on row slicing
            self.s2pc = self.s2pc.tocsr()
            # for k in self.tsg:
            #     assert(np.array(self.s2pc[k,:].todense())==self.seg2pts(k).T).all(),pdb.set_trace()
            #
            # This is wrong and asume a continuous indexation of points
            # TODO FIX : This problem cleanly
            # 
            # self.p2pc is only used in Gspos in outputGi_func only caled in case of 
            # multiprocessing 
            #
            # The temporary fix is to comment the 5 next lines
            #
            # mino = -min(self.Gs.nodes())+1
            # self.p2pc = sparse.lil_matrix((mino,2))
            # self.p2pc[-self.upnt,:]=self.pt.T
            # self.p2pc = self.p2pc.tocsr()
            # normal_ss = self.normal[:,self.tgs[self.lsss]]
            # self.normal = np.hstack((self.normal,normal_ss))
            # if problem here check file format 'z' should be a string
            lheight = array([v[1] for v in
                        nx.get_node_attributes(self.Gs, 'z').values()
                        if v[1] < 2000 ])
            #assert(len(lheight)>0),logger.error("no valid heights for segments")
            if len(lheight)>0:
                self.maxheight = np.max(lheight)
            else:
                self.maxheight = 3
            # self.maxheight=3.
            # calculate extremum of segments
            self.extrseg()

    def importshp(self, **kwargs):
        """ import layout from shape file

        Parameters
        ----------

        _fileshp :

        """

        pref = kwargs.pop('pref', [np.array([25481100, 6676890]), np.array([60.2043716, 24.6591147])])
        dist_m  = kwargs.pop('dist_m',250)
        latlon = kwargs.pop('latlon',True)
        bd = kwargs.pop('bd', [24, 60, 25, 61])



        fileshp = pyu.getlong(kwargs['_fileshp'], os.path.join('struc', 'shp'))
        polys = shp.Reader(fileshp)
        verts = []
        for poly in polys.iterShapes():
            verts.append(poly.points)

        npt = -1
        ns = 0
        xmin = 1e16
        ymin = 1e16
        xmax = -1e16
        ymax = -1e16

        self.name['WALL'] = []
        for p in verts:
            v = np.array(p) - pref[0][None, :]
            nv = np.sqrt(np.sum(v * v, axis=1))
            # if at least one point is in the radius the poygon is kept
            if (nv < dist_m).any():
                npoint = len(p)
                for k, point in enumerate(p):
                    # add a new node unless it is the last already existing
                    # point
                    if k != (npoint - 1):
                        if k == 0:
                            np0 = npt
                        self.Gs.add_node(npt)
                        x = point[0]
                        y = point[1]
                        xmin = min(x, xmin)
                        xmax = max(x, xmax)
                        ymin = min(y, ymin)
                        ymax = max(y, ymax)
                        self.Gs.pos[npt] = (x, y)
                        npt = npt - 1
                    # add a new segment from the second point
                    if (k > 0) & (k < npoint - 1):
                        ns = ns + 1
                        self.Gs.add_node(ns, name='WALL', z=[
                                         0, 10], offset=0, transition=False, connect=[npt + 1, npt + 2])
                        self.Gs.add_edge(npt + 1, ns)
                        self.Gs.add_edge(ns, npt + 2)
                        self.Gs.pos[ns] = tuple(
                            (np.array(self.Gs.pos[npt + 1]) + np.array(self.Gs.pos[npt + 2])) / 2.)
                    # add a new segment closing the polygon
                    if k == npoint - 1:
                        ns = ns + 1
                        self.Gs.add_node(ns, name='WALL', z=[
                                         0, 10], offset=0, transition=False, connect=[np0, npt + 1])
                        self.Gs.add_edge(np0, ns)
                        self.Gs.add_edge(ns, npt + 1)
                        self.Gs.pos[ns] = tuple(
                            (np.array(self.Gs.pos[npt + 1]) + np.array(self.Gs.pos[np0])) / 2.)
        #
        # TODO change lon_0 and lat_0 hard coded
        #
        self.m = Basemap(llcrnrlon=kwargs['bd'][0], llcrnrlat=kwargs['bd'][1],
                         urcrnrlon=kwargs['bd'][2], urcrnrlat=kwargs['bd'][3],
                         resolution='i', projection='cass', lon_0=24.5, lat_0=60.5)

        if latlon:
            lat_ref = pref[1][0]
            lon_ref = pref[1][1]
            x_ref, y_ref = self.m(lon_ref, lat_ref)
            Dx = pref[0][0] - x_ref
            Dy = pref[0][1] - y_ref
            pos = np.array(self.Gs.pos.values())
            for k, keys in enumerate(self.Gs.pos.keys()):
                self.Gs.pos[keys] = self.m( pos[k, 0] - Dx, pos[k, 1] - Dy, inverse=True)

            self.format = 'latlon'

    def importres(self,_fileres,**kwargs):
        """ import res format

        col1 : x1 coordinates
        col2 : y1 coordinates
        col3 : x2 coordinates
        col4 : y2 coordinates
        col5 : building height
        col6 : building number
        col7 : building class
        col8 : ground height

        Notes
        -----

        COST231 data Munich_buildings.res

        """
        fileres = pyu.getlong(_fileres, os.path.join('struc', 'res'))
        D  = np.fromfile(fileres,dtype='int',sep=' ')
        self.typ = 'outdoor'
        # number of integer
        N1 = len(D)
        # number of lines
        N2 = N1//8
        D = D.reshape(N2,8)
        # list of coordinates
        lcoords = []
        # list of ring
        lring = []
        # list of (z_ground, height_building)
        zring = []
        #
        bdg_old = 1
        self.zfloor=-1000
        self.zceil = 4000
        for e in range(N2):
            # p1 point coordinate
            p1 = ([D[e,0],D[e,1]])
            # p2 point coordinate
            p2 = ([D[e,2],D[e,3]])
            # (ground height,building height+ground_height)
            z  = (D[e,7],D[e,4]+D[e,7])
            # building number
            bdg =  D[e,5]
            # building class
            bdc =  D[e,6]
            # detect change of building
            if (bdg_old-bdg)!=0:
                ring = sh.LinearRing(lcoords)
                poly = sh.Polygon(ring)
                if poly.area>0:
                    lring.append(ring)
                    zring.append(z)
                    lcoords = []
            bdg_old=bdg
            # update lcoords
            if p1 not in lcoords:
                lcoords.append(p1)
            if p2 not in lcoords:
                lcoords.append(p2)

        npt = 1
        self.dpoly = {}
        for kpol,(r1,z1) in enumerate(zip(lring,zring)):
            x,y = r1.xy
            lseg = []

            for k2 in range(len(x)):
                new_pt = (x[k2],y[k2])
                kpos = list(self.Gs.pos.keys())
                vpos = list(self.Gs.pos.values())
                if new_pt not in vpos:
                    #
                    # add node point nde <0 and position
                    #
                    current_node_index = -npt
                    self.Gs.add_node(current_node_index)
                    self.Gs.pos[-npt] = new_pt
                    npt = npt + 1
                else:
                    u = [k for k in range(len(vpos)) if (vpos[k] == new_pt)]

                    current_node_index = kpos[u[0]]

                if k2>0: # at least already one point
                    ns = self.add_segment(current_node_index,
                                          previous_node_index,
                                          name='WALL',
                                          z=z1)
                    lseg.append(ns)
                else:
                    starting_node_index = current_node_index
                previous_node_index = current_node_index
            self.dpoly[kpol+1] = {'connect':lseg, 'z':z1, 'name':''}
            # last segment
            #ns = self.add_segment(previous_node_index, starting_node_index, name='WALL', z=z1)

    def importosm(self, **kwargs):
        """ import layout from osm file or osmapi

        Parameters
        ----------

        fileosm : string
        address : string
            address to be geocoded
        latlon : tuple
            (latitude,longitude) degrees
        dist_m : float
            distance in meter from the geocoded address (def 200 m )
        cart : boolean
            conversion in cartesian coordinates

        Notes
        -----

        The best and recommended manner to edit a layout is to use the
        josm editor in association with the piclayer plugin.
        This plugin allows to place a geo-adjusted image in the background
        which is very convenient for editing floorplan of buildings.

        In josm editor, nodes are numbered with negative indexes, while in
        pylayers they have a positive index.

        See Also
        --------

        pylayers.gis.osmparser.osmparse

        """

        self._fileosm = kwargs.pop('fileosm','')
        cart = kwargs.pop('cart',False)

        #
        #  zceil ansd zfloor are obtained from actual data
        #
        #  indoor default (0,3)
        #  outdoor default (0,3000)

        #if self.typ=='indoor':
        self.zceil = -1e10
        self.zfloor = 1e10

        if self._fileosm == '':  # by using osmapi address or latlon
            logger.info('load from osmapi')

            self.typ = kwargs.pop('typ','indoor')
            address = kwargs.pop('address','Rennes')
            latlon = kwargs.pop('latlon',0)

            if type(latlon) == 'str':
                 latlon = eval(latlon)

            dist_m = kwargs.pop('dist_m',200)

            coords, nodes, ways, m , latlon , dpoly = osm.getosm(address = address,
                                                          latlon = latlon,
                                                          dist_m = dist_m,
                                                          bcart = cart,
                                                          typ = self.typ)

            self.typ = 'outdoor'
            if cart:
                self.format = 'cart'
            else:
                self.format = 'latlon'

            if latlon == '0':
                self._filename = kwargs['address'].replace(' ', '_') + '.lay'
            else:
                lat = latlon[0]
                lon = latlon[1]

                self._filename = 'lat_' + \
                    str(lat).replace('.', '_') + '_lon_' + \
                    str(lon).replace('.', '_') + '.ini'
        else:  # by reading an osm file

            logger.info('load from osm file')
            # The osm file is supposed to be in $PROJECT/struc/osm directory
            fileosm = pyu.getlong(self._fileosm, os.path.join('struc', 'osm'))
            #coords, nodes, ways, relations, m = osm.osmparse(fileosm, typ=self.typ)
            # typ outdoor parse ways.buildings
            # typ indoor parse ways.ways
            # coords, nodes, ways, relations, m = osm.osmparse(fileosm)
            coords, nodes, ways, m , (lat,lon) , dpoly = osm.getosm(cart = cart,
                                                       filename = fileosm,
                                                       typ = self.typ,
                                                       bexcluded = self.bexcluded)
            if cart:
                self.format = 'cart'
            else:
                self.format = 'latlon'

            self._filename = self._fileosm.replace('osm', 'lay')

        self.dpoly = dpoly
        _np = 0  # _ to avoid name conflict with numpy alias
        _ns = 0
        ns = 0
        nss = 0
        # Reading points  (<0 index)

        # Reorganize points coordinates for detecting
        # duplicate nodes
        # duplicate nodes are saved in dict dup

        kp = [k for k in coords.xy]

        x = np.array([ coords.xy[x][0] for x in  kp ])
        y = np.array([ coords.xy[x][1] for x in  kp ])

        ux = np.argsort(x)
        x_prev = -100
        y_prev = -100
        dup = {}  # dictionnary of duplicate nodes
        for u in ux:
            # if node is not already a duplicate
            if x[u] == x_prev:
                # 2 consecutive points with same lon => check lat
                if y[u] == y_prev:
                    # node u is a duplicate
                    # udate dup dictionnary
                    # printu_prev ,k_prev, x_prev,y_prev
                    # print" ",u ,kp[u], x[u],y[u]
                    dup[kp[u]] = k_prev
            else:
                x_prev = x[u]
                y_prev = y[u]
                u_prev = u
                k_prev = kp[u]

        for npt in coords.xy:
            # if node is not duplicated add node
            if npt not in dup:
                self.Gs.add_node(npt)
                self.Gs.pos[npt] = tuple(coords.xy[npt])
                _np += 1

        # Reading segments
        #
        # ways of osm
        #
        for k, nseg in enumerate(ways.way):
            tahe = ways.way[nseg].refs
            for l in range(len(tahe) - 1):
                nta = tahe[l]
                nhe = tahe[l + 1]
                #
                # if a node is duplicated recover the original node
                #
                if nta in dup:
                    nta = dup[nta]
                if nhe in dup:
                    nhe = dup[nhe]

                d = ways.way[nseg].tags

                #
                # Convert string to integer if possible
                #

                for key in d:
                    try:
                        d[key] = eval(d[key])
                    except:
                        pass


                # getting segment slab information
                if 'slab' in d:
                        slab = d['slab']
                else:  # the default slab name is WALL
                        slab = "WALL"

                if 'z' in d:
                    z = d['z']
                else:
                    if self.typ == 'indoor':
                        z = (0, 3)
                    if self.typ == 'outdoor':
                        z = (0, 3000)
                if type(z[0])==str:
                    zmin = eval(z[0])
                else:
                    zmin = z[0]
                if type(z[1])==str:
                    zmax = eval(z[1])
                else:
                    zmax = z[1]

                #zmin = z[0]
                #zmax = z[1]

                if zmin < self.zfloor:
                    self.zfloor = zmin
                if zmax > self.zceil:
                    self.zceil = zmax

                if 'offset' in d:
                    offset = d['offset']
                else:
                    offset = 0
                #
                # get the common neighbor of nta and nhe if it exists
                #
                #v1.1 u1 = np.array(nx.neighbors(self.Gs, nta))
                #v1.1 u2 = np.array(nx.neighbors(self.Gs, nhe))
                # import ipdb
                # u1 = np.array(self.Gs.node[nta])
                # u2 = np.array(self.Gs.node[nhe])
                # inter_u1_u2 = np.intersect1d(u1, u2)
                #
                # Create  a new segment (iso segments are managed in add_segment)
                #
                ns = self.add_segment(nta, nhe, name=slab, z=z, offset=offset)

        self.Np = _np
        #self.Ns = _ns
        self.Nss = nss
        #
        #
        lon = array([self.Gs.pos[k][0] for k in self.Gs.pos])
        lat = array([self.Gs.pos[k][1] for k in self.Gs.pos])

        # bd = [lon.min(), lat.min(), lon.max(), lat.max()]
        # lon_0 = (bd[0] + bd[2]) / 2.
        # lat_0 = (bd[1] + bd[3]) / 2.

        # self.m = Basemap(llcrnrlon=bd[0], llcrnrlat=bd[1],
        #                  urcrnrlon=bd[2], urcrnrlat=bd[3],
        #                  resolution='i', projection='cass', lon_0=lon_0, lat_0=lat_0)


        self.m = m
        self.extent = (m.lonmin,m.lonmax,m.latmin,m.latmax)
        self.pll = self.m(self.extent[0],self.extent[2])
        self.pur = self.m(self.extent[1],self.extent[3])
        self.extent_c = (self.pll[0],self.pur[0],self.pll[1],self.pur[1])

        #
        # TODO use conversion function
        #
        if (cart and (self.format != 'cart')):
             x, y = self.m(lon, lat)
             self.Gs.pos = {k: (x[i], y[i]) for i, k in enumerate(self.Gs.pos)}
             self.format  = 'cart'

        # del coords
        # del nodes
        # del ways
        # del relations

        #
        # get slab and materials DataBase
        #
        # 1) create material database
        # 2) load materials database
        # 3) create slabs database
        # 4) add materials database to slab database
        # 5) load slabs database
        mat = sb.MatDB()
        mat.load(self._filematini)
        self.sl = sb.SlabDB()
        self.sl.mat = mat
        self.sl.load(self._fileslabini)

        #
        # update self.name with existing slabs database entries
        #
        for k in self.sl.keys():
            if k not in self.name:
                self.name[k] = []

        # convert graph Gs to numpy arrays for speed up post processing
        self.g2npy()

        #
        # add boundary
        #

        self.boundary()

        # save ini file
        self.save()

        #

    def exportosm(self):
        """  export layout in osm file format

        Parameters
        ----------

        _filename : string

        Notes
        -----

        See Also
        --------

        layout.loadosm
        layout.loadini
        layout.check

        """
        # export Layout in osm format
        # The osm filename basenam is the same as the _filename ini file

        _filename, ext = os.path.splitext(self._filename)
        filename = pyu.getlong(_filename + '.osm', 'struc/osm')

        if os.path.exists(filename):
            filename = pyu.getlong(_filename + '_.osm', 'struc/osm')

        fd = open(filename, "w")

        fd.write("<?xml version='1.0' encoding='UTF-8'?>\n")
        fd.write("<osm version='0.6' upload='false' generator='PyLayers'>\n")

        # creating points
        for n in self.Gs.pos:
            if n < 0:
                if n not in self.lboundary:
                    if self.format == 'latlon':
                        lon, lat = self.Gs.pos[n]
                    if self.format == 'cart':
                        x, y = self.Gs.pos[n]
                        lon, lat = self.m(x, y, inverse=True)
                    fd.write("<node id='" + str(n) + "' action='modify' visible='true' lat='" +
                             str(lat) + "' lon='" + str(lon) + "' />\n")

        for n in self.Gs.pos:
            if n > 0:
                #
                # Conditions pour ajout segments
                #
                # _AIR are not added
                #
                # outdoor AIR wall above buildings are not added
                # cond1 is wrong

                cond1 = (self.Gs.node[n]['name'] != '_AIR')
                cond2 = (self.Gs.node[n]['name'] == 'AIR')
                cond3 = (self.Gs.node[n]['z'][1] == self.zceil)
                cond4 = (self.Gs.node[n]['z'][0] == self.zfloor)
                cond5 = (cond2 and cond3)
                cond6 = (cond2 and cond4)
                cond7 = (cond2 and cond3 and cond4) 


                if (cond1 and (not cond5) and (not cond6)) or cond7: 
                    #v1.1 neigh = nx.neighbors(self.Gs, n)
                    neigh = self.Gs[n].keys() 
                    d = self.Gs.node[n]
                    #
                    noden = -10000000 - n
                    fd.write("<way id='" + str(noden) +
                             "' action='modify' visible='true'>\n")
                    fd.write("<nd ref='" + str(neigh[0]) + "' />\n")
                    fd.write("<nd ref='" + str(neigh[1]) + "' />\n")
                    fd.write("<tag k='name' v='" + str(d['name']) + "' />\n")
                    fd.write("<tag k='z' v=\"" + str(d['z']) + "\" />\n")
                    fd.write("<tag k='transition' v='" +
                             str(d['transition']) + "' />\n")
                    fd.write("</way>\n")

        fd.write("</osm>\n")
        fd.close()

    def savepickle(self):
        """ save graph in pickle format
        """
        pass

    def save(self):
        """ save Layout structure in a .lay file

        """

        current_version = 1.4

        if os.path.splitext(self._filename)[1]=='.ini':
            self._filename = self._filename.replace('.ini','.lay')
        #
        # version 1.3 : suppression of index in slab and materials
        #

        config = ConfigParser.RawConfigParser()
        config.optionxform = str
        config.add_section("info")
        config.add_section("points")
        config.add_section("segments")

        if hasattr(self,'dpoly'):
            config.add_section("polygons")

        config.add_section("files")
        config.add_section("slabs")
        config.add_section("materials")

        if self.format == 'latlon':
            config.set("info", "format", "latlon")
        else:
            config.set("info", "format", "cart")

        config.set("info", "version", current_version)
        config.set("info", "type", self.typ)

        if self.typ == 'indoor':
            config.add_section("indoor")
            config.set("indoor", "zceil", self.zceil)
            config.set("indoor", "zfloor", self.zfloor)

        if self.typ == 'outdoor':
            config.add_section("outdoor")

        #
        # save bounding box in latlon for reconstruction of self.m
        #
        if hasattr(self,"m"):
            config.add_section("latlon")
            config.set("latlon","llcrnrlon",self.m.llcrnrlon)
            config.set("latlon","llcrnrlat",self.m.llcrnrlat)
            config.set("latlon","urcrnrlon",self.m.urcrnrlon)
            config.set("latlon","urcrnrlat",self.m.urcrnrlat)
            config.set("latlon","projection",self.m.projection)

        # config.set("info",'Nsegments',self.Ns)
        # config.set("info",'Nsubsegments',self.Nss)

        #for k in self.display:
        #    config.set("display", k, self.display[k])

        # iterate on points
        # boundary nodes and air walls are not saved
        for n in self.Gs.pos:
            if n < 0:
                if n not in self.lboundary:
                    config.set("points", str(n), (self.Gs.pos[n][0], self.Gs.pos[n][1]))

        # iterate on segments
        for n in self.Gs.pos:
            if n > 0:
                cond1 = (self.Gs.node[n]['name'] != '_AIR')
                cond2 = (self.Gs.node[n]['name'] == 'AIR')
                cond3 = (self.Gs.node[n]['z'][1] == self.zceil)
                cond4 = (self.Gs.node[n]['z'][0] == self.zfloor)
                cond5 = (cond2 and cond3)
                cond6 = (cond2 and cond4)
                cond7 = (cond2 and cond3 and cond4)
                #
                # _AIR are not stored  (cond1)
                # AIR segment reaching zceil are not stored  (cond4)
                # AIR segment reaching zfloor are not stored (cond5)
                #
                if (cond1 and (not cond5) and (not cond6)) or cond7:
                    d = copy.deepcopy(self.Gs.node[n])
                    # v1.1 d['connect'] = nx.neighbors(self.Gs, n)
                    d['connect'] = list(self.Gs[n].keys())
                    try:
                        if d['transition']:
                            pass
                    except:
                        d['transition'] = False
                        try:
                            if 'DOOR' in d['ss_name']:
                                d['transition'] = True
                        except:
                            pass
                    # remove normal information from the strucure
                    try:
                        d.pop('norm')
                    except:
                        pass
                    # remove iso information from the strucure
                    try:
                        d.pop('iso')
                    except:
                        pass
                    # remove ncycles information from the strucure
                    try:
                        d.pop('ncycles')
                    except:
                        pass

                    # transition are saved only if True
                    if not d['transition']:
                        d.pop('transition')

                    # offset are saved only if not zero
                    if 'offset' in d:
                        if d['offset']==0:
                            d.pop('offset')


                    config.set("segments", str(n), d)
        #
        # [polygon]
        #   1 = { 'connect':[1,2,3,4], 'z':(100,115),name:''}
        #
        if hasattr(self,'dpoly'):
            for k in self.dpoly:
                config.set("polygons", str(k), self.dpoly[k])

        #
        # [ slabs ]
        #
        # get the list of used slabs

        lslab = [x for x in self.name if len(self.name[x]) > 0]
        lmat = []

        #
        # In case an osm file has been read; there is no .sl
        # By default all the available slabs and materials are provided
        #

        if not hasattr(self,'sl'):
            self.sl = sb.SlabDB(filemat='matDB.ini', fileslab='slabDB.ini')

        for s in lslab:
            ds = {}
            if s not in self.sl:
                if s not in self.sl.mat:
                    self.sl.mat.add(name=s,cval=6,sigma=0,typ='epsr')
                self.sl.add(s,[s],[0.1])

            #ds['index'] = self.sl[s]['index']
            ds['color'] = self.sl[s]['color']
            ds['lmatname'] = self.sl[s]['lmatname']
            for m in ds['lmatname']:
                if m not in lmat:
                    lmat.append(m)
            ds['lthick'] = self.sl[s]['lthick']
            ds['linewidth'] = self.sl[s]['linewidth']
            config.set("slabs", s, ds)

        if "_AIR" not in lslab:
            air = {'color': 'white', 'linewidth': 1,
                   'lthick': [0.1], 'lmatname': ['AIR']}
            config.set("slabs", "_AIR", air)

        if "AIR" not in lslab:
            air = {'color': 'white', 'linewidth': 1,
                   'lthick': [0.1], 'lmatname': ['AIR']}
            config.set("slabs", "AIR", air)

        if "CEIL" not in lslab:
            ceil = {'color': 'grey20', 'linewidth': 1,
                    'lthick': [0.1], 'lmatname': ['REINFORCED_CONCRETE']}
            config.set("slabs", "CEIL", ceil)

        if "FLOOR" not in lslab:
            floor = {'color': 'grey40', 'linewidth': 1,
                     'lthick': [0.1], 'lmatname': ['REINFORCED_CONCRETE']}
            config.set("slabs", "FLOOR", floor)

        #
        # [ materials ] 
        #
        for m in lmat:
            dm = self.sl.mat[m]
            try:
                dm.pop('name')
            except:
                pass
            # store UIT format only if it is used
            if 'a' in dm:
                if dm['a'] ==None:
                    dm.pop('a')
                    dm.pop('b')
                    dm.pop('c')
                    dm.pop('d')
            config.set("materials", m, dm)

        if "REINFORCED_CONCRETE" not in lmat:
            reic = {'mur': (
                1 + 0j), 'epr': (8.69999980927 + 0j), 'roughness': 0.0, 'sigma': 3.0}
            config.set("materials", "REINFORCED_CONCRETE", reic)
        # config.set("files",'materials',self.filematini)
        # config.set("files",'slab',self.fileslabini)
        #
        # [ furniture ]
        #
        config.set("files", 'furniture', self._filefur)
        #
        # handling olf format ( to be removed later) 
        #
        if os.path.splitext(self._filename)[1]=='.ini':
            fileout = self._filename.replace('.ini','.lay')
        else:
            fileout = self._filename

        filelay = pyu.getlong(fileout, pro.pstruc['DIRLAY'])
        print(filelay)
        fd = open(filelay, "w")
        config.write(fd)
        fd.close()
        # convert graph Gs to numpy arrays for speed up post processing
        # ideally an edited Layout should be locked while not saved.
        # self.g2npy()
        self._hash = hashlib.md5(open(filelay, 'rb').read()).hexdigest()

    def load_fast(self):
        """ load a layout from a .lay file

        The filename is in self._filename

        Format version 1.3
        ------------------

        [info]
        format = {cart | latlon}
        version =
        type = {indoor | outdoor}

        [points]
        -1 = (x,y)

        [segments]
        1 = {'slab':'',transition:boolean,'connect:[-1,-2],'z':(0,3)}

        [slabs]
        WALL = {'lthick':[,],'lmat':[,],'color:'','linewidth':float}

        [materials]
        BRICK = {'mur':complex,'epsr':complex,'sigma':float,'roughness':}

        [polygons]
        1 = {'connect':[1,2,3,4],'name':NAME,'z':(zmin,zmax)}

        [indoor]
        zceil =
        zfloor =

        [latlon]
        llcrnrlon = 
        llcrnrlat = 
        urcrnrlon =
        urcrnrlat = 
        projection = 



        """
        logger.info('load lay file in load fast')
        # di : dictionnary which reflects the content of ini file
        di = {}
        config = ConfigParser.RawConfigParser()
        config.optionxform = str
        filelay = pyu.getlong(self._filename, pro.pstruc['DIRLAY'])
        config.read(filelay)


        sections = config.sections()
        for section in sections:
            di[section] = {}
            options = config.options(section)
            for option in options:
                try:
                    di[section][option] = config.get(section, option)
                except:
                    print(section, option)

        self.Np = len(di['points'])
        self.Ns = len(di['segments'])

        #
        # Ng : number of polygons
        # polygons introduced in 1.4 format
        #
        if config.has_section('polygons'):
            self.Ng = len(di['polygons'])

        self.Gs = nx.Graph(name='Gs')
        self.Gs.pos = {}
        self.labels = {}

        #
        # [info]
        #    format     {cart,latlon}
        #    version    int
        #    type       {'indoor','outdoor'}
        if 'version' in di['info']:
            self.version = di['info']['version']

        if 'type' in di['info']:
            self.typ = di['info']['type']

        self.name = {}

        if ((self.typ!='indoor') &
            (self.typ!='outdoor') &
            (self.typ!='floorplan')):
            print("invalid file type in ",self._filename)
            return(None)

        #
        # [indoor]
        #   zceil
        #   zfloor
        #
        if self.typ == 'indoor':
            self.zceil = eval(di['indoor']['zceil'])
            self.zfloor = eval(di['indoor']['zfloor'])

        # old format
        if self.typ == 'floorplan':
            self.zceil = eval(di['floorplan']['zceil'])
            self.zfloor = eval(di['floorplan']['zfloor'])

        # from format 1.3 floorplan is call indoor
        if self.typ=='floorplan':
            self.typ = 'indoor'
        #
        # [outdoor]
        #   TODO add a DEM file
        #
        if self.typ == 'outdoor':
            if 'outdoor' in di:
                if 'zceil' in di['outdoor']:
                    self.zceil = eval(di['outdoor']['zceil'])
                else:
                    self.zceil  =  3000    # upper limit for AIR wall

                if 'zfloor' in di['outdoor']:
                    self.zfloor = eval(di['outdoor']['zfloor'])
                else:
                    self.zfloor = 0
            else:
                self.zfloor = 0
                self.zceil  =  3000    # upper limit for AIR walls

        #
        #
        # manage ini file with latlon coordinates
        #
        # if the format is latlon, coordinates are converted into
        # cartesian coordinates with the coords.cartesian method
        #

        logger.info('load latlon coordinates and convert')
        if 'format' in di['info']:
            if di['info']['format'] == 'latlon':
                or_coord_format = 'latlon'
                coords = osm.Coords()
                coords.clean()
                coords.latlon = {eval(i): np.array( eval(di['points'][i])) for i in di['points']}
                coords.boundary = np.hstack((np.min(np.vstack(coords.latlon.values()), axis=0),
                                             np.max(np.vstack(coords.latlon.values()), axis=0)))
                assert(len(coords.boundary)==4)
                coords.cartesian(cart=True)
            else:
                or_coord_format = 'cart'
        else:
            or_coord_format = 'cart'

        #
        # update display section
        #
        logger.info('load_fast : update display section')
        if 'display' in di:
            for k in di['display']:
                try:
                    self.display[k] = eval(di['display'][k])
                except:
                    self.display[k] = di['display'][k]

        # self.ax = self.display['box']
        #
        # [points]
        #
        # update points section

        logger.info('load_fast : reading points')
        lnodeindex = [ eval(x) for x in di['points']]
        self.Gs.add_nodes_from(lnodeindex)  # add point node
        if self.format=='cart':
            self.Gs.pos = { eval(i) : eval(di['points'][i]) for i in di['points']}
        else:
            self.Gs.pos = { i : coords.xy[i] for i in coords.xy}

        #    #
        #    # limitation of point precision is important for avoiding
        #    # topological problems in shapely.
        #    # Layout precision is hard limited to millimeter precision.
        #    #

        #self.Gs.pos[nodeindex] = (
        #    round(1000 * x) / 1000., round(1000 * y) / 1000.)
        #self.labels[nodeindex] = nn

        #
        # [segments]
        #
        # update segments section

        self.name['AIR'] = []
        self.name['_AIR'] = []
        #
        # get the segment maximum index
        #
        if len(di['segments'])>0:
            maxnum = max([eval(x) for x in di['segments'].keys()])
            logger.info('load_fast : reading segments')

            lsegnum = np.sort([ eval(x) for x in di['segments'].keys() ])
            maxnum = lsegnum[-1]

            for k in lsegnum:
                d = eval(di['segments'][str(k)])
                nta = d['connect'][0]
                nhe = d['connect'][1]
                name = d['name']

                z = d['z']

                if not 'transition' in d:
                    transition = False
                else:
                    transition = d['transition']

                if not 'offset' in d:
                    offset = 0
                else:
                    offset = d['offset']

                # add new segment
                #
                # The segment number is the same as in the .lay file
                #
                # Very useful feature
                #

                if self.typ=='indoor':
                    boutdoor = False
                else:
                    boutdoor = True

                num = self.add_segment(nta, nhe,
                                       num = k,
                                       name = name,
                                       transition = transition,
                                       offset = offset,
                                       z = z,
                                       maxnum = maxnum,
                                       boutdoor = boutdoor)

            # only for indoor
            # exploit iso for segment completion (AIR type)
            #  Complement single segment which do not reach zceil or zfloor with
            #  an iso segment with AIR property
            #

            if (self.typ == 'indoor'):
                segdone = []
                logger.info('segments completion')
                for key in di['segments']:
                    iseg = eval(key)
                    d = eval(di['segments'][key])
                    nta = d['connect'][0]
                    nhe = d['connect'][1]
                    # if not already done
                    if iseg not in segdone:
                        # get all the iso from the segment key
                        iso = copy.copy(self.Gs.node[iseg]['iso'])
                        # append key to iso
                        iso.append(iseg)
                        # stack all the intervals in increasing order
                        ziso = []
                        for ns in iso:
                            ziso.append(self.Gs.node[ns]['z'])

                        # get the complementary intervals
                        if self.typ == 'outdoor':
                            zmin = 1e6
                            zmax = -1e6
                            for iz in ziso:
                                zmin = np.minimum(zmin,min(iz))
                                zmax = np.maximum(zmax,max(iz))
                            ziso = [(zmin,zmax)]

                        # pyutil compint (get complementary interval)
                        zair = pyu.compint(ziso,self.zfloor,self.zceil)

                        # add AIR wall in the intervals
                        for za in zair:
                            num = self.add_segment(nta, nhe,
                                    name='AIR',
                                    offset=0,
                                    z=(za[0], za[1]))
                        segdone = segdone + iso

        #
        # add _AIR wall around the layout
        # This is not the place to add the boundary because
        # the boundary may depend on the position of tx and rx
        # boundary is call at the end of layout init 

        self.boundary(bg2npy = False)

        # load poygons
        if config.has_section('polygons'):
            logger.info("reading polygons")
            self.dpoly = { k : eval(di['polygons'][k]) for k in di['polygons']}
            #
            # TODO convert keys into int

        # compliant with config file without  material/slab information
        #
        # {latlon]
        #
        if config.has_section('latlon'):
            llcrnrlon = eval(config.get('latlon', 'llcrnrlon'))
            llcrnrlat = eval(config.get('latlon', 'llcrnrlat'))
            urcrnrlon = eval(config.get('latlon', 'urcrnrlon'))
            urcrnrlat = eval(config.get('latlon', 'urcrnrlat'))
            projection = config.get('latlon','projection')
            lon_0 = (llcrnrlon+urcrnrlon)/2.
            lat_0 = (llcrnrlat+urcrnrlat)/2.

            # Construction of Basemap for coordinates transformation
            self.m = Basemap(llcrnrlon=llcrnrlon,
                             llcrnrlat=llcrnrlat,
                             urcrnrlon=urcrnrlon,
                             urcrnrlat=urcrnrlat,
                             resolution='i',
                             projection=projection,
                             lon_0=lon_0,
                             lat_0=lat_0)

            self.extent = (llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat)
            self.pll = self.m(self.extent[0],self.extent[2])
            self.pur = self.m(self.extent[1],self.extent[3])
            self.extent_c = (self.pll[0],self.pur[0],self.pll[1],self.pur[1])


        if config.has_section('files'):
            # self.filematini=config.get('files','materials')
            # self.fileslabini=config.get('files','slab')
            self.filefur = config.get('files', 'furniture')

        if config.has_section('slabs'):
            #filemat = self._filename.replace('ini', 'mat')
            #fileslab = self._filename.replace('ini', 'slab')

            ds = di['slabs']
            dm = di['materials']

            for k in ds:
                ds[k] = eval(ds[k])
            for k in dm:
                dm[k] = eval(dm[k])

            self.sl = sb.SlabDB(ds=ds, dm=dm)

        # In this section we handle the ini file format evolution

        if 'fileoverlay' in self.display:
            self.display['overlay_file'] = self.display.pop('fileoverlay')
            self.display['overlay_axis'] = self.display['box']
            self.save()

        if 'inverse' in self.display:
            self.display['overlay_flip'] = ""
            self.display.pop('inverse')
            self.save()

        # convert graph Gs to numpy arrays for faster post processing
        logger.info('graph 2 numpy : g2npy')
        self.g2npy()

        # This hash should be saved with the gpickle file

        logger.info('save the hash')
        dirname = self._filename.replace('.lay','')
        path = os.path.join(pro.basename, 'struc', 'gpickle', dirname)
        if not os.path.exists(path):
            os.mkdir(path)
        fd = open(os.path.join(path,'.hash'),'w')
        fd.write(self._hash)
        fd.write('\n'+str(self.segboundary))
        fd.write('\n'+self.typ)
        fd.close()
        logger.info('dump a pickle of Gs')
        self.lbltg = ['s']
        self.dumpw()

    def load(self):
        """ load a layout from a .lay file

        The filename is in self._filename

        Format version 1.4
        ------------------

        [info]
        format = {cart | latlon}
        version =
        type = {indoor | outdoor}

        [points]
        -1 = (x,y)

        [segments]
        1 = {'slab':'',transition:boolean,'connect:[-1,-2],'z':(0,3)}

        [slabs]
        WALL = {'lthick':[,],'lmat':[,],'color:'','linewidth':float}

        [materials]
        BRICK = {'mur':complex,'epsr':complex,'sigma':float,'roughness':}

        [polygons]
        1 = {'connect':[1,2,3,4],'name':NAME,'z':(zmin,zmax)}

        [indoor]
        zceil =
        zfloor =

        [latlon]
        llcrnrlon =
        llcrnrlat =
        urcrnrlon =
        urcrnrlat =
        projection =


        """
        logger.info('load lay file')
        # di : dictionnary which reflects the content of ini file
        di = {}
        config = ConfigParser.RawConfigParser()
        config.optionxform = str
        filelay = pyu.getlong(self._filename, pro.pstruc['DIRLAY'])
        config.read(filelay)


        sections = config.sections()
        for section in sections:
            di[section] = {}
            options = config.options(section)
            for option in options:
                try:
                    di[section][option] = config.get(section, option)
                except:
                    print(section, option)

        self.Np = len(di['points'])
        self.Ns = len(di['segments'])
        #
        # Ng : number of polygons
        # polygons introduced in 1.4 format
        #
        if config.has_section('polygons'):
            self.Ng = len(di['polygons'])

        self.Gs = nx.Graph(name='Gs')
        self.Gs.pos = {}
        self.labels = {}

        #
        # [info]
        #    format     {cart,latlon}
        #    version    int
        #    type       {'indoor','outdoor'}
        if 'version' in di['info']:
            self.version = di['info']['version']

        if 'type' in di['info']:
            self.typ = di['info']['type']

        # floorplan is no longer a valid typ
        self.name = {}
        if ((self.typ!='indoor') 
                &   (self.typ!='outdoor') ):
            print("invalid file type in ",self._filename)
            return(None)

        #
        # [indoor]
        #   zceil
        #   zfloor
        #
        if self.typ == 'indoor':
            self.zceil = eval(di['indoor']['zceil'])
            self.zfloor = eval(di['indoor']['zfloor'])

        # [outdoor]
        #   TODO add a DEM file
        #
        if self.typ == 'outdoor':
            if 'outdoor' in di:
                if 'zceil' in di['outdoor']:
                    self.zceil = eval(di['outdoor']['zceil'])
                else:
                    self.zceil  =  3000    # upper limit for AIR wall

                if 'zfloor' in di['outdoor']:
                    self.zfloor = eval(di['outdoor']['zfloor'])
                else:
                    self.zfloor = 0
            else:
                self.zfloor = 0
                self.zceil  =  3000    # upper limit for AIR walls

        #
        #
        # manage ini file with latlon coordinates
        #
        # if the format is latlon, coordinates are converted into
        # cartesian coordinates with the coords.cartesian method
        #

        logger.info('load coordinates')
        if 'format' in di['info']:
            if di['info']['format'] == 'latlon':
                or_coord_format = 'latlon'
                coords = osm.Coords()
                coords.clean()
                coords.latlon = {i: np.array( eval(di['points'][i])) for i in di['points']}
                coords.boundary = np.hstack((np.min(np.array(coords.latlon.values()), axis=0),
                                             np.max(np.array(coords.latlon.values()), axis=0)))
                coords.cartesian(cart=True)
            else:
                or_coord_format = 'cart'
        else:
            or_coord_format = 'cart'

        #
        # update display section
        #

        logger.info('update display section')
        if 'display' in di:
            for k in di['display']:
                try:
                    self.display[k] = eval(di['display'][k])
                except:
                    self.display[k] = di['display'][k]

        # self.ax = self.display['box']
        #
        # [points]
        #
        # update points section
        logger.info('reading points')
        for nn in di['points']:
            nodeindex = eval(nn)
            if or_coord_format == 'latlon':
                x, y = coords.xy[nn]
            else:
                x, y = eval(di['points'][nn])

            #
            # limitation of point precision is important for avoiding
            # topological problems in shapely.
            # Layout precision is hard limited to millimeter precision.
            #

            self.Gs.add_node(nodeindex)  # add point node
            self.Gs.pos[nodeindex] = ( round(1000 * x) / 1000., round(1000 * y) / 1000.)
            self.labels[nodeindex] = nn

        #
        # [segments]
        #
        # update segments section

        self.name['AIR'] = []
        self.name['_AIR'] = []
        #
        # get the maximum index
        #
        maxnum = max([eval(x) for x in di['segments'].keys()])
        logger.info('reading segments')
        for key in di['segments']:

            d = eval(di['segments'][key])
            nta = d['connect'][0]
            nhe = d['connect'][1]
            #print(key,nta,nhe)


            name = d['name']
            z = d['z']

            if not 'transition' in d:
                transition = False
            else:
                transition = d['transition']

            if not 'offset' in d:
                offset = 0
            else:
                offset = d['offset']

            # add new segment
            #
            # The segment number is the same as in the .lay file
            #
            # Very useful feature
            #
            num = self.add_segment(nta, nhe,
                                   num = eval(key),
                                   name = name,
                                   transition = transition,
                                   offset = offset,
                                   z = z)

        # exploit iso for segment completion (AIR type)
        #
        #  Complement single segment which do not reach zceil or zfloor with
        #  an iso segment with AIR property
        #
        segdone = []
        logger.info('segments completion')
        for key in di['segments']:
            iseg = eval(key)
            d = eval(di['segments'][key])
            nta = d['connect'][0]
            nhe = d['connect'][1]
            # if not already done
            if iseg not in segdone:
                # get all the iso from the segment key
                iso = copy.copy(self.Gs.node[iseg]['iso'])
                # append key to iso
                iso.append(iseg)
                # stack all the intervals in increasing order
                ziso = []
                for ns in iso:
                    ziso.append(self.Gs.node[ns]['z'])
                # get the complementary intervals
                if self.typ == 'outdoor':
                    zmin = 1e6
                    zmax = -1e6
                    for iz in ziso:
                        zmin = np.minimum(zmin,min(iz))
                        zmax = np.maximum(zmax,max(iz))
                    ziso = [(zmin,zmax)]

                # pyutil compint (get complementary interval)
                zair = pyu.compint(ziso,self.zfloor,self.zceil)
                # add AIR wall in the intervals
                for za in zair:
                    num = self.add_segment(nta, nhe,
                            name='AIR',
                            offset=0,
                            z=(za[0], za[1]))
                segdone = segdone + iso

        #
        # add _AIR wall around the layout
        #
        self.boundary()
        if config.has_section('polygons'):
            logger.info("reading polygons")
            self.dpoly = di['polygons']

        # compliant with config file without  material/slab information
        #
        # {latlon]
        #
        if config.has_section('latlon'):
            llcrnrlon = eval(config.get('latlon', 'llcrnrlon'))
            llcrnrlat = eval(config.get('latlon', 'llcrnrlat'))
            urcrnrlon = eval(config.get('latlon', 'urcrnrlon'))
            urcrnrlat = eval(config.get('latlon', 'urcrnrlat'))
            projection = config.get('latlon','projection')
            lon_0 = (llcrnrlon+urcrnrlon)/2.
            lat_0 = (llcrnrlat+urcrnrlat)/2.

            # Construction of Basemap for coordinates transformation
            self.m = Basemap(llcrnrlon=llcrnrlon,
                             llcrnrlat=llcrnrlat,
                             urcrnrlon=urcrnrlon,
                             urcrnrlat=urcrnrlat,
                             resolution='i',
                             projection=projection,
                             lon_0=lon_0,
                             lat_0=lat_0)

            self.extent = (llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat)
            self.pll = self.m(self.extent[0],self.extent[2])
            self.pur = self.m(self.extent[1],self.extent[3])
            self.extent_c = (self.pll[0],self.pur[0],self.pll[1],self.pur[1])


        if config.has_section('files'):
            # self.filematini=config.get('files','materials')
            # self.fileslabini=config.get('files','slab')
            self._filefur = config.get('files', 'furniture')

        if config.has_section('slabs'):
            #filemat = self._filename.replace('ini', 'mat')
            #fileslab = self._filename.replace('ini', 'slab')

            ds = di['slabs']
            dm = di['materials']

            for k in ds:
                ds[k] = eval(ds[k])
            for k in dm:
                dm[k] = eval(dm[k])

            self.sl = sb.SlabDB(ds=ds, dm=dm)

        # In this section we handle the ini file format evolution

        if 'fileoverlay' in self.display:
            self.display['overlay_file'] = self.display.pop('fileoverlay')
            self.display['overlay_axis'] = self.display['box']
            self.save()

        if 'inverse' in self.display:
            self.display['overlay_flip'] = ""
            self.display.pop('inverse')
            self.save()

        # convert graph Gs to numpy arrays for faster post processing
        logger.info('g2npy')
        self.g2npy()
        #
        fd = open(filelay,'rb')
        self._hash = hashlib.md5(fd.read()).hexdigest()
        fd.close()

    def loadfur(self, _filefur):
        """ loadfur load a furniture file

        Parameters
        ----------

        _filefur  : string
            short name of the furniture ini file

        Notes
        -----

            Furniture objects are stored in self.lfur list

        Examples
        --------

        Load a Layout file and an associated furniture ini file



        .. plot::
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> from pylayers.gis.layout import *
            >>> L = Layout('WHERE1.lay')
            >>> L.loadfur('Furw1.ini')
            >>> fig = plt.figure()
            >>> ax = fig.gca()
            >>> fig,ax = L.showGs(fig=fig,ax=ax,furniture=True)
            >>> ti = plt.title('loadfur')
            >>> plt.show()


        """
        filefur = pyu.getlong(_filefur, pro.pstruc['DIRFUR'])
        config = ConfigParser.ConfigParser()
        config.read(filefur)
        furname = config.sections()
        self.lfur = []

        for name in furname:
            F = fur.Furniture()
            F.load(_filefur, name)
            self.lfur.append(F)
        self.filefur = _filefur

    def load_modif(self, _filename, build=True, cartesian=False, dist_m=400):
        """ load a Layout in different formats

        Parameters
        ----------

        _filename : string

        Notes
        -----
        +  .lay   : ini file format (natural one) DIRLAY


        """

        newfile = False
        filename = pyu.getlong(_filename, pro.pstruc['DIRLAY'])
        if os.path.exists(filename):  # which exists
            self.loadini(arg)
        else:  # which do not exist
            self._filename = _filename
            newfile = True
            print("new file", self._filename)

        #  construct geomfile (.off) for vizualisation with geomview
        self.subseg()
        if not newfile:
            try:
                self.geomfile()
            except:
                print("problem to construct geomfile")
        # if check:
        #     self.check()

        self.boundary(dx=10, dy=10)

        # create shapely polygons L._shseg

    def subseg(self):
        """ establishes the association : name <->  edgelist

        Returns
        -------

        dico : dict
               sub segment name as key and segment number as value

        """

        dico = {}
        listtransition = []

        for k in self.Gs.nodes.keys():
            dk = self.Gs.nodes[k]
            if 'transition' in dk:
                transition = dk['transition']
                if transition:
                    listtransition.append(k)

            if 'ss_name' in dk:
                lname = dk['ss_name']
                for j, name in enumerate(lname):
                    if name in dico:
                        dico[name].append((k, j))
                    else:
                        dico[name] = [(k, j)]

        self.dsseg = dico
        self.listtransition = listtransition
        return(dico)

    def add_pnod(self, p, e1, e2):
        """ Project point p on segment e1 along segment e2

        Parameters
        ----------

            p  : ndarray
                point
            e1 : int
                edge number 1
            e2 : int
                edge number 2

        ..todo
            This function is void
        """
        #p1 = p + alpha*ve2
        #p1 = pa + beta * (pb-pa)
        pass

    def add_fnod(self, p=(0.0, 0.0)):
        """ add free node  p

        Parameters
        ----------

        p :  (1x2) tuple

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('defstr.lay')
        >>> L.add_fnod((10.0,10.0))
        -13


        """
        # next free node
        if len(self.Gs.nodes)>0:
            num = -( -min(self.Gs.nodes) + 1 )
        else:
            num = -1
        self.Gs.add_node(num)
        self.Gs.pos[num] = p
        self.Np = self.Np + 1
        # update labels
        self.labels[num] = str(num)
        return(num)

    def add_nfpe(self, np0, s1, s2):
        """ Add node on s1 from projection of np0 along s2

        Parameters
        ----------

        np0  : point number
        s1   : edge number 1
        s2   : edge number 2

        """
        np1 = list(self.Gs[s1].keys())
        np2 = list(self.Gs[s2].keys())
        xA = self.Gs.pos[np1[0]][0]
        yA = self.Gs.pos[np1[0]][1]
        xB = self.Gs.pos[np1[1]][0]
        yB = self.Gs.pos[np1[1]][1]
        xC = self.Gs.pos[np2[0]][0]
        yC = self.Gs.pos[np2[0]][1]
        xD = self.Gs.pos[np2[1]][0]
        yD = self.Gs.pos[np2[1]][1]
        xP = self.Gs.pos[np0][0]
        yP = self.Gs.pos[np0][1]

        A = np.array([[xB - xA, xD - xC], [yB - yA, yD - yC]])
        b = np.array([xP - xA, yP - yA])
        x = sp.linalg.solve(A, b)
        if ((x[0] > 0.) & (x[0] < 1.0)):
            self.add_pons(s1, 1 - x[0])

    def add_pons(self, ns, alpha=0.5):
        """ add point on segment

        Parameters
        ----------

        ns  : int
            segment number

        alpha : parameterization of the point
            alpha = 0 (tail) alpha = 1 (head)

        Notes
        -----

        delete segment ns
        create 2 segments with same properties

        """
        # v1.1 nop = self.Gs.neighbors(ns)
        nop = list(self.Gs[ns])
        namens = self.Gs.nodes[ns]['name']
        zminns = self.Gs.nodes[ns]['z'][0]
        zmaxns = self.Gs.nodes[ns]['z'][1]
        p1 = np.array([self.Gs.pos[nop[0]][0], self.Gs.pos[nop[0]][1]])
        p2 = np.array([self.Gs.pos[nop[1]][0], self.Gs.pos[nop[1]][1]])
        p = tuple(alpha * p1 + (1 - alpha) * p2)
        num = self.add_fnod(p)
        # delete old edge ns
        self.del_segment(ns)
        # add new edge np[0] num
        self.add_segment(nop[0], num, name=namens, z=[
                         zminns, zmaxns], offset=0)
        # add new edge num np[1]
        self.add_segment(num, nop[1], name=namens, z=[
                         zminns, zmaxns], offset=0)


    def add_segment(self,n1,n2,**kwargs):
        """  add segment between node n1 and node n2

        Parameters
        ----------

        n1  : integer < 0
        n2  : integer < 0
        num : segment index (-1 default not given)
        maxnum : maximum number (-1 default not given)
        name : string
            layer name 'PARTITION'
        z : tuple of 2 floats
            default = (0,40000000)
        offset : float
            [-1,1] default (0)
        bootdoor : boolean
            if outdoor add an _AIR wall above the segment

        Returns
        -------

        num : segment number (>0)

        Notes
        -----

        A segment dictionnary has the following mandatory attributes

        name : slab name associated with segment
        z : list (zmin,zmax)   (meters)
        norm : array  (1x3)  segment normal
        transition : boolean
        ncycles : list of involved cycles
        connect : list of point number
        iso : list of isosegment

        If a segment is _AIR it cannnot be duplicated

        """

        num = kwargs.pop('num', -1)
        maxnum = kwargs.pop('maxnum', -1)
        transition = kwargs.pop('transition', False)
        name = kwargs.pop('name', 'PARTITION')
        offset = kwargs.pop('offset', 0)
        verbose = kwargs.pop('verbose', True)
        boutdoor = kwargs.pop('boutdoor', False)
        zmax = kwargs.pop('zmax', 3000)
        z = kwargs.pop('z', (0.0, zmax))
        # if 2 points are selected

        if ((n1 < 0) & (n2 < 0) & (n1 != n2)):
            nseg = [s for s in self.Gs.nodes if s > 0]
            if num==-1:
                if len(nseg) > 0:
                    num = max(maxnum+1, max(nseg) + 1)   # index not given 
                else: # first segment index not given
                    num = 1
            else:
                pass # segment index given
        else:
            if verbose:
                print("add_segment : error not a node", n1, n2)
            return

        # transition = False
        if (name == '_AIR'):
            # if name == 'AIR':
            transition = True

        #
        # unit vector along horizontal segment
        #

        p1 = np.array(self.Gs.pos[n1])
        p2 = np.array(self.Gs.pos[n2])
        p2mp1 = p2 - p1
        t = p2mp1 / np.sqrt(np.dot(p2mp1, p2mp1))

        #
        # n = t x z  (2D)
        #

        norm = np.array([t[1], -t[0], 0])

        #
        # Two segments with the same end points are iso segments
        #
        # Determine if there are existing segments with the same neighbors ?

        nbnta = self.Gs[n1].keys()
        nbnhe = self.Gs[n2].keys()

        #v1.1 nbnta = self.Gs.neighbors(n1)
        #nbnhe = self.Gs.neighbors(n2)
        #
        # add a segment node to Gs
        #

        if boutdoor == True:

            self.Gs.add_node(num, name=name,
                         z = z,
                         norm = norm,
                         transition = transition,
                         offset = offset,
                         connect = [n1, n2],
                         iso = [num + maxnum],
                         ncycles = []
                         )

            self.Gs.add_node(num + maxnum, name='AIR',
                         z = (z[1], zmax),
                         norm = norm,
                         transition = True,
                         offset = offset,
                         connect = [n1, n2],
                         iso = [num],
                         ncycles = [])
            self.Gs.pos[num] = tuple((p1 + p2) / 2.)
            self.Gs.pos[num+maxnum] = tuple((p1 + p2) / 2.)
        else:
            #
            # BAD IDEA : Not scalable
            #
            same_seg = list(set(nbnta).intersection(nbnhe))

            #
            # Impossible to have duplicated _AIR
            #
            # Warning : The 3 following lines are very important
            # it breaks buildGt if commented
            # Please do not comment them.
            #

            if (name == '_AIR'):
                if len(same_seg) > 0:
                    return None

            self.Gs.add_node(num, name=name,
                         z = z,
                         norm = norm,
                         transition = transition,
                         offset = offset,
                         connect = [n1, n2],
                         iso = [],
                         ncycles = []
                         )
            self.Gs.pos[num] = tuple((p1 + p2) / 2.)

            #
            # update iso of the 2 segments
            #
            for k in same_seg:
                if num not in self.Gs.nodes[k]['iso']:
                    self.Gs.nodes[k]['iso'].append(num)
                if k not in self.Gs.nodes[num]['iso']:
                    self.Gs.nodes[num]['iso'].append(k)

        #
        # Segment point position is placed at the middle of segment
        #


        #
        # Connectivity between segment node num and points nodes n1 and n2
        #

        self.Gs.add_edge(n1, num)
        self.Gs.add_edge(n2, num)

        if boutdoor:

            self.Gs.add_edge(n1, num + maxnum)
            self.Gs.add_edge(n2, num + maxnum)
        #
        # Update current total number of segments
        #
            self.Ns = self.Ns + 2
            try:
                self.name[name].append(num)
                self.name['AIR'].append(num + maxnum)
            except:
                self.name[name] = [num]
                self.name['AIR'] = [num+ maxnum]
            # update label
            self.labels[num] = str(num)
            self.labels[num + maxnum] = str(num + maxnum)

        else:
            self.Ns = self.Ns + 1
        # update slab name <-> edge number dictionnary
            try:
                self.name[name].append(num)
            except:
                self.name[name] = [num]
            # update label
            self.labels[num] = str(num)

        if name not in self.display['layers']:
            self.display['layers'].append(name)

        # update shseg
        self._shseg.update({num:sh.LineString((self.Gs.pos[n1], self.Gs.pos[n2]))})

        return(num)

    def merge_segment(self,n1,n2):
        """ merge segment n2 included in n1

        Parameters
        ----------

        n1 : int
            segment 1 (the larger) index
        n2 : int
            segment 2 (the smaller) index

        """
        # get height/slabname information from segment n1
        zn1 = self.Gs.nodes[n1]['z']
        namen1 = self.Gs.nodes[n1]['name']
        # get height/slabname information from segment n2
        zn2 = self.Gs.nodes[n2]['z']
        namen2 = self.Gs.nodes[n2]['name']
        if min(zn1)<min(zn2):
            znlow = (min(zn1),min(zn2))
        if max(zn1)>max(zn2):
            znhigh = (max(zn2),max(zn1))

        # get termination points of segment n1 (p1 -- p4)
        conn_n1 = self.Gs.nodes[n1]['connect']
        conn_n2 = self.Gs.nodes[n2]['connect']

        p1_index = conn_n1[0]
        p4_index = conn_n1[1]
        p2_index = conn_n2[0]
        p3_index = conn_n2[1]

        p1 = np.r_[self.Gs.pos[p1_index]]
        p2 = np.r_[self.Gs.pos[p2_index]]
        p3 = np.r_[self.Gs.pos[p3_index]]
        p4 = np.r_[self.Gs.pos[p4_index]]

        # determine point order p1 - p2 - p3 - p4

        v14 = p4 - p1
        v23 = p3 - p2

        if np.dot(v14,v23)<0:
            p2_index, p3_index = p3_index, p2_index
            p2, p3 = p3, p2


        # 1 delete segment n1

        self.del_segment([n1])


        # create new segment p1 - p2

        self.add_segment(p1_index, p2_index, z=zn1, name=namen1)

        # create new segment p3 - p4
        self.add_segment(p3_index, p4_index, z=zn1, name=namen1)

        # create new segment p2 - p3 with complementary heights

        if 'zlow' in locals():
            self.add_segment(p2_index, p3_index, z=znlow, name=namen1)

        if 'zhigh' in locals():
            self.add_segment(p2_index, p3_index, z=znhigh, name=namen1)

    def repair(self,dseg):
        """ repair layout

        Parameters
        ----------

        dseg : dict
            {ns : [np1,np2]}

        Notes
        -----

        Merge the superposed segments which has been determined by the check
        method.

        """
        for nseg in dseg:
            num_p = dseg[nseg]
            if len(num_p)==2:
                ns1 = np.r_[nx.neighbors(self.Gs,num_p[0])]
                ns2 = np.r_[nx.neighbors(self.Gs,num_p[1])]
                ns_inter = np.intersect1d(ns1,ns2)
                for nseg2 in ns_inter:
                    if ((self.Gs.nodes[nseg2]['name']!='AIR')
                        and ((self.Gs.nodes[nseg2]['name']!='_AIR'))):
                        self.merge_segment(nseg,nseg2)

    def wedge2(self, apnt):
        """ calculate wedge angle of a point

        Parameters
        ----------

        lpnt : array int
           list of point number


        """

        if isinstance(apnt, list):
            apnt = np.array(apnt)

        # 0. Find the position of diffraction point
        ptdiff = self.pt[:, self.iupnt[-apnt]]

        # 1. Find the associated segments and positions of a diff points

        #v1.1 aseg = map(lambda x: filter(lambda y: y not in self.name['AIR'],
        #                            nx.neighbors(self.Gs, x)),
        #           apnt)
        aseg = map(lambda x: filter(lambda y: y not in self.name['AIR'],
                                    self.Gs[x].keys()),apnt)
        # manage flat angle : diffraction by flat segment e.g. door limitation)
        [aseg[ix].extend(x) for ix, x in enumerate(aseg) if len(x) == 1]
        # get points positions
        pts = np.array(map(lambda x: self.seg2pts([x[0], x[1]]), aseg))

        pt1 = pts[:, 0:2, 0]  # tail seg1
        ph1 = pts[:, 2:4, 0]  # head seg1
        pt2 = pts[:, 0:2, 1]  # tail seg2
        ph2 = pts[:, 2:4, 1]  # head seg2

        # 2. Make the correct association
        # pts is (nb_diffraction_points x 4 x 2)
        # - The dimension 4 represent the 2x2 points: t1,h1 and t2,h2
        # tail and head of segemnt 1 and 2 respectively
        # a segment
        # - The dimension 2 is x,y
        #
        # The following aims to determine which tails and heads of
        # segments associated to a give diffraction point
        # are connected

        # point diff is pt1
        updpt1 = np.where(np.sum(ptdiff.T == pt1, axis=1) == 2)[0]
        # point diff is ph1
        updph1 = np.where(np.sum(ptdiff.T == ph1, axis=1) == 2)[0]

        # point diff is pt2
        updpt2 = np.where(np.sum(ptdiff.T == pt2, axis=1) == 2)[0]
        # point diff is ph2
        updph2 = np.where(np.sum(ptdiff.T == ph2, axis=1) == 2)[0]

        pa = np.empty((len(apnt), 2))
        pb = np.empty((len(apnt), 2))

        # seg 1 :
        # if pt1 diff point =>  ph1 is the other point
        pa[updpt1] = ph1[updpt1]
        # if ph1 diff point =>  pt1 is the other point
        pa[updph1] = pt1[updph1]
        # seg 2 :
        # if pt2 diff point =>  ph2 is the other point
        pb[updpt2] = ph2[updpt2]
        # if ph2 diff point =>  pt2 is the other point
        pb[updph2] = pt2[updph2]
        # pt is the diffraction point
        pt = ptdiff.T

        vptpa = pt - pa
        vptpan = vptpa.T / np.sqrt(np.sum((vptpa) * (vptpa), axis=1))
        vptpb = pt - pb
        vptpbn = vptpb.T / np.sqrt(np.sum((vptpb) * (vptpb), axis=1))
        v1 = vptpan
        v2 = vptpbn

        ang = geu.vecang(vptpbn, vptpan)
        ang[~uleft] = geu.vecang(vptpan, vptpan)

    def wedge(self, lpnt):
        """ calculate wedge angle of a point

        Parameters
        ----------

        lpnt : list of int
           list of point number


        """

        #v1.1 aseg = map(lambda x: filter(lambda y: y not in
        #                            self.name['AIR'],
        #                            nx.neighbors(self.Gs, x)),
        #           lpnt)
        aseg = map(lambda x: filter(lambda y: y not in
                                    self.name['AIR'],
                                    self.Gs[x]), lpnt)

        pts = np.array(map(lambda x: self.seg2pts([x[0], x[1]]).reshape(4, 2), aseg))
        #map(lambda x: pt ,pts)
        N = np.shape(pts)[0]
        sector = []
        for k in range(N):
            pt1 = pts[k, 0:2, 0]
            ph1 = pts[k, 2:4, 0]
            pt2 = pts[k, 0:2, 1]
            ph2 = pts[k, 2:4, 1]
            if (pt1 == pt2).all():
                pa = ph1
                pb = ph2
                pt = pt1
                ang = geu.sector(pa, pb, pt)
            if (pt1 == ph2).all():
                pa = ph1
                pb = pt2
                pt = pt1
                ang = geu.sector(pa, pb, pt)
            if (ph1 == pt2).all():
                pa = pt1
                pb = ph2
                pt = ph1
                ang = geu.sector(pa, pb, pt)
            if (ph1 == ph2).all():
                pa = pt1
                pb = pt2
                pt = ph1
                ang = geu.sector(pa, pb, pt)

            sector.append(ang)

        return(sector)

    def add_furniture(self, name='R1_C', matname='PARTITION', origin=(0., 0.),
                      zmin=0., height=0., width=0., length=0., angle=0.):
        """  add piece of furniture

        Parameters
        ----------

        name : string
            default = 'R1_C'
        matname : string
            default = 'PARTITION'
        origin : tuple of floats
        height : float
            default = 0
        width : float
            default = 0
        length : float
            default = 0
        angle : float
            default = 0

        """

        # compute the four points
        p0 = origin
        u = np.array([np.cos(angle * np.pi / 180),
                      np.sin(angle * np.pi / 180)])
        v = np.array([-np.sin(angle * np.pi / 180),
                      np.cos(angle * np.pi / 180)])
        p1 = p0 + u * length
        p2 = p1 + v * width
        p3 = p2 - u * length
        # adding free nodes
        n0 = self.add_fnod(p0)
        n1 = self.add_fnod(p1)
        n2 = self.add_fnod(p2)
        n3 = self.add_fnod(p3)
        # adding segments
        self.add_segment(n0, n1, name=matname, z=(zmin, zmin + height))
        self.add_segment(n1, n2, name=matname, z=(zmin, zmin + height))
        self.add_segment(n2, n3, name=matname, z=(zmin, zmin + height))
        self.add_segment(n3, n0, name=matname, z=(zmin, zmin + height))

    def add_furniture_file(self, _filefur, typ=''):
        """  add pieces of furniture from .ini files

        Parameters
        ----------

        _filefur : string

        """

        filefur = pyu.getlong(_filefur, pro.pstruc['DIRFUR'])
        config = ConfigParser.ConfigParser()
        config.read(filefur)
        furname = config.sections()
        for fur in furname:
            name = config.get(fur, "name")
            matname = config.get(fur, "matname")
            origin = tuple(ast.literal_eval(config.get(fur, "origin")))
            height = config.getfloat(fur, "height")
            width = config.getfloat(fur, "width")
            length = config.getfloat(fur, "length")
            angle = config.getfloat(fur, "angle")
            thickness = config.getfloat(fur, "thickness")
            # ~ if matname=='WOOD':
            # ~ zmin = height
            # ~ height=thickness
            # ~ else:
            # ~ zmin=0.0
            # .. todo: be more generic relate to floor level
            zmin = 0.0
            if typ == '':
                self.add_furniture(name, matname, origin,
                                   zmin, height, width, length, angle)
            else:
                try:
                    self.add_furniture(name, matname, origin,
                                       zmin, height, width, length, angle)
                except:
                    raise NameError('No such furniture type - ' + typ + '-')

    def del_points(self, lp):
        """ delete points in list lp

        Parameters
        ----------

        lp : list
            node list

        """

        # test if array
        if (type(lp) == np.ndarray):
            ln = list(ln)

        # test if list
        if (type(lp) != list):
            lp = [lp]

        #print("lp : ", lp)
        # get segments involved in points list
        ls = self.nd2seg(lp)

        #print("ls : ", ls)
        # 1) delete involved segments
        for k in ls:
            assert(k > 0)
            self.del_segment(k)
            #print('del ', k)
        # 2) delete involved points
        for n1 in lp:
            assert(n1 < 0)
            # v1.1 nbrs = self.Gs.neighbors(n1)
            nbrs = self.Gs[n1].keys()
            self.Gs.remove_node(n1)
            del self.Gs.pos[n1]
            self.labels.pop(n1)
            self.Np = self.Np - 1
        # 3) updating structures
        self.g2npy()

    def del_segment(self, le, verbose=True, g2npy=True):
        """ delete segments in le

        Parameters
        ----------

        le : list of segments number

        See Also
        --------

        pylayers.gis.layout.Layout.del_node

        Notes
        -----

        100% of time is in g2npy

        """
        if (type(le) == np.ndarray):
            le = list(le)

        if (type(le) != list):
            le = [le]

        for e in le:
            assert(e > 0)
            name = self.Gs.nodes[e]['name']
            iso = self.Gs.nodes[e]['iso']
            [self.Gs.nodes[i]['iso'].remove(e) for i in iso
             if e in self.Gs.nodes[i]['iso']]
            del self.Gs.pos[e]  # delete edge position
            self.Gs.remove_node(e)
            if e in self.labels:
                self.labels.pop(e)
            self.Ns = self.Ns - 1
            # update slab name <-> edge number dictionnary
            self.name[name].remove(e)
            # delete iso if required

            try:
                # remove shapely seg
                self._shseg.pop(e)
            except:
                pass
        if g2npy:
            self.g2npy()

    def point_touches_seg(self,pt,lseg=[],segtol=1e-2,tahetol=1e-2):
        """ determine if a point is touching a segment

            Parameters
            ----------

            pt : a point (2,)
            seg : a list of segments to test.
                 if [] => all Gs segments are tested

            segdtol : distance tolerance point to segment
            tahetol : distance tolerance point to segment extremeties
                        => a point on segment extremeties is considered
                           not touching the segseg



            Return
            ------

            ltseg : lsit of touched segments (by the point)

        """

        if lseg == []:
            lseg = self.Gs.nodes()

        ltseg = []
        allnodes = self.Gs.nodes()
        for s in lseg :
            if s > 0 and s in allnodes:
                n0,n1  = self.Gs.nodes[s]['connect']
                dta,dhe,h = geu.dptseg(np.array(pt)[:,None],
                            np.array(self.Gs.pos[n0])[:,None],
                            np.array(self.Gs.pos[n1])[:,None])
                if (h <= segtol) and ((dta > tahetol) and (dhe > tahetol)):
                    ltseg.append(s)
        return ltseg



    def seg_intersection(self,**kwargs):
        ''' determine if a segment intersects any other segment of the layout

            Parameters
            ----------

            shLine : a shapely LineString
            or
            ta,he : tail/head of a segment

            Returns
            -------

            llay_seg : list of layout's segments intersected
            lshP : list of shapely points of intersections.

            See Also
            --------

            editor.py

        '''

        if ('ta' in kwargs) and ('he' in kwargs):
            seg = sh.LineString((kwargs['ta'],kwargs['he']))
        elif 'shLine' in kwargs:
            seg = kwargs['shLine']

        # WARNING : use crosses instead of interesects
        # otherwise 2 segments connected to a same node
        # are considered as intersecting

        binter = [seg.crosses(x) for x in list(self._shseg.values())]
        if np.sum(binter) > 0:
            uinter = np.where(binter)[0]
            llay_seg = []
            lshP = []
            for k in uinter:
                # layout segment
                llay_seg.append(list(self._shseg.keys())[k])
                lay_shseg = self._shseg[llay_seg[-1]]
                # intersection shapely point
                lshP.append(seg.intersection(lay_shseg))

            return(llay_seg,lshP)
        else:
            return ([],[])


    def mask(self):
        """  returns the polygonal mask of the building

        Returns
        -------

        mask : geu.Polygon

        Notes
        -----

        This function assumes graph Gt has been generated

        """
        if hasattr(self,Gt):
            # takes the 1st cycle polygon
            p = self.Gt.node[1]['polyg']
            # get the exterior of the polygon
            ps = sh.Polygon(p.exterior)
            # make the union of the exterior of all the cycles
            #
            # cycle : -1 exterior
            #          0 ??
            #

            for k in self.Gt.node:
                if (k != 0) & (k != -1):
                    p = self.Gt.node[k]['polyg']
                    ps = ps.union(sh.Polygon(p.exterior))

            mask = geu.Polygon(ps)
            mask.setvnodes(self)
            return(mask)
        else:
            print("Gt not built")

    def translate(self, vec):
        """ translate layout

        Parameters
        ----------

     loa   vec :

        """
        for k in self.Gs.pos:
            pt = self.Gs.pos[k]
            self.Gs.pos[k] = (pt[0] + vec[0], pt[1] + vec[1])

    def rotate(self, angle=90):
        """ rotate the layout

        Parameters
        ----------

        angle : float
            (degrees)

        """

        a = angle * np.pi / 180

        for k in self.Gs.pos:
            pt = self.Gs.pos[k]
            ptr = np.dot(
                array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]), array(pt))
            self.Gs.pos[k] = (ptr[0], ptr[1])

        self.g2npy()

    def check2(self):
        """ Layout checking

        Returns
        -------

        tseg ; list of segment shapely

        """

        tseg = []

        for k in list(self.Gs.node.keys()):
            if k > 0:
                #v1.1 lnp = self.Gs.neighbors(k)
                lnp = list(self.Gs[k].keys())
                p1 = self.Gs.pos[lnp[0]]
                p2 = self.Gs.pos[lnp[1]]
                tseg.append(sh.LineString([(p1[0], p1[1]), (p2[0], p2[1])]))

        N = len(tseg)
        for k in combinations(range(N), 2):
            seg1 = tseg[k[0]]
            seg2 = tseg[k[1]]
            if seg1.crosses(seg2):
                print("crosses :", k[0], k[1])
            if seg1.contains(seg2):
                print("contains :", k[0], k[1])
            if seg2.contains(seg1):
                print("contains :", k[0], k[1])
            if seg1.overlaps(seg2):
                print("overlaps :", k[0], k[1])
            if seg2.overlaps(seg1):
                print("overlaps :", k[0], k[1])

        return(tseg)

    def cleanup(self):
        """ cleanup the Layout

        Notes
        -----

        1. Remove nodes which are not connected
        2. Remove supperimposed segments
        """
        lk = list(self.Gs.nodes.keys())
        for n in lk:
            if ((n < 0) & (self.Gs.degree(n) == 0)):
                self.Gs.remove_node(n)
                del self.Gs.pos[n]
                try:
                    self.Gv.remove_node(n)
                except:
                    pass
        self.Np = len(np.nonzero(np.array(list(self.Gs.nodes.keys())) < 0)[0])
        
        aseg_conn=[]
        for seg in self.Gs.nodes():
            if seg >0:
                n0,n1 = list(nx.neighbors(self.Gs,seg))
                aseg_conn.append([seg,n0,n1])
        aseg_conn = np.array(aseg_conn)

        # aseg_conn=np.array([[list(nx.neighbors(self.Gs,x))] for x in self.Gs.nodes() if x >0])
        uni,upos=np.unique(aseg_conn[:,1:],axis=0,return_index=True)
        utbd = [x for x in range(len(aseg_conn)) if not x in upos] 
        tbd = aseg_conn[utbd,0]
        for k in tbd:
            self.del_segment(k)

        self.g2npy()

    def info_segment(self, s1):
        """ information about segment

        Parameters
        ----------

        s1 : segment number

        """
        # v1.1 nebd = self.Gs.neighbors(s1)
        nebd = self.Gs[s1].keys()
        n1 = nebd[0]
        n2 = nebd[1]
        #v1.1 nns1 = self.Gs.neighbors(n1)
        #nns2 = self.Gs.neighbors(n2)
        nns1 = self.Gs[n1].keys()
        nns2 = self.Gs[n2].keys()
        ds1 = self.Gs.nodes[s1]
        print(n1, ' : ', nns1)
        print(n2, ' : ', nns2)
        print('------------')
        print('Slab     : ', ds1['name'])
        print('zmin (m) : ', ds1['z'][0])
        print('zmax (m) : ', ds1['z'][1])
        try:
            print('------------')
            a = ds1['ss_name']
            print('subseg Slabs  : ', ds1['ss_name'])
            print('subseg (zmin,zmax) (m) : ', ds1['ss_z'])
        except:
            pass

    def edit_seg(self, e1, data={}):
        """ edit segment

        Parameters
        ----------

        e1 : integer
            edge number
        data : dict
            dictionnary of value of seg or subseg

        Notes
        -----

        A segment has the following properties :
            + name  : string
            + z  :  tuple
            + transition : boolean (default FALSE)
            + offset : [-1,1]
        If a segment has subsegments attached the following properties are
        added :
            + ss_name : list of string
            + ss_z : list of subsegment e.q. [(min height (meters),max height (meters))]
            + ss_offset : list of offset in [0,1]
        """

        if data == {}:
            pass
        else:

            ename = self.Gs.nodes[e1]['name']
            # manage self.name
            self.name[ename].pop(self.name[ename].index(e1))
            # manage self.display['name']
            if len(self.name[ename]) == 0:
                try:
                    self.display['layers'].pop(
                        self.display['layers'].index(ename))
                except:
                    pass

            for k in data:
                self.Gs.nodes[e1][k] = data[k]

        if data['name'] in self.name:
            self.name[data['name']].append(e1)
        else:
            self.name[data['name']]=[e1]

        if data['name'] not in self.display['layers']:
            self.display['layers'].append(data['name'])

        return data

    def have_subseg(self, e1):
        """ check if edge e1 have subseg

        Parameters
        ----------

            e1 : int

        Returns
        -------

            have_subseg_bool : boolean


        """
        dk = self.Gs.nodes[e1]
        if len(dk['iso'])>0:
            return True
        else:
            return False

    def find_edgelist(self, edgelist, nodelist):
        """
        edgelist = find_edgelist(edgelist,nodelist)

        edgelist : input edgelist
        nodelist : input nodelist

        return the subset of edgelist

        Not Finished :

        """
        tail = self.tahe[0, edgelist]
        head = self.tahe[1, edgelist]

        nt = np.intersect1d_nu[tail, nodelist]
        nh = np.intersect1d_nu[head, nodelist]

        edgelist = edgelist[np.unique(ed_t, ed_h)]
        return(edgelist)

    def diag(self, p1, p2, l, al1, al2, quadsel=0):
        """ return edge list from a diagonal zone

        Parameters
        -----------

        p1  : np.array
        p2  : np.array
        tol :
        al1 :
        al2 :
        quadsel : 0   all quadrant
              2 1
              3 4

        Returns
        -------

        edgelist

        """
        x = self.pt[0, :]
        y = self.pt[1, :]

        #
        # selection du quadran
        #
        if (quadsel == 0):
            u0 = np.arange(self.Np)
        if (quadsel == 1):
            u0 = np.nonzero((y > p1[1]) & (x > p1[0]))[0]
        if (quadsel == 2):
            u0 = np.nonzero((y > p1[1]) & (x <= p1[0]))[0]
        if (quadsel == 3):
            u0 = np.nonzero((y <= p1[1]) & (x <= p1[0]))[0]
        if (quadsel == 4):
            u0 = np.nonzero((y <= p1[1]) & (x > p1[0]))[0]

        x_u0 = x[u0]
        y_u0 = y[u0]

        #
        # Permutation points
        #

        if (p1[0] > p2[0]):
            pt = p2
            p2 = p1
            p1 = pt
        #
        # Box length
        #

        Dx = p2[0] - p1[0]
        Dy = p2[1] - p1[1]

        L = np.sqrt(Dx ** 2 + Dy ** 2)

        #
        # p1 p2
        #
        if ((abs(Dx) > finfo(float).eps) & (abs(Dy) > finfo(float).eps)):
            a = Dy / Dx
            b = p1[1] - a * p1[0]
            b1 = p1[1] + p1[0] / a
            b2 = p2[1] + p2[0] / a

            delta_b = tol * L / abs(Dx)
            delta_b1 = al1 * L * L / abs(Dy)
            delta_b2 = al2 * L * L / abs(Dy)

            u1 = np.nonzero(y_u0 < a * x_u0 + b + delta_b / 2.)[0]
            x_u1 = x_u0[u1]
            y_u1 = y_u0[u1]
            u2 = np.nonzero(y_u1 > a * x_u1 + b - delta_b / 2.)[0]
            x_u2 = x_u1[u2]
            y_u2 = y_u1[u2]
            if (a > 0):
                u3 = np.nonzero(y_u2 > -x_u2 / a + b1 - delta_b1)[0]
                x_u3 = x_u2[u3]
                y_u3 = y_u2[u3]
                u4 = np.nonzero(y_u3 < -x_u3 / a + b2 + delta_b2)[0]
            else:
                u3 = np.nonzero(y_u2 < -x_u2 / a + b1 + delta_b1)[0]
                x_u3 = x_u2[u3]
                y_u3 = y_u2[u3]
                u4 = np.nonzero(y_u3 > -x_u3 / a + b2 - delta_b2)[0]
                x_u4 = x_u3[u4]
                y_u4 = y_u3[u4]
#
# p1 p2 vertical
#
        if (abs(Dx) <= finfo(float).eps):
            u1 = np.nonzero(x < p1[0] + tol / 2.)[0]
            x_u1 = x[u1]
            y_u1 = y[u1]
            u2 = np.nonzero(x_u1 > p1[0] - tol / 2.)[0]
            y_u2 = y[u2]
            if (p1[1] > p2[1]):
                u3 = np.nonzero(y_u2 < p1[1] + al1 * L)[0]
                y_u3 = y[u3]
                u4 = np.nonzero(y_u3 > p2[1] - al2 * L)[0]
            else:
                u3 = np.nonzero(y_u2 < p2[1] + al2 * L)[0]
                y_u3 = y[u3]
                u4 = np.nonzero(y_u3 > p1[1] - al1 * L)[0]
#
# p1 p2 horizontal
#
        if (abs(Dy) <= finfo(float).eps):
            u1 = np.nonzero(y < p1[1] + tol / 2.)[0]
            y_u1 = y[u1]
            u2 = np.nonzero(y_u1 > p1[1] - tol / 2.)[0]
            x_u2 = x[u2]
            if (p1(1) > p2(1)):
                u3 = np.nonzero(x_u2 < p1[0] + al1 * L)[0]
                x_u3 = x[u3]
                u4 = np.nonzero(x_u3 > p2[0] - al2 * L)[0]
            else:
                u3 = np.nonzero(x_u2 < p2[0] + al2 * L)[0]
                x_u3 = x[u3]
                u4 = np.nonzero(x > p1[0] - al1 * L)[0]

        nodelist = u0[u1[u2[u3[u4]]]]
        edgelist = np.arange(self.Ns)
        edgelist = self.find_edge_list(edgelist, nodelist)
        return(edgelist)

    def nd2seg(self, ndlist):
        """ convert node list to edge list

        Parameters
        ----------

        ndlist : list or ndarray
            node list

        Returns
        -------

        seglist : ndarray
            edge list


        Notes
        -----

        previously nd2ed

        """
        if isinstance(ndlist, np.ndarray):
            ndlist = ndlist.tolist()

        seglist = []

        # for n in ndlist:
        #    seglist = seglist + self.Gs.adj[n].keys()
        #l = map(lambda x: self.Gs.adj[x].keys(), ndlist)
        l = [ list(dict(self.Gs.adj[x]).keys()) for x in ndlist ]
        seglist = []
        for y in l :
            seglist = seglist + y
        #reduce(lambda x, y: x + y, l)

        return(np.unique(np.array(seglist)))

    def ed2nd(self, edlist):
        """ convert edgelist to nodelist

        Parameters
        ----------
        edlist : list or ndarray
            edge list

        Returns
        -------
        ndlist : ndarray
            node list

        """
        if isinstance(edlist, np.ndarray):
            edlist = edlist.tolist()
            # mecanisme de concatenation de listes
        ndlist = []
        for e in edlist:
            ndlist = ndlist + self.Gs.adj[e].keys()

        return(np.unique(ndlist))

    def get_zone(self, ax):
        """ get point list and segment list in a rectangular zone

        Parameters
        ----------
        ax  : list ot tuple
            [xmin,xmax,ymin,ymax]

        Returns
        -------
        ptlist,seglist

        """

        xmin = ax[0]
        xmax = ax[1]
        ymin = ax[2]
        ymax = ax[3]

        ptlist = []
        for n in self.Gs.nodes.keys():
            if n < 0:
                x = self.Gs.pos[n][0]
                y = self.Gs.pos[n][1]
                if ((x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)):
                    ptlist.append(n)

        seglist = self.nd2seg(ptlist)

        return ptlist, seglist

    def get_points(self, boxorpol , tol = 0.05):
        """ get points list and segments list in a polygonal zone

        Parameters
        ----------

        boxorpol  : list or tuple
            [xmin,xmax,ymin,ymax]
              or shapely Polygon

        Returns
        -------

        (pt,ke) : points coordinates and index
            pt : (2xn)
            ke : (,n)

        Notes
        -----

        This method returns all the existing Layout points inside a box zone or
        the boundary of a polygon

        """

        if type(boxorpol) == geu.Polygon:
            N = len(boxorpol.vnodes)/2
            eax  = boxorpol.bounds
            xmin = eax[0] - tol
            xmax = eax[2] + tol
            ymin = eax[1] - tol
            ymax = eax[3] + tol
        else:
            xmin = boxorpol[0]
            xmax = boxorpol[1]
            ymin = boxorpol[2]
            ymax = boxorpol[3]

        #
        # layout points
        #

        x = self.pt[0,:]
        y = self.pt[1,:]

        uxmin = (x>= xmin)
        uymin = (y>= ymin)
        uxmax = (x<= xmax)
        uymax = (y<= ymax)

        #
        # k True when all conditons are True simultaneously
        #

        k = np.where(uxmin*uymin*uxmax*uymax==1)[0]
        #pt = np.array(zip(x[k],y[k])).T
        # pt (2 x N )
        pt = np.vstack((x[k],y[k]))
        ke = self.upnt[k]

        # if(pt.shape[1]<N):
        #     plt.ion()
        #     fig,a=self.showG('s')
        #     a.plot(pt[0,:],pt[1,:],'or')
        #     a.plot(eax[0],eax[1],'or')
        #     plt.show()
        # ux = ((x>=xmin).all() and (x<=xmax).all())
        # uy = ((y>=ymin).all() and (y<=ymax).all())

        return((pt,ke))

    def angleonlink3(self, p1=np.array([0, 0, 1]), p2=np.array([10, 3, 1])):
        """ returns (seglist,angle) in retangular area defined by p1 and p2

        Parameters
        ----------

        p1 : np.array (3 x N) or (3,)
        p2 : np.array (3 x N) or (3,)

        Returns
        -------

        data : structured array x N
            'i' : index
            's' : slab
            'a' : angle (in radians)


        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('DLR2.lay')
        >>> p1 = np.array([0,0,1])
        >>> p2 = np.array([10,3,2])
        >>> data = L.angleonlink3(p1,p2)

        #array([(0, 141, 1.2793395519256592), (0, 62, 0.29145678877830505),
               (0, 65, 0.29145678877830505)],
              dtype=[('i', '<i8'), ('s', '<i8'), ('a', '<f4')])

        See Also
        --------

        antprop.loss.Losst
        geomutil.intersect3

        """

        sh1 = np.shape(p1)
        sh2 = np.shape(p2)

        assert sh1[0] == 3, pdb.set_trace()
        assert sh2[0] == 3, pdb.set_trace()

        #  p1 (1 axe) and p2 (2 or 3 axes) 
        if (len(sh1) < 2) & (len(sh2) > 1):
            p1 = np.outer(p1, np.ones(sh2[1]))

        #  p2 (1 axe) and p1 (2 or 3 axes) 
        if (len(sh2) < 2) & (len(sh1) > 1):
            p2 = np.outer(p2, np.ones(sh1[1]))

        if (len(sh2) < 2) & (len(sh1) < 2):
            p1 = np.outer(p1, np.ones(1))
            p2 = np.outer(p2, np.ones(1))

        # 3 x N
        u = p1 - p2
        # 1 x N
        nu = np.sqrt(np.sum(u * u, axis=0))
        # 3 x N
        un = u / nu[np.newaxis, :]

        #
        # warning : seglist contains the segment number in tahe not in Gs
        #
        seglist  =  self.seginframe2(p1[0:2], p2[0:2])
        useglist = np.unique(seglist)
        #seglist  = np.unique(self.seginframe(p1[0:2], p2[0:2]))

        upos = np.nonzero(useglist >= 0)[0]
        uneg = np.nonzero(useglist < 0)[0]


        # nNLOS = len(uneg) + 1
        # # retrieve the number of segments per link
        # if nNLOS > 1:
        #      llink = np.hstack(
        #         (uneg[0], np.hstack((uneg[1:], array([len(seglist)]))) - uneg - 1))
        # else:
        #     llink = np.array([len(seglist)])
        # [(link id,number of seg),...]
        # nl = zip(np.arange(nlink),llink)n

        useglist = useglist[upos]

        npta = self.tahe[0, useglist]
        nphe = self.tahe[1, useglist]

        Pta = self.pt[:, npta]
        Phe = self.pt[:, nphe]

        Nscreen = len(npta)
        # get segment height bounds
        zmin = np.array([self.Gs.nodes[x]['z'][0]
                         for x in self.tsg[useglist]])
        zmax = np.array([self.Gs.nodes[x]['z'][1]
                         for x in self.tsg[useglist]])
        # centroid of the screen
        Pg = np.vstack(((Phe + Pta) / 2., (zmax + zmin) / 2.))
        Ptahe = Phe - Pta
        L1 = np.sqrt(np.sum(Ptahe * Ptahe, axis=0))
        # 3 x Nscreen U1 is in plane xy
        U1 = np.vstack((Ptahe / L1, np.zeros(Nscreen)))
        L2 = zmax - zmin
        U2 = np.array([0, 0, 1])[:, None]  # 3 x 1  U2 is along z

        #
        # p1 : 3 x Ng
        # p2 : 3 x Ng
        # Pg : 3 x Nscreen
        # U1 : 3 x Nscreen
        # U2 : 3 x 1
        # L1 : ,Nscreen
        # L2 : ,Nscreen

        bo, pt = geu.intersect3(p1, p2, Pg, U1, U2, L1, L2)
        ubo = np.where(bo)

        Nseg = len(ubo[0])
        data = np.zeros(Nseg, dtype=[('i', 'i8'), ('s', 'i8'), ('a', np.float32)])

        data['i'] = ubo[0]
        data['s'] = self.tsg[useglist[ubo[1]]]

        #
        # Calculate angle of incidence refered from segment normal
        #

        norm = self.normal[:, useglist[ubo[1]]]
        # vector along the link
        uu = un[:, ubo[0]]
        unn = abs(np.sum(uu * norm, axis=0))
        angle = np.arccos(unn)

        data['a'] = angle

        return(data)

    def angleonlink(self, p1=np.array([0, 0]), p2=np.array([10, 3])):
        """ angleonlink(self,p1,p2) return (seglist,angle) between p1 and p2

        Parameters
        ----------

        p1 : np.array (2 x Np) or (2,)
        p2 : np.array (2 x Np) or (2,)

        Returns
        -------

        data['i']
        data['s'] : list of segment number
        data['a'] : angle (in radians) between segment and LOS axis

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('DLR.lay')
        >>> p1 = np.array([0,0])
        >>> p2 = np.array([10,3])
        >>> alpha = L.angleonlink(p1,p2)

        #array([(0, 141, 1.2793395519256592), (0, 62, 0.29145678877830505),
               (0, 65, 0.29145678877830505)],
              dtype=[('i', '<i8'), ('s', '<i8'), ('a', '<f4')])


        """

        sh1 = np.shape(p1)
        sh2 = np.shape(p2)
        assert sh1[0] == 2
        assert sh2[0] == 2

        if (len(sh1) < 2) & (len(sh2) > 1):
            p1 = np.outer(p1, np.ones(sh2[1]))

        if (len(sh2) < 2) & (len(sh1) > 1):
            p2 = np.outer(p2, np.ones(sh1[1]))

        if (len(sh2) < 2) & (len(sh1) < 2):
            p1 = np.outer(p1, np.ones(1))
            p2 = np.outer(p2, np.ones(1))

        # 2 x N
        u = p1 - p2
        # 1 x N
        nu = np.sqrt(np.sum(u * u, axis=0))
        # 2 x N
        un = u / nu[np.newaxis, :]

        #seglist = self.seginframe2(p1, p2)
        seglist = self.seginframe(p1, p2)
        upos = np.nonzero(seglist >= 0)[0]
        uneg = np.nonzero(seglist < 0)[0]

        nNLOS = len(uneg) + 1
        # retrieve the number of segments per link
        if nNLOS > 1:
            llink = np.hstack(
                (uneg[0], np.hstack((uneg[1:], array([len(seglist)]))) - uneg - 1))
        else:
            llink = np.array([len(seglist)])

        # llink : list of link length

        npta = self.tahe[0, seglist[upos]]
        nphe = self.tahe[1, seglist[upos]]

        Pta = self.pt[:, npta]
        Phe = self.pt[:, nphe]


        #
        # This part should possibly be improved
        #

        for i, nl in enumerate(llink):
            try:
                P1 = np.hstack((P1, np.outer(p1[:, i], np.ones(nl))))
                P2 = np.hstack((P2, np.outer(p2[:, i], np.ones(nl))))
                ilink = np.hstack(
                    (ilink, array([-1]), i * np.ones(nl, dtype='int')))
            except:
                P1 = np.outer(p1[:, i], np.ones(nl))
                P2 = np.outer(p2[:, i], np.ones(nl))
                ilink = i * np.ones(nl, dtype='int')

        bo = geu.intersect(P1, P2, Pta, Phe)

        upos_intersect = upos[bo]

        seglist2 = seglist[upos_intersect]

        idxlnk = ilink[upos_intersect]
        #
        # Calculate angle of incidence refered from segment normal
        #

        norm = self.normal[0:2, seglist2]
        # vector along the linkco
        uu = un[:,idxlnk]
        unn = abs(np.sum(uu * norm, axis=0))
        angle = np.arccos(unn)

        # seglist = seglist+1
        seglist = np.array([self.tsg[x] for x in  seglist2])
        data = np.zeros(len(seglist), dtype=[
                        ('i', 'i8'), ('s', 'i8'), ('a', np.float32)])

        #
        # update subsegment in seglist
        #
        # self.lsss

        data['i'] = idxlnk
        data['s'] = seglist
        data['a'] = angle

        return data

    def angleonlinkold(self, p1=np.array([0, 0]), p2=np.array([10, 3])):
        """ angleonlink(self,p1,p2) returns seglist between p1 and p2

        Parameters
        ----------

        p1 : (1 x 2 )
            [0,0]
        p2 : (1 x 2 )
            [10,3]

        Returns
        -------

        seglist : list
                  list of segment number on the link
        theta

        Examples
        --------

        #>>> from pylayers.gis.layout import *
        #>>> L = Layout('DLR.lay','matDB.ini','slabDB.ini')
        #>>> p1 = np.array([0,0])
        #>>> p2 = np.array([10,3])
        #>>> L.angleonlinkold(p1,p2)
        #(array([59, 62, 65]), array([ 1.27933953,  0.29145679,  0.29145679]))


        """

        logger.warning('This function is deprecated use')

        u = p1 - p2
        nu = np.sqrt(np.dot(u, u))
        un = u / nu

        seglist = self.seginframe(p1, p2)
        # new implementation of seginframe is faster
        #
        #seglist = self.seginframe2(p1, p2)

        npta = self.tahe[0, seglist]
        nphe = self.tahe[1, seglist]

        Pta = self.pt[:, npta]
        Phe = self.pt[:, nphe]

        P1 = np.outer(p1, np.ones(len(seglist)))
        P2 = np.outer(p2, np.ones(len(seglist)))

        bo = geu.intersect(P1, P2, Pta, Phe)

        seglist = seglist[bo]

        #
        # Calculate normal angle angle of incidence
        #
        tail = self.tahe[0, seglist]
        head = self.tahe[1, seglist]

        vn = np.vstack((self.pt[1, head] - self.pt[1, tail],
                        self.pt[0, head] - self.pt[0, tail]))
        mvn = np.outer(np.ones(2), np.sqrt(np.sum(vn * vn, axis=0)))

        n = vn / mvn
        uu = np.outer(un, np.ones(len(seglist)))
        unn = abs(np.sum(uu * n, axis=0))
        theta = np.arccos(unn)

        # printvn
        # printmvn
        # print'n :',n
        # print'un : ',unn
        # print'theta (deg)',the*180./pi

        # seglist = seglist+1
        seglist = np.array([self.tsg[x] for x in seglist])

        return(seglist, theta)

    def layeronlink(self, p1, p2):
        """

        layeronlink(self,p1,p2) return seglist between p1 and p2

        p1 : (1 x 2 )
        p2 : (1 x 2 )
        """
        seglist = self.seginframe(p1, p2)
        npta = self.tahe[0, seglist]
        nphe = self.tahe[1, seglist]

        Pta = self.pt[:, npta]
        Phe = self.pt[:, nphe]

        P1 = np.outer(p1, np.ones(len(seglist)))
        P2 = np.outer(p2, np.ones(len(seglist)))

        bool = np.intersect(P1, P2, Pta, Phe)

        seglist = seglist[bool]

        return seglist

    def seguv(self, iseg):
        """ returns unitary vector along segments

        Parameters
        ----------

        iseg : np.array
                index of segments

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('DLR.lay')
        >>> idx = np.array([1,2,3,17])
        >>> v1 = L.seguv(idx)
        >>> idx = np.array([1])
        >>> v2= L.seguv(idx)

        """
        # idx : npt
        idx = self.tgs[iseg]
        # tahe : 2 x npt
        tahe = self.tahe[:, idx]
        if len(iseg) > 1:
            ta = tahe[0, :]
            he = tahe[1, :]
        else:
            ta = tahe[0]
            he = tahe[1]
        pta = self.pt[:, ta]
        phe = self.pt[:, he]
        # v  : 2 x npt
        v = pta - phe
        # mv : npt
        mv = np.sqrt(np.sum(v * v, axis=0))
        # vn : 2 x npt
        if len(idx) > 1:
            vn = v / mv[np.newaxis, :]
        else:
            vn = (v / mv).reshape(2)
        return(vn)

    def seg2pts(self, aseg):
        """ convert segments array from Gs numerotation

        to corresponding termination points array in pt

        Parameters
        ----------

        aseg : np.array (,Ns) or int for single value:w
            array of segment number (>0)

        Returns
        -------

        pth : np.array (4 x Ns)
            pth is a vstacking of tail point (2,Ns) and head point (2,Ns)

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> import numpy as np
        >>> L = Layout('defstr.lay')
        >>> aseg = np.array([1,3,6])
        >>> pt =  L.seg2pts(aseg)

        Notes
        -----

        surprisingly self.s2pc is slower than this function

        """

        if not isinstance(aseg, np.ndarray):
            aseg = np.array([aseg])

        assert(len(np.where(aseg < 0)[0]) == 0)
        utahe = self.tgs[aseg]
        #if (utahe>=0).all():
        tahe = self.tahe[:, utahe]
        ptail = self.pt[:, tahe[0, :]]
        phead = self.pt[:, tahe[1, :]]
        pth = np.vstack((ptail, phead))
        pth = pth.reshape(pth.shape[0], pth.shape[-1])
        return pth
        #else:
        #    pdb.set_trace()

    def segpt(self, ptlist=np.array([0])):
        """ return the seg list of a sequence of point number

        Parameters
        ----------

        ptlist array(1xNp)
            point number array

        Returns
        -------

        seglist
            array seglist associated with ptlist

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('TA-Office.lay')
        >>> ptlist  = np.array([0,1])
        >>> seg = L.segpt(ptlist)

        Notes
        -----

        """
        seglist = np.array([], dtype=int)
        for i in ptlist:
            ut = np.nonzero(self.tahe[0, :] == i)[0]
            uv = np.nonzero(self.tahe[1, :] == i)[0]
            seglist = np.hstack((seglist, ut, uv))
        seglist = np.unique(seglist)
        return(seglist)

    def extrseg(self):
        """ calculate extremum of segments

        Notes
        -----

        update the following members
            `min_sx`
            `max_sx`
            `min_sy`
            `max_sy`
        Used in seginframe

        """
        # 2 x Np
        pt = self.pt
        # tahe 2 x Nseg
        #th = zip(self.tahe[0, :], self.tahe[1, :])
        ta = self.tahe[0,:]
        he = self.tahe[1,:]

        self.max_sx = np.maximum(pt[0,ta],pt[0,he])
        self.min_sx = np.minimum(pt[0,ta],pt[0,he])
        self.max_sy = np.maximum(pt[1,ta],pt[1,he])
        self.min_sy = np.minimum(pt[1,ta],pt[1,he])
        #self.max_sx = np.array([ np.maximum(pt[0, x[0]], pt[0, x[1]]) for x in th ])
        #self.min_sx = np.array([ np.minimum(pt[0, x[0]], pt[0, x[1]]) for x in th ])
        #self.max_sy = np.array([ np.maximum(pt[1, x[0]], pt[1, x[1]]) for x in th ])
        #self.min_sy = np.array([ np.minnimum(pt[1, x[0]], pt[1, x[1]]) for x in th ])

    def seginframe2(self, p1, p2):
        """ returns the seglist of a given zone defined by two points
        (vectorised version)

            Parameters
            ----------

            p1 array (2 x N)
                array of N 2D points
            p2 array (2 x N)
                array of N 2D points

            Returns
            -------

            seglist
                list of segment number inside a planar region defined by p1 an p2
                separated by -1 


            Examples
            --------

            .. plot::
                :include-source:

                >>> from pylayers.gis.layout import *
                >>> L = Layout('TA-Office.lay')
                >>> p1 = np.array([[0,0,0],[0,0,0]])
                >>> p2 = np.array([[10,10,10],[10,10,10]])
                >>> seglist = L.seginframe2(p1,p2)
                >>> edlist  = [ L.tsg[x] for x in  seglist ]
                >>> fig,ax = L.showG('s',edlist=edlist)

        """

        sh1 = np.shape(p1)
        sh2 = np.shape(p2)

        assert sh1[0] == 2
        assert sh2[0] == 2

        if (len(sh1) < 2) & (len(sh2) > 1):
            p1 = np.outer(p1, np.ones(sh2[1]))

        if (len(sh2) < 2) & (len(sh1) > 1):
            p2 = np.outer(p2, np.ones(sh1[1]))

        if (len(sh2) < 2) & (len(sh1) < 2):
            p1 = np.outer(p1, np.ones(1))
            p2 = np.outer(p2, np.ones(1))

        # clipping conditions to keep segment
        #
        # max_sx > min_x
        # min_sx < max_x
        # max_sy > min_y
        # min_sy < max_y

        # N x 1

        #max_x = [ max(x[1], x[0]) for x in  zip(p1[0, :], p2[0, :]) ]
        #min_x = [ min(x[1], x[0]) for x in  zip(p1[0, :], p2[0, :]) ]
        #max_y = [ max(x[1], x[0]) for x in  zip(p1[1, :], p2[1, :]) ]
        #min_y = [ min(x[1], x[0]) for x in  zip(p1[1, :], p2[1, :]) ]

        max_x = np.maximum(p1[0,:],p2[0,:])
        min_x = np.minimum(p1[0,:],p2[0,:])
        max_y = np.maximum(p1[1,:],p2[1,:])
        min_y = np.minimum(p1[1,:],p2[1,:])

        seglist = [ np.nonzero((self.max_sx > x[0]) &
                               (self.min_sx < x[1]) &
                               (self.max_sy > x[2]) &
                               (self.min_sy < x[3]))[0]
                   for x in zip(min_x, max_x, min_y, max_y) ]

        # np.array stacking
        # -1 acts as a deliminiter (not as a segment number)

        # seglist = reduce(lambda x, y: np.hstack((x, array([-1]), y)), seglist)
        x = np.array([]).astype(int)
        for y in seglist:
            x = np.hstack((x, np.array([-1]), y))

        return(x)

    def seginframe(self, p1, p2):
        """ return the seg list of a given zone defined by two points

            Parameters
            ----------

            p1
                array (1 x 2)
            p2
                array (1 x 2)

            Returns
            -------

            seglist
                list of segment number inside a planar region defined by p1 an p2


            Examples
            --------

            >>> from pylayers.gis.layout import *
            >>> L = Layout('TA-Office.lay')
            >>> p1 = np.array([0,0])
            >>> p2 = np.array([10,10])
            >>> L.seginframe(p1,p2)
            array([ 1,  3,  7,  8, 14, 15, 16, 17, 18, 20, 21, 23, 24, 26, 27, 29, 30,
                   32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 46, 47, 52, 53, 54,
                   55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
                   72, 73, 74, 75, 76, 77, 78, 81, 82, 85, 86])

        """
        #assert( (p1.shape==(1,2)) or (p1.shape==(2)))
        #assert( (p2.shape==(1,2)) or (p2.shape==(2)))
        max_x = max(p1[0], p2[0])
        min_x = min(p1[0], p2[0])
        max_y = max(p1[1], p2[1])
        min_y = min(p1[1], p2[1])

        Dx = max_x - min_x
        Dy = max_y - min_y

        if Dx < 0.5:
            max_x = max_x + 0.5
            min_x = min_x - 0.5

        if Dy < 0.5:
            max_y = max_y + 0.5
            min_y = min_y - 0.5

        if (Dy < Dx):
            up = np.nonzero((self.pt[0, :] < max_x) &
                            (self.pt[0, :] > min_x))[0]
        else:
            up = np.nonzero((self.pt[1, :] < max_y) &
                            (self.pt[1, :] > min_y))[0]

        seglist = self.segpt(up)

        return(seglist)

    def layerongrid(self, grid, Tx):
        """ grid Nx,Ny,2
        Tx   1x2
        .. todo:: layeron grid Not finished
        """
        Nx = grid.shape[0]
        Ny = grid.shape[1]

        for ix in range(Nx):
            for iy in range(Ny):
                p = grid[ix, iy, :]
                seglist, theta = self.layeronlink(p, Tx)

    def cycleinline(self, c1, c2):
        """ returns the intersection between a given line and all segments

        Parameters
        ----------

        c1 : int
            point
        c2 : int
            point

        Returns
        -------

        I : numpy.ndarray

        See Also
        --------

        pylayers.antprop.signature.Signatures.rays
        pylayers.gis.layout.Layout.seginframe2

        Notes
        -----

        This function is used to detect LOS conditions

        """

        I = np.array([]).reshape(3, 0)

        # polygon cycle 1
        poly1 = self.Gt.node[c1]['polyg']
        p1t = poly1.centroid.xy

        # polygon cycle 2
        poly2 = self.Gt.node[c2]['polyg']
        p2t = poly2.centroid.xy

        # centroid of cycle 1 and 2
        p1 = np.array([p1t[0][0], p1t[1][0]])
        p2 = np.array([p2t[0][0], p2t[1][0]])

        line = sh.LineString((p1, p2))

        # els = self.seginframe(p1,p2)
        # new implementation of seginframe is faster
        els = self.seginframe2(p1, p2)

        elg = self.tsg[els]

        lc = []
        ls = []
        I = np.array([]).reshape(2, 0)

        for seg in elg:
            #v1.1 ta, he = self.Gs.neighbors(seg)
            ta, he = self.Gs[seg]
            pa = np.array(self.Gs.pos[ta])
            pb = np.array(self.Gs.pos[he])

            segline = sh.LineString((pa, pb))

            if line.intersects(segline):
                lc.extend(self.Gs.nodes[seg]['ncycles'])
            # printseg,self.Gs.nodes[seg]['ncycles']
                ls.append(seg)
                psh = line.intersection(segline)
                I = np.hstack((I, np.array([[psh.x], [psh.y]])))
        v = (I - p1[:, np.newaxis])
        dv = np.sum(v * v, axis=0)
        u = np.argsort(dv)
        lss = np.array(ls)[u]

        lc = [c1]
        for s in lss:
            cy1, cy2 = self.Gs.nodes[s]['ncycles']
            if cy1 not in lc:
                lc.append(cy1)
            elif cy2 not in lc:
                lc.append(cy2)
            else:
                assert NameError('Bad transisiton in Layout.cycleinline')
        return lc

    def seginline(self, p1, p2):
        """ returns the intersection between a given line and all segments

        Parameters
        ----------
            p1 : numpy.ndarray
            p2 : numpy.ndarray

        Returns
        -------
            I : numpy.ndarray
        """
        I = np.array([]).reshape(3, 0)
        line = sh.LineString((p1, p2))
        for seg in self.Gs.nodes():
            if seg > 0:
                # v1.1 ta, he = self.Gs.neighbors(seg)
                ta, he = self.Gs[seg] 
                pa = np.array(self.Gs.pos[ta])
                pb = np.array(self.Gs.pos[he])
            else:
                pa = np.array(self.Gs.pos[seg])
                pb = pa

            segline = sh.LineString((pa, pb))
            if line.intersects(segline):
                psh = line.intersection(segline)
                liseg = np.array([[psh.x], [psh.y]])
                I = np.hstack((I, np.vstack(([[seg]], liseg))))
        return I

    def visilist(self, p):
        """ returns the list of nodes which are visible from point p

        Parameters
        ----------
        p
            np.array point

        Returns
        -------


        Notes
        -----

        AAS = [0:2pi]
        While (AAS != void set)
             1) Find segment ns either
                i)  the closest segment from p in AAS
                ii) neighbor of prec(ns)
             2) Find the edgelist visible from ns
            edgelist = vedgelist(ns)
             3) Check_occultation(p,ns,edgelist)
                Occultation 8  situations
                [p1,pM,p2] = [T,T,T]  : fully occulted
                         [     ]    partially visible
                         [F,F,F]  : fully visible
             4) Update Allowed Angular Sector  (AAS)

        """
        AAS = Intvl([0, 2 * pi])
        nsprev = np.inf
        edgelist = np.array([])

        while AAS.measure() != 0:
            if nsprev == np.inf:
                ns = self.closest(p, AAS)
            else:
                ns = self.neighbors(nsprev)
            edgelist = self.vedgelist(ns)
            [b1, bM, b2] = self.check - occultation(p, ns, edgelist)
            AAS = self.update(AAS,)

    def closest_edge(self, p, AAS):
        """ not implemented

        Parameters
        ----------

        This function return the closest segment from p which belong to
        the AAS (Allowed Angular Sector)

        [ns] = closest_edge(self,p,AAS)

        """
        pass
        # not implemented

    def visi_papb(self, pa, pb, edgelist=np.array([])):
        """
        visi_papb : determine if pa and pb are in visibility for the structure graph

        visi_papb(pa,pb,edgelist)

        pa       : 1x2
        pb       : 1x2
        edgelist : exclusion edge list

        """
        #
        # .. todo: avoid utilisation tahe
        #
        x = self.pt[0, :]
        y = self.pt[1, :]
        ta = self.tahe[0, :]
        he = self.tahe[1, :]

        x1 = x[ta]
        y1 = y[ta]
        x2 = x[he]
        y2 = y[he]

        den = (pb[1] - pa[1]) * (x2 - x1) - (pb[0] - pa[0]) * (y2 - y1)
        w = np.nonzero(abs(den) < 1e-12)[0]

        den[w] = 1e-12
        numa = (pb[0] - pa[0]) * (y1 - pa[1]) - (pb[1] - pa[1]) * \
            (x1 - pa[0])
        numb = (x2 - x1) * (y1 - pa[1]) - (y2 - y1) * (x1 - pa[0])

        ua = numa / den
        ub = numb / den

        #ua[edgelist] = 1000
        u = np.nonzero((ua >= 0) & (ua <= 1) & (ub >= 0) & (ub <= 1))[0]

    # Si le segment de droite pa-pb intercepte des paroies de la structure
        if (u != []):
            visi = 0
        else:
            visi = 1

        return(visi)

    def show_nodes(self, ndlist=[1e8], size=10, color='b', dlabels=False, font_size=15, alpha=1, node_shape='o', fig=[], ax=[]):
        """ show nodes

        Parameters
        ----------
        ndlist
        size   : int
            default 10
        color :  'b'
        dlabels : Boolean
            False
        font_size : int
            15
        alpha : float
            transparancy

        See Also
        --------

        show_segment
        showGs

        """
        if fig == []:
            fig = plt.figure()
        if ax == []:
            ax = fig.add_subplot(111)

        if type(ndlist) == np.ndarray:
            ndlist = list(ndlist)
        if len(ndlist) == 0:
            # ndlist.append(1e8)
            dlabels = False
        elif ndlist[0] == 1e8:
            ndlist = self.Gs.nodes.keys()
        # elif ndlist[0]==1e8:
        #    ndlist  = self.Gs.nodes.keys()

        # printndlist

        Z = nx.draw_networkx_nodes(self.Gs, self.Gs.pos, node_color=color,
                                   node_size=size, nodelist=ndlist, alpha=alpha,
                                   node_shape=node_shape, fig=fig, ax=ax)

        try:
            fig = Z.figure
            ax = Z.axes
        except:
            pass
        if dlabels:
            dicopos = {}
            dicolab = {}
            for n in ndlist:
                dicopos[n] = np.array(self.Gs.pos[n])
                dicolab[n] = self.labels[n]
            Z = nx.draw_networkx_labels(self.Gs, dicopos, dicolab,
                                        font_size=font_size, font_color=color, fig=fig, ax=ax)
            try:
                fig = Z.figure
                ax = Z.axes
            except:
                pass

        return fig, ax

    def show_seg1(self, edlist=[], alpha=1, width=1, size=2, color='black', font_size=15, dlabels=False):
        """ show segment

        Parameters
        ----------

        edlist
        alpha
        width
        size
        color
        font_size
        dlabels

        """
        if type(edlist) == 'ndarray':
            edlist = edlist.tolist()
        elif type(edlist) == int:
            edlist = [edlist]

        # printndlist
        nx.draw_networkx_nodes(
            self.Gs, self.Gs.pos, node_size=size, nodelist=edlist)
        if dlabels:
            dicopos = {}
            dicolab = {}
            for n in ndlist:
                # dicopos[n]=tuple(np.array(self.Gs.pos[n])+np.array((0.8,0.2)))
                dicopos[n] = np.array(self.Gs.pos[n])
                dicolab[n] = self.labels[n]
            nx.draw_networkx_labels(
                self.Gs, dicopos, dicolab, font_size=font_size)

    def show_segment(self, **kwargs):
        """ show segment

        Parameters
        ----------

        edlist : list
            segment list
        alpha : float
            transparency 0< alpha < 1
        width : float
            line width (default 1)
        color : string
            default 'black'
        dnodes : boolean
            display nodes ( Default False)
        dlabels : boolean
            display labels ( Default False)
        font_size : int
            Default 15

        See Also
        --------

        show_nodes

        """

        defaults = {'fig': [],
                    'ax': [],
                    'edlist': [],
                    'alpha': 1,
                    'width': 1,
                    'color': 'black',
                    'dnodes': False,
                    'dlabels': False,
                    'font_size': 15,
                    'node_shape': 'o'
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] == []:
            fig = plt.figure()
        else:
            fig = kwargs['fig']
        if kwargs['ax'] == []:
            ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']
        clrlist = []
        cold = pyu.coldict()

        # html color or string
        if kwargs['color'][0] != '#':
            clrlist.append(cold[kwargs['color']])
        else:
            if color == '#FFFFF0':
                color = '#00000F'
            clrlist.append(color)
        ecmap = clr.ListedColormap(clrlist)

        U = self.Gs.edges(kwargs['edlist'])
        # ue = (np.ones(2 * len(kwargs['edlist']))).astype('int').tolist()
        ue = np.ones(len(U),dtype='int').tolist()
        if len(U) > 0:
            Z = nx.draw_networkx_edges(self.Gs, self.Gs.pos, edgelist=U,
                                       edge_color=ue, edge_cmap=ecmap,
                                       alpha=kwargs['alpha'], width=kwargs['width'], fig=fig, ax=ax)
            try:
                fig = Z.figure
                ax = Z.axes
            except:
                pass

        if kwargs['dlabels']:
                # printedlist
                # nodelist = self.ed2nd(edlist)
            fig, ax = self.show_nodes(ndlist=kwargs['edlist'], dlabels=kwargs['dlabels'],
                                      color='b', font_size=kwargs['font_size'],
                                      node_shape=kwargs['node_shape'], fig=fig, ax=ax)
        if kwargs['dnodes']:
            fig, ax = self.show_nodes(ndlist=kwargs['edlist'], color='b', fig=fig, ax=ax)

        return fig, ax

    def show_layer(self, name, edlist=[], alpha=1, width=0,
                   color='black', dnodes=False, dthin=False,
                   dlabels=False, font_size=15, fGHz=[], fig=[], ax=[]):
        """ show layer

        Parameters
        ----------

        name :
        edlist : []
        alpha : float
            transparency
        width : int
            if width = 0 width depends on slab property
        color : string
            default black'
        dnodes :
            display nodes (False )
        dthin :
            display thin ( False )
        dlabels :
            display labels ( False )
        font_size

        See Also
        --------

        show_segment


        """

        if fig == []:
            fig = plt.figure()
        if ax == []:
            ax = fig.add_subplot(111)

        if edlist == []:
            edlist = self.name[name]
        else:
            # intersect layer edge list with local zone edge list (in function
            # argument)
            a1 = np.array(self.name[name])
            a2 = np.array(edlist)
            edlist = list(np.intersect1d(a1, a2))

        if self.display['thin']:
            fig, ax = self.show_segment(edlist=edlist,
                                        alpha=1,
                                        width=1,
                                        color=color,
                                        dlabels=dlabels,
                                        font_size=font_size, fig=fig, ax=ax)
        else:
            slab = self.sl[name]
            if width == 0:
                linewidth = slab['linewidth'] / 3.
            else:
                linewidth = width
            if fGHz == []:
                color = slab['color']
            else:
                if (name != 'METAL') & (name != 'METALIC'):
                    color = slab.tocolor
                else:
                    color = 'black'

            fig, ax = self.show_segment(edlist=edlist, alpha=1,
                                        width=linewidth, color=color, dnodes=dnodes,
                                        dlabels=dlabels, font_size=font_size, fig=fig, ax=ax)
        return fig, ax

    def _showGi(self, **kwargs):
        """  show graph of interactions Gi

        Parameters
        ----------
        seed : float
        alpha : float
            transparency
        sig : list of signatures (isequence of Gi nodes format)
        cycles : list
            [cystart,cyend]
        ninter : int
            interaction index
        inter : tuple
            interaction tuple

        See Also
        --------

        Signatures.siginter

        """
        defaults = {'seed':1,
                    'alpha':0.4,
                    'sig':[],
                    'cycles':[],
                    'ninter':0,
                    'node_size':30,
                    'fontsize':18,
                    'labels':False,
                    'inter':[]}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        edges = self.Gi.edges()
        cy = kwargs['cycles']
        if cy!=[]:
            pstart = self.Gt.pos[cy[0]]
            pstop = self.Gt.pos[cy[1]]

        if kwargs['sig']!=[]:
            lsig = kwargs['sig']
            edgelist = []
            startlist = []
            stoplist = []
            phe_start = np.array([])
            phe_stop = np.array([])
            phe_start.shape = (2,0)
            phe_stop.shape = (2,0)
            for sig in lsig:
                edgelist = edgelist + list(zip(sig[0:-1],sig[1:]))
                if cy!=[]:
                    p1 =  np.array(self.Gi.pos[sig[0]])[:,None]
                    p2 =  np.array(self.Gi.pos[sig[-1]])[:,None]
                    phe_start=np.hstack((phe_start,p1))
                    phe_stop=np.hstack((phe_stop,p2))
        elif kwargs['inter']!=[]:
            edinter = kwargs['inter']
            outlist = self.Gi[edinter[0]][edinter[1]]['output']
            outprob = outlist.values()
            edgelist = [(edinter[1],x) for x in outlist]
            dprob = dict(zip(edgelist,[str(x) for x in outprob]))
        elif kwargs['ninter']!=[]:
            edinter = [ e for e in edges][kwargs['ninter']]
            outlist = self.Gi[edinter[0]][edinter[1]]['output']
            outprob = outlist.values()
            edgelist = [(edinter[1],x) for x in outlist]
            dprob = dict(zip(edgelist,[str(x) for x in outprob]))
        else:
            pass


        ns = kwargs['node_size']
        np.random.seed(kwargs['seed'])
        fig = plt.figure(figsize=(20,10))
        ax1 = plt.subplot(121)
        pos = nx.spring_layout(self.Gi)
        nx.draw_networkx_nodes(self.Gi,pos,nodelist=[x for x in self.Gi.nodes() if len(x)==1],
                node_color='r',node_size=ns,ax=ax1,alpha=kwargs['alpha'])
        nx.draw_networkx_nodes(self.Gi,pos,nodelist=[x for x in self.Gi.nodes() if len(x)==2],
                node_color='b',node_size=ns,ax=ax1,alpha=kwargs['alpha'])
        nx.draw_networkx_nodes(self.Gi,pos,nodelist=[x for x in self.Gi.nodes() if len(x)==3],
                node_color='g',node_size=ns,ax=ax1,alpha=kwargs['alpha'])
        nx.draw_networkx_edges(self.Gi,pos,edgelist=self.Gi.edges(),width=.1,edge_color='k',arrow=False,ax=ax1)
        if (kwargs['sig']==[]):
            nx.draw_networkx_edges(self.Gi,pos,edgelist=[edinter],width=2,edge_color='g',arrow=False,ax=ax1)
        nx.draw_networkx_edges(self.Gi,pos,edgelist=edgelist,width=2,edge_color='r',arrow=False,ax=ax1)
        ax2 = plt.subplot(122)
        fig,ax2 = self.showG('s',aw=1,ax=ax2)
        nx.draw_networkx_nodes(self.Gi,self.Gi.pos,nodelist=[x for x in self.Gi.nodes() if len(x)==1],
                node_color='r',node_size=ns,ax=ax2,alpha=kwargs['alpha'])
        nx.draw_networkx_nodes(self.Gi,self.Gi.pos,nodelist=[x for x in self.Gi.nodes() if len(x)==2],
                node_color='b',node_size=ns,ax=ax2,alpha=kwargs['alpha'])
        nx.draw_networkx_nodes(self.Gi,self.Gi.pos,nodelist=[x for x in self.Gi.nodes() if len(x)==3],
                node_color='g',node_size=ns,ax=ax2,alpha=kwargs['alpha'])
        nx.draw_networkx_edges(self.Gi,self.Gi.pos,edgelist=self.Gi.edges(),width=.1,edge_color='k',arrow=False,ax=ax2)
        if kwargs['labels']:
            nx.draw_networkx_labels(self.Gi,self.Gi.pos,labels=[str(x) for x in self.Gi.nodes()],ax=ax2,fontsize=kwargs['fontsize'])
        if (kwargs['sig']==[]):
            nx.draw_networkx_edges(self.Gi,self.Gi.pos,edgelist=[edinter],width=2,edge_color='g',arrow=False,ax=ax2)
        nx.draw_networkx_edges(self.Gi,self.Gi.pos,edgelist=edgelist,width=2,edge_color='r',arrow=False,ax=ax2)
        if (kwargs['sig']==[]):
            nx.draw_networkx_edge_labels(self.Gi,self.Gi.pos,edge_labels=dprob,ax=ax2,fontsize=kwargs['fontsize'])
        if cy!=[]:
            ptstart = pstart[:,None]*np.ones(phe_start.shape[1])[None,:]
            ptstop = pstop[:,None]*np.ones(phe_start.shape[1])[None,:]
            plu.displot(ptstart,phe_start,ax=ax2,arrow=True)
            plu.displot(phe_stop,ptstop,ax=ax2,arrow=True)
        # interactions corresponding to edge en
#        int0, int1 = self.Gi.edges()[kwargs['en']]
#
#        print("int0 : ", int0)
#        print("int1 : ", int1)
#
#        # if interaction is tuple (R or T)
#        if ((len(int0) > 1) & (len(int1) > 1)):
#            nstr0 = int0[0]
#            nstr1 = int1[0]
#            e01 = self.Gi.edge[int0][int1]
#            lseg = []
#            if e01.has_key('output'):
#                output = e01['output']
#                print(" output ", output)
#                ltup = filter(lambda x: type(x) == tuple, output.keys())
#                lref = filter(lambda x: len(x) == 2, ltup)
#                ltran = filter(lambda x: len(x) == 3, ltup)
#                lseg = np.unique(np.array(map(lambda x: x[0], output.keys())))
#                probR = np.array(map(lambda x: output[x], lref))
#                segR = np.array(map(lambda x: x[0], lref))
#                probT = np.array(map(lambda x: output[x], ltran))
#                segT = np.array(map(lambda x: x[0], lref))
#                dprobR = dict(zip(segR, probR))
#                dprobT = dict(zip(segT, probT))
#            # print" Sum pR : ",sum(dprobR.values())
#            # print" Sum pT : ",sum(dprobT.values())
#            # print"lseg", lseg
#            # termination points from seg0 and seg1
#            pseg0 = self.s2pc[nstr0].toarray().reshape(2, 2).T
#            pseg1 = self.s2pc[nstr1].toarray().reshape(2, 2).T
#            #
#            # create the cone seg0 seg1
#            #
#            cn = cone.Cone()
#            cn.from2segs(pseg0, pseg1)
#            # show cone
#            # show Gt
#            self.display['thin'] = True
#            self.display['subseg'] = False
#            fig, ax = self.showG('s',aw=1,labels=True)
#            fig, ax = cn.show(fig=fig, ax=ax)
#            for nse in lseg:
#                ta, he = self.Gs.neighbors(nse)
#                pta = np.array(self.Gs.pos[ta])
#                phe = np.array(self.Gs.pos[he])
#
#                try:
#                    pR = dprobR[nse]
#                except:
#                    pR = 0
#
#                try:
#                    pT = dprobT[nse]
#                except:
#                    pT = 0
#
#                alpha = (pR + pT) / 2.
#                segment = ax.plot([pta[0], phe[0]],
#                                  [pta[1], phe[1]],
#                                  'g', linewidth=7, visible=True, alpha=alpha)
#
        return(fig, ax1)

    def _showGt(self, ax=[], roomlist=[], mode='indoor'):
        """ show topological graph Gt

        Parameters
        -----------
        ax : matlplotlib axes
        roomlist : list
            list of room numbers
        mode : string
            'indoor','open','area','start'

        """
        if not isinstance(ax, plt.Axes):
            fig = plt.gcf()
            ax = fig.gca()

        G = self.Gt

        for k, nc in enumerate(G.node.keys()):
            if nc!=0:
                poly = G.node[nc]['polyg']

                a = poly.signedarea()

                if mode == 'area':
                    if a < 0:
                        poly.plot(color='red', alpha=0.5, fig=fig, ax=ax)
                    else:
                        poly.plot(color='green', alpha=0.5, fig=fig, ax=ax)

                if mode == 'start':
                    if poly.vnodes[0] < 0:
                        poly.plot(color='blue', alpha=0.5, fig=fig, ax=ax)
                    else:
                        poly.plot(color='yellow', alpha=0.5, fig=fig, ax=ax)

                if mode == 'indoor':
                    if G.node[nc]['indoor']:
                        poly.plot(color='green', alpha=0.5, fig=fig, ax=ax)
                    else:
                        poly.plot(color='blue', alpha=0.5, fig=fig, ax=ax)

                if mode == 'open':
                    if G.node[nc]['isopen']:
                        poly.plot(color='green', alpha=0.5, fig=fig, ax=ax)
                    # else:
                    #     poly.plot(color='blue', alpha=0.5,fig=fig,ax=ax)

        ax.axis('scaled')

    def showGs(self, **kwargs):
        """ show structure graph Gs


        Parameters
        ----------

        ndlist  : np.array
            set of nodes to be displayed
        edlist  : np.array
            set of edges to be displayed
        roomlist : list
            default : []
        axis :
        width : int
            2
        fGHz : float
        show    : boolean
            default True
        furniture : boolean
            default False

        display parameters are defined in  display dictionnary

        Returns
        -------

        ax

        See Also
        --------

        pylayers.gis.layout.showG

        """

        defaults = {'ndlist': [],
                    'edlist': [],
                    'roomlist': [],
                    'axis': [],
                    'width': 2,
                    'fGHz': [],
                    'show': False,
                    'furniture': False,
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        args = {}
        for k in kwargs:
            if k not in defaults:
                args[k] = kwargs[k]

        if 'fig' not in kwargs:
            fig = plt.figure()
        else:
            fig = kwargs['fig']

        if 'ax' not in kwargs:
            ax = fig.add_subplot(111)
        else:
            ax = kwargs['ax']

        if self.display['clear']:
            ax.cla()

        # display overlay image
        if self.display['overlay']:
            # imok : Image is OK
            imok = False
            if len(self.display['overlay_file'].split('http:')) > 1:
                #img_file = urllib.urlopen(self.display['overlay_file'])
                img_file = urlopen(self.display['overlay_file'])
                #im = StringIO(img_file.read())
                image = Image.open(im)
                imok = True
            else:
                if self.display['overlay_file'] != '':
                    image = Image.open(os.path.join(
                        pro.basename, pro.pstruc['DIRIMAGE'], self.display['overlay_file']))
                    imok = True
            if imok:
                if 'v' in self.display['overlay_flip']:
                    image = image.transpose(Image.FLIP_LEFT_RIGHT)
                if 'h' in self.display['overlay_flip']:
                    image = image.transpose(Image.FLIP_TOP_BOTTOM)
                ax.imshow(image, extent=self.display[
                          'overlay_axis'], alpha=self.display['alpha'], origin='lower')

        if kwargs['ndlist'] == []:
            tn = np.array(list(self.Gs.nodes.keys()))
            u = np.nonzero(tn < 0)[0]
            ndlist = tn[u]

        if kwargs['edlist'] == []:
            tn = self.Gs.nodes.keys()
            #u  = np.nonzero(tn > 0)[0]
            #edlist = tn[u]
            edlist = filter(lambda x: (x > 0), tn)
            #& (not self.Gs.nodes[x].has_key('ss_name')),tn)
        else:
            edlist = kwargs['edlist']

        if self.display['nodes']:
            dlabels = self.display['ndlabel']
            fig, ax = self.show_nodes(
                ndlist, size=30, color='k', dlabels=dlabels, node_shape='s', fig=fig, ax=ax)

        if self.display['isonb']:
            if hasattr(self,'lsss'):
                seg = [x for x in self.Gs.nodes() if x >0]
                # psseg = np.array([[self.Gs.pos[x][0],self.Gs.pos[x][1]] for x in seg])
                # nbsseg = np.array([len(self.Gs.node[x]['iso']) for x in seg],dtype='int')
                try:
                    psseg = np.array([[self.Gs.pos[x][0],self.Gs.pos[x][1]] for x in seg 
                                   if len(self.Gs.nodes[x]['iso']) >1])  
                except:
                    import ipdb
                    ipdb.set_trace()

        #         [ax.text(psseg[x,0]+0.2,psseg[x,1]+0.2,str(nbsseg[x]),
        # fontdict={'size':8},ha='center') for x in range(len(seg))]
                [ax.text(psseg[x,0]+0.2,psseg[x,1]+0.2,'+',
                fontdict={'size':8},ha='center') for x in range(len(psseg))]

        if self.display['transition']:
            try:
                segwtrans = [y for y in [x for x in self.Gs.nodes() if x > 0]if
                        self.Gs.nodes[y]['transition']]
                posseg = np.array([self.Gs.pos[x] for x in segwtrans])
                normseg = np.array([self.Gs.nodes[x]['norm']
                                    for x in segwtrans])[:, :2]
                b1 = (posseg - normseg / 2)
                b2 = (posseg + normseg / 2)
                [ax.annotate('', xy=b1[x],
                             xycoords='data',
                             xytext=b2[x],
                             textcoords='data',
                             arrowprops={'arrowstyle': '<->'})
                 for x in range(len(segwtrans))]
            except:
                pass
        slablist = self.name.keys()
        if self.display['edges']:
            dlabels = self.display['edlabel']
            font_size = self.display['fontsize']
            dnodes = self.display['ednodes']
            dthin = self.display['thin']
            alpha = self.display['alpha']
            for nameslab in self.name:
                color = self.sl[nameslab]['color']
                edlist = self.name[nameslab]
                fig, ax = self.show_layer(nameslab, edlist=edlist, alpha=alpha,
                                          dthin=dthin, dnodes=dnodes, dlabels=dlabels,
                                          color=color,
                                          font_size=font_size,
                                          width=kwargs['width'],
                                          fGHz=kwargs['fGHz'],
                                          fig=fig, ax=ax)

        if self.display['subseg']:
            dico = self.subseg()
            for k in dico.keys():
                if kwargs['fGHz'] == []:
                    color = self.sl[k]['color']
                else:
                    if (k != 'METAL') & (k != 'METALIC'):
                        color = self.sl[k].tocolor(fGHz)
                        #color = 'red'
                    else:
                        color = 'black'
                        # printk,color
                edlist2 = []
                for ts in dico[k]:
                    edlist2.append(ts[0])
                    # edlist2.append(ts)
                edlist3 = list(set(edlist2).intersection(set(edlist)))
                # printk , color , edlist
                fig, ax = self.show_segment(
                    edlist=edlist3, color=color, alpha=1.0, width=2, fig=fig, ax=ax)

        if self.display['scaled']:
            ax.axis('scaled')
        ax.set_title(self.display['title'])
        #fig = plt.gcf()
        #ax  = fig.axes[0]
        #
        # TODO Not working in python 3
        #if self.display['ticksoff']:
        #    ax.xaxis.set_ticks([])
        #    for loc, spine in ax.spines.iteritems():
        #        spine.set_color('none')

        if kwargs['furniture']:
            if 'lfur' in self.__dict__:
                for fur1 in self.lfur:
                    if fur1.Matname == 'METAL':
                        fig, ax = fur1.show(fig, ax)
            else:
                print("Warning : no furniture file loaded")

        for nr in kwargs['roomlist']:
            ncy = self.Gr.node[nr]['cycle']
            fig, ax = self.Gt.node[ncy]['polyg'].plot(fig=fig, ax=ax)
        if kwargs['axis'] == []:
            ax.axis('scaled')
        else:
            ax.axis(kwargs['axis'])

        if kwargs['show']:
            plt.show()

        return fig, ax

    #@profile
    def build(self, graph='tvirw', verbose=False, difftol=0.15, multi=False):
        """ build graphs

        Parameters
        ----------

        graph : string composed of
            't' : Gt
            'v' : Gv
            'i' : Gi
            'r' : Gr
            'w" : Gw
        verbose : boolean
        difftol : diffraction tolerance
        multi : boolean
            enable multi processing

        Notes
        -----

        This function builds all the graph associated with the Layout. 

        Warning : by default the layout is saved (dumpw) after each build

        """
        # check that the graph and numpy table are aligned
        # if not call g2npy
        if self.pt.shape[1] != self.Np:
            self.g2npy()

        if not self.hasboundary:
            self.boundary()

        # to save graoh Gs
        # list of built graphs
        if 's' not in self.lbltg:
            self.lbltg.extend('s')

        Buildpbar = pbar(verbose,total=5,desc='Build Layout',position=0)

        if verbose:
            Buildpbar.update(1)
        if 't' in graph:
            logger.info('buildGt')
            self.buildGt(difftol=difftol, verbose=verbose, tqdmpos=1)
            self.lbltg.extend('t')
        if verbose:
            Buildpbar.update(1)
        if 'v' in graph:
            logger.info('buildGv')
            self.buildGv(verbose=verbose, tqdmpos=1)
            self.lbltg.extend('v')
        if verbose:
            Buildpbar.update(1)
        if 'i' in graph:
            logger.info('buildGi')
            self.buildGi(verbose=verbose, tqdmpos=1)
            if not multi:
                logger.info('outputGi')
                self.outputGi(verbose=verbose,tqdmpos=1)
            else:
                self.outputGi_mp()
            self.lbltg.extend('i')
        if verbose:
            Buildpbar.update(1)
        # if 'r' in graph:
        #     if verbose:
        #         print"Gr"
        #     self.buildGr()
        #     self.lbltg.extend('r')

        # if 'w' in graph and len(self.Gr.nodes())>1:
        #     self.buildGw()
        #     self.lbltg.extend('w')

        # add hash to node 0 of Gs

        filelay = pyu.getlong(self._filename, pro.pstruc['DIRLAY'])
        fd = open(filelay,'rb')
        _hash = hashlib.md5(fd.read()).hexdigest()
        fd.close()
        self.Gt.add_node(0, hash=_hash)

        # There is a dumpw after each build
        self.dumpw()
        self.isbuilt = True
        if verbose:
            Buildpbar.update(1)

    def dumpw(self):
        """ pickle dump of specified Graphs

        Notes
        -----

        graphs which are in lbltg are saved in pickle format

        't' : Gt
        's' : Gs
        'v' : Gv
        'i' : Gi
        'r' : Gr

        """
        # create layout directory
        if os.path.splitext(self._filename)[1]=='.ini':
            dirname = self._filename.replace('.ini','')
        if os.path.splitext(self._filename)[1]=='.lay':
            dirname = self._filename.replace('.lay','')
        path = os.path.join(pro.basename, 'struc', 'gpickle', dirname)

        if not os.path.isdir(path):
            os.mkdir(path)
        for g in self.lbltg:
            try:
                # if g in ['v','i']:
                #     gname1 ='G'+g
                #     write_gpickle(getattr(self,gname1),os.path.join(basename,'struc','gpickle','G'+g+'_'+self._filename+'.gpickle'))
                # else:
                gname = 'G' + g
                write_gpickle(getattr(self, gname), os.path.join(
                    path, 'G' + g + '.gpickle'))
            except:
                raise NameError(
                    'G' + g + ' graph cannot be saved, probably because it has not been built')

        if 's' in self.lbltg:
            if hasattr(self,'sl'):
                write_gpickle(getattr(self, 'sl'),
                          os.path.join(path, 'sl.gpickle'))
            if hasattr(self,'dpoly'):
                with open(os.path.join(path, 'dpoly.pickle'),'wb') as fd:
                    pickle.dump(getattr(self,'dpoly'),fd)

        # save dictionnary which maps string interaction to
        # [interaction node, interaction type]

        if 't' in self.lbltg:
            if hasattr(self,'ddiff'):
                write_gpickle(getattr(self, 'ddiff'),
                          os.path.join(path, 'ddiff.gpickle'))
            if hasattr(self,'lnss'):
                write_gpickle(getattr(self, 'lnss'),
                          os.path.join(path, 'lnss.gpickle'))
        if hasattr(self,'dca'):
            write_gpickle(getattr(self, 'dca'), os.path.join(path, 'dca.gpickle'))
        # write_gpickle(getattr(self,'sla'),os.path.join(path,'sla.gpickle'))
        if hasattr(self, 'm'):
            write_gpickle(getattr(self, 'm'), os.path.join(path, 'm.gpickle'))

    def dumpr(self, graphs='stvirw'):
        """ read of given graphs

        Notes
        -----

        graph : string
            's' : Gv
            't' : Gt
            'r' : Gr
            'v' : Gv
            'i' : Gi


        .gpickle files are store under the struc directory of the project
        specified by the $BASENAME environment variable

        """
        if os.path.splitext(self._filename)[1]=='.ini':
            dirname = self._filename.replace('.ini','')
        if os.path.splitext(self._filename)[1]=='.lay':
            dirname = self._filename.replace('.lay','')
        path = os.path.join(pro.basename, 'struc', 'gpickle', dirname)
        for g in graphs:
            try:
                # if g in ['v','i']:
                #     gname1 ='G'+g
                #     setattr(self, gname1, read_gpickle(os.path.join(basename,'struc','gpickle','G'+g+'_'+self._filename+'.gpickle')))
                # else:
                gname = 'G' + g
                filename = os.path.join(path, 'G' + g + '.gpickle')
                G = read_gpickle(filename)
                setattr(self, gname, G)
                self.lbltg.extend(g)
            except:
                print("Warning Unable to read graph G"+g)
                pass

        # retrieve md5 sum of the original ini file
        if 's' in graphs:
            #self._hash = self.Gs.node.pop(0)['hash']
            # self._hash = self.Gs.node[0]['hash']
            # update self.name
            lseg = [x for x in self.Gs.nodes if x > 0]
            for name in self.name:
                self.name[name] = [ x for x in lseg if self.Gs.nodes[x]['name'] == name]

            # TODO not necessary useful to call g2npy
            self.g2npy()
            # TODO use a pickle file instead of gpickle
            filesl = os.path.join(path, 'sl.gpickle')
            if os.path.isfile(filesl):
                sl = read_gpickle(filesl)
                setattr(self, 'sl', sl)

            filediff = os.path.join(path, 'ddiff.gpickle')
            if os.path.isfile(filediff):
                ddiff = read_gpickle(filediff)
                setattr(self, 'ddiff', ddiff)
            else:
                self.ddiff={}

            filelnss = os.path.join(path, 'lnss.gpickle')
            if os.path.isfile(filelnss):
                lnss = read_gpickle(filelnss)
                setattr(self, 'lnss', lnss)
            else :
                self.lnss=[]

            filedpoly = os.path.join(path, 'dpoly.pickle')
            if os.path.isfile(filedpoly):
                fd = open(filedpoly,'rb')
                dpoly = pickle.load(fd)
                setattr(self, 'dpoly', dpoly)
                #self.dpoly = {k:eval(self.dpoly[k]) for k in self.dpoly}
                self.dpoly = {k:self.dpoly[k] for k in self.dpoly}
            else :
                self.dpoly = {}

        filedca = os.path.join(path, 'dca.gpickle')
        if os.path.isfile(filedca):
            dca = read_gpickle(filedca)
            setattr(self, 'dca',dca)

        #
        # TODO Replace self.m by pyproj
        #
        filem = os.path.join(path, 'm.gpickle')
        if os.path.isfile(filem):
            setattr(self, 'm', read_gpickle(filem))
            self.extent = (self.m.lonmin,self.m.lonmax,self.m.latmin,self.m.latmax)
            self.pll = self.m(self.extent[0],self.extent[2])
            self.pur = self.m(self.extent[1],self.extent[3])
            self.extent_c = (self.pll[0],self.pur[0],self.pll[1],self.pur[1])


    def polysh2geu(self, poly):
        """ transform sh.Polygon into geu.Polygon
        """
        try:

            Gsnodes = np.array(self.Gs.nodes())
            # get node coordinates
            nodept = [self.Gs.pos[i] for i in Gsnodes]
            # transform into shapely points
            shpt = [sh.Point(pt) for pt in nodept]
            # IV 1 get nodes and vnodes
            # Create a ring to avoid taking points inside the polygon.
            # This helps to avoid polygon inside polygons
            # take exterior of polygon. embose it with buffer and find difference with original polygon*.
            # polye = poly.intersection((poly.exterior).buffer(1e-3))

            uvn = np.where([poly.exterior.buffer(1e-3).contains(p)
                            for p in shpt])[0]
            vnodes = Gsnodes[uvn]
            # IV 1.b transform vnodes to an ordered cycle with Cycle class
            # NOTE ! Using class cycle is MANDATORY
            # because, some extra vnodes can be pickup during the contain
            # process before
            S = nx.subgraph(self.Gs, vnodes)
            cycle = nx.cycle_basis(S)

            if len(cycle) > 1:
                lc = np.array([len(c) for c in cycle])
                dif = abs(lc - len(vnodes))
                ud = np.where(dif == min(dif))[0]
                cycle = cycle[ud]
            else:
                cycle = cycle[0]
            if cycle[0] > 0:
                cycle = np.roll(cycle, -1)
            pos = [self.Gs.pos[c] for c in cycle if c < 0]
            # IV 1.c create a new polygon with correct vnodes and correct
            # points
            P = geu.Polygon(p=pos, vnodes=cycle)

        except:
            import ipdb
            ipdb.set_trace()
        return P

    def getangles(self, poly, unit='rad', inside=True):
        """ find angles of a polygon

        Parameters
        ----------

        poly : geu.Polygon or sh.Polygon
        unit : str
            'deg' : degree values
            'rad' : radian values
        inside : boolean
            True :  compute the inside angles of the cycle.
                    (a.k.a. the interior of the polygon) 
            False : compute the outside angles of the cycle.
                    (a.k.a. the exterior of the polygon)

        Returns
        -------

        (u,a)
        u : int (Np)
            point number
        a : float (Np)
            associated angle to the point


        Notes
        -----

        http://www.mathopenref.com/polygonexteriorangles.html

        TODO : This function should be moved in geomutil.py  (NOT USED) 

        """

        if isinstance(poly, sh.Polygon):
            poly = polysh2geu(poly)

        cycle = poly.vnodes

        upt = cycle[cycle < 0]

        # rupt=np.roll(upt,1)         # for debug
        # rupt2=np.roll(upt,-1)         # for debug
        #
        # See OSM bug fix
        #
        pt = self.pt[:, self.iupnt[-upt]]
        if geu.SignedArea(pt) < 0:
            upt = upt[::-1]
            pt = pt[:, ::-1]

        ptroll = np.roll(pt, 1, axis=1)

        v = pt - ptroll
        v = np.hstack((v, v[:, 0][:, None]))
        vn = v / np.sqrt(np.sum((v) * (v), axis=0))
        v0 = vn[:, :-1]
        v1 = vn[:, 1:]
        cross = np.cross(v0.T, v1.T)
        dot = np.sum(v0 * v1, axis=0)
        ang = np.arctan2(cross, dot)
        uneg = ang < 0
        ang[uneg] = -ang[uneg] + np.pi
        ang[~uneg] = np.pi - ang[~uneg]

        if not inside:
            ang = 2 * np.pi - ang

        if unit == 'deg':
            return upt, ang * 180 / np.pi
        elif unit == 'rad':
            return upt, ang
            # atan2(cross(a,b)), dot(a,b))

    def pltlines(self, lines, fig=[], ax=[], color='r', alpha=1):
        """  plot a line with a specified color and transparency

        Parameters
        -----------

        lines : shapely lines
        fig : matplotlib figure
        ax : figure axis
        color  : string
        alpha  : float
            transparency

        See Also
        --------

        pylayers.gis.layout.Layout.plot

        """
        if fig == []:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        c = np.array([l.xy for l in lines])
        [ax.plot(x[0, :], x[1, :], color=color, alpha=alpha) for x in c]
        plt.axis(self.ax)
        plt.draw()

    def pltpoly(self, poly, fig=[], ax=[], color='r', alpha=0.2):
        """  plot a polygon with a specified color and transparency

        TODO : To be deplaced in an ither class
        """
        if fig == []:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()
        try:
            mpl = [ PolygonPatch(x, alpha=alpha, color=color) for x in poly]
        except:
            mpl = [ PolygonPatch(x, alpha=alpha, color=color) for x in [poly]]
        [ax.add_patch(x) for x in mpl]
        plt.axis(self.ax)
        plt.draw()

    def pltvnodes(self, vn, fig=[], ax=[]):
        """ plot vnodes

        Parameters
        ----------

        vn : list of nodes
        fig :
        ax :

        """
        if fig == []:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()

        if len(vn) > 0:
            X = np.array([self.Gs.pos[x] for x in vn])
            ax.plot(X[:, 0], X[:, 1], 'or')
            [ax.text(x[0], x[1], vn[xx]) for xx, x in enumerate(X)]

        return fig, ax

    def updateshseg(self):
        """ update shapely segment

        build a shapely object for all segments

        This function is called at the beginning of buildGt.

        See Also
        --------

        buildGt

        """

        seg_connect = {x: self.Gs.nodes[x]['connect']
                       for x in self.Gs.nodes() if x > 0}
        dpts = {x[0]: (self.Gs.pos[x[1][0]], self.Gs.pos[x[1][1]])
                for x in seg_connect.items()}

        self._shseg = {p[0]: sh.LineString(p[1]) for p in dpts.items()}

    def _triangle_old(self, poly_surround, poly_holes=[], mesh_holes=False):
        """
        perfome a delaunay partitioning on shapely polygons

        Parameters
        ----------

            poly_surround : sh.Polygon
                A single polygon to be partitionned
            poly_holes : list of sh.Polygon
                A list of polygon contained inside poly_surround. they are considered as holes
            mesh_holes : bool
                If True make the delaunay partition of poly_holes
                else : only partitioning poly_surround  and traits poly_holes as holes


        Returns
        -------
        T : dict
            dictionnary from triangle.triangulate library
            T.keys()
            ['segment_markers',
             'segments',
             'holes',
             'vertices',
             'vertex_markers',
             'triangles'
             ]


        Notes
        -----

        uses triangle library

        """

        if not isinstance(poly_surround, list):
            poly_surround = [poly_surround]

        lP = poly_surround + poly_holes

        vertices = np.ndarray(shape=(2, 0))
        segments = np.ndarray(shape=(2, 0), dtype='int')
        holes = np.ndarray(shape=(2, 0))
        segcpt = 0
        for p in lP:
            pts = np.array(p.exterior.xy)[:, :-1]
            vertices = np.hstack((vertices, pts))
            nbv = pts.shape[1]

            segments = np.hstack((segments, np.array(
                [np.arange(nbv), np.mod(range(1, nbv + 1), nbv)], dtype='int') + segcpt))
            segcpt = segcpt + nbv
            if not mesh_holes:
                holes = np.hstack((holes, np.array(p.centroid.xy)))

        if not mesh_holes:
            C = {'vertices': vertices.T, 'segments': segments.T, 'holes': holes.T}
        else:
            C = {'vertices': vertices.T, 'segments': segments.T}
        import ipdb
        ipdb.set_trace()
        T = triangle.triangulate(C, 'pa')

        # import triangle.plot as plot
        # ax=plt.gca()
        # plot.plot(ax,**T)

        return T

    def _merge_polygons(self, lP):
        """ merge triangle (polygon object) to cvx polygon

        Parameters
        ----------
            lP : list
                list of polygon to be merged

        Return
        ------

            lMP : list
                list of merged polygons

        """

        lMP = []
        # MERGE POLYGONS
        # move from delaunay triangles to convex polygons
        while lP != []:
            p = lP.pop(0)
            # restrict research to polygon that are touching themself
            restp = [(ix, x) for ix, x in enumerate(lP)
                     if isinstance(p.intersection(x), sh.LineString)]
            # self.pltpoly(p,ax=plt.gca())

            conv = False
            pold = p
            # for ip2,p2 in restp:
            for ip2, p2 in restp:
                # inter = p.intersection(p2)
                # if 2 triangles have a common segment
                p = p + p2
                if p.isconvex():
                    lP.pop(ip2)
                    lP.insert(0, p)
                    conv = True
                    break
                else:
                    # if pold not in cpolys:
                    #     cpolys.append(pold)
                    p = pold
            # if (ip2 >= len(polys)):# and (conv):
            # if conv :
            #     if p not in cpolys:
            #         cpolys.append(p)
            if restp == [] and conv == True:
                lMP.append(p)
            if not conv:  # else:
                if pold not in lMP:
                    lMP.append(pold)
            if len(lP) == 0:
                if p not in lMP:
                    lMP.append(p)

        return lMP

    def _triangle(self, holes=[], vnodes=[] ,bplot = False):
        """ Delaunay triangulation on shapely polygons

        Parameters
        ----------

        holes : ndarray
                if holes ==[] : it means the merge is applied on the interior of the layout (indoor)
                if holes == np.ndarray (centroid of polygon). indoor is discarded and delaunay
                        is applied on outdoor


        Returns
        -------

        T : dict
                dictionnary from triangle.triangulate library with the following keys
                ['segment_markers', 'segments', 'holes', 'vertices', 'vertex_markers', 'triangles']
        map_vertices : points index

        Notes
        -----

        This methods uses the `triangle` library

        """

        # this means Delaunay is applied on exterior
        # and internal polygons will be discarded
        segbounds = []
        ptbounds = []

        if holes == []:
            # remove air segments around layout
            pass
            # [segbounds.extend(nx.neighbors(L.Gs,x)) for x in L.lboundary]
            # ptbounds = L.lboundary

        if vnodes == []:
            vnodes = self.Gs.nodes()

        # find termination points of segments of layout
            if nx.__version__!='1.10':
                seg = np.array([self.Gs[x]  for x in vnodes
                        if x > 0
                        and x not in segbounds])
            else:
                seg = np.array([nx.neighbors(self.Gs, x) for x in vnodes
                        if x > 0
                        and x not in segbounds])

        # get vertices/points of layout
        ivertices = np.array([(x, self.Gs.pos[x][0], self.Gs.pos[x][1]) for x in vnodes
                              if x < 0
                              and x not in ptbounds])

        # map_vertices : points negative index  (Np,)

        map_vertices = ivertices[:, 0].astype('int')
        # vertices : coordinates (Np x 2)
        vertices = ivertices[:, 1:]

        sorter = np.argsort(map_vertices)

        # mapping between Gs graph segments and triangle segments
        segments = sorter[np.searchsorted(map_vertices, seg, sorter=sorter)]

        if holes == []:
            C = {'vertices': vertices, 'segments': segments}
        else:
            C = {'vertices': vertices, 'segments': segments, 'holes': holes}

        T = triangle.triangulate(C, 'pa')

        if bplot:
            import triangle.plot as plottri
            ax = plt.gca()
            plottri(ax,**T)
            ax = plt.gca()
            ax.get_xaxis().set_visible(True)
            ax.get_yaxis().set_visible(True)
            plt.show()

        return T, map_vertices

    def buildGt(self, check=True, difftol=0.01, verbose=False, tqdmpos=0):
        """ build graph of convex cycles

        Parameters
        ----------

        check : boolean
        difftol : float
        verbose : boolean
        tqdmpos : progressbar

        TODO :
        - add an option to only take outside polygon
            => pass to self._triangle a hole corresponding to centroid of
            polygon except those of boundary ( see buildGtold )

        """

        # 1. Do a Delaunay triangulation
        #       build a list of triangle polygons : lTP
        #       vnodes refers to the nodes of Gs
        #       if vnodes == 0 it means this is a created
        #       segment which is tagged as _AIR
        ###

        Gtpbar = pbar(verbose, total=100., desc ='BuildGt',position=tqdmpos)
        pbartmp = pbar(verbose, total=100., desc ='Triangulation',leave=True,position=tqdmpos+1)

        logger.info('buildGt : Triangulation')
        #
        #
        #
        # T  dict  ['segment_markers', 'segments', 'holes', 'vertices', 'vertex_markers', 'triangles']
        T, map_vertices = self._triangle()

        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        # coordinates of triangle points
        # ptri np.array Ntri x  3 x 2

        ptri = T['vertices'][T['triangles']]

        # check that any point of triangulation belong to Gs
        #
        # for k in range(ptri.shape[0]):
        #     for l in range(3):
        #         point_tri = ptri[k,l,:]
        #         if self.ispoint(point_tri,tol=0.01)==0:
        #             print(point_tri)
        # pdb.set_trace()

        # List of Triangle Polygons
        pbartmp = pbar(verbose,total=100.,
                        desc ='Transfer polygons list',
                        leave=True,
                        position=tqdmpos+1)

        logger.info('buildGt : create list of Triangle Polygons')
        lTP = [ geu.Polygon(x) for x in ptri ]

        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        # update vnodes of Polygons
        pbartmp = pbar(verbose,total=100.,
                        desc ='Update Polygons vnodes',
                        leave=True,
                        position=tqdmpos+1)
        #
        # p is a polygon
        # get_points(p) : get points from polygon
        # this is for limiting the search region for large Layout
        #

        logger.info('buildGt : apply setvnodes_new on each polygon')
        [ polygon.setvnodes_new(self.get_points(polygon), self) for polygon in lTP ]

        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)


        # 2.add air walls to triangle poly
        ###
        # luaw  : list of tuples
        # ( polygon , array of _AIR segments)
        pbartmp = pbar(verbose,total=100.,
                        desc ='Build list of airwalls',
                        leave=True,
                        position=tqdmpos+1)

        # segments of the polygon which do not belong to Gs are airwalls (code 0)
        #
        logger.info('buildGt : get list of airwalls')
        luaw = [(p, np.where(p.vnodes == 0)[0]) for p in lTP]

        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        #
        # For a triangle polygon
        # creates new _AIR segments for segments not belonging to Gs
        # vnodes == 0
        #

        cpt = 1./(len(luaw)+1)
        _airseg = []

        pbartmp = pbar(verbose,total=100., desc ='Add airwalls',leave=True,position=tqdmpos+1)

        logger.info('buildGt : add %d new airwalls segments ',len(luaw))
        for p, uaw in luaw:
            # for each vnodes == 0, add an _AIR
            if verbose :
                pbartmp.update(100.*cpt)
            for aw in uaw:
                modpt = len(p.vnodes)
                _airseg.append(self.add_segment(p.vnodes[np.mod(aw - 1, modpt)],
                                                p.vnodes[np.mod(aw + 1, modpt)],
                                                name = '_AIR',
                                                z = (0, 40000000),
                                                verbose = False))
            # update polygon segments with new added airwalls
            p.setvnodes_new(self.get_points(p),self)

        if verbose:
            Gtpbar.update(100./12.)


        pbartmp = pbar(verbose,total=100., desc ='Update Graph',leave=True,position=tqdmpos+1)
        logger.info('buildGt : create temporary graph')

        # tri : np.array (Ntri x 3 )
        tri = T['triangles']
        nbtri = len(T['triangles'])
        # temporary name/node_index of triangles
        MT = -np.arange(1, nbtri + 1)

        # 3. Create a temporary graph
        #
        # where : positive nodes (>0) are triangles segments
        #         negative nodes (<0) are triangles centroids
        # edges link triangle centroids to their respective segments

        # Ex represent list of points in Gs corresponging to segments
        # [pt_head pt_tail]

        E0 = map_vertices[tri[:, 1:]]
        E1 = map_vertices[tri[:, :2]]
        E2 = map_vertices[tri[:, 0::2]]

        # from [pt_tail pt_head] get segment id in Gs

        n0 = [self.numseg(e[0], e[1]) for e in E0]
        n1 = [self.numseg(e[0], e[1]) for e in E1]
        n2 = [self.numseg(e[0], e[1]) for e in E2]

        # creation of a temporary graph

        G = nx.Graph()
        G.add_edges_from(zip(n0, MT))
        G.add_edges_from(zip(n1, MT))
        G.add_edges_from(zip(n2, MT))

        # 4. search in the temporary graph
        ###
        # nodes of degree 2  :
        # - they correspond to Gs segments that link to triangle centroid
        # - their neighbors are the triangles centroids

        # find nodes of degree 2 (corresponding to segments linked to a
        # triangle centroid)

        rn = []
        rn.extend([un for un in n0 if nx.degree(G, un) == 2])
        rn.extend([un for un in n1 if nx.degree(G, un) == 2])
        rn.extend([un for un in n2 if nx.degree(G, un) == 2])
        rn = np.unique(rn)

        # determine the neighbors of those segments (the 2 connected triangles
        # centroids)
        # v1.1 neigh = [nx.neighbors(G, un) for un in rn]
        #neigh = [ dict(G[un]).keys() for un in rn ]

        neigh = [[n for n in nx.neighbors(G,un)] for un in rn ]

        # store into networkx compliant format
        # rlrn = range len rn
        rlrn = range(len(rn))
        uE = [(neigh[un][0], neigh[un][1], {'segment': [ rn[un]] + self.Gs.nodes[rn[un]]['iso']}) for un in rlrn]
        iuE = {rn[un]: [-neigh[un][0], -neigh[un][1]] for un in rlrn }

        # delete temporary graph G
        del G

        logger.info('buildGt : creates graph Gt')
        # create graph Gt
        self.Gt = nx.Graph(name='Gt')
        self.Gt.add_edges_from(uE)
        self.Gt = nx.relabel_nodes(self.Gt, lambda x: -x)


        # add nodes ro Gt 
        # add polyg  to nodes
        # add indoor to nodes
        # add isopen to nodes
        #
        # WARNING indoor: True is weird

        nno = [(n, {'polyg': lTP[n - 1], 'indoor':True, 'isopen':True})
               for n in self.Gt.nodes()]

        self.Gt.add_nodes_from(nno)
        self.Gt.pos = {}

        self.Gt.pos.update({n: np.array(
            self.Gt.nodes[n]['polyg'].centroid.xy).squeeze() for n in self.Gt.nodes()})

        #fig = plt.figure(figsize=(50,50))
        #self.showG('st',fig=fig)

        # self.Gtpos = {-MT[i]:pMT[i] for i in xrange(len(MT))}
        # plt.figure()
        # # G=nx.Graph()
        # # G.add_edges_from(E0)
        # # G.add_edges_from(E1)
        # # G.add_edges_from(E2)

        _airseg = np.array(_airseg)
        _airseg = _airseg[_airseg != np.array(None)].astype('int')
        _airseg = np.unique(_airseg)

        logger.info('buildGt : start Mikado like simplification ')
        #
        # Mikado like progression for simplification of a set of convex polygons
        #
        #    Loop over AIR segments
        #

        mapoldcy = {c: c for c in self.Gt.nodes()}

        # self.showG('st',aw=1)

        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        Nairseg = len(_airseg)
        cpt = 1./(Nairseg+1)
        pbartmp = pbar(verbose,total=100., desc ='Mikado',leave=True,position=tqdmpos+1)
        for a in _airseg:
            if verbose:
                pbartmp.update(100.*cpt)
            #
            # n0,n1 : cycle number
            #
            #debug = False

            n0, n1 = iuE[a]
            logger.debug(" segment %d : %d %d ",a,n0,n1)
            found = False
            while not found:
                nn0 = mapoldcy[n0]
                if n0 == nn0:
                    found = True
                else:
                    n0 = nn0

            found = False
            while not found:
                nn1 = mapoldcy[n1]
                if n1 == nn1:
                    found = True
                else:
                    n1 = nn1

            # if the 2 cycles are distincts
            if (n0 != n1):
                p0 = self.Gt.nodes[n0]['polyg']
                p1 = self.Gt.nodes[n1]['polyg']

                # Merge polygon
                P = p0 + p1
                # If the new Polygon is convex update Gt
                #
                if geu.isconvex(P):
                    logger.debug(" merge %d : %d %d ",a,n0,n1)
                    # updates vnodes of the new merged polygon
                    P.setvnodes_new(self.get_points(P),self)
                    # update edge
                    n0s = n0
                    n1s = n1
                    # get segments information from cycle n0
                    dne = dict(self.Gt[n0])
                    # remove connection to n0 to avoid a cycle being
                    # connected to itself
                    # v1.1 self.Gt[n1].pop(n0)
                    try:
                        dict(self.Gt[n1]).pop(n0)
                    except:
                        pdb.set_trace()
                    # add information from adjacent cycle n1
                    dne.update(dict(self.Gt[n1]))
                    # list of items of the merged dictionnary
                    ine = dne.items()
                    # update n0 with the new merged polygon
                    self.Gt.add_node(n0, polyg=P)
                    # connect new cycle n0 to neighbors
                    # for x in ine:
                    #     if x[0]!=n0:
                    #         ncy  = x[0]
                    #         dseg = x[1]
                    #         # a link between cycles already exists
                    #         if self.Gt.has_edge(n0,ncy):
                    #             dseg_prev = self.Gt.edge[n0][ncy]
                    #             dseg['segment']=list(set(dseg['segment']+dseg_prev['segment']))
                    #         printn0,ncy,dseg['segment']
                    #         self.Gt.add_edge(n0,ncy,segment=dseg['segment'])

                    self.Gt.add_edges_from([(n0, x[0], x[1])
                                            for x in ine if x[0] != n0])
                    # remove old cycle n1 n
                    self.Gt.remove_node(n1)
                    # update pos of the cycle with merged polygon centroid
                    self.Gt.pos[n0] = np.array((P.centroid.xy)).squeeze()
                    self.Gt.pos.pop(n1)
                    # delete _air segment a
                    # do not apply g2npy
                    self.del_segment(a, verbose=False, g2npy=False)
                    mapoldcy[n1] = n0
                    # fig,a=self.showG('st',aw=1)
                    # plt.show()
            else:
                    self.del_segment(a, verbose=False, g2npy=False)
        ######
        # fix renumbering Gt nodes

        if verbose:
            Gtpbar.update(100./12.)

        pbartmp = pbar(verbose,total=100., desc ='Update Gs ncy',leave=True,position=tqdmpos+1)

        pos = self.Gt.pos
        nl = {c: uc + 1 for uc, c in enumerate(self.Gt.nodes())}
        self.Gt = nx.relabel_nodes(self.Gt, nl)
        self.Gt.pos = {}
        self.Gt.pos = {nl[n]: pos[n] for n in nl}

        self._updGsncy()
        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)
        #
        # add cycle 0 to boundaries segments
        # cycle 0 is necessarily outdoor
        #
        self.Gt.add_node(0, indoor=False)
        for s in self.segboundary:
            self.Gs.nodes[s]['ncycles'].append(0)

        #
        # boundary adjascent cycles
        #
        #adjcyair = np.array(map(lambda x: filter(lambda y: y != 0,
        #                                         self.Gs.nodes[x]['ncycles'])[0], self.segboundary))
        adjcyair = np.array([[n for n in self.Gs.nodes[s]['ncycles'] if n!=0]
                             for s in self.segboundary]).ravel()
        # connect cycles separated by air wall to cycle 0
        for cy, seg in zip(adjcyair, self.segboundary):
            self.Gt.nodes[cy]['indoor'] = False
            self.Gt.nodes[cy]['isopen'] = True
            self.Gt.add_edge(0, cy, segment=[seg])

        #
        #
        #
        if check:
            # print("check len(ncycles) == 2",)
            nodes = [i for i in self.Gs.nodes() if i > 0]
            cncy = np.array([len(self.Gs.nodes[i]['ncycles']) for i in nodes])
            ucncyl = np.where(cncy < 2)[0]
            ucncym = np.where(cncy > 2)[0]
            assert len(ucncyl) == 0, "Some segments are connected to LESS than 2 cycles" + \
                str(np.array(nodes)[ucncyl])
            assert len(ucncym) == 0, "Some segments are connected to MORE than 2 cycles" + \
                str(np.array(nodes)[ucncym])
            # print("passed")

        # self.degree is updated in g2npy
        # self.degree has to be called before determination of diffraction points
        # which relies of the full determination of the degree of each point of Gs
        # including the corner point with degree 0 ( only connected to _AIR)

        self.g2npy()
        # find diffraction points : updating self.ddiff
        tqdmkwargs={'total':100.,'desc':'Find Diffractions','position':1}
        self._find_diffractions(difftol=difftol,verbose=verbose,tqdmkwargs=tqdmkwargs)
        if verbose:
            Gtpbar.update(100./12.)
            # print('find diffraction...Done 8/12')
            pbartmp = pbar(verbose,total=100., desc ='Diffraction on airwalls',leave=True,position=tqdmpos+1)

        # 
        # explanation of lnss
        #
        # list of diffraction point involving different segment 
        # list of diffraction point involving subsegment ( = iso segments)
        # needs checking height in rays.to3D for constructing the 3D ray
        #

        self.lnss = [x for x in self.ddiff if len(set(self.Gs[x]).intersection(set(self.lsss))) > 0]
        #set(nx.neighbors(self.Gs, x)).intersection(set(self.lsss))) > 0]


        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)
        #
        #   VIII -  Construct the list of interactions associated to each cycle
        #
        # Interaction labeling convention
        #
        #   tuple (npoint,)  : Diffraction on point npoint
        #   tuple (nseg,ncycle) : Reflection on nseg toward cycle ncycle
        #   tuple (nseg,cy0,cy1) : Transmission from cy0 to cy1 through nseg
        #
        #   At that stage the diffraction points are not included
        #   not enough information available.
        #   The diffraction points are not known yet
        tqdmkwargs={'total':100.,'desc':'List of interactions','position':1}

        self._interlist(verbose=verbose,tqdmkwargs=tqdmkwargs)
        if verbose:
            Gtpbar.update(100./12.)

        #
        # dca : dictionnary of cycles which have an air wall
        #
        pbartmp = pbar(verbose,total=100., desc ='Build dca',leave=True,position=tqdmpos+1)

        self.dca = {}
        for seg, d in self.Gs.node.items():
            if seg > 0:
                if ((d['name'] == 'AIR') or d['name'] == '_AIR'):
                    cy = d['ncycles']
                    try:
                        self.dca[cy[0]].append(cy[1])
                    except:
                        self.dca[cy[0]] = [cy[1]]
                    try:
                        self.dca[cy[1]].append(cy[0])
                    except:
                        self.dca[cy[1]] = [cy[0]]
        if verbose:
            # print('build dca...Done 11/12')
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        #
        # indoor property is spread by contagion
        #
        pbartmp = pbar(verbose,total=100., desc ='Indoor properties',leave=False,position=tqdmpos+1)


        visited = [0]
        #v1.1 to_visit = nx.neighbors(self.Gt, 0)
        to_visit = list(dict(self.Gt[0]).keys())

        law = self.name['_AIR'] + self.name['AIR']
        while len(to_visit) > 0:
            # get current cycle

            cur_cy = to_visit.pop()
            # get neighbors of current_cycle
            # v1.1 neighbors = nx.neighbors(self.Gt, cur_cy)
            neighbors = self.Gt[cur_cy].keys()
            # get neighbors separated by an air_wall
            neighbors_aw = [x for x in neighbors 
                            if (len(self.Gt[cur_cy][x]['segment'])==1 and
                                self.Gt[cur_cy][x]['segment'][0] in law
                                )
                            ]
            # get not visited neighbors_aw
            nv_neighbors_aw = [
                x for x in neighbors_aw if x not in (visited + to_visit)]
            # not visited neighbors air wall separated cycles are outdoor cycle
            for x in nv_neighbors_aw:
                self.Gt.nodes[x]['indoor'] = False
                self.Gt.nodes[x]['isopen'] = True
            # extend to_visit to not visited neighbors
            to_visit.extend(nv_neighbors_aw)
            visited.append(cur_cy)
        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        self.g2npy()

    def _visual_check(self,fontsize=18):
        """ visual checking of graphs

        Parameters
        ----------

        fontsize : int

        """
        fig, axs = plt.subplots(2, 2,figsize=(10,10))
        plt.subplots_adjust(left = 0 ,
                            right = 1.0, 
                            bottom = 0 ,
                            top = 1 ,
                            wspace = 0 , 
                            hspace =0)


        if hasattr(self,'Gs') and hasattr(self,'Gt'):
            ax = axs[0, 0]

            self.showG('s', aw=1, ax=ax, fig=fig)

            indoor = [self.Gt.nodes[p]['polyg']
                      for p in self.Gt.nodes() if p != 0 and self.Gt.nodes[p]['indoor']]
            outdoor = [self.Gt.nodes[p]['polyg']
                      for p in self.Gt.nodes()  if p != 0 and not self.Gt.nodes[p]['indoor']]

            self.pltpoly(indoor, color='r', ax=ax, fig=fig)
            self.pltpoly(outdoor, color='g', ax=ax, fig=fig)

            ax = axs[0, 1]
            f, ax = self.showG('s', aw=1, ax=ax, fig=fig)

            if hasattr(self,'ddiff'):
                diffpos = np.array([self.Gs.pos[x] for x in self.ddiff.keys()])
                ax.scatter(diffpos[:, 0], diffpos[:, 1],s=130)
            #ax.set_title('Diffraction points')

            ax = axs[1, 0]
            f, ax = self.showG('st', aw=1, ax=ax, fig=fig)
            #ax.set_title('$\mathcal{G}_t$',fontsize=fontsize)
            ax.set_axis_off

        if hasattr(self,'Gv'):
            ax = axs[1, 1]
            f, ax = self.showG('sv', aw=1, ax=ax, fig=fig)
            #ax.set_title('$\mathcal{G}_v$',fontsize=fontsize)
            ax.set_axis_off
        else:
            print('no Gv found. Yet computed ?')
        plt.savefig('visual_check.pdf')
        #plt.tight_layout()
        # axs[2,1].remove()

    def _delaunay(self, poly, polyholes=[]):
        """ make a Delaunay partitioning of a polygon

            If polyhole == []


                if a cycle is non convex

                1- find its polygon
                2- partition polygon into convex polygons (Delaunay)
                3- try to merge partitioned polygons in order to obtain
                   the minimal number of convex polygons


            If polyholes != []

                polygon poly contains holes (polyholes)

            This methods returns a partitioning of the polygon poly 
            into several convex polygons (voronoi). 

        Parameters
        ----------

            poly : sh.Polygon
            polyhole : list of sh.Polygon


        Returns
        -------
            ncpol : list
                list of new created geu.Polygons

        Notes
        -----

        The algorithm updates the Gt nodes and edges created into self.buildGt
        by adding new nodes and new _AIR segments.

        Called In 
        ---------

        pylayers.gis.layout.buildGt


        See Also
        --------

        pylayers.gis.layout.buildGt
        pylayers.gis.layout.add_segment
        pylayers.gis.layout.del_segment
        pylayers.util.geomutil.Polygon
        sp.spatial.Delaunay

        """

        pucs = np.array(poly.exterior.xy).T

        # keep all convex points (in + out) to build a Delaunay triangulation
        if polyholes != []:
            if not isinstance(polyholes, list):
                polyholes = [polyholes]
            for ph in polyholes:
                # sum up polyholes to their gathered polygones
                pucsh = np.array(ph.exterior.xy).T
                pucs = np.vstack((pucs, pucsh))

        if len(pucs) != 0:
            ####
            # perform a Delaunay Partioning
            ####

            trid = sp.spatial.Delaunay(pucs)
            tri = trid.simplices
            polys = []
            naw = []
            popo = []

            for t in tri:
                ts = geu.Polygon(pucs[t])
                # check if the new polygon is contained into
                # the original polygon (non guarantee by Delaunay)
                try:
                    C0 = poly.contains(ts)
                except:
                    from IPython.core.debugger import Tracer
                    Tracer()()
                if polyholes == []:
                    C = [False]
                    I = 0
                else:
                    C = [isinstance(ii.intersection(ts), sh.Polygon)
                         for ii in polyholes]

                popo.append(ts)
                # if poly contains triangle but not the polyholes
                # if polyholes !=[]:
                #     self.pltpoly([ts],color='b')
                #     import ipdb
                #     ipdb.set_trace()
                if C0 and (not np.any(C)):
                    # if polyholes!=[]:
                    #     self.pltpoly([ts],color='r')
                    #     plt.draw()

                    cp = ts
                    cp.setvnodes(self)
                    uaw = np.where(cp.vnodes == 0)[0]
                    lvn = len(cp.vnodes)
                    for i in uaw:
                        # keep track of created airwalls, because some
                        # of them will be destroyed in step 3.
                        naw.append(self.add_segment(
                                   cp.vnodes[np.mod(i - 1, lvn)],
                                   cp.vnodes[np.mod(i + 1, lvn)], name='_AIR'))
                    polys.append(cp)
            #
            # 3. merge Delaunay triangulation in order to obtain
            #   the larger convex polygons partitioning
            #
            diff = poly.difference(sh.MultiPolygon(polys))
            if isinstance(diff, sh.Polygon):
                diff = sh.MultiPolygon([diff])
            if isinstance(diff, sh.MultiPolygon):
                for d in diff:
                    extra = geu.Polygon(d)
                    extra.setvnodes(self)
                    polys.append(extra)

            cpolys = []
            nbpolys = len(polys)

            while polys != []:
                p = polys.pop(0)
                for ip2, p2 in enumerate(polys):
                    conv = False
                    inter = p.intersection(p2)
                    # if 2 triangles have a common segment
                    pold = p
                    if isinstance(inter, sh.LineString):
                        p = p + p2
                        if p.isconvex():
                            polys.pop(ip2)
                            polys.insert(0, p)
                            conv = True
                            break
                        else:
                            # if pold not in cpolys:
                            #     cpolys.append(pold)
                            p = pold
                # if (ip2 >= len(polys)):# and (conv):
                # if conv :
                #     if p not in cpolys:
                #         cpolys.append(p)
                if not conv:  # else:
                    if pold not in cpolys:
                        cpolys.append(pold)
                if len(polys) == 0:
                    cpolys.append(p)

            # 4. ensure the correct vnode numerotation of the polygons
            # and remove unecessary airwalls

            # ncpol : new created polygons
            ncpol = []
            vnodes = []
            for p in cpolys:
                interpoly = poly.intersection(p)
                if isinstance(interpoly, sh.MultiPolygon):
                    raise AttributeError('multi polygon encountered')
                else:
                    try:
                        ptmp = geu.Polygon(interpoly)
                        # ptmp = self.polysh2geu(interpoly)
                    except:
                        import ipdb
                        ipdb.set_trace()

                ptmp.setvnodes(self)
                ncpol.append(ptmp)
                vnodes.extend(ptmp.vnodes)

            # if no polyholes
            if polyholes == []:
                # 4bis
                # Check if all the original area is covered
                # sometimes, area surrounded by 2 new airwalls is not found
                # the following code re-add it.
                cpdiff = poly.difference(cascaded_union(cpolys))
                if isinstance(cpdiff, sh.Polygon):
                    cpdiff = sh.MultiPolygon([cpdiff])
                if isinstance(cpdiff, sh.MultiPolygon):
                    for cp in cpdiff:
                        ptmp = geu.Polygon(cp)
                        ptmp.setvnodes(self)
                        ncpol.append(ptmp)
                        vnodes.extend(ptmp.vnodes)

            daw = filter(lambda x: x not in vnodes, naw)

            for d in daw:
                self.del_segment(d, verbose=False, g2npy=False)
            self.g2npy()

        return ncpol

    def _updGsncy(self):
        """ update Gs ncycles using Gt information

        Update graph Gs segment with their 2 cycles information

        initialize a void list 'ncycles' for each segment of Gs

        See Also
        --------

        pylayers.gis.layout.buildGt
        pylayers.gis.layout.convexify

        """

        for k in self.Gs.node:
            self.Gs.nodes[k]['ncycles'] = []

        # filter out node 0
        Gtnodes = filter(lambda x: x != 0, self.Gt.nodes())

        # loop over all cycles
        for ncy in Gtnodes:
            # get vnodes : points and segments number
            vnodes = self.Gt.nodes[ncy]['polyg'].vnodes
            for n in vnodes:
                if n == 0:
                    pdb.set_trace()
                if ncy not in self.Gs.nodes[n]['ncycles']:
                    self.Gs.nodes[n]['ncycles'].append(ncy)
                    if n > 0:
                        if len(self.Gs.nodes[n]['ncycles']) > 2:
                            print(n, self.Gs.nodes[n]['ncycles'])
                            logger.warning(
                                'A segment cannot relate more than 2 cycles')

        for nseg in self.Gs.node:
            if nseg > 0:
                ncycles = self.Gs.nod[nseg]['ncycles']
                if len(ncycles) > 1:
                    #if nseg not in self.Gt.edge[ncycles[0]][ncycles[1]]['segment']:
                    #    self.Gt.edge[ncycles[0]][ncycles[1]][
                    #        'segment'].append(nseg)
                    if nseg not in self.Gt[ncycles[0]][ncycles[1]]['segment']:
                        self.Gt[ncycles[0]][ncycles[1]]['segment'].append(nseg)

    def _addoutcy(self, check=False):
        """ 
        Probably use in a future version of buildGt , managing the upcoming inifile
        add outside cycle (absorbant region index 0 )

        Parameters
        ----------

        check : Boolean 

        #   if ncycles is a list which has only one element then the adjascent
        #   cycle is the  outside region (cycle 0)

        """
        seg0 = []
        for macvx in self.macvx:
            seg = [i for i in macvx.vnodes if i > 0]
            seg0 = seg0 + seg
        [self.Gs.node[i]['ncycles'].append(0) for i in seg0]
        if check:
            print("check len(ncycles) == 2",)
            nodes = [i for i in self.Gs.nodes() if i > 0]
            cncy = np.array([len(self.Gs.nodes[i]['ncycles']) for i in nodes])
            ucncyl = np.where(cncy < 2)[0]
            ucncym = np.where(cncy > 2)[0]
            assert len(ucncyl) == 0, "Some segments are connected to LESS than 2 cycles" + \
                str(np.array(nodes)[ucncyl])
            assert len(ucncym) == 0, "Some segments are connected to MORE than 2 cycles" + \
                str(np.array(nodes)[ucncym])
            print("passed")

    def _interlist(self, nodelist=[],verbose = False,tqdmkwargs={}):
        """ Construct the list of interactions associated to each cycle


        Parameters
        ----------

        nodelist: list
            list of Gt nodes (cycles) for which interactions have to be found


        Notes
        -----

        if selfr.indoor==True , get list of interaction of Gt cycle with indoor =True
            else list of indoor interaction is skipped

         Interaction labeling convention

           tuple (npoint,)  : Diffraction on point npoint
           tuple (nseg,ncycle) : Reflection on nseg toward cycle ncycle
           tuple (nseg,cy0,cy1) : Transmission from cy0 to cy1 through nseg

           At that stage the diffraction points are not included
           not enough information available. The diffraction point are not
           known yet

        See Also
        --------

        pylayers.gis.layout.buildGt
        pylayers.gis.layout._convex_hull

        """

        if tqdmkwargs=={}:
            tqdmkwargs={'total':100.,
                        'desc':'list of interactions',
                        'position':0}

        if nodelist == []:
            nodelist = self.Gt.nodes()
        elif not isinstance(nodelist, list):
            nodelist = [nodelist]

        # for all cycles k (node of Gt)
        if verbose :
            cpt = 1./(len(nodelist)+1.)
            pbar =  tqdm.tqdm(tqdmkwargs)
        for k in nodelist:
            if verbose:
                pbar.update(100.*cpt)
            if k != 0:
                if self.typ=='indoor' or not self.Gt.nodes[k]['indoor']:
                    #vnodes = self.Gt.nodes[k]['vnodes']
                    vnodes = self.Gt.nodes[k]['polyg'].vnodes
                    ListInteractions = []
                    for inode in vnodes:
                        if inode > 0:   # segments
                            cy = set(self.Gs.nodes[inode]['ncycles'])
                            name = self.Gs.nodes[inode]['name']  # segment name
                            #
                            # Reflexion occurs on segment different
                            # from AIR and ABSORBENT  (segment number, cycle)
                            #
                            if ((name != '_AIR') & (name != 'AIR') & (name != 'ABSORBENT')):
                                ListInteractions.append((inode, k))
                            #
                            # Transmission requires 2 cycles separated by a
                            # segment which is different from METAL and ABSORBENT
                            #
                            # (segment number, cycle in , cycle out )
                            if len(cy) == 2:
                                if ('METAL' not in name) & ('ABSORBENT' not in name):
                                    ncy = list(cy.difference({k}))[0]
                                    ListInteractions.append((inode, k, ncy))
                                    ListInteractions.append((inode, ncy, k))
                        else:  # points
                            pass
                    # add list of interactions of a cycle
                    self.Gt.add_node(k, inter=ListInteractions)
                else:
                    self.Gt.add_node(k, inter=[])

    def _convex_hull(self, mask):
        """
        Add air walls to the layout enveloppe in self.Gs 
        in order the hull of the Layout to be convex.

        Parameters
        ----------

        mask : Polygon

        Returns
        -------

        polys : list of geu.Polygon
            nsew polygon of the convex hull

        self.macvx : convex mask of the layout

        Notes
        -----

        This is a post processing of BuildGt

        See Also
        --------

        pylayers.gis.layout._interlist


        """

        # 1 - Find differences between the convex hull and the Layout contour
        #     The result of the difference are polygons

        masku = cascaded_union(mask)
        ch = masku.convex_hull
        P = ch.difference(masku)
        polys = []
        if isinstance(P, sh.MultiPolygon):
            for p in P:
                if p.area > 1e-3:
                    polys.append(geu.Polygon(p))
                    polys[-1].setvnodes(self)

        lncy = []
        for p in polys:
            # p.coorddeter()
            uaw = np.where(p.vnodes == 0)
            for aw in uaw:
                # 2 - non existing segments are created as airwalls
                awid = self.add_segment(
                    p.vnodes[aw - 1][0], p.vnodes[aw + 1][0], name='AIR')
                p.vnodes[aw] = awid

        # U = cascaded_union([mask]+polys)

        # self.macvx = geu.Polygon(U)
        # self.macvx.setvnodes(self)

        return polys


    def buildGv(self, show=False,verbose=False,tqdmpos=0):
        """ build visibility graph

        Parameters
        ----------

        show : boolean
            default False
        verbose : boolean
        tqdmpos : progressbar

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('TA-Office.lay')
        >>> L.buildGt()
        >>> Ga = L.buildGr()
        >>> L.buildGv()

        Notes
        -----

        This method exploits cycles convexity.

        """
        if not hasattr(self,'ddiff'):
            self.ddiff={}
        Gvpbar = pbar(verbose,total=100., desc ='build Gv',position=tqdmpos)

        self.Gv = nx.Graph(name='Gv')
        #
        # loop over convex cycles (nodes of Gt)
        #
        self.dGv = {}  # dict of Gv graph

        cpt = 1./(len(self.Gt.node) + 1.)

        for icycle in self.Gt.node:
            if verbose:
                Gvpbar.update(100.*cpt)
            if icycle != 0:
                #if self.indoor or not self.Gt.nodes[icycle]['indoor']:
                    #print(icycle)
                #    pass
                #
                #  If indoor or outdoor all visibility are calculated
                #  If outdoor only visibility between iso = 'AIR' and '_AIR' are calculated 
                #
                #if self.indoor or not self.Gt.nodes[icycle]['indoor']:
                polyg = self.Gt.nodes[icycle]['polyg']

                # plt.show(polyg.plot(fig=plt.gcf(),ax=plt.gca())

                # take a single segment between 2 points 

                vnodes = polyg.vnodes

                # list of index of points in vodes
                unodes = np.where(vnodes<0)[0]

                # list of position of an incomplete list of segments 
                # used rule : after a point there is always a segment 
                useg = np.mod(unodes+1,len(vnodes))

                # list of points 
                #npt  = filter(lambda x: x < 0, vnodes)
                npt = [ x for x in vnodes if x <0 ]

                nseg_full = [x for x in vnodes if x > 0]
                # nseg : incomplete list of segments
                #
                # if mode outdoor and cycle is indoor only 
                # the part above the building (AIR and _AIR) is considered
                if ((self.typ=='outdoor') and (self.Gt.nodes[icycle]['indoor'])):
                    nseg = [ x for x in nseg_full if
                            ((self.Gs.nodes[x]['name']=='AIR') or
                                (self.Gs.nodes[x]['name']=='_AIR') ) ]
                else:
                    nseg = vnodes[useg]


                # # nseg_full : full list of segments
                # #nseg_full = filter(lambda x: x > 0, vnodes)

                # # keep only airwalls without iso single (_AIR)
                # nseg_single = filter(lambda x: len(self.Gs.nodes[x]['iso'])==0, nseg)

                # lair1 = self.name['AIR'] 
                # lair2 = self.name['_AIR']
                # lair  = lair1 + lair2

                # # list of airwalls in nseg_single

                # airwalls = filter(lambda x: x in lair, nseg_single)

                # diffraction points 

                ndiff = [x for x in npt if x in self.ddiff.keys()]
                #
                # Create a graph
                #

                Gv = nx.Graph(name='Gv')
                #
                # in convex case :
                #
                #    i)  every non aligned segments see each other
                #
                for nk in combinations(nseg, 2):
                    nk0 = self.tgs[nk[0]]
                    nk1 = self.tgs[nk[1]]
                    tahe0 = self.tahe[:, nk0]
                    tahe1 = self.tahe[:, nk1]

                    pta0 = self.pt[:, tahe0[0]]
                    phe0 = self.pt[:, tahe0[1]]
                    pta1 = self.pt[:, tahe1[0]]
                    phe1 = self.pt[:, tahe1[1]]

                    aligned = geu.is_aligned4(pta0,phe0,pta1,phe1)
                    # A0 = np.vstack((pta0, phe0, pta1))
                    # A0 = np.hstack((A0, np.ones((3, 1))))

                    # A1 = np.vstack((pta0, phe0, phe1))
                    # A1 = np.hstack((A1, np.ones((3, 1))))

                    # d0 = np.linalg.det(A0)
                    # d1 = np.linalg.det(A1)

                    #if not ((abs(d0) < 1e-1) & (abs(d1) < 1e-1)):
                    if not aligned:
                        if ((0 not in self.Gs.nodes[nk[0]]['ncycles']) and
                            (0 not in self.Gs.nodes[nk[1]]['ncycles'])):
                            # get the iso segments of both nk[0] and nk[1]
                            if ((self.typ=='indoor') or (not self.Gt.nodes[icycle]['indoor'])):
                                l0 = [nk[0]]+self.Gs.nodes[nk[0]]['iso']
                                l1 = [nk[1]]+self.Gs.nodes[nk[1]]['iso']
                            else:
                                l0 = [nk[0]]
                                l1 = [nk[1]]

                            for vlink in product(l0,l1):
                                #printicycle,vlink[0],vlink[1]
                                Gv.add_edge(vlink[0], vlink[1])

                #
                # Handle diffraction points
                #
                #    ii) all non adjascent valid diffraction points see each other
                #    iii) all valid diffraction points see segments non aligned
                #    with adjascent segments
                #
                #if diffraction:
                #
                # diffraction only if indoor or outdoor cycle if outdoor
                # 
                if ((self.typ=='indoor') or (not self.Gt.nodes[icycle]['indoor'])):
                    ndiffvalid = [ x for x in ndiff if icycle in self.ddiff[x][0]]

                        # non adjascent segment of vnodes see valid diffraction
                        # points
                    for idiff in ndiffvalid:
                        #
                        # segments voisins du point de diffraction valide
                        #
                        # v1.1 nsneigh = [x for x in 
                        #           nx.neighbors(self.Gs, idiff) 
                        #           if x in nseg_full]
                        nsneigh = [x for x in self.Gs[idiff] if x in nseg_full]
                        # segvalid : not adjascent segment
                        seen_from_neighbors = []

                        #
                        # point to point
                        #
                        for npoint in ndiffvalid:
                            if npoint != idiff:
                                Gv.add_edge(idiff, npoint)

                        #
                        # All the neighbors segment in visibility which are not connected to cycle 0
                        # and which are not neighbors of the point idiff
                        #
                        for x in nsneigh:
                            # v1.1 neighbx = [ y for y in nx.neighbors(Gv, x) 
                            #            if 0 not in self.Gs.nodes[y]['ncycles'] 
                            #            and y not in nsneigh]
                            neighbx = [ y for y in Gv[x] 
                                        if 0 not in self.Gs.nodes[y]['ncycles'] 
                                        and y not in nsneigh]
                            seen_from_neighbors += neighbx

                        for ns in seen_from_neighbors:
                            Gv.add_edge(idiff, ns)

                #
                # Graph Gv composition
                #

                self.Gv = nx.compose(self.Gv, Gv)
                self.dGv[icycle] = Gv

    def buildGi(self,verbose=False,tqdmpos=0):
        """ build graph of interactions

        Notes
        -----

        For each node of graph Gv creates
        5 different nodes associated to the same segment

        (np,) D
        (ns,cy0) R -> cy0
        (ns,cy1) R -> cy1
        (ns,cy0,cy1) T 0->1
        (ns,cy1,cy0) T 1->0

        Gi is an oriented Graph (DiGraph) 

        """

        Gipbar = pbar(verbose,total=100., desc ='Build Gi',position=tqdmpos)
        if verbose:
            Gipbar.update(0.)

        self.Gi = nx.DiGraph(name='Gi')
        self.Gi.pos = {}

        #
        # 1 ) Create nodes of Gi and their positions
        #
        # diffraction node  (D,)
        # reflexion node    (R,cy0)
        # transmission node (T,cy0,cy1)
        #

        cpt = 100./(len(self.Gv.node)+1)
        pbartmp = pbar(verbose,total=100., desc ='Create Gi nodes',position=tqdmpos+1)

        for n in self.Gv.node:
            if verbose:
                pbartmp.update(cpt)

            if n < 0:  # D
                self.Gi.add_node((n,))
                self.Gi.pos[(n,)] = self.Gs.pos[n]
            if n > 0:  # R | T
                cy = self.Gs.nodes[n]['ncycles']
                name = self.Gs.nodes[n]['name']
                assert(len(cy) == 2)
                cy0 = cy[0]
                cy1 = cy[1]

                #nei = self.Gs.neighbors(n)  # get neighbor
                nei = list(dict(self.Gs[n]).keys())  # get neighbor
                np1 = nei[0]
                np2 = nei[1]

                p1 = np.array(self.Gs.pos[np1])
                p2 = np.array(self.Gs.pos[np2])
                l = p1 - p2
                nl = np.dot(l, l)
                ln = l / nl

                delta = nl / 10.

                # On AIR or ABSORBENT there is no reflection

                if ((name != '_AIR') & (name != 'AIR') & (name != 'ABSORBENT')):
                    self.Gi.add_node((n, cy0))
                    self.Gi.pos[(n, cy0)] = tuple(self.Gs.pos[n] + ln * delta)
                    self.Gi.add_node((n, cy1))
                    self.Gi.pos[(n, cy1)] = tuple(self.Gs.pos[n] - ln * delta)

                # Through METAL or ABSORBENT there is no transmission
                # except if n has a subsegment

                if (name != 'METAL') & (name != 'ABSORBENT'):
                    self.Gi.add_node((n, cy0, cy1))
                    self.Gi.add_node((n, cy1, cy0))
                    self.Gi.pos[(n, cy0, cy1)] = tuple( self.Gs.pos[n] + ln * delta / 2.)
                    self.Gi.pos[(n, cy1, cy0)] = tuple( self.Gs.pos[n] - ln * delta / 2.)

        #
        # 2) Establishing link between interactions
        #
        # Loop over all Gt nodes cy 
        #
        #   if cy > 0 
        #     calculates vnodes of cycles
        #     for all node of vnodes
        #
        iprint = 0 
        if verbose :
            Gipbar.update(33.)

        cpt = 100./(len(self.Gt.node)+1)
        pbartmp = pbar(verbose,total=100., desc ='Create Gi nodes',position=tqdmpos+1)


        for cy in self.Gt.node:
            if verbose:
                pbartmp.update(cpt)
            # for all >0 convex cycles
            if cy > 0:
                vnodes = self.Gt.nodes[cy]['polyg'].vnodes
                npt = []
                #
                # find all diffraction points involved in the cycle cy 
                #
                for x in vnodes:
                    if x < 0:
                        if x in self.ddiff:
                            for y in self.ddiff[x][0]:
                                if y == cy:
                                    npt.append(x)

                nseg = [ k for k in vnodes if k>0 ]
                # all segments and diffraction points of the cycle
                vnodes = nseg + npt

                for nstr in vnodes:

                    if nstr in self.Gv.nodes():
                        # list 1 of interactions

                        li1 = []
                        if nstr > 0:
                            # output cycle 
                            # cy -> cyo1 
                            cyo1 = self.Gs.nodes[nstr]['ncycles']
                            cyo1 = [ x for x in cyo1 if x!= cy] [0]
                            #cyo1 = filter(lambda x: x != cy, cyo1)[0]

                            # R , Tin , Tout
                            if cyo1 > 0:
                                if (nstr, cy) in self.Gi.nodes():
                                    li1.append((nstr, cy))  # R 
                                if (nstr, cy, cyo1) in self.Gi.nodes():
                                    li1.append((nstr, cy, cyo1)) # T cy -> cyo1 
                                if (nstr, cyo1, cy) in self.Gi.nodes():
                                    li1.append((nstr, cyo1, cy)) # T : cyo1 -> cy 
                                # if (nstr,cy) in self.Gi.nodes():
                                #     li1 = [(nstr,cy),(nstr,cy,cyo1),(nstr,cyo1,cy)]
                                # else:# no reflection on airwall
                                #     li1 = [(nstr,cyo1,cy)]
                            else:
                                if (nstr, cy) in self.Gi.nodes():
                                    li1 = [(nstr, cy)]
                                # else:
                                #     li1 =[]
                        else:
                            # D
                            li1 = [(nstr,)]
                        # list of cycle entities in visibility of nstr
                        # v1.1 lneighb = nx.neighbors(self.Gv, nstr)
                        lneighb = list(dict(self.Gv[nstr]).keys())
                        #if (self.Gs.nodes[nstr]['name']=='AIR') or (
                        #        self.Gs.nodes[nstr]['name']=='_AIR'):
                        #    lneighcy = lneighb
                        #else:
                        # list of cycle entities in visibility of nstr in the same cycle 
                        lneighcy = [ x for x in lneighb if x in vnodes ] 
                        # lneighcy = filter(lambda x: x in vnodes, lneighb)

                        for nstrb in lneighcy:
                            if nstrb in self.Gv.nodes():
                                li2 = []
                                if nstrb > 0:
                                    cyo2 = self.Gs.nodes[nstrb]['ncycles']
                                    cyo2 = [ x for x in cyo2 if x!= cy] [0]
                                    #cyo2 = filter(lambda x: x != cy, cyo2)[0]
                                    if cyo2 > 0:
                                        if (nstrb, cy) in self.Gi.nodes():
                                            li2.append((nstrb, cy))
                                        if (nstrb, cy, cyo2) in self.Gi.nodes():
                                            li2.append((nstrb, cy, cyo2))
                                        if (nstrb, cyo2, cy) in self.Gi.nodes():
                                            li2.append((nstrb, cyo2, cy))
                                        # if (nstrb,cy) in self.Gi.nodes():
                                        #     li2 = [(nstrb,cy),(nstrb,cy,cyo2),(nstrb,cyo2,cy)]
                                        # else: #no reflection on airwall
                                        #     li2 = [(nstrb,cy,cyo2),(nstrb,cyo2,cy)]
                                    else:
                                        if (nstrb, cy) in self.Gi.nodes():
                                            li2 = [(nstrb, cy)]
                                else:
                                    li2 = [(nstrb,)]

                                # if cy==4:
                                #     printnstr,nstrb
                                #if iprint:
                                #     print("li1",li1)
                                #     print("li2",li2)
                                #if cy == 91:
                                #    print("     ",li2)

                                for i1 in li1:
                                    for i2 in li2:
                                        if (i1[0] != i2[0]):
                                            if ((len(i1) == 2) & (len(i2) == 2)):
                                                # print"RR"
                                                self.Gi.add_edge(i1, i2)
                                                self.Gi.add_edge(i2, i1)
                                            if ((len(i1) == 2) & (len(i2) == 3)):
                                                # print"RT"
                                                if i1[1] == i2[1]:
                                                    self.Gi.add_edge(i1, i2)
                                            if ((len(i1) == 3) & (len(i2) == 2)):
                                                # print"TR"
                                                if i1[2] == i2[1]:
                                                    self.Gi.add_edge(i1, i2)
                                            if ((len(i1) == 3) & (len(i2) == 3)):
                                                # print"TT"
                                                if i1[2] == i2[1]:
                                                    self.Gi.add_edge(i1, i2)
                                                if i2[2] == i1[1]:
                                                    self.Gi.add_edge(i2, i1)
                                            if ((len(i1) == 1) & (len(i2) == 3)):
                                                # print"DT"
                                                if i2[1] == cy:
                                                    self.Gi.add_edge(i1, i2)
                                            if ((len(i1) == 3) & (len(i2) == 1)):
                                                # print"TD"
                                                if i1[2] == cy:
                                                    self.Gi.add_edge(i1, i2)
                                            if ((len(i1) == 1) & (len(i2) == 2)):
                                                # print"DR"
                                                self.Gi.add_edge(i1, i2)
                                            if ((len(i1) == 2) & (len(i2) == 1)):
                                                # print"RD"
                                                self.Gi.add_edge(i1, i2)
                                            if ((len(i1) == 1) & (len(i2) == 1)):
                                                # print"DD"
                                                self.Gi.add_edge(i1, i2)
        if verbose :
            Gipbar.update(66.)
        # updating the list of interactions of a given cycle
        pbartmp = pbar(verbose,total=100.,
                       desc ='update interraction list',
                       leave=False,
                       position=tqdmpos+1)
        for c in self.Gt.node:
            if verbose:
                pbartmp.update(cpt)
            if c != 0:
                vnodes = self.Gt.nodes[c]['polyg'].vnodes
                for k in npt:
                    self.Gt.nodes[c]['inter'] += [(k,)]

        if verbose :
            Gipbar.update(100.)

        # cleaning deadend Gi 
        # if outdoor for all nodes of Gi 
        #   if not diffraction 
        #       if termination cycle is indoor 
        #           or if starting point is indoor 
        # then delete interaction 
        ldelete = []

        if self.typ=='outdoor':
            for k in list(dict(self.Gi.node).keys()):
                # R and T 
                if len(k)>1:
                    segtype = self.Gs.nodes[k[0]]['name']
                    if ((segtype!='AIR') and (segtype!='_AIR')):
                        cyend = k[-1] 
                        if self.Gt.nodes[cyend]['indoor']:
                            # if k[0]>0:
                            #     if self.Gs.nodes[k[0]]['name']!='AIR':
                            ldelete.append(k)
                        if len(k) == 3:
                            cystart = k[1]
                            if self.Gt.nodes[cystart]['indoor']:
                                # if k[0]>0:
                                #     if self.Gs.nodes[k[0]]['name']!='AIR':
                                ldelete.append(k)       

        self.Gi.remove_nodes_from(ldelete)
        # build adjacency matrix of Gi graph
        self.Gi_A = nx.adjacency_matrix(self.Gi)
        #store list of nodes of Gi ( for keeping order)
        self.Gi_no = self.Gi.nodes()

    def filterGi(self, situ='outdoor'):
        """ filter Gi to manage indoor/outdoor situations

        Not called

        """

        # get outdoor notes
        cy = np.array(self.Gt.nodes())
        uout = np.where([not self.Gt.nodes[i]['indoor'] for i in cy])
        cyout = cy[uout]

        inter = self.Gi.nodes()
        Ti = [i for i in inter if ((len(i) == 3) and i[0] > 0)]
        Ri = [i for i in inter if ((len(i) == 2) and i[0] > 0)]
        Di = [i for i in inter if i[0] < 0]

        Ti = [i for i in Ti if ((i[1] in cyout) and (i[2] in cyout))]
        Ri = [i for i in Ri if (i[1] in cyout)]
        Di = [i for i in Di if (i in self.ldiffout)]

        rinter = Ti + Ri + Di

        rGi = nx.subgraph(self.Gi, rinter)
        rGi.pos = {i: self.Gi.pos[i] for i in self.Gi.nodes()}

        self.Gi = rGi
        self.Gi.pos = rGi.pos

    #@profile
    def outputGi(self, verbose=False, tqdmpos=0.):
        """ filter output of Gi edges

        Parameters
        ----------

        L : Layout

        Notes
        -----

        Let assume a sequence (nstr0,nstr1,{nstr2A,nstr2B,...}) in a signature.
        This function checks whether this sequence is feasible or not
        , whatever the type of nstr0 and nstr1.
        The feasible outputs from nstr0 to nstr1 are stored in an output field of
        edge (nstr0,nstr1)

        See Also
        --------

        pylayers.util.cone.Cone.from2seg
        pylayers.util.cone.Cone.belong_seg


        """


        assert('Gi' in self.__dict__)

        oGipbar=pbar(verbose,total=100.,leave=False,desc='OutputGi',position=tqdmpos)
        # loop over all edges of Gi
        Nedges = len(self.Gi.edges())
        cpt = 100./Nedges
        # print "Gi Nedges :",Nedges
        cn = cone.Cone()

        for k, e in enumerate(self.Gi.edges()):

            # if (k%100)==0:
            # print"edge :  ",k
            # extract  both termination interactions nodes
            if verbose:
                oGipbar.update(cpt)

            i0 = e[0]
            i1 = e[1]

            nstr0 = i0[0]
            nstr1 = i1[0]

            # list of authorized outputs. Initialized void
            output = []

            # nstr1 : segment number of central interaction
            if nstr1 > 0:
                # central interaction is a segment
                pseg1 = self.seg2pts(nstr1).reshape(2, 2).T
                #pseg1 = self.s2pc[nstr1].toarray().reshape(2, 2).T
                # list all potential successors of interaction i1
                # v1.1 i2 = nx.neighbors(self.Gi, i1)
                i2 = list(dict(self.Gi[i1]).keys())
                # create a Cone object
                #cn = cone.Cone()
                # if starting from segment
                if nstr0 > 0:
                    pseg0 = self.seg2pts(nstr0).reshape(2, 2).T
                    #pseg0 = self.s2pc[nstr0].toarray().reshape(2, 2).T
                    # if nstr0 and nstr1 are connected segments
                    # v1.1 if (len(np.intersect1d(nx.neighbors(self.Gs, nstr0), nx.neighbors(self.Gs, nstr1))) == 0):
                    if (len(np.intersect1d(self.Gs[nstr0], self.Gs[nstr1])) == 0):
                        # from 2 not connected segment
                        cn.from2segs(pseg0, pseg1)
                    else:
                        # from 2 connected segments
                        cn.from2csegs(pseg0, pseg1)
                # if starting from a point
                else:
                    pt = np.array(self.Gs.pos[nstr0])
                    cn.fromptseg(pt, pseg1)
                    #

                ipoints = [x for x in i2 if len(x)==1 ]
                # i0      i1     i2[x]
                # Avoid to have the same diffaction point after reflection 
                # exemple :  (-10,),(245,12),(-10,) impossible 
                #                                 nstr0  nstr1 
                if nstr0<0:
                    ipoints = [x for x in ipoints if x[0]!=nstr0] 
                #ipoints = filter(lambda x: len(x) == 1, i2)
                pipoints = np.array([self.Gs.pos[ip[0]] for ip in ipoints]).T
                # filter tuple (R | T)
                #istup = filter(lambda x : type(eval(x))==tuple,i2)
                # map first argument segment number
                #isegments = np.unique(map(lambda x : eval(x)[0],istup))
                #isegments = np.unique(
                #    filter(lambda y: y > 0, map(lambda x: x[0], i2)))
                isegments = np.unique(np.array([ s for s in [ n[0] for n in i2]
                                                if s >0 ] ))
                # if nstr0 and nstr1 are adjescent segment remove nstr0 from
                # potential next interaction
                # Fix 01/2017
                # This is not always True if the angle between
                # the two adjascent segments is < pi/2
                # v1.1 nb_nstr0 = self.Gs.neighbors(nstr0)
                # v1.1 nb_nstr1 = self.Gs.neighbors(nstr1)
                nb_nstr0 = self.Gs[nstr0]
                nb_nstr1 = self.Gs[nstr1]

                #common_point = np.intersect1d(nb_nstr0,nb_nstr1)
                common_point = np.array([x for x in nb_nstr0 if x in nb_nstr1]).astype('int')

                if len(common_point) == 1:
                    num0 = [x for x in nb_nstr0 if x != common_point]
                    num1 = [x for x in nb_nstr1 if x != common_point]
                    p0 = np.array(self.Gs.pos[num0[0]])
                    p1 = np.array(self.Gs.pos[num1[0]])
                    pc = np.array(self.Gs.pos[common_point[0]])
                    v0 = p0 - pc
                    v1 = p1 - pc
                    v0n = v0/np.sqrt(np.sum(v0*v0))
                    v1n = v1/np.sqrt(np.sum(v1*v1))
                    if np.dot(v0n,v1n)<=0:
                        isegments = np.array([ x for x in isegments if x != nstr0 ])
                    #    filter(lambda x: x != nstr0, isegments))
                # there are one or more segments
                if len(isegments) > 0:
                    points = self.seg2pts(isegments)
                    pta = points[0:2, :]
                    phe = points[2:, :]
                    # add difraction points
                    # WARNING Diffraction points are added only if a segment is seen
                    # it should be the case in 99% of cases

                    if len(ipoints) > 0:
                        isegments = np.hstack((isegments, np.array(ipoints)[:, 0]))
                        pta = np.hstack((pta, pipoints))
                        phe = np.hstack((phe, pipoints))

                    # cn.show()

                    # if i0 == (38,79) and i1 == (135,79,23):
                    #     printi0,i1
                    #     import ipdb
                    #     ipdb.set_trace()
                    # i1 : interaction T
                    if len(i1) == 3:
                        #if ((e[0]==(53,17)) and (e[1]==(108,17,18))):
                        #    typ, prob = cn.belong_seg(pta, phe,visu=True)
                        #else:
                        #typ, prob = cn.belong_seg_old(pta, phe, prob=False)
                        typ = cn.belong_seg(pta, phe)
                        # if bs.any():
                        #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                        # if ~bs.any():
                        #    plu.displot(pta[:,~bs],phe[:,~bs],color='k')

                    # i1 : interaction R --> mirror
                    if len(i1) == 2:
                        Mpta = geu.mirror(pta, pseg1[:, 0], pseg1[:, 1])
                        Mphe = geu.mirror(phe, pseg1[:, 0], pseg1[:, 1])
                        typ = cn.belong_seg(Mpta, Mphe)
                        #typ, prob = cn.belong_seg_old(Mpta, Mphe, prob=False)
                    # keep segment with typ <> 0
                    utypseg = typ != 0
                    isegkeep = isegments[utypseg]
                    # dict   {numint : proba}
                    #dsegprob = {k: v for k, v in zip(isegkeep, prob[utypseg])}
                    #########
                    # output = filter(lambda x: x[0] in isegkeep, i2)
                    output = [x for x in i2 if x[0] in isegkeep]
                    # probint = map(lambda x: dsegprob[x[0]], output)
                    #probint = [dsegprob[x[0]] for x in output]
                    # dict interaction : proba
                    #dintprob = {k: v for k, v in zip(output, probint)}

                    # keep all segment above nstr1 and in Cone if T
                    # keep all segment below nstr1 and in Cone if R

            else:
                # central interaction is a point (nstr1 <0) 

                # 1) Simple approach
                #       output interaction are all visible interactions
                # 2) TO BE DONE
                #
                #       output of the diffraction points
                #       exploring
                # b
                #          + right of ISB
                #          + right of RSB
                #
                #  + using the wedge cone
                #  + using the incident cone
                #

                # v1.1 output = nx.neighbors(self.Gi, (nstr1,))
                output = self.Gi[(nstr1,)]
                #nout = len(output)
                #probint = np.ones(nout)  # temporarybns
                #dintprob = {k: v for k, v in zip(output, probint)}

            #self.Gi.add_edge(i0, i1, output=dintprob)
            self.Gi.add_edge(i0, i1, output=output)


    def outputGi_new(self,verbose=False,tqdmpos=0.):
        """ filter output of Gi edges

        this version of outputGi, uses sparses matrix instead of NetworkX for MP 
        purpose

        Parameters
        ----------

        L : Layout

        Notes
        -----

        Let assume a sequence (nstr0,nstr1,{nstr2A,nstr2B,...}) in a signature.
        This function checks whether this sequence is feasible or not
        , whatever the type of nstr0 and nstr1.
        The feasible outputs from nstr0 to nstr1 are stored in an output field of
        edge (nstr0,nstr1)

        See Also
        --------

        pylayers.util.cone.Cone.from2seg
        pylayers.util.cone.Cone.belong_seg


        """


        def Gspos(n):
            if n>0:
                return np.mean(self.s2pc[n].toarray().reshape(2,2),axis=0)
            else:
                return self.p2pc[-n].toarray()

        #s2pc = self.s2pc.toarray()
        #s2pu = self.s2pu.toarray()
        #p2pc = self.p2pc.toarray()
        #A = self.Gi_A.toarray()

        assert('Gi' in self.__dict__)

        oGipbar = pbar(verbose,total=100.,leave=False,desc='OutputGi',position=tqdmpos)
        # loop over all edges of Gi
        Nedges = len(self.Gi.edges())
        cpt = 100./Nedges
        # print "Gi Nedges :",Nedges
        for k, e in enumerate(self.Gi.edges()):

            # if (k%100)==0:
            # print"edge :  ",k
            # extract  both termination interactions nodes
            if verbose:
                oGipbar.update(cpt)
            i0 = e[0]  # first interaction 
            i1 = e[1]  # central interaction
            nstr0 = i0[0]
            nstr1 = i1[0]

            # list of authorized outputs. Initialized void
            output = []

            # nstr1 : segment number of central interaction
            if nstr1 > 0:
                # central interaction is a segment
                # pseg1 = self.s2pc[nstr1,:].toarray().reshape(2, 2).T
                pseg1 = self.s2pc[nstr1,:].toarray().reshape(2, 2).T
                # pseg1 = self.s2pc[nstr1,:].data.reshape(2, 2).T
                # pseg1o = self.seg2pts(nstr1).reshape(2, 2).T

                # create a Cone object
                cn = cone.Cone()
                # if starting from segment
                if nstr0 > 0:
                    # pseg0 = self.s2pc[nstr0,:].toarray().reshape(2, 2).T
                    pseg0 = self.s2pc[nstr0,:].toarray().reshape(2, 2).T
                    # pseg0 = self.s2pc[nstr0,:].data.reshape(2, 2).T
                    # pseg0o = self.seg2pts(nstr0).reshape(2, 2).T

                    # if nstr0 and nstr1 are connected segments
                    if self.sgsg[nstr0,nstr1] == 0:
                        # from 2 not connected segment
                        cn.from2segs(pseg0, pseg1)
                    else:
                        # from 2 connected segments
                        cn.from2csegs(pseg0, pseg1)
                # if starting from a point
                else:
                    pt = Gspos(nstr0)[0,:]
                    # pt = np.array(self.Gs.pos[nstr0])
                    cn.fromptseg(pt, pseg1)

                # list all potential successors of interaction i1
                ui2 = self.Gi_no.index(i1)
                ui = np.where(self.Gi_A[ui2,:].toarray()!=0)[1]
                i2 = [self.Gi_no[u] for u in ui]
                # i2 = nx.neighbors(self.Gi, i1)

                # how to find neighbors without network
                # ngi=L.Gi.nodes()
                # A=nx.adjacency_matrix(L.Gi)
                # inter = ngi[10]
                # u = ngi.index(inter)
                # ui = A[u,:].indices
                # neigh_inter = np.array([ngi[u] for u in ui])


                ipoints = [x for x in i2 if len(x)==1 ]

                #ipoints = filter(lambda x: len(x) == 1, i2)
                # pipoints = np.array([self.Gs.pos[ip[0]] for ip in ipoints]).T
                pipoints = np.array([Gspos(ip[0]) for ip in ipoints]).T
                # filter tuple (R | T)
                #istup = filter(lambda x : type(eval(x))==tuple,i2)
                # map first argument segment number
                #isegments = np.unique(map(lambda x : eval(x)[0],istup))
                # isegments = np.unique(
                #     filter(lambda y: y > 0, map(lambda x: x[0], i2)))
                isegments = np.unique([x[0] for x in i2 if x[0]>0])

                # if nstr0 and nstr1 are adjescent segment remove nstr0 from
                # potential next interaction
                # Fix 01/2017
                # This is not always True if the angle between 
                # the two adjascent segments is < pi/2
                # nb_nstr0 = self.Gs.neighbors(nstr0)
                # nb_nstr1 = self.Gs.neighbors(nstr1)
                # nb_nstr0 = np.array([self.s2pu[nstr0,0],self.s2pu[nstr0,1]])
                # nb_nstr1 = np.array([self.s2pu[nstr1,0],self.s2pu[nstr1,1]])
                # nb_nstr0 = self.s2pu[nstr0,:].toarray()[0]
                # nb_nstr1 = self.s2pu[nstr1,:].toarray()[0]

                # first interaction is a point
                if nstr0<0:
                    nb_nstr0 = [nstr0]
                else:
                    nb_nstr0 = self.s2pu[nstr0,:].toarray()[0,:]
                nb_nstr1 = self.s2pu[nstr1,:].toarray()[0,:]
                # common_point = np.intersect1d(nb_nstr0,nb_nstr1)
                common_point = np.array([x for x in nb_nstr0 if x in nb_nstr1])
                #print(common_point)

                # if len(common_point) == 1:
                #     pdb.set_trace()
                if common_point.any():
                    num0 = [x for x in nb_nstr0 if x != common_point]
                    num1 = [x for x in nb_nstr1 if x != common_point]
                    p0 = Gspos(num0[0])[0,:]
                    p1 = Gspos(num1[0])[0,:]
                    pc = Gspos(common_point[0])[0,:]

                    v0 = p0 - pc
                    v1 = p1 - pc

                    v0n = v0/np.sqrt(np.sum(v0*v0))
                    v1n = v1/np.sqrt(np.sum(v1*v1))
                    if np.dot(v0n,v1n)<=0:
                        isegments = np.array([ x for x in isegments if x != nstr0 ]) 
                    #    [ x for x in rle if  x != nstr0, isegments))
                # there are one or more segments
                # if len(isegments) > 0:
                if isegments.any():

                    li1 = len(i1)

                    # points = self.s2pc[isegments,:].toarray().T
                    points = self.s2pc[isegments,:].toarray().T
                    # points = self.s2pc[isegments,:].data.reshape(4,len(isegments))
                    # pointso = self.seg2pts(isegments)

                    pta = points[0:2, :]
                    phe = points[2:, :]
                    # add difraction points
                    # WARNING Diffraction points are added only if a segment is seen
                    # it should be the case in 99% of cases

                    if len(ipoints) > 0:
                        isegments = np.hstack(
                            (isegments, np.array(ipoints)[:, 0]))
                        pta = np.hstack((pta, pipoints))
                        phe = np.hstack((phe, pipoints))

                    # cn.show()

                    # if i0 == (38,79) and i1 == (135,79,23):
                    #     printi0,i1
                    #     import ipdb
                    #     ipdb.set_trace()
                    # i1 : interaction T
                    if li1 == 3:
                        typ, prob = cn.belong_seg(pta, phe, prob=False)
                        # if bs.any():
                        #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                        # if ~bs.any():
                        #    plu.displot(pta[:,~bs],phe[:,~bs],color='k')

                    # i1 : interaction R --> mirror
                    elif li1 == 2:
                        Mpta = geu.mirror(pta, pseg1[:, 0], pseg1[:, 1])
                        Mphe = geu.mirror(phe, pseg1[:, 0], pseg1[:, 1])
                        #typ, prob = cn.belong_seg_old(Mpta, Mphe, prob=False)
                        typ = cn.belong_seg(Mpta, Mphe)
                        # printi0,i1
                        # if ((i0 == (6, 0)) & (i1 == (7, 0))):
                        #    pdb.set_trace()
                        # if bs.any():
                        #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                        # if ~bs.any():
                        #    plu.displot(pta[:,~bs],phe[:,~bs],color='m')
                        #    plt.show()
                        #    pdb.set_trace())
                    ########
                    # SOMETIMES PROBA IS 0 WHEREAS SEG IS SEEN
                    ###########
                    # # keep segment with prob above a threshold
                    # isegkeep = isegments[prob>0]
                    # # dict   {numint : proba}
                    # dsegprob = {k:v for k,v in zip(isegkeep,prob[prob>0])}
                    # 4 lines are replaced by
                    # keep segment with prob above a threshold
                    utypseg = typ != 0
                    isegkeep = isegments[utypseg]
                    # dict   {numint : proba}
                    dsegprob = {k: v for k, v in zip(isegkeep, prob[utypseg])}
                    #########
                    # output = [ x for x in rle if  x[0] in isegkeep, i2)
                    output = [x for x in i2 if x[0] in isegkeep]
                    # probint = map(lambda x: dsegprob[x[0]], output)
                    probint = [dsegprob[x[0]] for x in output]
                    # dict interaction : proba
                    dintprob = {k: v for k, v in zip(output, probint)}

                    # keep all segment above nstr1 and in Cone if T
                    # keep all segment below nstr1 and in Cone if R

            else:
                # central interaction is a point

                # 1) Simple approach
                #       output interaction are all visible interactions
                # 2) TO BE DONE
                #
                #       output of the diffraction points
                #       exploring
                # b
                #          + right of ISB
                #          + right of RSB
                #
                #  + using the wedge cone
                #  + using the incident cone
                #

                # output = nx.neighbors(self.Gi, (nstr1,))
                uout = self.Gi_no.index((nstr1,))
                ui = np.where(self.Gi_A[uout,:].toarray()!=0)[1]
                output = [self.Gi_no[u] for u in ui]

                nout = len(output)
                probint = np.ones(nout)  # temporarybns
                dintprob = {k: v for k, v in zip(output,probint)}

            try:
                self.Gi.add_edge(i0, i1, output=dintprob)
            except:
                pass


    def outputGi_mp(self):
        """ filter output of Gi edges

        Parameters
        ----------

        L : Layout

        Notes
        -----

        Let assume a sequence (nstr0,nstr1,{nstr2A,nstr2B,...}) in a signature.
        This function checks whether this sequence is feasible or not
        , whatever the type of nstr0 and nstr1.
        The feasible outputs from nstr0 to nstr1 are stored in an output field of
        edge (nstr0,nstr1)

        See Also
        --------

        pylayers.util.cone.Cone.from2seg
        pylayers.util.cone.Cone.belong_seg


        """


        # assert('Gi' in self.__dict__)

        # oGipbar=pbar(verbose,total=100.,leave=False,desc='OutputGi',position=tqdmpos)
        # # loop over all edges of Gi
        # Nedges = len(self.Gi.edges())
        # cpt = 100./Nedges
        # print "Gi Nedges :",Nedges
        e = self.Gi.edges()
        #Gi_no = [self.Gi_no]*len(e)

        # densify sparse matrix
        #aGi_A = self.Gi_A.toarray()
        #ap2pc = self.p2pc.toarray()
        #asgsg = self.sgsg.toarray()
        #as2pc = self.s2pc.toarray()
        #as2pu = self.s2pu.toarray()

        global Gi_A
        global Gi_no
        global p2pc
        global sgsg
        global s2pc
        global s2pu

        Gi_A = self.Gi_A
        Gi_no = self.Gi_no
        p2pc = self.p2pc
        sgsg = self.sgsg
        s2pc = self.s2pc
        s2pu = self.s2pu


        #Gi_A = [aGi_A]*len(e)
        #p2pc = [ap2pc]*len(e)
        #s2pc = [as2pc]*len(e)
        #s2pu = [as2pu]*len(e)
        #sgsg = [asgsg]*len(e)

        pool = Pool(cpu_count())

        # multiprocessing style
        #Z=zip(e, Gi_no, Gi_A, p2pc, sgsg, s2pc, s2pu)
        #res = pool.map(outputGi_func,Z)
        Z = zip(e)
        res = pool.map(outputGi_func,Z)
        self.Gi.add_edges_from(res)




        # res = pool.map(outputGi_func_test,e)
        # print('e')
        # time.sleep(1)
        # res = pool.map(outputGi_func_test,Gi_no)
        # print('no')
        # time.sleep(1)
        # res = pool.map(outputGi_func_test,Gi_A)
        # print('A')
        # time.sleep(1)
        # res = pool.map(outputGi_func_test,Gspos)
        # print('pos')
        # time.sleep(1)
        # res = pool.map(outputGi_func_test,sgsg)
        # print('sgsg')
        # time.sleep(1)
        # res = pool.map(outputGi_func_test,s2pc)
        # print('s2pc')
        # time.sleep(1)
        # res = pool.map(outputGi_func_test,s2pu)
        # print('s2pu')
        # time.sleep(1)
        # res = pool.map(outputGi_func_test,Z)
        # print('Z')
        



    #def outputGi_func(arg):
           
        # if (k%100)==0:
        # print"edge :  ",k
        # extract  both termination interactions nodes

        #for k in arg:
        #    Z=arg*arg
        # e=arg[0]
        # s2pc=arg[1]
        # Gs=arg[2]
        # Gi=arg[3]

        # i0 = e[0]
        # i1 = e[1]
        # nstr0 = i0[0]
        # nstr1 = i1[0]
        # print(i0,i1)

        # for k in range(1000):
        #     y=k*k
        # # list of authorized outputs. Initialized void
        # output = []
        # # nstr1 : segment number of central interaction
        # if nstr1 > 0:
        #     # central interaction is a segment
        #     pseg1 = np.array(s2pc[nstr1,:].todense()).reshape(2, 2).T
        #     # create a Cone object
        #     cn = cone.Cone()
        #     # if starting from segment
        #     if nstr0 > 0:
        #         pseg0 = np.array(s2pc[nstr0,:].todense()).reshape(2, 2).T
        #         # if nstr0 and nstr1 are connected segments
        #         if (len(np.intersect1d(nx.neighbors(Gs, nstr0), nx.neighbors(Gs, nstr1))) == 0):
        #             # from 2 not connected segment
        #             cn.from2segs(pseg0, pseg1)
        #         else:
        #             # from 2 connected segments
        #             cn.from2csegs(pseg0, pseg1)
        #     # if starting from a point
        #     else:
        #         pt = np.array(Gs.pos[nstr0])
        #         cn.fromptseg(pt, pseg1)

        #     # list all potential successors of interaction i1
        #     i2 = nx.neighbors(Gi, i1)
        #     ipoints = [x for x in i2 if len(x)==1 ]
        #     #ipoints = [ x for x in rle if  len(x) == 1, i2)
        #     pipoints = np.array([Gs.pos[ip[0]] for ip in ipoints]).T
        #     # filter tuple (R | T)
        #     #istup = filter(lambda x : type(eval(x))==tuple,i2)
        #     # map first argument segment number
        #     #isegments = np.unique(map(lambda x : eval(x)[0],istup))
        #     isegments = np.unique(
        #         filter(lambda y: y > 0, map(lambda x: x[0], i2)))
        #     # if nstr0 and nstr1 are adjescent segment remove nstr0 from
        #     # potential next interaction
        #     # Fix 01/2017
        #     # This is not always True if the angle between 
        #     # the two adjascent segments is < pi/2
        #     nb_nstr0 = Gs.neighbors(nstr0)
        #     nb_nstr1 = Gs.neighbors(nstr1)
        #     common_point = np.intersect1d(nb_nstr0,nb_nstr1)
        #     if len(common_point) == 1:
        #         num0 = [x for x in nb_nstr0 if x != common_point]
        #         num1 = [x for x in nb_nstr1 if x != common_point]
        #         p0 = np.array(Gs.pos[num0[0]])
        #         p1 = np.array(Gs.pos[num1[0]])
        #         pc = np.array(Gs.pos[common_point[0]])
        #         v0 = p0-pc 
        #         v1 = p1-pc 
        #         v0n = v0/np.sqrt(np.sum(v0*v0))
        #         v1n = v1/np.sqrt(np.sum(v1*v1))
        #         if np.dot(v0n,v1n)<=0:
        #             isegments = np.array([ x for x in isegments if x != nstr0 ]) 
        #         #    [ x for x in rle if  x != nstr0, isegments))
        #     # there are one or more segments
        #     if len(isegments) > 0:
        #         points = np.array(s2pc[isegments,:].todense()).T
        #         pta = points[0:2, :]
        #         phe = points[2:, :]
        #         # add difraction points
        #         # WARNING Diffraction points are added only if a segment is seen
        #         # it should be the case in 99% of cases

        #         if len(ipoints) > 0:
        #             isegments = np.hstack(
        #                 (isegments, np.array(ipoints)[:, 0]))
        #             pta = np.hstack((pta, pipoints))
        #             phe = np.hstack((phe, pipoints))

        #         # cn.show()

        #         # if i0 == (38,79) and i1 == (135,79,23):
        #         #     printi0,i1
        #         #     import ipdb
        #         #     ipdb.set_trace()
        #         # i1 : interaction T
        #         if len(i1) == 3:
        #             typ, prob = cn.belong_seg(pta, phe)
        #             # if bs.any():
        #             #    plu.displot(pta[:,bs],phe[:,bs],color='g')
        #             # if ~bs.any():
        #             #    plu.displot(pta[:,~bs],phe[:,~bs],color='k')

        #         # i1 : interaction R --> mirror
        #         if len(i1) == 2:
        #             Mpta = geu.mirror(pta, pseg1[:, 0], pseg1[:, 1])
        #             Mphe = geu.mirror(phe, pseg1[:, 0], pseg1[:, 1])
        #             typ, prob = cn.belong_seg(Mpta, Mphe)
        #             # printi0,i1
        #             # if ((i0 == (6, 0)) & (i1 == (7, 0))):
        #             #    pdb.set_trace()
        #             # if bs.any():
        #             #    plu.displot(pta[:,bs],phe[:,bs],color='g')
        #             # if ~bs.any():
        #             #    plu.displot(pta[:,~bs],phe[:,~bs],color='m')
        #             #    plt.show()
        #             #    pdb.set_trace())
        #         ########
        #         # SOMETIMES PROBA IS 0 WHEREAS SEG IS SEEN
        #         ###########
        #         # # keep segment with prob above a threshold
        #         # isegkeep = isegments[prob>0]
        #         # # dict   {numint : proba}
        #         # dsegprob = {k:v for k,v in zip(isegkeep,prob[prob>0])}
        #         # 4 lines are replaced by
        #         # keep segment with prob above a threshold
        #         utypseg = typ != 0
        #         isegkeep = isegments[utypseg]
        #         # dict   {numint : proba}
        #         dsegprob = {k: v for k, v in zip(isegkeep, prob[utypseg])}
        #         #########
        #         # output = [ x for x in rle if  x[0] in isegkeep, i2)
        #         output = [x for x in i2 if x[0] in isegkeep]
        #         # probint = map(lambda x: dsegprob[x[0]], output)
        #         probint = [dsegprob[x[0]] for x in output]
        #         # dict interaction : proba
        #         dintprob = {k: v for k, v in zip(output, probint)}

        #         # keep all segment above nstr1 and in Cone if T
        #         # keep all segment below nstr1 and in Cone if R

        # else:
        #     # central interaction is a point

        #     # 1) Simple approach
        #     #       output interaction are all visible interactions
        #     # 2) TO BE DONE
        #     #
        #     #       output of the diffraction points
        #     #       exploring
        #     # b
        #     #          + right of ISB
        #     #          + right of RSB
        #     #
        #     #  + using the wedge cone
        #     #  + using the incident cone
        #     #

        #     output = nx.neighbors(Gi, (nstr1,))
        #     nout = len(output)
        #     probint = np.ones(nout)  # temporarybns
        #     dintprob = {k: v for k, v in zip(output, probint)}

        # return(i0,i1,dintprob)
        #self.Gi.add_edge(i0, i1, output=dintprob)

        

    def intercy(self, ncy, typ='source'):
        """ return the list of interactions seen from a cycle

        Parameters
        ----------

        ncy : cycle number( Project -> save project)
        typ : string
            if 'source' connect source cycle
            if 'target' connect target cycle

        Notes
        -----

        This method is called at the beginning of signature evaluation in order 
        to get the starting and ending interaction. It exploits the information 
        contained in teh graph Gi.

        """

        # list of interactions
        lint = self.Gi.nodes

        # list of tuple interactions (R|T)
        lD = [x for x in lint if len(x)==1]
        lR = [x for x in lint if len(x)==2]
        lT = [x for x in lint if len(x)==3]
        # lD = [ x for x in rle if  len(x) == 1, lint)
        # lR = [ x for x in rle if  len(x) == 2, lint)
        # lT = [ x for x in rle if  len(x) == 3, lint)

        # visible R|T source cycle is ncy

        lR = [ x for x in lR if  x[1] == ncy ]
        if typ == 'source':
            lT = [ x for x in lT if  x[1] == ncy ]
        if typ == 'target':
            lT = [ x for x in lT if  x[2] == ncy ]
        if typ == 'all':
            lT = lT
        # Finding the diffraction points
        # Diffraction points are different from indoor cycle and outdoor
        # cycles
        #
        # TODO check wedge validity.
        #

        vnodes = self.Gt.nodes[ncy]['polyg'].vnodes
        vpoints = [ x for x in vnodes if  x < 0 ]
        lD = []
        for x in vpoints:
            if x in self.ddiff:
                for y in self.ddiff[x][0]:
                    if y == ncy:
                        lD.append((x,))
        # indoor = self.Gt.nodes[ncy]['indoor']
        # if indoor:
        #     lD = map(lambda y : (y,),filter(lambda x : x in
        #                                     self.ldiffin,vpoints))
        # else:
        #     lD = map(lambda y : (y,),filter(lambda x : x in
        #                                     self.ldiffout,vpoints))

        return lR, lT, lD

    def show(self, **kwargs):
        """ show layout

        See also
        --------

        showG

        """
        defaults = {'show': True,
                    'fig': [],
                    'ax': [],
                    'nodes': False,
                    'edges': True,
                    'labels': False,
                    'alphan': 1.0,
                    'alphae': 1.0,
                    'width': 2,
                    'node_color': 'w',
                    'edge_color': 'k',
                    'node_size': 200,
                    'font_size': 30,
                    'nodelist': [],
                    'figsize': (5, 5),
                    'mode': 'cycle',
                    }
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        lair = []
        if 'AIR' in self.name:
            lair = self.name['AIR']
        if '_AIR' in self.name:
            lair = lair + self.name['_AIR']

        #
        # tsg : list of segment index for mapping with self.tahe

        segfilt = [ x for x in self.tsg if  x not in lair ]
        # get the association between segment and nx edges
        edges = self.Gs.edges()
        Ne = len(edges)

        # segments = np.array(edges)[:,0]
        # segments are >0 index so max in necesssarily
        # a segment number whatever the order
        segments = np.array([max(x) for x in edges])

        dse = {k: v for k, v in zip(segments, range(Ne))}

        edfilt = list(
            np.ravel(np.array(map(lambda x: [dse[x] - 1, dse[x]], segfilt))))

        # edgelist is to be understood as edges of Graph and not segments of
        # Layout
        if hasattr(self,'extent'):
            fig, ax = gkml.gearth_fig(self.extent,self.extent_c)
        else:
            fig = plt.gca()
            ax = plt.gca()

        fig, ax = self.showG('s', nodes=False, edgelist=edfilt,fig=fig,ax=ax)

        # display degree 1 nodes
        if 1 in self.degree:
            ldeg1 = list(self.degree[1])
            print(ldeg1)
            fig, ax = self.showG('s',
                                 fig=fig,
                                 ax=ax,
                                 nodelist=ldeg1,
                                 edges=kwargs['edges'],
                                 nodes=kwargs['nodes'],
                                 node_size=kwargs['node_size'],
                                 node_color='r')

        # display degree 4 nodes
        if 4 in self.degree:
            ldeg4 = list(self.degree[4])
            fig, ax = self.showG('s',
                                 fig=fig,
                                 ax=ax,
                                 nodelist=ldeg4,
                                 edges=kwargs['edges'],
                                 nodes=kwargs['nodes'],
                                 node_size=kwargs['node_size'],
                                 node_color='g')

        if hasattr(self,'extent'):
            pnglayout = 'kmllayout.png'
            kmzlayout = 'kmzlayout.kmz'
            fig.savefig(pnglayout,transparent=True,format='png')
            gkml.make_kml(self.extent,
                 figs = [pnglayout],
                 kmzfile = kmzlayout,
                 name = 'Layout')

        #     if k==1:
        #         fig,ax = self.showG('s',fig=fig,ax=ax,nodelist=ldeg,edges=False,nodes=True,node_size=50,node_color='c')
        #     if k==4:
        #         fig,ax = self.showG('s',fig=fig,ax=ax,nodelist=ldeg,nodes=False,node_size=50,node_color='b')

    def show_poly(self,**kwargs):
       """ show polygons

       Parameters
       ----------

       bmap : boolean (display smpoy map)

       Returns
       -------

       fig,ax,args,

       """

       import smopy
       bmap = kwargs.pop('bmap',True)
       figsize = kwargs.pop('figsize',(15,15))
       lm = self.extent[0]
       lM = self.extent[1]
       Lm = self.extent[2]
       LM = self.extent[3]
       zoom = 15
       self.map = smopy.Map((Lm,lm,LM,lM), z=zoom)
       bcond = bmap and hasattr(self,'map')
       ax = self.map.show_mpl(figsize=figsize)

       #ax = kwargs.pop('ax',plt.gca())

       BLUE='#6699cc'
       verts = []
       lbdg_height = []
       if hasattr(self.dpoly,'_xy'):
           for kpoly in self.dpoly:
               verts.append(self.dpoly[kpoly]._xy.T)
       else:
           for kpoly in self.dpoly:
               connect = self.dpoly[kpoly]['connect']
               z = self.dpoly[kpoly]['z']
               building_height = z[1] - z[0]
               lbdg_height.append(building_height)
               lpol = []
               lprev = 0
               #
               # Take care of the order of the sequence of points
               #
               for seg in connect:
                   ta = self.Gs.nodes[seg]['connect'][0]
                   he = self.Gs.nodes[seg]['connect'][1]
                   if ((he==lprev) or (lprev==0)):
                       lpol.append((self.Gs.pos[ta][0],self.Gs.pos[ta][1]))
                       lprev = ta
                   else:
                       lpol.append((self.Gs.pos[he][0],self.Gs.pos[he][1]))
                       lprev = he
               llL =  [ self.m(x[0],x[1],inverse=True) for x in lpol ]
               lpol2 = [ self.map.to_pixels(lL[1],lL[0]) for lL in llL]

               verts.append(lpol2)
       hmax = np.max(lbdg_height)
       hmin = np.min(lbdg_height)
       arg =  255*(np.array(lbdg_height)-hmin)/(hmax-hmin)
       facecolors = cm.jet(arg.astype(int))
       coll = PolyCollection(verts, edgecolors='k', facecolors=facecolors)
       ax.add_collection(coll)
       plt.axis('on')
       #ax.autoscale_view()
       #ax.axis('equal')
       fig = plt.gcf()
       return fig,ax,arg,facecolors

    def showG(self, graph='s', **kwargs):
        """ show the different graphs

        Parameters
        ----------

        graph : char
            't' : Gt 'r' : Gr 's' : Gs 'v' : Gv  'i' : Gi
        fig : matplotlib figure
            []
        ax : matplotlib figure
            []
        show : boolean
            False
        nodes : boolean
            alse
        edges : boolean
            True
        airwalls | aw: boolean
            display airwalls (False)
        subseg: boolean
            display subsegments (False)
        slab : boolean
            display color and width of slabs (False)
        labels : boolean |list
            display graph labels (False)
            if list precise label of which cycle to display
            (e.g. ['t'])
        alphan : float
            transparency of nodes (1.0)
        alphae : float
            transparency of edges (1.0)
        width : float
            line width (2)
        node_color: string
            w
        posnode_color: string
            positive node color (k)
        negnode_color: string
            negative node color (b)
        edge_color : string
            k
        node_size : float
            20
        font_size : float
            15,
        nodelist : list
            list of nodes to be displayed (all)
        edgelist : list
            list of edges to be displayed (all)
        mode : string
            'cycle' | 'none' | 'room'
        alphacy : string
            transparency of cycles (0.8)
        colorcy :
            '#abcdef'
        linter : list
            list of interaction for Gi
            ['RR','TT','RT','TR','RD','DR','TD','DT','DD']
        show0 : boolean
            If true display connection to cycle  0 of Gt (False)
        eded : boolean
            True
        ndnd : boolean
            True
        nded : boolean
            True
        width : int
            2
        nodelist : list
            []
        overlay : boolean

        diffraction :boolean
            False


        defaults = {'show': False,
                    'fig': [],
                    'ax': [],
                    'nodes': False,
                    'edges': True,
                    'sllist':[],
                    'airwalls': False,
                    'subseg': False,
                    'slab': True,
                    'labels': False,
                    'alphan': 1.0,
                    'alphae': 1.0,
                    'width': 2,
                    'node_color':'w',
                    'edge_color':'k',
                    'node_size':20,
                    'font_size':15,
                    'nodelist': [],
                    'edgelist': [],
                    'figsize': (5,5),
                    'mode':'nocycle',
                    'alphacy':0.8,
                    'colorcy':'abcdef',
                    'linter' : ['RR','TT','RT','TR','RD','DR','TD','DT','DD'],
                    'show0':False,
                    'axis':False,
                    'overlay':False,
                    'diffraction':False
                    }


        Examples
        --------

        .. plot::
            :include-source:

        >>> from pylayers.gis.layout import  *
        >>> import matplotlib.pyplot as plt
        >>> L = Layout('TA-Office.lay')
        >>> L.dumpr()
        >>> fig = plt.figure(figsize=(10,10))
        >>> ax = fig.add_subplot(221)
        >>> fig,ax = L.showG('s',fig=fig,ax=ax)
        >>> tis = plt.title("Gs")
        >>> ax = fig.add_subplot(222)
        >>> fig,ax = L.showG('t',fig=fig,ax=ax)
        >>> tit = plt.title("Gt")
        >>> ax = fig.add_subplot(223)
        >>> fig,ax = L.showG('r',fig=fig,ax=ax)
        >>> tic = plt.title("Gr")
        >>> ax = fig.add_subplot(224)
        >>> fig,ax = L.showG('v',fig=fig,ax=ax)
        >>> tiv = plt.title("Gv")
        >>> plt.show()

        See Also
        --------

        pylayers.util.graphutil.draw

        """
        defaults = {'show': False,
                    'fig': [],
                    'ax': [],
                    'nodes': [],
                    'edges': True,
                    'sllist': [],
                    'airwalls': False,
                    'aw': [],
                    'subseg': False,
                    'slab': True,
                    'labels': False,
                    'alphan': 1.0,
                    'alphae': 1.0,
                    'width': 2,
                    'node_color': 'w',
                    'edge_color': '',
                    'node_size': 20,
                    'font_size': 15,
                    'nodelist': [],
                    'edgelist': [],
                    'figsize': (8, 8),
                    'mode': 'nocycle',
                    'alphacy': 0.8,
                    'colorcy': '#abcdef',
                    'lvis': ['nn', 'ne', 'ee'],
                    'linter': ['RR', 'TT', 'RT', 'TR', 'RD', 'DR', 'TD', 'DT', 'DD'],
                    'show0': False,
                    'axis': True,
                    'overlay': False,
                    'diffraction': False
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['aw'] != []:
            kwargs['airwalls'] = kwargs['aw']

        if 'graph' in kwargs:
            graph = kwargs['graph']

        # get color dictionnary from pyutil

        cold = pyu.coldict()

        if isinstance(kwargs['labels'], list):
            labels = kwargs['labels']
        elif kwargs['labels'] == True:
            labels = ['s', 't', 'v', 'i', 'w']
        elif isinstance(kwargs['labels'], str):
            labels = kwargs['labels']
        else:
            labels = []

        if isinstance(kwargs['nodes'], list):
            dis_nodes = kwargs['nodes']
        elif kwargs['nodes'] == True:
            dis_nodes = ['s', 't', 'v', 'i', 'w']
        elif isinstance(kwargs['nodes'], str):
            dis_nodes = kwargs['nodes']
        else:
            dis_nodes = []
        #
        # s : structure graph
        #
        if 's' in graph:

            # not efficient
            G = self.Gs

            # lss = [ x for x in self.Gs.nodes if # self.Gs.nodes[x].has_key('ss_name')]
            # lss = [ x for x in lss if  len(self.Gs.nodes[x]['ss_name'])>0 ] 

            # keep track of segments already printed

            nodelistbkup = kwargs['nodelist']
            edgelistbkup = kwargs['edgelist']
            widthbkup = kwargs['width']
            nodecolbkup = kwargs['edge_color']
            # sllist  slab list
            try:
                sllist = [kwargs['sllist'].pop()]
            except:
                sllist = list(dict(self.name).keys())

            #
            # Draw segments slab per slab with proper linewidth and color
            #
            for lmat in sllist:
                lseg = self.name[lmat]
                if lseg != []:
                    lseg2 = [np.where(np.array(self.Gs.edges()) == i)[0] for i in lseg]
                    kwargs['edgelist'] = []
                    for y in lseg2:
                        kwargs['edgelist'] = kwargs['edgelist'] + list(y)
                    #kwargs['edgelist'] = list(reduce(lambda x, y: list(x) + list(y), lseg2))
                    if kwargs['slab']:
                        if self.sl[lmat]['color'][0]=="#":
                            kwargs['edge_color'] = self.sl[lmat]['color']
                        else:
                            kwargs['edge_color'] = cold[self.sl[lmat]['color']]

                        kwargs['width'] = self.sl[lmat]['linewidth']
                    else:
                        kwargs['edge_color'] = 'k'
                        kwargs['width'] = 1

                if 's' in labels:
                    kwargs['labels'] = True
                else:
                    kwargs['labels'] = False

                if 's' in dis_nodes:
                    kwargs['nodes'] = True
                else:
                    kwargs['nodes'] = False

            # plotting graph G
                kwargs['fig'], kwargs['ax'] = gru.draw(G, **kwargs)

            kwargs['nodelist'] = nodelistbkup
            kwargs['width'] = widthbkup
            kwargs['edge_color'] = nodecolbkup
            kwargs['edgelist'] = edgelistbkup

            if kwargs['subseg']:
                #
                # Display doors and windows subsegments with a slight offset
                #
                cold = pyu.coldict()
                d = self.subseg()
                for ss in list(dict(d).keys()):
                    color = cold[self.sl[ss]['color']]
                    for ns in d[ss]:
                        norm = self.Gs.nodes[ns[0]]['norm']
                        # v1.1 np1, np2 = self.Gs.neighbors(ns[0])
                        np1, np2 = self.Gs[ns[0]]
                        x = np.array(
                            [self.Gs.pos[np1][0], self.Gs.pos[np2][0]])
                        y = np.array(
                            [self.Gs.pos[np1][1], self.Gs.pos[np2][1]])
                        xoff = (1 + ns[1]) * 0.05 * norm[0]
                        yoff = (1 + ns[1]) * 0.05 * norm[1]
                        kwargs['ax'].plot(x + xoff, y + yoff,
                                          linewidth=2, color=color)

        #
        # t : graph of cycles
        #
        if 't' in graph:
            G = self.Gt
            if not kwargs['show0']:
                # filter out the 0 cycle
                nodes = list(G.nodes())
                edges = list(G.edges())
                nodf = [ x for x in nodes if x != 0 ]
                edf = [ x for x in np.arange(len(edges)) if ((edges[x][0]!=0) &
                                                             (edges[x][1]!=0))
                       ]
                kwargs['nodelist'] = nodf
                kwargs['edgelist'] = edf
            else:
                kwargs['nodelist'] = G.nodes()
                kwargs['edgelist'] = np.arange(len(G.edges()))

            if kwargs['edge_color'] == '':
                kwargs['edge_color'] = 'r'
            if 't' in labels:
                kwargs['labels'] = True
            else:
                kwargs['labels'] = False
            if 't' in dis_nodes:
                kwargs['nodes'] = True
            else:
                kwargs['nodes'] = False
            fig, ax = gru.draw(G, **kwargs)
            kwargs['fig'] = fig
            kwargs['ax'] = ax
        #
        # r : graph of rooms
        #
        if 'r' in graph:
            G = self.Gr
            if kwargs['edge_color'] == '':
                kwargs['edge_color'] = 'g'

            kwargs['fig'], kwargs['ax'] = gru.draw(self.Gs,
                                                   nodes=False, edges=True, alphacy=1.,
                                                   fig=kwargs['fig'], ax=kwargs['ax'], labels=False)
            if 'r' in labels:
                kwargs['labels'] = True
            else:
                kwargs['labels'] = False
            if 'r' in dis_nodes:
                kwargs['nodes'] = True
            else:
                kwargs['nodes'] = False
            fig, ax = gru.draw(G, **kwargs)
            kwargs['fig'] = fig
            kwargs['ax'] = ax
        #
        # v : visibility graph
        # In blue : segment segment
        # In red  : point point (Diffraction)
        # In green : point segment
        #
        if 'v' in graph:

            G = self.Gv
            G.pos = {}
            # nodes of Gv are nodes of Gs
            G.pos.update(self.Gs.pos)

            if kwargs['edge_color'] == '':
                kwargs['edge_color'] = 'm'

            edges = list(G.edges())
            rle = range(len(edges))
            eded = [ x for x in rle if (edges[x][0] > 0) & (edges[x][1] > 0)]
            ndnd = [ x for x in rle if (edges[x][0] < 0) & (edges[x][1] < 0)]
            nded = [ x for x in rle if (((edges[x][0] < 0) & (edges[x][1] > 0)) |
                                        ((edges[x][0] > 0) & (edges[x][1] < 0)))]
            if 'v' in labels:
                kwargs['labels'] = True
            else:
                kwargs['labels'] = False
            if 'v' in dis_nodes:
                kwargs['nodes'] = True
            else:
                kwargs['nodes'] = False

            if 'ee' in kwargs['lvis']:
                kwargs['edgelist'] = eded
                kwargs['edge_color'] = 'blue'
                kwargs['node_size'] = 200
                kwargs['fig'], kwargs['ax'] = gru.draw(G, **kwargs)
            if 'nn' in kwargs['lvis']:
                kwargs['edgelist'] = ndnd
                kwargs['edge_color'] = 'red'
                kwargs['fig'], kwargs['ax'] = gru.draw(G, **kwargs)
            if 'ne' in kwargs['lvis']:
                kwargs['edgelist'] = nded
                kwargs['edge_color'] = 'green'
                kwargs['fig'], kwargs['ax'] = gru.draw(G, **kwargs)
        #
        # i :  interaction graph
        #
        if 'i' in graph:

            G = self.Gi

            if kwargs['edge_color'] == '':
                kwargs['edge_color'] = 'k'

            #
            # Parsing the type of interactions
            #

            edges = list(G.edges())

            # range len edges

            rle = range(len(edges))

            DD = [ x for x in rle if   ((len(edges[x][0]) == 1) &
                                    (len(edges[x][1]) == 1))]

            RR = [ x for x in rle if  ((len(edges[x][0]) == 2) &
                                   (len(edges[x][1]) == 2))]

            TT = [ x for x in rle if  ((len(edges[x][0]) == 3) &
                                   (len(edges[x][1]) == 3))]

            RT = [ x for x in rle if  ((len(edges[x][0]) == 2) &
                                   (len(edges[x][1]) == 3))]

            TR = [ x for x in rle if  ((len(edges[x][0]) == 3) &
                                   (len(edges[x][1]) == 2))]

            RD = [ x for x in rle if   ((len(edges[x][0]) == 2) &
                                    (len(edges[x][1]) == 1))]

            TD = [ x for x in rle if   ((len(edges[x][0]) == 3) &
                                    (len(edges[x][1]) == 1))]

            DR = [ x for x in rle if   ((len(edges[x][0]) == 1) &
                                    (len(edges[x][1]) == 2))]

            DT = [ x for x in rle if   ((len(edges[x][0]) == 1) &
                                        (len(edges[x][1]) == 3))]

            tabcol = ['b', 'g', 'r', 'm', 'c', 'orange',
                      'purple', 'maroon', 'purple', 'k'][::-1]
            li = []
            if 'i' in labels:
                kwargs['labels'] = True
            else:
                kwargs['labels'] = False
            if 'v' in dis_nodes:
                kwargs['nodes'] = True
            else:
                kwargs['nodes'] = False
            for inter in kwargs['linter']:
                if len(eval(inter)) > 0:
                    li.append(inter)
                    kwargs['edgelist'] = eval(inter)
                    # ndlist = map(lambda x: edges[x][0],kwargs['edgelist'])+\
                    #          map(lambda x: edges[x][1],kwargs['edgelist'])
                    #ndlist = map(lambda x: edges[x][0], kwargs['edgelist']) +\
                    #    map(lambda x: edges[x][1], kwargs['edgelist'])
                    ndlist = [ edges[x][0]  for x in kwargs['edgelist']] + [edges[x][1]  for x in kwargs['edgelist']]
                    # keep only unique interaction
                    unique = []
                    [unique.append(it) for it in ndlist if it not in unique]

                    kwargs['nodelist'] = unique
                    kwargs['edge_color'] = tabcol.pop()
                    kwargs['fig'], kwargs['ax'] = gru.draw(G, **kwargs)
            legtxt = ['Gs'] + li
            # plt.legend(legtxt)
        #
        # w :  waypoint graph
        #
        if 'w' in graph:

            G = self.Gw

            if kwargs['edge_color'] == '':
                kwargs['edge_color'] = 'k'
            kwargs['fig'], kwargs['ax'] = gru.draw(self.Gs,
                                                   nodes=False, edges=True, alphacy=1.,
                                                   fig=kwargs['fig'], ax=kwargs['ax'], labels=False)
            if 'w' in labels:
                kwargs['labels'] = True
            else:
                kwargs['labels'] = False
            fig, ax = gru.draw(G, **kwargs)
            kwargs['fig'] = fig
            kwargs['ax'] = ax

        args = {'fig': kwargs['fig'], 'ax': kwargs['ax'], 'show': False}

        if len(kwargs['edgelist']) == 0:
            if kwargs['mode'] == 'cycle':
                for k, ncy in enumerate(list(dict(self.Gt.node).keys())):
                    if k != 0:
                        fig, ax = self.Gt.nodes[ncy]['polyg'].plot(
                            alpha=kwargs['alphacy'], color=kwargs['colorcy'], **args)
                        args['fig'] = fig
                        args['ax'] = ax
            if kwargs['mode'] == 'room':
                for k, nro in enumerate(list(dict(self.Gr.node.keys()))):
                    if k != 0:
                        fig, ax = self.Gr.nodes[nro]['cycle'].show(**args)
                        args['fig'] = fig
                        args['ax'] = ax

        kwargs['ax'].axis('scaled')
        if not kwargs['axis']:
            kwargs['ax'].axis('off')

        if kwargs['overlay']:
            imok = False
            if self.display['overlay_file'] != '':
                image = Image.open(os.path.join(
                    basename, pro.pstruc['DIRIMAGE'], self.display['overlay_file']))
                imok = True
            if imok:
                if 'v' in self.display['overlay_flip']:
                    print("flip v")
                    image = image.transpose(Image.FLIP_LEFT_RIGHT)
                if 'h' in self.display['overlay_flip']:
                    image = image.transpose(Image.FLIP_TOP_BOTTOM)
                    print("flip h")
                plt.axis()
                kwargs['ax'].imshow(np.array(image), extent=self.display[
                                    'overlay_axis'], alpha=self.display['alpha'], origin='lower')

        if kwargs['diffraction']:
            if len(self.ddiff.keys())>0:
                pt = np.array([self.Gs.pos[x] for x in self.ddiff.keys()])
                pta = np.array([self.Gs.pos[x] for x in self.lnss])
                kwargs['ax'].scatter(pt[:, 0], pt[:, 1], c='r', s=75)
                if len(self.lnss) > 0:
                    kwargs['ax'].scatter(pta[:, 0], pta[:, 1], c='b', s=20)
        if kwargs['show']:
            plt.show()

        return kwargs['fig'], kwargs['ax']

    def _showGv(self, **kwargs):
        """ show graph Gv (visibility)

        Parameters
        ----------
        display
        fig
        ax
        nodes    : boolean
            display nodes
        edges    : boolean
            display edges

        Returns
        -------
        fig : figure instance
        ax  : axes instance

        """
        defaults = {'show': False,
                    'ax': [],
                    'nodes': False,
                    'eded': True,
                    'ndnd': True,
                    'nded': True,
                    'linewidth': 2,
                    }

        for key, value in defaults.items():
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, value)
                kwargs[key] = value

        if kwargs['ax'] == []:
            fig = plt.figure()
            ax = fig.gca()
        else:
            ax = kwargs['ax']

        nodes = np.array(self.Gv.nodes())
        uneg = list(nodes[np.nonzero(nodes < 0)[0]])
        upos = list(nodes[np.nonzero(nodes > 0)[0]])
        if kwargs['nodes']:
            nx.draw_networkx_nodes(self.Gv, self.Gs.pos, nodelist=upos,
                                   node_color='blue', node_size=300, alpha=0.3)
            nx.draw_networkx_nodes(self.Gv, self.Gs.pos, nodelist=uneg,
                                   node_color='red', node_size=300, alpha=0.3)
            nx.draw_networkx_labels(self.Gv, self.Gs.pos)

        ndnd, nded, eded = gru.edgetype(self.Gv)

        if kwargs['eded']:
            nx.draw_networkx_edges(self.Gv, self.Gs.pos,
                                   edgelist=eded, edge_color='blue', width=2)
        if kwargs['ndnd']:
            nx.draw_networkx_edges(self.Gv, self.Gs.pos,
                                   edgelist=ndnd, edge_color='red', width=2)
        if kwargs['nded']:
            nx.draw_networkx_edges(self.Gv, self.Gs.pos,
                                   edgelist=nded, edge_color='green', width=2)

        if kwargs['show']:
            plt.show()

        return ax

    def waypointGw(self, nroom1, nroom2):
        """ get the waypoint between room1 and room2

        Parameters
        ----------

        nroom1
        nroom2

        Examples
        --------

            >>> from pylayers.gis.layout import *
            >>> L = Layout('TA-Office.lay')
            >>> L.build()

        Notes
        -----

        nodes of Gw are no longer room number

        """
        rooms = nx.dijkstra_path(self.Gw, nroom1, nroom2)
        return(rooms, [tuple(self.Gw.pos[i]) for i in rooms])

    def thwall(self, offx, offy):
        """ Create a list of wall tuples (Transit.world format )

        Parameters
        ----------

        offx
        offy

        Returns
        -------

        walls : list of wall tuples  (Transit format)

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('DLR.lay')
        >>> walls = L.thwall(0,0)

        """
        keyn = list(dict(self.Gs.node).keys())
        walls = []
        for nd in keyn:
            if nd > 0:
                #v1.1 nb = self.Gs.neighbors(nd)
                nb = list(dict(self.Gs[nd]).keys())
                pta = self.Gs.pos[nb[0]]
                phe = self.Gs.pos[nb[1]]
                pn = self.Gs.nodes[nd]['norm']
                name = self.Gs.nodes[nd]['name']
                transition = self.Gs.nodes[nd]['transition']
                sl = self.sl[name]
                thick = sum(sl['lthick'])

                p1 = np.array(pta) + \
                    np.array((pn[0], pn[1])) * thick / 2. + \
                    np.array([offx, offy])

                p2 = np.array(phe) + \
                    np.array((pn[0], pn[1])) * thick / 2.  + \
                    np.array([offx, offy])

                p3 = np.array(phe) - \
                    np.array((pn[0], pn[1])) * thick / 2.  + \
                    np.array([offx, offy])

                p4 = np.array(pta) - \
                    np.array((pn[0], pn[1])) * thick / 2.  + \
                    np.array([offx, offy])

                wall = (tuple(p1), tuple(p2), tuple(p3), tuple(p4))
                if not transition and name != 'AIR':
                    walls.append(wall)
        return(walls)

    def ptin(self, pt=np.array((0, 0, 0))):
        """ check if a point is in the Layout

        Parameters
        ----------

        pt : point (ndarray)

        Returns
        -------

        boolean : True if inside

        See Also
        --------

        ispoint

        """

        pt = pt[:2]

        x = np.array((self.ax[:2]))
        y = np.array((self.ax[2:]))

        # being   in [xmin xmax]
        c0 = pt[0] <= x[1] and pt[0] >= x[0]
        # being   in [ymin ymax]
        c1 = pt[1] <= y[1] and pt[1] >= y[0]

        return (c0 & c1)

    def ptGs2cy(self, n=-1):
        """ Gs node to cycle

        Parameters
        ----------
        upt : point (ndarray)

        Returns
        -------
        ncy : cycle number

        Notes
        -----
            If a cycle contains the Gs pointt this function returns the cycle(s) number
        """
        if n > 0:
            return self.Gs.nodes[n]['ncycles']
        else:
            nseg = list(dict(self.Gs[n]).keys())
            cy = []
            for nn in nseg:
                cy.extend(self.ptGs2cy(nn))
            cy = np.unique(cy).tolist()
            return cy

    def isindoor(self,pt=np.array([0,0])):
        """ test if a point is indoor

        Parameters
        ----------
        pt : np.array 1x2
            2d point

        Returns
        -------
        b1 : boolean
            True if indoor

        """
        cy = self.pt2cy(pt)
        b1 = self.Gt.nodes[cy]['indoor']
        return b1

    def pt2cy(self, pt=np.array((0, 0))):
        """ point to cycle

        Parameters
        ----------
        pt : point (ndarray)

        Returns
        -------
        ncy : cycle number

        Notes
        -----
            If a cycle contains point pt this function returns the cycle number
            Not scalable
            TODO : Find a way to accelerate this function

        See Also
        --------

        Layout.cy2pt

        """

        ptsh = sh.Point(pt[0], pt[1])
        cycle_exists = False

        for ncy in list(dict(self.Gt.nodes).keys()):
            if ncy > 0:
                criter1 = self.Gt.nodes[ncy]['polyg'].touches(ptsh)
                criter2 = self.Gt.nodes[ncy]['polyg'].contains(ptsh)
                if (criter1 or criter2):
                    cycle_exists = True
                    return(ncy)
        if not cycle_exists:
            raise NameError(str(pt) + " is not in any cycle")

    def cy2pt(self, cy=0, h=1.2):
        """returns a point into a given cycle

        Parameters
        ----------

        cy : int
            cycle number

        h : float
            point height

        Returns
        -------

        point  : nd.array
            3d point

        See Also
        --------

        Layout.pt2cy

        """

        if cy in self.Gt.nodes():
            pt = np.array((self.Gt.pos[cy]))
            pt = np.hstack((pt, h))
            return(pt)
        else:
            raise NameError("cycle " + str(cy) + " not in self.Gt")

    def pt2ro(self, pt=np.array((0, 0))):
        """ point to room

        Parameters
        ----------
        pt : point (ndarray)

        Returns
        -------
        nr : Room number

        Notes
        -----
            If a room contains point pt this function returns the room number

        """

        ptsh = sh.Point(pt[0], pt[1])
        ptshinroom = False
        for nr in list(dict(self.Gr.node.keys())):
            if self.Gr.nodes[nr]['polyg'].contains(ptsh)\
                    or self.Gr.nodes[nr]['polyg'].touches(ptsh):
                ptshinroom = True
                return(nr)
        if not ptshinroom:
            raise NameError(str(pt) + " is not in any room")

    def seg2ro(self, seg):
        """ return room number of a point

        Parameters
        ----------

        seg : int

        Returns
        -------
        nr : Room number

        Notes
        -----
        If a room contains point pt this function returns the room number

        """

        rooms = []
        for nr in list(dict(self.Gr.nodes.keys())):
            # if seg in self.Gt.nodes[self.Gr.nodes[nr]['cycle']]['vnodes']:
            ncy = self.Gr.nodes[nr]['cycle']
            if seg in self.Gt.nodes[ncy]['cycle'].cycle:
                rooms.append(nr)
        return rooms

    def room2segments(self, room):
        """ returns the segments of a room

        Parameters
        ----------

        room : int

        Returns
        -------

        seg : list

        """

        try:
            # old vnodes was there
            ncy = self.Gr.nodes[room]['cycle']
            seg = self.Gt.nodes[ncy].cycle
        except:
            raise NameError(str(room) + " is not in not on Gr")
        u = np.where(seg >= 0)
        seg = seg[u]
        return np.sort(seg.tolist())

    def room2nodes(self, room):
        """ returns the nodes of a room

        Parameters
        ----------

        room : int

        Returns
        -------

        nod : sorted list

        """

        try:
            ncy = self.Gr.nodes[room]['cycle']
            nod = self.Gt.nodes[ncy].cycle
            #nod = self.Gt.nodes[self.Gr.nodes[room]['cycle']]['vnodes']
        except:
            raise NameError(str(room) + " is not in not on Gr")
        u = np.where(nod < 0)
        nod = nod[u]

        return np.sort(nod.tolist())

    def get_diffslab(self,npt,lz):
        """ get the 2 slabs associated to a diffraction point

        Parameters
        ----------

        lnpt : diffraction point numbers (node of Gs)
        lz   : array of candidate heights of the diffraction point

        Notes
        -----

        As a diffraction point may involve iso segments the nature
        of the diffraction interaction depends on a height parameter
        This function extacts the couple of slab from this information

        Returns
        -------

        - a list of 2-segments . the length of this list == length of lz
        - a list of slab tuples.  the length of this list == length of lz

        [[443, 529], [444, 530]]
        [['WALL', 'WALL'], ['AIR', 'AIR']]

        """
        assert(npt in self.ddiff), logger.error('npt not a diffraction point')
        lcy = self.ddiff[npt][0]
        ls = []
        llz = len(lz)
        dz_seg= {z:[] for z in range(llz)}
        dz_sl= {z:[] for z in range(llz)}

        for cy in lcy:
            vn = set(self.Gt.nodes[cy]['polyg'].vnodes)
            # v1.1 lneig_pt = set(nx.neighbors(self.Gs,npt))
            lneig_pt = set(self.Gs[npt])
            lseg = lneig_pt.intersection(vn)
            lseg_valid = [ x for x in lseg if self.Gs.nodes[x]['name']!='_AIR']

            for x in lseg_valid:
                zsup = lz >self.Gs.nodes[x]['z'][0]
                zinf = lz <=self.Gs.nodes[x]['z'][1]
                z    = zsup & zinf
                uz = np.where(z)[0]
                # fill dz_seg at the correct height with a lseg_valid
                # and simulnaneously
                # fill dz_sl at the correct height with correspondong slab
                [(dz_seg[i].append(x),dz_sl[i].append(self.Gs.nodes[x]['name']))
                                                                    for i in uz]

        lseg  = list(dz_seg.values())
        lslab = list(dz_sl.values())
        return lseg, lslab

    def _find_diffractions(self, difftol=0.01,verbose = False,tqdmkwargs={}):
        """ find diffractions points of the Layout

        Parameters
        ----------

        difftol : float

            tolerance in radians

        Returns
        -------

        Update self.ddiff {nseg : ([ncy1,ncy2],wedge_angle)}

        """
        # dangles = self.get_Gt_angles()
        #
        # Problem here point number are converted into float64

        if tqdmkwargs=={}:
            tqdmkwargs={'total':100.,
                        'desc':'find_diffractions'}

        dangles = {cy: np.array(geu.get_pol_angles(self.Gt.nodes[cy]['polyg']))
                   for cy in self.Gt.nodes() if cy != 0}

        #
        # The candidate points for being diffraction points have degree 1 or 2
        # A point diffracts toward one or several cycles
        #
        #ldiff = list(np.hstack((self.degree[1],self.degree[2])).astype('int'))
        lpnt = [x for x in self.Gs.node if (x < 0 and x not in self.degree[0])]

        self.ddiff = {}

        if verbose :
            cpt = 1./(len(lpnt)+1)
            pbar = tqdm.tqdm(tqdmkwargs)
        for k in lpnt:
            if verbose :
                pbar.update(100.*cpt)
            # list of cycles associated with point k
            lcyk = self.Gs.nodes[k]['ncycles']
            if len(lcyk) > 2:
                # Subgraph of connected cycles around k
                Gtk = nx.subgraph(self.Gt, lcyk)
                # ordered list of connections between cycles
                try:
                    lccyk = nx.find_cycle(Gtk)
                except:
                    pdb.set_trace()

                # list of segment neighbours
                neigh = list(dict(self.Gs[k]).keys())
                # sega : list of air segment in neighors
                sega = [n for n in neigh if
                        (self.Gs.nodes[n]['name'] == 'AIR' or
                         self.Gs.node[n]['name'] == '_AIR')]

                sega_iso = [n for n in sega if len(self.Gs.node[n]['iso']) > 0]
                sega_eff = list(set(sega).difference(set(sega_iso)))
                nsector = len(neigh) - len(sega)

                dsector = {i: [] for i in range(nsector)}
                #
                # team building algo
                #
                ct = 0
                # if k ==-44:
                #     pdb.set_trace()
                for ccy in lccyk:

                    #segsep = self.Gt[ccy[0]][ccy[1]]['segment'][0]
                    segsep = self.Gt[ccy[0]][ccy[1]]['segment']
                    # filter only segments connected to point k (neigh)
                    lvseg = [x for x in segsep if x in neigh]
                    if len(lvseg) == 1 and (lvseg[0] in sega_eff):  # same sector
                        dsector[ct].append(ccy[1])
                    else:  # change sector
                        ct = (ct + 1) % nsector
                        dsector[ct].append(ccy[1])

                    # typslab = self.Gs.node[segsep]['name']
                    # if (typslab=='AIR' or typslab=='_AIR'): # same sector
                        # dsector[ct].append(ccy[1])
                    # else: # change sector
                        # ct=(ct+1)%nsector
                        # dsector[ct].append(ccy[1])
                        # lcy2.append(ccy[1])
                        # lcy1,lcy2 = lcy2,lcy1

                dagtot = {s: 0 for s in range(nsector)}
                save = []
                for s in dsector:
                    for cy in dsector[s]:
                        da = dangles[cy]
                        u = np.where(da[0, :].astype('int') == k)[0][0]
                        save.append((cy, da[1, u]))
                        dagtot[s] = dagtot[s] + da[1, u]
                for s in dagtot:
                    if dagtot[s] > (np.pi + difftol):
                        self.ddiff[k] = (dsector[s], dagtot[s])
                        break

                # if agtot1 > (np.pi+tol):
                #     self.ddiff[k]=(lcy1,agtot1)
                # elif 2*np.pi-agtot1 > (np.pi+tol):
                #     self.ddiff[k]=(lcy2,2*np.pi-agtot1)
            else:
                # diffraction by half-plane detected
                if k in self.degree[1]:
                    self.ddiff[k] = (lcyk, 2 * np.pi)


    def buildGr(self):
        """ build the graph of rooms Gr

        Notes
        -----

        adjascent rooms are connected

        Gr is at startup a deep copy of Gt

        The difficulty here is to take into account the AIR transition
        segments

        """

        self.Gr = copy.deepcopy(self.Gt)
        self.Gr.remove_node(0)
        self.Gr.remove_edges_from(self.Gt.edges())

        for e in list(self.Gt.edges()):
            if ((not 0 in e) and
                (self.Gt.node[e[0]]['indoor']) and
                (self.Gt.node[e[1]]['indoor']) ):

                seg = self.Gt[e[0]][e[1]]['segment']
                seg = np.unique(seg)
                trans_seg = [n for n in seg
                             if (self.Gs.node[n]['transition'])
                             and n not in self.segboundary]
                if trans_seg != []:
                    self.Gr.add_edge(e[0],e[1],segment=trans_seg)

        deg = dict(self.Gr.degree())
        #pdb.set_trace()
        self.Gr.remove_nodes_from([n for n in deg if deg[n] == 0])

    def buildGw(self):
        """ build Graph of waypaths


        See Also
        --------

        buildGr

        Notes
        -----

        for all edges of Gr (adjascent room)
            if room1 and room2 have a common transition


        """

        self.Gw = nx.Graph(name='Gw')
        self.Gw.pos = {}

        d_id = max(self.Gr.nodes())  # for numerotation of Gw nodes
        d_id_index = d_id + 1

        for e in self.Gr.edges():  # iterator on Gr edges

            self.Gw.add_node(e[0], room=e[0], door=False)
            self.Gw.add_node(e[1], room=e[1], door= False)

            # transitions of room e[0]
            # trans1 = self.Gr.node[e[0]]['segment']
            # # transitions of room e[1]
            # trans2 = self.Gr.node[e[1]]['segment']
            # Id = np.intersect1d(trans1, trans2)[0]  # list of common doors
            # import ipdb
            # ipdb.set_trace()
            Ids = self.Gr[e[0]][e[1]]['segment']
            # here is supposed that 2 room may have more than 1 door in common
            for Id in Ids:
                #v1.1 unode = self.Gs.neighbors(Id)  # get edge number of common doors
                unode = list(dict(self.Gs[Id]).keys())  # get edge number of common doors
                up0 = self.Gs.pos[unode[0]]
                up1 = self.Gs.pos[unode[1]]

                name = self.Gs.node[Id]['name']
                pn = self.Gs.node[Id]['norm']
                sl = self.sl[name]
                thick = (sum(sl['lthick']) / 2.) + 0.2

                # for ""doors"" extra waypoints points are added
                # in front and back of the aperture.
                # this is not done for AIR slabs
                if 'AIR' not in name:

                    # middle of the common door
                    pdoor0 = (np.array(up0) + pn[:2] * thick +
                              np.array(up1) + pn[:2] * thick) / 2.
                    pdoor1 = (np.array(up0) - pn[:2] * thick +
                              np.array(up1) - pn[:2] * thick) / 2.
                    P0 = sh.Point(pdoor0)
                    P1 = sh.Point(pdoor1)

                    ep0 = self.Gr.pos[e[0]]
                    ep1 = self.Gr.pos[e[1]]

                    if self.Gr.node[e[0]]['polyg'].contains(P0):
                        upd0 = d_id_index
                        self.Gw.pos[upd0] = pdoor0
                        self.Gw.add_node(upd0, room=e[0], door=True)
                        # if self.seginline(pdoor0,ep0).shape[1] <= 1:
                        self.Gw.add_edges_from([(e[0],upd0)])
                        d_id_index = d_id_index + 1

                        upd1 = d_id_index
                        self.Gw.pos[upd1] = pdoor1
                        self.Gw.add_node(upd1, room=e[1], door=True)
                        # if self.seginline(pdoor1,ep1).shape[1] <= 1:
                        self.Gw.add_edges_from([(e[1],upd1)])
                        d_id_index = d_id_index + 1
                    else:
                        upd0 = d_id_index
                        self.Gw.pos[upd0] = pdoor0
                        self.Gw.add_node(upd0, room=e[1], door=True)
                        # if self.seginline(pdoor0,ep1).shape[1] <= 1:
                        self.Gw.add_edges_from([(e[1],upd0)])
                        d_id_index = d_id_index + 1

                        upd1 = d_id_index
                        self.Gw.pos[upd1] = pdoor1
                        self.Gw.add_node(upd1, room=e[0], door=True)
                        # if self.seginline(pdoor1,ep0).shape[1] <= 1:
                        self.Gw.add_edges_from([(e[0],upd1)])
                        d_id_index = d_id_index + 1
                    self.Gw.add_edges_from([(upd0, upd1)])
                else:
                    self.Gw.add_edges_from([(e[0],e[1])])
        self.Gw.pos.update(self.Gr.pos)




    def info(self):
        """ gives information about the Layout
        """
        print("filestr : ", self._filename)
        # print("filematini : ", self.filematini)
        # print("fileslabini : ", self.fileslabini)
        try:
            print("filegeom : ", self.filegeom)
        except:
            print("geomfile (.off) has no been generated")

        # self.boundary()
        print("boundaries ", self.ax)
        print("number of Points :", self.Np)
        print("number of Segments :", self.Ns)
        print("number of Sub-Segments :", self.Nss)
        try:
            print("Gs Nodes : ", self.Gs.number_of_nodes())
            print("Gs Edges : ", self.Gs.number_of_edges())
        except:
            print("no Gs graph")

        try:
            print("Gt Nodes : ", self.Gt.number_of_nodes())
            print("Gt Edges : ", self.Gt.number_of_edges())
            print("vnodes = Gt.node[Nc]['polyg'].vnodes")
            print("poly = Gt.node[Nc]['polyg']")
        except:
            print("no Gt graph")

        try:
            print("Gr Nodes    :", self.Gr.number_of_nodes())
            print("Gr Edges    :", self.Gr.number_of_edges())
        except:
            print("no Gr graph")

    def facets3D(self, edlist, name='Layer', subseg=False):
        """ create facet 3D for geomview

        Parameters
        ----------

        edlist
        name : string
        subseg : boolean 


        """

        filename = name + '.list'
        filestruc = pyu.getlong(filename, pro.pstruc['DIRGEOM'])
        fos = open(filestruc, "w")
        fos.write("LIST{\n")
        for e in edlist:
            filename = self.facet3D(e, subseg)
            if filename == 'void':
                pass
            else:
                chaine = '{<' + filename + "}\n"
                fos.write(chaine)

        fos.write("}\n")
        fos.close()

    def numseg(self, ta, he, first=True):
        """ get segment number from 2 points index

        Parameters
        ----------

        ta  : int <0
        he  : int <0

        first : Boolean
            if True returns only one among the several iso segments
            else returns a np.array of iso segments

        Returns
        -------

        nseg : > 0
        if 0 not a segment

        See Also
        --------

        buildGt

        """
        nta = np.array(list(dict(self.Gs[ta]).keys()))
        nhe = np.array(list(dict(self.Gs[he]).keys()))
        nseg = np.intersect1d(nta, nhe)
        if len(nseg > 0):
            if first:
                return(nseg[0])
            else:
                return nseg
        else:
            return(0)

    def isseg(self, ta, he):
        """ test if ta<->he is a segment

        Parameters
        ----------

        ta  : int <0
        he  : int <0

        Returns
        -------

        boolean

        See Also
        --------

        editor.py

        """
        # transpose point numbering

        upnt = [ x for x in self.Gs.nodes() if x < 0 ]
        ta = np.nonzero(np.array(upnt) == ta)[0][0]
        he = np.nonzero(np.array(upnt) == he)[0][0]
        res = [x for x in zip(self.tahe[0], self.tahe[1])
              if (((x[0] == ta) & (x[1] == he)) |
                  ((x[0] == he) & (x[1] == ta)))  ]

        if len(res) > 0:
            return True
        else:
            return False

    def ispoint(self, pt, tol=0.05):
        """ check if pt is a point of the Layout

        Parameters
        ----------

        pt  : point (2,1)
        tol : float
            default (0.05 meters)

        if True the point number (<0) is returned
        else 0 is return

        Returns
        -------

        pt : point number if point exists 0 otherwise

        See Also
        --------

        pylayers.util.geomutil.Polygon.setvnodes

        """
        # print"ispoint : pt ", pt
        pts = np.array(list(self.Gs.pos.values())).T
        ke = np.array(list(self.Gs.pos.keys()))
        diff = pts - pt.reshape(2, 1)
        v = np.sqrt(np.sum(diff * diff, axis=0))
        nz = (v > tol)
        b = nz.prod()
        if b == 1:
            # if all layout points are different from pt
            #return(0,np.min(v))
            return(0)
        else:
            nup = np.where(nz == False)[0]
            if len(nup) == 1:
                return(ke[nup][0])
            else:
                mi = np.where(min(v[nup]) == v[nup])[0]
                return(ke[nup[mi]][0])

    def onseg(self, pt, tol=0.01):
        """ segment number from point (deprecated)

        return segment number which contains point pt

        Parameters
        ----------

        pt  np.array(1x2)
        tol = 0.01      tolerance

        """

        pts = np.array(self.Gs.pos.values()).T   # structure points
        ke = np.array(list(dict(self.Gs.pos).keys()))        # point keys
        n = np.shape(pts)[1]
        nbu = np.array([])
        if (n > 0):
            num = np.arange(n)                   #
            b = self.inbox(pt, tol)

            ta = self.tahe[0, b]
            he = self.tahe[1, b]

            nb = num[b]

            n = len(nb)
            p = np.outer(pt, np.ones(n))

            # printta
            v1 = p - pts[:, ta]
            v2 = pts[:, he] - p

            nv1 = np.sqrt(v1[0, :] * v1[0, :] + v1[1, :] * v1[1, :])
            nv2 = np.sqrt(v2[0, :] * v2[0, :] + v2[1, :] * v2[1, :])

            v1n = v1 / nv1
            v2n = v2 / nv2

            ps = v1n[0, :] * v2n[0, :] + v1n[1, :] * v2n[1, :]
            u = abs(1. - ps) < tol
            nbu = nb[u]

        return nbu

    def facet3D(self, e, subseg=False):
        """ calculate 3D facet from segment


        Parameters
        ---------

        s : int
            segment number
        subseg : boolean
            default False

        """
        P1 = np.array(np.zeros(3), dtype=np.float64)
        P2 = np.array(np.zeros(3), dtype=np.float64)
        P3 = np.array(np.zeros(3), dtype=np.float64)
        P4 = np.array(np.zeros(3), dtype=np.float64)
        # v1.1 nebr = self.Gs.neighbors(s)
        nebr = list(dict(self.Gs[s]).keys())
        n1 = nebr[0]
        n2 = nebr[1]

        P1[0:2] = np.array(self.Gs.pos[n1])
        P1[2] = self.Gs.node[s]['z'][0]

        P2[0:2] = np.array(self.Gs.pos[n2])
        P2[2] = self.Gs.node[s]['z'][0]

        P3[0:2] = np.array(self.Gs.pos[n2])
        P3[2] = self.Gs.node[s]['z'][1]

        P4[0:2] = np.array(self.Gs.pos[n1])
        P4[2] = self.Gs.node[s]['z'][1]

        cold = pyu.coldict()

        if subseg:
            nsseg = len(self.Gs.node[s]['ss_name'])
        else:
            nsseg = 0

        filename = 'fa' + str(s) + '.off'
        filestruc = pyu.getlong(filename, pro.pstruc['DIRGEOM'])
        fos = open(filestruc, "w")
        fos.write("OFF\n")
        fos.write("%d %d \n\n" % (1 + (nsseg + 1) * 4, nsseg + 1))
        fos.write("0.000 0.000 0.000\n")
        if subseg:
            try:
                for k, name in enumerate(self.Gs.node[s]['ss_name']):
                    P1[2] = self.Gs.node[s]['ss_z'][k][0]
                    P2[2] = self.Gs.node[s]['ss_z'][k][0]
                    P3[2] = self.Gs.node[s]['ss_z'][k][1]
                    P4[2] = self.Gs.node[s]['ss_z'][k][1]
                    fos.write("%6.3f %6.3f %6.3f \n" % (P1[0], P1[1], P1[2]))
                    fos.write("%6.3f %6.3f %6.3f \n" % (P2[0], P2[1], P2[2]))
                    fos.write("%6.3f %6.3f %6.3f \n" % (P3[0], P3[1], P3[2]))
                    fos.write("%6.3f %6.3f %6.3f \n" % (P4[0], P4[1], P4[2]))
            except:
                print('no subsegment on ', s)
                return('void')
        else:
            name = self.Gs.node[s]['name']
            fos.write("%6.3f %6.3f %6.3f \n" % (P1[0], P1[1], P1[2]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P2[0], P2[1], P2[2]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P3[0], P3[1], P3[2]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P4[0], P4[1], P4[2]))

        if subseg:
            for k, name in enumerate(self.Gs.node[s]['ss_name']):
                colname = sl[name]['color']
                colhex = cold[colname]
                col = pyu.rgb(colhex) / 255.
                fos.write("4 %i %i %i %i %6.3f %6.3f %6.3f 0.4\n" % (1 + 4 * k, 2 + 4 * k,
                                                                     3 + 4 * k, 4 + 4 * k, col[0], col[1], col[2]))
        else:
            name = self.Gs.node[s]['name']
            colname = sl[name]['color']
            colhex = cold[colname]
            col = pyu.rgb(colhex) / 255.
            fos.write("4 %i %i %i %i %6.3f %6.3f %6.3f 0.4\n" % (1, 2,
                                                                 3, 4, col[0], col[1], col[2]))

        return(filename)

    def geomfile(self, centered=False):
        """ create a .off geomview file

        Parameters
        ----------

        centered : Boolean
            if True the layout is centered around its center of gravity

        Notes
        -----

        The `.off` file can be vizualized through the show3 method

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('DLR.lay')
        >>> pg = L.geomfile()

        """

        # calculate center of gravity
        if centered:
            pg = np.sum(self.pt, axis=1) / np.shape(self.pt)[1]
        else:
            pg = np.array([0, 0])

        # en  = self.Ns # number of segments
        en = len(np.where(np.array(list(dict(self.Gs.node).keys())) > 0)[0])
        if en != self.Ns:
            logger.warning("wrong number of segments, consistency problem in layout")
        #cen = self.Nss
        # d : dictionnary of layout sub segments
        #
        d = self.subseg()
        cen = 0
        for k in d:
            lss = d[k]
            cen = cen + len(lss)

        if cen != self.Nss:
            logger.warning("wrong number of subsegments, consistency problem in layout")

        sl = self.sl
#
#        Create a polygon for each segment and subsegment
#
        P1 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P2 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P3 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P4 = np.array(np.zeros([3, en + cen], dtype=np.float64))

        ik = 0
        dikn = {}
        for i in list(dict(self.Gs.node).keys()):
            if i > 0:  # segment
                if ((self.Gs.node[i]['name'] != 'AIR') and
                        (self.Gs.node[i]['name'] != '_AIR')):
                    #v1.1 nebr = self.Gs.neighbors(i)
                    nebr = list(dict(self.Gs[i]).keys())
                    n1 = nebr[0]
                    n2 = nebr[1]
                    P1[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                    P1[2, ik] = self.Gs.node[i]['z'][0]

                    P2[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                    P2[2, ik] = self.Gs.node[i]['z'][0]

                    P3[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                    P3[2, ik] = self.Gs.node[i]['z'][1]

                    P4[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                    P4[2, ik] = self.Gs.node[i]['z'][1]
                    dikn[ik] = i
                    ik = ik + 1
                else:
                    en = en - 1
        # d = self.subseg()
        # k : ss_name v: seg number
        cpt = 0
        subseg = {}
        # pdb.set_trace()
        for k in d.keys():
            for l in d[k]:
                ids = l[0]
                subseg[cpt] = ids
                order = l[1]
                cpt = cpt + 1
                # v1.1 nebr = self.Gs.neighbors(l[0])
                nebr = list(dict(self.Gs[l[0]]).keys())
                n1 = nebr[0]
                n2 = nebr[1]
                # printik,n1,n2

                P1[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                P1[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                # printP1[:,ik]

                P2[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                P2[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                # printP2[:,ik]

                P3[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                P3[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                # printP3[:,ik]

                P4[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                P4[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                # printP4[:,ik]

                dikn[ik] = l
                ik = ik + 1

        npt = 4 * (en + cen)
        _filename, ext = os.path.splitext(self._filename)
        _filegeom = _filename + '.off'
        self.filegeom = _filegeom
        filegeom = pyu.getlong(_filegeom, pro.pstruc['DIRGEOM'])
        fos = open(filegeom, "w")
        fos.write("OFF\n")
        fos.write("%d %d \n\n" % (npt + 1, en + cen))
        fos.write("0.000 0.000 0.000\n")
        for i in range(en + cen):
            fos.write("%6.3f %6.3f %6.3f \n" % (P1[0, i], P1[1, i], P1[2, i]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P2[0, i], P2[1, i], P2[2, i]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P3[0, i], P3[1, i], P3[2, i]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P4[0, i], P4[1, i], P4[2, i]))

        cold = pyu.coldict()
#        ke   = cold.keys()
#

        for i in range(en + cen):
            q = 4 * i
            if i < en:
                #ne = i + 1
                ne = dikn[i]
                name = self.Gs.node[ne]['name']
            else:
                ne = dikn[i][0]
                order = dikn[i][1]
                #nss = i - en
                ##ne = subseg[nss]
                name = self.Gs.node[ne]['ss_name'][order]

#            if (i<en):
#                name = self.name[i]
#            else:
#                core = self.ce[subseg[i-en]][0]
#                name = sl.di[core]

            colname = sl[name]['color']

            colhex = cold[colname]
            col = pyu.rgb(colhex) / 255.
            fos.write("4 %i %i %i %i %6.3f %6.3f %6.3f 0.4\n" % (q +
                                                                 1, q + 2, q + 3, q + 4, col[0], col[1], col[2]))
        fos.close()
        return pg

    def _show3(self, centered=False, newfig=False, opacity=1., ceil_opacity=1., show_ceil=False, cyid=False, **kwargs):
        """ mayavi 3D vizualisation

        Parameters
        ----------

        newfig : Boolean
            create a new mayavi Figure
        opacity : float ([0,1])
            set slab opacity
        ceil_opacity : float
        centered : Boolean
            if True the layout is centered around its center of gravity
        cyid : boolean
            display cycle number
        show_ceil: boolean
            display ceil or not

        Notes
        -----

        The `.off` file can be vizualized through the show3 method

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.gis.layout import *
            >>> L = Layout()

        """

        #
        # calculate center of gravity of the layout
        #
        if centered:
            pg = np.sum(self.pt, axis=1) / np.shape(self.pt)[1]
        else:
            pg = np.array([0, 0])

        # en  = self.Ns # number of segments
        en = len(np.where(np.array(list(dict(self.Gs.node).keys())) > 0)[0])
        if en != self.Ns:
            logger.warning(
                "wrong number of segment consistency problem in layout")
        #cen = self.Nss
        # d : dictionnary of layout sub segments
        #
        d = self.subseg()
        cen = 0
        for k in d:
            lss = d[k]
            cen = cen + len(lss)

        if cen != self.Nss:
            logger.warning(
                "wrong number of subsegment consistency problem in layout")

        sl = self.sl
#
#        Create a 3D polygon for each segment and subsegment
#
        P1 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P2 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P3 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P4 = np.array(np.zeros([3, en + cen], dtype=np.float64))

        ik = 0
        dikn = {}

        #
        # segments which are not _AIR or AIR
        #
        for i in list(dict(self.Gs.node).keys()):
            if i > 0:  # segment
                if ((self.Gs.node[i]['name'] != 'AIR') and
                        (self.Gs.node[i]['name'] != '_AIR')):
                    #v1.1 nebr = self.Gs.neighbors(i)
                    nebr = list(dict(self.Gs[i]).keys())
                    n1 = nebr[0]
                    n2 = nebr[1]
                    P1[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                    P1[2, ik] = self.Gs.node[i]['z'][0]

                    P2[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                    P2[2, ik] = self.Gs.node[i]['z'][1]

                    P3[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                    P3[2, ik] = self.Gs.node[i]['z'][1]

                    P4[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                    P4[2, ik] = self.Gs.node[i]['z'][0]

                    dikn[ik] = i
                    ik = ik + 1

                else:

                    en = en - 1

        # d = self.subseg()
        # k : ss_name v: seg number
        cpt = 0
        subseg = {}
        for k in d.keys():
            for l in d[k]:
                ids = l[0]
                subseg[cpt] = ids
                order = l[1]
                cpt = cpt + 1
                # v1.1 nebr = self.Gs.neighbors(l[0])
                nebr = list(dict(self.Gs[l[0]]).keys())
                n1 = nebr[0]
                n2 = nebr[1]
                # printik,n1,n2

                P1[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                P1[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                # printP1[:,ik]

                P2[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                P2[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                # printP2[:,ik]

                P3[0:2, ik] = np.array(self.Gs.pos[n2]) - pg
                P3[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                # printP3[:,ik]

                P4[0:2, ik] = np.array(self.Gs.pos[n1]) - pg
                P4[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                # printP4[:,ik]

                dikn[ik] = l
                ik = ik + 1

        npt = 4 * (en + cen)
        npt_s = (en + cen)

        points = np.hstack((P1[:, 0:npt_s], P2[:, 0:npt_s]))
        points = np.hstack((points, P3[:, 0:npt_s]))
        points = np.hstack((points, P4[:, 0:npt_s]))
        points = points.T

        boxes = np.empty((int(npt / 4), 4), dtype='int')
        b = np.arange(int(npt / 4))

        boxes[:, 0] = b
        boxes[:, 1] = b + npt_s
        boxes[:, 2] = b + 2 * npt_s
        boxes[:, 3] = b + 3 * npt_s

#         _filename,ext = os.path.splitext(self._filename)
#         _filegeom = _filename+'.off'
#         self.filegeom=_filegeom
#         filegeom = pyu.getlong(_filegeom, pro.pstruc['DIRGEOM'])
#         fos = open(filegeom, "w")
#         fos.write("OFF\n")
#         fos.write("%d %d \n\n" % (npt + 1, en + cen))
#         fos.write("0.000 0.000 0.000\n")
#         for i in range(en + cen):
#             fos.write("%6.3f %6.3f %6.3f \n" % (P1[0, i], P1[1, i], P1[2, i]))
#             fos.write("%6.3f %6.3f %6.3f \n" % (P2[0, i], P2[1, i], P2[2, i]))
#             fos.write("%6.3f %6.3f %6.3f \n" % (P3[0, i], P3[1, i], P3[2, i]))
#             fos.write("%6.3f %6.3f %6.3f \n" % (P4[0, i], P4[1, i], P4[2, i]))

        cold = pyu.coldict()
        color = np.zeros((4 * (cen + en), 3))
        for i in range(en + cen):
            # q = 4 * i
            if i < en:
                ne = dikn[i]
                name = self.Gs.node[ne]['name']
            else:
                ne = dikn[i][0]
                order = dikn[i][1]
                name = self.Gs.node[ne]['ss_name'][order]

            colname = sl[name]['color']
            colhex = cold[colname]
            color[i, :] = pyu.rgb(colhex)
            color[i + npt_s, :] = pyu.rgb(colhex)
            color[i + 2 * npt_s, :] = pyu.rgb(colhex)
            color[i + 3 * npt_s, :] = pyu.rgb(colhex)

        colname = sl['FLOOR']['color']
        colhex = cold[colname]
        colf = np.repeat((pyu.rgb(colhex))[np.newaxis, :], 4, axis=0)
        color = np.vstack((color, colf))

        # trick for correcting  color assignement

        sc = tvtk.UnsignedCharArray()
        sc.from_array(color)
 
        # manage floor

        # if Gt doesn't exists

        try:
            self.ma.coorddeter()
            # z=np.ones(self.ma.xy.shape[1])
            z = np.zeros(self.ma.xy.shape[1])
            F = np.vstack((self.ma.xy, z))

            tri = np.arange(len(z))
            meshf = tvtk.PolyData(points=F.T, polys=np.array([tri]))

            meshf.point_data.scalars = sc
            meshf.point_data.scalars.name = 'scalars'

            surff = mlab.pipeline.surface(meshf, opacity=opacity)
            mlab.pipeline.surface(mlab.pipeline.extract_edges(surff),
                                  color=(0, 0, 0), )

        # otherwise
        except:

            floorx = np.array((points[:, 0].min(), points[:, 0].max()))
            floory = np.array((points[:, 1].min(), points[:, 1].max()))
            zmin = np.min(points[:, 2])
            Pf = np.array([floorx[0], floory[0], zmin])
            Pf = np.vstack((Pf, np.array([floorx[0], floory[1], zmin])))
            Pf = np.vstack((Pf, np.array([floorx[1], floory[1], zmin])))
            Pf = np.vstack((Pf, np.array([floorx[1], floory[0], zmin])))

            points = np.vstack((points, Pf))
            bf = np.arange(npt, npt + 4)
            boxes = np.vstack((boxes, bf))

        mesh = tvtk.PolyData(points=points, polys=boxes)
        mesh.point_data.scalars = sc
        mesh.point_data.scalars.name = 'scalars'

        if newfig:
            mlab.clf()
            f = mlab.figure(bgcolor=(1, 1, 1))
        else:
            f = mlab.gcf()
            f.scene.background = (1, 1, 1)

        f.scene.disable_render = True

        surf = mlab.pipeline.surface(mesh, opacity=opacity)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                              color=(0, 0, 0), )
        f.children[-1].name = 'Layout ' + self._filename

        if show_ceil == True:
            if len(self.Gt.nodes()) != 0:
                uin = [kn for kn in self.Gt.nodes() if self.Gt.node[kn]
                       ['indoor'] == True]

                ptc = np.ndarray(shape=(3, 0))
                boxc = np.ndarray(shape=(0, 3))
                cpt = 0
                for u in uin:
                    p = self.Gt.node[u]['polyg']
                    no = self.Gt.node[u]['polyg'].vnodes[
                        self.Gt.node[u]['polyg'].vnodes > 0]
                    for n in no:
                        if self.Gs.node[n]['z'][1] != 40000000:
                            h = self.Gs.node[n]['z'][1]
                            break
                    vert = {"vertices": np.array(p.exterior.xy).T}
                    dt = triangle.triangulate(vert)
                    nbpt = len(dt['vertices'])
                    pt = np.vstack((dt['vertices'].T, [h] * nbpt))
                    box = dt['triangles']

                    # if u == 114:
                    #     import ipdb
                    #     ipdb.set_trace()
                    #     box = np.roll(box,1,1)
                    ptc = np.hstack((ptc, pt))
                    boxc = np.vstack((boxc, box + cpt))
                    cpt = cpt + nbpt
                    # if box.shape[0] == 2 :
                    #     import ipdb
                    #     ipdb.set_trace()
                    #     print(cpt,nbpt)
                    #     print(box)
                    #     print(pt)
                    #     break

                # manage Ceil color

                colname = sl['CEIL']['color']
                colhex = cold[colname]
                colf = np.repeat((pyu.rgb(colhex))[np.newaxis, :], cpt, axis=0)
                # color = np.vstack((color, colf))
                color=colf

                # trick for correcting  color assignement

                sc = tvtk.UnsignedCharArray()
                sc.from_array(color)

                meshc = tvtk.PolyData(points=ptc.T, polys=boxc)
                meshc.point_data.scalars = sc
                meshc.point_data.scalars.name = 'scalars'
                mlab.pipeline.surface(
                    meshc, opacity=ceil_opacity, reset_zoom=False)

                # ptc =

                # ptcxy = np.array([self.Gt.node[u]['polyg'].exterior.xy[0],self.Gt.node[u]['polyg'].exterior.xy[1]])
                # ptcz = [self.Gs.node[self.Gt.node[u]['polyg'].vnodes[1]]['z'][1]]*len(self.Gt.node[u]['polyg'].exterior.xy[0])
                # ptc = np.vstack((ptcxy,ptcz))
                # nbpt = ptc.shape[1]
                # pdb
                # ceil = tvtk.PolyData(points=ptc.T, polys=np.arange(nbpt).reshape(1,nbpt))
                # surf2 = mlab.pipeline.surface(ceil, opacity=opacity)
                # import ipdb
                # ipdb.set_trace()

        if cyid:
            if len(self.Gt.nodes()) > 0:
                pk = self.Gt.pos.keys()
                v = np.array(self.Gt.pos.values())
                [mlab.text3d(v[ik, 0], v[ik, 1], 0.5, str(k))
                 for ik, k in enumerate(pk)]
        # if segpt:

        #     seg = dict(filter(lambda x: x[0]>0,self.Gs.pos.items()))
        #     pt = dict(filter(lambda x: x[0]<0,self.Gs.pos.items()))
        #     pseg = np.array(seg.values())
        #     ppt = np.array(pt.values())
        #     [mlab.text3d(pseg[ik,0],pseg[ik,1],0.5,str(k)) for ik,k in enumerate(seg)]
        #     [mlab.text3d(ppt[ik,0],ppt[ik,1],3.,str(k)) for ik,k in enumerate(pt)]

        f.scene.disable_render = False
        return(f)

    def show3(self, bdis=True, centered=True):
        """ geomview display of the indoor structure

        Parameters
        ----------

        bdis boolean (default True)
            boolean display (call geowview if True)
        centered : boolean
            if True center the layout before display


        """

        pg = self.geomfile(centered=centered)

        filename = pyu.getlong(self.filegeom, pro.pstruc['DIRGEOM'])
        if (bdis):
            #chaine = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
            chaine = "geomview  -b 1 1 1 " + filename + " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)

        return(pg)

    def signature(self, iTx, iRx):
        """ Determine signature between node iTx and node iRx

        Parameters
        ----------

        cy1  : int
            source cycle
        cy2  : int
            target cycle

        Returns
        -------

        sigarr    :
        signature :


        Warnings
        --------
        This a temporary function
            There is some algorithmic work to find the best way to determine signature
            T4 : limit the ndt to only edges and nodes in visibility from Tx

        """
        # Here we take all the vnodes >0  from the room
        #
        # Practically those list of nodes should depend on pTx , pRx
        #
        try:
            self.Gi
        except:
            raise NameError(
                'Interaction graph layout.Gi must be build before signature computation')
        if isinstance(iTx, np.ndarray):
            NroomTx = self.pt2ro(iTx)
        elif isinstance(iTx, int):
            NroomTx = iTx
        else:
            raise NameError('iTx must be an array or a room number')
        if isinstance(iRx, np.ndarray):
            NroomRx = self.pt2ro(iRx)
        elif isinstance(iRx, int):
            NroomRx = iRx
        else:
            raise NameError('iRx must be an array or a room number')

        if not self.Gr.has_node(NroomTx) or not self.Gr.has_node(NroomRx):
            raise AttributeError('Tx or Rx is not in Gr')

        #
        # .. todo:: modifier inter afin de ne pas retenir les points non diffractants
        #
        ndt = self.Gt.node[self.Gr.node[NroomTx]['cycle']]['inter']
        ndr = self.Gt.node[self.Gr.node[NroomRx]['cycle']]['inter']
        sigarr = np.array([]).reshape(2, 0)
        for nt in ndt:
            for nr in ndr:
                addpath = False
                if (type(nt) != type(nr)):
                    try:
                        path = nx.dijkstra_path(self.Gi, nt, nr)
                        addpath = True
                    except:
                        pass
                        # print'no path between ',nt,nr
                elif (nt != nr):
                    try:
                        path = nx.dijkstra_path(self.Gi, nt, nr)
                        addpath = True
                    except:
                        pass
                        # print'no path between ',nt,nr
                else:
                    addpath = True
                    path = [nt]
                if addpath:
                    sigarr = np.hstack((sigarr, np.array([[0], [0]])))
                    for interaction in path:
                        it = eval(interaction)
                        if type(it) == tuple:
                            sigarr = np.hstack((sigarr,
                                                np.array([[it[0]], [1]])))
                        elif it < 0:
                            sigarr = np.hstack((sigarr,
                                                np.array([[it], [-1]])))
                        else:
                            sigarr = np.hstack((sigarr, np.array([[it], [2]])))

        return sigarr

    def plot(self, **kwargs):
        """ plot the layout with shapely MultiLineString

        Parameters
        ---------

        show : boolean
        fig :figure
        ax  :
        labels : list
        nodes : boolean

        Returns
        -------

        fig, ax

        Examples
        --------

        >>> L= Layout('Munich.lay',bbuild=False)
        >>> L.plot(show=True)

        """

        #fig = kwargs.pop('fig', plt.gcf())
        #ax = kwargs.pop('ax', plt.gca())


        fig, ax = plt.subplots(facecolor='none')

        bnodes = kwargs.pop('bnodes', False)
        bsegs = kwargs.pop('bsegs', True)

        ax.axis(self.ax)
        fig.canvas.draw()

        xc = (self.ax[0] + self.ax[1])/2
        yc = (self.ax[2] + self.ax[3])/2

        # if isinstance(labels, bool):
        #     labels = ['s', 't', 'v', 'i', 'w']
        # elif isinstance(labels, str):
        #     labels = labels
        # else:
        #     labels = []

        k = list(self.Gs.pos.keys())
        v = list(self.Gs.pos.values())

        kk = np.array(k)
        vv = np.array(v)

        w = [str(x) for x in kk]

        #if 's' in labels:
        #    [ax.text(vv[i, 0], vv[i, 1], w[i]) for i in range(len(w))]

        if bnodes:
            point = ax.scatter([xc], [yc], color=[0,0,1], alpha=1)
            point.set_offsets(vv)
            ax.draw_artist(point)

            #ax.scatter(vv[:, 0], vv[:, 1])

        if bsegs:
            ML = sh.MultiLineString(list(self._shseg.values()))
            lseg = list(self._shseg.keys())
            lines_wall = [ [(l.xy[0][0],l.xy[1][0]),(l.xy[0][1],l.xy[1][1])] for k,l in
                    enumerate(ML) if lseg[k] in self.name['WALL']]
            lines_air = [ [(l.xy[0][0],l.xy[1][0]),(l.xy[0][1],l.xy[1][1])] for k,l in
                    enumerate(ML) if lseg[k] in self.name['_AIR']]
            lcwall = LineCollection(lines_wall,colors='k',linewidths=2)
            lcair = LineCollection(lines_air,colors='k',linewidths=1,linestyle='dotted')
            ax.add_collection(lcwall)
            ax.add_collection(lcair)

        plt.show()
        return fig, ax


    def get_Sg_pos(self, sigarr):
        """ return position of the signatures

        Parameters
        ----------

        sigarr : signature

        See Also
        --------

        showSig

        """

        signature = sigarr[0]
        sposfull = np.zeros((len(signature), 2))
        iz = np.nonzero(signature != 0)[0]
        spos = np.array([self.Gs.pos[i] for i in signature if i != 0])
        sposfull[iz, :] = spos

        return (sposfull)

    def plot_segments(self, lns, **kwargs):
        """"

        Parameters
        ----------
        lns
        *kwargs

        """
        defaults = {'show': False,
                    'fig': None,
                    'ax': None,
                    'color': 'b',
                    'linewidth': 1}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        elif kwargs['ax'] is None:
            ax = kwargs['fig'].add_subplot(111)
        else:
            fig = kwargs['fig']
            ax = kwargs['ax']

        # v1.1 nth = np.array(map(lambda n: nx.neighbors(self.Gs, n), lns))
        nth = np.array(map(lambda n: self.Gs[n], lns))
        nt = nth[:, 0]
        nh = nth[:, 1]
        # pt : 2 x Ns
        pt = np.array(map(lambda n:
                          [self.Gs.pos[n][0], self.Gs.pos[n][1]], nt)).T
        # ph : 2 x Ns
        ph = np.array(map(lambda n:
                          [self.Gs.pos[n][0], self.Gs.pos[n][1]], nh)).T

        fig, ax = plu.displot(pt, ph, fig=fig, ax=ax, color=kwargs['color'])

        return fig, ax

    def showSig(self, sigarr, Tx=None, Rx=None, fig=[], ax=None):
        """ Show signature

        Parameters
        ----------
        Tx  : np.array (2,1)
                Transmitter coordinates
        Rx  : np.array (2,1)
                Receipter coordinates
        sr  : boolean
                show room signature

        Returns
        -------
        fig   : figure instance
        ax    : axes instance
        lines : lines instance

        Examples
        --------

        """
        sig = sigarr[0]
        if fig == []:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        elif ax is None:
            ax = fig.add_subplot(111)
        lines = []
        ps = self.get_Sg_pos(sigarr)
        nz = np.nonzero(sig == 0)[0]
        mask = np.zeros((2, len(sig)))
        mask[:, nz] = 1
        vertices = np.ma.masked_array(ps.T, mask)
        lines.extend(ax.plot(vertices[0, :], vertices[1, :], color='k'))

        if Tx != []:
            itx = np.unique(sig[nz[1:-1] + 1], return_index=True)[1]
            itx2 = np.kron(itx, [1, 1])
            tx = ps[itx2]
            tx[range(0, len(tx), 2)] = Tx
            lines.extend(ax.plot(tx[:, 0], tx[:, 1], color='r'))
        if Rx != []:
            irx = np.unique(sig[nz[1:-1] - 1], return_index=True)[1]
            irx2 = np.kron(irx, [1, 1])
            rx = ps[irx2]
            rx[range(0, len(rx), 2)] = Rx
            lines.extend(ax.plot(rx[:, 0], rx[:, 1], color='b'))

        return (fig, ax, lines)
#        lines=[]
#        for s in sig:
#            l=[self.Gs.pos[s[ii]] for ii in xrange(len(s))]
#            if Tx!=None and Rx!=None:
#                l.insert(0,Tx)
#                l.insert(-1,Rx)
#            ls=sh.LineString(l)
#            x,y=ls.xy
#            lines.extend(ax.plot(x,y,'k',lw=0.1,alpha=0.2))
#        return (fig,ax,lines)

#    def distwall(self, p, nroom):
#        """ calculate distance to wall
#
#        Parameters
#        ----------
#
#        p : ndarray
#            point coordinate
#
#        nroom : int
#            room number of p
#
#        Returns
#        -------
#
#        dist
#                list of distances to walls of room nroom
#
#        Notes
#        -----
#
#        Return  dist a list of all the distances to the walls of a room
#
#
#        """
#        pp = sh.Point(p[0], p[1])
#
#        dist = []
#        p0_xy = []
#        p1_xy = []
#
#        vnode = self.Gr.node[nroom]['cycle'].cycle
#
#        # for j in range(len(Gr[nroom]['vnodes'])):
#        for j in range(len(vnodes)):
#            nn = self.b_Gr[5]['vnodes'][j]
#            nta = G1.tahe[0, nn - 1]
#            nhe = G1.tahe[1, nn - 1]
#            p0 = np.array([G1.pt[0, nta], G1.pt[1, nta]])
#            p1 = np.array([G1.pt[0, nhe], G1.pt[1, nhe]])
#            p0_xy.insert(j, p0)
#            p1_xy.insert(j, p1)
#
#        pstartwll = np.array(p0_xy)
#        pfinwll = np.array(p1_xy)
#
#        for i in range(len(self.b_Gr[nroom]['vnodes'])):
#            line_wall = sh.LineString([(pstartwll[i, 0],
#                                        pstartwll[i, 1]), (pfinwll[i, 0], pfinwll[i, 1])])
#            dist.insert(i, line_wall.distance(pp))
#        return(dist)

    def randTxRx(self):
        """returns random coordinates for Tx and Rx.

        Returns
        -------

        p_Tx : numpy.ndarray
             A point of the placement of the Tx
        p_Rx : numpy.ndarray
             A point of the placement of the Rx

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('defstr.lay')
        >>> p_Tx,p_Rx = L.randTxRx()

        Notes
        -----

        ex fn Tx_Rx_pos


        """

        Tx_x = rd.uniform(self.ax[0], self.ax[1])
        Tx_y = rd.uniform(self.ax[2], self.ax[3])
        Rx_x = rd.uniform(self.ax[0], self.ax[1])
        Rx_y = rd.uniform(self.ax[2], self.ax[3])

        p_Tx = np.array([Tx_x, Tx_y])
        p_Rx = np.array([Rx_x, Rx_y])

        return(p_Tx, p_Rx)

    def get_boundary(self):
        """ get and update Layout ax

        """
        xmax = max(p[0] for p in self.Gs.pos.values())
        xmin = min(p[0] for p in self.Gs.pos.values())
        ymax = max(p[1] for p in self.Gs.pos.values())
        ymin = min(p[1] for p in self.Gs.pos.values())
        self.ax = (xmin,xmax,ymin,ymax)

    def extend_boundary(self,pt,delta=1):
        """ extend boundary with a point

        Parameters
        ----------

        pt : np.array (,2)
        delta : offset distance

        """


        # determine which segment to delete
        # closest boundary segment to pt
        dmin = 1e15
        for ns in self.segboundary:
            ps = self.Gs.pos[ns]
            v = pt[0:2]-np.array(ps)
            dps_pt = np.linalg.norm(v)
            if dps_pt < dmin:
                dmin=dps_pt
                nsd = ns
                vn = v/dps_pt

        #
        n1,n2 = self.Gs.nodes[nsd]['connect']
        lcyn1 = self.Gs.nodes[n1]['ncycles']
        lcyn2 = self.Gs.nodes[n2]['ncycles']
        cyc = np.intersect1d(np.array(lcyn1),np.array(lcyn2))[0]
        polygon = self.Gt.nodes[cyc]['polyg']
        self.del_segment(nsd)
        self.segboundary.remove(nsd)
        # Create new layout point near pt
        npt = self.add_fnod((pt[0]+vn[0]*delta,pt[1]+vn[1]*delta))
        self.Gs.nodes[npt]['ncycles']=[cyc]
        # Create 2 new segment n1->npt and n2->npt
        ns1 = self.add_segment(n1, npt, name='_AIR')
        ns2 = self.add_segment(n2, npt, name='_AIR')
        self.segboundary.extend([ns1,ns2])
        self.lboundary.extend([npt])
        # Create a triangle polygon 
        p1 = np.array(self.Gs.pos[n1])
        p2 = np.array(self.Gs.pos[npt])
        p3 = np.array(self.Gs.pos[n2])
        pts = np.vstack((p1,p2,p3,p1)).T
        tri = geu.Polygon(pts)
        polysum = polygon+tri
        polysum.setvnodes(self)
         #tri.setvnodes(self)
        #pdb.set_trace()
        self.Gt.nodes[cyc]['polyg']=polysum
        self.Gt.pos[cyc]=(polysum.centroid.xy[0][0],polysum.centroid.xy[1][0])
        self.get_boundary()

    def boundary(self, **kwargs) :
        """ add a blank boundary around layout

        Parameters
        ----------

        percx : float
            percentage of Dx for x offset calculation (default 0.15)
        percy : float
           percentage of Dy for y offset calculation (default 0.15)
        xlim : tuple
        minD : minimum distance for boundary
        force : boolean
            force modification of boundaries even if one boundary already
            exists
        minD : int
            minimal distance over x and y

        self.lboundary is the list of the nodes of the added boundary
        self.axn is the zone without the boundary extension
        self.ax  is updated

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('defstr.lay')
        >>> L.boundary()

        Notes
        -----

        This function calls g2npy

        """

        percx = kwargs.pop('percx',0.15)
        percy = kwargs.pop('percy',0.15)
        xlim = kwargs.pop('xlim',())
        force = kwargs.pop('force', False)
        minD = kwargs.pop('minD', 10)
        bg2npy = kwargs.pop('bg2npy', True)

        if not self.hasboundary or force:

            if xlim != ():
                xmin = xlim[0]
                xmax = xlim[1]
                ymin = xlim[2]
                ymax = xlim[3]
            elif len(self.Gs.pos.values()) != 0:
                xmax = max(p[0] for p in self.Gs.pos.values())
                xmin = min(p[0] for p in self.Gs.pos.values())
                ymax = max(p[1] for p in self.Gs.pos.values())
                ymin = min(p[1] for p in self.Gs.pos.values())
            else:
                xmin = -20.
                xmax = 20.
                ymin = -10.
                ymax = 10.

            Dx = np.maximum(xmax - xmin, minD)
            Dy = np.maximum(ymax - ymin, minD)
            dx = Dx * percx
            dy = Dy * percy

            n1 = self.add_fnod((xmin - dx, ymin - dy))
            n2 = self.add_fnod((xmax + dx, ymin - dy))
            n3 = self.add_fnod((xmax + dx, ymax + dy))
            n4 = self.add_fnod((xmin - dx, ymax + dy))

            self.lboundary = [n1, n2, n3, n4]

            self.segboundary = []

            ns1 = self.add_segment(n1, n2, name='_AIR')
            ns2 = self.add_segment(n2, n3, name='_AIR')
            ns3 = self.add_segment(n3, n4, name='_AIR')
            ns4 = self.add_segment(n4, n1, name='_AIR')

            self.segboundary.append(ns1)
            self.segboundary.append(ns2)
            self.segboundary.append(ns3)
            self.segboundary.append(ns4)

            self.axn = (xmin, xmax, ymin, ymax)
            self.ax = (xmin - dx, xmax + dx, ymin - dy, ymax + dy)
            self.display['box'] = self.ax
            self.hasboundary = True

        elif xlim!=():

            # change points coordinates
            self.Gs.pos[self.lboundary[0]] = (xlim[0], xlim[2])
            self.Gs.pos[self.lboundary[1]] = (xlim[1], xlim[2])
            self.Gs.pos[self.lboundary[2]] = (xlim[1], xlim[3])
            self.Gs.pos[self.lboundary[3]] = (xlim[0], xlim[3])
            self.ax = xlim
            self.display['box'] = xlim

        if bg2npy:
            self.g2npy()


    def off_overlay(self, dx=0, dy=0):
        """ offset overlay image

        Parameters
        ----------

        dx : float
        dy : float

        """
        axis = (self.ax[0] + dx, self.ax[1] + dx,
                self.ax[2] + dy, self.ax[3] + dy)
        self.display['overlay_axis'] = axis

    def scl_overlay(self, ax=1.0, ay=1.0):
        """ scale overlay image

        Parameters
        ----------

        ax : float
        ay : float

        """

        axis = (self.ax[0] * ax, self.ax[1] * ax,
                self.ax[2] * ay, self.ax[3] * ay)
        self.display['overlay_axis'] = axis

    def get_paths(self, nd_in, nd_fin):
        """ returns the possible paths of graph Gs between two nodes.

        Parameters
        ----------
            nd_in: int
                initial graph node (segment or point)
            nd_fin: int
                final graph node (segment or point)

        Returns
        -------
            paths : list
                paths between nd_in and nd_fin
        """

        paths = gph.find_all_paths(self.Gs, nd_in, nd_fin)
        return paths


def outputGi_func_test(args):
    for k in range(10000):
        y = k*k+k*k
    return y

def outputGi_func(args):
# def outputGi_func(e, Gi_no, Gi_A, Gspos, sgsg, s2pc, s2pu):


    # for k in range(10000):
    #     y = k*k
    #     # time.sleep(0.01)
    # return y

    def Gspos(n):
        if n>0:
            #return np.mean(s2pc[n].reshape(2,2),axis=0)
            return np.mean(s2pc[n].toarray().reshape(2,2),axis=0)
        else:
            return p2pc[-n]

    e = args[0]
    #Gi_no = args[1]
    #Gi_A = args[2]
    #p2pc = args[3]
    #sgsg = args[4]
    #s2pc = args[5]
    #s2pu = args[6]

    print(e)



    i0 = e[0]
    i1 = e[1]
    nstr0 = i0[0]
    nstr1 = i1[0]

    # list of authorized outputs. Initialized void
    output = []

    # nstr1 : segment number of central interaction
    if nstr1 > 0:
        # central interaction is a segment
        # pseg1 = self.s2pc[nstr1,:].toarray().reshape(2, 2).T
        pseg1 = s2pc[nstr1,:].toarray().reshape(2, 2).T
        # pseg1 = self.s2pc[nstr1,:].data.reshape(2, 2).T
        # pseg1o = self.seg2pts(nstr1).reshape(2, 2).T

        # create a Cone object
        cn = cone.Cone()
        # if starting from segment
        if nstr0 > 0:
            # pseg0 = self.s2pc[nstr0,:].toarray().reshape(2, 2).T
            pseg0 = s2pc[nstr0,:].toarray().reshape(2, 2).T
            # pseg0 = self.s2pc[nstr0,:].data.reshape(2, 2).T
            # pseg0o = self.seg2pts(nstr0).reshape(2, 2).T

            # if nstr0 and nstr1 are connected segments
            if sgsg[nstr0,nstr1] == 0:
                # from 2 not connected segment
                cn.from2segs(pseg0, pseg1)
            else:
                # from 2 connected segments
                cn.from2csegs(pseg0, pseg1)
        # if starting from a point
        else:
            pt = Gspos(nstr0)
            cn.fromptseg(pt, pseg1)

        # list all potential successors of interaction i1
        ui2 = Gi_no.index(i1)
        ui = np.where(Gi_A[ui2,:]!=0)[0]
        i2 = [Gi_no[u] for u in ui]
        # i2 = nx.neighbors(self.Gi, i1)

        # how to find neighbors without network
        # ngi=L.Gi.nodes()
        # A=nx.adjacency_matrix(L.Gi)
        # inter = ngi[10]
        # u = ngi.index(inter)
        # ui = A[u,:].indices
        # neigh_inter = np.array([ngi[u] for u in ui])


        ipoints = [x for x in i2 if len(x)==1 ]
        #ipoints = filter(lambda x: len(x) == 1, i2)
        pipoints = np.array([Gspos(ip[0]) for ip in ipoints]).T
        # filter tuple (R | T)
        #istup = filter(lambda x : type(eval(x))==tuple,i2)
        # map first argument segment number
        #isegments = np.unique(map(lambda x : eval(x)[0],istup))
        # isegments = np.unique(
        #     filter(lambda y: y > 0, map(lambda x: x[0], i2)))
        isegments = np.unique([x[0] for x in i2 if x[0]>0])
        
        # if nstr0 and nstr1 are adjescent segment remove nstr0 from
        # potential next interaction
        # Fix 01/2017
        # This is not always True if the angle between 
        # the two adjascent segments is < pi/2
        # nb_nstr0 = self.Gs.neighbors(nstr0)
        # nb_nstr1 = self.Gs.neighbors(nstr1)
        # nb_nstr0 = np.array([self.s2pu[nstr0,0],self.s2pu[nstr0,1]])
        # nb_nstr1 = np.array([self.s2pu[nstr1,0],self.s2pu[nstr1,1]])
        nb_nstr0 = s2pu[nstr0,:].toarray()[0]
        nb_nstr1 = s2pu[nstr1,:].toarray()[0]
        print('nb_nstr0',nb_nstr0)
        #nb_nstr0 = s2pu[nstr0,:]
        #nb_nstr1 = s2pu[nstr1,:]
        # common_point = np.intersect1d(nb_nstr0,nb_nstr1)
        common_point = np.array([x for x in nb_nstr0 if x in nb_nstr1])
        # if len(common_point) == 1:
        if common_point.any():
            num0 = [x for x in nb_nstr0 if x != common_point]
            num1 = [x for x in nb_nstr1 if x != common_point]
            p0 = Gspos(num0[0])
            p1 = Gspos(num1[0])
            pc = Gspos(common_point[0])
            v0 = p0-pc 
            v1 = p1-pc 
            v0n = v0/np.sqrt(np.sum(v0*v0))
            v1n = v1/np.sqrt(np.sum(v1*v1))
            if np.dot(v0n,v1n)<=0:
                isegments = np.array([ x for x in isegments if x != nstr0 ]) 
            #    filter(lambda x: x != nstr0, isegments))
        # there are one or more segments
        # if len(isegments) > 0:
        if isegments.any():

            li1 = len(i1)

            points = self.s2pc[isegments,:].toarray().T
            #points = s2pc[isegments,:].T
            # points = self.s2pc[isegments,:].data.reshape(4,len(isegments))
            # pointso = self.seg2pts(isegments)

            pta = points[0:2, :]
            phe = points[2:, :]
            # add difraction points
            # WARNING Diffraction points are added only if a segment is seen
            # it should be the case in 99% of cases

            if len(ipoints) > 0:
                isegments = np.hstack(
                    (isegments, np.array(ipoints)[:, 0]))
                pta = np.hstack((pta, pipoints))
                phe = np.hstack((phe, pipoints))

            # cn.show()

            # if i0 == (38,79) and i1 == (135,79,23):
            #     printi0,i1
            #     import ipdb
            #     ipdb.set_trace()
            # i1 : interaction T
            if li1 == 3:
                typ, prob = cn.belong_seg(pta, phe)
                # if bs.any():
                #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                # if ~bs.any():
                #    plu.displot(pta[:,~bs],phe[:,~bs],color='k')

            # i1 : interaction R --> mirror
            elif li1 == 2:
                Mpta = geu.mirror(pta, pseg1[:, 0], pseg1[:, 1])
                Mphe = geu.mirror(phe, pseg1[:, 0], pseg1[:, 1])
                typ, prob = cn.belong_seg(Mpta, Mphe)
                # printi0,i1
                # if ((i0 == (6, 0)) & (i1 == (7, 0))):
                #    pdb.set_trace()
                # if bs.any():
                #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                # if ~bs.any():
                #    plu.displot(pta[:,~bs],phe[:,~bs],color='m')
                #    plt.show()
                #    pdb.set_trace())
            ########
            # SOMETIMES PROBA IS 0 WHEREAS SEG IS SEEN
            ###########
            # # keep segment with prob above a threshold
            # isegkeep = isegments[prob>0]
            # # dict   {numint : proba}
            # dsegprob = {k:v for k,v in zip(isegkeep,prob[prob>0])}
            # 4 lines are replaced by
            # keep segment with prob above a threshold
            utypseg = typ != 0
            isegkeep = isegments[utypseg]
            # dict   {numint : proba}
            dsegprob = {k: v for k, v in zip(isegkeep, prob[utypseg])}
            #########
            # output = filter(lambda x: x[0] in isegkeep, i2)
            output = [x for x in i2 if x[0] in isegkeep]
            # probint = map(lambda x: dsegprob[x[0]], output)
            probint = [dsegprob[x[0]] for x in output]
            # dict interaction : proba
            dintprob = {k: v for k, v in zip(output, probint)}

            # keep all segment above nstr1 and in Cone if T
            # keep all segment below nstr1 and in Cone if R

    else:
        # central interaction is a point

        # 1) Simple approach
        #       output interaction are all visible interactions
        # 2) TO BE DONE
        #
        #       output of the diffraction points
        #       exploring
        # b
        #          + right of ISB
        #          + right of RSB
        #
        #  + using the wedge cone
        #  + using the incident cone
        #

        # output = nx.neighbors(self.Gi, (nstr1,))
        uout = Gi_no.index((nstr1,))
        ui = np.where(Gi_A[uout,:]!=0)[0]
        output = [Gi_no[u] for u in ui]

        nout = len(output)
        probint = np.ones(nout)  # temporarybns
        dintprob = {k: v for k, v in zip(output, probint)}

    return (i0,i1, {'output':dintprob})
    # self.Gi.add_edge(i0, i1, output=dintprob)


if __name__ == "__main__":
    plt.ion()
    doctest.testmod()
    # L = Layout('Servon Sur Vilaine',verbose=True,dist_m=60)
    # L.build()    
