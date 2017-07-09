# -*- coding: utf-8 -*-
#
#
from __future__ import print_function
try:
    from tvtk.api import tvtk
    from mayavi import mlab
except:
    print('Layout:Mayavi is not installed')
import pdb
import sys
import os
import logging
import copy
import glob
import time
import tqdm
import numpy as np
import numpy.random as rd
import scipy as sp
import scipy.sparse as sparse
import doctest
import triangle
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import networkx as nx
from itertools import combinations, product
import ast
from networkx.readwrite import write_gpickle, read_gpickle
from mpl_toolkits.basemap import Basemap
import shapely.geometry as sh
import shapefile as shp
from shapely.ops import cascaded_union
from descartes.patch import PolygonPatch
from numpy import array
import PIL.Image as Image
import logging
if sys.version_info.major==2:
    from  urllib2 import urlopen
    import ConfigParser
else:
    from  urllib.request import urlopen
    import configparser
import hashlib
#from cStringIO import StringIO

from pathos.multiprocessing import ProcessingPool as Pool
from pathos.multiprocessing import cpu_count


# from multiprocessing import Pool
from functools import partial

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
import pylayers.util.easygui as eag
import pylayers.util.project as pro

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

    Graphs
    ------
    Gs : structure
    Gt : topology 
    Gv : visibility
    Gi : interaction
    Gr : room  
    Gm :  
    Gw : ways 

    Integer
    -------
    Np
    Ns 
    Nss 

    Tuple
    -----
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

    array
    -----
    pt : points coordinates  
    tahe : segment tail head 
    tgs : graph to segment
    tsg : segment to graph 
    upnt : array of point index

    sparse array
    ------------

    s2pc : segment to point coordinates
    s2pu : segment to point index
    sgsg 

    Slabs
    -----
    sl 


    String
    ------
    typ  : 'floorplan' | 'outdoor'
    coordinates : 'cart','lonlat'
    version
    _filename 
    _hash

    Dictionnaries
    -------------

    _shseg : keys / segment index 
             values / shapely LineString
    dca    : keys / Gt node 
             values / list of air wall 
    degree : keys / point degree
             values / array of index 
    display : dictionnary for controling various visualization
    dsseg : 

    boolean 
    -------

    indoor : if True allow indoor penetration 
    isbuilt 
    diffraction 

    heights
    -------
    maxheight
    zceil 
    zfloor 
    zmin

    """

    def __init__(self, string='',
                 _filematini='matDB.ini',
                 _fileslabini='slabDB.ini',
                 _filefur='',
                 bcheck=True,
                 bbuild=False,
                 bgraphs=True,
                 bindoor=False,
                 bdiffraction=False,
                 bverbose=False,
                 bcartesian=True,
                 dist_m=400,
                 typ='floorplan'):
        """ object constructor

        Parameters
        ----------

        arg : string
            layout file name, address or '(lat,lon)'
        _filematini :
            material dB file name
        _fileslabini :
            slab dB file name
        _filefur :
            furniture file name
        force : booleanlo
        check : boolean
        build : boolean 
        verbose : boolean 
        cartesian : boolean 
        dist_m : int 
        typ : string 
            'floorplan' | 'outdoor'


        """

        # mat = sb.MatDB()
        # mat.load(_filematini)

        # self.sl = sb.SlabDB()
        # self.sl.mat = mat
        # self.sl.load(_fileslabini)

        self.labels = {}

        self.Np = 0
        self.Ns = 0
        self.Nss = 0

        #
        # Initializing graphs
        #

        self.Gs = nx.Graph()
        self.Gr = nx.Graph()
        self.Gt = nx.Graph()
        self.Gm = nx.Graph()

        self.Gs.pos = {}
        self.tahe = np.zeros(([2, 0]), dtype=int)
        self.lbltg = []

        self.Gt.pos = {}

        #
        # related files
        #

        self.fileslabini = _fileslabini
        self.filematini = _filematini
        self.filefur = _filefur

        self.hasboundary = False
        self.coordinates = 'cart'
        self.version = '1.1'
        self.typ = typ
        # boolean 

        self.isbuilt = False
        self.loadosm = False
        # diffraction : activate diffraction 
        self.diffraction = bdiffraction
        # indoor : activate indoor propagation 
        self.indoor = bindoor

        #
        # setting display option
        #

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
        self.display['subsegnb'] = True
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
        self.display['box'] = (-50, 50, -50, 50)
        self.name = {}
        self.ax = self.display['box']
        self.zmin = 0
        self.maxheight = 3

        newfile = False
        loadini = False
        loadosm = False
        loadres = False
       
        #
        # Layout main argument
        #   If no .ini extension provided it is added
        #
        arg, ext = os.path.splitext(string)
        if arg != '':
            if ext == '.ini':
                self._filename = string
                loadini = True
            elif ext == '.osm':
                self._filename = arg + '.ini'
                loadosm = True
            elif ext == '.res':
                self._filename = arg + '.ini'
                loadres = True
            else:
                self.typ = 'outdoor'
        else:  # No argument
            self._filename = 'newfile.ini'
            newfile = True

        if not newfile:
            if loadini:
                filename = pyu.getlong(self._filename, pro.pstruc['DIRINI'])
                if os.path.exists(filename):  # which exists
                    self.load()
                else:  # which do not exist
                    newfile = True
                    print("new file - creating a void Layout", self._filename)
            elif loadosm:  # load .osm file
                self.importosm(_fileosm=string, cart=True,typ=self.typ)
                self.loadosm = True
            elif loadres:
                self.importres(_fileres=string)
                self.sl = sb.SlabDB()
            elif '(' in string:  # load from osmapi latlon in string
                self.importosm(latlon=string, dist_m=dist_m, cart=True,typ=self.typ)
                self.loadosm = True
            else:  # load from address geocoding
                self.importosm(address=string, dist_m=dist_m, cart=True,typ=self.typ)
                self.loadosm = True
            
            # add boundary if it not exist
            if not self.hasboundary:    
                self.boundary()
            self.subseg()
            self.updateshseg()
            try:
                self.geomfile()
            except:
                print("problem to construct geomfile")

            #
            # check layout 
            #
            if bcheck:
                self.check()

            # check if the graph gpickle files have been built
            
            if bgraphs:
                dirname = self._filename.replace('.ini','')
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
                    if self._hash != self.Gt.node[0]['hash']:
                        rebuild = True
                    else:  # reload
                        self.dumpr('stvirw')
                        self.isbuilt = True
                        bbuild = False

                else:
                    print("graphs have not been saved")
                    bbuild = True

            # build and save graphs 
            if bbuild:
                # ans = raw_input('Do you want to build the layout (y/N) ? ')
                # if ans.lower()=='y'
                self.build()
                self.lbltg.append('s')
                self.dumpw()

    def __repr__(self):
        st = '\n'
        st = st + "----------------\n"
        home = os.path.expanduser('~')
        with open(os.path.join(home,'.pylayers'),'r') as f:
            paths = f.readlines()
        uporj = paths.index('project\n')
        project = paths[uporj+1]
        st = st + "Project : " + project+'\n'
        if hasattr(self,'_hash'):
            st = st + self._filename + ' : ' + self._hash + "\n"
        else:
            st = st + self._filename + "\n"
        
        if self.isbuilt:
            st = st + 'Built with : ' + self.Gt.node[0]['hash'] + "\n"
        st = st + 'Type : '+ self.typ+'\n'
        if self.indoor:
            st = st + 'Indoor : Activated'+'\n'
        else:
            st = st + 'Indoor : Not activated'+'\n'

        if self.diffraction:
            st = st + 'Diffraction : Activated'+'\n'
        else:
            st = st + 'Diffraction : Not Activated'+'\n'
        if self.display['overlay_file'] != '':
            filename = pyu.getlong(
                self.display['overlay_file'], os.path.join('struc', 'images'))
            st = st + "Image('" + filename + "')\n"
        st = st + "Coordinates : " + self.coordinates + "\n"
        st = st + "----------------\n"
        if hasattr(self,'Gs'):
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
        if hasattr(self, 'degree'):
            for k in self.degree:
                if (k < 2) or (k > 3):
                    st = st + 'degree ' + \
                        str(k) + ' : ' + str(self.degree[k]) + "\n"
                else:
                    st = st + 'number of node point of degree ' + \
                        str(k) + ' : ' + str(len(self.degree[k])) + "\n"
        st = st + "\n"
        st = st + "xrange :" + str(self.ax[0:2]) + "\n"
        st = st + "yrange :" + str(self.ax[2:]) + "\n"
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
        st = st + "Segment s in Gs => s_ab coordinates \n"
        st = st + "s2pc : segment to point coordinates (sparse) [p1,p2] = L.s2pc.toarray().reshape(2,2).T \n"
        st = st + \
            "s -> u = self.tgs[s] -> v = self.tahe[:,u] -> s_ab = self.pt[:,v]\n\n"
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
            st = st + "s2pc : "+"from a Gs segment node to its 2 extremal points (tahe) index\n"
        if hasattr(self,'p2pu'):
            st = st + "p2pc : "+"from a Gs point node to its coordinates\n"
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

    def ls(self, typ='ini'):
        """ list the available file in dirstruc

        Parameters
        ----------

        typ : string optional
            {'ini'|'osm'|'wrl'}

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

        if typ == 'ini':
            pathname = os.path.join(pro.pstruc['DIRINI'], '*.' + typ)
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

        """

        newpoint = dict((k - offp, v)
                        for k, v in self.Gs.node.items() if k < 0)
        assert (np.array(newpoint.keys()) < 0).all()
        newseg = dict((k + offs, v) for k, v in self.Gs.node.items() if k > 0)
        assert (np.array(newseg.keys()) > 0).all()
        newpoint.update(newseg)
        self.Gs.node = newpoint

        newppoint = dict((k - offp, v)
                         for k, v in self.Gs.pos.items() if k < 0)
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
        self.Gs.adj = dseg
        self.Gs.edge = dseg

    def check(self, level=0):
        """ Check Layout consistency


        Parameters
        ----------

        level : int

        Returns
        -------

        consistent : Boolean
              True if consistent

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
        consistent = True

        nodes = self.Gs.nodes()
        if len(nodes) > 0:
            #
            # points
            # segments
            # degree of segments
            useg = filter(lambda x: x > 0, nodes)
            upnt = filter(lambda x: x < 0, nodes)
            degseg = map(lambda x: nx.degree(self.Gs, x), useg)

            #
            # 1)   all segments have degree 2
            #
            assert(np.all(array(degseg) == 2))

            #
            # degree of points
            # maximum degree of points
            #

            degpnt = map(lambda x: nx.degree(self.Gs, x),
                         upnt)  # points absolute degrees
            degmin = min(degpnt)
            degmax = max(degpnt)

            #
            #  No isolated points (degree 0)
            #  No points of degree 1
            #
            if (degmin <= 1):
                f, a = self.showG('s', aw=1)
                deg0 = filter(lambda x: nx.degree(self.Gs, x) == 0, upnt)
                deg1 = filter(lambda x: nx.degree(self.Gs, x) == 1, upnt)

                if len(deg0) > 0:
                    print("It exists degree 0 points :  %r" % deg0)
                    f, a = self.pltvnodes(deg0, fig=f, ax=a)
                if len(deg1) > 0:
                    print("It exists degree 1 points : %r" % deg1)
                    f, a = self.pltvnodes(deg1, fig=f, ax=a)

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
            ke = self.Gs.pos.keys()
            x = np.array(map(lambda x: x[0], self.Gs.pos.values()))
            y = np.array(map(lambda x: x[1], self.Gs.pos.values()))
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
            # Not very much scalable. Double for loop
            for s in useg:
                n1, n2 = np.array(self.Gs.neighbors(s))  # node s neighbors
                p1 = np.array(self.Gs.pos[n1])           # p1 --- p2
                p2 = np.array(self.Gs.pos[n2])  # s
                #
                # iterate on upnt : list of points
                for n in upnt:
                    if (n1 != n) & (n2 != n):
                        p = np.array(self.Gs.pos[n])
                        if geu.isBetween(p1, p2, p):
                            print(n1, p1)
                            print(n2, p2)
                            print(n, p)

                            logging.critical(
                                "segment %d contains point %d", s, n)
                            consistent = False
                if level > 0:
                    cycle = self.Gs.node[s]['ncycles']
                    if len(cycle) == 0:
                        logging.critical("segment %d has no cycle", s)
                    if len(cycle) == 3:
                        logging.critical(
                            "segment %d has cycle %s", s, str(cycle))
        #
        # check if Gs points are unique
        # segments can be duplicated
        #
        P = np.array([self.Gs.pos[k] for k in upnt])
        similar = geu.check_point_unicity(P)
        if len(similar) != 0:
            logging.critical(
                "points at index(es) %s in self.Gs.pos are similar", str(similar))
            consistent = False

        return(consistent)

    def clip(self, xmin, xmax, ymin, ymax):
        """ return the list of edges which cross or belong to the clipping zone

        Parameters
        ----------

        xmin : float
        xmax : float
        ymin : float
        ymax : float

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

        Notes
        -----

        This function updates the following arrays:

        self.pt   (2xNp)
        self.tahe (2xNs)
        self.tgs
        self.dca  : dictionnary of cycle with an airwall
        self.lsss : list of subsegments

        assert self.pt[self.iupnt[-1]] == self.pt[:,self.iupnt[-1]]

        """

        nodes = self.Gs.nodes()
        # nodes include points and segments

        # segment index
        # useg = filter(lambda x: x > 0, nodes)
        useg = [n for n in nodes if n >0]

        # points index
        # upnt = filter(lambda x: x < 0, nodes)
        upnt = [n for n in nodes if n < 0]


        # matrix segment-segment 
        # usage 
        # self.sgsg[seg1,seg2] => return common point
        mno = max(self.Gs.nodes())
        self.sgsg = sparse.lil_matrix((mno+1,mno+1),dtype='int')

        for s in useg:

            lpts = self.Gs.edge[s].keys()
            a = self.Gs.edge[lpts[0]].keys()
            b = self.Gs.edge[lpts[1]].keys()

            nsa = np.setdiff1d(a,b)
            nsb = np.setdiff1d(b,a)
            u = np.hstack((nsa,nsb))

            npta = [lpts[0]]*len(nsa)
            nptb = [lpts[1]]*len(nsb)
            ns = np.hstack((npta,nptb))

            self.sgsg[s,u]=ns




        # conversion in numpy array
        self.upnt = np.array((upnt))

        # association
        # pdb.set_trace()

        # utmp = np.array(zip(-self.upnt,np.arange(len(self.upnt))))
        # mutmp = max(utmp[:,0])
        # self.iupnt = -np.ones((mutmp+1),dtype='int')
        # self.iupnt[utmp[:,0]]=utmp[:,1]

        # degree of segment nodes
        degseg = map(lambda x: nx.degree(self.Gs, x), useg)
        
        assert(np.all(array(degseg) == 2))  # all segments must have degree 2

        #
        # self.degree : dictionnary (point degree : list of point index)
        #

        # points absolute degrees
        degpnt = np.array(map(lambda x: nx.degree(self.Gs, x), upnt))

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
            lseg = nx.neighbors(self.Gs, nupt)
            n = 0
            for ns in lseg:
                if ns in lairwall:
                    n = n + 1
            return n

        nairwall = np.array(map(nairwall, upnt))
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
        if verbose:
            print('Start node degree determination')
        for deg in range(degmax + 1):
            num = filter(lambda x: degpnt[x] == deg, range(
                len(degpnt)))  # position of degree 1 point
            # number of degree 1 points
            npt = np.array(map(lambda x: upnt[x], num))
            self.degree[deg] = npt

        if verbose:
            print('Node degree determination  : Done')
        #
        # convert geometric information in numpy array
        #

        self.pt = np.array(np.zeros([2, len(upnt)]), dtype=float)
        self.tahe = np.array(np.zeros([2, len(useg)]), dtype=int)

        self.Np = len(upnt)
        self.Ns = len(useg)

        self.pt[0, :] = np.array([self.Gs.pos[k][0] for k in upnt])
        self.pt[1, :] = np.array([self.Gs.pos[k][1] for k in upnt])

        if verbose:
            print('pt in np.array  : Done')

        self.pg = np.sum(self.pt, axis=1) / np.shape(self.pt)[1]
        self.pg = np.hstack((self.pg, 0.))

        # ntail = map(lambda x: nx.neighbors(self.Gs, x)[0], useg)
        # nhead = map(lambda x: nx.neighbors(self.Gs, x)[1], useg)
        ntahe = np.array([nx.neighbors(self.Gs, x) for x in useg])

        ntail = ntahe[:,0]
        nhead = ntahe[:,1]

        # create sparse matrix from a Gs segment node to its 2 extremal points (tahe) index
        mlgsn = max(self.Gs.nodes())+1
        self.s2pu = sparse.lil_matrix((mlgsn,2),dtype='int')
        self.s2pu[useg,:] = ntahe
        # convert to compressed row sparse matrix 
        # to be more efficient on row slicing
        self.s2pu = self.s2pu.tocsr()
        

        # tic = time.time()
        # self.tahe[0, :] = np.array(
        #      map(lambda x: np.nonzero(np.array(upnt) == x)[0][0], ntail))
        # self.tahe[1, :] = np.array(
        #    map(lambda x: np.nonzero(np.array(upnt) == x)[0][0], nhead))
        
        aupnt = np.array(upnt)
        self.tahe[0, :] = np.array([np.where(aupnt==x)[0][0] for x in ntail ])
        self.tahe[1, :] = np.array([np.where(aupnt==x)[0][0] for x in nhead ])
        
        if verbose:
            print('tahe in numpy array : Done')
        #
        # transcoding array between graph numbering (discontinuous) and numpy numbering (continuous)
        #
        Nsmax = 0
        self.tsg = np.array(useg)
        try:
            Nsmax = max(self.tsg)
        except:
            logging.warning("No segments in Layout yet")

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
            self.lsss = [x for x in useg if len(self.Gs.node[x]['iso']) > 0]

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
            #     #    pdb.set_trace()
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
        #pdb.set_trace()
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
        assert(len(lheight)>0),logging.error("no valid heights for segments")
        self.maxheight = np.max(lheight)
        # self.maxheight=3.
        # calculate extremum of segments
        self.extrseg()

    def importshp(self, **kwargs):
        """ import layout from shape file

        Parameters
        ----------

        _fileshp :

        """
        defaults = {'pref': [np.array([25481100, 6676890]), np.array([60.2043716, 24.6591147])],
                    'dist_m': 250,
                    'latlon': True,
                    'bd': [24, 60, 25, 61],
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

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
            v = np.array(p) - kwargs['pref'][0][None, :]
            nv = np.sqrt(np.sum(v * v, axis=1))
            # if at least one point is in the radius the poygon is kept
            if (nv < kwargs['dist_m']).any():
                # pdb.set_trace()
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

        if kwargs['latlon']:
            lat_ref = kwargs['pref'][1][0]
            lon_ref = kwargs['pref'][1][1]
            x_ref, y_ref = self.m(lon_ref, lat_ref)
            Dx = kwargs['pref'][0][0] - x_ref
            Dy = kwargs['pref'][0][1] - y_ref
            pos = np.array(self.Gs.pos.values())
            for k, keys in enumerate(self.Gs.pos.keys()):
                self.Gs.pos[keys] = self.m(
                    pos[k, 0] - Dx, pos[k, 1] - Dy, inverse=True)

            self.coordinates = 'latlon'

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

        """
        fileres = pyu.getlong(_fileres, os.path.join('struc', 'res'))
        D  = np.fromfile(fileres,dtype='int',sep=' ')
        # number of integer
        N1 = len(D)
        # number of lines
        N2 = N1/8
        D = D.reshape(N2,8)
        # list of coordinates
        lcoords = []
        # list of ring
        lring = [] 
        # list of (z_ground, height_building)
        zring = [] 
        # 
        bdg_old = 1
        for e in range(N2):
            # p1 point coordinate
            p1 = ([D[e,0],D[e,1]])
            # p2 point coordinate
            p2 = ([D[e,2],D[e,3]])
            # (ground height,building height) 
            z  = (D[e,7]-500,D[e,4])
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
        
        for r1,z1 in zip(lring,zring):
            x,y = r1.xy 
            
            for k2 in range(len(x)):
                new_pt = (x[k2],y[k2])
                kpos = self.Gs.pos.keys()
                vpos = self.Gs.pos.values()
                if new_pt not in vpos:
                    current_node_index = -npt
                    self.Gs.add_node(current_node_index)
                    self.Gs.pos[-npt] = new_pt
                    npt = npt + 1
                else:
                    u = [k for k in range(len(vpos)) if (vpos[k] == new_pt)]
                    
                    current_node_index = kpos[u[0]]

                if k2>0: # at least already one point
                    ns = self.add_segment(current_node_index, previous_node_index, name='WALL', z=z1)
                else:
                    starting_node_index  =   current_node_index
                previous_node_index = current_node_index
            # last segment    
            #ns = self.add_segment(previous_node_index, starting_node_index, name='WALL', z=z1)
        #pdb.set_trace()

    def importosm(self, **kwargs):
        """ import layout from osm file or osmapi

        Parameters
        ----------

        _fileosm : string
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
        which is very convenient for editing floorplan of buildings.. 

        In josm editor, nodes are numbered with negative indexes, while in 
        pylayers they have a positive index.

        See Also
        --------

        pylayers.gis.osmparser.osmparse

        """
        defaults = {'_fileosm': '',
                    'address': 'Rennes',
                    'typ': 'floorplan',
                    'latlon': '0',
                    'dist_m': 200,
                    'cart': False
                    }

        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        typ = kwargs['typ']
        address = kwargs['address']
        latlon = eval(kwargs['latlon'])
        dist_m = kwargs['dist_m']
        cart = kwargs['cart']
        
        #
        # TODO : Not clean get zceil from actual data
        #

        if self.typ=='floorplan':
            self.zceil = 3
            self.zfloor = 0
        
        if kwargs['_fileosm'] == '':  # by using osmapi address or latlon
            coords, nodes, ways, dpoly, m = osm.getosm(typ=typ,
                                                       address=address,
                                                       latlon=latlon,
                                                       dist_m=dist_m,
                                                       cart=cart)
            if cart:
                self.coordinates='cart'
            else:
                self.coordinates='latlon'
            if kwargs['latlon'] == '0':
                self._filename = kwargs['address'].replace(' ', '_') + '.ini'
            else:
                lat, lon = eval(kwargs['latlon'])
                self._filename = 'lat_' + \
                    str(lat).replace('.', '_') + '_lon_' + \
                    str(lon).replace('.', '_') + '.ini'
        else:  # by reading an osm file
            fileosm = pyu.getlong(kwargs['_fileosm'], os.path.join('struc', 'osm'))
            coords, nodes, ways, relations, m = osm.osmparse(fileosm, typ=typ)
            self.coordinates = 'latlon'
            self._filename = kwargs['_fileosm'].replace('osm', 'ini')
        
        # 2 valid typ : 'floorplan' and 'building'

        _np = 0  # _ to avoid name conflict with numpy alias
        _ns = 0
        ns = 0
        nss = 0

        # Reading points  (<0 index)

        # Reorganize points coordinates for detecting
        # duplicate nodes
        # duplicate nodes are saved in dict dup

        kp = [k for k in coords.xy]

        x = np.array(map(lambda x: coords.xy[x][0], kp))
        y = np.array(map(lambda x: coords.xy[x][1], kp))
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
        for k, nseg in enumerate(ways.way):
            tahe = ways.way[nseg].refs
            for l in range(len(tahe) - 1):
                nta = tahe[l]
                nhe = tahe[l + 1]
                #
                # if a node is duplicate recover the original node
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


                # getting segment information
                if 'name' in d:
                        slab = d['name']
                else:  # the default slab name is WALL
                        slab = "WALL"
                if 'z' in d:
                    z = d['z']
                else:
                    z = (0, 3)
                if 'offset' in d:
                    offset = d['offset']
                else:
                    offset = 0
                #
                # get the common neighbor of nta and nhe if it exists
                #
                u1 = np.array(nx.neighbors(self.Gs, nta))
                u2 = np.array(nx.neighbors(self.Gs, nhe))
                inter_u1_u2 = np.intersect1d(u1, u2)
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

        # pdb.set_trace()
        # self.m = Basemap(llcrnrlon=bd[0], llcrnrlat=bd[1],
        #                  urcrnrlon=bd[2], urcrnrlat=bd[3],
        #                  resolution='i', projection='cass', lon_0=lon_0, lat_0=lat_0)
        

        self.m = m
        if ((kwargs['cart']) and (self.coordinates!='cart')):
             x, y = self.m(lon, lat)
             self.Gs.pos = {k: (x[i], y[i]) for i, k in enumerate(self.Gs.pos)}
             self.coordinates = 'cart'

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
        mat.load(self.filematini)
        self.sl = sb.SlabDB()
        self.sl.mat = mat
        self.sl.load(self.fileslabini)

        #
        # update self.name with existing slabs database entries
        #
        for k in self.sl.keys():
            if k not in self.name:
                self.name[k] = []

        # convert graph Gs to numpy arrays for speed up post processing
        #pdb.set_trace()
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
                    if self.coordinates == 'latlon':
                        lon, lat = self.Gs.pos[n]
                    if self.coordinates == 'cart':
                        x, y = self.Gs.pos[n]
                        lon, lat = self.m(x, y, inverse=True)
                    fd.write("<node id='" + str(n) + "' action='modify' visible='true' lat='" +
                             str(lat) + "' lon='" + str(lon) + "' />\n")

        for n in self.Gs.pos:
            if n > 0:
                #
                # Conditions pour ajout segments
                # 
                cond1 = not ((not self.indoor)       and 
                         (self.Gs.node[n]['name']=='AIR')  and
                        (self.Gs.node[n][z][1]>2000)) 
                cond2 = (self.Gs.node[n]['name'] != '_AIR')
                if (cond1 and cond2):
                    neigh = nx.neighbors(self.Gs, n)
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

    def save(self):
        """ save structure in an ini file

        """
        print(self._filename)
        current_version = 1.2
        config = ConfigParser.RawConfigParser()
        config.optionxform = str
        config.add_section("info")
        config.add_section("points")
        config.add_section("segments")
        config.add_section("display")
        config.add_section("files")
        config.add_section("slabs")
        config.add_section("materials")
        if self.coordinates == 'latlon':
            config.set("info", "format", "latlon")
        else:
            config.set("info", "format", "cart")
        config.set("info", "version", current_version)
        config.set("info", "type", self.typ)

        if self.typ == 'floorplan':
            config.add_section("floorplan")
            config.set("floorplan", "zceil", self.zceil)
            config.set("floorplan", "zfloor", self.zfloor)

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

        # config.set("info",'Npoints',self.Np)
        # config.set("info",'Nsegments',self.Ns)
        # config.set("info",'Nsubsegments',self.Nss)

        for k in self.display:
            config.set("display", k, self.display[k])

        # iterate on points
        # boundary nodes and air walls are not saved
        for n in self.Gs.pos:
            if n < 0:
                if n not in self.lboundary:
                    config.set("points", str(
                        n), (self.Gs.pos[n][0], self.Gs.pos[n][1]))

        # iterate on segments
        for n in self.Gs.pos:
            if n > 0:
                if self.Gs.node[n]['name'] != '_AIR':
                    d = self.Gs.node[n]
                    # old format conversion
                    if 'ncycles' in d:
                        del d['ncycles']
                    if 'ss_ce1' in d:
                        del d['ss_ce1']
                    if 'ss_ce2' in d:
                        del d['ss_ce2']
                    if 'zmin' in d:
                        d['z'] = [d['zmin'], d['zmax']]
                        del(d['zmin'])
                        del(d['zmax'])
                    if 'ss_zmin' in d:
                        d['ss_z'] = [[d['ss_zmin'], d['ss_zmax']]]
                        d['ss_name'] = [d['ss_name']]
                        del(d['ss_zmin'])
                        del(d['ss_zmax'])

                    d['connect'] = nx.neighbors(self.Gs, n)
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
                    config.set("segments", str(n), d)

        # list of used slab
        lslab = [x for x in self.name if len(self.name[x]) > 0]
        lmat = []
        #
        # In case an osm file has been read; there is no .sl
        # By default all the available slab and materials are provided
        if not hasattr(self, 'sl'):
            self.sl = sb.SlabDB(filemat='matDB.ini', fileslab='slabDB.ini')

        for s in lslab:
            ds = {}
            if s not in self.sl:
                if s not in self.sl.mat:
                    self.sl.mat.add(name=s,cval=6,sigma=0,typ='epsr')
                self.sl.add(s,[s],[0.1])

            ds['index'] = self.sl[s]['index']
            ds['color'] = self.sl[s]['color']
            ds['lmatname'] = self.sl[s]['lmatname']
            for m in ds['lmatname']:
                if m not in lmat:
                    lmat.append(m)
            ds['lthick'] = self.sl[s]['lthick']
            ds['linewidth'] = self.sl[s]['linewidth']
            config.set("slabs", s, ds)

        if "_AIR" not in lslab:
            air = {'color': 'white', 'index': 1, 'linewidth': 1,
                   'lthick': [0.1], 'lmatname': ['AIR']}
            config.set("slabs", "_AIR", air)

        if "AIR" not in lslab:
            air = {'color': 'white', 'index': 1, 'linewidth': 1,
                   'lthick': [0.1], 'lmatname': ['AIR']}
            config.set("slabs", "AIR", air)

        if "CEIL" not in lslab:
            ceil = {'color': 'grey20', 'index': 6, 'linewidth': 1,
                    'lthick': [0.1], 'lmatname': ['REINFORCED_CONCRETE']}
            config.set("slabs", "CEIL", ceil)

        if "FLOOR" not in lslab:
            floor = {'color': 'grey40', 'index': 7, 'linewidth': 1,
                     'lthick': [0.1], 'lmatname': ['REINFORCED_CONCRETE']}
            config.set("slabs", "FLOOR", floor)

        for m in lmat:
            dm = self.sl.mat[m]
            config.set("materials", m, dm)

        if "REINFORCED_CONCRETE" not in lmat:
            reic = {'index': 6, 'name': 'REINFORCED_CONCRETE', 'mur': (
                1 + 0j), 'epr': (8.69999980927 + 0j), 'roughness': 0.0, 'sigma': 3.0}
            config.set("materials", "REINFORCED_CONCRETE", reic)
        # config.set("files",'materials',self.filematini)
        # config.set("files",'slab',self.fileslabini)
        config.set("files", 'furniture', self.filefur)
        fileini = pyu.getlong(self._filename, pro.pstruc['DIRINI'])
        fd = open(fileini, "w")
        config.write(fd)
        fd.close()
        # convert graph Gs to numpy arrays for speed up post processing
        # ideally an edited Layout should be locked while not saved.
        # self.g2npy()
        self._hash = hashlib.md5(open(fileini, 'rb').read()).hexdigest()

    def load(self):
        """ load a structure file from an .ini file

        The filename is self._filename

        """

        # di : dictionnary which reflects the content of ini file
        di = {}
        config = ConfigParser.RawConfigParser()
        config.optionxform = str
        fileini = pyu.getlong(self._filename, pro.pstruc['DIRINI'])
        config.read(fileini)
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
        self.Gs = nx.Graph()
        self.Gs.pos = {}
        self.labels = {}

        #
        # Check file version
        #
        if 'version' in di['info']:
            self.version = di['info']['version']
            self.name = {}
        else:
            self.version = 0.9
            mat = sb.MatDB()
            mat.load(self.filematini)
            self.sl = sb.SlabDB()
            self.sl.mat = mat
            self.sl.load(self.fileslabini)
            for k in self.sl.keys():
                self.name[k] = []

        if di['info'].has_key('type'):
            self.typ = di['info']['type']
            if self.typ == 'floorplan':
                self.zceil = eval(di['floorplan']['zceil'])
                self.zfloor = eval(di['floorplan']['zfloor'])
            if self.typ == 'outdoor':

                self.zceil = eval(di['outdoor']['zceil'])
                self.zfloor = eval(di['outdoor']['zfloor'])
        else:
            self.typ = 'floorplan'
            self.zfloor = 0
            self.zceil = 3
        # manage ini file with latlon coordinates
        #
        # if the format is latlon, coordinates are converted into
        # cartesian coordinates with the coords.cartesian method
        #

        if di['info'].has_key('format'):
            if di['info']['format'] == 'latlon':
                or_coord_format = 'latlon'
                coords = osm.Coords()
                coords.clean()
                coords.latlon = {i: np.array(
                    eval(di['points'][i])) for i in di['points']}
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
        for k in di['display']:
            try:
                self.display[k] = eval(di['display'][k])
            except:
                self.display[k] = di['display'][k]

        self.ax = self.display['box']

        #
        # POINTS
        #
        # update points section
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
            self.Gs.pos[nodeindex] = (
                round(1000 * x) / 1000., round(1000 * y) / 1000.)
            self.labels[nodeindex] = nn

        #
        # SEGMENTS
        #
        # update segments section

        self.name['AIR'] = []
        self.name['_AIR'] = []
        #
        # get the maximum index
        #
        maxnum = max([eval(x) for x in di['segments'].keys()])
        for k, key in enumerate(di['segments']):

            d = eval(di['segments'][key])
            nta = d['connect'][0]
            nhe = d['connect'][1]
            #print(key,nta,nhe)
            if not d.has_key('offset'):
                offset = 0
            else:
                offset = d['offset']
            
            name = d['name']
            z = d['z']
            num = self.add_segment(nta, nhe,
                                   num = eval(key), 
                                   name=name,
                                   offset=offset,
                                   z=z)

            # exploit iso for segment completion 
            #
            #  Complement single segment which do not reach zceil or zfloor with
            #  an iso segment with AIR property
            # 
            # if di['info']['type'] == 'outdoor':
            if z[1] < self.zceil:
                 num = self.add_segment(nta, nhe,
                                        name='AIR',
                                        maxnum = maxnum, 
                                        offset=offset,
                                        z=(z[1], self.zceil))

            if z[0] > self.zfloor:
                 num = self.add_segment(nta, nhe,
                                        name='AIR',
                                        maxnum = maxnum, 
                                        offset=offset,
                                        z=(self.zfloor,z[0]))

       
        self.boundary()
        
        # compliant with config file without  material/slab information
        if config.has_section('latlon'):
            llcrnrlon = eval(config.get('latlon', 'llcrnrlon'))
            llcrnrlat = eval(config.get('latlon', 'llcrnrlat'))
            urcrnrlon = eval(config.get('latlon', 'urcrnrlon'))
            urcrnrlat = eval(config.get('latlon', 'urcrnrlat'))
            lon_0 = (llcrnrlon+urcrnrlon)/2.
            lat_0 = (llcrnrlat+urcrnrlat)/2.

            # Construction of Basemap for coordinates transformation
            self.m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                    urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                resolution='i', projection='cass', lon_0=lon_0, lat_0=lat_0)
            

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

        if self.display.has_key('fileoverlay'):
            self.display['overlay_file'] = self.display.pop('fileoverlay')
            self.display['overlay_axis'] = self.display['box']
            self.save()

        if self.display.has_key('inverse'):
            self.display['overlay_flip'] = ""
            self.display.pop('inverse')
            self.save()

        # convert graph Gs to numpy arrays for faster post processing
        self.g2npy()
        #
        self._hash = hashlib.md5(open(fileini, 'rb').read()).hexdigest()
        

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
            >>> L = Layout('WHERE1.ini')
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
        +  .ini   : ini file format (natural one) DIRINI


        """

        newfile = False
        filename = pyu.getlong(_filename, pro.pstruc['DIRINI'])
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

        for k in self.Gs.node.keys():
            dk = self.Gs.node[k]
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

        Example
        -------

            >>> from pylayers.gis.layout import *
            >>> L = Layout('defstr.str')
            >>> L.add_fnod((10.0,10.0))
            -9


        """
        try:
            num = -(max(-np.array(self.Gs.node.keys())) + 1)
        except:
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
        np1 = self.Gs.neighbors(s1)
        np2 = self.Gs.neighbors(s2)
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
        # printxA,yA
        # printxB,yB
        # printxC,yC
        # printxD,yD
        # printxP,yP
        A = np.array([[xB - xA, xD - xC], [yB - yA, yD - yC]])
        b = np.array([xP - xA, yP - yA])
        x = sp.linalg.solve(A, b)
        if ((x[0] > 0.) & (x[0] < 1.0)):
            self.add_pons(e1, 1 - x[0])
        # printx

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
        nop = self.Gs.neighbors(ns)
        namens = self.Gs.node[ns]['name']
        zminns = self.Gs.node[ns]['z'][0]
        zmaxns = self.Gs.node[ns]['z'][1]
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

    def add_segment(self, 
                    n1,
                    n2,
                    num=-1,
                    maxnum=-1,
                    name='PARTITION', 
                    z=(0.0, 40000000), 
                    offset=0,
                    verbose=True):
        """  add segment between node n1 and node n2

        Parameters
        ----------

        n1  : integer < 0
        n2  : integer < 0
        name : string
            layer name 'PARTITION'
        z : list of float
            default = (0,40000000)
        offset : float
            [-1,1] default (0)

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

        # if 2 points are selected

        if ((n1 < 0) & (n2 < 0) & (n1 != n2)):
            nseg = [s for s in self.Gs.node if s > 0]
            if num==-1:
                if len(nseg) > 0:
                    num = max(maxnum+1,max(nseg) + 1)   # index not given 
                else: # first segment index not given
                    num = 1
            else:
                pass # segment index given  
        else:
            if verbose:
                print("add_segment : error not a node", n1, n2)
            return

        transition = False
        if (name == '_AIR'):
            # if name == 'AIR':
            transition = True

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
        # is there an other segments with the same neighbors ?

        nbnta = self.Gs.neighbors(n1)
        nbnhe = self.Gs.neighbors(n2)
        same_seg = list(set(nbnta).intersection(nbnhe))

        #
        # Impossible to have duplicated _AIR
        #
        if (name == '_AIR'):
            if len(same_seg) > 0:
                return None

        self.Gs.add_node(num, name=name,
                         z=z,
                         norm=norm,
                         transition=transition,
                         offset=offset,
                         connect=[n1, n2],
                         iso=[],
                         ncycles=[]
                         )

        for k in same_seg:
            self.Gs.node[k]['iso'].append(num)
            self.Gs.node[num]['iso'].append(k)

        #
        # Segment position in the middle
        #
        self.Gs.pos[num] = tuple((p1 + p2) / 2.)

        #
        # Connectivity
        #
        self.Gs.add_edge(n1, num)
        self.Gs.add_edge(n2, num)

        #
        # Update current total number of segments
        #
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
        return(num)

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

        aseg = map(lambda x: filter(lambda y: y not in self.name['AIR'],
                                    nx.neighbors(self.Gs, x)),
                   apnt)
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

        aseg = map(lambda x: filter(lambda y: y not in
                                    self.name['AIR'],
                                    nx.neighbors(self.Gs, x)),
                   lpnt)

        pts = np.array(map(lambda x: self.seg2pts(
            [x[0], x[1]]).reshape(4, 2), aseg))
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
        self.add_segment(n0, n1, matname, [zmin, zmin + height])
        self.add_segment(n1, n2, matname, [zmin, zmin + height])
        self.add_segment(n2, n3, matname, [zmin, zmin + height])
        self.add_segment(n3, n0, matname, [zmin, zmin + height])

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

        print("lp : ", lp)
        # get segments involved in points list
        ls = self.nd2seg(lp)

        print("ls : ", ls)
        # 1) delete involved segments
        for k in ls:
            assert(k > 0)
            self.del_segment(k)
            print('del ', k)
        # 2) delete involved points
        for n1 in lp:
            assert(n1 < 0)
            nbrs = self.Gs.neighbors(n1)
            self.Gs.remove_node(n1)
            del self.Gs.pos[n1]
            self.labels.pop(n1)
            self.Np = self.Np - 1
        # 3) updating structures
        self.g2npy()

    def del_segment(self, le, verbose=True, g2npy=True):
        """ delete segment e

        Parameters
        ----------

        le : list of segment number

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
            self.del_subseg(e, verbose=verbose)
            name = self.Gs.node[e]['name']
            del self.Gs.pos[e]  # delete edge position
            self.Gs.remove_node(e)
            self.labels.pop(e)
            self.Ns = self.Ns - 1
            # update slab name <-> edge number dictionnary
            self.name[name].remove(e)
            # delete subseg if required
            try:
                self.pop._shseg(e)
            except:
                pass
        if g2npy:
            self.g2npy()

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

        """

        tseg = []

        for k in self.Gs.node.keys():
            if k > 0:
                lnp = self.Gs.neighbors(k)
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

        """

        for n in self.Gs.node.keys():
            if ((n < 0) & (self.Gs.degree(n) == 0)):
                self.Gs.remove_node(n)
                del self.Gs.pos[n]
                try:
                    self.Gv.remove_node(n)
                except:
                    pass

        self.Np = len(np.nonzero(np.array(self.Gs.node.keys()) < 0)[0])
        self.g2npy()

    def displaygui(self):
        """ open a GUI for displaying configuration

        """

        displaygui = eag.multenterbox('', 'Display Parameters',
                                  ('filename',
                                   'nodes',
                                   'ednodes',
                                   'ndlabel',
                                   'edlabel',
                                   'edges',
                                   'subseg',
                                   'visu',
                                   'thin',
                                   'scaled',
                                   'overlay',
                                   'overlay_file',
                                   'overlay_flip',
                                   'alpha'),
                                  (self._filename,
                                   int(self.display['nodes']),
                                   int(self.display['ednodes']),
                                   int(self.display['ndlabel']),
                                   int(self.display['edlabel']),
                                   int(self.display['edges']),
                                   int(self.display['subseg']),
                                   int(self.display['visu']),
                                   int(self.display['thin']),
                                   int(self.display['scaled']),
                                   int(self.display['overlay']),
                                   self.display['overlay_file'],
                                   self.display['overlay_flip'],
                                   self.display['alpha']))
        if displaygui is not None:
            self._filename = displaygui[0]
            self.display['nodes'] = bool(eval(displaygui[1]))
            self.display['ednodes'] = bool(eval(displaygui[2]))
            self.display['ndlabel'] = bool(eval(displaygui[3]))
            self.display['edlabel'] = bool(eval(displaygui[4]))
            self.display['edges'] = bool(eval(displaygui[5]))
            self.display['subseg'] = bool(eval(displaygui[6]))
            self.display['visu'] = bool(eval(displaygui[7]))
            self.display['thin'] = bool(eval(displaygui[8]))
            self.display['scaled'] = bool(eval(displaygui[9]))
            self.display['overlay'] = bool(eval(displaygui[10]))
            self.display['overlay_file'] = displaygui[11]
            self.display['overlay_flip'] = eval(displaygui[12])
            self.display['alpha'] = eval(displaygui[14])

    def info_segment(self, s1):
        """ information about segment

        Parameters
        ----------

        s1 : segment number

        """
        nebd = self.Gs.neighbors(s1)
        n1 = nebd[0]
        n2 = nebd[1]
        nns1 = self.Gs.neighbors(n1)
        nns2 = self.Gs.neighbors(n2)
        ds1 = self.Gs.node[s1]
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

    def edit_point(self, np):
        """ edit point

        Parameters
        ----------

        np : integer
            point number

        """
        title = "Point (" + str(np) + ")"
        message = "Enter coordinates "
        pt = self.Gs.pos[np]
        data = eag.multenterbox(message, title, (('x', 'y')),
                            ((str(pt[0]), str(pt[1]))))
        self.Gs.pos[np] = tuple(eval(data[0]), eval(data[1]))

    def chgmss(self, ns, ss_name=[], ss_z=[], ss_offset=[], g2npy=True):
        """ change multi subsegments properties

        Parameters
        ----------

        ns : int
            segment number

        ss_name : list of Nss  string
            name of the different constitutive SLAB of the multi-segments

        ss_z : list of Nss tuple (zmin,zmax)

        ss_offset : list of subsegment offsets

        Examples
        --------

        See Also
        --------

        pylayers.gis.layout.g2npy

        """

        #
        # Check slabname in ss_name
        #
        #  update self.sl with new slab values

        SDB = sb.SlabDB()
        for sname in ss_name:
            if sname not in self.sl.keys():
                if sname not in SDB.keys():
                    print('')
                else:
                    slab = SDB[sname]
                    self.sl[slab['name']] = slab

        if ss_z != []:
            assert len(ss_name) == len(
                ss_z), 'Error incompatible size in chgmss'
        if ss_offset != []:
            assert len(ss_z) == len(
                ss_offset), 'Error incompatible size in chgmss'

        if ns in self.Gs.node.keys():
            if self.Gs.node[ns].has_key('ss_name'):
                if ss_name != []:
                    self.Gs.node[ns]['ss_name'] = ss_name
                if ss_z != []:
                    self.Gs.node[ns]['ss_z'] = ss_z
                if ss_offset != []:
                    self.Gs.node[ns]['ss_offset'] = ss_offset
                else:
                    self.Gs.node[ns]['ss_offset'] = [0] * len(ss_name)

                # update Layout information
                if g2npy:
                    self.g2npy()

    def edit_segment(self, e1, gui=True, outdata={}):
        """ edit segment WITH EasyGUI

        Parameters
        ----------

        e1 : integer
            edge number
        gui : boolean

        Notes
        -----

        A segment has the following properties :
            + name  : string
            + z  :  tuple
            + transition : boolean (default FALSE)

        If a segment has subsegments attached the following properties are
        added :
            + ss_name : string
            + ss_z : subsegment [(min height (meters),max height (meters))]


        """
        nebd = self.Gs.neighbors(e1)
        n1 = nebd[0]
        n2 = nebd[1]
        de1 = self.Gs.node[e1]
        title = "Segment (" + str(n1) + ',' + str(n2) + ")"
        message = str(self.sl.keys())
        if 'ss_name' not in de1.keys():
            de1k = ['name', 'z', 'transition', 'offset']
            de1v = [de1['name'], de1['z'], de1['transition'], de1['offset']]
        else:
            de1k = ['name', 'z', 'ss_name', 'ss_z', 'transition', 'ss_offset']
            de1v = [de1['name'], de1['z'], de1['ss_name'], de1['ss_z'],
                    de1['transition'], de1['ss_offset']]
        #de1v    = de1.values()
        if gui:
            outdata = {}
            data0 = choicebox('chose slab', title, self.sl.keys())
            try:
                data1 = eag.multenterbox(
                    'attribute for ' + data0, title, tuple(de1k[1:]), tuple(de1v[1:]))
                d1 = data1[0].split(' ')
                d1t = tuple((eval(d1[0]), eval(d1[1])))
                data1[0] = d1t
                data = [data0] + data1
                #data = eag.multenterbox(message, title, tuple(de1k), tuple(de1v))
                i = 0
                self.name[de1['name']].remove(e1)
                for k in de1k:
                    try:
                        self.Gs.node[e1][k] = eval(data[i])
                        outdata[k] = eval(data[i])
                    except:
                        self.Gs.node[e1][k] = data[i]
                        outdata[k] = data[i]
                        if k == 'name':
                            try:
                                self.name[data[i]].append(e1)
                            except:
                                self.name[data[i]] = [e1]
                    i = i + 1
            except:
                # if cancel
                pass
        else:
            if outdata == {}:
                pass
                # data = {}
                # val = '1'
                # while(val!='0'):
                #     clear
                #     print'0 : exit'
                #     for e,(k,v) in enumerate(zip(de1k,de1v)):
                #         printstr(e+1)+ ' '+k+': '+  str(v)+'\n'
                #     val = input('Your choice :')
                #     if val!='0':
                #         pass
            else:
                for k in de1k:
                    if k in ['z', 'name', 'transition', 'offset']:
                        self.Gs.node[e1][k] = outdata[k]
        return outdata

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

            ename = self.Gs.node[e1]['name']
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
                self.Gs.node[e1][k] = data[k]

        self.name[data['name']].append(e1)

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
        dk = self.Gs.node[e1]
        if 'ss_name' in dk:
            return True
        else:
            return False

    def del_subseg(self, e1, verbose=False):
        """ delete sub subsegent

        Parameters
        ----------
        e1 : integer
             segment number
        """
        assert (e1 > 0)
        if self.have_subseg(e1):
            self.Gs.node[e1].pop('ss_name')
            self.Gs.node[e1].pop('ss_z')
            self.Gs.node[e1].pop('ss_offset')
            try:
                self.Gs.node[e1].pop('ss_ce')
            except:
                pass
            self.Gs.node[e1].pop('transition')
            self.Nss -= 1
        elif verbose:
            print("no subseg to delete")

    def add_subseg(self, s1, name='DOOR', zmin=0, zmax=2.24, offset=0, transition=True):
        """ add a subsegment on a segment WITH EasyGUI

        Parameters
        ----------

        s1 : integer
            edge number > 0
        name : string
            slab name
        zmin : float
            default 0
        zmax : float
            default 2.24 m

        """
        self.info_segment(s1)
        message = str(self.sl.keys())
        title = 'Add a subsegment'
        data = eag.multenterbox(message, title, ('name', 'zmin', 'zmax', 'offset'),
                            (name, zmin, zmax, offset))
        try:
            self.Gs.node[s1]['ss_name'] = [data[0]]
            self.Nss += 1
            self.chgmss(s1, ss_name=[data[0]], ss_offset=[
                        eval(data[3])], ss_z=[(eval(data[1]), eval(data[2]))])
            # self.Gs.node[s1]['ss_z'] =
            # ce is for Pulsray compatibility
            #self.Gs.node[s1]['ss_ce'] = [ (0,0) ]
            self.Gs.node[s1]['transition'] = transition
            return True
        except:
            return False

    def update_sseg(self, s1, data={}, g2npy=True):
        """ update subsegment(s) on a segment

        Parameters
        ----------

        s1 : integer
            edge number > 0
        data = dict
            dictionnary of data

        """

        assert len(data['ss_name']) == len(
            data['ss_z']), 'Error incompatible size in chgmss'
        assert len(data['ss_z']) == len(data['ss_offset']
                                        ), 'Error incompatible size in chgmss'
        if s1 < 0:
            return False

        new_nbss = len(data['ss_name'])

        try:
            old_nbss = len(self.Gs.node[s1]['ss_name'])
        except:
            old_nbss = 0
        # update the number of subsegments for self.Nss
        deltaNss = new_nbss - old_nbss
        # printdeltaNss
        if new_nbss != 0:
            self.Gs.node[s1]['ss_name'] = [data['ss_name']]
            self.Nss += deltaNss
            self.chgmss(s1, ss_name=data['ss_name'],
                        ss_offset=data['ss_offset'],
                        ss_z=data['ss_z'], g2npy=g2npy)
            return True
        else:
            if self.Gs.node[s1].has_key('ss_name'):
                self.Gs.node[s1].pop('ss_name')
                self.Gs.node[s1].pop('ss_z')
                self.Gs.node[s1].pop('ss_offset')
                self.Nss += deltaNss
                if g2npy:
                    self.g2npy()
                return True
            else:
                return True

    def add_window(self, s1, z):
        """ add a window on segment

        Parameters
        ----------

        s1 : integer
            segment number
        z : tuple of float
            (zmin,zmax)

        """
        if (zmin > self.Gs.node[e1]['z'][0]) & (zmax < self.Gs.node[e1]['z'][1]):
            self.info_edge(s1)
            self.Gs.node[s1]['ss_name'].append('WINDOW')
            self.Gs.node[s1]['ss_z'].append((zmin, zmax))
            self.Gs.node[s1]['ss_ce'].append((0, 0))
            self.Gs.node[s1]['transition'] = False
            self.Nss += 1
        else:
            logging.warning('windows range is wrong')

    def add_door(self, s1, zmin, zmax):
        """ add a door on segment

        Parameters
        ----------

        s1 : integer
            segment number
        zmin : float
        zmax : float


        """
        if (zmin > self.Gs.node[s1]['z'][0]) & (zmax < self.Gs.node[s1]['z'][1]):
            self.info_segment(s1)
            self.Gs.node[s1]['ss_name'].append('DOOR')
            self.Gs.node[s1]['ss_zmin'].append((zmin, zmax))
            self.Gs.node[s1]['ss_ce'].append((0, 0))
            self.Gs.node[s1]['transition'] = True
            self.Nss += 1

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
        l = map(lambda x: self.Gs.adj[x].keys(), ndlist)
        seglist = reduce(lambda x, y: x + y, l)

        return(np.unique(seglist))

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
        for n in self.Gs.node.keys():
            if n < 0:
                x = self.Gs.pos[n][0]
                y = self.Gs.pos[n][1]
                if ((x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)):
                    ptlist.append(n)

        seglist = self.nd2seg(ptlist)

        return ptlist, seglist

    def get_points(self, ax):
        """ get points list and segments list in a rectangular zone

        Parameters
        ----------

        ax  : list ot tuple
            [xmin,xmax,ymin,ymax]
              or shapely Polygon 

        Returns
        -------

        (pt,ke)

        """


        if type(ax) == geu.Polygon:
            eax = ax.exterior.xy
            xmin = np.min(eax[0])
            xmax = np.max(eax[0])
            ymin = np.min(eax[1])
            ymax = np.max(eax[1])
        else:
            xmin = ax[0]
            xmax = ax[1]
            ymin = ax[2]
            ymax = ax[3]

        x = self.pt[0,:]
        y = self.pt[1,:]
        uxmin = (x>=xmin)
        uymin = (y>=ymin)
        uxmax = (x<=xmax)
        uymax = (y<=ymax)
        k  = np.where(uxmin*uymin*uxmax*uymax==1)[0]
        pt = np.array(zip(x[k],y[k])).T
        ke = self.upnt[k]
        # ux = ((x>=xmin).all() and (x<=xmax).all())
        # uy = ((y>=ymin).all() and (y<=ymax).all())
        return((pt,ke))

    def angleonlink3(self, p1=np.array([0, 0, 1]), p2=np.array([10, 3, 1])):
        """ angleonlink(self,p1,p2) return (seglist,angle) between p1 and p2

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
        >>> L = Layout('DLR.ini')
        >>> p1 = np.array([0,0,1])
        >>> p2 = np.array([10,3,2])
        >>> data = L.angleonlink3(p1,p2)

        #array([(0, 141, 1.2793395519256592), (0, 62, 0.29145678877830505),
               (0, 65, 0.29145678877830505)],
              dtype=[('i', '<i8'), ('s', '<i8'), ('a', '<f4')])

        See Also
        --------

        antprop.loss.Losst 
        
        """

        sh1 = np.shape(p1)
        sh2 = np.shape(p2)

        assert sh1[0] == 3
        assert sh2[0] == 3

        if (len(sh1) < 2) & (len(sh2) > 1):
            p1 = np.outer(p1, np.ones(sh2[1]))

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
        #
        
        seglist  = np.unique(self.seginframe2(p1[0:2], p2[0:2]))
    

        upos = np.nonzero(seglist >= 0)[0]
        uneg = np.nonzero(seglist < 0)[0]
        

        # nNLOS = len(uneg) + 1
        # # retrieve the number of segments per link
        # if nNLOS > 1:
        #     llink = np.hstack(
        #         (uneg[0], np.hstack((uneg[1:], array([len(seglist)]))) - uneg - 1))
        # else:
        #     llink = np.array([len(seglist)])
        # [(link id,number of seg),...]
        # nl = zip(np.arange(nlink),llink)n

        npta = self.tahe[0, seglist[upos]]
        nphe = self.tahe[1, seglist[upos]]

        Pta = self.pt[:, npta]
        Phe = self.pt[:, nphe]

        
        # #
        # # This part should possibly be improved
        # #

        # for i, nl in enumerate(llink):
        #     try:
        #         # P1 = np.hstack((P1,np.outer(p1[:,i],np.ones(nl))))
        #         # P2 = np.hstack((P2,np.outer(p2[:,i],np.ones(nl))))
        #         ilink = np.hstack(
        #             (ilink, array([-1]), i * np.ones(nl, dtype='int')))
        #     except:
        #         # P1 = np.outer(p1[:,i],np.ones(nl))
        #         # P2 = np.outer(p2[:,i],np.ones(nl))
        #         ilink = i * np.ones(nl, dtype='int')

        # check for intersection P1P2 PtaPhe
        # bo = geu.intersect(P1[0:-1], P2[0:-1], Pta, Phe)
        Nscreen = len(npta)
        # get segment height bounds
        zmin = np.array([self.Gs.node[x]['z'][0]
                         for x in self.tsg[seglist[upos]]])
        zmax = np.array([self.Gs.node[x]['z'][1]
                         for x in self.tsg[seglist[upos]]])
        # centroid of the screen
        Pg = np.vstack(((Phe + Pta) / 2., (zmax + zmin) / 2.))
        Ptahe = Phe - Pta
        L1 = np.sqrt(np.sum(Ptahe * Ptahe, axis=0))
        # 3 x Nscreen U1 is in plane xy
        U1 = np.vstack((Ptahe / L1, np.zeros(Nscreen)))
        L2 = zmax - zmin
        U2 = np.array([0, 0, 1])[:, None]  # 3 x 1  U2 is along z

        bo = geu.intersect3(p1, p2, Pg, U1, U2, L1, L2)

        ubo = np.where(bo)

        Nseg = len(ubo[0])
        data = np.zeros(Nseg, dtype=[
                        ('i', 'i8'), ('s', 'i8'), ('a', np.float32)])

        data['i']=ubo[0]
        data['s']=self.tsg[ubo[1]]
        
        #
        # Calculate angle of incidence refered from segment normal
        #

        norm = self.normal[:, ubo[1]]
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
        >>> L = Layout('DLR.ini')
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

        seglist = self.seginframe2(p1, p2)
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
        seglist = np.array(map(lambda x: self.tsg[x], seglist2))
        data = np.zeros(len(seglist), dtype=[
                        ('i', 'i8'), ('s', 'i8'), ('a', np.float32)])

        #
        # update subsegment in seglist
        #
        # self.lsss

        data['i'] = idxlnk
        data['s'] = seglist
        data['a'] = angle
        return(data)

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
        #>>> L = Layout('DLR.ini','matDB.ini','slabDB.ini')
        #>>> p1 = np.array([0,0])
        #>>> p2 = np.array([10,3])
        #>>> L.angleonlinkold(p1,p2)
        #(array([59, 62, 65]), array([ 1.27933953,  0.29145679,  0.29145679]))

        Notes
        -----


        """

        logging.warning('This function is deprecated use')

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
        seglist = np.array(map(lambda x: self.tsg[x], seglist))

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
        >>> L = Layout('DLR.ini')
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

        """ return the seg list of a sequence of point number

        Parameters
        ----------

        ptlist

            array(1xNp) Point number array

        Returns
        -------

        seglist
            array seglist associated with ptlist

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('TA-Office.ini')
        >>> ptlist  = np.array([0,1])
        >>> L.segpt(ptlist)
        array([44, 50, 84, 86])

        """

        #seglist = np.array([], dtype=int)
        ut = map(lambda x: np.hstack((np.nonzero(self.tahe[0, :] == x)[0],
                                      np.nonzero(self.tahe[1, :] == x)[0])), ptlist)
        utstack = reduce(lambda x, y: np.hstack((x, y)), ut)
        #uvstack = reduce(lambda x,y : np.hstack((x,y)),uv)
        # for i in ptlist:
        #    ut = np.nonzero(self.tahe[0, :] == i)[0]
        #    uv = np.nonzero(self.tahe[1, :] == i)[0]
        #    seglist = np.hstack((seglist, ut, uv))
        seglist = np.unique(utstack)

        return(seglist)

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
        >>> L = Layout('defstr.ini')
        >>> aseg = np.array([1,3,6])
        >>> pt =  L.seg2pts(aseg)

        OBSOLETE : Use self.s2pc instead

        """

        if not isinstance(aseg, np.ndarray):
            aseg = np.array([aseg])

        assert(len(np.where(aseg < 0)[0]) == 0)
        utahe = self.tgs[aseg]
        if (utahe>=0).all():
            tahe = self.tahe[:, utahe]
            ptail = self.pt[:, tahe[0, :]]
            phead = self.pt[:, tahe[1, :]]
            pth = np.vstack((ptail, phead))

            pth = pth.reshape(pth.shape[0], pth.shape[-1])
            return pth
        else:
            pdb.set_trace()



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
        >>> L = Layout('TA-Office.ini')
        >>> ptlist  = np.array([0,1])
        >>> seg = L.segpt(ptlist)

        Notes
        -----

        segpt is faster than segpt2

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
        th = zip(self.tahe[0, :], self.tahe[1, :])

        self.max_sx = np.array(
            map(lambda x: max(pt[0, x[0]], pt[0, x[1]]), th))
        self.min_sx = np.array(
            map(lambda x: min(pt[0, x[0]], pt[0, x[1]]), th))
        self.max_sy = np.array(
            map(lambda x: max(pt[1, x[0]], pt[1, x[1]]), th))
        self.min_sy = np.array(
            map(lambda x: min(pt[1, x[0]], pt[1, x[1]]), th))

    def seginframe2(self, p1, p2):
        """ returns the seg list of a given zone defined by two points

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


            Examples
            --------

            >>> from pylayers.gis.layout import *
            >>> L = Layout('TA-Office.ini')
            >>> p1 = np.array([[0,0,0],[0,0,0]])
            >>> p2 = np.array([[10,10,10],[10,10,10]])
            >>> seglist = L.seginframe2(p1,p2)
            >>> edlist  = map(lambda x: L.tsg[x],seglist)
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

        max_x = map(lambda x: max(x[1], x[0]), zip(p1[0, :], p2[0, :]))
        min_x = map(lambda x: min(x[1], x[0]), zip(p1[0, :], p2[0, :]))
        max_y = map(lambda x: max(x[1], x[0]), zip(p1[1, :], p2[1, :]))
        min_y = map(lambda x: min(x[1], x[0]), zip(p1[1, :], p2[1, :]))

        seglist = map(lambda x: np.nonzero((self.max_sx > x[0]) &
                                           (self.min_sx < x[1]) &
                                           (self.max_sy > x[2]) &
                                           (self.min_sy < x[3]))[0],
                      zip(min_x, max_x, min_y, max_y))

        # np.array stacking
        # -1 acts as a deliminiter (not as a segment number)

        seglist = reduce(lambda x, y: np.hstack((x, array([-1]), y)), seglist)

        return(seglist)

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
            >>> L = Layout('TA-Office.ini')
            >>> p1 = np.array([0,0])
            >>> p2 = np.array([10,10])
            >>> L.seginframe(p1,p2)
            array([ 1,  3,  7,  8, 14, 15, 16, 17, 18, 20, 21, 23, 24, 26, 27, 29, 30,
                   32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 46, 47, 52, 53, 54,
                   55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
                   72, 73, 74, 75, 76, 77, 78, 81, 82, 85, 86])

        """
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
            ta, he = self.Gs.neighbors(seg)
            pa = np.array(self.Gs.pos[ta])
            pb = np.array(self.Gs.pos[he])

            segline = sh.LineString((pa, pb))

            if line.intersects(segline):
                lc.extend(self.Gs.node[seg]['ncycles'])
            # printseg,self.Gs.node[seg]['ncycles']
                ls.append(seg)
                psh = line.intersection(segline)
                I = np.hstack((I, np.array([[psh.x], [psh.y]])))
        v = (I - p1[:, np.newaxis])
        dv = np.sum(v * v, axis=0)
        u = np.argsort(dv)
        lss = np.array(ls)[u]

        lc = [c1]
        for s in lss:
            cy1, cy2 = self.Gs.node[s]['ncycles']
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
                ta, he = self.Gs.neighbors(seg)
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
            ndlist = self.Gs.node.keys()
        # elif ndlist[0]==1e8:
        #    ndlist  = self.Gs.node.keys()

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
        ue = (np.ones(2 * len(kwargs['edlist']))).astype('int').tolist()
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
            fig, ax = self.show_nodes(
                ndlist=kwargs['edlist'], color='b', fig=fig, ax=ax)

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
                    color = slab.tocolor(fGHz)
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
            edinter = edges[kwargs['ninter']]
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
        #pdb.set_trace()
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

        # pdb.set_trace()
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
            tn = np.array(self.Gs.node.keys())
            u = np.nonzero(tn < 0)[0]
            ndlist = tn[u]

        if kwargs['edlist'] == []:
            tn = self.Gs.node.keys()
            #u  = np.nonzero(tn > 0)[0]
            #edlist = tn[u]
            edlist = filter(lambda x: (x > 0), tn)
            #& (not self.Gs.node[x].has_key('ss_name')),tn)
        else:
            edlist = kwargs['edlist']

        if self.display['nodes']:
            dlabels = self.display['ndlabel']
            fig, ax = self.show_nodes(
                ndlist, size=30, color='k', dlabels=dlabels, node_shape='s', fig=fig, ax=ax)

        # if self.display['subsegnb']:
        #     if hasattr(self,'lsss'):
        #         seg = self.lsss
        #         psseg = np.array([[self.Gs.pos[x][0],self.Gs.pos[x][1]] for x in seg])
        #         nbsseg = np.array([len(self.Gs.node[x]['ss_name']) for x in seg],dtype='int')

        #         [ax.text(psseg[x,0]+0.2,psseg[x,1]+0.2,str(nbsseg[x]),
        # fontdict={'size':8},ha='center') for x in range(len(seg))]

        if self.display['transition']:
            try:
                segwtrans = [y for y in [x for x in self.Gs.nodes() if x > 0]if self.Gs.node[
                    y]['transition']]
                posseg = np.array([self.Gs.pos[x] for x in segwtrans])
                normseg = np.array([self.Gs.node[x]['norm']
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
        if self.display['ticksoff']:
            ax.xaxis.set_ticks([])
            for loc, spine in ax.spines.iteritems():
                spine.set_color('none')

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

    def build(self, graph='tvirw',verbose=False,difftol=0.15,multi=False):
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

        This function build all the graph associated with the Layout. 

        Warning : by default the layout is saved (dumpw) after each build

        """
        # list of built graphs
        if not self.hasboundary:
            self.boundary()

        # to save graoh Gs
        self.lbltg.extend('s')

        Buildpbar = pbar(verbose,total=5,desc='Build Layout',position=0)

        if verbose:
            Buildpbar.update(1)
        if 't' in graph:
            self.buildGt(difftol=difftol,verbose=verbose,tqdmpos=1)
            self.lbltg.extend('t')
        if verbose:
            Buildpbar.update(1)
        if 'v' in graph:
            self.buildGv(verbose=verbose,tqdmpos=1)
            self.lbltg.extend('v')
        if verbose:
            Buildpbar.update(1)
        if 'i' in graph:
            self.buildGi(verbose=verbose,tqdmpos=1)
            #pdb.set_trace()
            if not multi:
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

        fileini = pyu.getlong(self._filename, pro.pstruc['DIRINI'])
        _hash = hashlib.md5(open(fileini, 'rb').read()).hexdigest()
        self.Gt.add_node(0, hash=_hash)

        # There is a dumpw after each build
        self.dumpw()
        self.isbuilt = True
        if verbose:
            Buildpbar.update(1)

    def dumpw(self):
        """ write a dump of given Graph

        Notes
        -----

        't' : Gt
        'r' : Gr
        's' : Gs
        'v' : Gv
        'i' : Gi

        """
        # create layout directory
        dirname = self._filename.replace('.ini','')
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
        # save dictionnary which maps string interaction to [interaction node,
        # interaction type]
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
        dirname = self._filename.replace('.ini','')
        path = os.path.join(pro.basename, 'struc', 'gpickle', dirname)
        for g in graphs:
            try:
                # if g in ['v','i']:
                #     gname1 ='G'+g
                #     setattr(self, gname1, read_gpickle(os.path.join(pro.basename,'struc','gpickle','G'+g+'_'+self._filename+'.gpickle')))
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
        # pdb.set_trace()
        if 's' in graphs:
            #self._hash = self.Gs.node.pop(0)['hash']
            # self._hash = self.Gs.node[0]['hash']
            # update self.name
            lseg = [x for x in self.Gs.node if x > 0]
            for name in self.name:
                self.name[name] = [
                    x for x in lseg if self.Gs.node[x]['name'] == name]

            self.g2npy()

        
            filediff = os.path.join(path, 'ddiff.gpickle')
            if os.path.isfile(filediff):
                ddiff = read_gpickle(filediff)
                setattr(self, 'ddiff', ddiff)
                self.diffraction = True
            else:
                self.ddiff={}
                self.diffraction=False
            filelnss = os.path.join(path, 'lnss.gpickle')
            if os.path.isfile(filelnss):
                lnss = read_gpickle(filelnss)
                setattr(self, 'lnss', lnss) 
            else : 
                self.lnss=[]

        filedca = os.path.join(path, 'dca.gpickle')
        if os.path.isfile(filedca):
            dca = read_gpickle(filedca)
            setattr(self, 'dca',dca) 

        filem = os.path.join(path, 'm.gpickle')
        if os.path.isfile(filem):
            setattr(self, 'm', read_gpickle(filem))

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
        fig
        ax 
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
        """
        if fig == []:
            fig = plt.gcf()
        if ax == []:
            ax = plt.gca()
        try:
            mpl = [PolygonPatch(x, alpha=alpha, color=color) for x in poly]
        except:
            mpl = [PolygonPatch(x, alpha=alpha, color=color) for x in [poly]]
        [ax.add_patch(x) for x in mpl]
        plt.axis(self.ax)
        plt.draw()

    def pltvnodes(self, vn, fig=[], ax=[],):
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

        seg_connect = {x: self.Gs.node[x]['connect']
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


        Return
        ------
            T : dict 
                dictionnary from triangle.triangulate library
                >>> T.keys()
                ['segment_markers', 'segments', 'holes', 'vertices', 'vertex_markers', 'triangles']


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

    def _triangle(self, holes=[], vnodes=[]):
        """
        perfom a Delaunay partitioning on shapely polygons

        Parameters
        ----------

            holes : ndarray
                if holes ==[] : it means the merge is applied on the interior of the layout (indoor)
                if holes == np.ndarray (centroid of polygon). indoor is discarded and delaunay
                        is applied on outdoor


        Return
        ------
            T : dict 
                dictionnary from triangle.triangulate library
                >>> T.keys()
                ['segment_markers', 'segments', 'holes', 'vertices', 'vertex_markers', 'triangles']


        Notes
        -----

        This methoc uses the triangle library

        """

        # this means Delaunay is applied on exterior
        # and inside polygon will be discarded
        segbounds = []
        ptbounds = []
        if holes == []:
            # remove air segments around layout
            pass
            # [segbounds.extend(nx.neighbors(L.Gs,x)) for x in L.lboundary]
            # ptbounds = L.lboundary

        if vnodes == []:
            vnodes = self.Gs.nodes()
        # find segments of layout
        seg = np.array([nx.neighbors(self.Gs, x) for x in vnodes
                        if x > 0
                        and x not in segbounds])
        # get vertices/points of layout
        ivertices = np.array([(x, self.Gs.pos[x][0], self.Gs.pos[x][1]) for x in vnodes
                              if x < 0
                              and x not in ptbounds])
        map_vertices = ivertices[:, 0]
        vertices = ivertices[:, 1:]
        sorter = np.argsort(map_vertices)

        # mapping between Gs graph segments and triangle segments
        segments = sorter[np.searchsorted(map_vertices, seg, sorter=sorter)]

        if holes == []:
            C = {'vertices': vertices, 'segments': segments}
        else:
            C = {'vertices': vertices, 'segments': segments, 'holes': holes}

        T = triangle.triangulate(C, 'pa')

        # import triangle.plot as plot
        # ax=plt.gca()
        # plot.plot(ax,**T)
        # plt.show()
        return T, map_vertices

    def buildGt(self, check=True,difftol=0.01,verbose=False,tqdmpos=0):
        """ build graph of convex cycle 

        Parameters
        ----------

        check : boolean 
        difftol : float 
        verbose : boolean 
        tqdmpos : progressbar 

        todo :
        - add an option to only take outside polygon 
            => pass to self._triangle a hole coreesponding to centroid of
            polygon except those of boundary ( see buildGtold )

        """

        # 1. Do a Delaunay triangulation
        # build a list of triangle polygons : lTP
        # vnodes refers to the nodes of Gs
        # if vnodes == 0 it means this is a created
        # segment which is tagged as _AIR
        ###



        # if verbose :
        #     Gtpbar = tqdm.tqdm(total=100., desc='BuildGt',position=0)
        #     pbar_awloop =  tqdm.tqdm(total=100., desc ='airwalls loop',leave=False,position=1)
        Gtpbar = pbar(verbose,total=100., desc ='BuildGt',position=tqdmpos)
        pbartmp = pbar(verbose,total=100., desc ='Triangulation',leave=True,position=tqdmpos+1)

        T, map_vertices = self._triangle()

        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)
        # point index are integer
        map_vertices = map_vertices.astype(int)
        ptri = T['vertices'][T['triangles']]

        # List of Triangle Polygons
        pbartmp = pbar(verbose,total=100., 
                        desc ='Transfer polygons list',
                        leave=True,
                        position=tqdmpos+1)

        lTP = [geu.Polygon(x) for x in ptri]
        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        # update vnodes of Polygons
        pbartmp = pbar(verbose,total=100., 
                        desc ='Update Polygons vnodes',
                        leave=True,
                        position=tqdmpos+1)

        # p is a polygon 
        # get_points(p) : get points from polygon
        # this is for limiting the search region for large Layout 

        [ p.setvnodes_new(self.get_points(p),self) for p in lTP ]

        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)


        # 2.add air walls to triangle poly
        ###
        # luaw  : list of tuples
        # ( polygon , array of _AIR segments)
        pbartmp = pbar(verbose,total=100., 
                        desc ='Buiild list of airwalls',
                        leave=True,
                        position=tqdmpos+1)
        luaw = [(p, np.where(p.vnodes == 0)[0]) for p in lTP]
        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)



        #
        # For a triangle polygon the number of vnodes
        # creates new _AIR segments
        #
        cpt = 1./(len(luaw)+1)
        _airseg = []

        pbartmp = pbar(verbose,total=100., desc ='Add airwalls',leave=True,position=tqdmpos+1)

        for p, uaw in luaw:
            # for each vnodes == 0, add an _AIR
            if verbose :
                pbartmp.update(100.*cpt)
            for aw in uaw:
                modpt = len(p.vnodes)
                _airseg.append(self.add_segment(p.vnodes[np.mod(aw - 1, modpt)],
                                                p.vnodes[
                                                    np.mod(aw + 1, modpt)], name='_AIR',
                                                z=(0, 40000000),
                                                verbose=False))
            # update polygon segments with new added airwalls
            p.setvnodes_new(self.get_points(p),self)
        if verbose:
            Gtpbar.update(100./12.)


        pbartmp = pbar(verbose,total=100., desc ='Update Graph',leave=True,position=tqdmpos+1)


        tri = T['triangles']
        nbtri = len(T['triangles'])
        # temporary name/node_index of triangles
        MT = -np.arange(1, nbtri + 1)

        # 3. Create a temporary graph
        # where : positive nodes (>0) are triangles segments
        # negative nodes (<0) are triangles centroids
        # edges link triangle centroids to their respective segments

        # Ex represent list of points in Gs corresponging to segments
        #[pt_head pt_tail]

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
        neigh = [nx.neighbors(G, un) for un in rn]

        # store into networkx compliant format

        uE = [(neigh[un][0], neigh[un][1], {'segment': [
               rn[un]] + self.Gs.node[rn[un]]['iso']}) for un in xrange(len(rn))]
        iuE = {rn[un]: [-neigh[un][0], -neigh[un][1]]
               for un in xrange(len(rn))}

        # delete temporary graph
        del G

        # pdb.set_trace()

        # create graph Gt
        self.Gt = nx.Graph()
        self.Gt.add_edges_from(uE)
        self.Gt = nx.relabel_nodes(self.Gt, lambda x: -x)

        # add polyg  to nodes
        # add indoor to nodes
        # add isopen to nodes

        nno = [(n, {'polyg': lTP[n - 1], 'indoor':True, 'isopen':True})
               for n in self.Gt.nodes()]

        self.Gt.add_nodes_from(nno)
        self.Gt.pos = {}
        self.Gt.pos.update({n: np.array(
            self.Gt.node[n]['polyg'].centroid.xy).squeeze() for n in self.Gt.nodes()})

        # self.Gtpos = {-MT[i]:pMT[i] for i in xrange(len(MT))}
        # plt.figure()
        # # G=nx.Graph()
        # # G.add_edges_from(E0)
        # # G.add_edges_from(E1)
        # # G.add_edges_from(E2)

        _airseg = np.unique(_airseg)
        _airseg = _airseg[_airseg != np.array(None)].astype('int')

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
            n0, n1 = iuE[a]
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

            p0 = self.Gt.node[n0]['polyg']
            p1 = self.Gt.node[n1]['polyg']

            # Merge polygon
            P = p0 + p1
            # If the new Polygon is convex update Gt
            #
            if geu.isconvex(P):
                # updates vnodes of the new merged polygon
                P.setvnodes_new(self.get_points(P),self)
                # update edge
                n0s = n0
                n1s = n1
                # get segments information from cycle n0
                dne = self.Gt[n0]
                # remove connection to n0 to avoid a cycle being
                # connected to itself
                self.Gt[n1].pop(n0)
                # add information from adjacent cycle n1
                dne.update(self.Gt[n1])
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
            self.Gs.node[s]['ncycles'].append(0)

        #
        # boundary adjascent cycles
        #
        adjcyair = np.array(map(lambda x: filter(lambda y: y != 0,
                                                 self.Gs.node[x]['ncycles'])[0], self.segboundary))
        # connect cycles separated by air wall to cycle 0
        for cy, seg in zip(adjcyair, self.segboundary):
            self.Gt.node[cy]['indoor'] = False
            self.Gt.node[cy]['isopen'] = True
            self.Gt.add_edge(0, cy, segment=[seg])
        
        # 
        #
        #
        if check:
            # print("check len(ncycles) == 2",)
            nodes = [i for i in self.Gs.nodes() if i > 0]
            cncy = np.array([len(self.Gs.node[i]['ncycles']) for i in nodes])
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
        if self.diffraction:
            tqdmkwargs={'total':100.,'desc':'Find Diffractions','position':1}
            self._find_diffractions(difftol=difftol,verbose=verbose,tqdmkwargs=tqdmkwargs)
            if verbose:
            # print('find diffraction...Done 8/12')
                Gtpbar.update(100./12.)
        # 
        # explanation of lnss
        #
        # list of diffraction point involving different segment 
        # list of diffraction point involving subsegment ( = iso segments)
        # needs checking height in rays.to3D for constructing the 3D ray
        #
            pbartmp = pbar(verbose,total=100., desc ='Diffraction on airwalls',leave=True,position=tqdmpos+1)

            self.lnss = [x for x in self.ddiff if len(
            set(nx.neighbors(self.Gs, x)).intersection(set(self.lsss))) > 0]


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
        to_visit = nx.neighbors(self.Gt, 0)
        law = self.name['_AIR'] + self.name['AIR']
        while len(to_visit) > 0:
            # get current cycle
            cur_cy = to_visit.pop()
            # get neighbors of current_cycle
            neighbors = nx.neighbors(self.Gt, cur_cy)
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
                self.Gt.node[x]['indoor'] = False
                self.Gt.node[x]['isopen'] = True
            # extend to_visit to not visited neighbors
            to_visit.extend(nv_neighbors_aw)
            visited.append(cur_cy)
        if verbose:
            pbartmp.update(100.)
            Gtpbar.update(100./12.)

        self.g2npy()

    def _visual_check(self):
        """ visual checking of graphs
        """
        fig, axs = plt.subplots(2, 2)

        try:
            ax = axs[0, 0]
            self.showG('s', aw=1, ax=ax, fig=fig)
            indoor = [self.Gt.node[p]['polyg']
                      for p in self.Gt.nodes() if p != 0 and self.Gt.node[p]['indoor']]
            outdoor = [self.Gt.node[p]['polyg'] for p in self.Gt.nodes(
            ) if p != 0 and not self.Gt.node[p]['indoor']]
            self.pltpoly(indoor, color='r', ax=ax, fig=fig)
            self.pltpoly(outdoor, color='g', ax=ax, fig=fig)
            ax.set_title('indoor red,outdoor green')
        except:
            print('error with polyg in Gt or indoor parameter')

        # try:

        #     ax = axs[0,1]
        #     f,ax=self.showG('s',aw=1,ax=ax,fig=fig)
        #     isopen = [self.Gt.node[p]['polyg'] for p in self.Gt.nodes() if p !=0 and self.Gt.node[p]['isopen'] ]
        #     self.pltpoly(isopen,ax=ax,fig=fig)
        #     ax.set_title('isopen red')
        # except:
        #     print('error with isopen parameter in Gt')

        try:
            ax = axs[0, 1]
            f, ax = self.showG('s', aw=1, ax=ax, fig=fig)
            diffpos = np.array([self.Gs.pos[x] for x in self.ddiff.keys()])
            ax.scatter(diffpos[:, 0], diffpos[:, 1])
            ax.set_title('diffraction points')
        except:
            print('no diffraction found. Yet computed ?')

        try:
            ax = axs[1, 0]
            f, ax = self.showG('st', labels='t', aw=1, ax=ax, fig=fig)
            ax.set_title('Gt')
        except:
            print('no Gt found. Yet computed ?')

        try:
            ax = axs[1, 1]
            f, ax = self.showG('sv', aw=1, ax=ax, fig=fig)
            ax.set_title('Gv')
        except:
            print('no Gv found. Yet computed ?')
        plt.tight_layout()
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

    def buildGt_old(self, check=True):
        """ 
        DEPRECATED
        build graph of convex cycles 


        Parameters
        ----------

        check : booolean


        """

        def pltpoly(poly, fig=[], ax=[]):
            if fig == []:
                fig = plt.gcf()
            if ax == []:
                ax = plt.gca()
            mpl = [PolygonPatch(x, alpha=0.2) for x in poly]
            [ax.add_patch(x) for x in mpl]
            plt.axis(self.ax)
            plt.draw()

        def pltGt(Gt, fig=[], ax=[]):
            if fig == []:
                fig = plt.gcf()
            if ax == []:
                ax = plt.gca()
            [Gt.node[x]['poly'].plot(fig=fig, ax=ax, alpha=0.2)
             for x in Gt.nodes()]
            plt.axis(self.ax)
            plt.draw()

        # I. get cycle basis
        C = nx.algorithms.cycles.cycle_basis(self.Gs)
        if C == []:
            C = [self.Gs]

        # pdb.set_trace()
        # II. create the hull of the layout by merging all polygons
        # corresponding to cycles basis
        poly = []
        for k, lnode in enumerate(C):
            npoints = filter(lambda x: x < 0, lnode)
            coords = map(lambda x: self.Gs.pos[x], npoints)
            poly.append(sh.Polygon(coords))
        # union all polygons
        # pdb.set_trace()
        ma = cascaded_union(poly)

        # transform into geomutil polygon
        # if  polygon is a layout
        if not isinstance(ma, sh.MultiPolygon):
            ma = geu.Polygon(ma)
            ma.setvnodes(self)

        else:
            # This is a fix for non enclosed layouts
            # with multiple non joint polygons (a.k.a. a city)
            # raise AttributeError('this is a city')
            macvx = ma.convex_hull

            streets = geu.Polygon(macvx.difference(ma))
            streets.setvnodes(self)

            # add air walls where to close the street & ma polygon
            uaw = np.where(streets.vnodes == 0)[0]
            lvn = len(streets.vnodes)
            for i in uaw:
                # keep trace of created airwalls, because some
                # of them will be destroyed in step 3.
                self.add_segment(
                    streets.vnodes[np.mod(i - 1, lvn)],
                    streets.vnodes[np.mod(i + 1, lvn)], name='_AIR')

            ma = self.polysh2geu(macvx)
            ma.setvnodes(self)
            streets.setvnodes(self)
            ma.setvnodes(self)

        self.ma = ma

        # III .FIND POLYGONS
        ###
        # polygons of each cycle are found by finding the interesection between
        # all segments of the layout and the layout hull.
        # The shapely diff return a multipolygon where all polygons corresponds to
        # a cycle
        #

        # get connected points from segments
        # connect is equivalent to self.tahe and lpos to self.pt
        #
        connect = [self.Gs.node[i]['connect']
                   for i in self.Gs.nodes() if i > 0]
        # get their coordinates
        lpos = np.array([(self.Gs.pos[i[0]], self.Gs.pos[i[1]])
                         for i in connect])
        pp = []
        lines = []

        for l in lpos:
            line = sh.LineString([l[0], l[1]])
            lines.append(line)
        # create associated multilines
        ml = sh.MultiLineString(lines)
        # increase buffer size ( width of polyline) to create a multipolygon
        R = []
        buffersize = 1e-9
        while not isinstance(R, sh.MultiPolygon) and buffersize < 1e-3:
            # create polygon from multiline by given a width to lines
            mlp = ml.buffer(buffersize)
            # the difference between the layout hull and polygons built from lines
            # returns the ndesired multipolygon
            R = self.ma.difference(mlp)
            # increase size of the buffer
            buffersize = buffersize * 10

        if isinstance(R, sh.Polygon):
            R = sh.MultiPolygon([R])

        # assert isinstance(R,sh.MultiPolygon), "Shapely.MultiPolygon decomposition Failed"

        ####################
        # Manage inner hole in polygons
        # ------------------------------
        # This part manages layout not correctly described, where
        # polygons remains in the middle of others
        ######

        # if !=0 it means some polygons are inside of others
        # which is not allowed. Some Layout modification will be performed

        Rgeu = []
        contain = {}

        for ur, r in enumerate(R):
            try:

                Rgeu.append(self.polysh2geu(r))
                # self.pltpoly([Rgeu[-1]],color='b')
                # plt.draw()
                # Rgeu.append(geu.Polygon(r))
                # Rgeu[-1].setvnodes(self)
            except:
                print("reject")
            # if area are not the same, it means that there is inner holes in r
            if not np.allclose(Rgeu[ur].area, r.area):
                # detect inclusion
                uc = np.where([Rgeu[ur].contains(p) for p in R])[0]
                contain[ur] = [c for c in uc if c != ur]

        # # split polygons with holes into several polygons without holes

        self.macvx = []
        for k in contain:

            polyholes = [Rgeu[i] for i in contain[k]]

            # 1 convexify polyholes
            polyg = self._convex_hull(polyholes)
            polyholes.extend(polyg)
            Rgeu.extend(polyg)
            # 2 delaunay on exterior
            ncpol = self._delaunay(Rgeu[k], polyholes=polyholes)

            Rgeu.pop(k)
            Rgeu.extend(ncpol)

            # add polyhole to convex mask ( macvx) list
            macvx = cascaded_union(polyholes)
            macvx = geu.Polygon(macvx)
            macvx.setvnodes(self)
            self.macvx.append(macvx)

        ####################
        # Manage  convex hull of the layout
        # -------------------

        # polys = self._convex_hull()
        # Rgeu.extend(polys)

        ####################
        # Manage Non convex polygons
        # -------------------
        # 1 . determine which polygons are not convex
        # 2 . apply a delaunay and tranform a single non convexpolygon
        # into several convex ( self.delaunay)
        # 3. remove old non convex polygon and readd new convex ones.

        ncpol = {}
        for ur, r in enumerate(Rgeu):
            if not r.isconvex():
                print('nt cvx')
                ncpol[ur] = self._delaunay(r)

                # plt.ion()
                # self.pltpoly(ncpol[ur],fig=plt.gcf(),ax=plt.gca())
                # plt.show()
                # plt.draw()
                # import ipdb
                # ipdb.set_trace()
        Rgeu = np.delete(Rgeu, ncpol.keys()).tolist()
        [Rgeu.extend(ncpol[k]) for k in ncpol]

        self.Gt = nx.Graph()
        self.Gt.pos = {}

        # IV Find Vnodes and Final polygons

        for n in self.Gs.node:
            if n > 0:
                self.Gs.node[n]['ncycles'] = []

        ncyid = -1

        sma = self.ma.vnodes[self.ma.vnodes > 0]
        # smac = self.macvx.vnodes[self.macvx.vnodes>0]
        # segma = np.unique(np.concatenate((sma,smac)))
        segma = sma
        # VI  add node 0
        #
        #   This shapely polygon has an interior
        #    Cycles = 0 exterior cycle (assumed outdoor)

        S = nx.subgraph(self.Gs, self.ma.vnodes)
        S.pos = {}
        S.pos.update({i: self.Gs.pos[i] for i in S.nodes()})
        cycle = cycl.Cycle(S, self.ma.vnodes)
        boundary = geu.Polygon(tuple(self.ax), delta=5)
        boundary.vnodes = self.ma.vnodes
        self.Gt.add_node(0, polyg=self.ma, cycle=cycle,
                         indoor=False, isopen=True)
        self.Gt.pos[0] = (self.ax[0], self.ax[2])

        # IV 1 get nodes and vnodes

        for ui, p in enumerate(Rgeu):
            cyid = ui + 1
            outdoor = False
            # IV 1.a get vnode associated to the polygon
            # get vnodes not in the correct order
            # uvn = np.where([r.buffer(1e-3).contains(p) for p in shpt])[0]
            # vnodes = Gsnodes[uvn]

            # IV 1.b transform vnodes to an ordered cycle with Cycle class
            # NOTE ! Using class cycle is MANDATORY
            # because, some extra vnodes can be picked up during the contain
            # process before

            S = nx.subgraph(self.Gs, p.vnodes)
            S.pos = {}
            S.pos.update({i: self.Gs.pos[i] for i in S.nodes()})

            cycle = cycl.Cycle(S)

            # IV 1.c create a new polygon with correct vnodes and correct points
            # P = geu.Polygon(p=cycle.p,vnodes=cycle.cycle)
            # import ipdb
            # ipdb.set_trace()
            # IV 1.d add node to Gt + position
            #

            seg = p.vnodes[p.vnodes > 0]
            lair = [x in self.name['AIR'] for x in seg]

            if sum(lair) > 0:
                isopen = True
            else:
                isopen = False

            # IV 1.e
            #   + add new node (convex cycle) to Gt
            #   + add centroid of cycle as position of cycle
            # if ((cyid==40) or (cyid==41)):
            #     pdb.set_trace()
            self.Gt.add_node(cyid, cycle=cycle, polyg=p,
                             isopen=isopen, indoor=True)
            self.Gt.pos.update({cyid: np.array(p.centroid.xy)[:, 0]})

        # IV 2. get edges
        for n1 in self.Gt.nodes():
            for n2 in self.Gt.nodes():
                if n1 != n2:
                    if self.Gt.node[n1]['polyg'].buffer(1e-3).touches(self.Gt.node[n2]['polyg']):
                        # find common segments
                        seg = np.array([n for n in self.Gt.node[n1]['cycle'].cycle if (
                            n in self.Gt.node[n2]['cycle'].cycle) and (n > 0)])
                        # if cycle are connected by at least a segmnet but not
                        # a point
                        if len(seg) > 0:
                            self.Gt.add_edge(n1, n2, segment=seg)

        # import ipdb
        # ipdb.set_trace()
        #  V update Gs
        #   V 1.Update graph Gs nodes with their cycles information
        #
        #   initialize a void list 'ncycles' for each node of Gs
        #
        self._updGsncy()
        # make a convex hull of layout
        # get segments of the mask ( equivalent to thoose connected to 0)
        # seg0 = [i for i in self.ma.vnodes if i >0]
        # [self.Gs.node[i]['ncycles'].append(0) for i in seg0]

        self._addoutcy(check)

        #   V 2. add outside cycle (absorbant region index 0 )
        #   if ncycles is a list which has only one element then the adjascent cycle is the
        #   outside region (cycle 0)
        #
        #   VII - Connect cycle 0 to each cycle connected to the layout
        #   boundary
        #

        # all segments of the Layout boundary
        nseg = filter(lambda x: x > 0, boundary.vnodes)
        # air segments of the Layout boundary
        nsegair = filter(lambda x: x in self.name['AIR'], nseg)
        # wall segments of the Layout boundary
        nsegwall = filter(lambda x: x not in self.name['AIR'], nseg)

        #
        # boundary adjascent cycles
        #

        adjcyair = np.array(map(lambda x: filter(lambda y: y != 0,
                                                 self.Gs.node[x]['ncycles'])[0], nsegair))
        adjcwall = np.array(map(lambda x: filter(lambda y: y != 0,
                                                 self.Gs.node[x]['ncycles'])[0], nsegwall))
        # pdb.set_trace()
        # adjcyair = np.unique(adjcyair)
        # adjcwall = np.unique(adjcwall)

        # connect cycles separated by air wall to cycle 0
        for cy, seg in zip(adjcyair, nsegair):
            self.Gt.node[cy]['indoor'] = False
            self.Gt.node[cy]['isopen'] = True
            self.Gt.add_edge(0, cy, segment=np.array([seg]))

        # connect cycles separated by wall to cycle 0
        for cy, seg in zip(adjcwall, nsegwall):
            self.Gt.add_edge(0, cy, segment=np.array([seg]))

        # IV Handle indoor/outdoor cycle
        #
        # Rule : A cycle is outdoor if it is separated from an outdoor cycle by an AIR segment
        #
        # An outdoor cycle has no associated ceil
        #
        for cy in self.Gt.nodes():
            lncy = nx.neighbors(self.Gt, cy)
            for ncy in lncy:
                segnum = self.Gt.edge[cy][ncy]['segment'][0]
                if self.Gs.node[segnum]['name'] == 'AIR':
                    if not self.Gt.node[cy]['indoor']:
                        self.Gt.node[ncy]['indoor'] = False

        self._find_diffractions()
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

        self._interlist()

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
            self.Gs.node[k]['ncycles'] = []

        # filter out node 0
        Gtnodes = filter(lambda x: x != 0, self.Gt.nodes())

        # loop over all cycles
        for ncy in Gtnodes:
            # get vnodes : points and segments number
            vnodes = self.Gt.node[ncy]['polyg'].vnodes
            for n in vnodes:
                if n == 0:
                    pdb.set_trace()
                if ncy not in self.Gs.node[n]['ncycles']:
                    self.Gs.node[n]['ncycles'].append(ncy)
                    if n > 0:
                        if len(self.Gs.node[n]['ncycles']) > 2:
                            print(n, self.Gs.node[n]['ncycles'])
                            logging.warning(
                                'A segment cannot relate more than 2 cycles')

        for nseg in self.Gs.node:
            if nseg > 0:
                ncycles = self.Gs.node[nseg]['ncycles']
                if len(ncycles) > 1:
                    if nseg not in self.Gt.edge[ncycles[0]][ncycles[1]]['segment']:
                        self.Gt.edge[ncycles[0]][ncycles[1]][
                            'segment'].append(nseg)

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
            cncy = np.array([len(self.Gs.node[i]['ncycles']) for i in nodes])
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
                if self.indoor or not self.Gt.node[k]['indoor']:
                    #vnodes = self.Gt.node[k]['vnodes']
                    vnodes = self.Gt.node[k]['polyg'].vnodes
                    ListInteractions = []
                    for inode in vnodes:
                        if inode > 0:   # segments
                            cy = set(self.Gs.node[inode]['ncycles'])
                            name = self.Gs.node[inode]['name']  # segment name
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
                                if (name != 'METAL') & (name != 'ABSORBENT'):
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

        self.Gw = nx.Graph()
        self.Gw.pos = {}

        d_id = max(self.Gr.nodes())  # for numerotation of Gw nodes
        d_id_index = d_id + 1

        for e in self.Gr.edges_iter():  # iterator on Gr edges

            self.Gw.add_node(e[0], {'room': e[0], 'door': False})
            self.Gw.add_node(e[1], {'room': e[1], 'door': False})

            # transitions of room e[0]
            trans1 = self.Gr.node[e[0]]['transition']
            # transitions of room e[1]
            trans2 = self.Gr.node[e[1]]['transition']
            Id = np.intersect1d(trans1, trans2)[0]  # list of common doors

            unode = self.Gs.neighbors(Id)  # get edge number of common doors
            up0 = self.Gs.pos[unode[0]]
            up1 = self.Gs.pos[unode[1]]

            name = self.Gs.node[Id]['name']

            if True:  # name != 'AIR':
                pn = self.Gs.node[Id]['norm']
                sl = self.sl[name]
                thick = (sum(sl['lthick']) / 2.) + 0.2

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
                    self.Gw.add_node(upd0, {'room': e[0], 'door': True})
                    # if self.seginline(pdoor0,ep0).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[0],upd0)])
                    d_id_index = d_id_index + 1

                    upd1 = d_id_index
                    self.Gw.pos[upd1] = pdoor1
                    self.Gw.add_node(upd1, {'room': e[1], 'door': True})
                    # if self.seginline(pdoor1,ep1).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[1],upd1)])
                    d_id_index = d_id_index + 1
                else:
                    upd0 = d_id_index
                    self.Gw.pos[upd0] = pdoor0
                    self.Gw.add_node(upd0, {'room': e[1], 'door': True})
                    # if self.seginline(pdoor0,ep1).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[1],upd0)])
                    d_id_index = d_id_index + 1

                    upd1 = d_id_index
                    self.Gw.pos[upd1] = pdoor1
                    self.Gw.add_node(upd1, {'room': e[0], 'door': True})
                    # if self.seginline(pdoor1,ep0).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[0],upd1)])
                    d_id_index = d_id_index + 1
                self.Gw.add_edges_from([(upd0, upd1)])
            # Airwalls case
            else:

                pdoor = (np.array(up0) + np.array(up1)) / 2
                self.Gw.pos[d_id_index] = pdoor
                self.Gw.add_edges_from(
                    [(e[0], d_id_index), (e[1], d_id_index)])
                d_id_index = d_id_index + 1

        self.Gw.pos.update(self.Gr.pos)

#####################
        # # SINGLE POINT SEGMENT TRANSITION
        # d_id = max(self.Gr.nodes())+1
        # # Add trnasistion segments in visibility
        # for e in self.Gr.edges_iter(): # iterator on Gr edges
        #     trans1 = self.Gr.node[e[0]]['transition']  # transitions of room e[0]
        #     trans2 = self.Gr.node[e[1]]['transition']  # transitions of room e[1]
        #     Id = np.intersect1d(trans1,trans2)[0]  # list of common doors

        #     unode = self.Gs.neighbors(Id) # get edge number of common doors
        #     p1 = self.Gs.pos[unode[0]]
        #     p2 = self.Gs.pos[unode[1]]
        # pdoor = (np.array(p1) + np.array(p2)) / 2  # middle of the common
        # door

        #     self.Gw.add_node(Id + d_id,door=True)     # new node
        #     self.Gw.pos[Id + d_id] = pdoor  # in the middle of the door|airwall
        #     pe0=self.Gr.pos[e[0]]
        #     pe1=self.Gr.pos[e[1]]

        #     # if only the trnasisiton intersects
        #     if self.seginline(pdoor,pe0).shape[1] == 1:
        #         self.Gw.add_edges_from([(e[0], Id + d_id)])
        #     if self.seginline(pdoor,pe1).shape[1] == 1:
        #         self.Gw.add_edges_from([(e[1], Id + d_id)])
        #     self.Gw.pos.update(self.Gr.pos)
####################
        # ADD CONVEX POINTS
        # d_id = max(self.Gw.nodes())+1
        # pcid = d_id
        # tmp = []
        # for n in self.Gr.nodes():
        #     # get segment number of the room
        #     tcc, nn = self.Gr.node[n]['polyg'].ptconvex()
        #     uconvex = np.nonzero(tcc == 1)[0]

        #     if len(uconvex) != 0 :
        #         lr = self.Gr.node[n]['polyg'].exterior
        #         x,y = lr.xy
        #         p = [np.array((x[i],y[i])) for i in uconvex]
        #         ln = self.Gr.node[n]['cycle'].cycle
        #         lp = filter(lambda x: x< 0, ln)
        #         pc = [lp[i] for i in uconvex]
        #         for uu,uc in enumerate(uconvex):
        #             # convex point position take into account wall width
        #             npc = nx.neighbors(self.Gs,pc[uu])
        #             # only corner are considered (point with 3 neighbors are not)
        #             if len(npc) <=2:
        #                 nname = [self.Gs.node[nd]['name'] for nd in npc]
        #                 npos = [self.Gs.pos[nd] for nd in npc]
        #                 nv = [p[uu]- i for i in npos]
        #                 nvn = sum(nv)/abs(sum(nv))
        #                 nsl = [self.sl[name] for name in nname]
        #                 thick = sum([sum(sl['lthick']) for sl in nsl])
        #                 # vector to add to the convex point position
        #                 v = nvn*thick
        #                 pid = uc+pcid
        #                 self.Gw.add_node(pid,{'diff':True,'room':n,'door':False})
        #                 self.Gw.pos.update({pid:p[uu]+v[:2]})
        #                 pcid = pcid +1

        # pdb.set_trace()
        for n in self.Gr.nodes():

            nr = nx.get_node_attributes(self.Gw, 'room').items()
            # get the Gw nodes in current room
            f = map(lambda y: y[0], filter(lambda x: x[1] == n, nr))
            for nw in combinations(f, 2):
                pf = map(lambda x: self.Gw.pos[x], nw)
                pf = np.array((pf))
                # if self.seginline(pf[0],pf[1]).shape[1] <= 1:
                d = np.sqrt(np.sum((pf[0] - pf[1])**2))
                self.Gw.add_edges_from([(nw[0], nw[1])], weight=d)

            # kudr = [kdr[u] for u in udr]
            # cdr = combinations(dr.keys()[udr],2)
            # for
            # import ipdb
            # ipdb.set_trace()

        # for n in self.Gr.nodes_iter():
        #     d = self.Gw.neighbors(n)   # neighbors of room n in Gw
        #     if len(d) > 1:
        #         self.Gw.add_edges_from(combinations(d, 2))

        # udn = nx.get_node_attributes(self.Gw,'diff').keys()
        # import ipdb
        # ipdb.set_trace()
        # for u in udn:
        #     [self.Gw[u].update({i:{'weigth':0.01}})for i in self.Gw.edge[u].keys()]
    # def buildGw(self):
    #     """ build Graph of waypaths

    #     See Also
    #     --------

    #     buildGr

    #     Notes
    #     -----

    #     for all edges of Gr (adjascent room)
    #         if room1 and room2 have a common transition
    #     """

    #     self.Gw = nx.Graph()
    #     self.Gw.pos = {}

    #     d_id = max(self.Gr.nodes()) # for numerotation of Gw nodes

    #     for e in self.Gr.edges_iter(): # iterator on Gr edges
    #         trans1 = self.Gr.node[e[0]]['transition']  # transitions of room e[0]
    #         trans2 = self.Gr.node[e[1]]['transition']  # transitions of room e[1]
    #         Id = np.intersect1d(trans1,trans2)[0]  # list of common doors

    #         unode = self.Gs.neighbors(Id) # get edge number of common doors
    #         p1 = self.Gs.pos[unode[0]]
    #         p2 = self.Gs.pos[unode[1]]
    # pdoor = (np.array(p1) + np.array(p2)) / 2  # middle of the common door

    #         self.Gw.add_node(Id + d_id,door=True)     # new node
    #         self.Gw.pos[Id + d_id] = pdoor  # in the middle of the door|airwall
    #         self.Gw.add_edges_from([(e[0], Id + d_id),
    #                                 (e[1], Id + d_id)])
    #         self.Gw.pos.update(self.Gr.pos)

    #     for n in self.Gr.nodes_iter():
    #         d = self.Gw.neighbors(n)   # neighbors of room n in Gw
    #         if len(d) > 1:
    #             self.Gw.add_edges_from(combinations(d, 2))

    # def buildGw2(self):
    #     """ build Graph of waypaths

    #     See Also
    #     --------

    #     buildGr

    #     Notes
    #     -----

    #     for all edges of Gr (adjascent room)
    #         if room1 and room2 have a common transition

    #     """

    #     self.Gw = nx.Graph()
    #     self.Gw.pos = {}

    #     d_id = max(self.Gr.nodes()) # for numerotation of Gw nodes
    #     d_id_index = d_id + 1
    #     for e in self.Gr.edges_iter(): # iterator on Gr edges
    #         trans1 = self.Gr.node[e[0]]['transition']  # transitions of room e[0]
    #         trans2 = self.Gr.node[e[1]]['transition']  # transitions of room e[1]
    #         Id = np.intersect1d(trans1,trans2)[0]  # list of common doors

    #         unode = self.Gs.neighbors(Id) # get edge number of common doors

    #         pn = self.Gs.node[Id]['norm']
    #         name = self.Gs.node[Id]['name']
    #         sl = self.sl[name]
    #         thick = sum(sl['lthick'])/2.

    #         p1 = self.Gs.pos[unode[0]]
    #         p2 = self.Gs.pos[unode[1]]

    #         # middle of the common door
    #         pdoor0 = (np.array(p1) + pn[:2] * thick +\
    #                   np.array(p2) + pn[:2] * thick) / 2
    #         pdoor1 = (np.array(p1) - pn[:2] * thick +\
    #                   np.array(p2) - pn[:2] * thick) / 2

    #         upd0 = d_id_index
    #         self.Gw.add_node(upd0)     # new node
    #         self.Gw.pos[upd0] = pdoor0
    #         d_id_index = d_id_index +1

    #         upd1 = d_id_index
    #         self.Gw.add_node(upd1)     # new node
    #         self.Gw.pos[upd1] = pdoor1
    #         d_id_index = d_id_index +1

    #         de0p0 = np.sqrt(np.sum((self.Gr.pos[e[0]]-pdoor0)**2))
    #         de0p1 = np.sqrt(np.sum((self.Gr.pos[e[0]]-pdoor1)**2))
    #         de1p0 = np.sqrt(np.sum((self.Gr.pos[e[1]]-pdoor0)**2))
    #         de1p1 = np.sqrt(np.sum((self.Gr.pos[e[1]]-pdoor1)**2))

    #         self.Gw.add_edges_from([(upd0, upd1)])
    #         if de0p0 < de0p1:
    #             self.Gw.add_edges_from([(e[0],upd0)])
    #         else:
    #             self.Gw.add_edges_from([(e[0],upd1)])
    #         if de1p0 < de1p1:
    #             self.Gw.add_edges_from([(e[1],upd0)])
    #         else:
    #             self.Gw.add_edges_from([(e[1],upd1)])

    #         self.Gw.pos.update(self.Gr.pos)
    #     for n in self.Gr.nodes_iter():
    #         d = self.Gw.neighbors(n)   # neighbors of room n in Gw
    #         if len(d) > 1:
    #             self.Gw.add_edges_from(combinations(d, 2))


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
        >>> L = Layout('TA-Office.ini')
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

        self.Gv = nx.Graph()
        #
        # loop over convex cycles (nodes of Gt)
        #
        self.dGv = {}  # dict of Gv graph

        cpt = 1./(len(self.Gt.node) + 1.)
        
        for icycle in self.Gt.node:
            if verbose:
                Gvpbar.update(100.*cpt)
            if icycle != 0:
                if self.indoor or not self.Gt.node[icycle]['indoor']:
                    #print(icycle)
                    pass
                #
                #  If indoor or outdoor all visibility are calculated
                #  If outdoor only visibility between iso = 'AIR' and '_AIR' are calculated 
                # 
                #if self.indoor or not self.Gt.node[icycle]['indoor']:
                polyg = self.Gt.node[icycle]['polyg']

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
                if ((not self.indoor) and (self.Gt.node[icycle]['indoor'])):
                    nseg = [ x for x in nseg_full if ((self.Gs.node[x]['name']=='AIR') or (self.Gs.node[x]['name']=='_AIR') ) ]
                else:
                    nseg = vnodes[useg]

                
                # # nseg_full : full list of segments
                # #nseg_full = filter(lambda x: x > 0, vnodes)

                # # keep only airwalls without iso single (_AIR)
                # nseg_single = filter(lambda x: len(self.Gs.node[x]['iso'])==0, nseg)

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

                Gv = nx.Graph()
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
                        if ((0 not in self.Gs.node[nk[0]]['ncycles']) and
                            (0 not in self.Gs.node[nk[1]]['ncycles'])):
                            # get the iso segments of both nk[0] and nk[1]
                            if ((self.indoor) or (not self.Gt.node[icycle]['indoor'])):
                                l0 = [nk[0]]+self.Gs.node[nk[0]]['iso']
                                l1 = [nk[1]]+self.Gs.node[nk[1]]['iso']
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
                if ((self.indoor) or (not self.Gt.node[icycle]['indoor'])):
                    ndiffvalid = [ x for x in ndiff if icycle in self.ddiff[x][0]]

                        # non adjascent segment of vnodes see valid diffraction
                        # points
                    for idiff in ndiffvalid:
                        #
                        # segments voisins du point de diffraction valide
                        #
                        nsneigh = [x for x in 
                                   nx.neighbors(self.Gs, idiff) 
                                   if x in nseg_full]
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
                        # and which are not neighbrs of the point idiff
                        #
                        for x in nsneigh:
                            neighbx = [ y for y in nx.neighbors(Gv, x) 
                                        if 0 not in self.Gs.node[y]['ncycles'] 
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

        self.Gi = nx.DiGraph()
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
            # espoo_journal debug
            #if n == 530:
            #    pdb.set_trace()
            if verbose:
                pbartmp.update(cpt)

            if n < 0:  # D
                self.Gi.add_node((n,))
                self.Gi.pos[(n,)] = self.Gs.pos[n]
            if n > 0:  # R | T
                cy = self.Gs.node[n]['ncycles']
                name = self.Gs.node[n]['name']
                assert(len(cy) == 2)
                cy0 = cy[0]
                cy1 = cy[1]

                nei = self.Gs.neighbors(n)  # get neighbor
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
                    self.Gi.pos[(n, cy0, cy1)] = tuple(
                        self.Gs.pos[n] + ln * delta / 2.)
                    self.Gi.pos[(n, cy1, cy0)] = tuple(
                        self.Gs.pos[n] - ln * delta / 2.)

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
                vnodes = self.Gt.node[cy]['polyg'].vnodes
                npt = []
                if self.diffraction:
                    #
                    # find all diffraction points involved in the cycle cy 
                    #
                    for x in vnodes:
                        if x < 0:
                            if self.ddiff.has_key(x):
                                for y in self.ddiff[x][0]:
                                    if y == cy:
                                        npt.append(x)
                    
                nseg = [ k for k in vnodes if k>0 ]
                if self.diffraction:
                # all segments and diffraction points of the cycle
                    vnodes = nseg + npt
                else:
                # only segments
                    vnodes = nseg

                for nstr in vnodes:

                    if nstr in self.Gv.nodes():
                        # list 1 of interactions
                        if nstr==108:
                            iprint = 1
                        else: 
                            iprint = 0
                        li1 = []
                        if nstr > 0:
                            # output cycle 
                            # cy -> cyo1 
                            cyo1 = self.Gs.node[nstr]['ncycles']
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
                        lneighb = nx.neighbors(self.Gv, nstr)
                        #if (self.Gs.node[nstr]['name']=='AIR') or (
                        #        self.Gs.node[nstr]['name']=='_AIR'):
                        #    lneighcy = lneighb
                        #else:
                        # list of cycle entities in visibility of nstr in the same cycle 
                        lneighcy = [ x for x in lneighb if x in vnodes ] 
                        # lneighcy = filter(lambda x: x in vnodes, lneighb)

                        for nstrb in lneighcy:
                            if nstrb in self.Gv.nodes():
                                li2 = []
                                if nstrb > 0:
                                    cyo2 = self.Gs.node[nstrb]['ncycles']
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
                                
                                for i1 in li1:
                                    # printli1
                                    for i2 in li2:
                                        # printli2
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
        # pdb.set_trace()
        pbartmp = pbar(verbose,total=100.,
                       desc ='update interraction list',
                       leave=False,
                       position=tqdmpos+1)

        for c in self.Gt.node:
            if verbose:
                pbartmp.update(cpt)
            if c != 0:
                vnodes = self.Gt.node[c]['polyg'].vnodes
                for k in npt:
                    self.Gt.node[c]['inter'] += [(k,)]

        if verbose :
            Gipbar.update(100.)

        # cleaning deadend Gi 
        # if not indoor for all nodes of Gi 
        # if not diffraction 
        # if termination cycle is indoor 
        # or if starting point is indoor 
        # then delte interaction 
        ldelete = []
        if not self.indoor:
            for k in self.Gi.node.keys():
                if len(k)>1:
                    segtype = self.Gs.node[k[0]]['name']
                    if ((segtype!='AIR') and (segtype!='_AIR')):
                        cyend = k[-1] 
                        if self.Gt.node[cyend]['indoor']:
                            # if k[0]>0:
                            #     if self.Gs.node[k[0]]['name']!='AIR':
                            ldelete.append(k)
                        if len(k) == 3:
                            cystart = k[1]
                            if self.Gt.node[cystart]['indoor']:
                                # if k[0]>0:
                                #     if self.Gs.node[k[0]]['name']!='AIR':
                                ldelete.append(k)       

        #print(ldelete)
        # pdb.set_trace()
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
        uout = np.where([not self.Gt.node[i]['indoor'] for i in cy])
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



    def outputGi(self,verbose=False,tqdmpos=0.):
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
                # list all potential successors of interaction i1
                i2 = nx.neighbors(self.Gi, i1)
                # create a Cone object
                cn = cone.Cone()
                # if starting from segment
                if nstr0 > 0:
                    pseg0 = self.seg2pts(nstr0).reshape(2, 2).T
                    # if nstr0 and nstr1 are connected segments
                    if (len(np.intersect1d(nx.neighbors(self.Gs, nstr0), nx.neighbors(self.Gs, nstr1))) == 0):
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

                ipoints = [x for x in i2 if len(x)==1 ]                               # i0      i1     i2[x]  
                # Avoid to have the same diffaction point after reflection exemple :  (-10,),(245,12),(-10,) impossible 
                #                                                                      nstr0  nstr1 
                if nstr0<0: 
                    ipoints = [x for x in ipoints if x[0]!=nstr0] 
                #ipoints = filter(lambda x: len(x) == 1, i2)
                pipoints = np.array([self.Gs.pos[ip[0]] for ip in ipoints]).T
                # filter tuple (R | T)
                #istup = filter(lambda x : type(eval(x))==tuple,i2)
                # map first argument segment number
                #isegments = np.unique(map(lambda x : eval(x)[0],istup))
                isegments = np.unique(
                    filter(lambda y: y > 0, map(lambda x: x[0], i2)))
                # if nstr0 and nstr1 are adjescent segment remove nstr0 from
                # potential next interaction
                # Fix 01/2017
                # This is not always True if the angle between 
                # the two adjascent segments is < pi/2
                nb_nstr0 = self.Gs.neighbors(nstr0)
                nb_nstr1 = self.Gs.neighbors(nstr1)
                common_point = np.intersect1d(nb_nstr0,nb_nstr1)
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
                    if len(i1) == 3:
                        #if ((e[0]==(53,17)) and (e[1]==(108,17,18))):
                        #    typ, prob = cn.belong_seg(pta, phe,visu=True)
                        #else:
                        typ, prob = cn.belong_seg(pta, phe)
                        # if bs.any():
                        #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                        # if ~bs.any():
                        #    plu.displot(pta[:,~bs],phe[:,~bs],color='k')

                    # i1 : interaction R --> mirror
                    if len(i1) == 2:
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

                output = nx.neighbors(self.Gi, (nstr1,))
                nout = len(output)
                probint = np.ones(nout)  # temporarybns
                dintprob = {k: v for k, v in zip(output, probint)}

            self.Gi.add_edge(i0, i1, output=dintprob)


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
        #     #ipoints = filter(lambda x: len(x) == 1, i2)
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
        #         #    filter(lambda x: x != nstr0, isegments))
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
        #         # output = filter(lambda x: x[0] in isegkeep, i2)
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
        lint = self.Gi.node

        # list of tuple interactions (R|T)
        lD = [x for x in lint if len(x)==1]
        lR = [x for x in lint if len(x)==2]
        lT = [x for x in lint if len(x)==3]
        # lD = filter(lambda x: len(x) == 1, lint)
        # lR = filter(lambda x: len(x) == 2, lint)
        # lT = filter(lambda x: len(x) == 3, lint)

        # visible R|T source cycle is ncy

        lR = filter(lambda x: x[1] == ncy, lR)
        if typ == 'source':
            lT = filter(lambda x: x[1] == ncy, lT)
        if typ == 'target':
            lT = filter(lambda x: x[2] == ncy, lT)
        if typ == 'all':
            lT = lT
        # Finding the diffraction points
        # Diffraction points are different from indoor cycle and outdoor
        # cycles
        #
        # TODO check wedge validity.
        #

        vnodes = self.Gt.node[ncy]['polyg'].vnodes
        vpoints = filter(lambda x: x < 0, vnodes)
        lD = []
        for x in vpoints:
            if self.ddiff.has_key(x):
                for y in self.ddiff[x][0]:
                    if y == ncy:
                        lD.append((x,))
        # indoor = self.Gt.node[ncy]['indoor']
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
        #
        segfilt = filter(lambda x: x not in lair, self.tsg)
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

        fig, ax = self.showG('s', nodes=False, edgelist=edfilt)

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

        #     if k==1:
        #         fig,ax = self.showG('s',fig=fig,ax=ax,nodelist=ldeg,edges=False,nodes=True,node_size=50,node_color='c')
        #     if k==4:
        #         fig,ax = self.showG('s',fig=fig,ax=ax,nodelist=ldeg,nodes=False,node_size=50,node_color='b')

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
            False
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
                    'slab': False,
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
            >>> L = Layout('TA-Office.ini')
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
                    'slab': False,
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
                    'figsize': (5, 5),
                    'mode': 'nocycle',
                    'alphacy': 0.8,
                    'colorcy': '#abcdef',
                    'lvis': ['nn', 'ne', 'ee'],
                    'linter': ['RR', 'TT', 'RT', 'TR', 'RD', 'DR', 'TD', 'DT', 'DD'],
                    'show0': False,
                    'axis': False,
                    'overlay': False,
                    'diffraction': False
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value
        if kwargs['aw'] != []:
            kwargs['airwalls'] = kwargs['aw']
        # overriding first argument graph
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

            # lss = filter(lambda x: self.Gs.node[x].has_key('ss_name'),self.Gs.nodes())
            # lss = filter(lambda x: len(self.Gs.node[x]['ss_name'])>0,lss)

            # keep track of segments already printed

            nodelistbkup = kwargs['nodelist']
            edgelistbkup = kwargs['edgelist']
            widthbkup = kwargs['width']
            nodecolbkup = kwargs['edge_color']

            try:
                sllist = [kwargs['sllist'].pop()]
            except:
                sllist = self.name.keys()

            for lmat in sllist:
                lseg = self.name[lmat]
                if lseg != []:
                    lseg2 = [np.where(np.array(self.Gs.edges()) == i)[0]
                             for i in lseg]
                    kwargs['edgelist'] = reduce(
                        lambda x, y: list(x) + list(y), lseg2)
                    if kwargs['slab']:
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
                for ss in d.keys():
                    color = cold[self.sl[ss]['color']]
                    for ns in d[ss]:
                        norm = self.Gs.node[ns[0]]['norm']
                        np1, np2 = self.Gs.neighbors(ns[0])
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
                nodes = G.nodes()
                edges = G.edges()
                nodf = filter(lambda x: x != 0, nodes)
                edf = filter(lambda x: ((edges[x][0] != 0) & (
                    edges[x][1] != 0)), np.arange(len(edges)))
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

            edges = G.edges()
            rle = range(len(edges))
            eded = filter(lambda x: (edges[x][0] > 0) & (edges[x][1] > 0), rle)
            ndnd = filter(lambda x: (edges[x][0] < 0) & (edges[x][1] < 0), rle)
            nded = filter(lambda x: (((edges[x][0] < 0) & (edges[x][1] > 0)) |
                                     ((edges[x][0] > 0) & (edges[x][1] < 0))), rle)
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

            edges = G.edges()

            # range len edges

            rle = range(len(edges))

            DD = filter(lambda x:  ((len(edges[x][0]) == 1) &
                                    (len(edges[x][1]) == 1)), rle)

            RR = filter(lambda x: ((len(edges[x][0]) == 2) &
                                   (len(edges[x][1]) == 2)), rle)

            TT = filter(lambda x: ((len(edges[x][0]) == 3) &
                                   (len(edges[x][1]) == 3)), rle)

            RT = filter(lambda x: ((len(edges[x][0]) == 2) &
                                   (len(edges[x][1]) == 3)), rle)

            TR = filter(lambda x: ((len(edges[x][0]) == 3) &
                                   (len(edges[x][1]) == 2)), rle)

            RD = filter(lambda x:  ((len(edges[x][0]) == 2) &
                                    (len(edges[x][1]) == 1)), rle)

            TD = filter(lambda x:  ((len(edges[x][0]) == 3) &
                                    (len(edges[x][1]) == 1)), rle)

            DR = filter(lambda x:  ((len(edges[x][0]) == 1) &
                                    (len(edges[x][1]) == 2)), rle)

            DT = filter(lambda x:  ((len(edges[x][0]) == 1) &
                                    (len(edges[x][1]) == 3)), rle)

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
                    ndlist = map(lambda x: edges[x][0], kwargs['edgelist']) +\
                        map(lambda x: edges[x][1], kwargs['edgelist'])
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
                for k, ncy in enumerate(self.Gt.node.keys()):
                    if k != 0:
                        fig, ax = self.Gt.node[ncy]['polyg'].plot(
                            alpha=kwargs['alphacy'], color=kwargs['colorcy'], **args)
                        args['fig'] = fig
                        args['ax'] = ax
            if kwargs['mode'] == 'room':
                for k, nro in enumerate(self.Gr.node.keys()):
                    if k != 0:
                        fig, ax = self.Gr.node[nro]['cycle'].show(**args)
                        args['fig'] = fig
                        args['ax'] = ax

        kwargs['ax'].axis('scaled')
        if not kwargs['axis']:
            kwargs['ax'].axis('off')

        if kwargs['overlay']:
            imok = False
            if self.display['overlay_file'] != '':
                image = Image.open(os.path.join(
                    pro.basename, pro.pstruc['DIRIMAGE'], self.display['overlay_file']))
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
            >>> L = Layout('WHERE1.ini')
            >>> L.dumpr()

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
        >>> L = Layout('DLR.ini')
        >>> walls = L.thwall(0,0)

        """
        keyn = self.Gs.node.keys()
        walls = []
        for nd in keyn:
            if nd > 0:
                nb = self.Gs.neighbors(nd)
                pta = self.Gs.pos[nb[0]]
                phe = self.Gs.pos[nb[1]]
                pn = self.Gs.node[nd]['norm']
                name = self.Gs.node[nd]['name']
                transition = self.Gs.node[nd]['transition']
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
            return self.Gs.node[n]['ncycles']
        else:
            nseg = self.Gs[n].keys()
            cy = []
            for nn in nseg:
                cy.extend(self.ptGs2cy(nn))
            cy = np.unique(cy).tolist()
            return cy

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

        See Also
        --------

        Layout.cy2pt

        """

        ptsh = sh.Point(pt[0], pt[1])
        cycle_exists = False

        for ncy in self.Gt.node.keys():
            if ncy > 0:
                criter1 = self.Gt.node[ncy]['polyg'].touches(ptsh)
                criter2 = self.Gt.node[ncy]['polyg'].contains(ptsh)
                if (criter1 or criter2):
                    cycle_exists = True
                    return(ncy)
        if not cycle_exists:
            raise NameError(str(pt) + " is not in any cycle")

    def cy2pt(self, cy=0, h=1.2):
        """return a point into a given cycle

        Parameters
        ----------

        cy : int
            cycle number

        h : float
            point height

        Returns
        -------

        point  : nd.array

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
        for nr in self.Gr.node.keys():
            if self.Gr.node[nr]['polyg'].contains(ptsh)\
                    or self.Gr.node[nr]['polyg'].touches(ptsh):
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
        for nr in self.Gr.node.keys():
            # if seg in self.Gt.node[self.Gr.node[nr]['cycle']]['vnodes']:
            ncy = self.Gr.node[nr]['cycle']
            if seg in self.Gt.node[ncy]['cycle'].cycle:
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
            ncy = self.Gr.node[room]['cycle']
            seg = self.Gt.node[ncy].cycle
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
            ncy = self.Gr.node[room]['cycle']
            nod = self.Gt.node[ncy].cycle
            #nod = self.Gt.node[self.Gr.node[room]['cycle']]['vnodes']
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

            Info
            ---- 
            As a diffraction point may involve iso segments the nature 
            of the diffraction interaction depends on a height parameter
            This function extact the couple of slab from this information

            Returns
            -------
            - a list of 2-segments list. the length of this list == length of lz
            - a list of slab tuples.  the length of this list == length of lz

            [[443, 529], [444, 530]]
            [['WALL', 'WALL'], ['AIR', 'AIR']]

        """
        assert(npt in self.ddiff), logging.error('npt not a diffraction point')
        lcy = self.ddiff[npt][0]
        ls = []
        llz = len(lz)
        dz_seg= {z:[] for z in range(llz)}
        dz_sl= {z:[] for z in range(llz)}

        for cy in lcy: 
            vn = set(self.Gt.node[cy]['polyg'].vnodes)   
            lneig_pt = set(nx.neighbors(self.Gs,npt))
            lseg = lneig_pt.intersection(vn)
            lseg_valid = [ x for x in lseg if self.Gs.node[x]['name']!='_AIR']

            for x in lseg_valid:
                zsup = lz >self.Gs.node[x]['z'][0]
                zinf = lz <=self.Gs.node[x]['z'][1]
                z    = zsup & zinf 
                uz = np.where(z)[0]
                # fill dz_seg at the correct height with a lseg_valid 
                # and simulnaneously 
                # fill dz_sl at the correct height with correspondong slab
                [(dz_seg[i].append(x),dz_sl[i].append(self.Gs.node[x]['name']))
                                                                    for i in uz]

        return dz_seg.values(),dz_sl.values()

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

        dangles = {cy: np.array(geu.get_pol_angles(self.Gt.node[cy]['polyg']))
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
            lcyk = self.Gs.node[k]['ncycles']
            if len(lcyk) > 2:
                # Subgraph of connected cycles around k
                Gtk = nx.subgraph(self.Gt, lcyk)
                # ordered list of connections between cycles
                try:
                    lccyk = nx.find_cycle(Gtk)
                except:
                    pdb.set_trace()

                # list of segment neighbours
                neigh = self.Gs[k].keys()
                # sega : list of air segment in neighors
                sega = [n for n in neigh if
                        (self.Gs.node[n]['name'] == 'AIR' or
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

    # def buildGr(self):
    #     """ build the graph of rooms Gr

    #     Notes
    #     -----

    #     + A room is a set of convex cycles connected together through _AIR or AIR segment
    #     + A room is not necessarily convex

    #     This graph is used for room electromagnetics evaluation

    #     """

    #     # list of all indoor convex cycle
    #     to_visit = [ x for x in self.Gt.node.keys() if self.Gt.node[x]['indoor']]
    #     # list of cycles already in a room
    #     visited = []
    #     # law : list of air walls
    #     law = self.name['AIR']+self.name['_AIR']
    #     room_cnt = 1
    #     r = {}
    #     cur_cy = to_visit[0]
    #     while len(to_visit)>0:
    #         cur_cy = to_visit.pop()
    #         #printcur_cy
    #         # get neighbors of current_cycle
    #         neighbors = nx.neighbors(self.Gt,cur_cy)
    #         # get neighbors separated from the current cycle by an air_wall
    #         neighbors_aw = [ x for x in neighbors if self.Gt[cur_cy][x]['segment'] in law ]
    #         # if at least one of the neighbors is already in a room --> same room
    #         # else --> new room
    #         neighbors_already_in_a_room = [ x for x in neighbors if x in visited ]
    #         for cy in neighbors:
    #             if len(already_in_a_room) >0 :
    #                 r[cy] = r[already_in_a_room[0]]
    #             else:
    #                 r[cy] = room_cnt
    #                 room_cnt += 1
    #             visited.append(cy)
    #         cur_cy = neighbors_aw[0]
    #     return r,to_visit,visited

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

        to_visit = [x for x in self.Gr.node.keys() if self.Gr.node[
            x]['indoor']]

        law = self.name['AIR'] + self.name['_AIR']
        while len(to_visit) > 0:
            remaining_nodes = self.Gr.node.keys()
            to_visit = [x for x in to_visit if x in remaining_nodes]
            cur_cy = to_visit.pop()
            ln = nx.neighbors(self.Gr, cur_cy)
            ln = [x for x in ln if self.Gr.node[x]['indoor']]
            # ls : list of segment
            ls = [self.Gr[cur_cy][x]['segment'] for x in ln]
            for x in zip(ln, ls):
                if x[1] in law:
                    p1 = self.Gr.node[cur_cy]['polyg']
                    p2 = self.Gr.node[x[0]]['polyg']
                    P = p1 + p2
                    P.setvnodes(self)
                    ine = nx.neighbors(self.Gr, x[0])
                    sne = [self.Gr[x[0]][y]['segment'] for y in ine]
                    self.Gr.add_node(cur_cy, polyg=P)
                    self.Gt.pos[cur_cy] = np.array((P.centroid.xy)).squeeze()
                    self.Gr.add_edges_from([(cur_cy, y)
                                            for y in ine if y != cur_cy])
                    for k in zip(ine, sne):
                        if k[0] != cur_cy:
                            self.Gr[cur_cy][k[0]]['segment'] = k[1]
                    self.Gr.remove_node(x[0])
                    del(self.Gr.pos[x[0]])

        #
        #  Connected components might not be all contiguous
        #  this is a problem because the concatenation of cycles
        #  operation requires cycles contiguity
        #
        # for n in self.Gr.nodes():
        #     self.Gr.node[n]['transition'] = []
        # ltrans = self.listtransition
        # ltmp = filter(lambda x:self.Gs.node[x]['name']!='AIR',ltrans)
        # ldoors = filter(lambda x:self.Gs.node[x]['name']!='_AIR',ltmp)
        # keys = self.Gr.node.keys()

        # for a in _airseg:
        #     n0,n1=iuE[a]
        #     found=False
        #     while not found:
        #         nn0 = mapoldcy[n0]
        #         if n0==nn0:
        #             found=True
        #         else:
        #             n0=nn0
        #     found=False
        #     while not found:
        #         nn1 = mapoldcy[n1]
        #         if n1==nn1:
        #             found=True
        #         else:
        #             n1=nn1

        #     p0=self.Gt.node[n0]['polyg']
        #     p1=self.Gt.node[n1]['polyg']

        #     # Merge polygon
        #     P = p0+p1
        #     # If the new Polygon is convex update Gt
        #     #
        #     if geu.isconvex(P):
        #         # updates vnodes of the new merged polygon
        #         P.setvnodes(self)
        #         # update edge
        #         n0s = n0
        #         n1s = n1
        #         # get segments information from cycle n0
        #         dne = self.Gt[n0]
        #         # remove connection to n0 to avoid a cycle being
        #         # connected to itself
        #         self.Gt[n1].pop(n0)
        #         # add information from adjacent cycle n1
        #         dne.update(self.Gt[n1])
        #         # list of items of the merged dictionnary
        #         ine = dne.items()
        #         #update n0 with the new merged polygon
        #         self.Gt.add_node(n0,polyg=P)
        #         # connect new cycle n0 to neighbors
        #         self.Gt.add_edges_from([(n0,x[0],x[1]) for x in ine if x[0] != n0])
        #         #remove old cycle n1
        #         self.Gt.remove_node(n1)
        #         # update pos of the cycle with merged polygon centroid
        #         self.Gt.pos[n0] = np.array((P.centroid.xy)).squeeze()
        #         # delete _air segment a
        #         # do not apply g2npy
        #         self.del_segment(a,verbose=False,g2npy=False)
        #         mapoldcy[n1]=n0
        # for cy in keys:
        #     p = self.Gr.node[cy]['polyg']
        #     lseg = [ x for x in p.vnodes if x > 0 ]
        #     #hasdoor = filter(lambda n : n in ldoors,lseg)
        #     hasdoor = [ x for x in lseg if x in ldoors]
        #     if len(hasdoor)>0:
        #         pass
        #     else:
        #         self.Gr.remove_node(cy)
        #         self.Gr.pos.pop(cy)

        # # Destroy edges which do not share a transition
        # for e in self.Gr.edges():
        #     keep = False
        #     if (e[0]>0) & (e[1]>0):
        #         cy1 = self.Gr.node[e[0]]['cycle']
        #         cy2 = self.Gr.node[e[1]]['cycle']
        #         f,b = cy1.intersect(cy2)
        #         for s in b:
        #             if s>0:
        #                 if self.Gs.node[s]['transition']:
        #                     keep = True
        #                     self.Gr.node[e[0]]['transition'].append(s)
        #                     self.Gr.node[e[1]]['transition'].append(s)

        #     if not keep:
        #         self.Gr.remove_edge(*e)

    def waypoint(self, nroom1, nroom2):
        """ get the waypoint between room1 and room2

        Parameters
        -----------

        nroom1 : int 
            room number 1
        nromm2 : int 
            room number 2

        Returns
        -------

        waypoint : 

        """
        rooms = nx.dijkstra_path(self.Gr, nroom1, nroom2)
        nroom = len(rooms)
        waypoint = []
        for k in np.arange(nroom - 1):
            room1 = rooms[k]
            proom1 = self.Gr.pos[room1]
            room2 = rooms[k + 1]
            doors1 = self.Gr.node[room1]['doors']
            doors2 = self.Gr.node[room2]['doors']
            doorId = np.intersect1d(doors1, doors2)[0]
            #
            # coord door
            #
            unode = self.Gs.neighbors(doorId)
            p1 = self.Gs.pos[unode[0]]
            p2 = self.Gs.pos[unode[1]]
            pdoor = (np.array(p1) + np.array(p2)) / 2
            waypoint.append((proom1[0], proom1[1]))
            waypoint.append((pdoor[0], pdoor[1]))

        proom2 = self.Gr.pos[nroom2]
        waypoint.append((proom2[0], proom2[1]))
        return(waypoint)

    def editor(self):
        """ invoke interactive layout graphical editor

        Notes
        -----

        point edition

            m  toggle point edition mode  (CP : Create Point)

                lclic same x
                rclic same y
                cclic free point

        segment edition

            [0-f] - display one of the 16 first layers
            x : save structure
            o : toggle overlay

        """

        fig = plt.gcf()
        ax = fig.add_subplot(111)
        self.display['nodes'] = True
        self.display['ednodes'] = True
        self.display['subsegnb'] = True
        self.display['transition'] = True
        self.display['ticksoff'] = True

        #self.af = SelectL2(self,fig=fig,ax=ax)
        self.af = SelectL(self, fig=fig, ax=ax)

        fig, ax = self.af.show(fig, ax, clear=True)

        self.cid1 = fig.canvas.mpl_connect('button_press_event',
                                           self.af.OnClick)
        self.cid2 = fig.canvas.mpl_connect('button_release_event',
                                           self.af.OnClickRelease)
        self.cid3 = fig.canvas.mpl_connect('motion_notify_event',
                                           self.af.OnMotion)
        self.cid4 = fig.canvas.mpl_connect('key_press_event',
                                           self.af.OnPress)
        self.cid5 = fig.canvas.mpl_connect('key_release_event',
                                           self.af.OnRelease)

        plt.draw()
        plt.axis('tight')
        plt.show()

#        """
#        """
#
#        # import gtk
#        from matplotlib.figure import Figure
#        from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#        from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
#        from matplotlib.backend_bases import key_press_handler
#
#        win = gtk.Window()
#        win.show_all()
#
#
#        win = gtk.Window()
#        win.connect("destroy", lambda x: gtk.main_quit())
#        win.set_default_size(400,300)
#        win.set_title("Embedding in GTK")
#
#        vbox = gtk.VBox()
#        win.add(vbox)
#
#        fig = Figure()
#        ax = fig.add_subplot(111)
#
#        fig,ax = self.showG('s',fig=fig,ax=ax)
#
#        canvas = FigureCanvas(fig)  # a gtk.DrawingArea
#        canvas.show()
#        vbox.pack_start(canvas)
#        toolbar = NavigationToolbar(canvas, win)
#        vbox.pack_start(toolbar, False, False)
#
#
#        def on_key_event(event):
#            print('you pressed %s'%event.key)
#            key_press_handler(event, canvas, toolbar)
#
#        canvas.mpl_connect('key_press_event', on_key_event)
#
#        win.show_all()
#        gtk.main()

    def editorTk(self):
        """ invoke interactive layout graphical editor

        Notes
        -----

        point edition

            m  toggle point edition mode  (CP : Create Point)

                lclic same x
                rclic same y
                cclic free point

        segment edition

            [0-f] - display one of the 16 first layers
            x : save structure
            o : toggle overlay

        """

        #import matplotlib
        # matplotlib.use('TkAgg')

        from matplotlib.backend_bases import key_press_handler
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
        from matplotlib.figure import Figure
        import Tkinter as Tk

        root = Tk.Tk()
        root.wm_title('Pylayers Layout Editor')

        fig = Figure()
        ax = fig.add_subplot(111)
        # ax.plot(np.arange(10))

        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        button = Tk.Button(master=root, text='Quit', command=sys.exit)
        button.pack(side=Tk.BOTTOM)

        self.display['nodes'] = True
        self.display['ednodes'] = True

        select = SelectL(self, canvas)

        # self.af.show(clear=True)

        self.cid1 = canvas.mpl_connect('button_press_event', select.OnClick)
        self.cid2 = canvas.mpl_connect('key_press_event', select.OnPress)
        # ax.axis('tight')
        canvas.show()
        Tk.mainloop()

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
        """
        facets3d(edlist,name)
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

        """
        nta = np.array(nx.neighbors(self.Gs, ta))
        nhe = np.array(nx.neighbors(self.Gs, he))
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

        """
        # transpose point numbering

        upnt = filter(lambda x: x < 0, self.Gs.nodes())
        ta = np.nonzero(np.array(upnt) == ta)[0][0]
        he = np.nonzero(np.array(upnt) == he)[0][0]
        res = filter(lambda x: (((x[0] == ta) & (x[1] == he))
                                | ((x[0] == he) & (x[1] == ta))), zip(self.tahe[0], self.tahe[1]))
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
        pts = np.array(self.Gs.pos.values()).T
        ke = np.array(self.Gs.pos.keys())
        u = pts - pt.reshape(2, 1)
        v = np.sqrt(np.sum(u * u, axis=0))
        nz = (v > tol)
        b = nz.prod()
        if b == 1:
            # if all layout points are different from pt
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
        ke = np.array(self.Gs.pos.keys())        # point keys
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
        nebr = self.Gs.neighbors(s)
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
        >>> L = Layout('DLR.ini')
        >>> pg = L.geomfile()

        """

        # calculate center of gravity
        if centered:
            pg = np.sum(self.pt, axis=1) / np.shape(self.pt)[1]
        else:
            pg = np.array([0, 0])

        # en  = self.Ns # number of segments
        en = len(np.where(np.array(self.Gs.node.keys()) > 0)[0])
        if en != self.Ns:
            logging.warning(
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
            logging.warning(
                "wrong number of subsegment consistency problem in layout")

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
        for i in self.Gs.node.keys():
            if i > 0:  # segment
                if ((self.Gs.node[i]['name'] != 'AIR') and
                        (self.Gs.node[i]['name'] != '_AIR')):
                    nebr = self.Gs.neighbors(i)
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
                nebr = self.Gs.neighbors(l[0])
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
        en = len(np.where(np.array(self.Gs.node.keys()) > 0)[0])
        if en != self.Ns:
            logging.warning(
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
            logging.warning(
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
        for i in self.Gs.node.keys():
            if i > 0:  # segment
                if ((self.Gs.node[i]['name'] != 'AIR') and
                        (self.Gs.node[i]['name'] != '_AIR')):
                    nebr = self.Gs.neighbors(i)
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
                nebr = self.Gs.neighbors(l[0])
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

        boxes = np.empty((npt / 4, 4), dtype='int')
        b = np.arange(npt / 4)

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
                    ptc = np.hstack((ptc, pt))
                    boxc = np.vstack((boxc, box + cpt))
                    cpt = cpt + nbpt

                # manage Ceil color
                colname = sl['CEIL']['color']
                colhex = cold[colname]
                colf = np.repeat((pyu.rgb(colhex))[np.newaxis, :], 4, axis=0)
                color = np.vstack((color, colf))

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
        """ plot the layout with shapely polygons

        Parameters
        ---------

        fig 
        ax 
        labels : list
        nodes : boolean

        """
        defaults = {'show': False,
                    'fig': [],
                    'ax': [],
                    'labels': [],
                    'nodes': False
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] == []:
            fig = plt.gcf()
        if kwargs['ax'] == []:
            ax = plt.gca()

        if isinstance(kwargs['labels'], list):
            labels = kwargs['labels']
        elif kwargs['labels'] == True:
            labels = ['s', 't', 'v', 'i', 'w']
        elif isinstance(kwargs['labels'], str):
            labels = kwargs['labels']
        else:
            labels = []

        k = self.Gs.pos.keys()
        v = self.Gs.pos.values()

        kk = np.array(k)
        vv = np.array(v)

        w = [str(x) for x in kk]

        if 's' in labels:
            [ax.text(vv[i, 0], vv[i, 1], w[i]) for i in range(len(w))]

        if kwargs['nodes']:
            ax.scatter(vv[:, 0], vv[:, 1])
        ML = sh.MultiLineString(self._shseg.values())

        self.pltlines(ML, color='k', fig=fig, ax=ax)

        return fig, ax

    def get_Sg_pos(self, sigarr):
        """ return position of the signatures

        Parameters
        ----------

        sigarr : signature
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

        nth = np.array(map(lambda n: nx.neighbors(self.Gs, n), lns))
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
        >>> L = Layout('defstr.ini')
        >>> p_Tx,p_Rx = L.randTxRx()

        Notes
        -----
            ex fn Tx_Rx_pos


        """

        # self.boundary()

        Tx_x = rd.uniform(self.ax[0], self.ax[1])
        Tx_y = rd.uniform(self.ax[2], self.ax[3])
        Rx_x = rd.uniform(self.ax[0], self.ax[1])
        Rx_y = rd.uniform(self.ax[2], self.ax[3])

        p_Tx = np.array([Tx_x, Tx_y])
        p_Rx = np.array([Rx_x, Rx_y])

        return(p_Tx, p_Rx)

    def boundary(self, percx=0.15, percy=0.15, xlim=(), force=False):
        """ add a blank boundary around layout

        Parameters
        ----------

        percx : float
            percentage of Dx for x offset calculation (default 0.15)
        percy : float
           percentage of Dy for y offset calculation (default 0.15)

        self.lboundary is the list of the nodes of the added boundary
        self.axn is the zone without the boundary extension
        self.ax  is updated

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('defstr.ini')
        >>> L.boundary()

        """
        if not self.hasboundary or force:
            if len(self.Gs.pos.values()) != 0:
                xmax = max(p[0] for p in self.Gs.pos.values())
                xmin = min(p[0] for p in self.Gs.pos.values())
                ymax = max(p[1] for p in self.Gs.pos.values())
                ymin = min(p[1] for p in self.Gs.pos.values())
            elif xlim == ():
                xmin = -20.
                xmax = 20.
                ymin = -10.
                ymax = 10.
            else:
                xmin = xlim[0]
                xmax = xlim[1]
                ymin = xlim[2]
                ymax = xlim[3]

            Dx = xmax - xmin
            Dy = ymax - ymin
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
