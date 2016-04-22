#-*- coding:Utf-8 -*-
"""
.. currentmodule:: pylayers.gis.layout

This class handle the description of an Indoor layout

Class Layout
============

.. autosummary::
    :toctree: generated

Utility functions
-----------------

.. autosummary::
    :toctree: generated/

    Layout.__init__
    Layout.__repr__
    Layout.__add__
    Layout.__mul__
    Layout.info
    Layout.ls
    Layout.delete
    Layout.check
    Layout.clip
    Layout.help
    Layout.g2npy
    Layout.check2
    Layout.cleanup
    Layout.boundary
    Layout.distwall
    Layout.randTxRx
    Layout.get_paths
    Layout.get_zone
    Layout.info_segment
    Layout.find_edgelist

Loading and Saving
------------------

.. autosummary::
    :toctree: generated/

    Layout.dumpw
    Layout.dumpr
    Layout.loadosm
    Layout.saveosm
    Layout.saveini
    Layout.loadini
    Layout.loadfur
    Layout.load
    Layout.loadstr
    Layout.savestr2
    Layout.save
    Layout.saveold

Layout editing
--------------

.. autosummary::
    :toctree: generated/

    Layout.editor
    Layout.editorGtk
    Layout.editorTk
    Layout.add_pnod
    Layout.add_fnod
    Layout.add_nfpe
    Layout.add_pons
    Layout.add_segment
    Layout.add_furniture
    Layout.add_furniture_file
    Layout.del_points
    Layout.del_segment
    Layout.thwall
    Layout.edit_point
    Layout.edit_segment
    Layout.add_window
    Layout.add_door

Layout transformation
---------------------

.. autosummary::
    :toctree: generated/

    Layout.translate
    Layout.rotate
    Layout.del_cycle
    Layout.chgmss
    Layout.diag
    Layout.nd2seg
    Layout.ed2nd
    Layout.seguv
    Layout.segpt2
    Layout.seg2pts
    Layout.segpt
    Layout.extrseg
    Layout.seginframe2
    Layout.seginframe
    Layout.layerongrid
    Layout.cycleinline
    Layout.seginline

Layout visibility
-----------------

.. autosummary::
    :toctree: generated/

    Layout.visilist
    Layout.closest_edge
    Layout.visi_papb
    Layout.angleonlink
    Layout.angleonlinkold
    Layout.layeronlink

Layout interactions
-------------------

    Layout.intercy

SubSegment Functions
--------------------

.. autosummary::
    :toctree: generated/

    Layout.subseg
    Layout.have_subseg
    Layout.del_subseg
    Layout.add_subseg

Vizualisation
--------------

.. autosummary::
    :toctree: generated/

    Layout.displaygui
    Layout.show_nodes
    Layout.show_seg1
    Layout.plot_segments
    Layout.show_segment
    Layout.show_layer
    Layout.facet3D
    Layout.facets3D
    Layout.geomfile
    Layout._show3
    Layout.show3

Showing Graphs
---------------

.. autosummary::
    :toctree: generated/

    Layout._showGi
    Layout._showGt
    Layout.showGs
    Layout.show
    Layout.showG
    Layout._showGv

Building Graphs
---------------

.. autosummary::
    :toctree: generated/

    Layout.build
    Layout.buildGt
    Layout.buildGw
    Layout.buildGv
    Layout.buildGi
    Layout.outputGi
    Layout.waypointGw
    Layout.builGr2
    Layout.buildGr
    Layout.buildGr3


Cycles and Rooms  Related Functions
-----------------------------------

.. autosummary::
    :toctree: generated/

    Layout.pt2cy
    Layout.cy2pt
    Layout.pt2ro
    Layout.seg2ro

    Layout.room2segments
    Layout.room2nodes

    Layout.waypoint

Testing Functions
------------------

.. autosummary::
    :toctree: generated/

    Layout.isseg
    Layout.numseg
    Layout.ispoint
    Layout.onseg


Signatures
----------

.. autosummary::
    :toctree: generated/

    Layout.signature
    Layout.showSig
    Layout.get_Sg_pos

"""
#
#
import pdb
import os
import copy
import glob
import pickle
import cPickle
import ConfigParser
import numpy as np
import numpy.random as rd
import scipy as sp
import struct as stru
from   scipy import io
import doctest
import random
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import networkx as nx
from networkx.readwrite import write_gpickle,read_gpickle
import shapely.geometry as sh
import shapely.ops as sho
from shapely.ops import cascaded_union
from descartes.patch import PolygonPatch
from numpy import array
import PIL.Image as Image
import logging
import urllib2 as urllib
from cStringIO import StringIO

from pylayers.antprop import slab as sb
from pylayers.util import geomutil as geu
from pylayers.util import plotutil as plu
from pylayers.util import pyutil as pyu
from pylayers.util import graphutil as gru
from pylayers.util import cone

#from  more_itertools import unique_everseen

# Handle furnitures
import pylayers.gis.furniture as fur
import pylayers.gis.osmparser as osm
#from pylayers.gis import cycles as Cycls
from pylayers.gis import cycles as cycl # new version of cycles
from pylayers.gis.selectl import SelectL
from pylayers.gui.editor_select import SelectL2

from pylayers.util.easygui import *
from pylayers.util.project import *
#from   PyUtil  import *
#from   interval import interval
from itertools import combinations
import pdb
import ast
import pylayers.util.graphutil as gph
from mpl_toolkits.basemap import Basemap
try:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi import mlab
except:
    print 'Layout:Mayavi is not installed'


class Layout(PyLayers):
    """ Handling Layout

    Attributes
    ----------

    Gs     : Garph of points and segment (structure)
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

    This class exploits `networkx` to store Layout information

    .. autosummary::

    """
    def __init__(self,_filename='defstr.ini',_filematini='matDB.ini',_fileslabini='slabDB.ini',_filefur='',force=False,check=True):
        """ object constructor

        Parameters
        ----------

        _filename : string
            layout file name
        _filematini :
            material dB file name
        _fileslabini :
            slab dB file name
        _filefur :
            furniture file name
        force : boolean
        check : boolean

        """

        mat = sb.MatDB()
        mat.load(_filematini)

        self.sl = sb.SlabDB()
        self.sl.mat = mat
        self.sl.load(_fileslabini)
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
        # related file names
        #

        self.filename = _filename
        self.fileslabini = _fileslabini
        self.filematini = _filematini
        self.filefur = _filefur


        self.hasboundary=False
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
        self.display['activelayer'] = self.sl.keys()[0]
        self.display['layers'] = []
        self.display['overlay'] = False
        self.display['overlay_flip'] = ""
        #self.display['overlay_file']="/home/buguen/Pyproject/data/image/"
        self.display['overlay_file'] = ""
        self.display['overlay_axis'] = ""
        self.display['layerset'] = self.sl.keys()
        self.display['box'] = (-50,50,-50,50)
        self.name = {}

        self.zmin = 0

        for k in self.sl.keys():
            self.name[k] = []

        self.load(_filename)
        # check layout integrity (default)
        if check:
            self.check()
        
        if not self.hasboundary:
            self.boundary(dx=10,dy=10)

        # If the layout has already been built then load the built structure
        if not force:
            try:
                self.dumpr()
            except:
                pass

    def __repr__(self):
        st = '\n'
        st = st + "----------------\n"
        st = st + self.filename + "\n"
        if self.display['overlay_file']!='':
            filename = pyu.getlong(self.display['overlay_file'],os.path.join('struc','images'))
            st = st + "Image('"+filename+"')\n"
        st = st + "----------------\n\n"
        st = st + "Number of points  : "+ str(self.Np)+"\n"
        st = st + "Number of segments  : "+str(self.Ns)+"\n"
        st = st + "Number of sub segments  : "+str(self.Nss)+"\n"
        st = st + "Number of cycles  : "+ str(len(self.Gt.node))+"\n"
        st = st + "Number of rooms  : "+ str(len(self.Gr.node))+"\n"
        if hasattr(self,'degree'):
            for k in self.degree:
                if  (k < 2) or (k>3):
                    st = st + 'degree '+str(k)+' : '+str(self.degree[k])+"\n"
                else:
                    st = st + 'number of node point of degree '+str(k)+' : '+str(len(self.degree[k]))+"\n"
        st = st + "\n"
        st = st + "xrange :"+ str(self.ax[0:2])+"\n"
        st = st + "yrange :"+ str(self.ax[2:])+"\n"
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
        if hasattr(self,'iupnt'):
            st = st + "iupnt : get point index in self.pt from point id  "+"\n"
        if hasattr(self,'lsss'):
            st = st + "lsss : list of segments with sub-segment"+"\n"
        if hasattr(self,'sridess'):
            st = st + "stridess : stride to calculate the index of a subsegment" +"\n"
        if hasattr(self,'sla'):
            st = st + "sla : list of all slab names (Nsmax+Nss+1)" +"\n"
        if hasattr(self,'degree'):
            st = st + "degree : degree of nodes " +"\n"

        return(st)

    def __add__(self, other):
        """ addition

        One can add a numpy array or an other layout
        """
        Ls = copy.deepcopy(self)
        if type(other)==np.ndarray:
            for k in Ls.Gs.pos:
                Ls.Gs.pos[k]=Ls.Gs.pos[k]+other[0:2]
        else:
            offp = -min(Ls.Gs.nodes())
            offs = max(Ls.Gs.nodes())
            other.offset_index(offp=offp,offs=offs)
            Ls.Gs.node.update(other.Gs.node)
            Ls.Gs.edge.update(other.Gs.edge)
            Ls.Gs.adj.update(other.Gs.adj)
            Ls.Gs.pos.update(other.Gs.pos)
            Ls.Np = Ls.Np+other.Np
            Ls.Ns = Ls.Ns+other.Ns
            Ls.Nss = Ls.Nss+other.Nss

        return(Ls)


    def __mul__(self,alpha):
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
            assert((type(alpha) == float) or (type(alpha) == int)), " not float"
            alpha=np.array([alpha, alpha, alpha])
        else:
            assert(len(alpha)==3), " not 3D"
        #
        # scaling x & y
        #
        x = np.array(Gs.pos.values())[:,0]
        x = x*alpha[0]

        y = np.array(Gs.pos.values())[:,1]
        y = y*alpha[1]

        xy = np.vstack((x,y)).T
        Ls.Gs.pos = dict(zip(Gs.pos.keys(),tuple(xy)))

        #
        # scaling z
        #

        nseg = filter(lambda x : x>0, Gs.nodes())
        for k in nseg:
            Ls.Gs.node[k]['z'] = tuple((np.array(Ls.Gs.node[k]['z'])-self.zmin)*alpha[2]+self.zmin)
            if Ls.Gs.node[k].has_key('ss_z'):
                Ls.Gs.node[k]['ss_z'] = list((np.array(Ls.Gs.node[k]['ss_z'])-self.zmin)*alpha[2]+self.zmin)

        #
        # updating numpy array from graph
        #

        Ls.g2npy()
        return Ls

    def ls(self, typ='ini'):
        """ list the available file in dirstruc

        Parameters
        ----------

        typ : string optional
            {'str'|'ini'|'osm'|'str2'|'wrl'}

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

        if typ=='str':
            pathname = os.path.join(pstruc['DIRSTRUC'],'*.'+ typ)
        if typ=='str2':
            pathname = os.path.join(pstruc['DIRSTRUC'],'*.' + typ)
        if typ=='ini':
            pathname = os.path.join(pstruc['DIRINI'],'*.' + typ)
        if typ=='osm':
            pathname = os.path.join(pstruc['DIROSM'],'*.' + typ)
        if typ=='wrl':
            pathname = os.path.join(pstruc['DIRWRL'],'*.' + typ)

        lfile_l = glob.glob(os.path.join(basename,pathname))
        lfile_s = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_s.append(fis)
        lfile_s.sort()

        return lfile_s


    def delete(self):
        """ delete Layout graphs

        delete  Gs

        """
        del self.Gs
        self.Gs = nx.Graph()
       
       

    def offset_index(self,offp=0,offs=0):
        """ offset points and segment index

        Parameters
        ----------
        offp : offset points
        offs : offset segments

        """

        newpoint = dict( (k-offp,v) for k,v in self.Gs.node.items() if k <0)
        assert (np.array(newpoint.keys())<0).all()
        newseg =   dict( (k+offs,v) for k,v in self.Gs.node.items() if k > 0)
        assert (np.array(newseg.keys())>0).all()
        newpoint.update(newseg)
        self.Gs.node = newpoint

        newppoint = dict( (k-offp,v) for k,v in self.Gs.pos.items() if k <0)
        newpseg =   dict( (k+offs,v) for k,v in self.Gs.pos.items() if k > 0)
        newppoint.update(newpseg)
        self.Gs.pos =  newppoint

        # adjascence list of segments
        ladjs = [self.Gs.adj[k] for k in self.Gs.adj.keys() if k > 0 ]
        # adjascence list of points
        ladjp = [self.Gs.adj[k] for k in self.Gs.adj.keys() if k < 0 ]

        nladjs = map(lambda x: dict((k-offp,v) for k,v in x.items()),ladjs)
        nladjp = map(lambda x: dict((k+offs,v) for k,v in x.items()),ladjp)

        lpt  = [k-offp for k in self.Gs.adj.keys() if k <0 ]
        lseg = [k+offs for k in self.Gs.adj.keys() if k >0 ]

        dpt  = dict(zip(lpt,nladjp))
        dseg = dict(zip(lseg,nladjs))
        dseg.update(dpt)
        self.Gs.adj =  dseg
        self.Gs.edge = dseg
        #pdb.set_trace()
        #dict(zip(self.Gs.keys(),))
        #self.Gs.adj = newapoint


    def check(self,level=0):
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

        """
        consistent = True
        nodes = self.Gs.nodes()
        if len(nodes)>0:
            #
            # points
            # segments
            # degree of segments
            useg  = filter(lambda x : x>0,nodes)
            upnt  = filter(lambda x : x<0,nodes)
            degseg  = map(lambda x : nx.degree(self.Gs,x),useg)

            #
            # 1)   all segments have degree 2
            #
            assert(np.all(array(degseg)==2))

            #
            # degree of points
            # maximum degree of points
            #
            degpnt = map(lambda x : nx.degree(self.Gs,x),upnt)  # points absolute degrees
            degmin = min(degpnt)
            degmax = max(degpnt)

            #
            #  No isolated points (degree 0)
            #  No points of degree 1
            #
            if (degmin<=1):
                deg0 = filter(lambda x: nx.degree(self.Gs,x)==0,upnt)
                deg1 = filter(lambda x: nx.degree(self.Gs,x)==1,upnt)
                assert (len(deg0)==0), "It exists degree 0 points :  %r" % deg0
                assert (len(deg1)==0), "It exists degree 1 points : %r" % deg1

            self.deg={}
            for deg in range(degmax+1):
                num = filter(lambda x : degpnt[x]==deg,range(len(degpnt))) # position of degree 1 point
                npt = map(lambda x : upnt[x],num)  # number of degree 1 points
                self.deg[deg] = npt


            #
            # check if there is duplicate points or segments
            #
            # TODO argsort x coordinate
            #


            ke = self.Gs.pos.keys()
            x = np.array(map(lambda x : x[0], self.Gs.pos.values()))
            y = np.array(map(lambda x : x[1], self.Gs.pos.values()))
            p   = np.vstack((x,y))
            d1  = p-np.roll(p,1,axis=1)
            sd1 = np.sum(np.abs(d1),axis=0)
            if not sd1.all()!=0:
               lu = np.where(sd1==0)[0]

               for u in lu:
                   if ke[u]>0:
                       self.del_segment(ke[u])
                   if ke[u]<0:
                       self.del_points(ke[u])

               nodes = self.Gs.nodes()
               useg  = filter(lambda x : x>0,nodes)
               upnt  = filter(lambda x : x<0,nodes)


            for s in useg:
                n1, n2 = np.array(self.Gs.neighbors(s))  # node s neighbors
                p1 = np.array(self.Gs.pos[n1])           # p1 --- p2
                p2 = np.array(self.Gs.pos[n2])           #     s
                for n in upnt:
                    if (n < 0) & (n1 != n) & (n2 != n):
                        p = np.array(self.Gs.pos[n])
                        if geu.isBetween(p1, p2, p):
                            print p1
                            print p
                            print p2
                            logging.critical("segment %d contains point %d",s,n)
                            consistent =False
                if level>0:
                    cycle = self.Gs.node[s]['ncycles']
                    if len(cycle)==0:
                        logging.critical("segment %d has no cycle",s)
                    if len(cycle)==3:
                        logging.critical("segment %d has cycle %s",s,str(cycle))
        #
        # check if Gs points are unique 
        #
        P = np.array(self.Gs.pos.values())
        similar = geu.check_point_unicity(P)
        if len(similar) !=0:
            logging.critical("points at index(es) %s in self.Gs.pos are similar",str(similar))
            consistent =False

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

    def g2npy(self):
        """ conversion from graphs to numpy arrays

        Notes
        -----

        This function updates the following arrays:

        self.pt (2xNp)
        self.tahe (2xNs)
        self.tgs
        self.dca  : dictionnary of cycle with an airwall
        self.lsss : list of subsegments

        assert self.pt[self.iupnt[-1]] == self.pt[:,self.iupnt[-1]]

        """


        nodes = self.Gs.nodes()
        # nodes include points and segments

        #segment index
        useg  = filter(lambda x : x>0,nodes)

        #points index
        upnt  = filter(lambda x : x<0,nodes)

        self.upnt = np.array((upnt))
        utmp = np.array(zip(-self.upnt,np.arange(len(self.upnt))))
        mutmp = max(utmp[:,0])
        self.iupnt = -np.ones((mutmp+1),dtype='int')
        self.iupnt[utmp[:,0]]=utmp[:,1]


        # degree of segment nodes
        degseg = map(lambda x : nx.degree(self.Gs,x),useg)

        assert(np.all(array(degseg)==2)) # all segments must have degree 2

        #
        # self.degree : dictionnary (point degree : list of point index)
        #

        degpnt = np.array(map(lambda x : nx.degree(self.Gs,x),upnt))  # points absolute degrees

        # lairwall : list of air wall segments

        lairwall = self.name['AIR']

        #
        #  function to couL.filename)nt airwall connected to a point
        #  probably this is not the faster solution
        #

        def nairwall(nupt):
            lseg = nx.neighbors(self.Gs,nupt)
            n = 0
            for ns in lseg:
                if ns in lairwall:
                    n = n+1
            return n

        nairwall = np.array(map(nairwall,upnt))

        #
        # if a node is connected to N air wall ==> deg = deg - N
        #

        degpnt = degpnt - nairwall

        try:
            degmax = max(degpnt)
        except:
            degmax = 1

        self.degree = {}
        for deg in range(degmax+1):
            num = filter(lambda x : degpnt[x]==deg,range(len(degpnt))) # position of degree 1 point
            npt = np.array(map(lambda x : upnt[x],num))  # number of degree 1 points
            self.degree[deg] = npt


        #
        # convert geometric information in numpy array
        #

        self.pt = np.array(np.zeros([2, len(upnt)]), dtype=float)
        self.tahe = np.array(np.zeros([2, len(useg)]), dtype=int)

        self.Np = len(upnt)
        self.Ns = len(useg)

        self.pt[0,:]= np.array([self.Gs.pos[k][0] for k in upnt])
        self.pt[1,:]= np.array([self.Gs.pos[k][1] for k in upnt])


        self.pg = np.sum(self.pt,axis=1)/np.shape(self.pt)[1]
        self.pg = np.hstack((self.pg,0.))

        ntail = map(lambda x : nx.neighbors(self.Gs,x)[0],useg)
        nhead = map(lambda x : nx.neighbors(self.Gs,x)[1],useg)

        self.tahe[0,:] = np.array(map(lambda x : np.nonzero(np.array(upnt)==x)[0][0],ntail))
        self.tahe[1,:] = np.array(map(lambda x : np.nonzero(np.array(upnt)==x)[0][0],nhead))

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
        #pdb.set_trace()
        if Nsmax >0:
            self.tgs = np.zeros(Nsmax+1,dtype=int)
            rag = np.arange(len(useg))
            self.tgs[self.tsg] = rag

            #
            # calculate normal to segment ta-he
            #
            # This could becomes obsolete once the normal will be calculated at
            # creation of the segment
            #

            X = np.vstack((self.pt[0,self.tahe[0,:]],self.pt[0,self.tahe[1,:]]))
            Y = np.vstack((self.pt[1,self.tahe[0,:]],self.pt[1,self.tahe[1,:]]))

            normx = Y[0,:]-Y[1,:]
            normy = X[1,:]-X[0,:]

            scale = np.sqrt(normx*normx+normy*normy)
            assert (scale.all()>0)
            self.normal = np.vstack((normx,normy,np.zeros(len(scale))))/scale


            #for ks in ds:
            #
            # lsss : list of subsegment
            #
            nsmax  = max(self.Gs.node.keys())
            # Warning
            # -------
            # nsmax can be different from the total number of segments
            # This means that the numerotation of segments do not need to be
            # contiguous.
            # stridess : length is equal to nsmax+1
            # sla is an array of string, index 0 is not used because there is
            # no such segment number.
            #
            self.lsss = []
            self.isss = []

            #self.stridess = np.array(np.zeros(nsmax+1),dtype=int)
            self.stridess = np.empty(nsmax+1,dtype=int)
            # +1 is for discarding index 0 (unused here)
            self.sla  = np.empty((nsmax+1+self.Nss), dtype='S20')
            self.offset = np.empty(nsmax+1+self.Nss,dtype=int)


            # Storing segment normals
            # Handling of subsegments
            #
            # index is for indexing subsegment after the nsmax value
            #
            index = nsmax+1
            for ks in useg:
                k = self.tgs[ks]                        # index numpy
                self.offset[k]=self.Gs.node[ks]['offset']
                self.Gs.node[ks]['norm'] = self.normal[:,k]  # update normal
                nameslab  = self.Gs.node[ks]['name']   # update sla array
                assert nameslab!='', "segment "+str(ks)+ " is not defined"
                self.sla[ks] = nameslab
                # stridess is different from 0 only for subsegments
                self.stridess[ks] = 0                   # initialize stridess[ks]
                #if index==155:
                #    pdb.set_trace()
                if self.Gs.node[ks].has_key('ss_name'): # if segment has sub segment
                    nss = len(self.Gs.node[ks]['ss_name'])  # retrieve number of sseg
                    self.stridess[ks]=index-1           # update stridess[ks] dict
                    for uk,slabname in enumerate(self.Gs.node[ks]['ss_name']):
                        self.lsss.append(ks)
                        self.sla[index] = slabname
                        self.isss.append(index)
                        self.offset[index] = self.Gs.node[ks]['ss_offset'][uk]
                        index = index+1

        # append sub segment normal to normal

        normal_ss = self.normal[:,self.tgs[self.lsss]]
        self.normal = np.hstack((self.normal,normal_ss))


        # if problem here check file format 'z' should be a string
        self.maxheight = np.max([v[1] for v in nx.get_node_attributes(self.Gs,'z').values()])
        # calculate extremum of segments
        self.extrseg()


    def loadosm(self, _fileosm):
        """ load layout from an osm file 

        Parameters
        ----------

        _fileosm : string


        Notes
        -----

        In JOSM nodes are numbered with negative indexes. It is not valid to
        have a positive node number. To remain compliant with the PyLayers
        convention which assumes that <0 nodes are points and >0 nodes are segments,
        in the osm format, segments are numbered negatively with a known offset
        of 1e7=10000000. The convention is set back when loading the osm file.

        """

        self.filename = _fileosm
        fileosm = pyu.getlong(_fileosm,os.path.join('struc','osm'))
        coords,nodes,ways,relations,m = osm.osmparse(fileosm,typ='floorplan')
        _np = 0 # _ to avoid name conflict with numpy alias
        _ns = 0
        ns  = 0
        nss  = 0
        for npt in coords.xy:
            self.Gs.add_node(npt)
            self.Gs.pos[npt] = tuple(coords.xy[npt])
            _np+=1

        for k,nseg in enumerate(ways.way):
            tahe = ways.way[nseg].refs
            if len(tahe)==2:
                nta = tahe[0]
                nhe = tahe[1]
                d  = ways.way[nseg].tags

                # old format conversion
                if d.has_key('zmin'):
                    d['z']=[d['zmin'],d['zmax']]
                    del(d['zmin'])
                    del(d['zmax'])
                if d.has_key('ss_zmin'):
                    d['ss_z']=[[d['ss_zmin'],d['ss_zmax']]]
                    d['ss_name']=[d['ss_name']]
                    del(d['ss_zmin'])
                    del(d['ss_zmax'])
                #/old format conversion
                for key in d:
                    try:
                        d[key]=eval(d[key])
                    except:
                        pass
                # avoid  0 value (not a segment number)
#                ns = k+1
#                # transcode segment index
#                if d.has_key('name'):
#                    name = d['name']
#                else:
#                    name = 'AIR'
#                    d['name'] = 'AIR'
#                self.Gs.add_node(ns)
#                self.Gs.add_edge(nta,ns)
#                self.Gs.add_edge(ns,nhe)
#                self.Gs.node[ns] = d
#                self.Gs.pos[ns] = tuple((np.array(self.Gs.pos[nta])+np.array(self.Gs.pos[nhe]))/2.)
#                if name not in self.display['layers']:
#                    self.display['layers'].append(name)
#                self.labels[ns] = str(ns)
#                if d.has_key('ss_name'):
#                    nss+=len(d['ss_name'])
#                    for n in d['ss_name']:
#                        if n in self.name:
#                            self.name[n].append(ns)
#                        else:
#                            self.name[n]=[ns]
#                if name in self.name:
#                    self.name[name].append(ns)
#                else:
#                    self.name[name] = [ns]
                ns = self.add_segment(nta,nhe,name=d['name'],z=[eval(u) for u in d['z']],offset=0)
               # self.chgmss(s1,ss_name=d['ss_name'],ss_offset=d['ss_offset'],ss_z=d['ss_z'])
                if d.has_key('ss_name'):
                    nss+=len(d['ss_name'])
                    self.chgmss(ns,ss_name=d['ss_name'],ss_z=[[eval(u) for u in v ] for v in d['ss_z']])

#                _ns+=1

        self.Np = _np
        #self.Ns = _ns
        self.Nss = nss
        #del coords
        #del nodes
        #del ways
        #del relations
        # convert graph Gs to numpy arrays for speed up post processing
        self.g2npy()

    def saveosm(self, _fileosm):
        """  save layout in osm file format

        Parameters
        ----------

        _fileosm : string

        """
        fileosm = pyu.getlong(_fileosm,'struc')
        #
        #
        #
        lonmin = -2
        lonmax = -1
        latmin = 47
        latmax = 48
        lon_0 = -1.5
        lat_0 = 47.5
        m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,
            resolution='i',projection='cass',lon_0=lon_0,lat_0=lat_0)
        fd = open(fileosm,"w")
        fd.write("<?xml version='1.0' encoding='UTF-8'?>\n")
        fd.write("<osm version='0.6' upload='false' generator='PyLayers'>\n")

        for n in self.Gs.pos:
            if n <0:
                x,y = self.Gs.pos[n]
                lon,lat = m(x,y,inverse=True)
                fd.write("<node id='"+str(n)+"' action='modify' visible='true' lat='"+str(lat)+"' lon='"+str(lon)+"' />\n")

        for n in self.Gs.pos:
            if n >0:
                neigh = nx.neighbors(self.Gs,n)
                d = self.Gs.node[n]
                noden = -10000000-n
                fd.write("<way id='"+str(noden)+"' action='modify' visible='true'>\n")
                fd.write("<nd ref='"+str(neigh[0])+"' />\n")
                fd.write("<nd ref='"+str(neigh[1])+"' />\n")
                fd.write("<tag k='name' v='"+str(d['name'])+"' />\n")
                fd.write("<tag k='z' v=\""+str(d['z'])+"\" />\n")
                fd.write("<tag k='transition' v='"+str(d['transition'])+"' />\n")
                if d.has_key('ss_name'):
                    ch = str(d['ss_name'])
                    fd.write("<tag k='ss_name' v=\""+ch+"\" />\n")
                    fd.write("<tag k='ss_z' v=\""+str(d['ss_z'])+"\" />\n")
                fd.write("</way>\n")


        fd.write("</osm>\n")

    def saveini(self, _fileini):
        """ save structure in an ini file

        Parameters
        ----------
        _fileini : string
                   short filemame with extension

        """
        config = ConfigParser.ConfigParser()
        config.add_section("info")
        config.add_section("points")
        config.add_section("segments")
        config.add_section("display")
        config.add_section("files")
        config.set("info",'Npoints',self.Np)
        config.set("info",'Nsegments',self.Ns)
        config.set("info",'Nsubsegments',self.Nss)

        for k in self.display:
            config.set("display",k,self.display[k])

        # boundary nodes and air walls are not stored
        for n in self.Gs.pos:
            if n <0:
                if n not in self.lboundary:
                    config.set("points",str(n),self.Gs.pos[n])

        for n in self.Gs.pos:
            if n >0:
                if self.Gs.node[n]['name']!='AIR':
                    d = self.Gs.node[n]
                    # old format conversion
                    if d.has_key('ss_ce1'):
                        del d['ss_ce1']
                    if d.has_key('ss_ce2'):
                        del d['ss_ce2']
                    if d.has_key('zmin'):
                        d['z']=[d['zmin'],d['zmax']]
                        del(d['zmin'])
                        del(d['zmax'])
                    if d.has_key('ss_zmin'):
                        d['ss_z']=[[d['ss_zmin'],d['ss_zmax']]]
                        d['ss_name']=[d['ss_name']]
                        del(d['ss_zmin'])
                        del(d['ss_zmax'])

                    d['connect'] = nx.neighbors(self.Gs,n)
                    try:
                        if d['transition']:
                            pass
                    except:
                        d['transition']=False
                        try:
                            if 'DOOR' in d['ss_name']:
                                d['transition']=True
                        except:
                            pass
                    config.set("segments",str(n),d)
        config.set("files",'materials',self.filematini)
        config.set("files",'slab',self.fileslabini)
        config.set("files",'furniture',self.filefur)
        fileini = pyu.getlong(_fileini,pstruc['DIRINI'])
        fd = open(fileini,"w")
        config.write(fd)
        fd.close()

        # convert graph Gs to numpy arrays for speed up post processing
        # ideally an edited Layout should be locked while not saved.

        self.g2npy()


    def loadini(self, _fileini):
        """ load a structure file from an .ini file

        Parameters
        ----------

        _fileini : string
            file name extension .ini


        """
        self.filename=_fileini
        di = {}
        config = ConfigParser.ConfigParser()
        fileini = pyu.getlong(_fileini,pstruc['DIRINI'])
        config.read(fileini)
        sections = config.sections()
        for section in sections:
            di[section]={}
            options = config.options(section)
            for option in options:
                try:
                    di[section][option] = config.get(section,option)
                except:
                    print section, option

        self.Np = len(di['points'])
        self.Ns = len(di['segments'])
        self.Gs = nx.Graph()
        self.Gs.pos = {}
        self.labels = {}

        # manage ini file with latlon coordinates
        if di['info'].has_key('format'):
            if di['info']['format']=='latlon':
                or_coord_format = 'latlon'
                coords = osm.Coords()
                coords.clean()
                coords.latlon={i:np.array(eval(di['points'][i])) for i in di['points']}
                coords.boundary=np.hstack((np.min(np.array(coords.latlon.values()),axis=0),
                                           np.max(np.array(coords.latlon.values()),axis=0)))
                coords.cartesian(cart=True)
        else :
            or_coord_format = 'cart'
        #
        # update display section
        #
        for k in di['display']:
            try:
                self.display[k]=eval(di['display'][k])
            except:
                self.display[k]=di['display'][k]

        self.ax = self.display['box']

        # update points section
        for nn in di['points']:
            nodeindex = eval(nn)
            if or_coord_format=='latlon':
                x,y = coords.xy[nn]
            else :
                x,y = eval(di['points'][nn])

            #
            # limitation of point precision is important for avoiding
            # topological problems in shapely.
            # Layout precision is hard limited to millimeter precision.
            #

            self.Gs.add_node(nodeindex)  # add point node
            self.Gs.pos[nodeindex] = (round(1000*x)/1000.,round(1000*y)/1000.)
            self.labels[nodeindex] = nn

        # update segments section
        Nss = 0
        for ns in di['segments']:
            self.Gs.add_node(eval(ns)) # add segment node
            d = eval(di['segments'][ns])
            if d.has_key('ss_name'):
                Nss = Nss + len(d['ss_name'])
                ss_offset=[]
                for n in d['ss_name']:
                    if n in self.name:
                        self.name[n].append(eval(ns))
                    else:
                        self.name[n]=[eval(ns)]
                    ss_offset.append(0)
                if not d.has_key('ss_offset'):
                    d['ss_offset']=ss_offset

            name = d['name']
            nta = d['connect'][0]
            nhe = d['connect'][1]
            if not d.has_key('offset'):
                d['offset']=0


            x,y = tuple((np.array(self.Gs.pos[nta])+np.array(self.Gs.pos[nhe]))/2.)
            # round to mm
            self.Gs.pos[eval(ns)] = (round(1000*x)/1000.,round(1000*y)/1000.)
            self.Gs.node[eval(ns)] = d
            self.Gs.add_edge(nta,eval(ns))
            self.Gs.add_edge(eval(ns),nhe)
            if name not in self.display['layers']:
                self.display['layers'].append(name)
            self.labels[eval(ns)] = ns
            if name in self.name:
                self.name[name].append(eval(ns))
            else:
                self.name[name] = [eval(ns)]
        self.Nss = Nss

        # compliant with config file without  material/slab information

        if config.has_section('files'):
            self.filematini=config.get('files','materials')
            self.fileslabini=config.get('files','slab')
            self.filefur=config.get('files','furniture')

        # In this section we handle the ini file format evolution

        if self.display.has_key('fileoverlay'):
            self.display['overlay_file'] = self.display.pop('fileoverlay')
            self.display['overlay_axis'] = self.display['box'] 
            self.saveini(_fileini)

        if self.display.has_key('inverse'):
            self.display['overlay_flip'] = ""
            self.display.pop('inverse')

            self.saveini(_fileini)
        # convert graph Gs to numpy arrays for faster post processing
        #pdb.set_trace()
        self.g2npy()
        # 


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
        filefur = pyu.getlong(_filefur, pstruc['DIRFUR'])
        config = ConfigParser.ConfigParser()
        config.read(filefur)
        furname = config.sections()
        self.lfur = []
        for name in furname:
            F = fur.Furniture()
            F.load(_filefur, name)
            self.lfur.append(F)
        self.filefur=_filefur

    def load(self,_filename):
        """ load a Layout in different formats

        Parameters
        ----------

        _filename : string

        Notes
        -----

        Available formats are :

        +  .ini   : ini file format (natural one) DIRINI
        +  .str2  : native Pyray (C implementation) DIRSTRUC
        +  .str   : binary file with visibility DIRSTRUC
        +  .osm   : opens street map format  DIROSM


        layout files are stored in the directory pstruc['DIRxxx']

        """
        filename,ext=os.path.splitext(_filename)
        if ext=='.osm':
            filename = pyu.getlong(_filename,pstruc['DIROSM'])
            if os.path.exists(filename):
                self.loadosm(_filename)
            else:
                self.filename = _filename
                print "new file",self.filename
        elif ext=='.str':
            filename = pyu.getlong(_filename,pstruc['DIRSTRUC'])
            if os.path.exists(filename):
                self.loadstr(_filename,self.filematini,self.fileslabini)
            else:
                self.filename = _filename
                print "new file",self.filename
        elif ext=='.str2':
            filename = pyu.getlong(_filename,pstruc['DIRSTRUC'])
            if os.path.exists(filename):
                self.loadstr2(_filename,self.filematini,self.fileslabini)
            else:
                self.filename = _filename
                print "new file",self.filename
        elif ext=='.ini':
            filename = pyu.getlong(_filename,pstruc['DIRINI'])
            if os.path.exists(filename):
                self.loadini(_filename)
            else:
                self.filename = _filename
                print "new file",self.filename
        else:
            raise NameError('layout filename extension not recognized')

        self.lbltg=['s']
        #  construct geomfile (.off) for vizualisation with geomview
        self.subseg()
        if os.path.exists(filename):
            try:
                self.geomfile()
            except:
                print "problem to construct geomfile"

    # def loadstr(self, _filename, _filematini='matDB.ini', _fileslabini='slabDB.ini'):
    #     """ loadstr load a .str de PulsRay

    #     Parameters
    #     ----------

    #     _filename : string
    #     _filematini  : string
    #         default 'matDB.ini'
    #     _fileslabini : string
    #         default 'slabDB.ini'

    #     Examples
    #     --------

    #     >>> from pylayers.gis.layout import *
    #     >>> L = Layout()
    #     >>> L.loadstr('defstr.str')

    #     """

    #     self.filename = _filename
    #     self.delete()
    #     mat = sb.MatDB()
    #     mat.load(_filematini)
    #     self.sl = sb.SlabDB()
    #     self.sl.mat = mat
    #     self.sl.load(_fileslabini)
    #     self.labels = {}
    #     self.name = {}
    #     self.Gs.pos = {}
    #     lname = []
    #     filename = pyu.getlong(_filename, pstruc['DIRSTRUC'])
    #     fo = open(filename, "rb")
    #     data = fo.read()
    #     fo.close()

    #     #
    #     # Read : Np Ns Nss
    #     #        Number of Nodes           nn
    #     #        Number of Edges           en
    #     #        Number of Sub Segments    cen
    #     #
    #     data_nn = data[0:4]
    #     Np = stru.unpack('i', data_nn)[0]
    #     data_en = data[4:8]
    #     Ns = stru.unpack('i', data_en)[0]
    #     data_cen = data[8:12]
    #     Nss = stru.unpack('i', data_cen)[0]
    #     self.Np = Np
    #     self.Ns = Ns
    #     self.Nss = Nss

    #     codesl = np.array(np.zeros(Ns), dtype=int)
    #     codes = np.array(np.zeros(Ns), dtype=int)

    #     # tahe : segment tail and head point index
    #     tahe = np.array(np.zeros([2, Ns]), dtype=int)
    #     ini = 12
    #     for i in range(Ns):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         #self.tahe[0,i]= stru.unpack('i',dt)[0]-1
    #         tahe[0, i] = stru.unpack('i', dt)[0] - 1

    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         #self.tahe[1,i]= stru.unpack('i',dt)[0] -1
    #         tahe[1, i] = stru.unpack('i', dt)[0] - 1

    #     # x : tableau des coordonnees x des noeuds
    #     pt = np.array(np.zeros([2, Np], dtype=np.float64))
    #     ini = stop
    #     for i in range(Np):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         pt[0, i] = stru.unpack('d', dt)[0]
    #     # y : tableau des coordinates y des noeuds
    #     ini = stop
    #     for i in range(Np):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         pt[1, i] = stru.unpack('d', dt)[0]
    #     #--------------------------------------------
    #     # Node labelling (structure nodes)
    #     #--------------------------------------------
    #     for k in range(Np):
    #         self.Gs.add_node(-(k + 1))
    #         self.Gs.pos[-(k + 1)] = (pt[0, k], pt[1, k])
    #         self.labels[-(k + 1)] = str(-(k + 1))

    #     #
    #     # y : type de noeud
    #     #
    #     typ = np.array(np.zeros(Np), dtype=int)
    #     codep = np.array(np.zeros(Np), dtype=int)
    #     ini = stop
    #     for i in range(Np):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         typ[i] = stru.unpack('i', dt)[0]
    #         codep[i] = stru.unpack('i', dt)[0]
    #     #
    #     # agi : tableau des angles initiaux des noeuds de type 2
    #     #
    #     ag = np.array(np.zeros([3, Np], dtype=np.float64))
    #     ini = stop
    #     for i in range(Np):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ag[0, i] = stru.unpack('d', dt)[0]
    #     # agf : tableau des angles finaux des noeuds de type 2
    #     ini = stop
    #     for i in range(Np):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ag[1, i] = stru.unpack('d', dt)[0]
    #     # nN : tableau des parametres d'ouverture de diedre des noeuds de type 2
    #     nN = np.array(1.0 * np.zeros(Np))
    #     ini = stop
    #     for i in range(Np):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ag[2, i] = stru.unpack('d', dt)[0]
    #     #eml  =
    #     em = np.array(np.zeros([3, Ns]), dtype=int)
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         em[0, i] = stru.unpack('i', dt)[0]
    #     #emr  =
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         em[1, i] = stru.unpack('i', dt)[0]
    #     #emc  =
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         em[2, i] = stru.unpack('i', dt)[0]
    #         codes[i] = -2
    #         codesl[i] = em[2, i]
    #         name = self.sl.di[codesl[i]]
    #         lname.append(name)
    #     #thickness =
    #     thick = np.array(1.0 * np.zeros(Ns))
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         thick[i] = stru.unpack('d', dt)[0]
    #     #ehmin =
    #     z = np.array(np.zeros([2, Ns], dtype=np.float64))
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         z[0, i] = stru.unpack('d', dt)[0]
    #     #ehmax =
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         z[1, i] = stru.unpack('d', dt)[0]

    #     norm = np.array(np.zeros([2, Ns], dtype=np.float64))
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 16 * i
    #         stop = ini + 16 * (i + 1)
    #         dt1 = data[start:start + 8]
    #         norm[0, i] = stru.unpack('d', dt1)[0]
    #         dt2 = data[start + 8:stop]
    #         norm[1, i] = stru.unpack('d', dt2)[0]
    #     #
    #     # read matrice node-node
    #     #
    #     ini = stop
    #     nd_nd = np.zeros([Np, Np], dtype=int)
    #     for i in range(Np):
    #         for j in range(Np):
    #             k = Np * i + j
    #             start = ini + 4 * k
    #             stop = ini + 4 * (k + 1)
    #             dt = data[start:stop]
    #             nd_nd[i][j] = stru.unpack('i', dt)[0]
    #     #
    #     # read matrice node-edge
    #     #
    #     ini = stop
    #     nd_ed = np.zeros([Ns, Np], dtype=int)
    #     for i in range(Ns):
    #         for j in range(Np):
    #             k = Np * i + j
    #             start = ini + 4 * k
    #             stop = ini + 4 * (k + 1)
    #             dt = data[start:stop]
    #             nd_ed[i][j] = stru.unpack('i', dt)[0]
    #     #
    #     # read mat_i
    #     #
    #     mat_i = np.array(np.zeros(Ns), dtype=int)
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         mat_i[i] = stru.unpack('i', dt)[0]
    #     #
    #     # read mat_d
    #     #
    #     mat_d = np.array(1.0 * np.zeros(Ns))
    #     ini = stop
    #     for i in range(Ns):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         mat_d[i] = stru.unpack('d', dt)[0]
    #     #
    #     # read matrice ed-ed
    #     #
    #     ini = stop
    #     ed_ed = np.zeros([Ns, Ns], dtype=int)
    #     for i in range(Ns):
    #         for j in range(Ns):
    #             k = Ns * i + j
    #             start = ini + 4 * k
    #             stop = ini + 4 * (k + 1)
    #             dt = data[start:stop]
    #             ed_ed[i][j] = stru.unpack('i', dt)[0]

    #     # Sous segments
    #     #
    #     # read ce_core  (A COMPLETER)
    #     #
    #     ce_core = np.array(np.zeros(Nss), dtype=int)
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         ce_core[i] = stru.unpack('i', dt)[0]
    #     #
    #     # read ce_thick
    #     #
    #     ce_thick = np.array(np.zeros(Nss))
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_thick[i] = stru.unpack('d', dt)[0]
    #     #
    #     # read ce_prop_i
    #     #
    #     ce_prop = np.array(np.zeros([2, Nss]), dtype=int)
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         ce_prop[0, i] = stru.unpack('i', dt)[0]
    #     #
    #     # read ce_wall_floor_ceil
    #     #
    #     ce_wall_floor_ceil = np.array(np.zeros(Nss), dtype=int)
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         ce_wall_floor_ceil[i] = stru.unpack('i', dt)[0]
    #     #
    #     # read ce_ed
    #     #
    #     ce_ed = np.array(np.zeros(Nss), dtype=int)
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 4 * i
    #         stop = ini + 4 * (i + 1)
    #         dt = data[start:stop]
    #         ce_ed[i] = stru.unpack('i', dt)[0]
    #     #
    #     # read ce_prop_d
    #     #
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_prop[1, i] = stru.unpack('d', dt)[0]
    #         #   self.ce_prop[i]= stru.unpack('d',dt)[0]
    #     #
    #     # read ce_xmin
    #     #
    #     ce_xmin = np.array(np.zeros(Nss))
    #     ce_xmax = np.array(np.zeros(Nss))
    #     ce_ymin = np.array(np.zeros(Nss))
    #     ce_ymax = np.array(np.zeros(Nss))
    #     ce_zmin = np.array(np.zeros(Nss))
    #     ce_zmax = np.array(np.zeros(Nss))
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_xmin[i] = stru.unpack('d', dt)[0]
    #     #
    #     # read ce_xmax
    #     #
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_xmax[i] = stru.unpack('d', dt)[0]
    #     #
    #     # read ce_ymin
    #     #
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_ymin[i] = stru.unpack('d', dt)[0]
    #     #
    #     # read ce_ymax
    #     #
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_ymax[i] = stru.unpack('d', dt)[0]
    #     #
    #     # read ce_zmin
    #     #
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_zmin[i] = stru.unpack('d', dt)[0]
    #     #
    #     # read ce_zmax
    #     #
    #     ini = stop
    #     for i in range(Nss):
    #         start = ini + 8 * i
    #         stop = ini + 8 * (i + 1)
    #         dt = data[start:stop]
    #         ce_zmax[i] = stru.unpack('d', dt)[0]

    #     ce = {}
    #     for i in range(Nss):
    #         ce[ce_ed[i] - 1] = (ce_core[i],
    #                             ce_wall_floor_ceil[i],
    #                             ce_prop[0, i],
    #                             ce_zmin[i],
    #                             ce_zmax[i],
    #                             ce_xmin[i],
    #                             ce_xmax[i],
    #                             ce_ymin[i],
    #                             ce_ymax[i])
    #     #self.udbox()
    #     #self.laylist()
    #     #for i in self.layl:
    #     #    self.display['Layer'].append(i)
    #     #    self.display['ActiveLayer'].append(i)

    #     #----------------------------------------
    #     # Node labelling (structure edges)
    #     #----------------------------------------
    #     self.display['layers']=[]
    #     for k in range(Ns):
    #         self.Gs.add_node(k + 1, name=lname[k])
    #         self.Gs.add_node(k + 1, z=(z[0, k],z[1, k]))
    #         self.Gs.add_node(k + 1, norm=np.array([norm[0, k],
    #                                                norm[1, k], 0.]))
    #         nta = tahe[0, k]
    #         nhe = tahe[1, k]
    #         self.Gs.pos[k + 1] = ((pt[0, nta] + pt[0, nhe]) / 2.,
    #                               (pt[1, nta] + pt[1, nhe]) / 2.)
    #         self.Gs.add_edge(-(nta + 1), k + 1)
    #         self.Gs.add_edge(k + 1, -(nhe + 1))
    #         self.labels[k + 1] = str(k + 1)
    #         if lname[k] not in self.display['layers']:
    #             self.display['layers'].append(lname[k])

    #         if lname[k] in self.name:
    #             self.name[lname[k]].append(k + 1)
    #         else:
    #             self.name[lname[k]] = [k + 1]
    #     #
    #     # Update sub-segment
    #     #
    #     for k in ce:
    #         self.Gs.add_node(k + 1, ss_name=[self.sl.di[ce[k][0]]])
    #         self.Gs.add_node(k + 1, ss_ce=[(ce[k][1],ce[k][2])])
    #         self.Gs.add_node(k + 1, ss_z=[(ce[k][3],ce[k][4])])

    #     self.ndnd = nd_nd
    #     self.eded = ed_ed
    #     self.nded = nd_ed
    #     #
    #     # Create connectivity graph Gc
    #     #   update Gc with nd_nd ed_ed
    #     #
    #     self.Gc = nx.Graph()
    #     self.Gc.add_nodes_from(self.Gs.nodes())
    #     pos = self.Gs.pos
    #     #
    #     # !! Incomplet  (To Do nd_ed)
    #     #
    #     Np = np.shape(nd_nd)[0]
    #     for k in range(Np):
    #         nnp = -(k + 1)
    #         kvu = np.nonzero(nd_nd[k] == 3)
    #         nc = -kvu[0] - 1
    #         for l in nc:
    #             self.Gc.add_edge(nnp, l)

    #     Ns = np.shape(ed_ed)[0]
    #     for k in range(Ns):
    #         ne = k + 1
    #         kvu = np.nonzero(ed_ed[k] != 0)
    #         nc = kvu[0] + 1
    #         for l in nc:
    #             self.Gc.add_edge(ne, l)


    #     self.Gc.pos = pos
    #     #
    #     # The numpy format is conserved for acceleration
    #     #
    #     self.pt = pt
    #     self.tahe = tahe
    #     self.display['activelayer'] = self.sl.keys()[0]
    #     #
    #     # update boundary
    #     #
    #     self.boundary(1, 1)

    def loadstr2(self, _filename, _filematini='matDB.ini', _fileslabini='slabDB.ini'):
        """ load a Graph from a str2 file

            Parameters
            ----------
            _filename : string
                str2 filename
            _filematini
                mat filename
            _fileslabini
                slab filename

            Notes
            -----

                str2 format is as follow

                Np Ns Nss
                xp_1 yp_1 zp_1 codep_1
                ...
                xp_Np yp_Np zp_Np codep_Np
                tail_1 head_1 left_1 core_1 right_1 zmin_1 zmax_1
                ...
                tail_Ns head_Ns left_Ns core_Ns right_Ns zmin_Ns zmax_Ns
                segId_1 SlabCode_1 wall_floor_ceil_1 prop_d_1 zmin_1 zmax_1
                ...
                segId_Nss SlabCode_Nss wall_floor_ceil_Nss prop_d_Nss zmin_Nss zmax_Nss

        """

        self.delete()
        self.filename = _filename
        mat = sb.MatDB()
        mat.load(_filematini)

        self.sl = sb.SlabDB()
        self.sl.mat = mat
        self.sl.load(_fileslabini)

        self.labels = {}
        self.name = {}
        self.Gs.pos = {}

        self.Np = 0
        self.Ns = 0
        self.Nss = 0

        filename = pyu.getlong(_filename, pstruc['DIRSTRUC'])
        try:
            fo = open(filename)
        except:
            print "no file named ",filename
            return

        lines = fo.readlines()
        fo.close()
        l1 = lines[0].split()

        #
        # Parse the .str2 header NP NSEG NCOSEG
        #

        Np = int(l1[0])
        Ns = int(l1[1])
        Nss = int(l1[2])

        self.Np = Np
        self.Ns = Ns
        self.Nss = Nss

        lname = []


        pt = np.array(np.zeros([2, Np], dtype=np.float64))
        codep = np.array(np.zeros(Np, dtype=int))
        ag = np.array(np.zeros([3, Np], dtype=np.float64))
        tahe = np.array(np.zeros([2, Ns], dtype=int))
        codesl = np.array(np.zeros(Ns), dtype=int)
        codes = np.array(np.zeros(Ns), dtype=int)

        em = np.array(np.zeros([3, Ns]), dtype=int)
        thick = np.array(np.zeros(Ns))

        ed_mat_prop_d = np.array(np.zeros(Ns, dtype=float))
        height = np.array(np.zeros([2, Ns], dtype=float))
        z = np.array(np.zeros([2, Ns], dtype=float))

        ce_ed = np.array(np.zeros(Nss), dtype=int)
        ce_core = np.array(np.zeros(Nss), dtype=int)
        ce_wall_floor_ceil = np.array(np.zeros(Nss))
        ce_prop_d = np.array(np.zeros(Nss, dtype=np.float64))
        ce_zmin = np.array(np.zeros(Nss, dtype=np.float64))
        ce_zmax = np.array(np.zeros(Nss, dtype=np.float64))
        #
        # Read points
        #
        for i in range(Np):
            dt = lines[i + 1].split()
            pt[0, i] = float(dt[0])
            pt[1, i] = float(dt[1])
            codep[i] = int(dt[2])
            #ag[0:i]=float(dt[3])
            #ag[1:i]=float(dt[4])
            #ag[2:i]=float(dt[5])
        ind1 = i + 2
        #--------------------------------------------
        # Node labelling (structure nodes)
        #--------------------------------------------
        for k in range(Np):
            self.Gs.add_node(-(k + 1))
            self.Gs.pos[-(k + 1)] = (pt[0, k], pt[1, k])
            self.labels[-(k + 1)] = str(-(k + 1))
        #
        # Read segments
        #
        for i in range(Ns):
            dt = lines[i + ind1].split()
            tahe[0, i] = int(dt[0])
            tahe[1, i] = int(dt[1])
            # em : [ left right core ]
            em[1, i] = int(dt[2])
            em[2, i] = int(dt[3])
            em[0, i] = int(dt[4])
            codes[i] = -2
            codesl[i] = em[2, i]
            lname.append(self.sl.di[em[2, i]])
            ed_mat_prop_d[i] = float(dt[5])
            z[0, i] = float(dt[6])
            z[1, i] = float(dt[7])

        ind2 = i + ind1 + 1
        #
        # Read co-segments   ( Transposer dans loadstr)
        #
        ce = {}
        for i in range(Nss):
            dt = lines[i + ind2].split()
            ce_ed[i] = int(dt[0])
            ce_core[i] = int(dt[1])
            ce_wall_floor_ceil[i] = int(dt[2])
            ce_prop_d[i] = float(dt[3])
            ce_zmin[i] = float(dt[4])
            ce_zmax[i] = float(dt[5])
            ce[int(dt[0]) - 1] = (int(dt[1]),
                                  int(dt[2]),
                                  float(dt[3]),
                                  float(dt[4]),
                                  float(dt[5]))

        #----------------------------------------
        # Node labelling (structure edges)
        #----------------------------------------
        self.display['layers']=[]
        for k in range(Ns):
            #print k, lname[k]
            self.Gs.add_node(k + 1, name=lname[k])
            self.Gs.add_node(k + 1, z=(z[0, k],z[1, k]))
            #self.Gs.add_node(k+1,norm=np.array([norm[0,k],norm[1,k],0.]))
            nta = tahe[0, k] - 1
            nhe = tahe[1, k] - 1
            self.Gs.pos[k + 1] = ((pt[0, nta] + pt[0, nhe]) / 2.,
                                  (pt[1, nta] + pt[1, nhe]) / 2.)
            self.Gs.add_edge(-(nta + 1), k + 1)
            self.Gs.add_edge(k + 1, -(nhe + 1))
            self.labels[k + 1] = str(k + 1)
            # update list of layers
            if lname[k] not in self.display['layers']:
                self.display['layers'].append(lname[k])

            if lname[k] in self.name:
                self.name[lname[k]].append(k + 1)
            else:
                self.name[lname[k]] = [k + 1]
        #
        # Update sub-segment
        #
        for k in ce:
            self.Gs.add_node(k + 1, ss_name=[self.sl.di[ce[k][0]]])
            self.Gs.add_node(k + 1, ss_ce= [(ce[k][1],ce[k][2])])
            self.Gs.add_node(k + 1, ss_z = [(ce[k][3],ce[k][4])])

        #
        # Nodes are numbered from 1 in .str2
        # Nodes are numbered from 0 in Graph
        #
        #self.udbox()
        #self.boundary()
        #self.laylist()
        #for i in self.layl:
        #    self.display['Layer'].append(i)
        #    self.display['ActiveLayer'].append(i)

        self.pt = pt
        self.tahe = tahe
        self.display['activelayer'] = self.sl.keys()[0]
        #self.boundary(1,1)

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
                for j,name in enumerate(lname):
                    if name in dico:
                        dico[name].append((k,j))
                    else:
                        dico[name] = [(k,j)]

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
        #print xA,yA
        #print xB,yB
        #print xC,yC
        #print xD,yD
        #print xP,yP
        A = np.array([[xB - xA, xD - xC], [yB - yA, yD - yC]])
        b = np.array([xP - xA, yP - yA])
        x = sp.linalg.solve(A, b)
        if ((x[0] > 0.) & (x[0] < 1.0)):
            self.add_pons(e1, 1 - x[0])
        #print x

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
        self.add_segment(nop[0], num, name=namens, z = [zminns,zmaxns], offset=0)
        # add new edge num np[1]
        self.add_segment(num, nop[1], name=namens, z = [zminns,zmaxns], offset=0)

    def add_segment(self, n1, n2, name='PARTITION',z=[0.0,3.0],offset=0):
        """  add segment between node n1 and node n2

        Parameters
        ----------

        n1  : integer < 0
        n2  : integer < 0
        name : string
            layer name 'PARTITION'
        z : list of float
            default = [0,3.0]
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

        """

        # if 2 points are selected
        if ((n1 < 0) & (n2 < 0) & (n1 != n2)):
            nn = np.array(self.Gs.node.keys())  ## nn : node list array     (can be empty)
            up = np.nonzero(nn > 0)[0]          ## up : segment index (>O)  (can be empty)
            lp = len(up)                        ## lp : number of segments  (can be zero)
            if lp>0:
                e1 = np.arange(lp) + 1          ## e1 : ordered list of segment number
            else:
                e1 = np.array([1])
            e2 = nn[up]                         ## e2 : current list of segment number
            c = ~np.in1d(e1, e2)                ## c  : e1 not in e2 (free segment number)
            tn = e1[c]                          ## tn[c] free segment number
            #print tn
            try:
                num = tn[0]
            except:
                num = max(self.Gs.node.keys()) + 1
                if num == 0:
                    num = 1
        else:
            print "add_segment : error not a node", n1, n2
            return
        transition = False
        if name == 'AIR':
            transition=True

        p1 = np.array(self.Gs.pos[n1])
        p2 = np.array(self.Gs.pos[n2])
        p2mp1 = p2 - p1
        t = p2mp1 / np.sqrt(np.dot(p2mp1, p2mp1))

        #
        # n = t x z
        #

        norm = np.array([t[1], -t[0], 0])
        self.Gs.add_node(num, name=name, 
                              z=z,
                              norm=norm,
                              transition=transition,
                              offset=offset,
                              connect=[n1,n2])

        self.Gs.pos[num] = tuple((p1 + p2) / 2.)
        self.Gs.add_edge(n1, num)
        self.Gs.add_edge(n2, num)
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

    def wedge2(self,apnt):
        """ calculate wedge angle of a point

        Parameters
        ----------

        lpnt : array int
           list of point number


        """

        if isinstance(apnt,list):
            apnt = np.array(apnt)

        ##0. Find the position of diffraction point
        ptdiff = self.pt[:,self.iupnt[-apnt]]

        ##1. Find the associated segments and positions of a diff points

        aseg = map(lambda x : filter(lambda y : y not in self.name['AIR'],
                                         nx.neighbors(self.Gs,x)),
                                         apnt)
        #manage flat angle : diffraction by flat segment e.g. door limitation)
        [aseg[ix].extend(x) for ix,x in enumerate(aseg) if len(x)==1]
        # get points positions
        pts = np.array(map(lambda x : self.seg2pts([x[0],x[1]]),aseg))


        pt1 = pts[:,0:2,0]#tail seg1
        ph1 = pts[:,2:4,0]#head seg1
        pt2 = pts[:,0:2,1]#tail seg2
        ph2 = pts[:,2:4,1]#head seg2


        ##2. Make the correct association
        #pts is (nb_diffraction_points x 4 x 2)
        #- The dimension 4 represent the 2x2 points: t1,h1 and t2,h2
        # tail and head of segemnt 1 and 2 respectively
        #a segment 
        #- The dimension 2 is x,y
        #
        # The following aims to determine which tails and heads of 
        # segments associated to a give diffraction point 
        # are connected

        # point diff is pt1
        updpt1 = np.where(np.sum(ptdiff.T==pt1,axis=1)==2)[0]
        # point diff is ph1
        updph1 = np.where(np.sum(ptdiff.T==ph1,axis=1)==2)[0]

        # point diff is pt2
        updpt2 = np.where(np.sum(ptdiff.T==pt2,axis=1)==2)[0]
        # point diff is ph2
        updph2 = np.where(np.sum(ptdiff.T==ph2,axis=1)==2)[0]

        pa = np.empty((len(apnt),2))
        pb = np.empty((len(apnt),2))

        ####seg 1 :
        #if pt1 diff point =>  ph1 is the other point
        pa[updpt1]= ph1[updpt1]
        #if ph1 diff point =>  pt1 is the other point
        pa[updph1]= pt1[updph1]
        ####seg 2 :
        #if pt2 diff point =>  ph2 is the other point
        pb[updpt2]= ph2[updpt2]
        #if ph2 diff point =>  pt2 is the other point
        pb[updph2]= pt2[updph2]
        #pt is the diffraction point
        pt = ptdiff.T



        vptpa = pt-pa
        vptpan = vptpa.T / np.sqrt(np.sum((vptpa)*(vptpa),axis=1))
        vptpb = pt-pb
        vptpbn = vptpb.T / np.sqrt(np.sum((vptpb)*(vptpb),axis=1))
        v1=vptpan
        v2=vptpbn

        import ipdb
        ipdb.set_trace()
        ang = geu.vecang(vptpbn,vptpan)
        ang[~uleft] = geu.vecang(vptpan,vptpan)
        

    def get_singlGt_angles(self, cy, unit= 'rad', inside=True):

        """ find angles of a single Gt cycle of the layout. 

        Parameters
        ----------
        cy : int
            Gt cyle number
        unit : str
            'deg' : degree values
            'rad' : radian values
        inside : bollean
            True :  compute the inside angles of the cycle.
                    (a.k.a. in regard of the interior of the polygon) 
            False : compute the outside angles of the cycle.
                    (a.k.a.  in regard of the exterior of the polygon)

        Return
        ------

        (u,a)
        u : int (Np)
            point number
        a : float (Np)
            associated angle to the point


        Notes
        -----

        http://www.mathopenref.com/polygonexteriorangles.html

        """

        if cy > 0:
            cycle = self.Gt.node[cy]['cycle'].cycle
        else: #handle outdoorcycle 0
            try:
                cycle = self.ma.vnodes
            except:
                self.ma = self.mask()
                cycle = self.ma.vnodes
        upt = cycle[cycle<0]

        # rupt=np.roll(upt,1)         # for debug
        # rupt2=np.roll(upt,-1)         # for debug
        pt = self.pt[:,self.iupnt[-upt]]
        if geu.SignedArea(pt)<0:
            upt = upt[::-1]
            pt = pt [:,::-1]


        ptroll = np.roll(pt,1,axis=1)

        v = pt-ptroll
        v = np.hstack((v,v[:,0][:,None]))
        vn = v / np.sqrt(np.sum((v)*(v),axis=0))
        v0 = vn[:,:-1]
        v1 = vn[:,1:]
        cross = np.cross(v0.T,v1.T)
        dot = np.sum(v0*v1,axis=0)
        ang = np.arctan2(cross,dot)
        uneg = ang <0
        ang[uneg] = -ang[uneg]+np.pi
        ang[~uneg] = np.pi-ang[~uneg]

        if not inside : 
            ang = 2*np.pi-ang


        if unit == 'deg':
            return upt,ang*180/np.pi
        elif unit == 'rad':
            return upt,ang
        # atan2(cross(a,b)), dot(a,b))


    def get_Gt_angles(self):
        """ find angles of all Gt cycles of the layout. 

        Parameters
        ----------
        None

        Return
        ------

        (u,a)
        u : int (Np)
            point number
        a : float (Np)
            associated angle to the point


        See Also
        --------

        pylayer.layout.get_singlGt_angles
        pylayer.layout.g2npy
        
        """

        dangles = {}
        for cy in self.Gt.nodes():
            uc,ac = self.get_singlGt_angles(cy,inside=True)
            dangles[cy]=np.array(([uc,ac]))
        return dangles

    def wedge(self,lpnt):
        """ calculate wedge angle of a point

        Parameters
        ----------

        lpnt : list of int
           list of point number


        """

        aseg = map(lambda x : filter(lambda y : y not in
                                     self.name['AIR'],
                                     nx.neighbors(self.Gs,x)),
                                     lpnt)

        pts = np.array(map(lambda x : self.seg2pts([x[0],x[1]]).reshape(4,2),aseg))
        #map(lambda x: pt ,pts)
        N = np.shape(pts)[0]
        sector = []
        for k in range(N):
            pt1 = pts[k,0:2,0]
            ph1 = pts[k,2:4,0]
            pt2 = pts[k,0:2,1]
            ph2 = pts[k,2:4,1]
            if (pt1==pt2).all():
                pa = ph1
                pb = ph2
                pt = pt1
                ang = geu.sector(pa,pb,pt)
            if (pt1==ph2).all():
                pa = ph1
                pb = pt2
                pt = pt1
                ang = geu.sector(pa,pb,pt)
            if (ph1==pt2).all():
                pa = pt1
                pb = ph2
                pt = ph1
                ang = geu.sector(pa,pb,pt)
            if (ph1==ph2).all():
                pa = pt1
                pb = pt2
                pt = ph1
                ang = geu.sector(pa,pb,pt)

            sector.append(ang)


        return(sector)


    def add_furniture(self, name='R1_C', matname='PARTITION', origin=(0.,0.),
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
        self.add_segment(n0, n1, matname, [zmin, zmin+height])
        self.add_segment(n1, n2, matname, [zmin, zmin+height])
        self.add_segment(n2, n3, matname, [zmin, zmin+height])
        self.add_segment(n3, n0, matname, [zmin, zmin+height])

    def add_furniture_file(self, _filefur, typ=''):
        """  add pieces of furniture from .ini files

        Parameters
        ----------

        _filefur : string

        """

        filefur = pyu.getlong(_filefur, pstruc['DIRFUR'])
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
            thickness=config.getfloat(fur, "thickness")
            #~ if matname=='WOOD':
                #~ zmin = height
                #~ height=thickness
            #~ else:
                #~ zmin=0.0
            # .. todo: be more generic relate to floor level
            zmin = 0.0
            if typ=='':
                self.add_furniture(name, matname, origin, zmin, height, width, length, angle)
            else:
                try:
                    self.add_furniture(name, matname, origin, zmin, height, width, length, angle)
                except:
                    raise NameError('No such furniture type - '+typ+'-')

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

        print "lp : ",lp
        # get segments involved in points list
        ls = self.nd2seg(lp)

        print "ls : ",ls
        # 1) delete involved segments
        for k in ls:
            assert(k>0)
            self.del_segment(k)

        # 2) delete involved points
        for n1 in lp:
            assert(n1<0)
            nbrs = self.Gs.neighbors(n1)
            self.Gs.remove_node(n1)
            del self.Gs.pos[n1]
            self.labels.pop(n1)
            self.Np = self.Np - 1
        # 3) updating structures
        self.g2npy()

    def del_segment(self,le,verbose=True):
        """ delete segment e

        Parameters
        ----------

        le : list of segment number

        See Also
        --------

        pylayers.gis.layout.Layout.del_node

        """
        if (type(le) == np.ndarray):
            le = list(le)

        if (type(le) != list):
            le = [le]

        for e in le:
            assert(e>0)
            self.del_subseg(e,verbose=verbose)
            name = self.Gs.node[e]['name']
            del self.Gs.pos[e] # delete edge position
            self.Gs.remove_node(e)
            self.labels.pop(e)
            self.Ns = self.Ns - 1
            # update slab name <-> edge number dictionnary
            self.name[name].remove(e)
            # delete subseg if required
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
            if (k!=0) & (k!=-1):
                p = self.Gt.node[k]['polyg']
                ps = ps.union(sh.Polygon(p.exterior))

        mask = geu.Polygon(ps)
        mask.setvnodes(self)
        return(mask)


    def translate(self,vec):
        """ translate layout

        Parameters
        ----------

        vec :

        """
        for k in self.Gs.pos:
            pt=self.Gs.pos[k]
            self.Gs.pos[k]=(pt[0]+vec[0],pt[1]+vec[1])

    def rotate(self,angle=90):
        """ rotate the layout

        Parameters
        ----------

        angle : float
            (degrees)

        """

        a = angle*np.pi/180

        for k in self.Gs.pos:
            pt  = self.Gs.pos[k]
            ptr = np.dot(array([[np.cos(a), -np.sin(a)],[np.sin(a),np.cos(a)]]),array(pt))
            self.Gs.pos[k]=(ptr[0],ptr[1])

        self.g2npy()

    def del_cycle(self, lnc):
        """ delete a cycle

        Parameters
        ----------

        lnc :  list of cycle number

        """

        if (type(lnc) == np.ndarray):
            lnc = list(lnc)

        if (type(lnc) == int):
            lnc = [lnc]

        # for all cycles
        for nc in lnc:
            # get nodes of the cycles
            vnodes = np.array(self.Gt.node[nc]['cycle'].cycle)
            # get neighbors cycles
            neigh = self.Gt.neighbors(nc)
            # array of nodes of neighbors
            tvn = np.array([])
            # for all neihbor cycles
            for ncy in neigh:
                #vn = np.array(self.Gt.node[ncy]['vnodes'])
                # get nodes of neighbor cycle
                vn = np.array(self.Gt.node[ncy]['cycle'].cycle)
                # append nodes in tvn
                try:
                    tvn = np.hstack((tvn, vn))
                except:
                    tvn = vn
            # remove multiple values
            utvn = np.unique(tvn)
            # destroy nodes which are not involved with neighbors
            udel = vnodes[~np.in1d(vnodes, utvn)]

            # delete cycle
            self.Gt.remove_node(nc)
            # delete nodes in udel
            self.del_segment(udel)

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
                print "crosses :", k[0], k[1]
            if seg1.contains(seg2):
                print "contains :", k[0], k[1]
            if seg2.contains(seg1):
                print "contains :", k[0], k[1]
            if seg1.overlaps(seg2):
                print "overlaps :", k[0], k[1]
            if seg2.overlaps(seg1):
                print "overlaps :", k[0], k[1]

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

        displaygui = multenterbox('', 'Display Parameters',
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
                                  (self.filename,
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
            self.filename = displaygui[0]
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
        print n1, ' : ', nns1
        print n2, ' : ', nns2
        print '------------'
        print 'Slab     : ', ds1['name']
        print 'zmin (m) : ', ds1['z'][0]
        print 'zmax (m) : ', ds1['z'][1]
        try:
            print '------------'
            a = ds1['ss_name']
            print 'subseg Slabs  : ', ds1['ss_name']
            print 'subseg (zmin,zmax) (m) : ', ds1['ss_z']
        except:
            pass

    def edit_point(self, np):
        """ edit point

        Parameters
        ----------

        np : integer
            point number

        """
        title = "Point (" + str(np)  + ")"
        message = "Enter coordinates "
        pt = self.Gs.pos[np]
        data = multenterbox(message, title, (('x','y')),
                            ((str(pt[0]),str(pt[1]))))
        self.Gs.pos[np]=tuple(eval(data[0]),eval(data[1]))


    def chgmss(self,ns,ss_name=[],ss_z=[],ss_offset=[]):
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

        if ss_z!=[]:
            assert len(ss_name)==len(ss_z),'Error incompatible size in chgmss'
        if ss_offset!=[]:
            assert len(ss_z)==len(ss_offset),'Error incompatible size in chgmss'

        if ns in self.Gs.node.keys():
            if self.Gs.node[ns].has_key('ss_name'):
                if ss_name!=[]:
                    self.Gs.node[ns]['ss_name']=ss_name
                if ss_z!=[]:
                    self.Gs.node[ns]['ss_z']=ss_z
                if ss_offset!=[]:
                    self.Gs.node[ns]['ss_offset']=ss_offset
                else:
                    self.Gs.node[ns]['ss_offset']=[0]*len(ss_name)


                # update Layout information
                self.g2npy()

    def edit_segment(self, e1 , gui=True,outdata={}):
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
            de1k = ['name', 'z','transition','offset']
            de1v = [de1['name'],de1['z'],de1['transition'],de1['offset']]
        else:
            de1k = ['name', 'z', 'ss_name', 'ss_z','transition','ss_offset']
            de1v = [de1['name'], de1['z'], de1['ss_name'], de1['ss_z'],
                    de1['transition'],de1['ss_offset']]
        #de1v    = de1.values()
        if gui:
            outdata={}
            data0 = choicebox('chose slab',title,self.sl.keys())
            try:
                data1 = multenterbox('attribute for ' + data0, title, tuple(de1k[1:]), tuple(de1v[1:]))
                d1 = data1[0].split(' ')
                d1t = tuple((eval(d1[0]),eval(d1[1])))
                data1[0]=d1t
                data = [data0]+data1
                #data = multenterbox(message, title, tuple(de1k), tuple(de1v))
                i = 0
                self.name[de1['name']].remove(e1)
                for k in de1k:
                    try:
                        self.Gs.node[e1][k] = eval(data[i])
                        outdata[k]=eval(data[i])
                    except:
                        self.Gs.node[e1][k] = data[i]
                        outdata[k]=data[i]
                        if k == 'name':
                            try:
                                self.name[data[i]].append(e1)
                            except:
                                self.name[data[i]] = [e1]
                    i = i + 1
            except:
                #if cancel
                pass
        else:
            if outdata=={}:
                pass
                # data = {}
                # val = '1'
                # while(val!='0'):
                #     clear
                #     print '0 : exit'
                #     for e,(k,v) in enumerate(zip(de1k,de1v)):
                #         print str(e+1)+ ' '+k+': '+  str(v)+'\n'
                #     val = input('Your choice :')
                #     if val!='0':
                #         pass
            else:
                for k in de1k:
                    if k in ['z','name','transition','offset']:
                        self.Gs.node[e1][k] = outdata[k]
        return outdata

    def edit_seg(self, e1 ,data={}):
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

        if data=={}:
            pass
        else:

            ename = self.Gs.node[e1]['name']
            # manage self.name
            self.name[ename].pop(self.name[ename].index(e1))
            # manage self.display['name']
            if len(self.name[ename]) == 0:
                try:
                    self.display['layers'].pop(self.display['layers'].index(ename))
                except:
                    pass

            for k in data:
                self.Gs.node[e1][k] = data[k]


        self.name[data['name']].append(e1)

        if data['name'] not in self.display['layers']:
            self.display['layers'].append(data['name'])

        return data


    def have_subseg(self, e1):
        """
        have_subseg
        """
        dk = self.Gs.node[e1]
        if 'ss_name' in dk:
            return True
        else:
            return False

    def del_subseg(self, e1,verbose=False):
        """ delete sub subsegent

        Parameters
        ----------
        e1 : integer
             segment number
        """
        assert (e1>0)
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
            print "no subseg to delete"

    def add_subseg(self,s1,name='DOOR',zmin=0,zmax=2.24,offset=0,transition=True):
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
        data = multenterbox(message, title, ('name', 'zmin', 'zmax','offset'),
                                            (name, zmin, zmax,offset))
        try:
            self.Gs.node[s1]['ss_name'] = [data[0]]
            self.Nss += 1
            self.chgmss(s1,ss_name=[data[0]],ss_offset=[eval(data[3])],ss_z=[(eval(data[1]),eval(data[2]))])
            #self.Gs.node[s1]['ss_z'] =
            # ce is for Pulsray compatibility
            #self.Gs.node[s1]['ss_ce'] = [ (0,0) ]
            self.Gs.node[s1]['transition'] = transition
            return True
        except:
            return False


    def update_sseg(self,s1,data={}):
        """ update subsegment(s) on a segment

        Parameters
        ----------

        s1 : integer
            edge number > 0
        data = dict
            dictionnary of data

        """

        assert len(data['ss_name'])==len(data['ss_z']),'Error incompatible size in chgmss'
        assert len(data['ss_z'])==len(data['ss_offset']),'Error incompatible size in chgmss'
        if s1<0:
            return False

        new_nbss = len(data['ss_name'])

        try:
            old_nbss = len(self.Gs.node[s1]['ss_name'])
        except:
            old_nbss = 0
        # update the number of subsegments for self.Nss
        deltaNss= new_nbss - old_nbss
        print deltaNss
        if new_nbss != 0:
            self.Gs.node[s1]['ss_name'] = [data['ss_name']]
            self.Nss += deltaNss
            self.chgmss(s1,ss_name=data['ss_name'],
                ss_offset=data['ss_offset'],
                ss_z=data['ss_z'])
            return True
        else:
            if self.Gs.node[s1].has_key('ss_name'):
                self.Gs.node[s1].pop('ss_name')
                self.Gs.node[s1].pop('ss_z')
                self.Gs.node[s1].pop('ss_offset')
                self.Nss += deltaNss
                self.g2npy()
                return True
            else :
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
        if (zmin>self.Gs.node[e1]['z'][0])&(zmax<self.Gs.node[e1]['z'][1]):
            self.info_edge(s1)
            self.Gs.node[s1]['ss_name'].append('WINDOW')
            self.Gs.node[s1]['ss_z'].append((zmin,zmax))
            self.Gs.node[s1]['ss_ce'].append((0,0))
            self.Gs.node[s1]['transition'] =False
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
        if (zmin>self.Gs.node[s1]['z'][0])&(zmax<self.Gs.node[s1]['z'][1]):
            self.info_segment(s1)
            self.Gs.node[s1]['ss_name'].append('DOOR')
            self.Gs.node[s1]['ss_zmin'].append((zmin,zmax))
            self.Gs.node[s1]['ss_ce'].append((0,0))
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
        if isinstance(ndlist,np.ndarray):
            ndlist = ndlist.tolist()

        seglist = []

        #for n in ndlist:
        #    seglist = seglist + self.Gs.adj[n].keys()
        l = map(lambda x :self.Gs.adj[x].keys(),ndlist)
        seglist = reduce(lambda x ,y : x+y,l)

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
        if isinstance(edlist,np.ndarray):
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

    def savestr2(self, _filename='default.str2', furniture=False):
        """ save Layout in .str2 format

        Parameters
        ----------

        _filename : string
            file is  written in the struc directory of the current Project
            directory which is defined through the environment variable $BASENAME
            furniture :  boolean

        Notes
        -----

        To produce the .str file

                > newstruc -str2 file.str2 -conf ../project.conf

        .. todo:: Create a savestr from the Layout Class requires Gv

        """
        if furniture:
            #
            # Create node list
            #

            #
            # Create edge list
            #

            #
            # Create subseg
            #
            pass

        sl = self.sl
        filename = pyu.getlong(_filename,pstruc['DIRSTRUC'])

        nn = self.Np
        ne = self.Ns
        nss = self.Nss

        cnn = str(nn)
        cne = str(ne)
        cnss = str(nss)

        fo = open(filename, 'w')
        #
        # Write in .str2 file
        #
        chaine = cnn + " " + cne + " " + cnss + "\n"
        fo.write(chaine)

        dnode = {}
        ni = 1

        #
        # Reorder segments and points
        #
        # ..todo:: do a spatial reordering
        #

        nodes = np.array(self.Gs.node.keys())
        useg = np.nonzero(nodes > 0)
        upoint = np.nonzero(nodes < 0)
        npoint = nodes[upoint]
        nseg = nodes[useg]

        for i in npoint:
            #
            # points
            #
            x = str(self.Gs.pos[i][0]).replace(',', '.')
            y = str(self.Gs.pos[i][1]).replace(',', '.')
            deg = self.Gs.degree(i)
            if deg > 2:
                deg = 0
            codep = str(deg)
            chaine = x + " " + y + " " + codep + " 0.0 0.0 0.0\n"
            dnode[i] = ni
            ni = ni + 1
            fo.write(chaine)

        for i in nseg:
            #
            # segments
            #
            ta = dnode[self.Gs.neighbors(i)[0]]
            he = dnode[self.Gs.neighbors(i)[1]]
            cta = str(ta)
            che = str(he)
            name = self.Gs.node[i]['name']
            core = str(sl[name]['index'])
            zmin = str(self.Gs.node[i]['z'][0])
            zmax = str(self.Gs.node[i]['z'][1])
            chaine = cta + " " + che + " 1 " + core + " " + " 1 " + \
                " " + " 0 " + " " + zmin + " " + zmax + "\n"
            fo.write(chaine)

        for k, i in enumerate(nseg):
            #
            # sub-segment
            #

            if 'ss_name' in self.Gs.node[i]:

                name = str(self.Gs.node[i]['ss_name'])
                try:
                    core = str(sl[name]['index'])
                except:
                    core = str(sl[eval(name)[0]]['index'])

                if self.Gs.node[i].has_key('ss_ce1'):
                    ce1 = str(self.Gs.node[i]['ss_ce1'][0][0])
                else:
                    ce1 = str(0)

                if self.Gs.node[i].has_key('ss_ce2'):
                    ce2 = str(self.Gs.node[i]['ss_ce2'][0][1])
                else:
                    ce2 = str(0)

                ss_zmin = str(self.Gs.node[i]['ss_z'][0][0])
                ss_zmax = str(self.Gs.node[i]['ss_z'][0][1])
                chaine = str(k + 1) + " " + core + " " + ce1 + \
                    " " + ce2 + " " + ss_zmin + " " + ss_zmax +  "\n"
                fo.write(chaine)

        fo.close()

    def angleonlink3(self, p1=np.array([0,0,1]), p2=np.array([10, 3,1])):
        """ angleonlink(self,p1,p2) return (seglist,angle) between p1 and p2

        Parameters
        ----------

        p1 : np.array (2 x Np) or (2,)
        p2 : np.array (2 x Np) or (2,)

        Returns
        -------

        seglist : list
                  list of segment number on the link
        angle   : angle (in radians) between segment and LOS axis

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('DLR.ini')
        >>> p1 = np.array([0,0,1])
        >>> p2 = np.array([10,3,2])
        >>> alpha = L.angleonlink3(p1,p2)

        #array([(0, 141, 1.2793395519256592), (0, 62, 0.29145678877830505),
               (0, 65, 0.29145678877830505)],
              dtype=[('i', '<i8'), ('s', '<i8'), ('a', '<f4')])


        """

        sh1 = np.shape(p1)
        sh2 = np.shape(p2)

        assert sh1[0]==3
        assert sh2[0]==3

        if (len(sh1)<2) & (len(sh2)>1):
            p1 = np.outer(p1,np.ones(sh2[1]))

        if (len(sh2)<2) & (len(sh1)>1):
            p2 = np.outer(p2,np.ones(sh1[1]))

        if (len(sh2)<2) & (len(sh1)<2):
            p1 = np.outer(p1,np.ones(1))
            p2 = np.outer(p2,np.ones(1))

        # 3 x N
        u = p1 - p2
        # 1 x N
        nu = np.sqrt(np.sum(u*u,axis=0))
        # 3 x N
        un = u / nu[np.newaxis,:]

        seglist = self.seginframe2(p1[0:2], p2[0:2])

        upos = np.nonzero(seglist>=0)[0]
        uneg = np.nonzero(seglist<0)[0]

        nNLOS = len(uneg)+1
        # retrieve the number of segments per link
        if nNLOS>1:
            llink = np.hstack((uneg[0],np.hstack((uneg[1:],array([len(seglist)])))-uneg-1))
        else:
            llink = np.array([len(seglist)])
        # [(link id,number of seg),...]
        #nl = zip(np.arange(nlink),llink)

        npta = self.tahe[0, seglist[upos]]
        nphe = self.tahe[1, seglist[upos]]

        Pta = self.pt[:, npta]
        Phe = self.pt[:, nphe]

        #
        # This part should possibly be improved
        #

        for i,nl in enumerate(llink):
            try:
                P1 = np.hstack((P1,np.outer(p1[:,i],np.ones(nl))))
                P2 = np.hstack((P2,np.outer(p2[:,i],np.ones(nl))))
                ilink = np.hstack((ilink,array([-1]),i*np.ones(nl,dtype='int')))
            except:
                P1 = np.outer(p1[:,i],np.ones(nl))
                P2 = np.outer(p2[:,i],np.ones(nl))
                ilink = i*np.ones(nl,dtype='int')

        # check for intersection P1P2 PtaPhe
        bo = geu.intersect(P1, P2, Pta, Phe)

        upos_intersect = upos[bo]

        seglist2 = seglist[upos_intersect]
        idxlnk = ilink[upos_intersect]

        #
        # Calculate angle of incidence refered from segment normal
        #

        norm  = self.normal[:,seglist2]
        # vector along the link
        uu = un[:,idxlnk]
        unn = abs(np.sum(uu * norm, axis=0))
        angle = np.arccos(unn)

        # seglist = seglist+1
        seglist = np.array(map(lambda x : self.tsg[x],seglist2))
        data = np.zeros(len(seglist),dtype=[('i','i8'),('s','i8'),('a',np.float32)])

        #
        # update subsegment in seglist
        #
        # self.sla
        # self.lsss
        # self.stridess
        #
        sseglist = map(lambda x: self.stridess[x]+1 if x in self.lsss else x,seglist)

        data['i'] = idxlnk
        data['s'] = sseglist
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

        seglist : list
                  list of segment number on the link
        angle   : angle (in radians) between segment and LOS axis

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

        assert sh1[0]==2
        assert sh2[0]==2

        if (len(sh1)<2) & (len(sh2)>1):
            p1 = np.outer(p1,np.ones(sh2[1]))

        if (len(sh2)<2) & (len(sh1)>1):
            p2 = np.outer(p2,np.ones(sh1[1]))

        if (len(sh2)<2) & (len(sh1)<2):
            p1 = np.outer(p1,np.ones(1))
            p2 = np.outer(p2,np.ones(1))

        # 2 x N
        u = p1 - p2
        # 1 x N
        nu = np.sqrt(np.sum(u*u,axis=0))
        # 2 x N
        un = u / nu[np.newaxis,:]

        seglist = self.seginframe2(p1, p2)

        upos = np.nonzero(seglist>=0)[0]
        uneg = np.nonzero(seglist<0)[0]

        nNLOS = len(uneg)+1
        # retrieve the number of segments per link
        if nNLOS>1:
            llink = np.hstack((uneg[0],np.hstack((uneg[1:],array([len(seglist)])))-uneg-1))
        else:
            llink = np.array([len(seglist)])
        # [(link id,number of seg),...]
        #nl = zip(np.arange(nlink),llink)

        npta = self.tahe[0, seglist[upos]]
        nphe = self.tahe[1, seglist[upos]]

        Pta = self.pt[:, npta]
        Phe = self.pt[:, nphe]

        #
        # This part should possibly be improved
        #

        for i,nl in enumerate(llink):
            try:
                P1 = np.hstack((P1,np.outer(p1[:,i],np.ones(nl))))
                P2 = np.hstack((P2,np.outer(p2[:,i],np.ones(nl))))
                ilink = np.hstack((ilink,array([-1]),i*np.ones(nl,dtype='int')))
            except:
                P1 = np.outer(p1[:,i],np.ones(nl))
                P2 = np.outer(p2[:,i],np.ones(nl))
                ilink = i*np.ones(nl,dtype='int')

        bo = geu.intersect(P1, P2, Pta, Phe)

        upos_intersect = upos[bo]

        seglist2 = seglist[upos_intersect]
        idxlnk = ilink[upos_intersect]

        #
        # Calculate angle of incidence refered from segment normal
        #

        norm  = self.normal[0:2,seglist2]
        # vector along the link
        uu = un[:,idxlnk]
        unn = abs(np.sum(uu * norm, axis=0))
        angle = np.arccos(unn)

        # seglist = seglist+1
        seglist = np.array(map(lambda x : self.tsg[x],seglist2))
        data = np.zeros(len(seglist),dtype=[('i','i8'),('s','i8'),('a',np.float32)])

        #
        # update subsegment in seglist
        #
        # self.sla
        # self.lsss
        # self.stridess
        #
        sseglist = map(lambda x: self.stridess[x]+1 if x in self.lsss else x,seglist)

        data['i'] = idxlnk
        data['s'] = sseglist
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

        Pta  = self.pt[:, npta]
        Phe  = self.pt[:, nphe]

        P1 = np.outer(p1, np.ones(len(seglist)))
        P2 = np.outer(p2, np.ones(len(seglist)))

        bo = geu.intersect(P1, P2, Pta, Phe)

        seglist = seglist[bo]

        #
        # Calculate normal angle angle of incidence
        #
        tail = self.tahe[0, seglist]
        head = self.tahe[1, seglist]

        vn  = np.vstack((self.pt[1, head] - self.pt[1, tail],
                        self.pt[0, head] - self.pt[0, tail]))
        mvn = np.outer(np.ones(2), np.sqrt(np.sum(vn * vn, axis=0)))

        n = vn / mvn
        uu = np.outer(un, np.ones(len(seglist)))
        unn = abs(np.sum(uu * n, axis=0))
        theta = np.arccos(unn)

        #print vn
        #print mvn
        #print 'n :',n
        #print 'un : ',unn
        #print 'theta (deg)',the*180./pi

        # seglist = seglist+1
        seglist = np.array(map(lambda x : self.tsg[x],seglist))

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
        idx  = self.tgs[iseg]
        # tahe : 2 x npt
        tahe = self.tahe[:,idx]
        if  len(iseg)>1:
            ta   = tahe[0,:]
            he   = tahe[1,:]
        else:
            ta = tahe[0]
            he = tahe[1]
        pta  = self.pt[:,ta]
        phe  = self.pt[:,he]
        # v  : 2 x npt
        v    = pta-phe
        # mv : npt
        mv   = np.sqrt(np.sum(v*v,axis=0))
        # vn : 2 x npt
        if len(idx)>1:
            vn  = v/mv[np.newaxis,:]
        else:
            vn  = (v/mv).reshape(2)
        return(vn)



    def segpt2(self, ptlist=np.array([0])):
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
        ut = map(lambda x : np.hstack((np.nonzero(self.tahe[0,:]==x)[0],
                                       np.nonzero(self.tahe[1,:]==x)[0])), ptlist)
        utstack = reduce(lambda x,y : np.hstack((x,y)),ut)
        #uvstack = reduce(lambda x,y : np.hstack((x,y)),uv)
        #for i in ptlist:
        #    ut = np.nonzero(self.tahe[0, :] == i)[0]
        #    uv = np.nonzero(self.tahe[1, :] == i)[0]
        #    seglist = np.hstack((seglist, ut, uv))
        seglist = np.unique(utstack)

        return(seglist)

    def seg2pts(self,aseg):
        """ convert segments array to corresponding termination points array

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

        """

        if not isinstance(aseg,np.ndarray):
            aseg = np.array([aseg])

        assert(len(np.where(aseg<0)[0])==0)
        utahe = self.tgs[aseg]
        tahe =  self.tahe[:,utahe]
        ptail =  self.pt[:,tahe[0,:]]
        phead = self.pt[:,tahe[1,:]]
        pth = np.vstack((ptail,phead))

        pth = pth.reshape(pth.shape[0],pth.shape[-1])
        return pth

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
        pt   = self.pt
        # tahe 2 x Nseg
        th  = zip(self.tahe[0,:],self.tahe[1,:])

        self.max_sx = np.array(map(lambda x : max(pt[0,x[0]],pt[0,x[1]]), th))
        self.min_sx = np.array(map(lambda x : min(pt[0,x[0]],pt[0,x[1]]), th))
        self.max_sy = np.array(map(lambda x : max(pt[1,x[0]],pt[1,x[1]]), th))
        self.min_sy = np.array(map(lambda x : min(pt[1,x[0]],pt[1,x[1]]), th))

    def seginframe2(self, p1, p2):
        """ returns the seg list of a given zone defined by two points

            Parameters
            ----------

            p1 array (2 x N)
            p2 array (2 x N)

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

        assert sh1[0]==2
        assert sh2[0]==2

        if (len(sh1)<2) & (len(sh2)>1):
            p1 = np.outer(p1,np.ones(sh2[1]))

        if (len(sh2)<2) & (len(sh1)>1):
            p2 = np.outer(p2,np.ones(sh1[1]))

        if (len(sh2)<2) & (len(sh1)<2):
            p1 = np.outer(p1,np.ones(1))
            p2 = np.outer(p2,np.ones(1))

        # clipping conditions to keep segment
        #
        # max_sx > min_x
        # min_sx < max_x
        # max_sy > min_y
        # min_sy < max_y


        # N x 1

        max_x = map(lambda x : max(x[1],x[0]),zip(p1[0,:], p2[0,:]))
        min_x = map(lambda x : min(x[1],x[0]),zip(p1[0,:], p2[0,:]))
        max_y = map(lambda x : max(x[1],x[0]),zip(p1[1,:], p2[1,:]))
        min_y = map(lambda x : min(x[1],x[0]),zip(p1[1,:], p2[1,:]))

        seglist = map(lambda x : np.nonzero( (self.max_sx > x[0]) &
                                         (self.min_sx < x[1]) &
                                         (self.max_sy > x[2]) &
                                         (self.min_sy < x[3]) )[0],
                                         zip(min_x,max_x,min_y,max_y))

        # np.array stacking
        # -1 acts as a deliminiter (not a segment number)

        seglist = reduce(lambda x,y : np.hstack((x,array([-1]),y)),seglist)

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
            max_x=max_x+0.5
            min_x=min_x-0.5

        if Dy < 0.5:
            max_y=max_y+0.5
            min_y=min_y-0.5



        if (Dy < Dx):
            up = np.nonzero((self.pt[0, :] < max_x) & (self.pt[ 0, :] > min_x))[0]
        else:
            up = np.nonzero((self.pt[1, :] < max_y) & (self.pt[ 1, :] > min_y))[0]

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

        I = np.array([]).reshape(3,0)

        # polygon cycle 1
        poly1 = self.Gt.node[c1]['polyg']
        p1t = poly1.centroid.xy

        # polygon cycle 2
        poly2 = self.Gt.node[c2]['polyg']
        p2t = poly2.centroid.xy

        # centroid of cycle 1 and 2
        p1 = np.array([p1t[0][0],p1t[1][0]])
        p2 = np.array([p2t[0][0],p2t[1][0]])

        line = sh.LineString((p1,p2))

        # els = self.seginframe(p1,p2)
        # new implementation of seginframe is faster
        els = self.seginframe2(p1,p2)

        elg = self.tsg[els]

        lc = []
        ls = []
        I  = np.array([]).reshape(2,0)

        for seg in elg:
            ta, he = self.Gs.neighbors(seg)
            pa = np.array(self.Gs.pos[ta])
            pb = np.array(self.Gs.pos[he])

            segline = sh.LineString((pa,pb))


            if line.intersects(segline):
                lc.extend(self.Gs.node[seg]['ncycles'])
            #print seg,self.Gs.node[seg]['ncycles']
                ls.append(seg)
                psh = line.intersection(segline)
                I = np.hstack((I, np.array([[psh.x],[psh.y]])))
        v = (I-p1[:,np.newaxis])
        dv = np.sum(v*v,axis=0)
        u = np.argsort(dv)
        lss = np.array(ls)[u]

        lc=[c1]
        for s in lss:
            cy1,cy2 = self.Gs.node[s]['ncycles']
            if cy1 not in lc :
                lc.append(cy1)
            elif cy2 not in lc :
                lc.append(cy2)
            else :
                assert NameError('Bad transisiton in Layout.cycleinline')
        return lc




    def seginline(self, p1, p2):
        """
        Returns the intersection between a given line and all segments
        Parameters
        ----------
            p1 : numpy.ndarray
            p2 : numpy.ndarray
        Returns
        -------
            I : numpy.ndarray
        """
        I = np.array([]).reshape(3,0)
        line = sh.LineString((p1,p2))
        for seg in self.Gs.nodes():
            if seg>0:
                ta, he = self.Gs.neighbors(seg)
                pa = np.array(self.Gs.pos[ta])
                pb = np.array(self.Gs.pos[he])
            else:
                pa = np.array(self.Gs.pos[seg])
                pb = pa

            segline = sh.LineString((pa,pb))
            if line.intersects(segline):
                psh = line.intersection(segline)
                liseg = np.array([[psh.x],[psh.y]])
                I = np.hstack((I, np.vstack(([[seg]],liseg))))
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

    def save(self,_filename=[]):
        """ save layout

        Parameters
        ----------

        _filename : short file name (without path)

        """
        if _filename==[]:
            racine, ext = os.path.splitext(self.filename)
            _filename = racine + '.str2'
            _fileini = racine + '.ini'
            self.savestr2(_filename)
            self.saveini(_fileini)
            print "structure saved in ", _filename
            print "structure saved in ", _fileini
        else:
            racine, ext = os.path.splitext(_filename)
            if ext == '.str2':
                self.savestr2(_filename)
                print "structure saved in ", _filename
            if ext == '.ini':
                self.saveini(_filename)
                print "structure saved in ", _filename

   

    def show_nodes(self, ndlist=[1e8], size=10, color='b', dlabels=False, font_size=15, alpha=1,node_shape='o',fig=[],ax=[]):
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
        if fig ==[]:
            fig = plt.figure()
        if ax ==[]:
            ax = fig.add_subplot(111)

        if type(ndlist) == np.ndarray:
            ndlist = list(ndlist)
        if len(ndlist) == 0:
            # ndlist.append(1e8)
            dlabels = False
        elif ndlist[0] == 1e8:
            ndlist = self.Gs.node.keys()
        #elif ndlist[0]==1e8:
        #    ndlist  = self.Gs.node.keys()

        #print ndlist
        Z=nx.draw_networkx_nodes(self.Gs, self.Gs.pos, node_color=color,
                               node_size=size, nodelist=ndlist, alpha=alpha,
                               node_shape=node_shape,fig=fig,ax=ax)
        try:
            fig=Z.figure
            ax=Z.axes
        except:
            pass
        if dlabels:
            dicopos = {}
            dicolab = {}
            for n in ndlist:
                dicopos[n] = np.array(self.Gs.pos[n])
                dicolab[n] = self.labels[n]
            Z=nx.draw_networkx_labels(self.Gs, dicopos, dicolab,
                                    font_size=font_size, font_color=color,fig=fig,ax=ax)
            try:
                fig=Z.figure
                ax=Z.axes
            except:
                pass


        return fig,ax


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

        #print ndlist
        nx.draw_networkx_nodes(
            self.Gs, self.Gs.pos, node_size=size, nodelist=edlist)
        if dlabels:
            dicopos = {}
            dicolab = {}
            for n in ndlist:
                #dicopos[n]=tuple(np.array(self.Gs.pos[n])+np.array((0.8,0.2)))
                dicopos[n] = np.array(self.Gs.pos[n])
                dicolab[n] = self.labels[n]
            nx.draw_networkx_labels(
                self.Gs, dicopos, dicolab, font_size=font_size)

    def show_segment(self,**kwargs):
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

        defaults = { 'fig':[],
                    'ax':[],
                    'edlist': [],
                    'alpha':1,
                    'width':1,
                    'color':'black',
                    'dnodes':False,
                    'dlabels':False,
                    'font_size':15,
                    'node_shape':'o'
                   }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] ==[]:
            fig = plt.figure()
        else :
            fig = kwargs['fig']
        if kwargs['ax'] ==[]:
            ax = fig.add_subplot(111)
        else :
            ax=kwargs['ax']
        clrlist = []
        cold = pyu.coldict()

        # html color or string
        if kwargs['color'][0]!='#':
            clrlist.append(cold[kwargs['color']])
        else:
            if color=='#FFFFF0':
                color = '#00000F'
            clrlist.append(color)
        ecmap = clr.ListedColormap(clrlist)

        U = self.Gs.edges(kwargs['edlist'])
        ue = (np.ones(2 * len(kwargs['edlist']))).astype('int').tolist()
        if len(U) >0:
            Z=nx.draw_networkx_edges(self.Gs, self.Gs.pos, edgelist=U,
                               edge_color=ue, edge_cmap=ecmap,
                               alpha=kwargs['alpha'], width=kwargs['width'],fig=fig,ax=ax)
            try:
                fig=Z.figure
                ax=Z.axes
            except:
                pass

        if kwargs['dlabels']:
               # print edlist
               # nodelist = self.ed2nd(edlist)
            fig,ax=self.show_nodes(ndlist=kwargs['edlist'], dlabels=kwargs['dlabels'],
                            color='b', font_size=kwargs['font_size'],
                            node_shape=kwargs['node_shape'],fig=fig,ax=ax)
        if kwargs['dnodes']:
            fig,ax=self.show_nodes(ndlist=kwargs['edlist'], color='b',fig=fig,ax=ax)

        return fig,ax


    def show_layer(self, name, edlist=[], alpha=1, width=0,
                   color='black', dnodes=False, dthin=False,
                   dlabels=False, font_size=15,fGHz=[],fig=[],ax=[]):
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

        if fig ==[]:
            fig = plt.figure()
        if ax ==[]:
            ax = fig.add_subplot(111)


        if edlist == []:
            edlist = self.name[name]
        else:
            # intersect layer edge list with local zone edge list (in function argument)
            a1 = np.array(self.name[name])
            a2 = np.array(edlist)
            edlist = list(np.intersect1d(a1, a2))

        if self.display['thin']:
            fig,ax=self.show_segment(edlist=edlist,
                              alpha=1,
                              width=1,
                              color=color,
                              dlabels=dlabels,
                              font_size=font_size,fig=fig,ax=ax)
        else:
            slab = self.sl[name]
            if width==0:
                linewidth = slab['linewidth'] / 3.
            else:
                linewidth = width
            if fGHz==[]:
                color = slab['color']
            else:
                if (name!='METAL') & (name!='METALIC'):
                    color = slab.tocolor(fGHz)
                else:
                    color = 'black'

            fig,ax=self.show_segment(edlist=edlist, alpha=1,
                            width=linewidth, color=color, dnodes=dnodes,
                            dlabels=dlabels, font_size=font_size,fig=fig,ax=ax)
        return fig,ax

    def _showGi(self, **kwargs):
        """  show graph of interactions Gi

        Parameters
        ----------

        en  : int
            edge number

        """

        # interactions corresponding to edge en
        int0,int1 = self.Gi.edges()[kwargs['en']]

        print "int0 : ",int0
        print "int1 : ",int1

        # if interaction is tuple (R or T)
        if ((len(int0)>1) & (len(int1)>1)):
            nstr0 = int0[0]
            nstr1 = int1[0]
            output = self.Gi.edge[int0][int1]['output']
            print " output ", output
            ltup = filter(lambda x : type(x)==tuple,output.keys())
            lref = filter(lambda x : len(x)==2,ltup)
            ltran =filter(lambda x : len(x)==3,ltup)
            lseg = np.unique(np.array(map(lambda x : x[0],output.keys())))
            probR = np.array(map(lambda x : output[x],lref))
            segR = np.array(map(lambda x : x[0],lref))
            probT = np.array(map(lambda x : output[x],ltran))
            segT = np.array(map(lambda x : x[0],lref))
            dprobR = dict(zip(segR,probR))
            dprobT = dict(zip(segT,probT))
            #print " Sum pR : ",sum(dprobR.values())
            #print " Sum pT : ",sum(dprobT.values())
            #print "lseg", lseg
            # termination points from seg0 and seg1
            pseg0 = self.seg2pts(nstr0).reshape(2,2).T
            pseg1 = self.seg2pts(nstr1).reshape(2,2).T
            #
            # create the cone seg0 seg1
            #
            cn = cone.Cone()
            cn.from2segs(pseg0,pseg1)
            # show cone
            # show Gt
            self.display['thin']=True
            self.display['subseg']=False
            fig,ax = self.showGs()
            fig,ax = cn.show(fig = fig,ax = ax)
            for nse in lseg:
                ta, he = self.Gs.neighbors(nse)
                pta = np.array(self.Gs.pos[ta])
                phe = np.array(self.Gs.pos[he])

                try:
                    pR= dprobR[nse]
                except:
                    pR = 0

                try:
                    pT = dprobT[nse]
                except:
                    pT = 0

                alpha = (pR+pT)/2.
                segment = ax.plot([pta[0],phe[0]],
                                  [pta[1],phe[1]],
                                   'g',linewidth=7, visible=True,alpha=alpha)

            return(fig,ax)

    def _showGt(self, ax=[], roomlist=[],mode='indoor'):
        """ show topological graph Gt

        Parameters
        -----------
        ax : matlplotlib axes
        roomlist : list
            list of room numbers

        """
        if not isinstance(ax, plt.Axes):
            fig = plt.gcf()
            ax = fig.gca()

        
        G = self.Gt
        
        #pdb.set_trace()
        for k, nc in enumerate(G.node.keys()):
            poly = G.node[nc]['polyg']
            

            a = poly.signedarea()

            if mode == 'area':
                if a < 0:
                    poly.plot(color='red',alpha=0.5,fig=fig,ax=ax)
                else:
                    poly.plot(color='green', alpha=0.5,fig=fig,ax=ax)

            if mode == 'start':
                if poly.vnodes[0] < 0:
                    poly.plot(color='blue',alpha=0.5,fig=fig,ax=ax)
                else:
                    poly.plot(color='yellow', alpha=0.5,fig=fig,ax=ax)

            if mode == 'indoor':   
                if G.node[nc]['indoor']:
                    poly.plot(color='green',alpha=0.5,fig=fig,ax=ax)
                else:
                    poly.plot(color='blue', alpha=0.5,fig=fig,ax=ax) 

            if mode == 'open':   
                if G.node[nc]['isopen']:
                    poly.plot(color='green',alpha=0.5,fig=fig,ax=ax)
                # else:
                #     poly.plot(color='blue', alpha=0.5,fig=fig,ax=ax)         

        ax.axis('scaled')


    def showGs(self,**kwargs):
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


        defaults = {'ndlist' : [],
                    'edlist': [],
                    'roomlist' : [],
                    'axis' : [],
                    'width': 2,
                    'fGHz' : [],
                    'show':False,
                    'furniture':False,
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
            if len(self.display['overlay_file'].split('http:'))>1:
                img_file = urllib.urlopen(self.display['overlay_file'])
                im = StringIO(img_file.read())
                image = Image.open(im)
                imok =True
            else:
                if self.display['overlay_file']!='':
                    image = Image.open(os.path.join(basename,pstruc['DIRIMAGE'],self.display['overlay_file']))
                    imok =True
            if imok:
                if 'v' in self.display['overlay_flip']:
                    image = image.transpose(Image.FLIP_LEFT_RIGHT)
                if 'h' in self.display['overlay_flip']:
                    image = image.transpose(Image.FLIP_TOP_BOTTOM)
                ax.imshow(image, extent=self.display['overlay_axis'],alpha=self.display['alpha'],origin='lower')

        if kwargs['ndlist'] == []:
            tn = np.array(self.Gs.node.keys())
            u = np.nonzero(tn < 0)[0]
            ndlist = tn[u]

        if kwargs['edlist'] == []:
            tn = self.Gs.node.keys()
            #u  = np.nonzero(tn > 0)[0]
            #edlist = tn[u]
            edlist = filter(lambda x: (x>0) ,tn)
            #& (not self.Gs.node[x].has_key('ss_name')),tn)

        if self.display['nodes']:
            dlabels = self.display['ndlabel']
            fig,ax=self.show_nodes(ndlist, size=30, color='k', dlabels=dlabels,node_shape='s',fig=fig,ax=ax)

        if self.display['subsegnb']:
            if hasattr(self,'lsss'):
                seg = self.lsss
                psseg = np.array([[self.Gs.pos[x][0],self.Gs.pos[x][1]] for x in seg])
                nbsseg = np.array([len(self.Gs.node[x]['ss_name']) for x in seg],dtype='int')

                [ax.text(psseg[x,0]+0.2,psseg[x,1]+0.2,str(nbsseg[x]),
                    fontdict={'size':8},ha='center') for x in range(len(seg))]

        if self.display['transition']:
            try:
                segwtrans = [y for y in [x for x in self.Gs.nodes() if x>0 ]if self.Gs.node[y]['transition']]
                posseg = np.array([self.Gs.pos[x] for x in segwtrans])
                normseg = np.array([self.Gs.node[x]['norm'] for x in segwtrans])[:,:2]
                b1 = (posseg-normseg/2)
                b2 = (posseg+normseg/2)
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
            for nameslab in self.display['layers']:
                color = self.sl[nameslab]['color']
                edlist = self.name[nameslab]
                fig,ax=self.show_layer(nameslab, edlist=edlist, alpha=alpha,
                                dthin=dthin, dnodes=dnodes, dlabels=dlabels,
                                color =color,
                                font_size=font_size,
                                width=kwargs['width'],
                                fGHz=kwargs['fGHz'],
                                fig=fig,ax=ax)

        if self.display['subseg']:
            dico = self.subseg()
            for k in dico.keys():
                if kwargs['fGHz']==[]:
                    color = self.sl[k]['color']
                else:
                    if (k!='METAL') & (k!='METALIC'):
                        color = self.sl[k].tocolor(fGHz)
                        #color = 'red'
                    else:
                        color='black'
                        #print k,color
                edlist2 = []
                for ts in dico[k]:
                    edlist2.append(ts[0])
                    #edlist2.append(ts)
                edlist3 = list(set(edlist2).intersection(set(edlist)))
                #print k , color , edlist
                fig,ax=self.show_segment(edlist=edlist3, color=color, alpha=1.0,width=2,fig=fig,ax=ax)

        if self.display['scaled']:
            ax.axis('scaled')
        ax.set_title(self.display['title'])
        #fig = plt.gcf()
        #ax  = fig.axes[0]
        if self.display['ticksoff']:
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            for loc, spine in ax.spines.iteritems():
                spine.set_color('none')

        if kwargs['furniture']:
            if 'lfur' in self.__dict__:
                for fur1 in self.lfur:
                    if fur1.Matname == 'METAL':
                        fig,ax = fur1.show(fig, ax)
            else:
                print "Warning : no furniture file loaded"


        for nr in kwargs['roomlist']:
            ncy = self.Gr.node[nr]['cycle']
            fig,ax=self.Gt.node[ncy]['polyg'].plot(fig=fig,ax=ax)
        if kwargs['axis']==[]:
            ax.axis('scaled')
        else:
            ax.axis(kwargs['axis'])

        if kwargs['show']:
            plt.show()


        return fig,ax

    def build(self, graph='tcvirw',verbose=False):
        """ build graphs

        Parameters
        ----------

        graph : string composed of
            't' : Gt
            'v' : Gv
            'i' : Gi
            'r' : Gr
            'w" : Gw

        Notes
        -----

        This function can build all the graph associated with the Layout. 
        Warning : by default the layout is saved after each build

        """
        # list of built graphs

        if 't' in graph:
            if verbose:
                print "Gt"
            self.buildGt()
            self.lbltg.extend('t')

        if 'r' in graph:
            if verbose:
                print "Gr"
            self.buildGr()
            self.lbltg.extend('r')

        if 'w' in graph and len(self.Gr.nodes())>1:
            self.buildGw()
            self.lbltg.extend('w')
      
        if 'v' in graph:
            if verbose:
                print "Gv"
            self.buildGv()
            self.lbltg.extend('v')

        if 'i' in graph:
            if verbose:
                print "Gi"
            self.buildGi()
            self.outputGi()
            self.lbltg.extend('i')

        if 'w' in graph and len(self.Gr.nodes())>1:
            if verbose:
                print "Gw"
            self.buildGw()
            self.lbltg.extend('w')
            pass

        # dca : dictionnary of cycles which have an air wall
        

        self.dca={}
        for seg,d in self.Gs.node.items():
            if seg >0 :
                if d['name'] == 'AIR':
                    cy = d['ncycles']
                    try:
                        self.dca[cy[0]].append(cy[1])
                    except:
                        self.dca[cy[0]]=[cy[1]]
                    try:
                        self.dca[cy[1]].append(cy[0])
                    except:
                        self.dca[cy[1]]=[cy[0]]

        # f=os.path.splitext(self.filename)
        # if f[1] =='.ini':
        #     self.saveini(self.filename)
        # else :
        #     self.saveini(f[0] +'.ini')

    


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
        path = os.path.join(basename,'struc','gpickle',self.filename)
        if not os.path.isdir(path):
           os.mkdir(path)
        for g in self.lbltg:
            try:
                if g in ['v','i']:
                    gname1 ='G'+g
                    write_gpickle(getattr(self,gname1),os.path.join(basename,'struc','gpickle','G'+g+'_'+self.filename+'.gpickle'))
                else:
                    gname='G'+g
                    write_gpickle(getattr(self,gname),os.path.join(path,'G'+g+'.gpickle'))
            except:
                raise NameError('G'+g+' graph cannot be saved, probably because it has not been built')
        # save dictionnary which maps string interaction to [interactionnode, interaction type]
        if 't' in self.lbltg:
            write_gpickle(getattr(self,'ldiffin'),os.path.join(path,'ldiffin.gpickle'))
            write_gpickle(getattr(self,'ldiffout'),os.path.join(path,'ldiffout.gpickle'))
        write_gpickle(getattr(self,'dca'),os.path.join(path,'dca.gpickle'))


        root,ext = os.path.splitext(self.filename)
        if ext == '.ini':
            self.saveini(self.filename)

    def dumpr(self):
        """ read a dump of given Graph

        Notes
        -----

        graph : string
            't' : Gt
            'r' : Gr
            'v' : Gv
            'i' : Gi


        .gpickle files are store under the struc directory of the project
        specified by the $BASENAME environment variable

        """
        graphs=['t','v','i','r','w']
        path = os.path.join(basename,'struc','gpickle',self.filename)
        for g in graphs:
            try:
                if g in ['v','i']:
                    gname1 ='G'+g
                    setattr(self, gname1, read_gpickle(os.path.join(basename,'struc','gpickle','G'+g+'_'+self.filename+'.gpickle')))
                else:
                    gname='G'+g
                    setattr(self, gname,read_gpickle(os.path.join(path,'G'+g+'.gpickle')))
                self.lbltg.extend(g)
            except:
                pass
                #print 'G',g,' not saved'

        #
        # fixing bug #136
        # update ncycles attributes of Gs from information in Gt
        #
        for k in self.Gs.node:
            if k>0:
                self.Gs.node[k]['ncycles']=[]

        for k in self.Gt.node:
            if k != 0:
                vnodes = self.Gt.node[k]['cycle'].cycle
                if vnodes[0]<0:
                    self.Gt.node[k]['polyg'].vnodes = vnodes
                else:
                    self.Gt.node[k]['polyg'].vnodes = np.roll(vnodes,-1)
                for inode in vnodes:
                    if inode > 0:   # segments
                        if k not in self.Gs.node[inode]['ncycles']:
                            self.Gs.node[inode]['ncycles'].append(k)
                            if len(self.Gs.node[inode]['ncycles'])>2:
                                print inode,self.Gs.node[inode]['ncycles']
                                logging.warning('dumpr : a segment cannot relate more than 2 cycles')
        # if ncycles is a list with only one element the other cycle is the
        # outside region (cycle -1)
        for k in self.Gs.node:
            if k>0:
                if len(self.Gs.node[k]['ncycles'])==1:
                    self.Gs.node[k]['ncycles'].append(-1)
        # load dictionnary which maps string interaction to [interactionnode, interaction type]
        if 't' in graphs :
            setattr(self,'ldiffin', read_gpickle(os.path.join(path,'ldiffin.gpickle')))
            setattr(self,'ldiffout', read_gpickle(os.path.join(path,'ldiffout.gpickle')))
        setattr(self,'dca', read_gpickle(os.path.join(path,'dca.gpickle')))



    def polysh2geu(self,poly):
        """ transform sh.Polygon into geu.Polygon
        """
        try:

            Gsnodes = np.array(self.Gs.nodes())
            # get node coordinates
            nodept = [self.Gs.pos[i] for i in Gsnodes]
            # transform into shapely points
            shpt  = [sh.Point(pt) for pt in nodept]
            # IV 1 get nodes and vnodes
            # Create a ring to avoid taking points inside the polygon.
            # This helps to avoid polygon inside polygons
            # take exterior of polygon. embose it with buffer and find difference with original polygon*.
            # polye = poly.intersection((poly.exterior).buffer(1e-3))

            uvn = np.where([poly.exterior.buffer(1e-3).contains(p) for p in shpt])[0]
            vnodes = Gsnodes[uvn]
            # IV 1.b transform vnodes to an ordered cycle with Cycle class 
            # NOTE ! Using class cycle is MANDATORY
            # because, some extra vnodes can be pickup during the contain 
            # process before
            S = nx.subgraph(self.Gs,vnodes)
            cycle = nx.cycle_basis(S)
            
            if len(cycle) > 1:
                lc = np.array([len(c) for c in cycle])
                dif = abs(lc - len(vnodes))
                ud = np.where(dif == min(dif))[0]
                cycle = cycle[ud]
            else: 
                cycle = cycle[0]
            if cycle[0]>0:
                cycle=np.roll(cycle,-1)
            pos = [self.Gs.pos[c] for c in cycle if c <0]
            # IV 1.c create a new polygon with correct vnodes and correct points
            P = geu.Polygon(p=pos,vnodes=cycle)

        except:
            import ipdb
            ipdb.set_trace()
        return P

    def getangles(self,poly, unit= 'rad', inside=True):
        """ find angles of a polygo_n

        Parameters
        ----------

        poly : geu.Polygon or sh.Polygon
        unit : str
            'deg' : degree values
            'rad' : radian values
        inside : bollean
            True :  compute the inside angles of the cycle.
                    (a.k.a. in regard of the interior of the polygon) 
            False : compute the outside angles of the cycle.
                    (a.k.a.  in regard of the exterior of the polygon)

        Return
        ------

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

        if isinstance(poly,sh.Polygon):
            poly = polysh2geu(poly)
        
        cycle=poly.vnodes



        upt = cycle[cycle<0]

        # rupt=np.roll(upt,1)         # for debug
        # rupt2=np.roll(upt,-1)         # for debug
        pt = self.pt[:,self.iupnt[-upt]]
        if geu.SignedArea(pt)<0:
            upt = upt[::-1]
            pt = pt [:,::-1]


        ptroll = np.roll(pt,1,axis=1)

        v = pt-ptroll
        v = np.hstack((v,v[:,0][:,None]))
        vn = v / np.sqrt(np.sum((v)*(v),axis=0))
        v0 = vn[:,:-1]
        v1 = vn[:,1:]
        cross = np.cross(v0.T,v1.T)
        dot = np.sum(v0*v1,axis=0)
        ang = np.arctan2(cross,dot)
        uneg = ang <0
        ang[uneg] = -ang[uneg]+np.pi
        ang[~uneg] = np.pi-ang[~uneg]

        if not inside : 
            ang = 2*np.pi-ang


        if unit == 'deg':
            return upt,ang*180/np.pi
        elif unit == 'rad':
            return upt,ang
            # atan2(cross(a,b)), dot(a,b))


    # def _delaunay(self,poly,polyholes=[]):
    #     """ make a Delaunay partitioning of a polygon

    #         If polyhole == []

                       
    #             if a cycle is non convex

    #             1- find its polygon
    #             2- partition polygon into convex polygons (Delaunay)
    #             3- try to merge partitioned polygons in order to obtain
    #                the minimal number of convex polygons


    #         If polyholes != []
    
    #             polygon poly contains holes (polyholes)

    #         This methods returns a partitioning of the polygon poly 
    #         into several convex polygons (voronoi). 

    #     Parameters
    #     ----------

    #         poly : geu.Polygon
    #         polyhole : list of geu.Polygon


    #     Return
    #     ------
    #         ncpol : list
    #             list of new created polygons

    #     Notes
    #     -----

    #     The algorithm updates the Gt nodes and edges created into self.buildGt
    #     by adding new nodes and new AIR segments.

    #     Called In 
    #     ---------

    #     pylayers.gis.layout.buildGt
    

    #     See Also
    #     --------

    #     pylayers.gis.layout.buildGt
    #     pylayers.gis.layout.add_segment
    #     pylayers.gis.layout.del_segment
    #     pylayers.util.geomutil.Polygon
    #     sp.spatial.Delaunay

    #     """

    #     cvex,ccve = poly.ptconvex2()
    #     ucs = cvex+ccve

    #     # keep all convex points (in + out) to build a Delaunay triangulation

    #     if polyholes != []:
    #         # sum up polyholes to their gathered polygone
    #         cp = cascaded_union(polyholes)
    #         if isinstance(cp,sh.Polygon):
    #             cp=[cp]
    #         cp=[geu.Polygon(c) for c in cp]
    #         [c.setvnodes(self) for c in cp]
    #         tmp=[]
    #         for c in cp :
    #             cvexh,ccveh = c.ptconvex2()
    #             tmp = tmp + cvexh +ccveh
    #         # tmp=[]
    #         # for p in polyholes:
    #         #     cvexh,ccveh = p.ptconvex2()
    #         #     tmp = tmp + cvexh +ccveh

    #         ucs=ucs+tmp


    #     if len(ucs) !=0 :
    #         try:
    #             pucs = array(map(lambda x: self.Gs.pos[x], ucs))
    #         except:
    #             import ipdb
    #             ipdb.set_trace()
    #         pucs = np.vstack((pucs,pucs[-1]))
    #         ####
    #         #### perform a Delaunay Partioning
    #         ####
    #         if len(ucs) >2:
    #             trid=sp.spatial.Delaunay(pucs)
    #             tri =trid.simplices
    #             polys = []
    #             naw = []
    #             for t in tri:
    #                 ts = geu.Polygon(pucs[t])
    #                 # check if the new polygon is contained into
    #                 #the original polygon (non guarantee by Delaunay)
    #                 C0 = poly.contains(ts)
    #                 if polyholes == []:
    #                     C=[False]
    #                     I=0
    #                 else: 
    #                     C = [isinstance(ii.intersection(ts),sh.Polygon) for ii in polyholes]

    #                 # if poly contains triangle but not the polyholes
    #                 # self.pltpoly([ts],color='b')
    #                 if C0 and (not np.any(C) ):
    #                     # self.pltpoly([ts],color='r')
    #                     # plt.draw()
    #                     # import ipdb
    #                     # ipdb.set_trace()
    #                     cp =ts
    #                     cp.setvnodes(self)
    #                     uaw = np.where(cp.vnodes == 0)[0]
    #                     lvn = len(cp.vnodes)
    #                     for i in uaw:
    #                         #keep track of created airwalls, because some
    #                         #of them will be destroyed in step 3.
    #                         naw.append(self.add_segment(
    #                                    cp.vnodes[np.mod(i-1,lvn)],
    #                                    cp.vnodes[np.mod(i+1,lvn)]
    #                                    ,name='AIR'))
    #                     polys.append(cp)
    #         #
    #         # 3. merge delaunay triangulation in order to obtain
    #         #   the larger convex polygons partioning
    #         #

    #         cpolys = []
    #         nbpolys = len(polys)

    #         while polys !=[]:
    #             p = polys.pop(0)
    #             for ip2,p2 in enumerate(polys):
    #                 conv=False
    #                 inter = p.intersection(p2)
    #                 # if 2 triangles have a common segment
    #                 pold = p
    #                 if isinstance(inter,sh.LineString):
    #                     p = p + p2
    #                     if p.isconvex():
    #                         polys.pop(ip2)
    #                         polys.insert(0, p)
    #                         conv=True
    #                         break
    #                     else:
    #                         # if pold not in cpolys:
    #                         #     cpolys.append(pold)
    #                         p = pold
    #             # if (ip2 >= len(polys)):# and (conv):
    #             # if conv :
    #             #     if p not in cpolys:
    #             #         cpolys.append(p)
    #             if not conv:#else:
    #                 if pold not in cpolys:
    #                     cpolys.append(pold)
    #             if len(polys) == 0:
    #                 cpolys.append(p)

    #         #### 4. ensure the correct vnode numerotation of the polygons
    #         #### and remove unecessary airwalls

    #         # ncpol : new created polygons
    #         ncpol = []
    #         vnodes = []
    #         for p in cpolys:
    #             interpoly = poly.intersection(p)
    #             if isinstance(interpoly,sh.MultiPolygon):
    #                 raise AttributeError('multi polygon encountered')
    #             else :
    #                 try:
    #                     ptmp = geu.Polygon(interpoly)
    #                     # ptmp = self.polysh2geu(interpoly)
    #                 except:
    #                     import ipdb
    #                     ipdb.set_trace()

    #             ptmp.setvnodes(self)
    #             ncpol.append(ptmp)
    #             vnodes.extend(ptmp.vnodes)


    #         # if no polyholes
    #         if polyholes == []:
    #             ### 4bis 
    #             # Check if all the original area is covered 
    #             # sometimes, area surrounded by 2 new airwalls is not found
    #             # the following code re-add it.
    #             cpdiff=poly.difference(cascaded_union(cpolys))
    #             if isinstance(cpdiff,sh.Polygon):
    #                 cpdiff=sh.MultiPolygon([cpdiff])
    #             if isinstance(cpdiff,sh.MultiPolygon):
    #                 for cp in cpdiff:
    #                     ptmp = geu.Polygon(cp)
    #                     ptmp.setvnodes(self)
    #                     ncpol.append(ptmp)
    #                     vnodes.extend(ptmp.vnodes)



    #         daw = filter(lambda x: x not in vnodes,naw)

    #         for d in daw:
    #             self.del_segment(d,verbose=False)
    #     return ncpol

    def pltlines(self,lines,fig=[],ax=[],color='r'):
        """  plot a line with a specified color and transparency
        """
        if fig == []:
            fig=plt.gcf()
        if ax == []:
            ax=plt.gca()

        c = np.array([l.xy for l in lines])
        [ax.plot(x[0,:],x[1,:],color=color) for x in c]
        plt.axis(self.ax)
        plt.draw()




    def pltpoly(self,poly,fig=[],ax=[],color='r',alpha=0.2):
        """  plot a polygon with a specified color and transparency
        """
        if fig == []:
            fig=plt.gcf()
        if ax == []:
            ax=plt.gca()
        try:
            mpl = [PolygonPatch(x,alpha=alpha,color=color) for x in poly]
        except:
            mpl = [PolygonPatch(x,alpha=alpha,color=color) for x in [poly]]
        [ax.add_patch(x) for x in mpl]
        plt.axis(self.ax)
        plt.draw()

    def pltvnodes(self,vn,fig=[],ax=[],):
        if fig == []:
            fig=plt.gcf()
        if ax == []:
            ax=plt.gca()
        X=np.array([self.Gs.pos[x] for x in vn])
        [ax.text(x[0],x[1],vn[xx]) for xx,x in enumerate(X)]


    def updateshseg(self):
        seg_connect = {x:self.Gs.node[x]['connect'] for x in self.Gs.nodes() if x >0}
        dpts = {x[0]:(self.Gs.pos[x[1][0]],self.Gs.pos[x[1][1]]) for x in seg_connect.items() }
        self._shseg = {p[0]:sh.LineString(p[1]) for p in dpts.items()}

    def buildGt(self,check=False):
        """ 
        build graph of convex cycles 



        Parameters
        ----------

        check : booolean


        """


        # I . get polygon from segments
        # seg_connect = {x:self.Gs.node[x]['connect'] for x in self.Gs.nodes() if x >0}
        # pdb.set_trace()
        # dpts = {x[0]:(self.Gs.pos[x[1][0]],self.Gs.pos[x[1][1]]) for x in seg_connect.items() }
        # self._shseg = {p[0]:sh.LineString(p[1]) for p in dpts.items()}
        self.updateshseg()
        X=sho.polygonize(self._shseg.values())
        P=[x for x in X]
        NP=[]

        # II . make polygon convex
        # for each polygon :
        # - is polygon convex: 
        #   -yes : you're done !
        #   - no :
        #       -has polygon an inner hole ? 
        #           - yes : delaunay on polygon, excluing the polygon inside ( polyhole)
        #           - no : delaunay on the polygon
        #
        for p in P:
            pin = [z for z in sho.polygonize(p.interiors)]
            # no holes in polygon
            if pin == []:
                # delaunay only if polygon not convex
                if not geu.isconvex(p):
                    A=self._delaunay(p)
                    NP.extend(A)
                else:

                    A=geu.Polygon(p)
                    A.setvnodes(self)
                    NP.extend([A])
            # hole in polygon
            else:
                # print 'hole'
                A=self._delaunay(p,pin)
                NP.extend(A)


        # III . create Gt nodes
        for ui,p in enumerate(NP):
            cyid = ui+1
            outdoor = False
            # III 1.a get vnode associated to the polygon
            # get vnodes not in the correct order
            # uvn = np.where([r.buffer(1e-3).contains(p) for p in shpt])[0]
            # vnodes = Gsnodes[uvn]

            # III 1.b transform vnodes to an ordered cycle with Cycle class 
            # NOTE ! Using class cycle is MANDATORY
            # because, some extra vnodes can be picked up during the contain 
            # process before


            S = nx.subgraph(self.Gs,p.vnodes)
            S.pos={}
            S.pos.update({i:self.Gs.pos[i] for i in S.nodes()})

            cycle = cycl.Cycle(S)

            # III 1.c create a new polygon with correct vnodes and correct points
            # P = geu.Polygon(p=cycle.p,vnodes=cycle.cycle)
            # import ipdb
            # ipdb.set_trace()
            # IV 1.d add node to Gt + position
            #
            

            seg = p.vnodes[p.vnodes>0]
            lair = [x in self.name['AIR'] for x in seg]
            
            if sum(lair)>0:
                isopen = True
            else:
                isopen = False

            # III 1.e add nodes
            self.Gt.add_node(cyid,cycle=cycle,polyg=p,isopen=isopen,indoor=True)
            self.Gt.pos.update({cyid:np.array(p.centroid.xy)[:,0]})

        # IV  create Gt edges
        for n1 in self.Gt.nodes():
            for n2 in self.Gt.nodes():
                if n1!= n2:
                    if self.Gt.node[n1]['polyg'].touches(self.Gt.node[n2]['polyg']):
                        # find common segments
                        seg = np.array([n for n in self.Gt.node[n1]['cycle'].cycle if (n in self.Gt.node[n2]['cycle'].cycle) and (n>0)])
                        # if cycle are connected by at least a segmnet but not a point
                        if len(seg)>0:
                            self.Gt.add_edge(n1,n2,segment=seg)


        #  V update Gs
        #   V Update graph Gs nodes with their cycles information
        #
        #   initialize a void list 'ncycles' for each node of Gs
        #
        
        self._updGsncy()

        # add cycle 0 to boundaries segments
        for s in self.segboundary:
            self.Gs.node[s]['ncycles'].append(0)

        #   VI - Connect cycle 0 to each cycle connected to the layout
        #   boundary
        #

        bp = sho.polygonize([self._shseg[x] for x in self.segboundary])
        bp = [p for p in bp]
        boundary=geu.Polygon(bp[0])
        boundary.setvnodes(self)
        # all segments of the Layout boundary
        nseg = filter(lambda x : x >0 , boundary.vnodes)
        # air segments of the Layout boundary
        nsegair = filter(lambda x : x in self.name['AIR'],nseg)
        # wall segments of the Layout boundary
        nsegwall = filter(lambda x : x not in self.name['AIR'],nseg)

        #
        # ldiffin  : list of indoor diffraction points
        # ldiffout : list of outdoor diffraction points (belong to layout boundary)
        #

        # self.ldiffin  = filter(lambda x : x not in boundary.vnodes,self.ldiff)
        # self.ldiffout = filter(lambda x : x in boundary.vnodes,self.ldiff)

        #
        # boundary adjascent cycles
        #
        adjcyair = np.array(map(lambda x : filter(lambda y: y!=0,
                                      self.Gs.node[x]['ncycles'])[0],nsegair))
        adjcwall = np.array(map(lambda x : filter(lambda y: y!=0,
                                      self.Gs.node[x]['ncycles'])[0],nsegwall))
        # pdb.set_trace()
        # adjcyair = np.unique(adjcyair)
        # adjcwall = np.unique(adjcwall)

        # connect cycles separated by air wall to cycle 0  
        for cy,seg in zip(adjcyair,nsegair):
            self.Gt.node[cy]['indoor'] = False
            self.Gt.node[cy]['isopen'] = True
            self.Gt.add_edge(0,cy,segment=np.array([seg]))

        # connect cycles separated by wall to cycle 0      
        for cy,seg in zip(adjcwall,nsegwall):
            self.Gt.add_edge(0,cy,segment=np.array([seg]))
        
        # VII Handle indoor/outdoor cycle 
        #
        # Rule : A cycle is outdoor if it is separated from an outdoor cycle by an AIR segment 
        #
        # An outdoor cycle has no associated ceil
        #
        for cy in self.Gt.nodes():
            if cy != 0:
                lncy = nx.neighbors(self.Gt,cy)  
                for ncy in lncy:
                    segnum = self.Gt.edge[cy][ncy]['segment'][0]
                    if self.Gs.node[segnum]['name']=='AIR':  
                        if not self.Gt.node[cy]['indoor']:
                            self.Gt.node[ncy]['indoor']=False



        
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


        # TO DO : remove this one and manage add/deletion of seg lines at
        # the end of  self._delaunay
        self.updateshseg()



        if check : 
            print "check len(ncycles) == 2",
            nodes = [i for i in self.Gs.nodes() if i>0]
            cncy = np.array([len(self.Gs.node[i]['ncycles']) for i in nodes])
            ucncyl = np.where(cncy<2)[0]
            ucncym = np.where(cncy>2)[0]
            assert len(ucncyl)==0,"Some segments are connected to LESS than 2 cycles" + str(np.array(nodes)[ucncyl])
            assert len(ucncym)==0,"Some segments are connected to MORE than 2 cycles" + str(np.array(nodes)[ucncym])
            print "passed"


    def _delaunay(self,poly,polyholes=[]):
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


        Return
        ------
            ncpol : list
                list of new created geu.Polygons

        Notes
        -----

        The algorithm updates the Gt nodes and edges created into self.buildGt
        by adding new nodes and new AIR segments.

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
            if not isinstance(polyholes,list):
                polyholes=[polyholes]
            for ph in polyholes:
            # sum up polyholes to their gathered polygone
                pucsh  = np.array(ph.exterior.xy).T
                pucs = np.vstack((pucs,pucsh))


        if len(pucs) !=0 :
            ####
            #### perform a Delaunay Partioning
            ####

            trid=sp.spatial.Delaunay(pucs)
            tri =trid.simplices
            polys = []
            naw = []
            popo=[]
            for t in tri:
                ts = geu.Polygon(pucs[t])
                # check if the new polygon is contained into
                #the original polygon (non guarantee by Delaunay)
                C0 = poly.contains(ts)
                if polyholes == []:
                    C=[False]
                    I=0
                else: 
                    C = [isinstance(ii.intersection(ts),sh.Polygon) for ii in polyholes]

                popo.append(ts)
                # if poly contains triangle but not the polyholes
                # if polyholes !=[]:
                #     self.pltpoly([ts],color='b')
                #     import ipdb
                #     ipdb.set_trace()
                if C0 and (not np.any(C) ):
                    # if polyholes!=[]:
                    #     self.pltpoly([ts],color='r')
                    #     plt.draw()

                    cp =ts
                    cp.setvnodes(self)
                    uaw = np.where(cp.vnodes == 0)[0]
                    lvn = len(cp.vnodes)
                    for i in uaw:
                        #keep track of created airwalls, because some
                        #of them will be destroyed in step 3.
                        naw.append(self.add_segment(
                                   cp.vnodes[np.mod(i-1,lvn)],
                                   cp.vnodes[np.mod(i+1,lvn)]
                                   ,name='AIR'))
                    polys.append(cp)
        #
            # 3. merge delaunay triangulation in order to obtain
            #   the larger convex polygons partioning
            #
            diff = poly.difference(sh.MultiPolygon(polys))
            if isinstance(diff,sh.Polygon):
                diff=sh.MultiPolygon([diff])
            if isinstance(diff,sh.MultiPolygon):
                for d in diff:
                    extra = geu.Polygon(d)
                    extra.setvnodes(self)
                    polys.append(extra)

            cpolys = []
            nbpolys = len(polys)

            while polys !=[]:
                p = polys.pop(0)
                for ip2,p2 in enumerate(polys):
                    conv=False
                    inter = p.intersection(p2)
                    # if 2 triangles have a common segment
                    pold = p
                    if isinstance(inter,sh.LineString):
                        p = p + p2
                        if p.isconvex():
                            polys.pop(ip2)
                            polys.insert(0, p)
                            conv=True
                            break
                        else:
                            # if pold not in cpolys:
                            #     cpolys.append(pold)
                            p = pold
                # if (ip2 >= len(polys)):# and (conv):
                # if conv :
                #     if p not in cpolys:
                #         cpolys.append(p)
                if not conv:#else:
                    if pold not in cpolys:
                        cpolys.append(pold)
                if len(polys) == 0:
                    cpolys.append(p)

            #### 4. ensure the correct vnode numerotation of the polygons
            #### and remove unecessary airwalls

            # ncpol : new created polygons
            ncpol = []
            vnodes = []
            for p in cpolys:
                interpoly = poly.intersection(p)
                if isinstance(interpoly,sh.MultiPolygon):
                    raise AttributeError('multi polygon encountered')
                else :
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
                ### 4bis 
                # Check if all the original area is covered 
                # sometimes, area surrounded by 2 new airwalls is not found
                # the following code re-add it.
                cpdiff=poly.difference(cascaded_union(cpolys))
                if isinstance(cpdiff,sh.Polygon):
                    cpdiff=sh.MultiPolygon([cpdiff])
                if isinstance(cpdiff,sh.MultiPolygon):
                    for cp in cpdiff:
                        ptmp = geu.Polygon(cp)
                        ptmp.setvnodes(self)
                        ncpol.append(ptmp)
                        vnodes.extend(ptmp.vnodes)



            daw = filter(lambda x: x not in vnodes,naw)

            for d in daw:
                self.del_segment(d,verbose=False)


        return ncpol








    def buildGt_old(self,check=True):

        """ 
        DEPRECATED
        build graph of convex cycles 


        Parameters
        ----------

        check : booolean


        """

        def pltpoly(poly,fig=[],ax=[]):
            if fig == []:
                fig=plt.gcf()
            if ax == []:
                ax=plt.gca()
            mpl = [PolygonPatch(x,alpha=0.2) for x in poly]
            [ax.add_patch(x) for x in mpl]
            plt.axis(self.ax)
            plt.draw()

        def pltGt(Gt,fig=[],ax=[]):
            if fig == []:
                fig=plt.gcf()
            if ax == []:
                ax=plt.gca()
            [Gt.node[x]['poly'].plot(fig=fig,ax=ax,alpha=0.2) for x in Gt.nodes()]
            plt.axis(self.ax)
            plt.draw()


        # I. get cycle basis
        C = nx.algorithms.cycles.cycle_basis(self.Gs)
        if C==[]:
            C = [self.Gs]

        #pdb.set_trace()
        # II. create the hull of the layout by merging all polygons 
        # corresponding to cycles basis
        poly=[]
        for k,lnode in enumerate(C):
            npoints = filter(lambda x : x <0 ,lnode)
            coords  = map(lambda x : self.Gs.pos[x],npoints)
            poly.append(sh.Polygon(coords))
        # union all polygons
        #pdb.set_trace()
        ma = cascaded_union(poly)

        # transform into geomutil polygon
        # if  polygon is a layout
        if not isinstance(ma,sh.MultiPolygon):
            ma = geu.Polygon(ma)
            ma.setvnodes(self)

        else :
            # This is a fix for non enclosed layouts 
            # with multiple non joint polygons (a.k.a. a city)
            # raise AttributeError('this is a city')
            macvx = ma.convex_hull

            streets = geu.Polygon(macvx.difference(ma))
            streets.setvnodes(self)


            # add air walls where to close the street & ma polygon
            uaw = np.where(streets.vnodes == 0)[0]
            lvn = len(streets.vnodes)
            for i in uaw:
                #keep trace of created airwalls, because some
                #of them will be destroyed in step 3.
                self.add_segment(
                           streets.vnodes[np.mod(i-1,lvn)],
                           streets.vnodes[np.mod(i+1,lvn)]
                           ,name='AIR')


            ma= self.polysh2geu(macvx)
            ma.setvnodes(self)
            streets.setvnodes(self)
            ma.setvnodes(self)


        self.ma = ma

        ###### III .FIND POLYGONS
        ###
        # polygons of each cycle are found by finding the interesection between 
        # all segments of the layout and the layout hull.
        # The shapely diff return a multipolygon where all polygons corresponds to 
        # a cycle
        #

        # get connected points from segments
        # connect is equivalent to self.tahe and lpos to self.pt
        #
        connect = [self.Gs.node[i]['connect'] for i in self.Gs.nodes() if i>0]
        # get their coordinates
        lpos = np.array([(self.Gs.pos[i[0]],self.Gs.pos[i[1]]) for i in connect])
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
        while  not isinstance(R ,sh.MultiPolygon) and buffersize<1e-3  :
            # create polygon from multiline by given a width to lines
            mlp = ml.buffer(buffersize)
            # the difference between the layout hull and polygons built from lines 
            # returns the ndesired multipolygon
            R = self.ma.difference(mlp)
            # increase size of the buffer
            buffersize = buffersize*10

        if isinstance(R,sh.Polygon):
            R=sh.MultiPolygon([R])

        # assert isinstance(R,sh.MultiPolygon), "Shapely.MultiPolygon decomposition Failed"

        ####################
        #### Manage inner hole in polygons
        #### ------------------------------
        ### This part manages layout not correctly described, where 
        ### polygons remains in the middle of others
        ######

        # if !=0 it means some polygons are inside of others
        # which is not allowed. Some Layout modification will be performed


        Rgeu = []
        contain = {}

        for ur,r in enumerate(R):
            try:

                Rgeu.append(self.polysh2geu(r))
                # self.pltpoly([Rgeu[-1]],color='b')
                # plt.draw()
                # Rgeu.append(geu.Polygon(r))
                # Rgeu[-1].setvnodes(self)
            except:
                print "reject"
            # if area are not the same, it means that there is inner holes in r
            if not np.allclose(Rgeu[ur].area,r.area):
                # detect inclusion
                uc = np.where([Rgeu[ur].contains(p) for p in R])[0]
                contain[ur] = [c for c in uc if c != ur]

        # # split polygons with holes into several polygons without holes

        self.macvx=[]
        for k in contain:

            polyholes = [ Rgeu[i] for i in contain[k] ]
            
            # 1 convexify polyholes
            polyg = self._convex_hull(polyholes)
            polyholes.extend(polyg)
            Rgeu.extend(polyg)
            # 2 delaunay on exterior
            ncpol = self._delaunay(Rgeu[k],polyholes=polyholes)

            Rgeu.pop(k)
            Rgeu.extend(ncpol)

            # add polyhole to convex mask ( macvx) list 
            macvx = cascaded_union(polyholes)
            macvx = geu.Polygon(macvx)
            macvx.setvnodes(self)
            self.macvx.append(macvx)


        ####################
        #### Manage  convex hull of the layout
        #### -------------------


        # polys = self._convex_hull()
        # Rgeu.extend(polys)

        ####################
        #### Manage Non convex polygons
        #### -------------------
        #### 1 . determine which polygons are not convex
        #### 2 . apply a delaunay and tranform a single non convexpolygon
        ###      into several convex ( self.delaunay)
        ###  3. remove old non convex polygon and readd new convex ones.
    
        ncpol={}
        for ur,r in enumerate(Rgeu):
            if not r.isconvex():
                print 'nt cvx'
                ncpol[ur] = self._delaunay(r)

                # plt.ion()
                # self.pltpoly(ncpol[ur],fig=plt.gcf(),ax=plt.gca())
                # plt.show()
                # plt.draw()
                # import ipdb
                # ipdb.set_trace()
        Rgeu = np.delete(Rgeu,ncpol.keys()).tolist()
        [Rgeu.extend(ncpol[k]) for k in ncpol]

       
        self.Gt=nx.Graph()
        self.Gt.pos={}

        #### IV Find Vnodes and Final polygons

        for n in self.Gs.node:
            if n>0:
                self.Gs.node[n]['ncycles']=[]


        ncyid = -1

        sma  = self.ma.vnodes[self.ma.vnodes>0]
        # smac = self.macvx.vnodes[self.macvx.vnodes>0]
        # segma = np.unique(np.concatenate((sma,smac)))
        segma=sma
        # VI  add node 0
        #
        #   This shapely polygon has an interior 
        #    Cycles = 0 exterior cycle (assumed outdoor)

        S = nx.subgraph(self.Gs,self.ma.vnodes)
        S.pos={}
        S.pos.update({i:self.Gs.pos[i] for i in S.nodes()})
        cycle = cycl.Cycle(S,self.ma.vnodes)
        boundary = geu.Polygon(tuple(self.ax),delta=5)
        boundary.vnodes = self.ma.vnodes
        self.Gt.add_node(0,polyg=self.ma,cycle=cycle,indoor=False,isopen=True)
        self.Gt.pos[0]=(self.ax[0],self.ax[2])
       


        # IV 1 get nodes and vnodes
        
        for ui,p in enumerate(Rgeu):
            cyid = ui+1
            outdoor = False
            # IV 1.a get vnode associated to the polygon
            # get vnodes not in the correct order
            # uvn = np.where([r.buffer(1e-3).contains(p) for p in shpt])[0]
            # vnodes = Gsnodes[uvn]

            # IV 1.b transform vnodes to an ordered cycle with Cycle class 
            # NOTE ! Using class cycle is MANDATORY
            # because, some extra vnodes can be picked up during the contain 
            # process before


            S = nx.subgraph(self.Gs,p.vnodes)
            S.pos={}
            S.pos.update({i:self.Gs.pos[i] for i in S.nodes()})

            cycle = cycl.Cycle(S)

            # IV 1.c create a new polygon with correct vnodes and correct points
            # P = geu.Polygon(p=cycle.p,vnodes=cycle.cycle)
            # import ipdb
            # ipdb.set_trace()
            # IV 1.d add node to Gt + position
            #
            

            seg = p.vnodes[p.vnodes>0]
            lair = [x in self.name['AIR'] for x in seg]
            
            if sum(lair)>0:
                isopen = True
            else:
                isopen = False

            # IV 1.e 
            #   + add new node (convex cycle) to Gt 
            #   + add centroid of cycle as position of cycle
            # if ((cyid==40) or (cyid==41)):
            #     pdb.set_trace()
            self.Gt.add_node(cyid,cycle=cycle,polyg=p,isopen=isopen,indoor=True)
            self.Gt.pos.update({cyid:np.array(p.centroid.xy)[:,0]})

        # IV 2. get edges
        for n1 in self.Gt.nodes():
            for n2 in self.Gt.nodes():
                if n1!= n2:
                    if self.Gt.node[n1]['polyg'].buffer(1e-3).touches(self.Gt.node[n2]['polyg']):
                        # find common segments
                        seg = np.array([n for n in self.Gt.node[n1]['cycle'].cycle if (n in self.Gt.node[n2]['cycle'].cycle) and (n>0)])
                        # if cycle are connected by at least a segmnet but not a point
                        if len(seg)>0:
                            self.Gt.add_edge(n1,n2,segment=seg)

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
        
        import ipdb
        ipdb.set_trace()
        self._addoutcy(check)

        #   V 2. add outside cycle (absorbant region index 0 )
        #   if ncycles is a list which has only one element then the adjascent cycle is the
        #   outside region (cycle 0)
         #
        #   VII - Connect cycle 0 to each cycle connected to the layout
        #   boundary
        #

        # all segments of the Layout boundary
        nseg = filter(lambda x : x >0 , boundary.vnodes)
        # air segments of the Layout boundary
        nsegair = filter(lambda x : x in self.name['AIR'],nseg)
        # wall segments of the Layout boundary
        nsegwall = filter(lambda x : x not in self.name['AIR'],nseg)


        #
        # boundary adjascent cycles
        #

        adjcyair = np.array(map(lambda x : filter(lambda y: y!=0,
                                      self.Gs.node[x]['ncycles'])[0],nsegair))
        adjcwall = np.array(map(lambda x : filter(lambda y: y!=0,
                                      self.Gs.node[x]['ncycles'])[0],nsegwall))
        # pdb.set_trace()
        # adjcyair = np.unique(adjcyair)
        # adjcwall = np.unique(adjcwall)

        # connect cycles separated by air wall to cycle 0  
        for cy,seg in zip(adjcyair,nsegair):
            self.Gt.node[cy]['indoor'] = False
            self.Gt.node[cy]['isopen'] = True
            self.Gt.add_edge(0,cy,segment=np.array([seg]))

        # connect cycles separated by wall to cycle 0      
        for cy,seg in zip(adjcwall,nsegwall):
            self.Gt.add_edge(0,cy,segment=np.array([seg]))
        
        # IV Handle indoor/outdoor cycle 
        #
        # Rule : A cycle is outdoor if it is separated from an outdoor cycle by an AIR segment 
        #
        # An outdoor cycle has no associated ceil
        #
        for cy in self.Gt.nodes():
            lncy = nx.neighbors(self.Gt,cy)  
            for ncy in lncy:
                segnum = self.Gt.edge[cy][ncy]['segment'][0]
                if self.Gs.node[segnum]['name']=='AIR':  
                    if not self.Gt.node[cy]['indoor']:
                        self.Gt.node[ncy]['indoor']=False


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

        Update graph Gs segment with their 2 cycles  information
        initialize a void list 'ncycles' for each segment of Gs

        See Also
        --------

        pylayers.gis.layout.buildGt
        pylayers.gis.layout.convexify

        """

        for k in self.Gs.node:
            self.Gs.node[k]['ncycles']=[]
       
        # remove node 0
        Gtnodes = filter(lambda x: x != 0, self.Gt.nodes())
        
        for ncy in Gtnodes:
            #vnodes = np.array(self.Gt.node[k]['vnodes'])

            vnodes = np.array(self.Gt.node[ncy]['cycle'].cycle)
            
            #vnodes = self.Gt.node[k]['polyg'].vnodes
            for n in vnodes:
                if ncy not in self.Gs.node[n]['ncycles']:
                    self.Gs.node[n]['ncycles'].append(ncy)
                    if n>0:
                        if len(self.Gs.node[n]['ncycles'])>2:
                            print n,self.Gs.node[n]['ncycles']
                            logging.warning('A segment cannot relate more than 2 cycles')
             



    def _addoutcy(self,check=False):
        """ 
        Probably use in a future version of buildGt , managing the upcoming inifile
        add outside cycle (absorbant region index 0 )

        Parameters
        ----------

        check : Boolean 

        #   if ncycles is a list which has only one element then the adjascent
        #   cycle is the  outside region (cycle 0)

        """
        seg0=[]
        for macvx in self.macvx:
            seg = [i for i in macvx.vnodes if i >0]
            seg0 = seg0 + seg
        [self.Gs.node[i]['ncycles'].append(0) for i in seg0]
        if check :
            print "check len(ncycles) == 2",
            nodes = [i for i in self.Gs.nodes() if i>0]
            cncy = np.array([len(self.Gs.node[i]['ncycles']) for i in nodes])
            ucncyl = np.where(cncy<2)[0]
            ucncym = np.where(cncy>2)[0]
            assert len(ucncyl)==0,"Some segments are connected to LESS than 2 cycles" + str(np.array(nodes)[ucncyl])
            assert len(ucncym)==0,"Some segments are connected to MORE than 2 cycles" + str(np.array(nodes)[ucncym])
            print "passed"


    def _interlist(self,nodelist=[]):
        """ Construct the list of interactions associated to each cycle


        Parameters
        ----------

        nodelist: list
        list of Gt nodes (cycles) for which interactions have to be found

        Notes
        -----


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

        if nodelist == []:
            nodelist = self.Gt.nodes()
        elif not isinstance(nodelist,list):
            nodelist=[nodelist]

        # for all cycles k (node of Gt)
        for k in nodelist:
            if k !=0:
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
                        if (name!='AIR') & (name!='ABSORBENT'):
                            ListInteractions.append((inode, k))
                        #
                        # Transmission requires 2 cycles separated by a
                        # segment which is different from METAL and ABSORBENT
                        #
                        # (segment number, cycle in , cycle out )
                        if len(cy) == 2:
                            if (name!='METAL') & (name!='ABSORBENT'):
                                ncy = list(cy.difference({k}))[0]
                                ListInteractions.append((inode, k, ncy))
                                ListInteractions.append((inode, ncy, k))
                    else:  # points
                        pass
                # add list of interactions of a cycle
                self.Gt.add_node(k, inter=ListInteractions)

    def _convex_hull(self,mask):
        """
        Add air walls to the layout enveloppe in self.Gs 
        in order the hull of the Layout to be convex.

        Return
        ------

        polys : list of geu.Polygon
            new polygon of the convex hull

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
        if isinstance(P,sh.MultiPolygon):
            for p in P:
                if p.area > 1e-3:
                    polys.append(geu.Polygon(p))
                    polys[-1].setvnodes(self)


        lncy=[]
        for p in polys:
            # p.coorddeter()
            uaw = np.where(p.vnodes == 0)
            for aw in uaw :
                #2 - non existing segments are created as airwalls
                awid = self.add_segment(p.vnodes[aw-1][0], p.vnodes[aw+1][0], name='AIR')
                p.vnodes[aw] = awid


        # U = cascaded_union([mask]+polys)

        # self.macvx = geu.Polygon(U)
        # self.macvx.setvnodes(self)

        return polys


    def _convexify_old(self,check=False):
        """ determine which cycles are not convex and convexify it.

            if a cycle is non convex

            1- find its polygon
            2- partition polygon into convex polygons (Delaunay)
            3- try to merge partitioned polygons in order to obtain
               the minimal number of convex polygons
            4- re-number/re-create Gt by creating  new 'convex' cycles
               from those convex polygons

        Notes
        -----

        The algorithm updates the Gt nodes and edges created into self.buildGt
        by adding new nodes and connnections

        See Also
        --------

        pylayers.gis.layout.buildGt
        pylayers.gis.layout.add_segment
        pylayers.gis.layout.del_segment
        pylayers.util.geomutil.Polygon
        sp.spatial.Delaunay

        """
        #function for debug purpose
        def polyplot(poly,pgs=True):
            if pgs:
                fig,ax=self.showG('s')

            color=['r','b','g']*10
            for ip, p in enumerate(poly):
                fig,ax = p.plot(fig=fig,ax=ax,color=color[ip],alpha =0.5)
        # lacy : list of added cycles
        lacy =[]
        Gt=copy.deepcopy(self.Gt)
        for n in self.Gt.nodes():
            #if indoor cycle
            if n > 0:

                ncy=max(self.Gt.nodes())

                ####
                #### 1 Determine if pt convex in cycle
                ####
                import ipdb
                ipdb.set_trace()
                if not self.Gt.node[n]['polyg'].isconvex():#self.Gt.node[n]['indoor']:
                    no = self.Gt.node[n]['cycle'].cycle
                    cvex,ccve= self.Gt.node[n]['polyg'].ptconvex2()
                    # keep all convex point (in + out) to build teh delaunay triangulation
                    ucs= cvex+ccve

                    if len(ucs) !=0:
                        pucs = array(map(lambda x: self.Gs.pos[x], ucs))
                        pucs = np.vstack((pucs,pucs[-1]))
                        ####
                        #### 2 perform a Delaunay Partioning
                        ####
                        if len(ucs) >2:
                            trid = sp.spatial.Delaunay(pucs)
                            tri = trid.simplices
                            polys = []
                            naw = []
                            for t in tri:
                                ts = geu.Polygon(pucs[t])
                                # check if the new polygon is contained into
                                #the original polygon (non guaratee by Delaunay)
                                C = self.Gt.node[n]['polyg'].contains(ts)
                                if C:
                                    cp =ts
                                    cp.setvnodes(self)
                                    uaw = np.where(cp.vnodes == 0)[0]
                                    lvn = len(cp.vnodes)
                                    for i in uaw:
                                        #keep trace of created airwalls, because some
                                        #of them will be destroyed in step 3.
                                        naw.append(self.add_segment(
                                                   cp.vnodes[np.mod(i-1,lvn)],
                                                   cp.vnodes[np.mod(i+1,lvn)]
                                                   ,name='AIR'))
                                    polys.append(cp)
                        #
                        # 3. merge delaunay triangulation in order to obtain
                        #   the larger convex polygons partioning
                        #
                        cpolys = []
                        nbpolys = len(polys)

                        while polys !=[]:
                            p = polys.pop(0)
                            for ip2,p2 in enumerate(polys):
                                conv=False
                                inter = p.intersection(p2)
                                #if 2 triangles have a common segment
                                pold = p
                                if isinstance(inter,sh.LineString):
                                    p = p + p2
                                    if p.isconvex():
                                        polys.pop(ip2)
                                        polys.insert(0, p)
                                        conv=True
                                        break
                                    else:
                                        # if pold not in cpolys:
                                        #     cpolys.append(pold)
                                        p = pold
                            # if (ip2 >= len(polys)):# and (conv):
                            # if conv :
                            #     if p not in cpolys:
                            #         cpolys.append(p)
                            if not conv:#else:
                                if pold not in cpolys:
                                    cpolys.append(pold)
                            if len(polys) == 0:
                                cpolys.append(p)

                        #### 4. ensure the correct vnode numerotaion of the polygons
                        #### and remove unecessary airwalls

                        # ncpol : new created polygons
                        ncpol = []
                        vnodes=[]

                        for p in cpolys:
                            interpoly = self.Gt.node[n]['polyg'].intersection(p)
                            if isinstance(interpoly,sh.MultiPolygon):
                                raise AttributeError('multi polygon encountered')
                            else :
                                try:
                                    ptmp = geu.Polygon(interpoly)
                                except:
                                    import ipdb
                                    ipdb.set_trace()

                            ptmp.setvnodes(self)
                            ncpol.append(ptmp)
                            vnodes.extend(ptmp.vnodes)

                        ### 4bis 
                        # Check if all the original area is covered 
                        # sometimes, area surrounded by 2 new airwalls is not found
                        # the following code re-add it.
                        cpdiff=self.Gt.node[n]['polyg'].difference(cascaded_union(cpolys))
                        if isinstance(cpdiff,sh.Polygon):
                            cpdiff=sh.MultiPolygon([cpdiff])
                        if isinstance(cpdiff,sh.MultiPolygon):
                            for cp in cpdiff:
                                ptmp = geu.Polygon(cp)
                                ptmp.setvnodes(self)
                                ncpol.append(ptmp)
                                vnodes.extend(ptmp.vnodes)




                        # [i.plot(fig=plt.gcf(),ax=plt.gca(),color='g') for i in ncpol]
                        # plt.draw()

                        #air walls to be deleted (because origin Delaunay triangle
                        # has been merged )
                        daw = filter(lambda x: x not in vnodes,naw)

                        for d in daw:
                            self.del_segment(d,verbose=False)
                        #remove old cycle
                        self.Gt.remove_node(n)
                        nbpolys=len(ncpol)

                        # lcyid: (new) list of cycle id
                        lcyid = [n] + range(ncy+1,ncy+(nbpolys))
                        lacy.extend(lcyid)
                        for ip,p in enumerate(ncpol):

                            # in order to obtain p.xy
                            cyid = lcyid[ip]
                            # replace by new ones
                            lnode = p.vnodes
                            G = nx.subgraph(self.Gs,lnode)
                            G.pos = {}

                            G.pos.update({l: self.Gs.pos[l] for l in lnode})

                            if not len (p.vnodes) == len(G.nodes()):
                                # some point are counted twice during
                                # the "self.Gt.node[n]['polyg'].intersection(p)"
                                #them make a unique with preserved order
                                _, idx = np.unique(p.vnodes, return_index=True)
                                pv = p.vnodes[np.sort(idx)]
                                pp = p.ndarray()[:,np.sort(idx)[::2]/2]
                                p = geu.Polygon(p=pp,vnodes=pv)

                            cy  = cycl.Cycle(G,lnode=p.vnodes)
                            self.Gt.add_node(cyid,cycle=cy)
                            self.Gt.node[cyid]['polyg'] = p#geu.Polygon(p.xy,cy.cycle)
                            self.Gt.node[cyid]['indoor']=True
                            self.Gt.node[cyid]['isopen']=True
                            self.Gt.pos[cyid] = tuple(cy.g)

        for k in combinations(self.Gt.nodes(), 2):
            if not 0 in k:
                vnodes0 = np.array(self.Gt.node[k[0]]['cycle'].cycle)
                vnodes1 = np.array(self.Gt.node[k[1]]['cycle'].cycle)
                #
                # Connect Cycles if they share at least one segments
                #
                intersection_vnodes = np.intersect1d(vnodes0, vnodes1)

                if len(intersection_vnodes) > 1:
                    segment = intersection_vnodes[np.where(intersection_vnodes>0)]
                    self.Gt.add_edge(k[0], k[1],segment= segment)
            else:

                vnodes0 = self.Gt.node[k[0]]['polyg'].vnodes
                vnodes1 = self.Gt.node[k[1]]['polyg'].vnodes
                intersection_vnodes = np.intersect1d(vnodes0, vnodes1)

                if len(intersection_vnodes) > 1:
                    segment = intersection_vnodes[np.where(intersection_vnodes>0)]
                    self.Gt.add_edge(k[0], k[1],segment= segment)


        self.mask()
        # #update self.Gs.node[x]['ncycles']
        self._updGsncy()
        # #add outside cycle to Gs.node[x]['ncycles']
        self._addoutcy(check)
        # #update interaction list into Gt.nodes (cycles)
        self._interlist(nodelist=lacy)


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

        d_id = max(self.Gr.nodes()) # for numerotation of Gw nodes
        d_id_index = d_id + 1

        for e in self.Gr.edges_iter(): # iterator on Gr edges

            self.Gw.add_node(e[0],{'room':e[0],'door':False})
            self.Gw.add_node(e[1],{'room':e[1],'door':False})

            trans1 = self.Gr.node[e[0]]['transition']  # transitions of room e[0]
            trans2 = self.Gr.node[e[1]]['transition']  # transitions of room e[1]
            Id = np.intersect1d(trans1,trans2)[0]  # list of common doors

            unode = self.Gs.neighbors(Id) # get edge number of common doors
            up0 = self.Gs.pos[unode[0]]
            up1 = self.Gs.pos[unode[1]]

            name = self.Gs.node[Id]['name']

            if True:#name != 'AIR':
                pn = self.Gs.node[Id]['norm']
                sl = self.sl[name]
                thick = (sum(sl['lthick'])/2.)+0.2



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
                    self.Gw.add_node(upd0,{'room':e[0],'door':True})
                    # if self.seginline(pdoor0,ep0).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[0],upd0)])
                    d_id_index = d_id_index +1

                    upd1 = d_id_index
                    self.Gw.pos[upd1] = pdoor1
                    self.Gw.add_node(upd1,{'room':e[1],'door':True})
                    # if self.seginline(pdoor1,ep1).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[1],upd1)])
                    d_id_index = d_id_index +1
                else :
                    upd0 = d_id_index
                    self.Gw.pos[upd0] = pdoor0
                    self.Gw.add_node(upd0,{'room':e[1],'door':True})
                    # if self.seginline(pdoor0,ep1).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[1],upd0)])
                    d_id_index = d_id_index +1

                    upd1 = d_id_index
                    self.Gw.pos[upd1] = pdoor1
                    self.Gw.add_node(upd1,{'room':e[0],'door':True})
                    # if self.seginline(pdoor1,ep0).shape[1] <= 1:
                    #     self.Gw.add_edges_from([(e[0],upd1)])
                    d_id_index = d_id_index +1
                self.Gw.add_edges_from([(upd0,upd1)])
            # Airwalls case
            else :

                pdoor = (np.array(up0)+np.array(up1)) / 2
                self.Gw.pos[d_id_index] = pdoor
                self.Gw.add_edges_from([(e[0],d_id_index),(e[1],d_id_index)])
                d_id_index = d_id_index +1

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
        #     pdoor = (np.array(p1) + np.array(p2)) / 2  # middle of the common door

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


        #pdb.set_trace()
        for n in self.Gr.nodes():

            nr = nx.get_node_attributes(self.Gw,'room').items()
            # get the Gw nodes in current room
            f = map(lambda y: y[0],filter(lambda x: x[1] == n,nr))
            for nw in combinations(f,2):
                pf = map(lambda x: self.Gw.pos[x],nw)
                pf =  np.array((pf))
                #if self.seginline(pf[0],pf[1]).shape[1] <= 1:
                d = np.sqrt(np.sum((pf[0]-pf[1])**2))
                self.Gw.add_edges_from([(nw[0],nw[1])],weight=d)


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
    #         pdoor = (np.array(p1) + np.array(p2)) / 2  # middle of the common door

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

    def buildGv(self):

        Gv = nx.Graph()
        for c in self.Gt.nodes():
            if c != 0:

                ext =  self.Gt.node[c]['polyg'].exterior
                vn  = self.Gt.node[c]['polyg'].vnodes
                seg = list(vn[vn>0])
                ptdiff = [ v  for v in vn if v in self.ddiff and c in self.ddiff[v][0]]
                entities = seg+ptdiff
                combit = combinations(entities,2)
                comb =  [(x[0],x[1]) for x in combit]
                lines = [sh.LineString([self.Gs.pos[x[0]],self.Gs.pos[x[1]]]) for x in comb]
                within = np.array([l.within(ext) for l in lines])
                # looking for lines not superimposed with polygon exterior
                ut = np.where(~within)[0]
                comb= np.array(comb)
                [Gv.add_edge(v[0],v[1]) for v in comb[ut]]


                # diff diff processing
                dd = [(x[0],x[1]) for x in comb if x[0]<0 and x[1]<0]
                lined = [sh.LineString([self.Gs.pos[x[0]],self.Gs.pos[x[1]]]) for x in dd]
                uair = [self.Gs.node[x]['name']=='AIR' for x in seg]
                segair = [s for us,s in enumerate(seg) if uair[us]]
                lineair = sh.MultiLineString([self._shseg[x] for x in segair])
                ext2 = ext.difference(lineair)
                withind = np.array([l.within(ext2) for l in lined])
                if len(withind)!=0:
                    utd = np.where(~withind)[0]
                    dd= np.array(dd)
                    [Gv.add_edge(v[0],v[1]) for v in dd[utd]]
                
        self.Gv = Gv 

    # def buildGv(self, show=False):
    #     """ build visibility graph

    #     Parameters
    #     ----------

    #     show : boolean
    #         default False

    #     Examples
    #     --------

    #     >>> from pylayers.gis.layout import *
    #     >>> L = Layout('TA-Office.ini')
    #     >>> L.buildGt()
    #     >>> Ga = L.buildGr()
    #     >>> L.buildGv()

    #     Notes
    #     -----

    #     This method exploits cycles convexity.

    #     """

    #     self.Gv = nx.Graph()
    #     #
    #     # loop over cycles
    #     #
    #     self.dGv = {}  # dict of Gv graph

    #     for icycle in self.Gt.node:
    #         if icycle!=0:
    #             polyg = self.Gt.node[icycle]['polyg']
    #             vnodes = polyg.vnodes

    #             npt  = filter(lambda x : x<0,vnodes)
    #             nseg = filter(lambda x : x>0,vnodes)

    #             airwalls = filter(lambda x : x in self.name['AIR'],nseg)
    #             ndiff = [ x for x in npt if x in self.ddiff.keys()]
    #             #
    #             # Create a graph
    #             #
    #             Gv = nx.Graph()
    #             #
    #             # in convex case :
    #             #
    #             #    i)  every non aligned segments see each other
    #             #
    #             for nk in combinations(nseg, 2):
    #                 nk0 = self.tgs[nk[0]]
    #                 nk1 = self.tgs[nk[1]]
    #                 tahe0 = self.tahe[:,nk0]
    #                 tahe1 = self.tahe[:,nk1]

    #                 pta0 = self.pt[:,tahe0[0]]
    #                 phe0 = self.pt[:,tahe0[1]]
    #                 pta1 = self.pt[:,tahe1[0]]
    #                 phe1 = self.pt[:,tahe1[1]]

    #                 A0 = np.vstack((pta0,phe0,pta1))
    #                 A0 = np.hstack((A0,np.ones((3,1))))

    #                 A1 = np.vstack((pta0,phe0,phe1))
    #                 A1 = np.hstack((A1,np.ones((3,1))))

    #                 d0 = np.linalg.det(A0)
    #                 d1 = np.linalg.det(A1)

    #                 if not ((abs(d0)<1e-3) & (abs(d1)<1e-3)):
    #                     Gv.add_edge(nk[0],nk[1])

    #             #
    #             # Handle diffraction points
    #             #
    #             #    ii) all non adjascent valid diffraction points see each other
    #             #    iii) all valid diffraction points see segments non aligned
    #             #    with adjascent segments

    #             ndiffvalid =  [ x for x in ndiff if icycle in self.ddiff[x][0] ]

    #             # non adjascent segment of vnodes see valid diffraction points
    #             for idiff in ndiffvalid:
    #                 # idiff segment neighbors
    #                 nsneigh = [ x for x in nx.neighbors(self.Gs,idiff) if x in nseg and x not in airwalls]
    #                 # segvalid : not adjascent segment
    #                 # seen from neighbors
    #                 for x in nsneigh:
    #                     seen_from_neighbors.append(nx.neighbors(Gv,x))
    #                 #segvalid = [ x for x in nseg if x in nx.neighbors(Gv,nsneigh)]
    #                 # nbidiff valid
    #                 nsneigh =  [ x for x in nsneigh if x not in airwalls]
    #                 # excluded diffraction points

    #                 for ns in segvalid:
    #                     Gv.add_edge(idiff,ns)

    #                 for npoint in ndiffvalid:
    #                     if npoint !=idiff:
    #                         Gv.add_edge(idiff,npoint)
    #             #
    #             # Graph Gv composition
    #             #

    #             self.Gv  = nx.compose(self.Gv, Gv)
    #           self.dGv[icycle] = Gv


    def buildGi(self):
        """ build graph of interactions

        Notes
        -----

        For each node > of graph Gs creates
        4 different nodes associated to the same segment

        (ns,cy0) R -> cy0
        (ns,cy1) R -> cy1
        (ns,cy0,cy1) T 0->1
        (ns,cy1,cy0) T 1->0

        """
        self.Gi = nx.DiGraph()
        self.Gi.pos = {}
        #
        # 1 ) Create nodes of Gi and their positions
        #
        # (D,)
        # (R,cy0)
        # (T,cy0,cy1)
        #
        for n in self.Gv.node:
            if n < 0: # D
                self.Gi.add_node((n,))
                self.Gi.pos[(n,)] = self.Gs.pos[n]
            if n > 0: # R | T
                cy = self.Gs.node[n]['ncycles']
                name = self.Gs.node[n]['name']
                try:
                    cy0 = cy[0]
                    cy1 = cy[1]
                except:
                    import ipdb
                    ipdb.set_trace()
                nei = self.Gs.neighbors(n)  # get neighbor
                np1 = nei[0]
                np2 = nei[1]

                p1 = np.array(self.Gs.pos[np1])
                p2 = np.array(self.Gs.pos[np2])
                l = p1 - p2
                nl = np.dot(l, l)
                ln = l / nl

                delta = nl / 10
                # On AIR or ABSORBENT there is no reflection
                # except if n is a subsegment
                if ((name!='AIR') & (name!='ABSORBENT')) or (n in self.lsss):
                    self.Gi.add_node((n,cy0))
                    self.Gi.add_node((n,cy1))
                    self.Gi.pos[(n, cy0)] = tuple(self.Gs.pos[n] + ln * delta)
                    self.Gi.pos[(n, cy1)] = tuple(self.Gs.pos[n] - ln * delta)

                # Through METAL or ABSORBENT there is no transmission
                # except if n is a subsegment
                if (name!='METAL') & (name!='ABSORBENT') or (n in self.lsss):
                    self.Gi.add_node((n,cy0,cy1))
                    self.Gi.add_node((n,cy1,cy0))
                    self.Gi.pos[(n, cy0, cy1)] = tuple(self.Gs.pos[n]+ln*delta/2.)
                    self.Gi.pos[(n, cy1, cy0)] = tuple(self.Gs.pos[n]-ln*delta/2.)

                else :
                    self.Gi.pos[(n, cy0, cy1)] = tuple(self.Gs.pos[n]+ln*delta/2.)
                    self.Gi.pos[(n, cy1, cy0)] = tuple(self.Gs.pos[n]-ln*delta/2.)

        
        #
        # 2) Establishing link between interactions
        #
        for cy in self.Gt.node:
            if cy >0:
                vnodes = self.Gt.node[cy]['polyg'].vnodes
                # indoor = self.Gt.node[cy]['indoor']
                npt=[]
                for x in vnodes:
                    if x < 0:
                        if self.ddiff.has_key(x):
                            for y in self.ddiff[x][0]:
                                if y == cy:
                                    npt.append(x)
                

                nseg = filter(lambda x : x>0,vnodes)
                vnodes = nseg+npt
                for nstr in vnodes:
                    if nstr in self.Gv.nodes():
                        if nstr>0:
                            cyo1 = self.Gs.node[nstr]['ncycles']
                            cyo1 = filter(lambda x : x!=cy,cyo1)[0]
                            # R , Tin , Tout
                            if cyo1>0:
                                if (nstr,cy) in self.Gi.nodes():
                                    li1 = [(nstr,cy),(nstr,cy,cyo1),(nstr,cyo1,cy)]
                                else:# no reflection on airwall
                                    li1 = [(nstr,cyo1,cy)]
                            else:
                                if (nstr,cy) in self.Gi.nodes():
                                    li1 = [(nstr,cy)]
                                else:
                                    li1 =[]
                        else:
                            # D
                            li1 =[(nstr,)]
                        # list of cycle entities in visibility of nstr
                        lneighb  = nx.neighbors(self.Gv,nstr)
                        lneighcy = filter(lambda x: x in vnodes,lneighb)

                        for nstrb in lneighcy:
                            if nstrb in self.Gv.nodes():
                                if nstrb>0:
                                    cyo2 = self.Gs.node[nstrb]['ncycles']
                                    cyo2 = filter(lambda x : x!=cy,cyo2)[0]
                                    if cyo2>0:
                                        if (nstrb,cy) in self.Gi.nodes():
                                            li2 = [(nstrb,cy),(nstrb,cy,cyo2),(nstrb,cyo2,cy)]
                                        else: #no reflection on airwall
                                            li2 = [(nstrb,cy,cyo2),(nstrb,cyo2,cy)]
                                    else:
                                        if (nstrb,cy) in self.Gi.nodes():
                                            li2 = [(nstrb,cy)]
                                else:
                                    li2 = [(nstrb,)]

                                for i1 in li1:
                                    #print li1
                                    for i2 in li2:
                                        #print li2
                                        if (i1[0]!=i2[0]):
                                            if ((len(i1)==2) & (len(i2)==2)):
                                                #print "RR"
                                                self.Gi.add_edge(i1,i2)
                                                self.Gi.add_edge(i2,i1)
                                            if ((len(i1)==2) & (len(i2)==3)):
                                                #print "RT"
                                                if i1[1]==i2[1]:
                                                    self.Gi.add_edge(i1,i2)
                                            if ((len(i1)==3) & (len(i2)==2)):
                                                #print "TR"
                                                if i1[2]==i2[1]:
                                                    self.Gi.add_edge(i1,i2)
                                            if ((len(i1)==3) & (len(i2)==3)):
                                                #print "TT"
                                                if i1[2]==i2[1]:
                                                    self.Gi.add_edge(i1,i2)
                                                if i2[1]==i1[2]:
                                                    self.Gi.add_edge(i2,i1)
                                            if ((len(i1)==1) & (len(i2)==3)):
                                                #print "DT"
                                                if  i2[1]==cy:
                                                    self.Gi.add_edge(i1,i2)
                                            if ((len(i1)==3) & (len(i2)==1)):
                                                #print "TD"
                                                if  i1[2]==cy:
                                                    self.Gi.add_edge(i1,i2)
                                            if ((len(i1)==1) & (len(i2)==2)):
                                                #print "DR"
                                                self.Gi.add_edge(i1,i2)
                                            if ((len(i1)==2) & (len(i2)==1)):
                                                #print "RD"
                                                self.Gi.add_edge(i1,i2)
                                            if ((len(i1)==1) & (len(i2)==1)):
                                                #print "DD"
                                                self.Gi.add_edge(i1,i2)




        # updating the list of interactions of a given cycle
        for c in self.Gt.node:
            if c != 0:
                vnodes = self.Gt.node[c]['polyg'].vnodes
                for k in npt:
                    self.Gt.node[c]['inter']+= [(k,)]

    def filterGi(self, situ = 'outdoor'):
        """ Filter Gi to manage indoor/ outdoor situations

        Not called

        """

        # get outdoor notes
        cy = np.array(self.Gt.nodes())
        uout = np.where([not self.Gt.node[i]['indoor'] for i in cy])
        cyout = cy[uout]

        inter = self.Gi.nodes()
        Ti = [i for i in inter if ((len(i) == 3) and i[0]>0)]
        Ri = [i for i in inter if ((len(i) == 2) and i[0]>0)]
        Di = [i for i in inter if i[0]<0]

        Ti = [i for i in Ti if ((i[1] in cyout) and (i[2] in cyout))]
        Ri = [i for i in Ri if (i[1] in cyout) ]
        Di = [i for i in Di if (i in self.ldiffout) ]

        rinter = Ti + Ri + Di

        rGi = nx.subgraph(self.Gi,rinter)
        rGi.pos={i:self.Gi.pos[i] for i in self.Gi.nodes()}

        self.Gi = rGi
        self.Gi.pos = rGi.pos





    def outputGi(self):
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
        # loop over all edges of Gi
        Nedges = len(self.Gi.edges())
        #print "Gi Nedges :",Nedges
        for k,e in enumerate(self.Gi.edges()):
            #if (k%100)==0:
            #print "edge :  ",k
            # extract  both termination interactions nodes

            i0 = e[0]
            i1 = e[1]

            nstr0 = i0[0]
            nstr1 = i1[0]

            # list of authorized outputs. Initialized void
            output = []

            # nstr1 : segment number of central interaction
            if nstr1>0:
                # central interaction is a segment
                pseg1 = self.seg2pts(nstr1).reshape(2,2).T
                # create a Cone object
                cn = cone.Cone()
                # if starting from segment
                if nstr0>0:
                    pseg0 = self.seg2pts(nstr0).reshape(2,2).T
                    # if nstr0 and nstr1 are connected segments
                    if (len(np.intersect1d(nx.neighbors(self.Gs,nstr0),nx.neighbors(self.Gs,nstr1)))==0):
                        # from 2 not connected segment
                        cn.from2segs(pseg0,pseg1)
                    else:
                        # from 2 connected segments
                        cn.from2csegs(pseg0,pseg1)
                # if starting from a point
                else:
                    pt = np.array(self.Gs.pos[nstr0])
                    cn.fromptseg(pt,pseg1)

                # list all potential successors of interaction i1
                i2 = nx.neighbors(self.Gi,i1)
                ipoints = filter(lambda x: len(x)==1,i2)
                pipoints = np.array([self.Gs.pos[ip[0]] for ip in ipoints]).T
                # filter tuple (R | T)
                #istup = filter(lambda x : type(eval(x))==tuple,i2)
                # map first argument segment number
                #isegments = np.unique(map(lambda x : eval(x)[0],istup))
                isegments = np.unique(filter(lambda y: y>0,map(lambda x: x[0],i2)))
                # if nstr0 and nstr1 are adjescent segment remove nstr0 from
                # potential next interaction
                if len(np.intersect1d(self.Gs.neighbors(nstr0),self.Gs.neighbors(nstr1)))>0:
                       isegments = np.array(filter(lambda x : x!=nstr0,isegments))
                # there are one or more segments
                if len(isegments)>0:
                    points = self.seg2pts(isegments)
                    pta = points[0:2,:]
                    phe = points[2:,:]
                    # add difraction points
                    # WARNING Diffraction points are added only if a segment is seen
                    # it should be the case in 99% of cases

                    if len(ipoints)>0:
                        isegments = np.hstack((isegments,np.array(ipoints)[:,0]))
                        pta = np.hstack((pta,pipoints))
                        phe = np.hstack((phe,pipoints))

                    #cn.show()

                    # if i0 == (38,79) and i1 == (135,79,23):
                    #     print i0,i1
                    #     import ipdb
                    #     ipdb.set_trace()
                    # i1 : interaction T
                    if len(i1)==3:
                        typ,prob = cn.belong_seg(pta,phe)
                        #if bs.any():
                        #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                        #if ~bs.any():
                        #    plu.displot(pta[:,~bs],phe[:,~bs],color='k')

                    # i1 : interaction R --> mirror
                    if len(i1)==2:
                        Mpta = geu.mirror(pta,pseg1[:,0],pseg1[:,1])
                        Mphe = geu.mirror(phe,pseg1[:,0],pseg1[:,1])
                        typ,prob = cn.belong_seg(Mpta,Mphe)
                        #print i0,i1
                        #if ((i0 == (6, 0)) & (i1 == (7, 0))):
                        #    pdb.set_trace()
                        #if bs.any():
                        #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                        #if ~bs.any():
                        #    plu.displot(pta[:,~bs],phe[:,~bs],color='m')
                        #    plt.show()
                        #    pdb.set_trace())
                    ########
                    #SOMETIMES PROBA IS 0 WHEREAS SEG IS SEEN
                    ###########
                    # # keep segment with prob above a threshold
                    # isegkeep = isegments[prob>0]
                    # # dict   {numint : proba}
                    # dsegprob = {k:v for k,v in zip(isegkeep,prob[prob>0])}
                    ########### 4 lines are replaced by
                    # keep segment with prob above a threshold
                    utypseg = typ!=0
                    isegkeep = isegments[utypseg]
                    # dict   {numint : proba}
                    dsegprob = {k:v for k,v in zip(isegkeep,prob[utypseg])}
                    #########
                    output = filter(lambda x : x[0] in isegkeep, i2)
                    probint = map(lambda x: dsegprob[x[0]],output)
                    # dict interaction : proba
                    dintprob = {k:v for k,v in zip(output,probint)}

                    # keep all segment above nstr1 and in Cone if T
                    # keep all segment below nstr1 and in Cone if R

            else:
                # central interaction is a point

                # 1) Simple approach
                #       output is all visible interaction
                # 2) TO BE DONE
                #
                #       output of the diffraction points
                #       exploring
                #b
                #          + right of ISB
                #          + right of RSB
                #
                #  + using the wedge cone
                #  + using the incident cone
                #

                output = nx.neighbors(self.Gi,(nstr1,))
                nout = len(output)
                probint = np.ones(nout) # temporarybns
                dintprob = {k:v for k,v in zip(output,probint)}


            self.Gi.add_edge(i0,i1,output=dintprob)

    def intercy(self,ncy,typ='source'):
        """ return the list of interactions seen from a cycle

        Parameters
        ----------

        ncy : cycle number( Project -> save proejct)
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
        lD = filter(lambda x: len(x)==1,lint)
        lR = filter(lambda x: len(x)==2,lint)
        lT = filter(lambda x: len(x)==3,lint)

        # visible R|T source cycle is ncy

        lR = filter(lambda x : x[1]==ncy,lR)
        if typ=='source':
            lT = filter(lambda x: x[1]==ncy,lT)
        if typ=='target':
            lT = filter(lambda x: x[2]==ncy,lT)
        if typ=='all':
            lT=lT
        # Finding the diffraction points
        # Diffraction points are different from indoor cycle and outdoor
        # cycles
        #
        # TODO check wedge validity.
        #
        vnodes = self.Gt.node[ncy]['polyg'].vnodes
        vpoints = filter(lambda x: x<0,vnodes)
        indoor = self.Gt.node[ncy]['indoor']
        if indoor:
            lD = map(lambda y : (y,),filter(lambda x : x in
                                            self.ldiffin,vpoints))
        else:
            lD = map(lambda y : (y,),filter(lambda x : x in
                                            self.ldiffout,vpoints))

        return lR,lT,lD


    # def intercyGc2Gt(self,ncy,typ='source'):
    #     """ return the list of interactions in Gt from a Gc cycle

    #     Parameters
    #     ----------

    #     ncy : cycle number from Gc
    #     typ : string
    #         if 'source' connect source cycle
    #         if 'target' connect target cycle

    #     """

    #     # list of interactions

    #     lint = self.Gi.node
    #     lTa=[]
    #     lRa=[]
    #     lDa=[]
    #     if self.Gt.node[ncy].has_key('merged'):
    #         cym = self.Gt.node[ncy]['merged']
    #         lcy = self.Gc.node[cym]['merged']

    #     else :
    #         lcy=[ncy]
    #     # lint = self.Gi.node
    #     for c in lcy:
    #         # list of tuple interactions (R|T)
    #         lD = filter(lambda x: len(x)==1,lint)
    #         lR = filter(lambda x: len(x)==2,lint)
    #         lT = filter(lambda x: len(x)==3,lint)

    #         # visible R|T source cycle is ncy

    #         lR = filter(lambda x : x[1]==c,lR)
    #         if typ=='source':
    #             lT = filter(lambda x: x[1]==c,lT)
    #         if typ=='target':
    #             lT = filter(lambda x: x[2]==c,lT)


    #         # Finding the diffraction points
    #         # Diffraction points are different from indoor cycle and outdoor
    #         # cycles
    #         #
    #         # TODO check wedge validity.
    #         #
    #         vnodes = self.Gt.node[ncy]['polyg'].vnodes
    #         vpoints = filter(lambda x: x<0,vnodes)
    #         indoor = self.Gt.node[ncy]['indoor']
    #         if indoor:
    #             lD = map(lambda y : (y,),filter(lambda x : x in
    #                                             self.ldiffin,vpoints))
    #         else:
    #             lD = map(lambda y : (y,),filter(lambda x : x in
    #                                             self.ldiffout,vpoints))
    #         lTa.extend(lT)
    #         lRa.extend(lR)
    #         lDa.extend(lD)
    #     return lRa,lTa,lDa


#    def showGraph(self,**kwargs):
#        """
#        Parameters
#        ----------
#        print n,nb
#        show : boolean
#        fig
#        nodes
#        eded
#        ndnd
#        nded
#        linewidth
#        roomlis
#        """
#        defaults={'show':False,
#                  'fig':[],
#                  'ax':[],
#                  'nodes':False,
#                  'eded':True,
#                  'ndnd':True,
#                  'nded':True,
#                  'linewidth':2,
#                  'nodelist':[]
#                  }
#        for key, value in defaults.items():
#            if kwargs.has_key(key):
#                setattr(self, key, kwargs[key])
#            else:
#                setattr(self, key, value)
#                kwargs[key]=value
#        if kwargs['fig']==[]:
#            fig = plt.figure()
#            fig.set_frameon(True)
#        else:
#            fig = kwargs['fig']
#        if kwargs['ax']==[]:
#            ax = fig.gca()
#        else:
#            ax = kwargs['ax']
#
#        if graph=='t':
#            G = self.Gt
#            nx.draw(G,G.pos)
#        if graph=='r':
#            G = self.Gr
#            nx.draw(G,G.pos)
#        if graph=='s':
#            G = self.Gs
#            nx.draw(G,G.pos)
#        if graph=='v':
#            G = self.Gv
#            nx.draw(G,self.Gs.pos)
#
#        for k,ncy in enumerate(self.Gt.node.keys()):
#            self.Gt.node[ncy]['polyg'].plot()
#        ax.axis('scaled')
#        # Display doors and windows
#        d = self.subseg()
#        for ss in d.keys():
#            if ss=='DOOR':
#                color='red'
#            if ss=='3D_WINDOW_GLASS':
#                color='blue'
#            if ss=='WINDOW_GLASS':
#                color='cyan'
#            for ns in d[ss]:
#                np1,np2 = self.Gs.neighbors(ns)
#                x  = [self.Gs.pos[np1][0],self.Gs.pos[np2][0]]
#                y  = [self.Gs.pos[np1][1],self.Gs.pos[np2][1]]
#                ax.plot(x,y,linewidth=2,color=color)
#        if kwargs['show']:
#            plt.show()

    def show(self,**kwargs):
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
                    'node_color':'w',
                    'edge_color':'k',
                    'node_size':20,
                    'font_size':30,
                    'nodelist': [],
                    'figsize': (5,5),
                    'mode':'cycle',
                    }
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        segfilt = filter(lambda x : x not in self.name['AIR'], self.tsg)
        # get the association between segment and nx edges
        edges = self.Gs.edges()
        Ne = len(edges)
        segments = np.array(edges)[:,0]
        dse = {k:v for k,v in zip(segments,range(Ne))}

        edfilt = list(np.ravel(np.array(map(lambda x : [dse[x]-1,dse[x]],segfilt))))

        # edgelist is to be understood as edge of graph and not segments of layout

        fig,ax = self.showG('s',nodes=False,edgelist=edfilt)

        # display degree 1 nodes
        if 1 in self.degree:
            ldeg1  = list(self.degree[1])
            print ldeg1
            fig,ax = self.showG('s',
                fig=fig,
                ax=ax,
                nodelist=ldeg1,
                edges=kwargs['edges'],
                nodes=kwargs['nodes'],
                node_size=kwargs['node_size'],
                node_color='r')

        # display degree 4 nodes
        if 4 in self.degree:
            ldeg4  = list(self.degree[4])
            fig,ax = self.showG('s',
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
        u""" show the different graphs

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
                    'overlay':False
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
                    'nodes': False,
                    'edges': True,
                    'sllist':[],
                    'airwalls': False,
                    'aw': [],
                    'subseg': False,
                    'slab': False,
                    'labels': False,
                    'alphan': 1.0,
                    'alphae': 1.0,
                    'width': 2,
                    'node_color':'w',
                    'edge_color':'',
                    'node_size':20,
                    'font_size':15,
                    'nodelist': [],
                    'edgelist': [],
                    'figsize': (5,5),
                    'mode':'nocycle',
                    'alphacy':0.8,
                    'colorcy':'#abcdef',
                    'lvis' : ['nn','ne','ee'],
                    'linter' : ['RR','TT','RT','TR','RD','DR','TD','DT','DD'],
                    'show0':False,
                    'axis':False,
                    'overlay':False
                    }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value
        if kwargs['aw'] != []:
            kwargs['airwalls']=kwargs['aw']
        # overriding first argument graph
        if 'graph' in kwargs:
            graph = kwargs['graph']

        # get color dictionnary from pyutil

        cold = pyu.coldict()

        if isinstance(kwargs['labels'],list):
            labels = kwargs['labels']
        elif kwargs['labels'] == True:
            labels=['s','t','v','i','w']
        elif isinstance(kwargs['labels'],str):
            labels=kwargs['labels']
        else:
            labels=[]
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

            try :
                sllist = [kwargs['sllist'].pop()]
            except:
                sllist = self.name.keys()

            for lmat in sllist:
                lseg = self.name[lmat]
                if lseg != []:
                    lseg2 = [np.where(np.array(self.Gs.edges()) == i)[0] for i in lseg]
                    kwargs['edgelist'] = reduce(lambda x,y:list(x)+list(y),lseg2)
                    if kwargs['slab']:
                        kwargs['edge_color']=cold[self.sl[lmat]['color']]
                        kwargs['width']=self.sl[lmat]['linewidth']
                    else:
                        kwargs['edge_color']='k'
                        kwargs['width']=1

                if 's' in labels:
                   kwargs['labels']=True
                else:
                   kwargs['labels']=False
                kwargs['fig'],kwargs['ax'] = gru.draw(G,**kwargs)

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
                        x = np.array([self.Gs.pos[np1][0], self.Gs.pos[np2][0]])
                        y = np.array([self.Gs.pos[np1][1], self.Gs.pos[np2][1]])
                        xoff = (1+ns[1])*0.05*norm[0]
                        yoff = (1+ns[1])*0.05*norm[1]
                        kwargs['ax'].plot(x+xoff, y+yoff, linewidth=2, color=color)


        #
        # t : graph of cycles
        #
        if 't' in graph:
            G = self.Gt
            if not kwargs['show0']:
                # filter out the 0 cycle
                nodes = G.nodes()
                edges = G.edges()
                nodf = filter(lambda x : x!=0,nodes)
                edf  = filter(lambda x: ((edges[x][0]!=0) & (edges[x][1]!=0)),np.arange(len(edges)))
                kwargs['nodelist']=nodf
                kwargs['edgelist']=edf
            else :
                kwargs['nodelist'] = G.nodes()
                kwargs['edgelist'] = np.arange(len(G.edges()))

            if kwargs['edge_color']=='':
                kwargs['edge_color'] ='r'
            if 't' in labels:
                kwargs['labels']=True
            else:
                kwargs['labels']=False
            fig,ax = gru.draw(G,**kwargs)
            kwargs['fig']=fig
            kwargs['ax']=ax
        #
        # r : graph of rooms
        #
        if 'r' in graph:
            G = self.Gr
            if kwargs['edge_color']=='':
                kwargs['edge_color'] ='g'

            kwargs['fig'],kwargs['ax'] = gru.draw(self.Gs,
                              nodes=False,edges=True,alphacy=1.,
                              fig=kwargs['fig'],ax=kwargs['ax'],labels=False)
            if 'r' in labels:
                kwargs['labels']=True
            else:
                kwargs['labels']=False
            fig,ax = gru.draw(G,**kwargs)
            kwargs['fig']=fig
            kwargs['ax']=ax
                #
        # v : visibility graph
        # In blue : segment segment
        # In red  : point point (Diffraction)
        # In green : point segment
        #
        if 'v' in graph:

            G = self.Gv
            G.pos={}
            # nodes of Gv are nodes of Gs
            G.pos.update(self.Gs.pos)

            if kwargs['edge_color']=='':
                kwargs['edge_color'] ='m'

            edges = G.edges()
            rle = range(len(edges))
            eded  = filter(lambda x: (edges[x][0]>0) & (edges[x][1]>0),rle)
            ndnd  = filter(lambda x: (edges[x][0]<0) & (edges[x][1]<0),rle)
            nded  = filter(lambda x: (((edges[x][0]<0) & (edges[x][1]>0)) |
                                     ((edges[x][0]>0) & (edges[x][1]<0))),rle)
            if 'v' in labels:
                kwargs['labels']=True
            else:
                kwargs['labels']=False

            if 'ee' in kwargs['lvis']:
                kwargs['edgelist'] = eded
                kwargs['edge_color']='blue'
                kwargs['node_size'] = 200
                kwargs['fig'],kwargs['ax'] = gru.draw(G,**kwargs)
            if 'nn' in kwargs['lvis']:
                kwargs['edgelist'] = ndnd
                kwargs['edge_color'] ='red'
                kwargs['fig'],kwargs['ax'] = gru.draw(G,**kwargs)
            if 'ne' in kwargs['lvis']:
                kwargs['edgelist'] = nded
                kwargs['edge_color']='green'
                kwargs['fig'],kwargs['ax'] = gru.draw(G,**kwargs)
        #
        # i :  interaction graph
        #
        if 'i' in graph:

            G = self.Gi

            if kwargs['edge_color']=='':
                kwargs['edge_color'] ='k'

            #
            # Parsing the type of interactions
            #

            edges = G.edges()

            # range len edges

            rle = range(len(edges))

            DD = filter(lambda x:  ((len(edges[x][0])==1) &
                                    (len(edges[x][1])==1)),rle)

            RR = filter(lambda x: ((len(edges[x][0])==2) &
                                   (len(edges[x][1])==2)),rle)

            TT = filter(lambda x: ((len(edges[x][0])==3) &
                                   (len(edges[x][1])==3)),rle)

            RT = filter(lambda x: ((len(edges[x][0])==2) &
                                   (len(edges[x][1])==3)),rle)

            TR = filter(lambda x: ((len(edges[x][0])==3) &
                                   (len(edges[x][1])==2)),rle)

            RD = filter(lambda x:  ((len(edges[x][0])==2) &
                                    (len(edges[x][1])==1)),rle)

            TD = filter(lambda x:  ((len(edges[x][0]) ==3) &
                                    (len(edges[x][1]) ==1)),rle)

            DR = filter(lambda x:  ((len(edges[x][0]) ==1) &
                                    (len(edges[x][1])==2)),rle )

            DT = filter(lambda x:  ((len(edges[x][0]) ==1) &
                                    (len(edges[x][1]) ==3)),rle )

            tabcol = ['b','g','r','m','c','orange','purple','maroon','purple','k'][::-1]
            li = []
            if 'i' in labels:
                kwargs['labels']=True
            else:
                kwargs['labels']=False
            for inter in kwargs['linter']:
                if len(eval(inter))>0:
                    li.append(inter)
                    kwargs['edgelist'] = eval(inter)
                    # ndlist = map(lambda x: edges[x][0],kwargs['edgelist'])+\
                    #          map(lambda x: edges[x][1],kwargs['edgelist'])
                    ndlist = map(lambda x: edges[x][0],kwargs['edgelist'])+\
                             map(lambda x: edges[x][1],kwargs['edgelist'])
                    # keep only unique interaction
                    unique = []
                    [unique.append(it) for it in ndlist if it not in unique]

                    kwargs['nodelist'] = unique
                    kwargs['edge_color'] = tabcol.pop()
                    kwargs['fig'],kwargs['ax'] = gru.draw(G,**kwargs)
            legtxt = ['Gs'] + li
            #plt.legend(legtxt)
        #
        # w :  waypoint graph
        #
        if 'w' in graph:

            G = self.Gw

            if kwargs['edge_color']=='':
                kwargs['edge_color'] ='k'
            kwargs['fig'],kwargs['ax'] = gru.draw(self.Gs,
                              nodes=False,edges=True,alphacy=1.,
                              fig=kwargs['fig'],ax=kwargs['ax'],labels=False)
            if 'w' in labels:
                kwargs['labels']=True
            else:
                kwargs['labels']=False
            fig,ax = gru.draw(G,**kwargs)
            kwargs['fig']=fig
            kwargs['ax']=ax

        #args = {'fig':kwargs['fig'],'ax':kwargs['ax'],'show':False}

        if len(kwargs['edgelist'])==0:
            if kwargs['mode']=='cycle':
                for k, ncy in enumerate(self.Gt.node.keys()):
                    fig,ax = self.Gt.node[ncy]['polyg'].plot(alpha=kwargs['alphacy'],color=kwargs['colorcy'],**args)
                    args['fig']=fig
                    args['ax']=ax
            if kwargs['mode']=='room':
                for k, nro in enumerate(self.Gr.node.keys()):
                    fig,ax = self.Gr.node[nro]['cycle'].show(**args)
                    args['fig']=fig
                    args['ax']=ax

        kwargs['ax'].axis('scaled')
        if not kwargs['axis']:
            kwargs['ax'].axis('off')


        if kwargs['overlay']:
            imok = False
            if self.display['overlay_file']!='':
                image = Image.open(os.path.join(basename,pstruc['DIRIMAGE'],self.display['overlay_file']))
                imok =True
            if imok:
                if 'v' in self.display['overlay_flip']:
                    print "flip v"
                    image = image.transpose(Image.FLIP_LEFT_RIGHT)
                if 'h' in self.display['overlay_flip']:
                    image = image.transpose(Image.FLIP_TOP_BOTTOM)
                    print "flip h"
                kwargs['ax'].imshow(image, extent=self.display['overlay_axis'],alpha=self.display['alpha'],origin='lower')


        if kwargs['show']:
            plt.show()

        return kwargs['fig'],kwargs['ax']

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
            ax  = fig.gca()
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
        return(rooms,[tuple(self.Gw.pos[i]) for i in rooms])

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

    def ptin(self,pt=np.array((0, 0, 0))):
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

        x= np.array((self.ax[:2]))
        y= np.array((self.ax[2:]))

        # being   in [xmin xmax]
        c0 = pt[0]<x[1] and  pt[0]>x[0]
        # being   in [ymin ymax]
        c1 = pt[1]<y[1] and  pt[1]>y[0]


        return (c0 & c1)


    def ptGs2cy(self,n=-1):
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
        if n >0:
            return self.Gs.node[n]['ncycles']
        else:
            nseg = self.Gs[n].keys()
            cy=[]
            for nn in nseg :
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
            if ncy>0:
                criter1 = self.Gt.node[ncy]['polyg'].touches(ptsh)
                criter2 = self.Gt.node[ncy]['polyg'].contains(ptsh)
                if (criter1 or criter2):
                    cycle_exists = True
                    return(ncy)
        if not cycle_exists:
            raise NameError(str(pt)+" is not in any cycle")

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
            pt=np.hstack((pt,h))
            return(pt)
        else:
            raise NameError("cycle "+str(cy)+" not in self.Gt")




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
            raise NameError(str(pt)+" is not in any room")

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
            #if seg in self.Gt.node[self.Gr.node[nr]['cycle']]['vnodes']:
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
            raise NameError(str(room)+" is not in not on Gr")
        u = np.where(seg>=0)
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
            raise NameError(str(room)+" is not in not on Gr")
        u = np.where(nod<0)
        nod = nod[u]

        return np.sort(nod.tolist())


    def builGr2(self,v):
        """ alternative buildGr method (deprecated)


        parameters
        ----------
        v

        """
        # 1 : Cycles which are connected with an airwall are merged
        # 2 : The remaining cycles which contains at least one common door are
        # connected

        self.Gr = nx.Graph()
        self.Gr.pos = {}
        for ncy in self.Gt.nodes:
            pass
        a = map(lambda x: self.Gs.node[x]['ncycles'],v)
        b = reduce(lambda x,y: x+y,a)
        involvedcycles = np.unique(np.array(b))
        # list of cycles which are already involved in rooms
        alreadythere = filter(lambda x: x in cycleroom.keys(),involvedcycles)


    # def buildGc(self):
    #     """ build the graph of cycles

    #     Returns
    #     -------

    #     Ga : graph of adjascent rooms

    #     Notes
    #     -----

    #     An adjascent graph (connecting cycle sharing airwalls) is created.
    #     This graph contains multiples sub-graph which are not connected with each others

    #     lGa is a list of each of thoose subgraphs

    #     For each of those sub graphs, all their nodes are merged until they only have a single node.
    #     The main difficulty here , is to merge the polygons of those nodes.
    #     Because GeomUtil.Polygon, does not support multiple Polygons, nodes of a graph Ga must be 
    #     merged into a specific order.
    #     This exaplains the dfs_successors and following in the code !:w



    #     """

    #     #
    #     # Create a graph of adjascent cycles
    #     #

    #     Ga = nx.Graph()
    #     Ga.pos ={}
    #     for k in self.Gt.edge:
    #         dk = self.Gt.edge[k]
    #         for cy in dk:
    #             try:
    #                 segs = dk[cy]['segment']
    #             except:
    #                 segs=[]
    #             for s in segs:
    #                 if self.Gs.node[s]['name']=='AIR':
    #                     if k not in Ga.node:
    #                         Ga.add_node(k)
    #                         Ga.pos[k]=self.Gt.pos[k]
    #                     if cy not in Ga.node:
    #                         Ga.add_node(cy)
    #                         Ga.pos[cy]=self.Gt.pos[cy]
    #                     Ga.add_edge(k,cy)

    #     # deep copy of Gt
    #     self.Gc = copy.deepcopy(self.Gt)

    #     # list of connected subgraphs of Gt
    #     lGa = nx.connected_component_subgraphs(Ga)

    #     #dictionnary mapping
    #     # key = old cycle ID :  value =new cycle ID
    #     dmap={x:x for x in self.Gc.nodes()}


    #     def merge_cycle(root,child):
    #         """ method to merge a child cycle into a root cycle
    #         """
    #         try:
    #             self.Gc.node[root]['polyg']+=self.Gc.node[child]['polyg'] # here the merging
    #             self.Gc.node[root]['cycle']+=self.Gc.node[child]['cycle'] # here the merging
    #             # import ipdb
    #             # ipdb.set_trace()
    #             # self.Gc.node[root]['polyg'].plot(fig=plt.gcf(),ax=plt.gca(),color='g')
    #             # self.Gc.node[root]['cycle'].show(fig=plt.gcf(),ax=plt.gca())
    #             # plt.draw()
    #         except:
    #             import ipdb
    #             ipdb.set_trace()
    #         try:
    #             self.Gc.node[root]['merged'].append(child)
    #         except:
    #             self.Gc.node[root]['merged']=[child]

    #         dmap[child] = root

            
    #         self.Gc.add_edges_from([(root,x) for x in self.Gc.edge[child].keys()])
    #         # update Gt
    #         self.Gt.node[child]['merged'] = root
    #         # a merged cycle is open
    #         self.Gt.node[child]['isopen'] = True
    #         self.Gt.node[root]['merged'] = root
    #         self.Gc.remove_node(child)
    #         # self.Gc.pos[root]=tuple(self.Gc.node[root]['cycle'].g)
    #         self.Gc.pos[root]=tuple(np.array(self.Gc.node[root]['polyg'].centroid.xy)[:,0])


    #     for Ga in lGa:


    #         # root node of subgraph
    #         r = Ga.nodes()[0]
    #         cnctd = [r]
    #         # depth first search successors tree rooted on r
    #         # this return dict where keys are cycles and values
    #         # are list of cycle (a.k.a. child ) connected to the key cycle
    #         dn = nx.dfs_successors(Ga,r)

    #         # merge 

    #         # loop on each root cycles
    #         for uk,k in enumerate(dn.keys()):
    #             # loop on child cycles
    #             for v in dn[k]:
    #                 # 
    #                 root =dmap[k]
    #                 while dmap[root] !=root :
    #                     root = dmap[root]
    #                 # if the child cycle has already been merged in a previous
    #                 # loop; its ID has changed but thoose of its childs too. 
    #                 # This update the ID of childs of v
    #                 if v in dn.keys()[:uk]:
    #                     # list of cycle to be modified
    #                     lcy = dn[v]
    #                     [dmap.update({c:root}) for c in lcy]

    #                 merge_cycle(root,v)
    #     # find diffractions of the layout
    #     self._find_diffractions()





    # def buildGc_old(self):
    #     """ 
    #     DEPRECATED

    #     This version of buildGc do not works well with
    #     large layout. 
    #     Indeed, the algorithm cannot merge large room
    #     with multiple airwalls

    #     build the graph of cycles



    #     Returns
    #     -------

    #     Ga : graph of adjascent rooms

    #     Notes
    #     -----

    #     adjascent rooms are connected
    #     Gr is at first a deep copy of Gt

    #     The difficulty here is to take into account the AIR transition
    #     segments

    #     """
    #     #
    #     # Create a graph of adjascent cycles
    #     #

    #     Ga = nx.Graph()
    #     Ga.pos ={}
    #     for k in self.Gt.edge:
    #         dk = self.Gt.edge[k]
    #         for cy in dk:
    #             try:
    #                 segs = dk[cy]['segment']
    #             except:
    #                 segs=[]
    #             for s in segs:
    #                 if self.Gs.node[s]['name']=='AIR':
    #                     if k not in Ga.node:
    #                         Ga.add_node(k)
    #                         Ga.pos[k]=self.Gt.pos[k]
    #                     if cy not in Ga.node:
    #                         Ga.add_node(cy)
    #                         Ga.pos[cy]=self.Gt.pos[cy]
    #                     Ga.add_edge(k,cy)

    #     # deep copy of Gt
    #     self.Gc = copy.deepcopy(self.Gt)

    #     # list of connected subgraphs of Gt
    #     lGa = nx.connected_component_subgraphs(Ga)
    #     connected = []
    #     for Ga in lGa:
    #         # root node of subgraph
    #         r = Ga.nodes()[0]
    #         cnctd = [r]
    #         # depth first search successors tree rooted on r
    #         dn = nx.dfs_successors(Ga,r)
    #         Nlevel = len(dn)

    #         #
    #         # Probably it exists a simpler manner to obtain
    #         # the sequence of connected nodes
    #         #
    #         succ =[]

    #         while dn.keys()!=[]:
                
    #             try:
    #                 succ = succ+ dn.pop(r)
    #             except:
    #                 r = dn.pop()
    #             if r ==102:
    #                 import ipdb
    #                 ipdb.set_trace()


    #             for i in range(len(succ)):
    #                 try:
    #                     n = succ.pop()
    #                 except:
    #                     r = succ.pop()
    #                     cnctd.append(r)
    #                     break
    #                 if n in dn.keys():
    #                     r = n
    #                     cnctd.append(r)
    #                 else:
    #                     cnctd.append(n)
    #                     # try:
    #                     #     r = succ.pop()
    #                     #     cnctd.append(r)
    #                     # except:
    #                     #     break
    #         #for i in range(Nlevel):
    #         #    succ = dn.pop(r)
    #         #    cnctd = cnctd + succ
    #         #    for k in dn:
    #         #        if k in succ:
    #         #            r = k
    #         #            break
    #         connected.append(cnctd)
    #     #
    #     # Merge all air-connected cycles
    #     #  for all conected components
    #     #  example licy = [[22,78,5],[3,4]] 2 cycles are connected
    #     #
    #     import ipdb
    #     ipdb.set_trace()
    #     merge2c=[]
    #     for licy in connected:
    #         root = licy[0]      # pick the first cycle as root
    #         merged = [root]     # merged cycle is void
    #         tomerge = licy[-1:0:-1]  # pick the inverse remaining part as tomerge list
    #         #for cy in tomerge: #
    #         while tomerge!=[]:

    #             ncy = tomerge.pop()
    #             #print "ncy = ",ncy
    #             # testing cycle contiguity before merging

    #             try:
    #                 croot = self.Gc.node[root]['cycle']
    #             except:
    #                 import ipdb
    #                 ipdb.set_trace()
    #             cy = self.Gc.node[ncy]['cycle']
    #             flip,path = croot.intersect(cy)

    #             if len(path) < 1 and [root,ncy] not in merge2c:
    #                 print tomerge,root,ncy,merge2c
    #                 tomerge.insert(0,ncy)
    #             else:
    #                 print 'merge' 

    #                 neigh = nx.neighbors(self.Gc,ncy) # all neighbors of 5
    #                 self.Gc.node[root]['polyg']+=self.Gc.node[ncy]['polyg'] # here the merging
    #                 self.Gc.node[root]['cycle']+=self.Gc.node[ncy]['cycle'] # here the merging
    #                 merged.append(ncy)
    #                 merge2c.append([root,ncy])
    #                 #print self.Gc.node[root]['polyg'].exterior.xy
    #                 for k in neigh:
    #                     if k!= root:
    #                         self.Gc.add_edge(root,k)

    #         # keep track of merged convex cycles
    #         self.Gc.node[root]['merged'] = merged
    #         for cy in merged:
    #             self.Gt.node[cy]['merged'] = root
    #             # a merged cycle is open
    #             self.Gt.node[cy]['isopen'] = True
    #         self.Gt.node[root]['merged'] = root
    #         # remove merged cycles
    #         for cy in merged:
    #             if cy != root:
    #                 self.Gc.remove_node(cy)
    #         # update pos of root cycle with new center of gravity
    #         self.Gc.pos[root]=tuple(self.Gc.node[root]['cycle'].g)

    #     self._find_diffractions()




    #     return(Ga)



    def _find_diffractions(self,tol=0.01):
        """ Find diffractions points of the Layout

        Parameters
        ----------
        tol : float
            telerance in radians

        Returns
        -------

        Update self.ddiff {nseg : ([ncy1,ncy2],wedge_angle)}

        """
        dangles = self.get_Gt_angles()

        #
        # The canditate points for being diffraction points have degree 1 or 2
        # A point diffracts toward one or several cycles
        #
        #ldiff = list(np.hstack((self.degree[1],self.degree[2])).astype('int'))
        lpnt = [x for x in self.Gs.node if (x <0 and x not in self.degree[0]) ]
        self.ddiff = {}

        for k in lpnt:
            # list of cycles associated with point k
            lcyk = self.Gs.node[k]['ncycles']
            if len(lcyk)>2:
                # Subgraph of connected cycles around k 
                Gtk = nx.subgraph(self.Gt,lcyk)
                # ordered list of connections between cycles
                try:
                    lccyk = nx.find_cycle(Gtk)
                except:
                    pdb.set_trace()
                
                neigh = self.Gs[k].keys()
                sega  = [ n for n in neigh if self.Gs.node[n]['name']=='AIR' ]

                nsector = len(neigh)-len(sega)


                dsector = {i:[] for i in range(nsector)}
                #
                # team building algo
                # 
                ct = 0
                for ccy in lccyk:

                    segsep = self.Gt[ccy[0]][ccy[1]]['segment'][0]
                    typslab = self.Gs.node[segsep]['name']
                    if typslab=='AIR': # same sector
                        dsector[ct].append(ccy[1])
                    else: # change sector
                        ct=(ct+1)%nsector
                        dsector[ct].append(ccy[1])
                        # lcy2.append(ccy[1])
                        # lcy1,lcy2 = lcy2,lcy1


                dagtot = {s:0 for s in range(nsector)} 
                for s in dsector:
                    for cy in dsector[s]:
                        da = dangles[cy]
                        u  = np.where(da[0,:].astype('int')==k)[0][0]
                        dagtot[s] = dagtot[s] + da[1,u]
                
                for s in dagtot:
                    if dagtot[s]>(np.pi+tol):
                        self.ddiff[k]=(dsector[s],dagtot[s])
                        break

                # if agtot1 > (np.pi+tol):
                #     self.ddiff[k]=(lcy1,agtot1)
                # elif 2*np.pi-agtot1 > (np.pi+tol):            
                #     self.ddiff[k]=(lcy2,2*np.pi-agtot1)
            else:
                # diffraction by half-plane detected
                if k in self.degree[1]:
                    self.ddiff[k]=(lcy,2*np.pi)


    def buildGr(self):
        """ build the graph of rooms Gr


        Returns
        -------

        Ga : graph of adjascent rooms

        Notes
        -----

        adjascent rooms are connected
        Gr is at first a deep copy of Gt

        The difficulty here is to take into account the AIR transition
        segments

        """
        
        self.Gr = copy.deepcopy(self.Gt)
        
        # delete node 0 which cannot be a room 

        try:
            del(self.Gr.node[0])
        except:
            pass
        #
        #  Connected components might not be all contiguous
        #  this a problem because the concatenation of cycles
        #  operation below requires cycles contiguity
        #
        for n in self.Gr.nodes():
            self.Gr.node[n]['transition'] = []
        ltrans = self.listtransition
        ldoors = filter(lambda x:self.Gs.node[x]['name']!='AIR',ltrans)

        # Destroy cycles which have no doors

        keys = self.Gr.node.keys()
        for cy in keys:
            if cy>0:
                lseg = self.Gr.node[cy]['cycle'].cycle
                hasdoor = filter(lambda n : n in ldoors,lseg)
                if len(hasdoor)>0:
                    pass
                else:
                    self.Gr.remove_node(cy)
                    self.Gr.pos.pop(cy)

        # Destroy edges which do not share a door
        for e in self.Gr.edges():
            keep = False
            if (e[0]>0) & (e[1]>0):
                cy1 = self.Gr.node[e[0]]['cycle']
                cy2 = self.Gr.node[e[1]]['cycle']
                f,b = cy1.intersect(cy2)
                for s in b:
                    if s>0:
                        if self.Gs.node[s]['transition']:
                            keep = True
                            self.Gr.node[e[0]]['transition'].append(s)
                            self.Gr.node[e[1]]['transition'].append(s)

            if not keep:
                self.Gr.remove_edge(*e)


    # def buildGr3(self):
    #     """ build Graph of rooms

    #     Summary
    #     -------

    #         A room is a set of cycles which contains at least one door

    #         This function requires Gt

    #     """
    #     self.Gr = nx.Graph()
    #     self.Gr.pos = {}
    #     #self.doors ={}
    #     self.transition = {}
    #     self.airwall = {}
    #     d = self.subseg()
    #     # rcpt : rooms counter
    #     rcpt = 0
    #     # ltrans : list of transition segment
    #     ltrans = np.array(self.listtransition)
    #     # lairwalls : list of air walls
    #     lairwalls = filter(lambda x:self.Gs.node[x]['name']=='AIR',ltrans)
    #     # ldoors : list of doors segment number
    #     ldoors = filter(lambda x:self.Gs.node[x]['name']!='AIR',ltrans)
    #     #
    #     # For all cycles
    #     #
    #     # Rule : add a new room if :
    #     #       + the cycle has a transition segment which is not an air wall
    #     #       unless
    #     #       + there already exists a created room which is separated from the
    #     #       current cycle by an airwall
    #     #
    #     #
    #     # roomcycles dict room : list of cycles number involved in room
    #     # cycleroom dict cycle : room number
    #     roomcycles = {}
    #     cycleroom = {}
    #     for k in self.Gt.node:
    #         #if k==5:
    #         #    pdb.set_trace()
    #         # list of segments from the cycle
    #         # which have:
    #         #  ldoors
    #         #  lairwalls
    #         #
    #         lseg = self.Gt.node[k]['cycle'].cycle
    #         u = np.intersect1d(lseg, ldoors)
    #         v = np.intersect1d(lseg, lairwalls)
    #         alreadythere =[]
    #         #
    #         # Analysis of cycles which are connected via an air-wall
    #         #
    #         # cyclehasdoor is True if cycle k has a door segment
    #         # hasdoors is True if at least one adjascent cycle from the same
    #         # room has a door
    #         # doors : np array with doors associated to room k
    #         cyclehasdoor = False
    #         hasdoors = False
    #         doors = np.array([])

    #         if len(u)>0:
    #             cyclehasdoor = True
    #             doors = np.array(u)
    #         # this should be a recursive function
    #         if len(v)>0:
    #             # b : list of cycles which involve an airwall
    #             a = map(lambda x: self.Gs.node[x]['ncycles'],v)
    #             b = reduce(lambda x,y: x+y,a)
    #             involvedcycles = np.unique(np.array(b))
    #             # list of cycles which are already involved in rooms
    #             alreadythere = filter(lambda x: x in cycleroom.keys(),involvedcycles)
    #             notyet = filter(lambda x: x not in cycleroom.keys(),involvedcycles)
    #             for cy1 in involvedcycles:
    #                 lseg1 = self.Gt.node[cy1]['cycle'].cycle
    #                 u1 = np.intersect1d(lseg1, ldoors)
    #                 if len(u1)>0:
    #                     hasdoors = True
    #                     doors = np.unique(np.hstack((doors,u1)))

    #             print "cycle "+str(k)
    #             #print "airwall segment "+str(v)
    #             #print "involved cycles "+str(involvedcycles)
    #             print "already there "+str(alreadythere)
    #             print "not there yet "+ str(notyet)
    #             #print "cycles involved ",cycleroom.keys()
    #         #
    #         # If cycle has a door (transition which is not an air wall)
    #         # Then create a new room
    #         #
    #         if (cyclehasdoor|hasdoors) & (len(alreadythere)==0):
    #             #self.Gr.add_node(j, cycle=k, doors=u)
    #             self.Gr.add_node(rcpt, cycle=[k], transitions=doors)
    #             self.Gr.pos[rcpt] = self.Gt.pos[k]
    #             self.Gr.node[rcpt]['polyg']=self.Gt.node[k]['polyg']
    #             #roomcycles[rcpt].append(k)
    #             cycleroom[k]=rcpt
    #             # add transitions
    #             for ku in u:
    #                 try:
    #                     self.transition[ku].append(rcpt)
    #                 except:
    #                     self.transition[ku] = [rcpt]

    #             # Merge cycles which are separated by an airwall
    #             if len(v) > 0:
    #                 for kv in v:
    #                     ncy  = filter(lambda x : x !=k,self.Gs.node[kv]['ncycles'])[0]
    #                     self.Gr.node[rcpt]['cycle'].append(ncy)
    #             # increment room counter
    #             rcpt += 1
    #         if (len(alreadythere)>0):
    #             ncy = alreadythere[0]
    #             roomn = cycleroom[ncy]
    #             for ncy in notyet:
    #                 print "merging cycle "+str(ncy)+"in room "+str(roomn)
    #                 self.Gr.node[roomn]['polyg'] += self.Gt.node[ncy]['polyg']


    #     # add connection between rooms
    #     for k in self.transition:
    #         room1room2 = self.transition[k]
    #         if len(room1room2) == 2:
    #             self.Gr.add_edge(room1room2[0], room1room2[1])


    def waypoint(self, nroom1, nroom2):
        """
        get the waypoint between room1 and room2
        waypoint = waypoint(nroom1,nroom2)
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
         ax  = fig.add_subplot(111)
         self.display['nodes']=True
         self.display['ednodes']=True
         self.display['subsegnb']=True
         self.display['transition']=True
         self.display['ticksoff']=True

         #self.af = SelectL2(self,fig=fig,ax=ax)
         self.af = SelectL(self,fig=fig,ax=ax)

         fig,ax = self.af.show(fig,ax,clear=True)

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
        #matplotlib.use('TkAgg')

        from matplotlib.backend_bases import key_press_handler
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
        from matplotlib.figure import Figure
        import Tkinter as Tk

        root = Tk.Tk()
        root.wm_title('Pylayers Layout Editor')

        fig = Figure()
        ax  = fig.add_subplot(111)
        #ax.plot(np.arange(10))

        canvas = FigureCanvasTkAgg(fig,master=root)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas,root)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)

        button = Tk.Button(master=root,text='Quit',command=sys.exit)
        button.pack(side=Tk.BOTTOM)

        self.display['nodes']=True
        self.display['ednodes']=True

        select = SelectL(self,canvas)

        #self.af.show(clear=True)

        self.cid1 = canvas.mpl_connect('button_press_event', select.OnClick)
        self.cid2 = canvas.mpl_connect('key_press_event', select.OnPress)
        #ax.axis('tight')
        canvas.show()
        Tk.mainloop()

    def info(self):
        """ gives information about the Layout
        """
        print "filestr : ", self.filename
        print "filematini : ", self.filematini
        print "fileslabini : ", self.fileslabini
        try:
            print "filegeom : ", self.filegeom
        except:
            print "geomfile (.off) has no been generated"

        #self.boundary()
        print "boundaries ",self.ax
        print "number of Points :", self.Np
        print "number of Segments :", self.Ns
        print "number of Sub-Segments :", self.Nss
        try:
            print "Gs Nodes : ", self.Gs.number_of_nodes()
            print "Gs Edges : ", self.Gs.number_of_edges()
        except:
            print "no Gs graph"

        try:
            print "Gt Nodes : ", self.Gt.number_of_nodes()
            print "Gt Edges : ", self.Gt.number_of_edges()
            print "vnodes = Gt.node[Nc]['cycles'].cycle "
            print "poly = Gt.node[Nc]['cycle'].polyg "
        except:
            print "no Gt graph"

        try:
            print "Gr Nodes    :", self.Gr.number_of_nodes()
            print "Gr Edges    :", self.Gr.number_of_edges()
            print "Nc  = Gr.node[nroom]['cycles']  "
        except:
            print "no Gr graph"

    def facets3D(self, edlist, name='Layer', subseg=False):
        """
        facets3d(edlist,name)
        """

        filename = name + '.list'
        filestruc = pyu.getlong(filename, pstruc['DIRGEOM'])
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

    def numseg(self,ta,he):
        """ get segment number from 2 points index

        Parameters
        ----------

        ta  : int <0
        he  : int <0

        Returns
        -------

        nseg : > 0
        if 0 not a segment

        """
        nta = np.array(nx.neighbors(self.Gs,ta))
        nhe = np.array(nx.neighbors(self.Gs,he))
        nseg = np.intersect1d(nta,nhe)
        if len(nseg>0):
            return(nseg[0])
        else:
            return(0)

    def isseg(self,ta,he):
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

        upnt  = filter(lambda x : x<0,self.Gs.nodes())
        ta = np.nonzero(np.array(upnt)==ta)[0][0]
        he = np.nonzero(np.array(upnt)==he)[0][0]
        res = filter(lambda x : (((x[0] == ta) &   (x[1] == he))
                     |           ((x[0] == he) &   (x[1] == ta)))
                     , zip(self.tahe[0],self.tahe[1]))
        if len(res)>0:
            return True
        else:
            return False



    def ispoint(self, pt, tol=0.05):
        """ verify if pt is a point of Layout

        Parameters
        ----------

        pt  : point (2,1)
        tol : float
            default (0.05 meters)

        if True the point number (<0) is returned
        else 0 is return

        Returns
        -------

        pt : point number

        See Also
        --------

        pylayers.util.geomutil.Polygon.setvnodes

        """
        # print "ispoint : pt ", pt
        pts = np.array(self.Gs.pos.values()).T
        ke  = np.array(self.Gs.pos.keys())
        u = pts - pt.reshape(2, 1)
        v = np.sqrt(np.sum(u * u, axis=0))
        nz = (v > tol)
        b = nz.prod()
        if b == 1:
            # if all layout points are different from pt
            return(0)
        else:
            nup = np.where(nz == False)[0]
            if len(nup)==1:
                return(ke[nup][0])
            else :
                mi=np.where(min(v[nup])==v[nup])[0]
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

            #print ta
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
        ----------

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
        filestruc = pyu.getlong(filename, pstruc['DIRGEOM'])
        fos = open(filestruc, "w")
        fos.write("OFF\n")
        fos.write("%d %d \n\n" % (1+(nsseg+1)*4, nsseg+1))
        fos.write("0.000 0.000 0.000\n")
        if subseg:
            try:
                for k,name in enumerate(self.Gs.node[s]['ss_name']):
                    P1[2] = self.Gs.node[s]['ss_z'][k][0]
                    P2[2] = self.Gs.node[s]['ss_z'][k][0]
                    P3[2] = self.Gs.node[s]['ss_z'][k][1]
                    P4[2] = self.Gs.node[s]['ss_z'][k][1]
                    fos.write("%6.3f %6.3f %6.3f \n" % (P1[0], P1[1], P1[2]))
                    fos.write("%6.3f %6.3f %6.3f \n" % (P2[0], P2[1], P2[2]))
                    fos.write("%6.3f %6.3f %6.3f \n" % (P3[0], P3[1], P3[2]))
                    fos.write("%6.3f %6.3f %6.3f \n" % (P4[0], P4[1], P4[2]))
            except:
                print 'no subsegment on ', s
                return('void')
        else:
            name = self.Gs.node[s]['name']
            fos.write("%6.3f %6.3f %6.3f \n" % (P1[0], P1[1], P1[2]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P2[0], P2[1], P2[2]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P3[0], P3[1], P3[2]))
            fos.write("%6.3f %6.3f %6.3f \n" % (P4[0], P4[1], P4[2]))

        if subseg:
            for k,name in enumerate(self.Gs.node[s]['ss_name']):
                colname = sl[name]['color']
                colhex = cold[colname]
                col = pyu.rgb(colhex) / 255.
                fos.write("4 %i %i %i %i %6.3f %6.3f %6.3f 0.4\n" % (1+4*k, 2+4*k,
                3+4*k, 4+4*k, col[0], col[1], col[2]))
        else:
            name = self.Gs.node[s]['name']
            colname = sl[name]['color']
            colhex = cold[colname]
            col = pyu.rgb(colhex) / 255.
            fos.write("4 %i %i %i %i %6.3f %6.3f %6.3f 0.4\n" % (1, 2,
            3, 4, col[0], col[1], col[2]))

        return(filename)

    def geomfile(self,centered=False):
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
            pg = np.sum(self.pt,axis=1)/np.shape(self.pt)[1]
        else:
            pg = np.array([0,0])

        #en  = self.Ns # number of segments
        en  = len(np.where(np.array(self.Gs.node.keys())>0)[0])
        if en != self.Ns:
            logging.warning("wrong number of segment consistency problem in layout")
        #cen = self.Nss
        # d : dictionnary of layout sub segments
        #
        d = self.subseg()
        cen = 0
        for k in d:
            lss = d[k]
            cen = cen + len(lss)

        if cen != self.Nss:
            logging.warning("wrong number of subsegment consistency problem in layout")

        sl = self.sl
#
#        Create a polygon for each segment and subsegment
#
        P1 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P2 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P3 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P4 = np.array(np.zeros([3, en + cen], dtype=np.float64))

        ik   = 0
        dikn = {}
        for i in self.Gs.node.keys():
            if i > 0:  # segment
                if self.Gs.node[i]['name']!='AIR':
                    nebr = self.Gs.neighbors(i)
                    n1 = nebr[0]
                    n2 = nebr[1]
                    P1[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                    P1[2, ik] = self.Gs.node[i]['z'][0]

                    P2[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                    P2[2, ik] = self.Gs.node[i]['z'][0]

                    P3[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                    P3[2, ik] = self.Gs.node[i]['z'][1]

                    P4[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                    P4[2, ik] = self.Gs.node[i]['z'][1]
                    dikn[ik]=i
                    ik = ik + 1
                else:
                    en = en-1
        # d = self.subseg()
        # k : ss_name v: seg number
        cpt = 0
        subseg = {}
        #pdb.set_trace()
        for k in d.keys():
            for l in d[k]:
                ids = l[0]
                subseg[cpt] = ids
                order = l[1]
                cpt = cpt + 1
                nebr = self.Gs.neighbors(l[0])
                n1 = nebr[0]
                n2 = nebr[1]
                #print ik,n1,n2

                P1[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                P1[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                #print P1[:,ik]

                P2[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                P2[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                #print P2[:,ik]

                P3[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                P3[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                #print P3[:,ik]

                P4[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                P4[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                #print P4[:,ik]

                dikn[ik] = l
                ik = ik + 1

        npt = 4 * (en + cen)
        _filename,ext = os.path.splitext(self.filename)
        _filegeom = _filename+'.off'
        self.filegeom=_filegeom
        filegeom = pyu.getlong(_filegeom, pstruc['DIRGEOM'])
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

    def _show3(self,centered=False,newfig=False,opacity=1.,cyid=False,**kwargs):
        """ create a .off geomview file

        Parameters
        ----------

        newfig : Boolean
            create a new mayavi Figure
        opacity : float ([0,1])
            set slab opacity
        centered : Boolean
            if True the layout is centered around its center of gravity
        cyid : boolean
            display cycle number

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

        # calculate center of gravity
        if centered:
            pg = np.sum(self.pt,axis=1)/np.shape(self.pt)[1]
        else:
            pg = np.array([0,0])

        #en  = self.Ns # number of segments
        en  = len(np.where(np.array(self.Gs.node.keys())>0)[0])
        if en != self.Ns:
            logging.warning("wrong number of segment consistency problem in layout")
        #cen = self.Nss
        # d : dictionnary of layout sub segments
        #
        d = self.subseg()
        cen = 0
        for k in d:
            lss = d[k]
            cen = cen + len(lss)

        if cen != self.Nss:
            logging.warning("wrong number of subsegment consistency problem in layout")

        sl = self.sl
#
#        Create a polygon for each segment and subsegment
#
        P1 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P2 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P3 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P4 = np.array(np.zeros([3, en + cen], dtype=np.float64))

        ik   = 0
        dikn = {}

        for i in self.Gs.node.keys():
            if i > 0:  # segment
                if self.Gs.node[i]['name']!='AIR':
                    nebr = self.Gs.neighbors(i)
                    n1 = nebr[0]
                    n2 = nebr[1]
                    P1[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                    P1[2, ik] = self.Gs.node[i]['z'][0]

                    P2[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                    P2[2, ik] = self.Gs.node[i]['z'][1]

                    P3[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                    P3[2, ik] = self.Gs.node[i]['z'][1]

                    P4[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                    P4[2, ik] = self.Gs.node[i]['z'][0]

                    dikn[ik]=i
                    ik = ik + 1

                else:

                    en = en-1


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
                #print ik,n1,n2

                P1[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                P1[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                #print P1[:,ik]

                P2[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                P2[2, ik] = self.Gs.node[ids]['ss_z'][order][0]
                #print P2[:,ik]

                P3[0:2, ik] = np.array(self.Gs.pos[n2])-pg
                P3[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                #print P3[:,ik]

                P4[0:2, ik] = np.array(self.Gs.pos[n1])-pg
                P4[2, ik] = self.Gs.node[ids]['ss_z'][order][1]
                #print P4[:,ik]

                dikn[ik] = l
                ik = ik + 1

        npt = 4 * (en + cen)
        npt_s = (en + cen)

        points = np.hstack((P1[:,0:npt_s],P2[:,0:npt_s]))
        points = np.hstack((points,P3[:,0:npt_s]))
        points = np.hstack((points,P4[:,0:npt_s]))
        points=points.T
        boxes=np.empty((npt/4,4),dtype='int')
        b = np.arange(npt/4)
        boxes[:,0]=b
        boxes[:,1]=b+npt_s
        boxes[:,2]=b+2*npt_s
        boxes[:,3]=b+3*npt_s







#         _filename,ext = os.path.splitext(self.filename)
#         _filegeom = _filename+'.off'
#         self.filegeom=_filegeom
#         filegeom = pyu.getlong(_filegeom, pstruc['DIRGEOM'])
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
        color=np.zeros((4*(cen+en),3))
        for i in range(en + cen):
            q = 4 * i
            if i < en:
                ne = dikn[i]
                name = self.Gs.node[ne]['name']
            else:
                ne = dikn[i][0]
                order = dikn[i][1]
                name = self.Gs.node[ne]['ss_name'][order]

            colname = sl[name]['color']
            colhex = cold[colname]
            color[i,:] = pyu.rgb(colhex)
            color[i+npt_s,:] = pyu.rgb(colhex)
            color[i+2*npt_s,:] = pyu.rgb(colhex)
            color[i+3*npt_s,:] = pyu.rgb(colhex)



        colname = sl['FLOOR']['color']
        colhex = cold[colname]
        colf = np.repeat((pyu.rgb(colhex))[np.newaxis,:],4,axis=0)
        color = np.vstack((color,colf))

        # trick for correcting  color assignement

        sc=tvtk.UnsignedCharArray()
        sc.from_array(color)


        # manage floor

        #if Gt doesn't exists

        try:
            self.ma.coorddeter()
            # z=np.ones(self.ma.xy.shape[1])
            z=np.zeros(self.ma.xy.shape[1])
            F=np.vstack((self.ma.xy,z))
            tri = np.arange(len(z))
            meshf = tvtk.PolyData(points=F.T, polys=np.array([tri]))
            meshf.point_data.scalars = sc
            meshf.point_data.scalars.name = 'scalars'
            surff = mlab.pipeline.surface(meshf, opacity=opacity)
            mlab.pipeline.surface(mlab.pipeline.extract_edges(surff),
                                        color=(0, 0, 0), )

        # otherwise
        except:

            floorx= np.array((points[:,0].min(),points[:,0].max()))
            floory= np.array((points[:,1].min(),points[:,1].max()))
            zmin= np.min(points[:,2])
            Pf = np.array([floorx[0],floory[0],zmin])
            Pf = np.vstack((Pf,np.array([floorx[0],floory[1],zmin])))
            Pf = np.vstack((Pf,np.array([floorx[1],floory[1],zmin])))
            Pf = np.vstack((Pf,np.array([floorx[1],floory[0],zmin])))

            points = np.vstack((points,Pf))
            bf =np.arange(npt,npt+4)
            boxes = np.vstack((boxes,bf))




        mesh = tvtk.PolyData(points=points, polys=boxes)
        mesh.point_data.scalars = sc
        mesh.point_data.scalars.name = 'scalars'


        if newfig:
            mlab.clf()
            f = mlab.figure(bgcolor=(1,1,1))
        else :
            f = mlab.gcf()
            f.scene.background=(1,1,1)

        f.scene.disable_render = True

        surf = mlab.pipeline.surface(mesh, opacity=opacity)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                    color=(0, 0, 0), )
        f.children[-1].name='Layout ' + self.filename

        if cyid:
            if len(self.Gt.nodes())>0:
                pk=self.Gt.pos.keys()
                v=np.array(self.Gt.pos.values())
                [mlab.text3d(v[ik,0],v[ik,1],0.5,str(k)) for ik,k in enumerate(pk)]
        # if segpt:

        #     seg = dict(filter(lambda x: x[0]>0,self.Gs.pos.items()))
        #     pt = dict(filter(lambda x: x[0]<0,self.Gs.pos.items()))
        #     pseg = np.array(seg.values())
        #     ppt = np.array(pt.values())
        #     [mlab.text3d(pseg[ik,0],pseg[ik,1],0.5,str(k)) for ik,k in enumerate(seg)]
        #     [mlab.text3d(ppt[ik,0],ppt[ik,1],3.,str(k)) for ik,k in enumerate(pt)]

        f.scene.disable_render = False

    def show3(self, bdis=True,centered=True):
        """ geomview display of the indoor structure

        Parameters
        ----------

        bdis boolean (default True)
            boolean display (call geowview if True)
        centered : boolean
            if True center the layout before display


        """

        pg = self.geomfile(centered=centered)

        filename = pyu.getlong(self.filegeom, pstruc['DIRGEOM'])
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
            raise NameError('Interaction graph layout.Gi must be build before signature computation')
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
                        #print 'no path between ',nt,nr
                elif (nt != nr):
                    try:
                        path = nx.dijkstra_path(self.Gi, nt, nr)
                        addpath = True
                    except:
                        pass
                        #print 'no path between ',nt,nr
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

    def plot_segments(self,lns,**kwargs):
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
            fig =  kwargs['fig']
            ax = kwargs['ax']

        nth = np.array(map(lambda n: nx.neighbors(self.Gs,n),lns))
        nt = nth[:,0]
        nh = nth[:,1]
        # pt : 2 x Ns
        pt  = np.array(map(lambda n :
                           [self.Gs.pos[n][0],self.Gs.pos[n][1]],nt)).T
        # ph : 2 x Ns
        ph  = np.array(map(lambda n :
                           [self.Gs.pos[n][0],self.Gs.pos[n][1]],nh)).T

        fig,ax = plu.displot(pt,ph,fig=fig,ax=ax,color=kwargs['color'])

        return fig,ax
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
        sig =sigarr[0]
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

    def distwall(self, p, nroom):
        """ calculate distance to wall

        Parameters
        ----------

        p : ndarray
            point coordinate

        nroom : int
            room number of p

        Returns
        -------

        dist
                list of distances to walls of room nroom

        Notes
        -----

        Return  dist a list of all the distances to the walls of a room


        """
        pp = sh.Point(p[0], p[1])

        dist = []
        p0_xy = []
        p1_xy = []

    
        vnode = self.Gr.node[nroom]['cycle'].cycle

        #for j in range(len(Gr[nroom]['vnodes'])):
        for j in range(len(vnodes)):
            nn = self.b_Gr[5]['vnodes'][j]
            nta = G1.tahe[0, nn - 1]
            nhe = G1.tahe[1, nn - 1]
            p0 = np.array([G1.pt[0, nta], G1.pt[1, nta]])
            p1 = np.array([G1.pt[0, nhe], G1.pt[1, nhe]])
            p0_xy.insert(j, p0)
            p1_xy.insert(j, p1)

        pstartwll = np.array(p0_xy)
        pfinwll = np.array(p1_xy)

        for i in range(len(self.b_Gr[nroom]['vnodes'])):
            line_wall = sh.LineString([(pstartwll[i, 0],
                pstartwll[i, 1]), (pfinwll[i, 0], pfinwll[i, 1])])
            dist.insert(i, line_wall.distance(pp))
        return(dist)

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

        #self.boundary()

        Tx_x = rd.uniform(self.ax[0], self.ax[1])
        Tx_y = rd.uniform(self.ax[2], self.ax[3])
        Rx_x = rd.uniform(self.ax[0], self.ax[1])
        Rx_y = rd.uniform(self.ax[2], self.ax[3])

        p_Tx = np.array([Tx_x, Tx_y])
        p_Rx = np.array([Rx_x, Rx_y])

        return(p_Tx, p_Rx)

    def boundary(self, dx=0, dy=0,xlim=()):
        """ add a blank boundary around layout

        Parameters
        ----------

        dx : float
            x offset (default 0)
        dy : float
            y offset (default 0 )

        self.ax is updated

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout('defstr.str')
        >>> L.boundary()

        """
        if not self.hasboundary:
            if len(self.Gs.pos.values())!=0:
                xmax = max(p[0] for p in self.Gs.pos.values())
                xmin = min(p[0] for p in self.Gs.pos.values())
                ymax = max(p[1] for p in self.Gs.pos.values())
                ymin = min(p[1] for p in self.Gs.pos.values())
            elif xlim==():
                xmin = -20.
                xmax = 20.
                ymin = -10.
                ymax = 10.
            else:
                xmin = xlim[0]
                xmax = xlim[1]
                ymin = xlim[2]
                ymax = xlim[3]

            n1 = self.add_fnod((xmin-dx,ymin-dy))
            n2 = self.add_fnod((xmax+dx,ymin-dy))
            n3 = self.add_fnod((xmax+dx,ymax+dy))
            n4 = self.add_fnod((xmin-dx,ymax+dy))

            self.lboundary=[n1,n2,n3,n4]

            self.segboundary = []
            self.segboundary.append(self.add_segment(n1, n2, name='AIR'))
            self.segboundary.append(self.add_segment(n2, n3, name='AIR'))
            self.segboundary.append(self.add_segment(n3, n4, name='AIR'))
            self.segboundary.append(self.add_segment(n4, n1, name='AIR'))

            self.ax = (xmin - dx, xmax + dx, ymin - dy, ymax + dy)
            self.display['box'] = self.ax
            self.hasboundary = True


    def off_overlay(self,dx=0,dy=0):
        """ offset overlay image

        Paramaters
        ----------

        dx : float
        dy : float

        """
        axis = (self.ax[0]+dx,self.ax[1]+dx,self.ax[2]+dy,self.ax[3]+dy)
        self.display['overlay_axis'] = axis

    def scl_overlay(self,ax=1.0,ay=1.0):
        """ scale overlay image

        Paramaters
        ----------

        ax : float
        ay : float

        """
        axis = (self.ax[0]*ax,self.ax[1]*ax,self.ax[2]*ay,self.ax[3]*ay)
        self.display['overlay_axis'] = axis

    def get_paths(self,nd_in, nd_fin):
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


if __name__ == "__main__":
    #plt.ion()
    doctest.testmod()
    #L = Layout('defstr3.ini')
