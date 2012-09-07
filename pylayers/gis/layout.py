# -*- coding:Utf-8 -*-
#
# Class Layout
# Class SelectL
#
#
import pdb
import os
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
import shapely.geometry as sh
from   shapely.ops import cascaded_union
from   descartes.patch import PolygonPatch
import Image

from pylayers.antprop import slab as sb
from pylayers.util import geomutil as geu
from pylayers.util import pyutil as pyu
from pylayers.util import graphutil as gru
# Handle furnitures
import pylayers.gis.furniture as fur
from pylayers.gis import cycles as Cycls
import pylayers.gis.selectl
from pylayers.util.easygui import *
from pylayers.util.project import *
#from   PyUtil  import *
#from   interval import interval
from itertools import combinations
import pdb
#
#


class Layout(object):
    """ Handling of Layout

    Attributes
    ----------
    Gs     : Structure graph
    Gt     : Topological graph  (indicates topological relationships between rooms)
    Gr     : Graph of room
    Gv     : Graph of visibility
    Gc     : Connection graph (indicates visbility relationships)
    Nnode  : Number of nodes of Gs
    Nedge  : Number of edges of Gs
    pt     : points sequence
    tahe   : tail head

    Notes
    ------
     This class exploits `networkx` to store Layout information

    Methods
    -------

        dir_ray_sign(self,Tx_x,Tx_y,Rx_x,Rx_y)
            to be explained
        v_color(self,ob)
            to be explained
        help(self)
            give information about the Layout
        ss-dico = L.subseg(self)
            extract the dictionnary of subseg
            interact(self)  : edit the layout interctively
        add_pnod(self,p,e1,e2)
            tbd
        add_fnod(self,p=(0.0,0.0))
            add free node p
        add_nfpe(self,np0,e1,e2)
            tbd
        add_none(self,ns,alpha=0.5)
           add node on edge
        add_edge(self,n1,n2,name='PARTITION',zmin=0,zmax=3.0)
            add edge between n1 and n2
        del_node(self,n1)
            delete node n1
        del_edge(self,e1)
            add edge between n1 and n2
        displaygui(self)
        info_edge(self,e1)
        edit_edge(self,e1)
        have_subseg(self,e1)
        del_subseg(self,e1)
        add_subseg(self,e1)
        add_window(self,e1)
        add_door(self,e1)
        check(self)
        ispoint(self,pt,tol=0.45)
        onseg(self,pt,tol=0.01)
        checkvis(self,p,edgelist,nodelist)
        visilist(self,p)
        visi_papb(self,pa,pb,edgelist=np.array([]))
        boundary(self)

        clip(self,xmin,xmax,ymin,ymax)
        find_edgelist(self,edgelist,nodelist)
        diag(self,p1,p2,l,al1,al2,quadsel=0)
        nd2ed(self,ndlist)
        ed2nd(self,edlist)
        nl,el = L.zone(self,xmin=0,xmax=10,ymin=0,ymax=10)
        closest_edge(self,p,AAS)
        walls = L.thwall(self,offx,offy)
            extract walls as polygon --> Simpy World
        save(self,filename)
        load(self,filename)
        geomfile(self,filename='Lstruc.off')
        savestr2(self)
            export Layout as an  .str2 file

        show_nodes(self,ndlist=[1e8],size=10,color='b',
            dlabels=False,font_size=15,alpha=1)
        show_edge(self,edlist=[],alpha=1,width=1,size=2,
            color='black',font_size=15,dlabels=False)
        show_edges(self,edlist=[],alpha=1,width=1,color='black',
            dnodes=False,dlabels=False,font_size=15)
        show_layer(self,name,edlist=[],alpha=1,width=1,color='black',
            dnodes=False,dthin=False,dlabels=False,font_size=15)
        showGs(self,ndlist=[],edlist=[])
        show3(self,bdis=True)
        buildGc(self)
        buildGt(self)
        buildGr(self)
        buildGi(self)
        waypoint(self,nroom1,nroom2)
        info(self)
        facets3D(self,edlist,name='Layer',subseg=False)
        facet3D(self,e,subseg=False)
        signature(self,pTx,pRx,NroomTx,NroomRx)
        distwall(self,p,nroom)

    """
    def __init__(self, _fileslab='simul8.slab', _filemat='simul8.mat'):

        mat = sb.MatDB()
        mat.load(_filemat)

        self.sl = sb.SlabDB()
        self.sl.mat = mat
        self.sl.load(_fileslab)

        self.Gs = nx.Graph()
        self.Gc = nx.Graph()
        #self.Gt=nx.Graph()
        self.Gm = nx.Graph()
        self.labels = {}
        self.Gs.pos = {}
        self.Nn = 0
        self.Ne = 0
        self.Nss = 0
        self.filename = 'Lstruc.str2'
        self.display = {}
        self.display['title'] = ''
        self.display['ticksoff'] = True
        self.display['nodes'] = False
        self.display['ndsize'] = 10
        self.display['ndlabel'] = False
        self.display['ndlblsize'] = 10
        self.display['edlblsize'] = 10
        self.display['fontsize'] = 10
        self.display['edlabel'] = False
        self.display['edges'] = True
        self.display['ednodes'] = False
        self.display['subseg'] = True
        self.display['visu'] = False
        self.display['thin'] = False
        self.display['scaled'] = True
        self.display['alpha'] = 0.5
        self.display['layer'] = []
        self.display['clear'] = False
        self.display['activelayer'] = self.sl.keys()
        self.display['overlay'] = False
        #self.display['fileoverlay']="/home/buguen/Pyproject/data/image/"
        self.display['fileoverlay'] = "TA-Office.png"
        self.display['box'] = (-11.4, 19.525, -8.58, 23.41)
        self.display['subLayer'] = self.sl.keys()
        self.name = {}
        for k in self.sl.keys():
            self.name[k] = []

        self.ax = (-10, 10, -10, 10)

    def dir(self, typ='str'):
        """ list the available file in dirstruc

        Parameters
        ----------
        typ : string optional
            {'str'|}

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

        .. plot::
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> from pylayers.gis.layout import *
            >>> L = Layout()
            >>> for _filename in L.dir():
            >>>    L.load(_filename)
            >>>    L.showGs()
            >>>    plt.title(_filename)
            >>>    plt.figure()
            >>>    plt.show()

        """

        pathname = strdir + '/*.' + typ
        lfile_l = glob.glob(pathname)
        lfile_s = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_s.append(fis)
        lfile_s.sort()
        return lfile_s

    def delete(self):
        """ delete Layout graphs

        release  Gs,Gc,Gm

        """
        del self.Gs
        del self.Gc
        del self.Gm
        self.Gs = nx.Graph()
        self.Gc = nx.Graph()
        self.Gm = nx.Graph()

    def check(self):
        """ Check Layout consistency

        See Also
        --------
        GeomUtil.isBetween

        Notes
        -----
        For all edges
            get the 2 vertices
                for all the other vertices
                    check if it belongs to segment
        """
        consistent = True
        for e in self.Gs.node.keys():
            if e > 0:
                n1, n2 = np.array(self.Gs.neighbors(e))
                p1 = np.array(self.Gs.pos[n1])
                p2 = np.array(self.Gs.pos[n2])
                for n in self.Gs.node.keys():
                    if (n < 0) & (n1 != n) & (n2 != n):
                        p = np.array(self.Gs.pos[n])
                        if geu.isBetween(p1, p2, p):
                            print "Warning segment ", e, "contains point ", n
                            return(consistent)

    def clip(self, xmin, xmax, ymin, ymax):
        """ return the list of edges which cross or belong to the clipping zone

         .. todo this is wrong

         Parameters
         ----------
            xmin
            xmax
            ymin
            ymax

         Algorithm :
              1) Determine all segments outside the clipping zone for condition to test
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

        iseg = np.arange(self.Ne)

        return np.setdiff1d(iseg, u)

#    def dir_ray_sign(self,Tx_x,Tx_y,Rx_x,Rx_y):
#        """
#            L.dir_ray_sign(Tx_x,Tx_y,Rx_x,Rx_y): takes transmitter and receiver coordinates
#            and returns the shortest ray and his signature
#        """
#
#        xmin=min(Tx_x,Rx_x)
#        xmax=max(Tx_x,Rx_x)
#        ymin=min(Tx_y,Rx_y)
#        ymax=max(Tx_y,Rx_y)
#
#        seg=self.clip(xmin,xmax,ymin,ymax)
#
#        # creation of direct signature and ray
#        sign=[]
#        sign_Tx_Rx=[]
#        line = sh.LineString([(Tx_x, Tx_y), (Rx_x, Rx_y)])
#        for i in range(len(seg)):
#            n1,n2=self.Gs.neighbors(seg[i]+1)
#            line_i = sh.LineString([(self.Gs.pos[n1][0], self.Gs.pos[n1][1]), (self.Gs.pos[n2][0], self.Gs.pos[n2][1])])
#            if line.crosses(line_i)==True:
#                sign.append(seg[i]+1)
#        fig = plt.figure(1, dpi=90)
#        ax = fig.add_subplot(1,1,1)
#
#        geu.plot_coords2(ax, line)
#        geu.plot_bounds2(ax, line)
#        geu.plot_line2(ax, line)
#        sign_Tx_Rx.append('Tx')
#
#        for i in range(len(sign)):
#                 sign_Tx_Rx.append(sign[i])
#
#        sign_Tx_Rx.append('Rx')
#
#
#
#        self.showGs()
#        plt.draw()
#        plt.show()
#        return(sign_Tx_Rx)

    def help(self):
        print "L.showGs(clear=True)"
        print "L.showGs(edlist=L.subseg()['WOOD'],dthin=False,dlabels=True)"

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
            >>> L = Layout()
            >>> L.load('sircut.str')
            >>> L.loadfur('Furw1.ini')
            >>> f,a = L.showGs()
            >>> plt.show()


        """
        filefur = pyu.getlong(_filefur, 'struc')
        config = ConfigParser.ConfigParser()
        config.read(filefur)
        furname = config.sections()
        self.lfur = []
        for name in furname:
            F = fur.Furniture()
            F.load(_filefur, name)
            self.lfur.append(F)

    def load(self,_filename,_filemat='simul9.mat',_fileslab='simul9.slab'):
        """ load a Layout in different formats

        Parameters
        ----------
        _filename
        _filemat
        _fileslab
        """
        filename,ext=os.path.splitext(_filename)
        if ext=='.str':
            self.loadstr(_filename,_filemat,_fileslab)
        elif ext=='.str2':
            self.loadstr2(_filename,_filemat,_fileslab)
        else:
            raise NameError('layout filename extension not recognized')


    def loadstr(self, _filename, _filemat='simul9.mat', _fileslab='simul9.slab'):
        """ loadstr load a .str de PulsRay

        Parameters
        ----------
        _filename : string
        _filemat  : string
        _fileslab : string

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L=Layout()
        >>> L.loadstr('exemple.str')

        """
        self.filename = _filename
        self.delete()
        mat = sb.MatDB()
        mat.load(_filemat)
        self.sl = sb.SlabDB()
        self.sl.mat = mat
        self.sl.load(_fileslab)

        self.labels = {}
        self.name = {}
        self.Gs.pos = {}

        lname = []
        filename = pyu.getlong(_filename, 'struc')
        fo = open(filename, "rb")
        data = fo.read()
        fo.close()

        #
        # Read : Nn Ns Nss
        #        Number of Nodes           nn
        #        Number of Edges           en
        #        Number of Sub Segments    cen
        #
        data_nn = data[0:4]
        Nn = stru.unpack('i', data_nn)[0]
        data_en = data[4:8]
        Ne = stru.unpack('i', data_en)[0]
        data_cen = data[8:12]
        Nss = stru.unpack('i', data_cen)[0]
        self.Nn = Nn
        self.Ne = Ne
        self.Nss = Nss

        codesl = np.array(np.zeros(Ne), dtype=int)
        codes = np.array(np.zeros(Ne), dtype=int)

        # tahe : tableau des noeuds initiaux et finaux des sgements
        tahe = np.array(np.zeros([2, Ne]), dtype=int)
        ini = 12
        for i in range(Ne):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            #self.tahe[0,i]= stru.unpack('i',dt)[0]-1
            tahe[0, i] = stru.unpack('i', dt)[0] - 1

        ini = stop
        for i in range(Ne):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            #self.tahe[1,i]= stru.unpack('i',dt)[0] -1
            tahe[1, i] = stru.unpack('i', dt)[0] - 1

        # x : tableau des coordonnees x des noeuds
        pt = np.array(np.zeros([2, Nn], dtype=np.float64))
        ini = stop
        for i in range(Nn):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            pt[0, i] = stru.unpack('d', dt)[0]
        # y : tableau des coordinates y des noeuds
        ini = stop
        for i in range(Nn):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            pt[1, i] = stru.unpack('d', dt)[0]
        #--------------------------------------------
        # Node labelling (structure nodes)
        #--------------------------------------------
        for k in range(Nn):
            self.Gs.add_node(-(k + 1))
            self.Gs.pos[-(k + 1)] = (pt[0, k], pt[1, k])
            self.labels[-(k + 1)] = str(-(k + 1))

        #
        # y : type de noeud
        #
        typ = np.array(np.zeros(Nn), dtype=int)
        codep = np.array(np.zeros(Nn), dtype=int)
        ini = stop
        for i in range(Nn):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            typ[i] = stru.unpack('i', dt)[0]
            codep[i] = stru.unpack('i', dt)[0]
        #
        # agi : tableau des angles initiaux des noeuds de type 2
        #
        ag = np.array(np.zeros([3, Nn], dtype=np.float64))
        ini = stop
        for i in range(Nn):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ag[0, i] = stru.unpack('d', dt)[0]
        # agf : tableau des angles finaux des noeuds de type 2
        ini = stop
        for i in range(Nn):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ag[1, i] = stru.unpack('d', dt)[0]
        # nN : tableau des parametres d'ouverture de diedre des noeuds de type 2
        nN = np.array(1.0 * np.zeros(Nn))
        ini = stop
        for i in range(Nn):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ag[2, i] = stru.unpack('d', dt)[0]
        #eml  =
        em = np.array(np.zeros([3, Ne]), dtype=int)
        ini = stop
        for i in range(Ne):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            em[0, i] = stru.unpack('i', dt)[0]
        #emr  =
        ini = stop
        for i in range(Ne):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            em[1, i] = stru.unpack('i', dt)[0]
        #emc  =
        ini = stop
        for i in range(Ne):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            em[2, i] = stru.unpack('i', dt)[0]
            codes[i] = -2
            codesl[i] = em[2, i]
            name = self.sl.di[codesl[i]]
            lname.append(name)
        #thickness =
        thick = np.array(1.0 * np.zeros(Ne))
        ini = stop
        for i in range(Ne):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            thick[i] = stru.unpack('d', dt)[0]
        #ehmin =
        z = np.array(np.zeros([2, Ne], dtype=np.float64))
        ini = stop
        for i in range(Ne):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            z[0, i] = stru.unpack('d', dt)[0]
        #ehmax =
        ini = stop
        for i in range(Ne):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            z[1, i] = stru.unpack('d', dt)[0]

        norm = np.array(np.zeros([2, Ne], dtype=np.float64))
        ini = stop
        for i in range(Ne):
            start = ini + 16 * i
            stop = ini + 16 * (i + 1)
            dt1 = data[start:start + 8]
            norm[0, i] = stru.unpack('d', dt1)[0]
            dt2 = data[start + 8:stop]
            norm[1, i] = stru.unpack('d', dt2)[0]
        #
        # read matrice node-node
        #
        ini = stop
        nd_nd = np.zeros([Nn, Nn], dtype=int)
        for i in range(Nn):
            for j in range(Nn):
                k = Nn * i + j
                start = ini + 4 * k
                stop = ini + 4 * (k + 1)
                dt = data[start:stop]
                nd_nd[i][j] = stru.unpack('i', dt)[0]
        #
        # read matrice node-edge
        #
        ini = stop
        nd_ed = np.zeros([Ne, Nn], dtype=int)
        for i in range(Ne):
            for j in range(Nn):
                k = Nn * i + j
                start = ini + 4 * k
                stop = ini + 4 * (k + 1)
                dt = data[start:stop]
                nd_ed[i][j] = stru.unpack('i', dt)[0]
        #
        # read mat_i
        #
        mat_i = np.array(np.zeros(Ne), dtype=int)
        ini = stop
        for i in range(Ne):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            mat_i[i] = stru.unpack('i', dt)[0]
        #
        # read mat_d
        #
        mat_d = np.array(1.0 * np.zeros(Ne))
        ini = stop
        for i in range(Ne):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            mat_d[i] = stru.unpack('d', dt)[0]
        #
        # read matrice ed-ed
        #
        ini = stop
        ed_ed = np.zeros([Ne, Ne], dtype=int)
        for i in range(Ne):
            for j in range(Ne):
                k = Ne * i + j
                start = ini + 4 * k
                stop = ini + 4 * (k + 1)
                dt = data[start:stop]
                ed_ed[i][j] = stru.unpack('i', dt)[0]

        # Sous segments
        #
        # read ce_core  (A COMPLETER)
        #
        ce_core = np.array(np.zeros(Nss), dtype=int)
        ini = stop
        for i in range(Nss):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            ce_core[i] = stru.unpack('i', dt)[0]
        #
        # read ce_thick
        #
        ce_thick = np.array(np.zeros(Nss))
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_thick[i] = stru.unpack('d', dt)[0]
        #
        # read ce_prop_i
        #
        ce_prop = np.array(np.zeros([2, Nss]), dtype=int)
        ini = stop
        for i in range(Nss):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            ce_prop[0, i] = stru.unpack('i', dt)[0]
        #
        # read ce_wall_floor_ceil
        #
        ce_wall_floor_ceil = np.array(np.zeros(Nss), dtype=int)
        ini = stop
        for i in range(Nss):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            ce_wall_floor_ceil[i] = stru.unpack('i', dt)[0]
        #
        # read ce_ed
        #
        ce_ed = np.array(np.zeros(Nss), dtype=int)
        ini = stop
        for i in range(Nss):
            start = ini + 4 * i
            stop = ini + 4 * (i + 1)
            dt = data[start:stop]
            ce_ed[i] = stru.unpack('i', dt)[0]
        #
        # read ce_prop_d
        #
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_prop[1, i] = stru.unpack('d', dt)[0]
            #   self.ce_prop[i]= stru.unpack('d',dt)[0]
        #
        # read ce_xmin
        #
        ce_xmin = np.array(np.zeros(Nss))
        ce_xmax = np.array(np.zeros(Nss))
        ce_ymin = np.array(np.zeros(Nss))
        ce_ymax = np.array(np.zeros(Nss))
        ce_zmin = np.array(np.zeros(Nss))
        ce_zmax = np.array(np.zeros(Nss))
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_xmin[i] = stru.unpack('d', dt)[0]
        #
        # read ce_xmax
        #
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_xmax[i] = stru.unpack('d', dt)[0]
        #
        # read ce_ymin
        #
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_ymin[i] = stru.unpack('d', dt)[0]
        #
        # read ce_ymax
        #
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_ymax[i] = stru.unpack('d', dt)[0]
        #
        # read ce_zmin
        #
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_zmin[i] = stru.unpack('d', dt)[0]
        #
        # read ce_zmax
        #
        ini = stop
        for i in range(Nss):
            start = ini + 8 * i
            stop = ini + 8 * (i + 1)
            dt = data[start:stop]
            ce_zmax[i] = stru.unpack('d', dt)[0]

        ce = {}
        for i in range(Nss):
            ce[ce_ed[i] - 1] = (ce_core[i],
                                ce_wall_floor_ceil[i],
                                ce_prop[0, i],
                                ce_zmin[i],
                                ce_zmax[i],
                                ce_xmin[i],
                                ce_xmax[i],
                                ce_ymin[i],
                                ce_ymax[i])
        #self.udbox()
        #self.laylist()
        #for i in self.layl:
        #    self.display['Layer'].append(i)
        #    self.display['ActiveLayer'].append(i)

        #----------------------------------------
        # Node labelling (structure edges)
        #----------------------------------------
        for k in range(Ne):
            self.Gs.add_node(k + 1, name=lname[k])
            self.Gs.add_node(k + 1, zmin=z[0, k])
            self.Gs.add_node(k + 1, zmax=z[1, k])
            self.Gs.add_node(k + 1, norm=np.array([norm[0, k],
                                                   norm[1, k], 0.]))
            nta = tahe[0, k]
            nhe = tahe[1, k]
            self.Gs.pos[k + 1] = ((pt[0, nta] + pt[0, nhe]) /
                                  2., (pt[1, nta] + pt[1, nhe]) / 2.)
            self.Gs.add_edge(-(nta + 1), k + 1)
            self.Gs.add_edge(k + 1, -(nhe + 1))
            self.labels[k + 1] = str(k + 1)
            if lname[k] in self.name:
                self.name[lname[k]].append(k + 1)
            else:
                self.name[lname[k]] = [k + 1]
        #
        # Update sub-segment
        #
        for k in ce:
            self.Gs.add_node(k + 1, ss_name=self.sl.di[ce[k][0]])
            self.Gs.add_node(k + 1, ss_ce1=ce[k][1])
            self.Gs.add_node(k + 1, ss_ce2=ce[k][2])
            self.Gs.add_node(k + 1, ss_zmin=ce[k][3])
            self.Gs.add_node(k + 1, ss_zmax=ce[k][4])

        #
        # Create connectivity graph Gc
        #   update Gc with nd_nd ed_ed
        #
        self.Gc = nx.Graph()
        self.Gc.add_nodes_from(self.Gs.nodes())
        pos = self.Gs.pos
        #
        # !! Incomplet  (To Do nd_ed)
        #
        Nn = np.shape(nd_nd)[0]
        for k in range(Nn):
            nnp = -(k + 1)
            kvu = sp.nonzero(nd_nd[k] == 3)
            nc = -kvu[0] - 1
            for l in nc:
                self.Gc.add_edge(nnp, l)

        Ne = np.shape(ed_ed)[0]
        for k in range(Ne):
            ne = k + 1
            kvu = sp.nonzero(ed_ed[k] != 0)
            nc = kvu[0] + 1
            for l in nc:
                self.Gc.add_edge(ne, l)
        self.Gc.pos = pos
        #
        # The numpy format is conserved for acceleration
        #
        self.pt = pt
        self.tahe = tahe
        self.display['activelayer'] = self.name.keys()
        #
        # update boundary
        #
        self.boundary(1, 1)

    def loadstr2(self, _filename, _filemat='simul9.mat', _fileslab='simul9.slab'):
        """ load a Graph from a str2 file

            Parameters
            ----------
            _filename : string
                str2 filename
            _filemat
                mat filename
            _fileslab
                slab filename

            Notes
            -----
            .. sourcecode::

                str2 format is as follow

                Np Ns Nss
                xp_1 yp_1 codep_1
                ...
                xp_Np yp_Np codep_Np
                tail_1 head_1 left_1 core_1 right_1 zmin_1 zmax_1
                ...
                tail_Ns head_Ns left_Ns core_Ns right_Ns zmin_Ns zmax_Ns
                segId_1 SlabCode_1 wall_floor_ceil_1 prop_d_1 zmin_1 zmax_1
                ...
                segId_Nss SlabCode_Nss wall_floor_ceil_Nss prop_d_Nss zmin_Nss zmax_Nss

            Examples
            --------
            >>> from pylayers.gis.layout import *
            >>> L = Layout()
            >>> L.loadstr2('Lstruc.str2')

        """

        self.delete()
        self.filename = _filename
        mat = sb.MatDB()
        mat.load(_filemat)

        self.sl = sb.SlabDB()
        self.sl.mat = mat
        self.sl.load(_fileslab)

        filename = pyu.getlong(_filename, 'struc')
        fo = open(filename)
        lines = fo.readlines()
        fo.close()
        l1 = lines[0].split()
        #
        # Parse the .str2 header NP NSEG NCOSEG
        #
        Nn = int(l1[0])
        Ne = int(l1[1])
        Nss = int(l1[2])
        self.Nn = Nn
        self.Ne = Ne
        self.Nss = Nss

        lname = []

        self.labels = {}
        self.name = {}
        self.Gs.pos = {}

        pt = np.array(np.zeros([2, Nn], dtype=np.float64))
        codep = np.array(np.zeros(Nn, dtype=int))
        ag = np.array(np.zeros([3, Nn], dtype=np.float64))
        tahe = np.array(np.zeros([2, Ne], dtype=int))
        codesl = np.array(np.zeros(Ne), dtype=int)
        codes = np.array(np.zeros(Ne), dtype=int)

        em = np.array(np.zeros([3, Ne]), dtype=int)
        thick = np.array(np.zeros(Ne))

        ed_mat_prop_d = np.array(np.zeros(Ne, dtype=float))
        height = np.array(np.zeros([2, Ne], dtype=float))
        z = np.array(np.zeros([2, Ne], dtype=float))

        ce_ed = np.array(np.zeros(Nss), dtype=int)
        ce_core = np.array(np.zeros(Nss), dtype=int)
        ce_wall_floor_ceil = np.array(np.zeros(Nss))
        ce_prop_d = np.array(np.zeros(Nss, dtype=np.float64))
        ce_zmin = np.array(np.zeros(Nss, dtype=np.float64))
        ce_zmax = np.array(np.zeros(Nss, dtype=np.float64))
        #
        # Read points
        #
        for i in range(Nn):
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
        for k in range(Nn):
            self.Gs.add_node(-(k + 1))
            self.Gs.pos[-(k + 1)] = (pt[0, k], pt[1, k])
            self.labels[-(k + 1)] = str(-(k + 1))
        #
        # Read segments
        #
        for i in range(Ne):
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

        for k in range(Ne):
            #print k, lname[k]
            self.Gs.add_node(k + 1, name=lname[k])
            self.Gs.add_node(k + 1, zmin=z[0, k])
            self.Gs.add_node(k + 1, zmax=z[1, k])
            #self.Gs.add_node(k+1,norm=np.array([norm[0,k],norm[1,k],0.]))
            nta = tahe[0, k] - 1
            nhe = tahe[1, k] - 1
            self.Gs.pos[k + 1] = ((pt[0, nta] + pt[0, nhe]) /
                                  2., (pt[1, nta] + pt[1, nhe]) / 2.)
            self.Gs.add_edge(-(nta + 1), k + 1)
            self.Gs.add_edge(k + 1, -(nhe + 1))
            self.labels[k + 1] = str(k + 1)
            if lname[k] in self.name:
                self.name[lname[k]].append(k + 1)
            else:
                self.name[lname[k]] = [k + 1]
        #
        # Update sub-segment
        #
        for k in ce:
            self.Gs.add_node(k + 1, ss_name=self.sl.di[ce[k][0]])
            self.Gs.add_node(k + 1, ss_ce1=ce[k][1])
            self.Gs.add_node(k + 1, ss_ce2=ce[k][2])
            self.Gs.add_node(k + 1, ss_zmin=ce[k][3])
            self.Gs.add_node(k + 1, ss_zmax=ce[k][4])

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
        self.display['activelayer'] = self.name.keys()
        #self.boundary(1,1)

    def subseg(self):
        """ establish the association : name <->  edgelist

        Returns
        -------
        A dictionnary with sub seg name as key  and edge number as value
        """
        dico = {}
        for k in self.Gs.node.keys():
            dk = self.Gs.node[k]
            if 'ss_name' in dk:
                name = dk['ss_name']
                if name in dico:
                    dico[name].append(k)
                else:
                    dico[name] = [k]
        self.dsseg = dico
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
        p is a tuple

        >>> from pylayers.gis.layout import *
        >>> L = Layout()
        >>> L.loadstr('exemple.str')
        >>> L.add_fnod((10.0,10.0))
        -9


        """
        try:
            num = -(max(-np.array(self.Gs.node.keys())) + 1)
        except:
            num = -1
        self.Gs.add_node(num)
        self.Gc.add_node(num)
        self.Gs.pos[num] = p
        self.Nn = self.Nn + 1
        # update labels
        self.labels[num] = str(num)
        return(num)

    def add_nfpe(self, np0, e1, e2):
        """ Add node on e1 from projection of np0 along e2

        Parameters
        ----------
            np0  : point number
            e1   : edge number 1
            e2   : edge number 2
        """
        np1 = self.Gs.neighbors(e1)
        np2 = self.Gs.neighbors(e2)
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
            self.add_none(e1, 1 - x[0])
        #print x

    def add_none(self, ns, alpha=0.5):
        """
        L.add_none(ns,alpha=0.5) : add node on edge ns
        alpha is a parameterisation on the  edge
        """
        nop = self.Gs.neighbors(ns)
        namens = self.Gs.node[ns]['name']
        zminns = self.Gs.node[ns]['zmin']
        zmaxns = self.Gs.node[ns]['zmax']
        p1 = np.array([self.Gs.pos[nop[0]][0], self.Gs.pos[nop[0]][1]])
        p2 = np.array([self.Gs.pos[nop[1]][0], self.Gs.pos[nop[1]][1]])
        p = tuple(alpha * p1 + (1 - alpha) * p2)
        num = self.add_fnod(p)
        # delete old edge ns
        self.del_edge(ns)
        # add new edge np[0] num
        self.add_edge(nop[0], num, name=namens, zmin=zminns, zmax=zmaxns)
        # add new edge num np[1]
        self.add_edge(num, nop[1], name=namens, zmin=zminns, zmax=zmaxns)

    def add_edge(self, n1, n2, name='PARTITION', zmin=0, zmax=3.0):
        """
        L.add_edge(n1,n2,name='PARTITION',zmin=0,zmax=3.0)     : add edge between n1 and n2
        """
        if ((n1 < 0) & (n2 < 0)):
            nn = np.array(self.Gs.node.keys())
            up = np.nonzero(nn > 0)[0]
            lp = len(up)
            e1 = np.arange(lp) + 1
            e2 = nn[up]
            c = ~np.in1d(e1, e2)
            tn = e1[c]
            #print tn
            try:
                num = tn[0]
            except:
                num = max(self.Gs.node.keys()) + 1
                if num == 0:
                    num = 1
        else:
            print "add_edge : error not a node", n1, n2
            return
        p1 = np.array(self.Gs.pos[n1])
        p2 = np.array(self.Gs.pos[n2])
        p2mp1 = p2 - p1
        t = p2mp1 / np.sqrt(np.dot(p2mp1, p2mp1))
        #
        # n = t x z
        norm = np.array([t[1], -t[0], 0])
        self.Gs.add_node(num, name=name)
        self.Gs.add_node(num, zmin=zmin)
        self.Gs.add_node(num, zmax=zmax)
        self.Gs.add_node(num, norm=norm)
        self.Gs.pos[num] = tuple((p1 + p2) / 2.)
        self.Gs.add_edge(n1, num)
        self.Gs.add_edge(n2, num)
        self.Ne = self.Ne + 1
        # update slab name <-> edge number dictionnary
        try:
            self.name[name].append(num)
        except:
            self.name[name] = [num]
        # update label
        self.labels[num] = str(num)
        if name not in self.display['activelayer']:
            self.display['activelayer'].append(name)
        return(num)

    def del_node(self, ln):
        """ delete node in list ln
        """
        if (type(ln) == np.ndarray):
            ln = list(ln)

        if (type(ln) == int):
            ln = [ln]

        for n1 in ln:
            nbrs = self.Gs.neighbors(n1)
            nbrc = self.Gc.neighbors(n1)
            self.Gs.remove_node(n1)
            self.Gc.remove_node(n1)
            for k in nbrs:
                self.del_edge(k)
            #
            # .. todo :: del_node Layout.py :  Attention Graph Gc non mis a jour
            #
            self.labels.pop(n1)
            self.Nn = self.Nn - 1

    def del_edge(self, le):
        """ delete edge e

        Parameters
        ----------
        le : list of segment number

        Notes
        -----

        """
        if (type(le) == np.ndarray):
            le = list(le)

        if (type(le) == int):
            le = [le]

        for e in le:
            if e > 0:
                name = self.Gs.node[e]['name']
                self.Gs.remove_node(e)
                self.labels.pop(e)
                self.Ne = self.Ne - 1
                # update slab name <-> edge number dictionnary
                self.name[name].remove(e)
                # delete subseg if required

    def del_cycle(self, lnc):
        """ delete a cycle

        Parameters
        ----------

        nc :  cycle number

        """
        if (type(lnc) == np.ndarray):
            lnc = list(lnc)

        if (type(lnc) == int):
            lnc = [lnc]

        for nc in lnc:
            vnodes = np.array(self.Gt.node[nc]['vnodes'])
            neigh = self.Gt.neighbors(nc)
            tvn = np.array([])

            for ncy in neigh:
                vn = np.array(self.Gt.node[ncy]['vnodes'])
                try:
                    tvn = np.hstack((tvn, vn))
                except:
                    tvn = vn

            utvn = np.unique(tvn)
            udel = vnodes[~np.in1d(vnodes, utvn)]

            # delete cycle
            self.Gt.remove_node(nc)
            # delete nodes in udel
            self.del_edge(udel)

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
                if (k != 19) and (k != 54) and (k != 79) and (k != 254):
                    print k
                    cy = self.Gs.node[k]['ncycles']
                    print cy
                    if len(cy) > 2:
                        print "more than 2 cycles in segment ", k

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
                try:
                    self.Gc.remove_node(n)
                except:
                    pass
                try:
                    self.Gv.remove_node(n)
                except:
                    pass

        self.Nn = len(np.nonzero(np.array(self.Gs.node.keys()) < 0)[0])

    def displaygui(self):
        """
        displaygui() : open a GUI for display configuration
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
                                   'fileoverlay',
                                   'box',
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
                                   self.display['fileoverlay'],
                                   str(self.display['box']),
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
            self.display['fileoverlay'] = displaygui[11]
            self.display['box'] = eval(displaygui[12])
            self.display['alpha'] = eval(displaygui[13])

    def info_edge(self, e1):
        """
        info_edge(e1)
        """
        nebd = self.Gs.neighbors(e1)
        n1 = nebd[0]
        n2 = nebd[1]
        nns1 = self.Gs.neighbors(n1)
        nns2 = self.Gs.neighbors(n2)
        de1 = self.Gs.node[e1]
        print n1, ' : ', nns1
        print n2, ' : ', nns2
        print '------------'
        print 'Slab     : ', de1['name']
        print 'zmin (m) : ', de1['zmin']
        print 'zmax (m) : ', de1['zmax']
        try:
            print '------------'
            a = de1['ss_name']
            print 'subseg Slab     : ', de1['ss_name']
            print 'subseg zmin (m) : ', de1['ss_zmin']
            print 'subseg zmax (m) : ', de1['ss_zmax']
        except:
            pass

    def edit_edge(self, e1):
        """
        L.edit_edge(e1)       : edit edge e1
        """
        nebd = self.Gs.neighbors(e1)
        n1 = nebd[0]
        n2 = nebd[1]
        de1 = self.Gs.node[e1]
        title = "Segment (" + str(n1) + ',' + str(n2) + ")"
        message = str(self.sl.keys())
        N = len(de1.keys())
        if N == 4:
            de1k = ['name', 'zmin', 'zmax']
            de1v = [de1['name'], de1['zmin'], de1['zmax']]
        else:
            de1k = ['name', 'zmin', 'zmax', 'ss_name', 'ss_zmin', 'ss_zmax']
            de1v = [de1['name'], de1['zmin'], de1['zmax'], de1[
                'ss_name'], de1['ss_zmin'], de1['ss_zmax']]
        #de1v    = de1.values()
        data = multenterbox(message, title, tuple(de1k), tuple(de1v))
        i = 0
        self.name[de1['name']].remove(e1)
        for k in de1k:
            try:
                self.Gs.node[e1][k] = eval(data[i])
            except:
                self.Gs.node[e1][k] = data[i]
                if k == 'name':
                    try:
                        self.name[data[i]].append(e1)
                    except:
                        self.name[data[i]] = [e1]
            i = i + 1
        #
        # Update L.name
        #

    def have_subseg(self, e1):
        """
        have_subseg
        """
        dk = self.Gs.node[e1]
        if 'ss_name' in dk:
            return True
        else:
            return False

    def del_subseg(self, e1):
        """
        del_subseg(e1)

        del subseg information on e1

        """
        if self.have_subseg(e1):
            self.Gs.node[e1].pop('ss_name')
            self.Gs.node[e1].pop('ss_zmin')
            self.Gs.node[e1].pop('ss_zmax')
            self.Gs.node[e1].pop('ss_ce1')
            self.Gs.node[e1].pop('ss_ce2')
            self.Nss -= 1
        else:
            print "no subseg to delete"

    def add_subseg(self, e1):
        """
        add_subseg(e1)

        add a sub segment on segment e1

        """
        if self.have_subseg(e1):
            print "a subseg already exists"
        else:
            self.info_edge(e1)
            message = str(self.sl.keys())
            data = multenterbox(message, title, ('name',
                                                 'zmin', 'zmax'), ('DOOR', 0, 2.24))

            self.Gs.node[e1]['ss_name'] = data[0]
            self.Gs.node[e1]['ss_zmin'] = eval(data[1])
            self.Gs.node[e1]['ss_zmax'] = eval(data[2])
            self.Gs.node[e1]['ss_ce1'] = 0
            self.Gs.node[e1]['ss_ce2'] = 0
            self.Nss += 1

    def add_window(self, e1):
        pass

    def add_door(self, e1):
        pass

    def find_edgelist(self, edgelist, nodelist):
        """
        edgelist = find_edgelist(edgelist,nodelist)

        edgelist : input edgelist
        nodelist : input nodelist

        return the subset of edgelist

        Not Finished :

        """
        #
        # TODO : eviter utilisation de self.tahe
        #
        tail = self.tahe[0, edgelist]
        head = self.tahe[1, edgelist]

        nt = np.intersect1d_nu[tail, nodelist]
        nh = np.intersect1d_nu[head, nodelist]

        edgelist = edgelist[np.unique(ed_t, ed_h)]
        return(edgelist)

    def diag(self, p1, p2, l, al1, al2, quadsel=0):
        """

        diag (p1,p2,l,al1,al2,quadsel)

        p1  :
        p2  :
        al1
        al2

        quadsel : 0   all quadrant
              2 1
              3 4

        WARNING : Not Tested
        """
        x = self.pt[0, :]
        y = self.pt[1, :]

        #
        # selection du quadran
        #
        if (quadsel == 0):
            u0 = np.arange(self.Nn)
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
        # Parametre de la droite p1 p2 (cas general)
        #
        if ((abs(Dx) > finfo(float).eps) & (abs(Dy) > finfo(float).eps)):
            a = Dy / Dx
            b = p1[1] - a * p1[0]
            b1 = p1[1] + p1[0] / a
            b2 = p2[1] + p2[0] / a

            delta_b = l * L / abs(Dx)
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
            u1 = np.nonzero(x < p1[0] + l / 2.)[0]
            x_u1 = x[u1]
            y_u1 = y[u1]
            u2 = np.nonzero(x_u1 > p1[0] - l / 2.)[0]
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
            u1 = np.nonzero(y < p1[1] + l / 2.)[0]
            y_u1 = y[u1]
            u2 = np.nonzero(y_u1 > p1[1] - l / 2.)[0]
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
        edgelist = np.arange(self.Ne)
        edgelist = self.find_edge_list(edgelist, nodelist)
        return(edgelist)

    def nd2ed(self, ndlist):
        """
        convert nodelist to edgelist
        """
        if type(ndlist) == ndarray:
            ndlist = ndlist.tolist()
            #mecanisme puissant de concatenation de listes
        edlist = []
        for n in ndlist:
            edlist = edlist + self.Gs.adj[n].keys()

        return(unique(edlist))

    def ed2nd(self, edlist):
        """
        convert edgelist to nodelist
        """
        if type(edlist) == ndarray:
            edlist = edlist.tolist()
            # mecanisme de concatenation de listes
        ndlist = []
        for e in edlist:
            ndlist = ndlist + self.Gs.adj[e].keys()

        return(unique(ndlist))

    def zone(self, xmin=0, xmax=10, ymin=0, ymax=10):
        """
        nl,el = zone(xmin=0,xmax=10,ymin=0,,ymax=10):

        """
        ndlist = []
        for n in self.Gs.node.keys():
            if n < 0:
                x = self.Gs.pos[n][0]
                y = self.Gs.pos[n][1]
                if ((x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)):
                    ndlist.append(n)
        edlist = self.nd2ed(ndlist)
        return ndlist, edlist

    def savestr2(self, _filename, furniture=False):
        """ save Layout in .str2 format

        Parameters
        ----------
        _filename : string
            file is  written in the struc directory of the current Project
            directory which is defined through the environment variable $BASENAME
            extension .str2 is added to short name filename
            furniture :  boolean

        Notes
        -----

        To produce the .str file

                > newstruc -str2 file.str2 -conf ../project.conf

        .. todo:: Create a savestr from the Layout Class require Gv

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
        filename = pyu.getlong(_filename, 'struc')
        nn = self.Nn
        ne = self.Ne
        nss = self.Nss

        cnn = str(nn)
        cne = str(ne)
        cnss = str(nss)

        fo = open(filename + '.str2', 'w')
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
            zmin = str(self.Gs.node[i]['zmin'])
            zmax = str(self.Gs.node[i]['zmax'])
            chaine = cta + " " + che + " 1 " + core + " " + " 1 " + \
                " " + " 0 " + " " + zmin + " " + zmax + "\n"
            fo.write(chaine)

        for k, i in enumerate(nseg):
            #
            # sub-segment
            #
            if 'ss_name' in self.Gs.node[i]:
                name = str(self.Gs.node[i]['ss_name'])
                core = str(sl[name]['index'])
                ce1 = str(self.Gs.node[i]['ss_ce1'])
                ce2 = str(self.Gs.node[i]['ss_ce2'])
                zmin = str(self.Gs.node[i]['ss_zmin'])
                zmax = str(self.Gs.node[i]['ss_zmax'])
                chaine = str(k + 1) + " " + core + " " + ce1 + \
                    " " + ce2 + " " + zmin + " " + zmax + "\n"
                fo.write(chaine)

        fo.close()

    def angleonlink(self, p1=np.array([0, 0]), p2=np.array([10, 3])):
        """
            angleonlink(self,p1,p2) return seglist between p1 and p2
            p1 : (1 x 2 )
            p2 : (1 x 2 )
            return (seglist , theta)

            >>> L=Layout()
            >>> L.loadstr('sircut.str','simul9.mat','simul9.slab')
            >>> p1 = np.array([0,0])
            >>> p2 = np.array([10,3])
            >>> seglist,theta = L.angleonlink(p1,p2)

        """
        u = p1 - p2
        nu = np.sqrt(np.dot(u, u))
        un = u / nu

        seglist = self.seginframe(p1, p2)
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
        #print vn
        #print mvn
        #print 'n :',n
        #print 'un : ',unn
        #print 'theta (deg)',the*180./pi

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

    def segpt(self, ptlist=np.array([0])):
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

        >>> L = Layout()
        >>> L.loadstr('exemple.str')
        >>> ptlist  = np.array([0,1])
        >>> seglist = L.segpt(ptlist)

        """
        seglist = np.array([], dtype=int)
        for i in ptlist:
            ut = np.nonzero(self.tahe[0, :] == i)[0]
            uv = np.nonzero(self.tahe[1, :] == i)[0]
            seglist = np.hstack((seglist, ut, uv))
        seglist = np.unique(seglist)

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

            >>> L = Layout()
            >>> L.loadstr('sircut.str')
            >>> p1 = np.array([0,0])
            >>> p2 = np.array([10,10])
            >>> seglist = L.seginframe(p1,p2)
            >>> assert len(seglist)==97,"something has changed in sircut.str"

        """

        max_x = max(p1[0], p2[0])
        min_x = min(p1[0], p2[0])
        max_y = max(p1[1], p2[1])
        min_y = min(p1[1], p2[1])

        Dx = max_x - min_x
        Dy = max_y - min_y

        if (Dy < Dx):
            up = np.nonzero((self.pt[0, :] < max_x) & (self.pt[
                0, :] > min_x))[0]
        else:
            up = np.nonzero((self.pt[1, :] < max_y) & (self.pt[
                1, :] > min_y))[0]

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

    def checkvis(self, p, edgelist, nodelist):
        pass

    def visilist(self, p):
        """ returns the list of nodes from Gc which are visible from point p

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
        """

        This function return the closest segment from p which belong to
        the AAS (Allowed Angular Sector)

        [ns] = closest_edge(self,p,AAS)

        """
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
        w = sp.nonzero(abs(den) < 1e-12)[0]

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

    def save(self, filename):
        """ save Layout

        Parameters
        ----------
        filename : string

        Notes
        -----
            File extension is .gml

        """
        fileGs = filename + 'Gs' + '.gml'
        nx.write_gml(self.Gs, fileGs)
        fileGc = filename + 'Gc' + '.gml'
        nx.write_gml(self.Gc, fileGc)

    def loadG(self, filename):
        """ load Layout

        Parameters
        ----------
        filename : string

        Notes
        -----
            File extension is .gml

        """
        fileGs = filename + 'Gs' + '.gml'
        self.Gs = nx.read_gml(fileGs)
        fileGc = filename + 'Gc' + '.gml'
        self.Gc = nx.read_gml(fileGc)

    def show_nodes(self, ndlist=[1e8], size=10, color='b', dlabels=False, font_size=15, alpha=1):
        """ show nodes

        show_nodes(self,ndlist=[],size=10,color='b'):

        """
        if type(ndlist) == np.ndarray:
            ndlist = ndlist.tolist()
        if len(ndlist) == 0:
            ndlist.append(1e8)
            dlabels = False
        if ndlist[0] == 1e8:
            ndlist = self.Gs.node.keys()
        #elif ndlist[0]==1e8:
        #    ndlist  = self.Gs.node.keys()

        #print ndlist
        nx.draw_networkx_nodes(self.Gs, self.Gs.pos, node_color=color,
                               node_size=size, nodelist=ndlist, alpha=alpha)
        if dlabels:
            dicopos = {}
            dicolab = {}
            for n in ndlist:
                dicopos[n] = np.array(self.Gs.pos[n])
                dicolab[n] = self.labels[n]
            nx.draw_networkx_labels(self.Gs, dicopos, dicolab,
                                    font_size=font_size, font_color=color)

    def show_edge(self, edlist=[], alpha=1, width=1, size=2, color='black', font_size=15, dlabels=False):
        """
        show_edges

        show_edges(self,edlist=[],alpha=1,width=1,color='black')

        """
        if type(ndlist) == ndarray:
            ndlist = ndlist.tolist()

        #print ndlist
        nx.draw_networkx_nodes(
            self.Gs, self.Gs.pos, node_size=size, nodelist=ndlist)
        if dlabels:
            dicopos = {}
            dicolab = {}
            for n in ndlist:
                #dicopos[n]=tuple(np.array(self.Gs.pos[n])+np.array((0.8,0.2)))
                dicopos[n] = np.array(self.Gs.pos[n])
                dicolab[n] = self.labels[n]
            nx.draw_networkx_labels(
                self.Gs, dicopos, dicolab, font_size=font_size)

    def show_edges(self, edlist=[], alpha=1, width=1, color='black', dnodes=False, dlabels=False, font_size=15):
        """
        show_edges

        Parameters
        ----------
            edlist
            alpha
            width
            color
            dnodes
            dlabels
            font_size


        """
        clrlist = []
        cold = pyu.coldict()
        clrlist.append(cold[color])
        ecmap = clr.ListedColormap(clrlist)
        U = self.Gs.edges(edlist)
        ue = (np.ones(2 * len(edlist))).astype('int').tolist()
        nx.draw_networkx_edges(self.Gs, self.Gs.pos, edgelist=U,
                               edge_color=ue, edge_cmap=ecmap, alpha=alpha, width=width)
        if dlabels:
               # print edlist
               # nodelist = self.ed2nd(edlist)
            self.show_nodes(ndlist=edlist, dlabels=dlabels,
                            color='b', font_size=font_size)
        if dnodes:
            self.show_nodes(ndlist=edlist, color='b')

    def show_layer(self, name, edlist=[], alpha=1, width=1, color='black', dnodes=False, dthin=False, dlabels=False, font_size=15):
        """
        show_layer

        show_layer(name,alpha=1,width=1,color='black',dthin=False,dlabels=False):


        """
        if edlist == []:
            edlist = self.name[name]
        else:
            # intersect layer edge list with local zone edge list (in function argument)
            a1 = np.array(self.name[name])
            a2 = np.array(edlist)
            edlist = list(np.intersect1d(a1, a2))

        if self.display['thin']:
            self.show_edges(edlist, alpha=1, width=1,
                            color=color, dlabels=dlabels, font_size=font_size)
        else:
            slab = self.sl[name]
            linewidth = slab['linewidth'] / 3.
            color = slab['color']
            self.show_edges(edlist, alpha=1, width=linewidth, color=color, dnodes=dnodes, dlabels=dlabels, font_size=font_size)

    def showGt(self, ax=[], roomlist=[]):
        """
        """
        if not isinstance(ax, plt.Axes):
            fig = plt.gcf()
            ax = fig.gca()

        for k, nc in enumerate(self.Gt.node.keys()):
            poly = self.Gt.node[nc]['polyg']
            a = poly.signedarea()
            if poly.vnodes[0] < 0:
                if a < 0:
                    poly.plot(color='red')
                else:
                    poly.plot(color='red', alpha=0.5)
            else:
                if a < 0:
                    poly.plot(color='blue')
                else:
                    poly.plot(color='blue', alpha=0.5)
        ax.axis('scaled')

    def showGs(self, ax=[], ndlist=[], edlist=[], show=False, furniture=False, roomlist=[]):
        """
        show structure graph Gs

        Parameters
        ----------
        ndlist  : np.array
            set of nodes to be displayed
        edlist  : np.array
            set of edges to be displayed
        show    : boolean
            default True
        furniture : boolean
            default False

        display parameters are defined in  L.display dictionnary

        """

        if not isinstance(ax, plt.Axes):
            fig = plt.gcf()
            ax = fig.add_subplot(111)

        if furniture:
            if 'lfur' in self.__dict__:
                for fur1 in self.lfur:
                    if fur1.Matname == 'METAL':
                        fur1.show(fig, ax)
            else:
                print "Warning : no furniture file loaded"

        if self.display['clear']:
            ax.cla()
        if self.display['overlay']:
            image = Image.open(strdir + '/' + self.display['fileoverlay'])
            ax.imshow(image, origin='lower', extent=(0, 40, 0, 15), alpha=0.8)
        if ndlist == []:
            tn = np.array(self.Gs.node.keys())
            u = np.nonzero(tn < 0)[0]
            ndlist = tn[u]
        if edlist == []:
            tn = np.array(self.Gs.node.keys())
            u = np.nonzero(tn > 0)[0]
            edlist = tn[u]
        if self.display['nodes']:
            dlabels = self.display['ndlabel']
            self.show_nodes(ndlist, size=10, color='r', dlabels=dlabels)
        slablist = self.name.keys()
        if self.display['edges']:
            dlabels = self.display['edlabel']
            font_size = self.display['fontsize']
            dnodes = self.display['ednodes']
            dthin = self.display['thin']
            alpha = self.display['alpha']
            for nameslab in self.display['activelayer']:
                #print len(edlist)
                self.show_layer(nameslab, edlist=edlist, alpha=alpha, dthin=dthin, dnodes=dnodes, dlabels=dlabels, font_size=font_size)
        if self.display['subseg']:
            dico = self.subseg()
            for k in dico.keys():
                edlist = dico[k]
                color = self.sl[k]['color']
                self.show_edges(edlist=edlist, color='black', alpha=1)

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

        for nr in roomlist:
            ncy = self.Gr.node[nr]['cycle']
            self.Gt.node[ncy]['polyg'].plot()

        if show:
            plt.show()

        fig = plt.gcf()
        return fig, ax

    def buildGc(self):
        """ build the connectivity graph

        nd_nd  : node to node only convex to convex visibility is taken into account
        nd_ed  : node to edge
        ed_ed  : edge to edge

        .. todo: To be Continued
        Faire ce travail piece par piece
        Ce code implemnte une condition necessaire mais non suffisante

        Il existe des points de degre <=2 qui ne sont pas diffractant
        si en zone concave

                __________
                |
                |
                |
                |
        Return
        ------
        ncoin , ndiff

        """
        #
        # First step
        #
        #for nr in self.Gr.node()
        ncoin = np.array([])
        ndiff = np.array([])
        # first step
        # find all node (<0) with  degree < 3
        #
        for n in self.Gs.nodes():
            deg = self.Gs.degree(n)
            if deg > 2:
                ncoin = np.hstack((ncoin, n)).astype('int')
            else:
                if n < 0:
                    ndiff = np.hstack((ndiff, n)).astype('int')
        return ncoin, ndiff

    def buildGt(self):
        """ Built topological graph Gt

        Notes
        -----
        1. Exploit `cycle_basis` function of NetworkX
        2. Each discovered cycle in the graph Gs is transform in a Cycles.Cycles object
        3. LC (List of Cycles) contains the list of all these Cycle object
        4. Create graph G_t each cycle of Gs is a node of Gt
        5. Seek for Cycle inter connectivity
                Algorithm :
                    For k in cycles :
                        vnodesk = get vnodes(k)
                            For l in cycles > k
                                vnodesl = get vnodes(l)
                                    nkinnl = vnodesk :math:\cap vnodesl



        See Also
        --------
            nx.algorithms.cycles.cycle_basis

        """
        C = nx.algorithms.cycles.cycle_basis(self.Gs)
        LC = []
        for c in C:
            Cy = Cycls.Cycle(self.Gs, c)
            LC.append(Cy)

        Cys = Cycls.Cycles(LC, self.Gs)
        self.Gt = Cys.Gt

        #N = len(self.Gt.nodes())
        #for k in self.Gt.nodes():
        #    nk = np.array(self.Gt.node[k]['vnodes'])
        #    for l in np.arange(k+1,N):
        #        nl = np.array(self.Gt.node[l]['vnodes'])
        #        nkinl = np.intersect1d(nk,nl)
        #        if len(nkinl!=0):
        #            self.Gt.add_edge(k,l)
        Ncycles = len(self.Gt.nodes())

        #
        #  Update graph Gs with cycle number information
        #
        for k in range(Ncycles):
            vnodes = np.array(self.Gt.node[k]['vnodes'])
            for n in vnodes:
                try:
                    self.Gs.node[n]['ncycles'].append(k)
                except:
                    self.Gs.node[n]['ncycles'] = [k]

        #
        #  Seek for Cycle inter connectivity
        #
        for k in combinations(range(Ncycles), 2):
            vnodes0 = np.array(self.Gt.node[k[0]]['vnodes'])
            vnodes1 = np.array(self.Gt.node[k[1]]['vnodes'])
            #
            # Connect Cycles if they share nodes
            #
            intersection_vnodes = np.intersect1d(vnodes0, vnodes1)

            if len(intersection_vnodes != 0):
                self.Gt.add_edge(k[0], k[1])

        #
        # Construct the polygon associated to each cycle
        #
        for k in self.Gt.nodes():
            vnodes = self.Gt.node[k]['vnodes']
            u_neg = np.nonzero(vnodes < 0)[0]
            npoints = vnodes[u_neg]
            coords = []
            #
            # Loop over points
            #
            for ind in npoints:
                coords.append(self.Gs.pos[ind])
            polk = geu.Polygon(sh.MultiPoint(tuple(coords)), vnodes)

            self.Gt.add_node(k, polyg=polk)
        #
        # Construct the list of interactions associated to each cycle
        #
        # Interaction labeling convention
        #
        #   negative integer : Diffraction on point |ni|
        #   positive integer : Transmission through segment ni
        #   tuple (nseg,ncycle) : Reflection on nseg toward cycle ncycle
        #
        #
        #    At that stage the diffraction point are not included
        #    not enough information available
        #
        for k in self.Gt.nodes():
            vnodes = self.Gt.node[k]['vnodes']
            ListInteractions = []
            for inode in vnodes:
                if inode > 0:
                    ListInteractions.append(str(inode))
                    ListInteractions.append(str((inode, k)))
            self.Gt.add_node(k, inter=ListInteractions)

    def buildGw(self):
        """ build Graph of waypaths

        See Also
        --------
        buildGr

        """
        self.Gw = nx.Graph()
        self.Gw.pos = {}
        d_id = max(self.Gr.nodes())
        for e in self.Gr.edges_iter():
            doors1 = self.Gr.node[e[0]]['doors']
            doors2 = self.Gr.node[e[1]]['doors']
#            if len(doors1) > 1:
#                pdb.set_trace()
            doorId = np.intersect1d(doors1, doors2)[0]
            unode = self.Gs.neighbors(doorId)
            p1 = self.Gs.pos[unode[0]]
            p2 = self.Gs.pos[unode[1]]
            pdoor = (np.array(p1) + np.array(p2)) / 2
            self.Gw.add_node(doorId + d_id)
            self.Gw.pos[doorId + d_id] = pdoor
            self.Gw.add_edges_from(
                [(e[0], doorId + d_id), (e[1], doorId + d_id)])
            self.Gw.pos.update(self.Gr.pos)
        for n in self.Gr.nodes_iter():
            d = self.Gw.neighbors(n)
            if len(d) > 1:
                self.Gw.add_edges_from(combinations(d, 2))

    def buildGv(self, show=False):
        """ build global visibility graph

        Parameters
        ----------
        display : boolean
            default False

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout()
        >>> L.loadstr('exemple.str')
        >>> L.buildGt()
        >>> L.buildGr()
        >>> L.buildGv()

        """

        self.Gv = nx.Graph()
        #
        # loop over rooms
        #

        for nr in self.Gr.node:
            udeg2 = []
            udeg1 = []
            icycle = self.Gr.node[nr]['cycle']  # id of cycle
            room = self.Gt.node[icycle]      # cycle of the room
            polyg = room['polyg']             # pol
            vnodes = room['vnodes']
            #
            # seek node of degree 2
            #
            # udeg2 is the index of the deg 2 point in the sequence of points
            for ik, inode in enumerate(vnodes):
                deg = self.Gs.degree(inode)
                if vnodes[0] < 0:
                    index = ik / 2
                else:
                    index = (ik - 1) / 2
                if inode < 0:
                    if deg == 2:
                        udeg2.append(index)
                    if deg == 1:
                        udeg1.append(index)    # warning not used
            Gv = polyg.buildGv(show=show, udeg2=udeg2)
            #
            # Graph Gv aggregation
            #
            self.Gv = nx.compose(self.Gv, Gv)

    def buildGi(self):
        """ Build graph of interaction

        Notes
        -----

        For each node > of graph Gs creates
            3 different nodes associated to the same segment
            R+ T R-

        """
        self.Gi = nx.Graph()
        self.Gi.pos = {}
        for n in self.Gv.node:
            if n < 0:
                self.Gi.add_node(str(n))
                self.Gi.pos[n] = self.Gs.pos[n]
            if n > 0:
                cy = self.Gs.node[n]['ncycles']
                if n == 129:
                    pdb.set_trace()
                if len(cy) == 2:
                    cy0 = cy[0]
                    cy1 = cy[1]
                    self.Gi.add_node(str((n, cy0)))
                    self.Gi.add_node(str((n, cy1)))
                    self.Gi.add_node(str(n))
                    nei = self.Gs.neighbors(n)
                    np1 = nei[0]
                    np2 = nei[1]
                    p1 = np.array(self.Gs.pos[np1])
                    p2 = np.array(self.Gs.pos[np2])
                    l = p1 - p2
                    nl = np.dot(l, l)
                    ln = l / nl
                    delta = nl / 10
                    self.Gi.pos[n] = self.Gs.pos[n]
                    self.Gi.pos[(n, cy0)] = self.Gs.pos[n] + ln * delta
                    self.Gi.pos[(n, cy1)] = self.Gs.pos[n] - ln * delta

                if len(cy) == 1:
                    self.Gi.add_node(str((n, cy[0])))
                    self.Gi.pos[(n, cy[0])] = self.Gs.pos[n]

        #
        # Loop over interactions
        #
        for sn in self.Gi.node:
            n = eval(sn)
            print n
            if n == (161, 53):
                pass
                #pdb.set_trace()
            if isinstance(n, tuple):  # reflection
                ns = n[0]
                nc = n[1]
                vnodes = self.Gt.node[nc]['vnodes']
                neigh = self.Gv.neighbors(ns)  # find neighbors
                for nb in neigh:
                    if nb in vnodes:           # Si Voisin dans cycle reflexion
                        if nb > 0:
                            node1 = str(n)
                            node2 = str((nb, nc))
                            if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                self.Gi.add_edge(node1, node2)
                                    # R-R
                            else:
                                print node1, node2
                                pdb.set_trace()

                            if len(self.Gs.node[nb]['ncycles']) == 2:  # R-T
                                node1 = str(n)
                                node2 = str(nb)
                                if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                    self.Gi.add_edge(node1, node2)
                                else:
                                    print node1, node2
                                    pdb_set_trace()
                        else:                                    # R-D
                            node1 = str(n)
                            node2 = str(nb)
                            if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                self.Gi.add_edge(node1, node2)
                            else:
                                print node1, node2
                                pdb_set_trace()
            else:  # transmission or diffraction
                if n > 0:  # transmission
                    ns = n
                    cy = self.Gs.node[n]['ncycles']
                    cy0 = cy[0]
                    cy1 = cy[1]
                    vnodes0 = self.Gt.node[cy0]['vnodes']
                    vnodes1 = self.Gt.node[cy1]['vnodes']
                    neigh = self.Gv.neighbors(ns)  # find neighbors
                    for nb in neigh:
                        if nb in vnodes0:    # Si Voisin dans cycle reflexion 0
                            if nb > 0:
                                node1 = str(n)
                                node2 = str((nb, cy0))
                                if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                    self.Gi.add_edge(node1, node2)
                                else:
                                    pdb.set_trace()
                                if len(self.Gs.node[nb]['ncycles']) == 2:  # R-T
                                    node1 = str(n)
                                    node2 = str(nb)
                                    if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                        self.Gi.add_edge(node1, node2)
                                    else:
                                        print node1, node2
                                        pdb.set_trace()
                            else:
                                node1 = str(n)
                                node2 = str(nb)
                                if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                    self.Gi.add_edge(str(n), str(nb))
                                else:
                                    print node1, node2
                                    pdb.set_trace()
                        if nb in vnodes1:    # Si Voisin dans cycle reflexion 1
                            if nb > 0:
                                node1 = str(n)
                                node2 = str((nb, cy1))
                                if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                    self.Gi.add_edge(node1, node2)
                                else:
                                    print node1, node2
                                    pdb.set_trace()
                                if len(self.Gs.node[nb]['ncycles']) == 2:  # R-T
                                    node1 = str(n)
                                    node2 = str(nb)
                                    if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                        self.Gi.add_edge(node1, node2)
                                    else:
                                        print node1, node2
                                        pdb.set_trace()
                            else:
                                node1 = str(n)
                                node2 = str(nb)
                                if (node1 in self.Gi.node.keys()) & (node2 in self.Gi.node.keys()):
                                    self.Gi.add_edge(node1, node2)
                                else:
                                    print node1, node2
                                    pdb.set_trace()

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
    def showG(self, graph='r', **kwargs):
        """
        Parameters
        ----------
        show : boolean
        fig
        nodes
        eded
        ndnd
        nded
        linewidth
        nodelist
        """

        defaults = {'show': False,
                    'fig': [],
                    'ax': [],
                    'nodes': False,
                    'eded': True,
                    'ndnd': True,
                    'nded': True,
                    'linewidth': 2,
                    'nodelist': []
                    }

        for key, value in defaults.items():
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, value)
                kwargs[key] = value

        if kwargs['fig'] == []:
            fig = plt.figure()
            fig.set_frameon(True)
        else:
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            ax = fig.gca()
        else:
            ax = kwargs['ax']

        if 't' in graph:
            G = self.Gt
            nx.draw(G, G.pos, node_color='r', edge_color='r')
        if 'r' in graph:
            G = self.Gr
            nx.draw(G, G.pos, node_color='g', edge_color='g')
        if 's' in graph:
            G = self.Gs
            nx.draw(G, G.pos, node_color='b', edge_color='b')
        if 'v' in graph:
            G = self.Gv
            nx.draw(G, self.Gs.pos, node_color='m', edge_color='m')

        for k, ncy in enumerate(self.Gt.node.keys()):
            self.Gt.node[ncy]['polyg'].plot()

        ax.axis('scaled')
        # Display doors and windows
        d = self.subseg()
        for ss in d.keys():
            if ss == 'DOOR':
                color = 'red'
            if ss == '3D_WINDOW_GLASS':
                color = 'blue'
            if ss == 'WINDOW_GLASS':
                color = 'cyan'
            for ns in d[ss]:
                np1, np2 = self.Gs.neighbors(ns)
                x = [self.Gs.pos[np1][0], self.Gs.pos[np2][0]]
                y = [self.Gs.pos[np1][1], self.Gs.pos[np2][1]]
                ax.plot(x, y, linewidth=2, color=color)

        if kwargs['show']:
            plt.show()

    def showGv(self, **kwargs):
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

        Examples
        --------

        .. plot::
           :include-source:

            >>> from pylayers.gis.Layout import *
            >>> L = Layout()
            >>> L.loadstr('exemple.str')
            >>> L.buildGt()
            >>> L.buildGr()
            >>> L.buildGv()
            >>> fig,ax=L.showGs()
            >>> L.showGv(fig=fig,ax=ax)
            >>> plt.axis('off')
            >>> plt.show()

        """
        defaults = {'show': False,
                    'fig': [],
                    'ax': [],
                    'nodes': False,
                    'eded': True,
                    'ndnd': True,
                    'nded': True,
                    'linewidth': 2
                    }

        for key, value in defaults.items():
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, value)
                kwargs[key] = value

        if kwargs['fig'] == []:
            fig = plt.figure()
            fig.set_frameon(True)
        else:
            fig = kwargs['fig']

        if kwargs['ax'] == []:
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
                                   edgelist=eded, edge_color='blue', linewidth=2)
        if kwargs['ndnd']:
            nx.draw_networkx_edges(self.Gv, self.Gs.pos,
                                   edgelist=ndnd, edge_color='red', linewidth=2)
        if kwargs['nded']:
            nx.draw_networkx_edges(self.Gv, self.Gs.pos,
                                   edgelist=nded, edge_color='green', linewidth=2)

        if kwargs['show']:
            plt.show()

        return(fig, ax)

    def waypointGw(self, nroom1, nroom2):
        """ get the waypoint between room1 and room2

        Parameters
        ----------
            nroom1
            nroom2

        Examples
        --------
            >>> from pylayers.gis.layout import *
            >>> L = Layout()
            >>> L.loadstr('Lstruc.str')
            >>> L.buildGt()
            >>> L.buildGr()
            >>> L.buildGw()
            >>> nroom1 = 1
            >>> nroom2 = 6
            >>> waypoint = L.waypointGw(nroom1,nroom2)

        """
        rooms = nx.dijkstra_path(self.Gw, nroom1, nroom2)
        return([tuple(self.Gw.pos[i]) for i in rooms])

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
        >>> L = Layout()
        >>> L.loadstr('sircut.str2')
        >>> L.thwall(0,0)

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
                try:
                    ss_name = self.Gs.node[nd]['ss_name']
                except:
                    ss_name = ''
                l = self.sl[name]
                thick = sum(l['lthick'])

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
                #if ss_name!="WOOD":
                if ss_name != "DOOR":
                    walls.append(wall)
        return(walls)

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
        for nr in self.Gr.node.keys():
            if self.Gt.node[self.Gr.node[nr]['cycle']]['polyg'].contains(ptsh):
                return(nr)

    def buildGr(self):
        """ build Graph of room

        Summary
        -------
            A room is a cycle with at least one door.
            This function requires graph Gt
        """
        self.Gr = nx.Graph()
        self.Gr.pos = {}
        self.Door = {}
        d = self.subseg()
        #ldoorseg    = np.array(d['WOOD'])
        ldoorseg = np.array(d['DOOR'])
        j = 0
        #
        # Pour tous les cycles
        #
        for k in self.Gt.node:
            lseg = self.Gt.node[k]['vnodes']
            u = np.intersect1d(lseg, ldoorseg)
            #
            # Si le cycle contient une porte, c'est une piece
            #
            if len(u) > 0:
                self.Gr.add_node(j, cycle=k, doors=u)
                self.Gr.pos[j] = self.Gt.pos[k]
                for ku in u:
                    try:
                        self.Door[ku].append(j)
                    except:
                        self.Door[ku] = [j]
                j = j + 1
        for k in self.Door:
            room1room2 = self.Door[k]
            #print room1room2
            # Une porte d'entree separe une piece de l'exterieur du
            # batiment
            # On pourra ulterieurement creer un noeud exterieur pour
            # l'inter indoor
            if len(room1room2) == 2:
                self.Gr.add_edge(room1room2[0], room1room2[1])

    def waypoint(self, nroom1, nroom2):
        """
        get the waypoint between room1 and room2
        waypoint = L.waypoint(nroom1,nroom2)
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

    def interact(self):
        """ interact

        Notes
        -----
        point edition

            p create point

                lclic same x
                rclic same y
                cclic free point

        segment edition

            [0-f] - display one of the 16 first layers
            x : save structure
            o : toggle overlay

        """
        fig = plt.figure()
        plt.axis(self.ax)
        self.af = SelectL.SelectL(self, fig)
        self.af.show()
        self.cid1 = fig.canvas.mpl_connect(
            'button_press_event', self.af.OnClick)
        self.cid2 = fig.canvas.mpl_connect('key_press_event', self.af.OnPress)
        #fig.canvas.mpl_connect('button_press_event',af.OnClick)
        #fig.canvas.mpl_connect('key_press_event',af.OnPress)
        plt.draw()
        plt.show()

    def info(self):
        """ gives information about the Layout
        """
        print "filestruc : ", self.filename
        try:
            print "filegeom : ", self.filegeom
        except:
            print "geomfile (.off) has no been generated"

        print "number of Nodes :", self.Nn
        print "number of Segments :", self.Ne
        print "number of Sub-Segments :", self.Nss
        try:
            print "Gs Nodes : ", self.Gs.number_of_nodes()
            print "Gs Edges : ", self.Gs.number_of_edges()
        except:
            print "no Gs graph"

        try:
            print "Gc Nodes : ", self.Gc.number_of_nodes()
            print "Gc Edges : ", self.Gc.number_of_edges()
        except:
            print "no Gc graph"

        try:
            print "Gt Nodes : ", self.Gt.number_of_nodes()
            print "Gt Edges : ", self.Gt.number_of_edges()
            print "vnodes = L.Gt.node[Nc]['vnodes'] "
            print "poly = L.Gt.node[Nc]['polyg'] "
        except:
            print "no Gt graph"

        try:
            print "Gr Nodes    :", self.Gr.number_of_nodes()
            print "Gr Edges    :", self.Gr.number_of_edges()
            print "Nc  = L.Gr.node[nroom]['cycles']  "
        except:
            print "no Gr graph"

    def facets3D(self, edlist, name='Layer', subseg=False):
        """
        facets3d(edlist,name)
        """

        filename = name + '.list'
        filestruc = pyu.getlong(filename, "geom")
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

    def ispoint(self, pt, tol=0.45):
        """
        ispoint(pt,tol) : verify if pt is a point of Layout

        if True the point number (numbered from 1) is return
        else -1 is return

        """
        print "ispoint : pt ", pt
        pts = np.array(self.Gs.pos.values()).T
        ke = np.array(self.Gs.pos.keys())
        u = pts - pt.reshape(2, 1)
        v = np.sqrt(np.sum(u * u, axis=0))
        nz = (v > tol)
        b = nz.prod()
        if b == 1:
            return(0)
        else:
            nup = np.nonzero(nz == False)
            return(ke[nup[0]][0])

    def onseg(self, pt, tol=0.01):
        """
        onseg(pt,tol)

        return the segment number which contains point pt

        pt  np.array(1x2)  it is a 2D point
        tol = 0.01      tolerance

        """

        pts = np.array(self.Gs.pos.values()).T
        ke = np.array(self.Gs.pos.keys())
        n = np.shape(pts)[1]
        nbu = np.array([])
        if (n > 0):
            num = np.arange(n)
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
        """ facet3D

        Parameters
        ----------
        e
        subseg : boolean
            default False
        """
        P1 = np.array(np.zeros(3), dtype=np.float64)
        P2 = np.array(np.zeros(3), dtype=np.float64)
        P3 = np.array(np.zeros(3), dtype=np.float64)
        P4 = np.array(np.zeros(3), dtype=np.float64)
        nebr = self.Gs.neighbors(e)
        n1 = nebr[0]
        n2 = nebr[1]
        P1[0:2] = np.array(self.Gs.pos[n1])
        P1[2] = self.Gs.node[e]['zmin']

        P2[0:2] = np.array(self.Gs.pos[n2])
        P2[2] = self.Gs.node[e]['zmin']

        P3[0:2] = np.array(self.Gs.pos[n2])
        P3[2] = self.Gs.node[e]['zmax']

        P4[0:2] = np.array(self.Gs.pos[n1])
        P4[2] = self.Gs.node[e]['zmax']

        cold = pyu.coldict()
        if subseg:
            try:
                name = self.Gs.node[e]['ss_name']
                P1[2] = self.Gs.node[e]['ss_zmin']
                P2[2] = self.Gs.node[e]['ss_zmin']
                P3[2] = self.Gs.node[e]['ss_zmax']
                P4[2] = self.Gs.node[e]['ss_zmax']
            except:
                print 'no subsegment on ', e
                return('void')
        else:
            name = self.Gs.node[e]['name']
        colname = sl[name]['color']
        colhex = cold[colname]
        col = pyu.rgb(colhex) / 255.

        filename = 'fa' + str(e) + '.off'
        filestruc = pyu.getlong(filename, "geom")
        fos = open(filestruc, "w")
        fos.write("OFF\n")
        fos.write("%d %d \n\n" % (5, 1))
        fos.write("0.000 0.000 0.000\n")
        fos.write("%6.3f %6.3f %6.3f \n" % (P1[0], P1[1], P1[2]))
        fos.write("%6.3f %6.3f %6.3f \n" % (P2[0], P2[1], P2[2]))
        fos.write("%6.3f %6.3f %6.3f \n" % (P3[0], P3[1], P3[2]))
        fos.write("%6.3f %6.3f %6.3f \n" % (P4[0], P4[1], P4[2]))
        fos.write("4 %i %i %i %i %6.3f %6.3f %6.3f 0.4\n" % (1, 2,
            3, 4, col[0], col[1], col[2]))
        return(filename)

    def geomfile(self):
        """ create a geomview file of the layout

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> L = Layout()
        >>> L.geomfile()

        """
        en = self.Ne
        cen = self.Nss

        sl = self.sl
#
#        Creation d'un plan pour chaque segment et chaque sous-segment
#
        P1 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P2 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P3 = np.array(np.zeros([3, en + cen], dtype=np.float64))
        P4 = np.array(np.zeros([3, en + cen], dtype=np.float64))

        ik = 0
        for i in self.Gs.node.keys():
            if i > 0:
                nebr = self.Gs.neighbors(i)
                n1 = nebr[0]
                n2 = nebr[1]
                P1[0:2, ik] = np.array(self.Gs.pos[n1])
                P1[2, ik] = self.Gs.node[i]['zmin']

                P2[0:2, ik] = np.array(self.Gs.pos[n2])
                P2[2, ik] = self.Gs.node[i]['zmin']

                P3[0:2, ik] = np.array(self.Gs.pos[n2])
                P3[2, ik] = self.Gs.node[i]['zmax']

                P4[0:2, ik] = np.array(self.Gs.pos[n1])
                P4[2, ik] = self.Gs.node[i]['zmax']
                ik = ik + 1

        d = self.subseg()
        cpt = 0
        subseg = {}
        for k in d.keys():
            for l in d[k]:
                subseg[cpt] = l
                cpt = cpt + 1
                nebr = self.Gs.neighbors(l)
                n1 = nebr[0]
                n2 = nebr[1]
                #print ik,n1,n2

                P1[0:2, ik] = np.array(self.Gs.pos[n1])
                P1[2, ik] = self.Gs.node[l]['ss_zmin']
                #print P1[:,ik]

                P2[0:2, ik] = np.array(self.Gs.pos[n2])
                P2[2, ik] = self.Gs.node[l]['ss_zmin']
                #print P2[:,ik]

                P3[0:2, ik] = np.array(self.Gs.pos[n2])
                P3[2, ik] = self.Gs.node[l]['ss_zmax']
                #print P3[:,ik]

                P4[0:2, ik] = np.array(self.Gs.pos[n1])
                P4[2, ik] = self.Gs.node[l]['ss_zmax']
                #print P4[:,ik]

                ik = ik + 1

#        subseg = self.ce.keys()
#        for j in range(cen):
#            ie = subseg[j]
#            t  = self.tahe[0,ie]
#            h  = self.tahe[1,ie]
#            k  = ii+j
#
#            P1[0:2,k]=self.Gs.pos[t]
#            P1[2,k]=self.Gs.node[ie]['zmin']
#
#            P2[0:2,k]=self.Gs.pos[h]
#            P2[2,k]=self.Gs.node[ie]['zmin']
#
#            P3[0:2,k]=self.Gs.pos[h]
#            P3[2,k]=self.Gs.node[ie]['zmax']
#
#            P4[0:2,k]=self.Gs.pos[t]
#            P4[2,k]=self.Gs.node[ie]['zmax']
#
        npt = 4 * (en + cen)
        _filename,ext = os.path.splitext(self.filename)
        _filegeom = _filename+'.off'
        self.filegeom=_filegeom
        filegeom = pyu.getlong(_filegeom, "geom")
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
        for i in range(en + cen):
            q = 4 * i
            if i < en:
                ne = i + 1
                name = self.Gs.node[ne]['name']
            else:
                nss = i - en
                ne = subseg[nss]
                name = self.Gs.node[ne]['ss_name']

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

    def show3(self, bdis=True):
        """ geomview display of the indoor structure

        Parameters
        ----------
            bdis
                boolean (default True)
        """
        self.geomfile()
        filename = pyu.getlong(self.filegeom, "geom")
        if (bdis):
            #chaine = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
            chaine = "geomview  -b 1 1 1 " + filename + " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)

    def signature(self, iTx, iRx):
        """ Determine signature between node iTx and node iRx

        Parameters
        ----------
        pTx  : np.array (2, 1)
                ransmitter coordinates
        pRx  : np.array (2, 1)
                Receiver coordinates

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

        signature = []

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
                    signature.append(path)
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

        return(sigarr, signature)

    def get_Sg_pos(self, signature):
        """ return position of the signatures
        """
        signature = signature[0][0]
        sposfull = np.zeros((len(signature), 2))
        iz = np.nonzero(signature != 0)[0]
        spos = np.array([self.Gs.pos[i] for i in signature if i != 0])
        sposfull[iz, :] = spos
        return (sposfull)

    def showSig(self, sig, Tx=None, Rx=None, fig=plt.figure(), ax=None):
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

        """

        if fig is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        elif ax is None:
            ax = fig.add_subplot(111)
        lines = []
        ps = self.get_Sg_pos(sig)
        nz = np.nonzero(sig == 0)[0]
        if sr != []:
            mask = np.zeros((2, len(sig)))
            mask[:, nz] = 1
            vertices = np.ma.masked_array(ps.T, mask)
            lines.extend(ax.plot(vertices[0, :], vertices[1, :], color='k'))
        pdb.set_trace()
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
            p
                point
            nroom
                room number of p
        Return
        ------
            dist
                list of distances to walls of room nroom

        Notes
        -----

            Return  dist list which is a list of all the distances to the walls of the room

        """
        pp = Point(p[0], p[1])

        dist = []
        p0_xy = []
        p1_xy = []

        Nc = L.Gr.node[nroom]['cycle']
        vnode = L.Gt.node[Nc]['vnodes']

        for j in range(len(L.Gr[nroom]['vnodes'])):
            nn = L.b_Gr[5]['vnodes'][j]
            nta = G1.tahe[0, nn - 1]
            nhe = G1.tahe[1, nn - 1]
            p0 = np.array([G1.pt[0, nta], G1.pt[1, nta]])
            p1 = np.array([G1.pt[0, nhe], G1.pt[1, nhe]])
            p0_xy.insert(j, p0)
            p1_xy.insert(j, p1)

        pstartwll = np.array(p0_xy)
        pfinwll = np.array(p1_xy)

        for i in range(len(L.b_Gr[nroom]['vnodes'])):
            line_wall = LineString([(pstartwll[i, 0],
                pstartwll[i, 1]), (pfinwll[i, 0], pfinwll[i, 1])])
            dist.insert(i, line_wall.distance(pp))
        return(dist)

    def randTxRx(self):
        """Returns random coordinates for Tx and Rx.

        Returns
        -------
        p_Tx : numpy.ndarray
             A point of the placement of the Tx
        p_Rx : numpy.ndarray
             A point of the placement of the Rx

        Examples
        --------

        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> p_Tx,p_Rx = L.randTxRx()

        Notes
        -----
            ex fn Tx_Rx_pos


        """

        self.boundary()

        Tx_x = rd.uniform(self.ax[0], self.ax[1])
        Tx_y = rd.uniform(self.ax[2], self.ax[3])
        Rx_x = rd.uniform(self.ax[0], self.ax[1])
        Rx_y = rd.uniform(self.ax[2], self.ax[3])

        p_Tx = np.array([Tx_x, Tx_y])
        p_Rx = np.array([Rx_x, Rx_y])

        return(p_Tx, p_Rx)

    def boundary(self, dx=0, dy=0):
        """ set a boundary around layout

        Parameters
        ----------
        dx : float
            x offset (default 0)
        dy : float
            y offset (default 0 )

        self.ax is updated

        Examples
        --------

        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> L.boundary()

        """
        xmax = max(p[0] for p in self.Gs.pos.values())
        xmin = min(p[0] for p in self.Gs.pos.values())
        ymax = max(p[1] for p in self.Gs.pos.values())
        ymin = min(p[1] for p in self.Gs.pos.values())

        self.ax = (xmin - dx, xmax + dx, ymin - dy, ymax + dy)

    def loadGv(self, _fileGv):
        """Load Layout Gv file which contains Gv graph Values and reconstruct Gv graph.

        Parameters
        ----------
        _fileGv : str
        _fileGv is an *.Gv file, which contains all information about graph of visibility (Gv)

        Returns
        -------
        Gv_re : networkx.classes.graph.Graph
            The graph of visibility reconstructed
        pos : dict
            Explicitly set positions of the nodes
        labels : dict
            Node labels in a dictionary keyed by edge two-tuple of text labels
        (default=None), Only labels for the keys in the dictionary are drawn.

        Notes
        -----
         Extension of Gv file is lo (to be changed)
         Those file are placed in the layout directory


        Examples
        --------

        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> ncoin,ndiff = L.buildGc()
        >>> L.buildGt()
        >>> L.buildGr()
        >>> _fileGv = 'exemple.lo'
        >>> Gv_re,pos,labels=L.loadGv( _fileGv)
        >>> assert Gv_re.nodes()[5]== 6,'Mistake'
        >>> a=plt.title('Test Gv loadGv')
        >>> fig,ax = L.showGs()
        >>> nx.draw(Gv_re,pos,node_color='r')
        >>> plt.show()
        >>> plt.clf()


        """
        Gv_re = nx.Graph()
        pos = {}
        labels = {}
        fileGv = pyu.getlong(_fileGv, 'layout')
        if os.path.isfile(fileGv):
            data = io.loadmat(fileGv, appendmat=False)
            # reconstruct Gv from data
            Gv_lab = data['lab']
            pos_val = data['p_val']
            pos_keys = data['p_keys']
            Gv_edges = data['edges']
            Gv_nodes = data['node']

            for i in range(len(pos_keys)):
                Gv_re.add_node(Gv_nodes[i][0])
                pos[pos_keys[i][0]] = (pos_val[i][0], pos_val[i][1])
                labels[Gv_lab[i][0]] = str(Gv_lab[i][0])

            for j in range(len(Gv_edges)):
                Gv_re.add_edge(Gv_edges[j, 0], Gv_edges[j, 1])

        return(Gv_re, pos, labels)

#    def GvF(self, _fileGv):
#        """build GvF graph
#        
#        GvF Graph of visibility between the walls and
#        the furnitures of indoor environment.
#
#        Parameters
#        ----------
#        _fileGv : .lo file which contains graph of visibility (Gv).
#
#        Returns
#        -------
#        GvF : networkx.classes.graph.Graph
#            The graph of visibility between the walls and the furnitures(GvF).
#        posF : dict
#            Explicitly set positions of the nodes
#        labelsF : dict
#            Node labels in a dictionary keyed by edge two-tuple of text labels
#            (default=None), Only labels for the keys in the dictionary are drawn.
#
#        Examples
#        --------
#        >>> L = Layout()
#        >>> L.loadstr('tag.str','simul8.mat','simul8.slab')
#        >>> L.buildGt()
#        >>> L.buildGr()
#        >>> GvF,posF,labelsF=L.GvF('tag.lo')
#        >>> assert posF[-14][0]==8.5,'Mistake'
#
#        """
#
#        visik1, visik2 = self.visi_list_F()
#        posF = {}
#        labelsF = {}
#        GF = nx.Graph()
#        Gv_re, pos_re, labels_re = self.loadGv(_fileGv)
#
#        for i in range(len(visik1)):
#            GF.add_node(visik1[i])
#            posF[visik1[i]] = (self.Gs.pos[visik1[i]][0],
#                self.Gs.pos[visik1[i]][1])
#            labelsF[visik1[i]] = str(visik1[i])
#            GF.add_node(visik2[i])
#            posF[visik2[i]] = (self.Gs.pos[visik2[i]][0],
#                self.Gs.pos[visik2[i]][1])
#            labelsF[visik2[i]] = str(visik2[i])
#            GF.add_edge(visik1[i], visik2[i])
#
#        posF.update(pos_re, **posF)
#        labelsF.update(labels_re, **labelsF)
#        GvF = nx.compose(Gv_re, GF)
#nx.draw(GvF,posF)
#L.showGs()
#show()

        return(GvF, posF, labelsF)

    def saveGv(self, _fileGv):
        """Creates Layout's Gv file which contains Gv graph values.

        Parameters
        ----------
        _fileGv : str
        _fileGv is an *.lo file, in which saves all information about graph of
        visibility (Gv).


        """
        # create lo file
        fileGv = pyu.getlong(_fileGv, 'layout')

        # writing into lo file
        if os.path.isfile(fileGv):
            print fileGv, 'Warning fileGv already exist'
        else:
            print 'create ', fileGv, ' file'
            data = dict(lab=self.labels_new.keys(),
                      node=self.Gv_new.nodes(),
                      edges=self.Gv_new.edges(),
                      p_val=self.pos_new.values(),
                      p_keys=self.pos_new.keys())

            io.savemat(fileGv, data, appendmat=False)

    def savelay(self, _fileGv, _filelay):
        """ Creates Layout's lay file which contains graphs values.

        Parameters
        ----------
        _fileGv : str
            _fileGv is an *.lo file, in which saves all information about graph of
        visibility (Gv).
        _filelay : str
            _fileGv is an *.lay file, in which all information about graphs of
        structure (Gs), topological (Gt), rooms (Gr) and visibility (Gv) saves.

        Examples
        --------
        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> ncoin,ndiff = L.buildGc()
        >>> L.buildGt()
        >>> L.buildGr()
        >>> _fileGv  = 'exemple.lo'
        >>> _filelay = 'example.lay'
        >>> filelay= pyu.getlong(_filelay,'layout')
        >>> assert os.path.isfile(filelay)==True,'already exist'


        """
        Gv, pos, labels = self.loadGv(_fileGv)

        filelay = pyu.getlong(_filelay, 'layout')
        if os.path.isfile(filelay):
            print filelay, ' already exist'
        else:
            print 'create ', filelay, ' file'
            data_Gs = dict(node=self.Gs.nodes(), edges=self.Gs.edges(), p_val=self.Gs.pos.values(), p_keys=self.Gs.pos.keys())
            data_Gt = dict(node=self.Gt.nodes(), edges=self.Gt.edges(), p_val=self.Gt.pos.values(), p_keys=self.Gt.pos.keys())
            data_Gr = dict(node=self.Gr.nodes(), edges=self.Gr.edges(), p_val=self.Gr.pos.values(), p_keys=self.Gr.pos.keys())
            data_Gv = dict(node=Gv.nodes(), edges=Gv.edges(
                ), p_val=pos.values(), p_keys=pos.keys())
            data_graph = dict(Gs=data_Gs, Gt=data_Gt, Gr=data_Gr, Gv=data_Gv)
            cPickle.dump(data_graph, open(filelay, "wb"))

    def loadlay(self, _filelay):
        """Loads Layout's lay file which contains graphs values.

        Parameters
        ----------
        _filelay : str
            _fileGv is an *.lay file, in which contains all information about graphs of
        structure (Gs), topological (Gt), rooms (Gr) and visibility (Gv).

        Returns
        -------
        data_graph : dict
            A dictionary which contains all information about graphs Gs, Gt, Gr, Gv.

        Examples
        --------
        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> ncoin,ndiff = L.buildGc()
        >>> L.buildGt()
        >>> L.buildGr()
        >>> _filelay = 'exemple.lay'
        >>> data_graph=L.loadlay(_filelay)
        """
        filelay = pyu.getlong(_filelay, 'layout')
        data_graph = cPickle.load(open(filelay))

        return(data_graph)

#    def poly_door(self, nnode):
#        """Creates a polygon with the space between a door.
#
#        Decides the position of the doors in indoor environment, than moves the
#        door extremities points from the extremities d = 5cm and builds a polygon.
#        The polygon includes the length of the doors thus a part of not allowed
#        zone of rooms is also included in required part of the Allowed Zone.
#
#        Parameters
#        ----------
#        nnode : int
#            The nnode number of a node of Gt graph.
#
#        Returns
#        -------
#        describe : shapely.geometry.polygon.Polygon
#            A polygon between a door.
#
#        References
#        ----------
#        http://gispython.org/shapely/docs/1.2/manual.html#object.buffer
#
#        ..  [1] Report "Simulation based on ray tracing for the modeling of
#        dynamic wideband channel" - Localization application, p 31, 2011.
#
#        .. image:: _static/AZ_2.png
#            :width: 200
#            :height: 200
#            :scale: 60 %
#            :alt: alternate text
#
#        Examples
#        --------
#        >>> L = Layout()
#        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
#        >>> ncoin,ndiff = L.buildGc()
#        >>> L.buildGt()
#        >>> L.buildGr()
#        >>> nnode=0
#        >>> dilated=L.poly_door(nnode)
#        >>> assert dilated.area==0.78413712263648416,'Mistake'
#
#        """
#
#        nd = self.Gr.node[nnode]['doors'][0]
#        n1, n2 = self.Gs.neighbors(nd)
#        x1, y1 = self.Gs.pos[n1]
#        x2, y2 = self.Gs.pos[n2]
#        deltx = self.delta(x1, x2)
#        delty = self.delta(y1, y2)
#
#        if deltx > delty:
#            x1, x2 = self.new_coor_d(x1, x2)
#            y1 = y1
#            y2 = y2
#        if deltx < delty:
#            y1, y2 = self.new_coor_d(y1, y2)
#            x1 = x1
#            x2 = x2
#
#        line = sh.LineString([(x1, y1), (x2, y2)])
#        dilated = line.buffer(0.5)
#        return(dilated)

#    def moby_allowed_zone(self):
#        """Defineds Allowed Mobility Zones for an entire indoor environment.
#
#        By adding the Allowed Mobility Zones in the rooms to the Allowed
#        Mobility Zones within the doors, there is the complete required Allowed
#        Mobility Zones within entire indoor environment.
#
#        Returns
#        -------
#        describe : shapely.geometry.multipolygon.MultiPolygon
#            Allowed Mobility Zones for an entire indoor environment
#
#        References
#        ----------
#        http://gispython.org/shapely/docs/1.2/manual.html#object.union
#
#        ..  [1] Report "Simulation based on ray tracing for the modeling of
#        dynamic wideband channel" - Localization application, p 31, 2011.
#
#        .. image:: _static/AZ_3.png
#            :width: 200
#            :height: 200
#            :scale: 60 %
#            :alt: alternate text
#
#        Examples
#        --------
#        >>> L = Layout()
#        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
#        >>> ncoin,ndiff = L.buildGc()
#        >>> L.buildGt()
#        >>> L.buildGr()
#        >>> u=L.moby_allowed_zone()
#        >>> assert u.area== 39.799999999999997,'Mistake'
#
#        """
#
#        polygon1 = []
#        dilated = []
#        for i in range(len(self.Gt.node)):
#            polygon1.append(self.polygon_delta(i))
#
#        for j in range(len(self.Gr.node)):
#            dilated.append(self.poly_door(j))
#
#        polygons = polygon1 + dilated
#        u = cascaded_union(polygons)
#        return(u)

    def segs_cross(self, line):
        """Gives segment/subsegment numbers crossed by a line.

        Parameters
        ----------
        line : shapely.geometry.linestring.LineString
            The line formed by the Tx and Rx coordinates.

        Returns
        -------
        segs_cross : list
        Segment/subsegment numbers crossed by a line.

        References
        ----------
        ..  [1] Report "Simulation based on ray tracing for the modeling of
        dynamic wideband channel" - Localization application, pp 44 - 45, 2011.

        Examples
        --------
        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> ncoin,ndiff = L.buildGc()
        >>> L.buildGt()
        >>> L.buildGr()
        >>> _fileGv='exemple.lo'
        >>> p_Tx=np.array([2,0])
        >>> p_Rx=np.array([7.5,0])
        >>> Nc=10
        >>> color_all=pyu.randcol(Nc)
        >>> color=color_all[0]
        >>> signature,num_seg_val=L.filtre_sign(p_Tx,p_Rx,_fileGv,color)
        >>> kuku=num_seg_val[0]
        >>> Tx_x,Tx_y = p_Tx
        >>> Rx_x,Rx_y = p_Rx
        >>> ptx=sh.Point(Tx_x, Tx_y)
        >>> prx=sh.Point(Rx_x, Rx_y)
        >>> Tx_im_x, Tx_im_y=L.image_M_Tx(p_Tx,kuku)
        >>> line=sh.LineString([(Tx_im_x, Tx_im_y),(Rx_x, Rx_y)])
        >>> p_int=L.cross_shortest(line)
        >>> line1=sh.LineString([(Tx_x, Tx_y),(p_int[0].x,p_int[0].y),(Rx_x, Rx_y)])
        >>> sign=L.segs_cross(line1)
        >>> segs_cross=L.segs_cross(line)

        """
        segs_cross = []

        nodes = self.Gs.nodes()

        for i in range(len(nodes)):
            if nodes[i] >= 0:
                n1, n2 = self.Gs.neighbors(nodes[i])
                    #+1 for coherancy with PyRay result
                line_i = sh.LineString([(self.Gs.pos[n1][0],
                                         self.Gs.pos[n1][1]),
                                        (self.Gs.pos[n2][0],
                                         self.Gs.pos[n2][1])])
                if line.crosses(line_i):
                    segs_cross.append(nodes[i])

        return(segs_cross)

    def sign_cross(self, pTx, pRx, pt, seg_num):
        """Returns some possible signatures which have at last one reflexion.

        Parameters
        ----------
        pTx : numpy.ndarray
            The point of the placement of the Tx
        pRx : numpy.ndarray
            The point of the placement of the Rx
        pt : numpy.ndarray
            The point of the placement of a Reflection
        seg_num : list
            Segment/subsegment numbers crossed by a line formed by Tx and Rx coordinates.

        Returns
        -------
        sign : list
            Some possible signature at last with one reflexion.

        Raises
        ------
        Need to be developed more.

        References
        ----------
        ..  [1] Report "Simulation based on ray tracing for the modeling of
                dynamic wideband channel" - Localization application, p 45, 2011.

        Examples
        --------
        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> ncoin,ndiff = L.buildGc()
        >>> L.buildGt()
        >>> L.buildGr()
        >>> _fileGv='exemple.lo'
        >>> pTx=np.array([2,0])
        >>> pRx=np.array([7.5,0])
        >>> p_Tx=pTx
        >>> p_Rx=pRx
        >>> Nc=10
        >>> color_all=pyu.randcol(Nc)
        >>> color=color_all[0]
        >>> signature,num_seg_val=L.filtre_sign(p_Tx,p_Rx,_fileGv,color)
        >>> kuku=num_seg_val[0]
        >>> Tx_x,Tx_y = p_Tx
        >>> Rx_x,Rx_y = p_Rx
        >>> ptx=sh.Point(Tx_x, Tx_y)
        >>> prx=sh.Point(Rx_x, Rx_y)
        >>> Tx_im_x, Tx_im_y=L.image_M_Tx(p_Tx,kuku)
        >>> line=sh.LineString([(Tx_im_x, Tx_im_y),(Rx_x, Rx_y)])
        >>> p_int=L.cross_shortest(line)
        >>> pt=p_int[0]
        >>> line1=sh.LineString([(Tx_x, Tx_y),(pt.x,pt.y),(Rx_x, Rx_y)])
        >>> sign=L.segs_cross(line1)
        >>> segs_cross=L.segs_cross(line)
        >>> seg_num=segs_cross[0]
        >>> sign=L.sign_cross(pTx,pRx, pt,seg_num)

        """

        Tx_x, Tx_y = pTx
        Rx_x, Rx_y = pRx

        px = pt.x
        py = pt.y

        line_Tx_p = sh.LineString([(Tx_x, Tx_y), (px, py)])
        line_Rx_p = sh.LineString([(px, py), (Rx_x, Rx_y)])
        segs_Tx_p = self.segs_cross(line_Tx_p)
        segs_Rx_p = self.segs_cross(line_Rx_p)
        dist_Tx = []
        dist_Rx = []
        sign = []
        for i in range(len(segs_Tx_p)):
            d = sh.Point(Tx_x, Tx_y).distance(sh.Point(self.Gs.pos[
                segs_Tx_p[i]][0], self.Gs.pos[segs_Tx_p[i]][1]))
            dist_Tx.append(d)
        distTx = np.argsort(dist_Tx)
        for j in range(len(segs_Rx_p)):
            la = sh.Point(self.Gs.pos[segs_Rx_p[j]][0], self.Gs.pos[
                segs_Rx_p[j]][1]).distance(sh.Point(Tx_x, Tx_y))
            dist_Rx.append(la)
        distRx = np.argsort(dist_Rx)
        #sign.append(0) # Modif 28/10/11
        for k in range(len(segs_Tx_p)):
            sign.append(segs_Tx_p[distTx[k]])
        #sign.append(segs_cross[0]) # Modif 28/10/11
        sign.append(seg_num)
        for h in range(len(segs_Rx_p)):
            sign.append(segs_Rx_p[distRx[h]])
        #sign.append(0)
        return (sign)

#    def cross_shortest(self, line):
#        """Returns all points crossed by the line.
#
#        Returns all points crossed by the line which cross the shortest ray
#        in the middle of the shortest ray and they form an angle of 90 degree.
#
#        Parameters
#        ----------
#        line : shapely.geometry.linestring.LineString
#            A line which forms an angle of 90 degree with the line formed by
#        the Tx and Rx coordinates.
#
#        Raises
#        ------
#        Need to be developed more.
#
#        Returns
#        -------
#        p_int : type
#            Points crossed by the line
#
#        References
#        ----------
#        ..  [1] Report "Simulation based on ray tracing for the modeling of
#        dynamic wideband channel" - Localization application, p 45, 2011.
#
#        Examples
#        --------
#        >>> L = Layout()
#        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
#        >>> ncoin,ndiff = L.buildGc()
#        >>> L.buildGt()
#        >>> L.buildGr()
#        >>> _fileGv='exemple.lo'
#        >>> p_Tx=np.array([2,0])
#        >>> p_Rx=np.array([7.5,0])
#        >>> Nc=10
#        >>> color_all=pyu.randcol(Nc)
#        >>> color=color_all[0]
#        >>> signature,num_seg_val=L.filtre_sign(p_Tx,p_Rx,_fileGv,color)
#        >>> kuku = num_seg_val[0]
#        >>> Tx_x,Tx_y = p_Tx
#        >>> Rx_x,Rx_y = p_Rx
#        >>> ptx=sh.Point(Tx_x, Tx_y)
#        >>> prx=sh.Point(Rx_x, Rx_y)
#        >>> Tx_im_x, Tx_im_y=L.image_M_Tx(p_Tx,kuku)
#        >>> line=sh.LineString([(Tx_im_x, Tx_im_y),(Rx_x, Rx_y)])
#        >>> p_int=L.cross_shortest(line)
#        >>> line1=sh.LineString([(Tx_x, Tx_y),(p_int[0].x,p_int[0].y),(Rx_x, Rx_y)])
#        >>> sign=L.segs_cross(line1)
#        >>> p_int=L.cross_shortest(line1)
#        >>> assert p_int[0].x== 4.75,'Mistake'
#
#        """
#
#        segs_cross = self.segs_cross(line)
#
#        p_int = []
#        for k in range(len(segs_cross)):
#            n1, n2 = self.Gs.neighbors(segs_cross[k])
#            line_k = sh.LineString([(self.Gs.pos[n1][0], self.Gs.pos[n1]
#                [1]), (self.Gs.pos[n2][0], self.Gs.pos[n2][1])])
#            p_int.append(line.intersection(line_k))
#
#        return(p_int)

    def plot_sign2(self, p_Tx, p_Rx):
        """Displays a lot of possible signatures.

        Displays a lot of possible signatures but not all.

        Parameters
        ----------
        p_Tx : numpy.ndarray
            The point of the placement of the Tx
        p_Rx : numpy.ndarray
            The point of the placement of the Rx

        References
        ----------
        ..  [1] Report "Simulation based on ray tracing for the modeling of
        dynamic wideband channel" - Localization application, pp 44 - 45, 2011.

        Examples
        --------
        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> p_Tx=np.array([2,0])
        >>> p_Rx=np.array([7.5,0])
        >>> Nc=10
        >>> color_all=pyu.randcol(Nc)
        >>> color=color_all[0]
        >>> L.plot_sign2(p_Tx,p_Rx)

        """
        Nc = 10
        color_all = pyu.randcol(Nc)
        color = color_all[0]
        fig = plt.figure(1, dpi=90)
        ax = fig.add_subplot(1, 1, 1)
        Tx_x, Tx_y = p_Tx
        Rx_x, Rx_y = p_Rx

        xmin = min(Tx_x, Rx_x)
        xmax = max(Tx_x, Rx_x)
        ymin = min(Tx_y, Rx_y)
        ymax = max(Tx_y, Rx_y)
        seg = self.clip(xmin, xmax, ymin, ymax)
        line_sh = sh.LineString([(Tx_x, Tx_y), (Rx_x, Rx_y)])
        for i in range(len(seg)):
            p0_i = self.Gs.pos[seg[i] + 1]
            line = sh.LineString(
                [(Tx_x, Tx_y), (p0_i[0], p0_i[1]), (Rx_x, Rx_y)])
            geu.plot_coords3(ax, line, color)
            geu.plot_bounds3(ax, line, color)
            geu.plot_line3(ax, line, color)

        fig = plt.figure(1, dpi=90)
        ax = fig.add_subplot(1, 1, 1)

    def points_image(self, p_Tx):
        """Returns all image points of the bases

         Parameters
         ----------
         p_Tx : numpy.ndarray
             A point of the placement of the Tx

         Returns
         -------
         Tx_im : list
              A list of tuples: containing the coordinates of image points.
         nodes_posi : list
              A list of positive numbers of the graph structure.

        Examples
        --------
        >>> L = Layout()
        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
        >>> L.buildGt()
        >>> L.buildGr()
        >>> p_Tx=np.array([2,0])
        >>> Tx_im,nodes_posi=L.points_image(p_Tx)
        >>> assert Tx_im[0][0]==2,'Mistake'

        """

        nodes_posi = self.nodes_posi()
        Tx_im = []
        for mk in range(len(nodes_posi)):
            num_seg = nodes_posi[mk]
            Tx_x, Tx_y = p_Tx
            n1, n2 = self.Gs.neighbors(num_seg)
            pa_x, pa_y = self.Gs.pos[n1]
            pb_x, pb_y = self.Gs.pos[n2]
            delx = abs(abs(pa_x) - abs(pb_x))
            dely = abs(abs(pa_y) - abs(pb_y))
            dist = sh.Point(Tx_x, Tx_y).distance(
                sh.LineString([(pa_x, pa_y), (pb_x, pb_y)]))
            if delx <= dely:
                if max(pa_x, pb_x) <= Tx_x:
                    Tx_x_im = Tx_x - 2 * dist
                else:
                    Tx_x_im = Tx_x + 2 * dist
                Tx_y_im = Tx_y
                Tx_im.append((Tx_x_im, Tx_y_im))
            if delx > dely:
                if max(pa_y, pb_y) <= Tx_y:
                    Tx_y_im = Tx_y - 2 * dist
                else:
                    Tx_y_im = Tx_y + 2 * dist
                Tx_x_im = Tx_x
                Tx_im.append((Tx_x_im, Tx_y_im))

        return(Tx_im, nodes_posi)

#    def filtre_sign(self, p_Tx, p_Rx, _filelo, color):
#        """Filtering and validating of signatures.
#
#        Parameters
#        ----------
#        p_Tx : numpy.ndarray
#            Tx coordinates
#        p_Rx : numpy.ndarray
#            Rx coordinates
#        _filelo : str
#            The str file, which contains the graph of visibility.
#
#        Returns
#        -------
#        signature : list
#                A list which contains found signatures
#        num_seg_val : list
#                The wall numbers the ray has the first interaction with
#
#        Examples
#        --------
#        >>> L=Layout()
#        >>> L.loadstr('exemple.str','simul8.mat','simul8.slab')
#        >>> ncoin,ndiff=L.buildGc()
#        >>> L.buildGt()
#        >>> L.buildGr()
#        >>> Nc=10
#        >>> color_all = pyu.randcol(Nc)
#        >>> color = color_all[0]
#        >>> p_Tx  = np.array([2,0])
#        >>> p_Rx  = np.array([7.5,0])
#        >>> _filelo = 'exemple.lo'
#        >>> signature,num_seg_val=L.filtre_sign(p_Tx,p_Rx,_filelo,color)
#
#
#        """
#
#        Nc = 10
#        color_all = pyu.randcol(Nc)
#        color = color_all[0]
#
#        Gsnodes = self.Gs.nodes()
#        nodes_posi = []
#        for i in range(len(Gsnodes)):
#            if Gsnodes[i] >= 0:
#                nodes_posi.append(Gsnodes[i])
#            else:
#                pass
#        all_sign = []  # validation of signachure
#        signature = []  # correct siractures
#        num_seg_val = []  # validated first interaction
#        for i in range(len(nodes_posi)):
#            num_seg = nodes_posi[i]
#            sign = self.find_more_sign(p_Tx, p_Rx, num_seg, color)
#            all_sign.append(sign)
#        for i in range(len(all_sign)):
#            if all_sign[i] != 0:
#                signature.append(all_sign[i])
#                num_seg_val.append(nodes_posi[i])
#            else:
#                pass
#
#        return(signature, num_seg_val)

#    def image_M_Tx(self, p_Tx, num_seg):
#        """Calculation of point image.
#
#        Calculates point image for a segment situated between a transmitter and a receiver.
#
#        Parameters
#        ----------
#        p_Tx : numpy.ndarray
#            The point of the placement of the Tx.
#        num_seg : int
#            The the segment number by which relationship the image point is built.
#
#        Returns
#        -------
#        Tx_im : tuples
#            The coordinates of images point of the Tx calculated for a segment.
#
#        Examples
#        --------
#        >>> L = Layout()
#        >>> L.loadstr('exemple.str')
#        >>> ndiff,ncoin = L.buildGc()
#        >>> L.buildGt()
#        >>> L.buildGr()
#        >>> _filelo='exemple.lo'
#        >>> p_Tx=np.array([2,0])
#        >>> p_Rx=np.array([7.5,0])
#        >>> Nc=10
#        >>> color_all=pyu.randcol(Nc)
#        >>> color=color_all[0]
#        >>> signature,num_seg_val=L.filtre_sign(p_Tx,p_Rx,_filelo,color)
#        >>> Tx_im=L.image_M_Tx( p_Tx,num_seg_val[0])
#
#        """
#
#        Tx_im = []
#        Tx_x, Tx_y = p_Tx
#        n1, n2 = self.Gs.neighbors(num_seg)
#        pa_x, pa_y = self.Gs.pos[n1]
#        pb_x, pb_y = self.Gs.pos[n2]
#        a1 = np.sqrt((Tx_x - pa_x) ** 2 + (Tx_y - pa_y) ** 2)
#        a2 = np.sqrt((pb_x - Tx_x) ** 2 + (pb_y - Tx_y) ** 2)
#        seg_line = sh.LineString([(pa_x, pa_y), (pb_x, pb_y)])
#        m1 = 2 * (pa_x - pb_x)
#        m2 = a2 ** 2 - a1 ** 2 + pa_x ** 2 - pb_x ** 2 + pa_y ** 2 - pb_y ** 2
#        m3 = 2 * (pb_y - pa_y)
#        #x=(m2/m1+m3*y/m1)
#        k1 = m2 / (1. * m1)
#        k2 = m3 / (1. * m1)
#        #x=k1+k2*y
#        f1 = 1 + k2 ** 2
#        f2 = 2 * (k1 * k2 - k2 * pb_x - pb_y)
#        f3 = k1 ** 2 - a2 ** 2 - 2 * k1 * pb_x + pb_x ** 2 + pb_y ** 2
#        #f1*y**2+f2*y+f3=0
#        b = f2 / (1. * f1)
#        ck = f3 / (1. * f1)
#        #y**2+b*y+c=0
#        di = b ** 2 - 4 * ck
#        delx = pb_x - pa_x
#        dely = pb_y - pa_y
#        if delx > dely:
#            y = (-b + np.sqrt(di)) / 2
#            x = k1 + k2 * y
#        else:
#            y = (-b - np.sqrt(di)) / 2
#            x = k1 + k2 * y
#        Tx_im.append((x, y))
#
#        return(Tx_im[0])


if __name__ == "__main__":
    doctest.testmod()
