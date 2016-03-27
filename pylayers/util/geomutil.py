# -*- coding:Utf-8 -*-
"""
.. currentmodule:: pylayers.util.geomutil

=================================================
Geometry Module (:mod:`pylayers.util.geomutil`)
================================================

Geomview Class
==============

.. autosummary::
    :toctree: generated/

Geomview.__init__
Geomview.show3

Geomlist Class
==============

.. autosummary::
    :toctree: generated/

    Geomlist.__init__
    Geomlist.append

GeomVect Class
==============

.. autosummary::
    :toctree: generated/

    GeomVect.__init__
    GeomVect.segments
    GeomVect.geomBase
    GeomVect.points

Geomoff Class
=============

.. autosummary::
    :toctree: generated/

    Geomoff.__init__
    Geomoff.loadpt
    Geomoff.savept
    Geomoff.polygon
    Geomoff.polygons
    Geomoff.cylinder
    Geomoff.box
    Geomoff.pattern

Plot_Shapely Class
==================

.. autosummary::
    :toctree: generated/

    Plot_shapely.__init__
    Plot_shapely.plot_coords
    Plot_shapely.plot_ligne
    Plot_shapely.plot_polygon
    Plot_shapely.plot_multi

LineString Class
==================

.. autosummary::
    :toctree: generated/

    LineString.__init__
    LineString.plot

PolyGon Class
=============

.. autosummary::
    :toctree: generated/

    Polygon.__init__
    Polygon.plot
    Polygon.__init__
    Polygon.__add__
    Polygon.__repr__
    Polygon.ndarray
    Polygon.signedarea
    Polygon.plot
    Polygon.simplify
    Polygon.buildGv
    Polygon.showGv
    Polygon.ptconvex


Utility Functions
=================

.. autosummary::
    :toctree: generated/

     angular
     SignedArea
     Centroid
     Lr2n
     isBetween
     pvec
     pvecn
     onb
     vec_sph
     ellipse
     normalize
     ptonseg
     dptseg
     linet
     ccw
     intersect
     isaligned
     isleft
     isleftorequal
     affine
     cylmap
     mul3
     MRot3
     MEulerAngle
     SphericalBasis
     angledir
     BTB_rx
     BTB_tx

     plot_coords
     plot_bounds
     plot_line
     v_color

     plotPolygon
     shrinkPolygon
     shrinkPolygon2
     simplifyPolygon
     wall_delta
     plot_coords2
     plot_bounds2
     plot_line2
     plot_coords3
     plot_bounds3
     plot_line3
     valid_wedge
     sector
     dist
     line_intersection
     linepoly_intersection
     mirror
     distseg
     dmin3d



"""
from pylayers.util import easygui
from pylayers.antprop.slab import Slab, SlabDB, Mat, MatDB
import shapely.geometry as sh
import scipy.linalg as la
import pdb
import logging
import networkx as nx
import doctest
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from   matplotlib.patches import Circle
import scipy as sp
import numpy as np
from   scipy.linalg import toeplitz, triu
from   pylayers.util.project import *
import pylayers.util.pyutil as pyu
import pylayers.util.graphutil as gru
#from antenna import *
import shapely.geometry as shg
import shapely.ops as sho
from   descartes.patch import PolygonPatch
from   shapely.wkt import loads
from   itertools import combinations,permutations


COLOR = {
    True: '#6699cc',
    False: '#ff3333'
}

class Plot_shapely(PyLayers):
    """draw Shapely with matplotlib - pylab
     Plot_shapely.py
     Author : Martin Laloux 2010

    """

    def __init__(self, obj, ax, coul=None, alph=1):
        """ object constructor

         Parameters
         ----------
         ax :
             pylab Axes
         obj  : geometric object
         coul : matplotlib color
         alph : transparency

         Examples
         --------

         >>> from shapely.wkt import loads
         >>> import matplotlib.pylab as plt
         >>> ax = plt.gca()
         >>> ligne = loads('LINESTRING (3 1, 4 4, 5 5, 5 6)')
         >>> a  = Plot_shapely(ligne,ax,'r', 0.5)
         >>> a.plot
         >>> Plot_shapely(ligne,ax,'#FFEC00').plot
         >>> plt.show()
        """
        self.obj = obj
        self.type = obj.geom_type
        self.ax = ax
        self.coul = coul
        self.alph = alph

    def plot_coords(self):
        """ points
        """
        x, y = self.obj.xy
        self.ax.plot(x, y, 'o', color=self.coul)

    def plot_ligne(self):
        """lines"""
        x, y = self.obj.xy
        self.ax.plot(x, y, color=self.coul, alpha=self.alph, linewidth=3)

    def plot_polygon(self):
        """polygons"""
        patch = PolygonPatch(self.obj, facecolor=self.coul,
                             edgecolor='#000000', alpha=self.alph)
        self.ax.add_patch(patch)

    def plot_multi(self):
        """multipoints, multilignes,multipolygones + GeometryCollection"""
        for elem in self.obj:
            Plot_shapely(elem, self.ax, self.coul, self.alph).plot

    @property
    def plot(self):
        """draw w.r.t geometrical type"""
        if self.type == 'Point':
            self.plot_coords()
        elif self.type == 'Polygon':
            self.plot_polygon()
        elif self.type == 'LineString':
            self.plot_ligne()
        elif "Multi" in self.type:
            """ex. MultiPolygon"""
            self.plot_multi()
        elif self.type == 'GeometryCollection':
            self.plot_multi()
        elif self.type == 'LinearRing':
            self.plot_line()
        else:
            raise ValueError("unknown: %s" % self.type)


class LineString(PyLayers,shg.LineString):
    """ Overloaded shapely LineString class
    """
    def __init__(self,p):

        if type(p) == shg.polygon.Polygon:
            self.Np = shape(p.exterior.xy)[1] - 1
            shg.LineString.__init__(self, p)

        if type(p) == shg.multipoint.MultiPoint:
            self.Np = np.shape(p)[0]
            shg.LineString.__init__(self, p)

        if type(p) == list:
            p = np.array(p)

        if type(p) == np.ndarray:
            self.Np = np.shape(p)[1]
            tp = []
            for k in range(self.Np):
                tp.append(p[:, k])
            tp.append(tp[0])
            tu = tuple(tp)
            shg.LineString.__init__(self, tu)

    def plot(self,**kwargs):
        """ plot LineString

        Parameters
        ----------

        show : boolean
        fig : figure object
        ax  : axes object
        linewidth : int 
        color :  string
            default #abcdef"
        alpha :  float
            transparency   (default 0.8)
        figsize : tuple

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.geomutil import *
            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> l1 = np.array([[0,1,1,0],[0,0,1,1]])
            >>> L1 = LineString(l1)
            >>> l2 = [[3,4,4,3],[1,1,2,2]]
            >>> L2 = LineString(l2)
            >>> fig,ax = L1.plot(color='red',alpha=0.3,linewidth=3)
            >>> fig,ax = L2.plot(fig=fig,ax=ax,color='blue',alpha=0.7,linewidth=2)
            >>> title = plt.title('test plotting LineString')

        """

        defaults = {'show': False,
                'fig': [],
                'ax': [],
                'color':'#abcdef',
                'linewidth':1,
                'alpha':0.8 ,
                'figsize':(10,10)
                 }
        #
        # update default values
        #
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value
        #
        # getting fig and ax
        #
        if kwargs['fig'] == []:
            fig = plt.figure(figsize=kwargs['figsize'])
            fig.set_frameon(True)
        else:
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            ax = fig.gca()
        else:
            ax = kwargs['ax']

        x, y = self.xy

        ax.plot(x, y,
                color = kwargs['color'],
                alpha=kwargs['alpha'],
                linewidth = kwargs['linewidth'])

        if kwargs['show']:
            plt.show()

        return fig,ax

#-----------------------------------------------------------
#   Functions used for calculation of visibility graph Gv
#-----------------------------------------------------------

class Polygon(PyLayers,shg.Polygon):
    """ Overloaded shapely Polygon class

    Attributes
    ----------

    Methods
    -------
    plot
    ptconvex
    buildGv
    ndarray :
        get a ndarray from a Polygon
    signedarea :
        get the signed area of the polygon

    """
    def __init__(self, p=[[3, 4, 4, 3], [1, 1, 2, 2]], vnodes=[],delta=0):
        """ object constructor

        Parameters
        ----------

        p : list
            2xNp np.array
            shg.MultiPoint
            shg.Polygon
            tuple : self.ax
        vnodes : list of alternating points and segments numbers
            default = [] in this case a regular ordered sequence
            is generated.

        Notes
        -----

        Convention : a Polygon as an equal number of points and segments
        There is an implicit closure between first and last point

        """

        if type(p) == sh.multipolygon.MultiPolygon:
            raise AttributeError('MultiPolygon are not allowed')

        if type(p) == shg.polygon.Polygon:
            self.Np = np.shape(p.exterior.xy)[1] - 1
            p = np.vstack((p.exterior.xy[0][0:-1],p.exterior.xy[1][0:-1]))
            #shg.Polygon.__init__(self, pt)
            #

        if type(p) == tuple:
            xmin = p[0]-delta
            xmax = p[1]+delta
            ymin = p[2]-delta
            ymax = p[3]+delta
            p = [[xmin,xmin,xmax,xmax],[ymin,ymax,ymax,ymin]]

        if type(p) == shg.multipoint.MultiPoint:
            self.Np = np.shape(p)[0]
            shg.Polygon.__init__(self, p)

        if type(p) == list:
            p = np.array(p)

        if type(p) == np.ndarray:
            if np.shape(p)[1] == 2:
                p = p.T
            self.Np = np.shape(p)[1]
            tp = []
            for k in range(self.Np):
                tp.append(p[:, k])
            tp.append(tp[0])
            tu = tuple(tp)
            shg.Polygon.__init__(self, tu)

        self.Np = np.shape(self.exterior.xy)[1] - 1

        if vnodes != []:
            self.vnodes = np.array(vnodes)
            # check if always True
            # very important fic for buildGv
            # now vnodes starts always with <0
            if self.vnodes[0]>0:
                self.vnodes = np.roll(self.vnodes,-1)
        else:
            # create sequence
            #
            # -1 1 -2 2 -3 3 ... -(Np-1) (Np-1)
            #
            u = np.array([-1, 1])
            v = np.arange(self.Np) + 1
            self.vnodes = np.kron(v, u)

    def __add__(self,p):
        """ add 2 polygons

        Parameters
        ----------

        p : Polygon

        Returns
        -------

        pm : merged polygon or unchanged polygon

        """
        pnew = self.union(p)
        # v0   = self.vnodes
        #v1   = p.vnodes
        #nseg0 = filter(lambda x:x>0,v0)
        #nseg1 = filter(lambda x:x>0,v1)

        #commseg = np.intersect1d(nseg0,nseg1)[0]

        #is0 = np.where(nseg0==commseg)[0][0]
        #is1 = np.where(nseg1==commseg)[0][0]

        #rs0 = np.roll(v0,2*is0-1)[1:]
        #rs1 = np.roll(v1,2*is1-1)[1:]
        #if rs1[0]==rs0[0]:
        #    rs1=rs1[::-1]

        #print rs0
        #print rs1
        #assert(rs0[0]==rs1[-1])
        #assert(rs0[-1]==rs1[0])
        #vnodes = np.hstack((rs0,rs1[1:-1]))
        #self.vnodes = vnodes
        #p2 = Polygon(pnew,vnodes=vnodes)
        p2 = Polygon(pnew)
        #
        # Not finished
        #
        return(p2)

        #p1 = np.vstack((pnew.exterior.xy[0],pnew.exterior.xy[1]))
        #p2 = Polygon(p1)
        #return(p2)
        #if isinstance(pnew,sh.polygon.Polygon):
        #    p1 = np.vstack((pnew.exterior.xy[0],pnew.exterior.xy[1]))
        #    return(p2)
        #else:
        #    pdb.set_trace()
        #    return(self)

    def __repr__(self):
        st = ''
        p = self.ndarray()
        sh = np.shape(p)
        for k in range(sh[1]):
            st = st + '('+str(p[0,k])+','+str(p[1,k])+')\n'

        # vnodes to link with external nodes numerotation
        st = st + '\nvnodes : ('
        for k in range(len(self.vnodes)):
            st = st + str(self.vnodes[k])+' '
        st = st+')\n'

        return(st)


    @property
    def xy(self):
        return self._xy

    @xy.setter
    def xy(self,xy):
        self._xy = xy

    @xy.getter
    def xy(self):
        return self._xy


    def setvnodes(self,L):
        """ update vnodes member from Layout

        Parameters
        ----------

        L : pylayers.layout.Layout

        See Also
        --------

        pylayers.layout.Layout.ispoint

        """
        x,y = self.exterior.xy
        # npts = map(lambda x :
        #            L.ispoint(np.array(x),tol=0.01),zip(x[0:-1],y[0:-1]))
        npts = [L.ispoint(np.array(xx),tol=0.01) for xx in zip(x[0:-1],y[0:-1])]
        seg = zip(npts,np.roll(npts,-1))
        try:
            nseg = map(lambda x : L.numseg(x[0],x[1]),seg)
        except:
            import ipdb
            ipdb.set_trace()
        vnodes = np.kron(npts,np.array([1,0]))+np.kron(nseg,np.array([0,1]))
        self.vnodes = vnodes

    def ndarray(self):
        """ get a ndarray from a Polygon

        Returns
        -------
            p : ndarray (2xNp)

        Examples
        --------
        >>> from pylayers.util.geomutil import *
        >>> p1 = np.array([[0,1,1,0],[0,0,1,1]])
        >>> P1 = Polygon(p1)

        """
        lring = self.exterior
        x, y = lring.xy
        p = np.array([x[0:-1], y[0:-1]])
        return(p)


    def signedarea(self):
        """ get the signed area of the polygon

        """
        p = self.ndarray()
        return sum(np.hstack((p[0, 1::], p[0, 0:1])) * (np.hstack((p[1, 2::], p[1, 0:2])) - p[1, :])) / 2.



    def coorddeter(self):
        """ determine polygon coordinates
        """

        self.xy = np.array([self.exterior.xy[0],self.exterior.xy[1]])

    def isconvex(self,tol = 1e-2):
        """ Determine if a polygon is convex

        Parameters
        ----------
        tol : tolerence on aligned point

        Returns
        -------
         True if convex

        Notes
        -----

        the algorithm tests all triplet of point and L.determine 
        if the third point is left to the 2 first.
        a tolerance can be introduce in cases where the polygon is 
        almost convex.

        """
        self.coorddeter()
        p = self.xy[:,:-1]
        a = p
        b = np.roll(p,1,axis=1)
        c = np.roll(p,2,axis=1)
        return ( np.sum(isleft(a,b,c,tol=tol)) == 0 ) or \
                (np.sum(isleft(c,b,a,tol=tol)) == 0)

    def coorddeter(self):
        """ determine polygon coordinates
        """

        self.xy = np.array([self.exterior.xy[0],self.exterior.xy[1]])

    def isconvex(self,tol = 1e-2):
        """ Determine if a polygon is convex

        Parameters
        ----------
        tol : tolerence on aligned point

        Returns
        -------
        True if convex

        Notes
        -----

        the algorithm tests all triplet of point and L.determine
        if the third point is left to the 2 first.
        a tolerance can be introduce in cases where the polygon is
        almost convex.

        """
        self.coorddeter()
        p = self.xy[:,:-1]
        a = p
        b = np.roll(p,1,axis=1)
        c = np.roll(p,2,axis=1)
        return ( np.sum(isleft(a,b,c,tol=tol)) == 0 ) or \
                (np.sum(isleft(c,b,a,tol=tol)) == 0)


    def plot(self,**kwargs):
        """ plot function

        Parameters
        ----------

        color :  string
            default #abcdef"
        alpha :  float
            transparency   (default 0.8)
        vnodes : bool
            display vnodes
        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.geomutil import *
            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> p1 = np.array([[0,1,1,0],[0,0,1,1]])
            >>> P1 = Polygon(p1)
            >>> p2 = [[3,4,4,3],[1,1,2,2]]
            >>> P2 = Polygon(p2)
            >>> p3 = [np.array([10,10]),np.array([11,10]),np.array([11,11]),np.array([10,11])]
            >>> P3 = Polygon(p3)
            >>> fig,ax = P1.plot(color='red',alpha=0.3)
            >>> fig,ax = P2.plot(fig=fig,ax=ax,color='blue',alpha=0.7)
            >>> fig,ax = P3.plot(fig=fig,ax=ax,color='green',alpha=1)
            >>> title = plt.title('test plotting polygons')

        """

        defaults = {'show': False,
                'fig': [],
                'ax': [],
                'vnodes':False,
                'color':'#abcdef',
                'edgecolor':'#000000',
                'alpha':0.8 ,
                'figsize':(10,10)
                 }
        #
        # update default values
        #
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value
        #
        # getting fig and ax
        #
        if kwargs['fig'] == []:
            fig = plt.figure(figsize=kwargs['figsize'])
            fig.set_frameon(True)
        else:
            fig = kwargs['fig']

        if kwargs['ax'] == []:
            ax = fig.gca()
        else:
            ax = kwargs['ax']

        x, y  = self.exterior.xy
        numpt = filter(lambda x: x <0,self.vnodes)

        ax.fill(x, y,
                color = kwargs['color'],
                alpha = kwargs['alpha'],
                ec = kwargs['edgecolor'])
        if kwargs['vnodes']:
            for k in range(len(numpt)):
                ax.text(x[k],y[k],numpt[k])

        if kwargs['show']:
            plt.show()

        return fig,ax

    def coorddeter(self):
        """ determine polygon coordinates
        """

        self.xy = np.array([self.exterior.xy[0],self.exterior.xy[1]])

    def simplify(self):
        """ simplify polygon - suppress adjacent colinear segments

        Returns
        -------
            poly2 : simplified polygon

        Examples
        --------

        Before


        After


        """
        p = np.array(self.exterior.xy)
        N = np.shape(p)[1]
        q = p[:, 0].reshape(2, 1)
        for k in range(N - 2):
            v1 = p[:, k + 1] - p[:, k]
            v2 = p[:, k + 2] - p[:, k + 1]
            v1n = v1 / np.sqrt(np.dot(v1, v1))
            v2n = v2 / np.sqrt(np.dot(v2, v2))
            u = np.dot(v1n, v2n)
            if u < 0.98:
                q = np.hstack((q, p[:, k + 1].reshape(2, 1)))
        vini = q[:, 1] - q[:, 0]
        vin = vini / np.sqrt(np.dot(vini, vini))
        v = np.dot(v2n, vin)
        if v > 0.98:
            q = q[:, 1:]
        y = q.T.copy()
        ls = shg.asLineString(y)
        poly2 = shg.Polygon(ls)
        return(poly2)


    def buildGvc(self, **kwargs):
        """ Create visibility graph for a convex polygon

        Parameters
        ----------

        display   : boolean
            default : False
        fig       : matplotlib.figure.pyplot
        ax        : axes
        udeg2     : np.array indexes of points of degree 2
            default = []
        eded   : boolean
            default True
        indoor : boolean
            default True

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.geomutil import *
            >>> import shapely.geometry as shg
            >>> import matplotlib.pyplot as plt
            >>> points  = shg.MultiPoint([(0, 0), (0, 1), (2.5,1), (2.5, 2), \
                                          (2.8,2), (2.8, 1.1), (3.2, 1.1), \
                                          (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
            >>> polyg   = Polygon(points)
            >>> Gv      = polyg.buildGv(show=True)
            >>> plt.axis('off')
            (-0.5, 4.0, -0.5, 2.5)
            >>> title = plt.title('Testing buildGv')


        Notes
        -----

        Segment k and (k+1)%N share segment (k+1)%N
        The degree of a point is dependent from other polygons around

        Topological error can be raised if the point coordinates accuracy
        is not limited.

        Nodes of polygon are numbered in the global graph in vnodes member.

        See Also
        --------

        pylayers.gis.layout.Layout.buildGv

        """

        defaults = {'udeg2': np.array([]),
                    'eded':True,
                    'open':True,
                    'indoor':True
                    }

##       initialize function attributes

        for key, value in defaults.items():
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, value)
                kwargs[key] = value

        Gv = nx.Graph()
        Gv.pos = {}
        if kwargs['open']:
            pass
        else:
            pass

        lring = self.exterior
        #
        # Calculate interior normals
        #
        x, y = lring.xy
        p = np.array([x[0:-1], y[0:-1]])
        #
        # determine convex points
        #
        #pdb.set_trace()
        tcc, n = self.ptconvex()
        # Np = self.Np
        Np = np.shape(self.exterior.xy)[1] -1
        #
        # retrieve
        #  npt points label sequence
        #  nseg segments label sequence
        #
        # vnodes do not necessarily start with a point
        #

        npt  = filter(lambda x : x<0,self.vnodes)
        nseg = filter(lambda x : x>0,self.vnodes)

        #
        # in convex case all segments see all segments
        #
        for nk in combinations(nseg, 2):
            Gv.add_edge(nk[0],nk[1],weight=0)

        #
        # Update position of points in Gv
        #
        for nk in Gv.node:
            Gv.pos[nk] = (p[0, nk], p[1, nk])

        xr, yr = lring.xy

        #
        # Determine diffraction points
        #
        # deg2 : if null:
        #           the point is kept
        #        if convex:
        #           the point is kept
        #        else:
        #           the point is not kept
        #
        if indoor:
            uconvex = np.nonzero(tcc == 1)[0] # convex point position
        else:
            uconvex = np.nonzero(tcc == -1)[0] # convex point position
        uzero = np.nonzero(tcc == 0)[0]   # planar point (joining two parallel segment)
        udiffdoor = np.intersect1d(uzero, udeg2)  # degree 2 paralell points are often doors and windows
        udiff = np.hstack((uconvex, udiffdoor)).astype('int') # diffracting point
        #
        # 1) Calculate node-node visibility
        #
        #
        # Between all combinations of diffracting points
        # create a segment and check whether it is fully included in the
        # polygon.
        # If verified then there is a visibility between the 2 points.
        #
        for nk in combinations(udiff, 2):
            p1 = p[:, nk[0]]
            p2 = p[:, nk[1]]
            seg = shg.LineString(((p1[0], p1[1]), (p2[0], p2[1])))
            if self.contains(seg):
                Gv.add_edge(npt[nk[0]], npt[nk[1]], weight=0)

        #
        #  2) Calculate edge-edge and node-edge visibility
        #

        for nk in range(Np):   # loop on range of number of points
            ptk = p[:, nk]     # tail point
            phk = p[:, (nk + 1) % Np] # head point (%Np to get 0 as last point)

            # lnk : unitary vector on segment nk
            lk = phk - ptk
            nlk = np.sqrt(np.dot(lk, lk))
            lnk = lk / nlk

            # the epsilon is (1/1000) of the segment length
            epsilonk = nlk / 1000.  # this can be dangerous (epsilon can be large)

            # x--o----------------------o--x
            #    +eps                  -eps
            pcornert = ptk + lnk * epsilonk  # + n[:,nk]*epsilon
            pcornerh = phk - lnk * epsilonk  # + n[:,nk]*epsilon

        #
        # in any case no ray towark nk
        # if nk is convex no ray toward (nk-1)%Np
        #
        # start from the two extremity of the segment
            for i, pcorner in enumerate([pcornert, pcornerh]):
                #
                #  if tail point
                #           remove nk segment
                #  and if the point is convex
                #          remove previous segment
                #
                #  si point head
                #
                listpoint = range(Np)
                listpoint.remove(nk)   # remove current point
                if i == 0:  # first iteration pcornert
                    if nk in uconvex:  # == 1
                        listpoint.remove((nk - 1) % Np)
                if i == 1:  # second iteration pcornerh
                    if (nk + 1) % Np in uconvex: # ==1
                        listpoint.remove((nk + 1) % Np)

                for ns in listpoint:
                    pts = p[:, ns]
                    phs = p[:, (ns + 1) % Np]
                    # Add B.Uguen 2/01/2014 no possible visibility relation between aligned segments
                    if (not (isaligned(pts,phs,ptk) & isaligned(pts,phs,phk))):
                        ls = phs - pts
                        nls = np.sqrt(np.dot(ls, ls))
                        lns = ls / nls
                        epsilons = nls / 1000.
                        pte = pts + lns * epsilons  # + n[:,ns]*epsilon
                        phe = phs - lns * epsilons  # + n[:,ns]*epsilon
                        tbr = pyu.bitreverse(16, 5) / 16.
                        for alpha in tbr:
                            pa = pte + alpha * (phe - pte)
                            seg = shg.LineString((pcorner, pa))
                            #print "seg: ",seg.xy
                            #if npt[nk] == -3:
                            #    plt.plot(np.array([pcorner[0],pa[0]]),np.array([pcorner[1],pa[1]]),linewidth=0.2,color='k')
                            #    plt.draw()
                            # topological error can be raised here
                            seg2 = self.intersection(seg)
                            #if self.contains(seg):
                            if seg2.almost_equals(seg, decimal=4):
                                #print alpha,nk,ns
                                #plt.plot(np.array([pcorner[0],pa[0]]),np.array([pcorner[1],pa[1]]),linewidth=2,color='r')
                                #Gv.add_edge(-(uconvex[nk]+1),ns+1,weight=10)
                                if i == 0:
                                    if nk in udiff:
                                        Gv.add_edge(npt[nk], nseg[ns], weight=1)
                                        #plt.plot(np.array([Gv.pos[npt[nk]][0],Gv.pos[nseg[ns]][0]]),np.array([Gv.pos[npt[nk]][1],Gv.pos[nseg[ns]][1]]),'r')
                                if i == 1:
                                    if (nk + 1) % Np in udiff:
                                        Gv.add_edge(npt[(nk + 1) % Np], nseg[ns], weight=1)
                                        #plt.plot(np.array([Gv.pos[npt[(nk+1)%Np]][0],Gv.pos[nseg[ns]][0]]),np.array([Gv.pos[npt[(nk+1)%Np]][1],Gv.pos[nseg[ns]][1]]),'g')
                                    #plt.draw()
                                #if i==1:
                                #if (((nseg[nk]==10) & (nseg[ns]==7)) or
                                #    ((nseg[nk]==7) & (nseg[ns]==10))):
                                #    pdb.set_trace()
                                if nseg[nk] != nseg[ns]:
                                    if kwargs['eded']:
                                        Gv.add_edge(nseg[nk], nseg[ns], weight=1)
                                    #else:
                                    #    print nseg[nk],nseg[ns]
                                    #    print pts,phs
                                    #    print ptk,phk
                                    #if (((nseg[nk]==10) & (nseg[ns]==7)) or
                                    #    ((nseg[nk]==7) & (nseg[ns]==10))):
                                    #    plt.plot(np.array([Gv.pos[nseg[nk]][0],Gv.pos[nseg[ns]][0]]),np.array([Gv.pos[nseg[nk]][1],Gv.pos[nseg[ns]][1]]),'b')
                                    #    plt.plot(np.array([pcorner[0],pa[0]]),np.array([pcorner[1],pa[1]]),'b')
                                    #    print "seg: ",seg.xy
                                    #    print "seg2: ",seg2.xy
                                    #    print nseg[nk],nseg[ns]
                                    #    print pcorner , ptk
                                    #    print  alpha , pa ,pte
                                    #    plt.draw()
                                    #    raw_input()
                                break
                    #else:
                        #print p
                        #print ns
                        #print nk
                        #print 'nsegnk : ',nseg[nk]
                        #print 'nsegns', nseg[ns]
                        #print 'ptk : ',ptk
                        #print 'phk : ',phk
                        #print 'pts : ',pts
                        #print 'phs : ',phs
                        #print "aligne :",nseg[nk],nseg[ns]
                        #pdb.set_trace()

        if kwargs['show']:
            nodes = np.array(Gv.nodes())
            uneg = list(nodes[np.nonzero(nodes < 0)[0]])
            upos = list(nodes[np.nonzero(nodes > 0)[0]])
            nx.draw_networkx_nodes(Gv, Gv.pos, nodelist=upos,
                                   node_color='blue', node_size=300, alpha=0.3)
            nx.draw_networkx_nodes(Gv, Gv.pos, nodelist=uneg,
                                   node_color='red', node_size=300, alpha=0.3)
            nx.draw_networkx_labels(Gv, Gv.pos)

            ndnd, nded, eded = gru.edgetype(Gv)

            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=eded,
                                   edge_color='blue', width=2)
            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=ndnd,
                                   edge_color='red', width=2)
            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=nded,
                                   edge_color='green', width=2)

            #label = {}
            #for (u,v) in Gv.edges():
            #    d = Gv.get_edge_data(u,v)
            #    label[(u,v)]=d['weight']

            #edge_label=nx.draw_networkx_edge_labels(Gv,Gv.pos,edge_labels=label)

        return(Gv)
    def buildGv(self, **kwargs):
        """ Create  visibility graph for a polygon

        Parameters
        ----------

        display   : boolean
            default : False
        fig       : matplotlib.figure.pyplot
        ax        : axes
        udeg2     : np.array indexes of points of degree 2
            default = []
        eded   : boolean
            default True
        indoor : boolean
            default True

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.geomutil import *
            >>> import shapely.geometry as shg
            >>> import matplotlib.pyplot as plt
            >>> points  = shg.MultiPoint([(0, 0), (0, 1), (2.5,1), (2.5, 2), \
                                          (2.8,2), (2.8, 1.1), (3.2, 1.1), \
                                          (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
            >>> polyg   = Polygon(points)
            >>> Gv      = polyg.buildGv(show=True)
            >>> plt.axis('off')
            (-0.5, 4.0, -0.5, 2.5)
            >>> title = plt.title('Testing buildGv')


        Notes
        -----

        Segment k and (k+1)%N share segment (k+1)%N
        The degree of a point is dependent from other polygons around

        Topological error can be raised if the point coordinates accuracy
        is not limited.

        Nodes of polygon are numbered in the global graph in vnodes member.

        See Also
        --------

        pylayers.gis.layout.Layout.buildGv

        """

        defaults = {'show': False,
                    'fig': [],
                    'ax': [],
                    'udeg2': np.array([]),
                    'eded':True,
                    'indoor':True
                    }

##       initialize function attributes

        for key, value in defaults.items():
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, value)
                kwargs[key] = value

        #self.args=args
        if kwargs['show']:
            if kwargs['fig'] == []:
                fig = plt.figure(figsize=(20,20))
                fig.set_frameon(True)
            else:
                fig = kwargs['fig']

            if kwargs['ax'] == []:
                ax = fig.gca()
            else:
                ax = kwargs['ax']
            plt.ion()

        udeg2 = kwargs['udeg2']

        GRAY = '#999999'
        Gv = nx.Graph()
        Gv.pos = {}

        #pdb.set_trace()
        lring = self.exterior
        #
        # Calculate interior normals
        #
        x, y = lring.xy
        p = np.array([x[0:-1], y[0:-1]])
        #
        # determine convex points
        #
        #pdb.set_trace()
        tcc, n = self.ptconvex()
        # Np = self.Np
        Np = np.shape(self.exterior.xy)[1] -1
        #
        # retrieve
        #  npt points label sequence
        #  nseg segments label sequence
        #
        # vnodes do not necessarily start with a point
        #
        if self.vnodes[0] < 0:
            ipt  = 2 * np.arange(Np)
            iseg = 2 * np.arange(Np) + 1
        else:
            ipt = 2 * np.arange(Np) + 1
            iseg = 2 * np.arange(Np)

        npt = self.vnodes[ipt]
        nseg = self.vnodes[iseg]
        #print "npt : ",npt
        #print "nseg : ",nseg

        assert  np.all(npt < 0), "something wrong with points"
        assert  np.all(nseg > 0), "something wrong with segments"
        #
        #
        # Create middle point on lring
        #
        # Warning lring recopy the node at the end of the sequence
        #
        # A problem arises from the fact that a vnodes sequence
        # do no necessarily starts with a point (negative node)
        #
        #
        tpm = []
        for ik, k in enumerate(lring.coords):
            pt = np.array(k)
            try:
                pm = (pt + pm1) / 2.
                if self.vnodes[0] < 0:
                    Gv.pos[nseg[ik - 1]] = (pm[0], pm[1])
                else:
                    Gv.pos[nseg[ik % Np]] = (pm[0], pm[1])
                tpm.append(pm)
                pm1 = pt
            except:
                pm1 = pt
        #
        # Update position of points in Gv
        #
        for nk in range(Np):
            #nnode = -(nk+1)
            Gv.pos[npt[nk]] = (p[0, nk], p[1, nk])

        xr, yr = lring.xy

        #
        # Determine diffraction points
        #
        # deg2 : if null:
        #           the point is kept
        #        if convex:
        #           the point is kept
        #        else:
        #           the point is not kept
        #
        if kwargs['indoor']:
            uconvex = np.nonzero(tcc == 1)[0] # convex point position
        else:
            uconvex = np.nonzero(tcc == -1)[0] # convex point position
        uzero = np.nonzero(tcc == 0)[0]   # planar point (joining two parallel segment)
        udiffdoor = np.intersect1d(uzero, udeg2)  # degree 2 paralell points are often doors and windows
        udiff = np.hstack((uconvex, udiffdoor)).astype('int') # diffracting point
        #print "vnodes",self.vnodes
        #print "tcc : ",tcc
        #print "uzero : ",uzero
        #print "udiffdoor : ",udiffdoor
        #print "udiff",udiff
        #print "udeg2",udeg2
        #print "npt",npt
        #if udiff!=[]:
        #    print "diff : ",npt[udiff]
        #if udeg2!=[]:
        #    print "deg2 : ",npt[udeg2]
        #if uzero!=[]:
        #    print "zero :",npt[uzero]
        #
        # if show == True display points and polygon
        #
        if kwargs['show']:
            points1 = shg.MultiPoint(lring)
            for k, pt in enumerate(points1):
                if k in uconvex:
                    ax.plot(pt.x, pt.y, 'o', color='red')
                elif k in udiffdoor:
                    ax.plot(pt.x, pt.y, 'o', color='blue')
                else:
                    ax.plot(pt.x, pt.y, 'o', color=GRAY)

            patch = PolygonPatch(self, facecolor='#6699cc',
                                 edgecolor='#000000', alpha=0.5, zorder=2)
            ax.add_patch(patch)
        #pdb.set_trace()
        #
        #  1) Calculate node-node visibility
        #
        # The algorithm exploits definition of convexity.
        #
        # Between all combinations of diffracting points
        # create a segment and check whether it is fully included in the
        # polygon.
        # If verified then there is a visibility between the 2 points.
        #
        for nk in combinations(udiff, 2):
            p1 = p[:, nk[0]]
            p2 = p[:, nk[1]]
            seg = shg.LineString(((p1[0], p1[1]), (p2[0], p2[1])))
            if self.contains(seg):
                Gv.add_edge(npt[nk[0]], npt[nk[1]], weight=0)

        #
        #  2) Calculate edge-edge and node-edge visibility
        #

        for nk in range(Np):   # loop on range of number of points
            ptk = p[:, nk]     # tail point
            phk = p[:, (nk + 1) % Np] # head point (%Np to get 0 as last point)

            # lnk : unitary vector on segment nk
            lk = phk - ptk
            nlk = np.sqrt(np.dot(lk, lk))
            lnk = lk / nlk

            # the epsilon is (1/1000) of the segment length
            epsilonk = nlk / 1000.  # this can be dangerous (epsilon can be large)

            # x--o----------------------o--x
            #    +eps                  -eps
            pcornert = ptk + lnk * epsilonk  # + n[:,nk]*epsilon
            pcornerh = phk - lnk * epsilonk  # + n[:,nk]*epsilon

        #
        # in any case no ray towark nk
        # if nk is convex no ray toward (nk-1)%Np
        #
        # start from the two extremity of the segment
            for i, pcorner in enumerate([pcornert, pcornerh]):
                #
                #  if tail point
                #           remove nk segment
                #  and if the point is convex
                #          remove previous segment
                #
                #  si point head
                #
                listpoint = range(Np)
                listpoint.remove(nk)   # remove current point
                if i == 0:  # first iteration pcornert
                    if nk in uconvex:  # == 1
                        listpoint.remove((nk - 1) % Np)
                if i == 1:  # second iteration pcornerh
                    if (nk + 1) % Np in uconvex: # ==1
                        listpoint.remove((nk + 1) % Np)

                for ns in listpoint:
                    pts = p[:, ns]
                    phs = p[:, (ns + 1) % Np]
                    # Add B.Uguen 2/01/2014 no possible visibility relation between aligned segments
                    if (not (isaligned(pts,phs,ptk) & isaligned(pts,phs,phk))):
                        ls = phs - pts
                        nls = np.sqrt(np.dot(ls, ls))
                        lns = ls / nls
                        epsilons = nls / 1000.
                        pte = pts + lns * epsilons  # + n[:,ns]*epsilon
                        phe = phs - lns * epsilons  # + n[:,ns]*epsilon
                        tbr = pyu.bitreverse(16, 5) / 16.
                        for alpha in tbr:
                            pa = pte + alpha * (phe - pte)
                            seg = shg.LineString((pcorner, pa))
                            #print "seg: ",seg.xy
                            #if npt[nk] == -3:
                            #    plt.plot(np.array([pcorner[0],pa[0]]),np.array([pcorner[1],pa[1]]),linewidth=0.2,color='k')
                            #    plt.draw()
                            # topological error can be raised here
                            seg2 = self.intersection(seg)
                            #if self.contains(seg):
                            if seg2.almost_equals(seg, decimal=4):
                                #print alpha,nk,ns
                                #plt.plot(np.array([pcorner[0],pa[0]]),np.array([pcorner[1],pa[1]]),linewidth=2,color='r')
                                #Gv.add_edge(-(uconvex[nk]+1),ns+1,weight=10)
                                if i == 0:
                                    if nk in udiff:
                                        Gv.add_edge(npt[nk], nseg[ns], weight=1)
                                        #plt.plot(np.array([Gv.pos[npt[nk]][0],Gv.pos[nseg[ns]][0]]),np.array([Gv.pos[npt[nk]][1],Gv.pos[nseg[ns]][1]]),'r')
                                if i == 1:
                                    if (nk + 1) % Np in udiff:
                                        Gv.add_edge(npt[(nk + 1) % Np], nseg[ns], weight=1)
                                        #plt.plot(np.array([Gv.pos[npt[(nk+1)%Np]][0],Gv.pos[nseg[ns]][0]]),np.array([Gv.pos[npt[(nk+1)%Np]][1],Gv.pos[nseg[ns]][1]]),'g')
                                    #plt.draw()
                                #if i==1:
                                #if (((nseg[nk]==10) & (nseg[ns]==7)) or
                                #    ((nseg[nk]==7) & (nseg[ns]==10))):
                                #    pdb.set_trace()
                                if nseg[nk] != nseg[ns]:
                                    if kwargs['eded']:
                                        Gv.add_edge(nseg[nk], nseg[ns], weight=1)
                                    #else:
                                    #    print nseg[nk],nseg[ns]
                                    #    print pts,phs
                                    #    print ptk,phk
                                    #if (((nseg[nk]==10) & (nseg[ns]==7)) or
                                    #    ((nseg[nk]==7) & (nseg[ns]==10))):
                                    #    plt.plot(np.array([Gv.pos[nseg[nk]][0],Gv.pos[nseg[ns]][0]]),np.array([Gv.pos[nseg[nk]][1],Gv.pos[nseg[ns]][1]]),'b')
                                    #    plt.plot(np.array([pcorner[0],pa[0]]),np.array([pcorner[1],pa[1]]),'b')
                                    #    print "seg: ",seg.xy
                                    #    print "seg2: ",seg2.xy
                                    #    print nseg[nk],nseg[ns]
                                    #    print pcorner , ptk
                                    #    print  alpha , pa ,pte
                                    #    plt.draw()
                                    #    raw_input()
                                break
                    #else:
                        #print p
                        #print ns
                        #print nk
                        #print 'nsegnk : ',nseg[nk]
                        #print 'nsegns', nseg[ns]
                        #print 'ptk : ',ptk
                        #print 'phk : ',phk
                        #print 'pts : ',pts
                        #print 'phs : ',phs
                        #print "aligne :",nseg[nk],nseg[ns]
                        #pdb.set_trace()

        if kwargs['show']:
            nodes = np.array(Gv.nodes())
            uneg = list(nodes[np.nonzero(nodes < 0)[0]])
            upos = list(nodes[np.nonzero(nodes > 0)[0]])
            nx.draw_networkx_nodes(Gv, Gv.pos, nodelist=upos,
                                   node_color='blue', node_size=300, alpha=0.3)
            nx.draw_networkx_nodes(Gv, Gv.pos, nodelist=uneg,
                                   node_color='red', node_size=300, alpha=0.3)
            nx.draw_networkx_labels(Gv, Gv.pos)

            ndnd, nded, eded = gru.edgetype(Gv)

            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=eded,
                                   edge_color='blue', width=2)
            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=ndnd,
                                   edge_color='red', width=2)
            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=nded,
                                   edge_color='green', width=2)

            #label = {}
            #for (u,v) in Gv.edges():
            #    d = Gv.get_edge_data(u,v)
            #    label[(u,v)]=d['weight']

            #edge_label=nx.draw_networkx_edge_labels(Gv,Gv.pos,edge_labels=label)

        return(Gv)
    def showGv(self, **kwargs):
        """ show graph Gv

        Parameters
        ----------
        display
        fig
        ax
        ndnd    : boolean
            display node/node
        nded    : boolean
            display node/edge
        eded    : boolean
            display edge/edge
        linewidth: float
            default 2

        """
        defaults = {'display': False,
                    'fig': [],
                    'ax': [],
                    'ndnd': True,
                    'nded': False,
                    'ndnd': False,
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

        lring = self.exterior
        points = shg.MultiPoint(lring)
        for k, pt in enumerate(points):
            if tcc[k % Np] == 1:
                ax.plot(pt.x, pt.y, 'o', color='red')
            else:
                ax.plot(pt.x, pt.y, 'o', color=GRAY)
            k = k + 1

        patch = PolygonPatch(self, facecolor='#6699cc',
                             edgecolor='#6699cc', alpha=0.5, zorder=2)
        ax.add_patch(patch)
        nodes = np.array(Gv.nodes())

        uneg = list(nodes[np.nonzero(nodes < 0)[0]])
        upos = list(nodes[np.nonzero(nodes > 0)[0]])
        if kwargs['nodes']:
            nx.draw_networkx_nodes(Gv, Gv.pos, nodelist=upos,
                                   node_color='blue', node_size=300, alpha=0.3)
            nx.draw_networkx_nodes(Gv, Gv.pos, nodelist=uneg,
                                   node_color='red', node_size=300, alpha=0.3)
            nx.draw_networkx_labels(Gv, Gv.pos)

            ndnd, nded, eded = gru.edgetype(Gv)

        if kwargs['eded']:
            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=eded,
                                   edge_color='blue', width=2)
        if kwargs['ndnd']:
            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=ndnd,
                                   edge_color='red', width=2)
        if kwargs['nded']:
            nx.draw_networkx_edges(Gv, Gv.pos, edgelist=nded,
                                   edge_color='green', width=2)

        return(fig, ax)

    def ptconvex2(self):
        """ Determine convex / concave points in the Polygon

            !!! Warning !!! cvex and ccve can be switched
            depends on the Polygon direction of travel
        
        Returns
        -------
        cvex : list of convex points
        ccve : list of concave points

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.geomutil import *
            >>> import shapely.geometry as shg
            >>> import matplotlib.pyplot as plt
            >>> points  = shg.MultiPoint([(0, 0), (0, 1), (3.2, 1), (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
            >>> polyg1   = Polygon(points)
            >>> cvex,ccave   = polyg.ptconvex2() 
            >>> points  = shg.MultiPoint([(0, 0), (0, 1), (-3.2, 1), (-3.2, 0.7), (-0.4, 0.7), (-0.4, 0)])
            >>> polyg1   = Polygon(points)
            >>> cvex,ccave   = polyg.ptconvex2() 

        """
        if not hasattr(self,'xy'):
            self.coorddeter()
        
        pts = filter(lambda x: x<0,self.vnodes)
        A=self.xy[:,:-1]
        B=np.roll(A,-1)
        C=np.roll(B,-1)
        if self.signedarea()>0:
            cw = ccw(C,B,A)
        else :
            cw = ccw(A,B,C)
        cvex = np.array(pts)[np.roll(cw,+1)]
        ccve = np.array(pts)[np.roll(~cw,+1)]
        
        return cvex.tolist(),ccve.tolist()

    def ptconvex(self, display=False):
        """ Return a list of booleans indicating points convexity

        Parameters
        ----------

        display : boolean
            default False


        Returns
        -------

        tcc     : np.array (1x Nseg)
            1 if convex , -1 if concav , 0 if plane
        n       :  array(2xNseg)
            segments normals

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.geomutil import *
            >>> import shapely.geometry as shg
            >>> import matplotlib.pyplot as plt
            >>> points  = shg.MultiPoint([(0, 0), (0, 1), (3.2, 1), (3.2, 0.7), (0.4, 0.7), (0.4, 0)])
            >>> N = len(points)
            >>> polyg   = Polygon(points)
            >>> tcc,n   = polyg.ptconvex()
            >>> #k = 0
            >>> #for p in points:
            >>> #  if tcc[k] == 1 :
            >>> #      plt.plot(p.x, p.y, 'o', color='red',alpha=1)
            >>> #  else:
            >>> #      plt.plot(p.x, p.y, 'o', color='blue',alpha=0.3)
            >>> #  k = k+1
            >>> #polyg.plot()
            >>> #plt.figure()
            >>> #points  = shg.MultiPoint([(0, 0), (1, 1), (2, 0), (1, 0)])
            >>> #poly    = Polygon(points)
            >>> #tcc,n   = polyg.ptconvex()
            >>> #poly.plot()


        Notes
        ------

            This function determines the convex and concav points of a polygon.
            As there is no orientation convention for the polygon the sign of the cross
            product can't be directly interpreted. So we exploit the following
            property :

            Let N be the number of points of the Polygon. N  = Nx + Nc where
            Nx is the number of convex points and Nc the number of concav points

            We have Nx >= Nc

            If a point is common to two parallel segments, the cross product is = 0

        See Also
        --------

        Lr2n

        """

        lring = self.exterior

        #
        # Calculate interior normals
        #

        x, y = lring.xy
        Np = len(x) - 1
        Nseg = Np
        p = np.array([x[0:-1], y[0:-1]])

        n = Lr2n(p)

        tcc = np.zeros(Np)

        #
        # cross product between two adjascent normals
        #
        for k in range(Nseg):
            nk = n[:, (k - 1) % Nseg]
            nkp1 = n[:, k]
            v = np.cross(nk, nkp1)
            tcc[k] = v

        #
        # warning this test is fragile
        #
        # debug : print tcc
        #
        # The purpose here is to remove flat transition
        #
        upos = np.nonzero(tcc > 1e-2)[0]
        uneg = np.nonzero(tcc < -1e-2)[0]

        if len(upos) > len(uneg):
            nconvex = uneg
            nconcav = upos

        if len(upos) < len(uneg):
            nconvex = upos
            nconcav = uneg

        if len(upos) == len(uneg):
            logging.warning("polygon is a star")
            # self.plot()
            # pdb.set_trace()


        tcc = np.zeros(Np)
        tcc[nconvex] = 1
        tcc[nconcav] = -1
        #print "ptseg tcc ",tcc
        upos = np.nonzero(tcc > 1e-4)[0]
        return(tcc, n)

class Geomview(PyLayers):
    """ Geomview file class

    This class is parent of  GeomVect Geomlist Geomoff

    Methods
    -------

    show3

    """
    def __init__(self, _filename,clear=False):
        filename = pyu.getlong(_filename, "geom")
        self.filename = filename
        if clear:
            fd = open(self.filename,'w')
            fd.close()

    def show3(self):
        """
        .. todo:
             change  background
             look for other geomview options
        """
        chaine = "geomview  -b 1 1 1 " + self.filename + " 2>/dev/null &"
        os.system(chaine)


class Geomlist(Geomview):
    """

    """
    def __init__(self, _filename,clear=False):
        _filename = _filename + '.list'
        Geomview.__init__(self, _filename,clear=clear)

    def append(self, strg):
        """
           append a line in .list file
        """
        fd = open(self.filename, 'a')
        fd.write(strg)
        fd.close()


class GeomVect(Geomview):
    """ Geomview VECT file class

       + NPolylines  NVertices  NColors
       + Nv[0] ... Nv[NPolylines-1]     # number of vertices in each polyline
       + Nc[0] ... Nc[NPolylines-1]     # number of colors supplied in each polyline
       + Vert[0] ... Vert[NVertices-1]  # All the vertices (3*NVertices floats)
       + Color[0] ... Color[NColors-1]  # All the colors   (4*NColors floats, RGBA)

        VECT objects represent lists of polylines (strings of connected line segments, possibly closed).

       A degenerate polyline can be used to represent a point:
       A VECT file begins with the key word VECT or 4VECT and three integers:

       NLines, NVertices, and NColors.
           Here NLines is the number of polylines in the file,
           NVertices the total number of vertices, and NColors the number of
           colors as explained below.
           Next come NLines 16-bit integers
       Nv[0] Nv[1] Nv[2] ... Nv[NLines-1]
           giving the number of vertices in each polyline.
           A negative number indicates a closed polyline; 1 denotes a single-pixel point.
           The sum (of absolute values) of the Nv[i] must equal NVertices.
           Next come NLines more 16-bit integers
           Nc[i]: the number of colors in each polyline.
           Normally one of three values:
                0 : No color is specified for this polyline.
                    It's drawn in the same color as the previous polyline.
                1 : A single color is specified.
                    The entire polyline is drawn in that color.
                abs(Nv[i]) : Each vertex has a color.
                    Either each segment is drawn in the corresponding color,
                    or the colors are smoothly interpolated along the line segments,
                    depending on the implementation.

            Next come NVertices groups of 3 or 4 floating-point numbers:
                            the coordinates of all the vertices.
            If the keyword is 4VECT then there are 4 values per vertex.
            The first abs(Nv[0]) of them form the first polyline,
                    the next abs(Nv[1]) form the second and so on.
            Finally NColors groups of 4 floating-point numbers give red,
                    green, blue and alpha (opacity) values.
            The first Nc[0] of them apply to the first polyline, and so on.

    Methods
    -------

    geomBase
        display a frame
    ellipse
        display an ellipse
    points
        display a set of points

    """
    def __init__(self, _filename='geomdef',clear=False):
        _filename = _filename + '.vect'
        Geomview.__init__(self, _filename,clear=clear)

    def segments(self, ds, i2d=True, linewidth=2):
        """ display segments

        Parameters
        ----------

        ds   : dictionnary
            len ds
        i2d  : boolean (defaut True)
            2d indicator
        linewidth : float
            default 2

        """
        fo = open(self.filename, "w")
        fo.write("appearance { linewidth %d }\n" % linewidth)
        fo.write("VECT\n")
        Ns = len(ds)
        fo.write("%d %d %d\n" % (Ns, 2 * Ns, 0))
            # 3 Lines 6 Vertices 3 colors
        for k in range(Ns):
            fo.write("2 ")
        fo.write("\n")
        for k in range(Ns):
            fo.write("0 ")
        fo.write("\n")
        for k in ds:
            (pta, phe) = ds[k]
            if i2d:
                fo.write("%6.3f %6.3f %6.3f\n" % (pta[0], pta[1], 0.0))
                fo.write("%6.3f %6.3f %6.3f\n" % (phe[0], phe[1], 0.0))
            else:
                fo.write("%6.3f %6.3f %6.3f\n" % (pta[0], pta[1], pta[2]))
                fo.write("%6.3f %6.3f %6.3f\n" % (phe[0], phe[1], phe[2]))
        fo.close()

    def geomBase(self, M, pt=np.array([0., 0., 0.]), col=np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]]),
                 linewidth=3,scale=1):
        """ Construct a geomview vect file for vizualisation of a frame

        Notes
        -----

        by default the geomview filename is base0.vect


        Parameters
        ----------
        M     : ndarray (3 x 3 )
            [ v1, v2, v3 ]
        pt    : np.array
            origin point  (default (0,0,0))
        col   :
            color (3x3)
        linewidth :
            linewidth (default 3)

        Examples
        --------

        >>> from pylayers.util.geomutil import *
        >>> import numpy as np
        >>> v1 = np.array([1,0,0])
        >>> v2 = np.array([0,1,0])
        >>> v3 = np.array([0,0,1])
        >>> M  = np.vstack((v1,v2,v3))
        >>> #gv = GeomVect('test')
        >>> #gv.geomBase(M)
        >>> #gv.show3()

        """
        M = M * scale
        fo = open(self.filename, "w")
        fo.write("appearance { linewidth %d }\n" % linewidth)
        fo.write("VECT\n")
        fo.write("3 6 3\n")   # 3 Lines 6 Vertices 3 colors
        fo.write("2 2 2\n")   # 2 points per lines
        fo.write("1 1 1\n")   # 1 color per line
        fo.write("%6.3f %6.3f %6.3f\n" % (pt[0], pt[1], pt[2]))
        fo.write("%6.3f %6.3f %6.3f\n" % (pt[0] + M[0, 0], pt[1] +
                                          M[1, 0], pt[2] + M[2, 0]))
        fo.write("%6.3f %6.3f %6.3f\n" % (pt[0], pt[1], pt[2]))
        fo.write("%6.3f %6.3f %6.3f\n" % (pt[0] + M[0, 1], pt[1] +
                                          M[1, 1], pt[2] + M[2, 1]))
        fo.write("%6.3f %6.3f %6.3f\n" % (pt[0], pt[1], pt[2]))
        fo.write("%6.3f %6.3f %6.3f\n" % (pt[0] + M[0, 2], pt[1] +
                                          M[1, 2], pt[2] + M[2, 2]))
        fo.write("%6.3f %6.3f %6.3f  0.\n" % (col[0, 0], col[0, 1], col[0, 2]))
        fo.write("%6.3f %6.3f %6.3f  0.\n" % (col[1, 0], col[1, 1], col[1, 2]))
        fo.write("%6.3f %6.3f %6.3f  0.\n" % (col[2, 0], col[2, 1], col[2, 2]))
        #fo.write("{<}\n")
        fo.close()



    def points(self, pt, colorname='blue'):
        """ Geomview display a set of points with color

        Parameters
        ----------
        pt
           sequence of points np.ndarray or dictionnary whose value is a tuple (x,y,z)
        colorname
            a colorname from coldict keys

        Examples
        --------

        >>> import numpy as np
        >>> from pylayers.util.geomutil import *
        >>> import scipy as sp
        >>> pt1 = sp.rand(3,10)
        >>> pt2 = { 1:(0,0,0),2:(10,10,10),3:(0,10,0),4:(10,0,0)}
        >>> gv1 = GeomVect('test1')
        >>> gv1.points(pt1)
        >>> #gv1.show3()
        >>> gv2 = GeomVect('test2')
        >>> gv2.points(pt2)
        >>> #gv2.show3()

        .. todo::
            colorbar depending of a value associated with point
        """
        fo = open(self.filename, "w")
        if type(pt) == list:
            pt = np.array(pt).reshape(3,1)
        if type(pt) == dict:
            npt = len(pt.keys())
        if type(pt) == np.ndarray:
            npt = np.shape(pt)[1]
        snpt = str(npt) + "\n"
        snpt2 = str(npt) + " " + str(npt) + " " + str(npt) + "\n"
        if npt > 1:
            fo.write("appearance{\n")
            fo.write("linewidth 8}\n")
            fo.write("VECT\n")
            fo.write(snpt2)
            fo.write("1 " * npt + "\n")
            fo.write("1 " * npt + "\n")
        else:
            fo.write("ESPHERE\n")
            fo.write("0.2\n")

        if type(pt) == dict:
            for i in range(npt):
                x = str(pt[pt.keys()[i]][0])
                y = str(pt[pt.keys()[i]][1])
                try:
                    z = str(pt[pt.keys()[i]][2])
                except:
                    z = str(0.0)
                chaine = x + " " + y + " " + z + "\n"
                fo.write(chaine)
        if type(pt) == np.ndarray:
            for i in range(npt):
                x = str(pt[0, i]).replace(',', '.')
                y = str(pt[1, i]).replace(',', '.')
                try:
                    z = str(pt[2, i]).replace(',', '.')
                except:
                    z = str(0.0)
                chaine = x + " " + y + " " + z + "\n"
                fo.write(chaine)

        coldic = pyu.coldict()
        col = pyu.rgb(coldic[colorname], 'float')
        if npt > 1:
            for i in range(npt):
                fo.write("%6.3f %6.3f %6.3f 1\n" % (col[0], col[1], col[2]))
        fo.close()


class Geomoff(Geomview):
    """

    Notes
    -----

    Class Geomview OFF File (Object File Format)
    [ST][C][N][4][n]OFF  #header keyword
    [Ndim] # spac dimension of vertices, present only if nOFF
    NVertices NFaces NEdges

    x[0],y[0] z[0]

    # Vertices,possibly with normals
    #colors, and/or texture coordinates, in that order, if the
    # prefixes N , C , ST are present
    # If 4OFF , each vertex has 4 components
    # including a final homogeneous component
    # If nOFF, each vertex has Ndim components
    # If 4nOFF , each vertex has Ndim+1 components
    ....
    x[NVertices-1],y[NVertices-1],z[NVertices-1]

    # Faces
    # Nv = # vertices on this face
    # v[0] ... v[Nv-1] : vertex indices
    # in range 0... NVertices -1

    Nv v[0] v[1] ....v[Nv-1] colorspec

    # colorspec continues past v[Nv-1]
    # to end-of-line may be 0 to 4 numbers
    # nothing default
    # integer : colormap index (read from the file cmap.fmap)
    # 3 or 4 integers RGB[A] values 0..255
    #
    """
    def __init__(self, _filename= 'geomoff'):
        _filename = _filename + '.off'
        Geomview.__init__(self, _filename)

    def loadpt(self):
        """ load points

        """
        fo = open(self.filename,'r')
        lis = fo.readlines()
        typ,nv,nf,ne=lis[0].split(' ')
        if typ<>'OFF':
            logging.critical('not an off file')
        nv = eval(nv)
        nf = eval(nf)
        ne = eval(ne)
        for k in range(nv):
            #x,y,z = lis[k+1].split(' ')
            pt = np.fromstring(lis[k+1],dtype=float,sep=' ')
            try:
                t = np.vstack((t,pt))
            except:
                t = pt
        return(t)

    def savept(self,ptnew,_fileoff):
        """
        """
        fo = open(self.filename,'r')
        lis = fo.readlines()
        typ,nv,nf,ne=lis[0].split(' ')
        if typ<>'OFF':
            logging.critical('not an off file')
        else:
            try:
                nv = eval(nv)
                nf = eval(nf)
                ne = eval(ne)
            except:
                logging.critical('load off wrong number of values')
        fo.close()

        fileoff = pyu.getlong(_fileoff, "geom")
        fo = open(fileoff,'w')
        fo.write(lis[0])
        for k in range(nv):
            fo.write(str(ptnew[k,0])+' '+str(ptnew[k,1])+' '+str(ptnew[k,2])+' '+'\n')
        for li in lis[k+2:]:
            fo.write(li)
        fo.close()

    def polygon(self, p, poly):
        """  create geomview off for polygon

        Parameters
        ----------
        p    : nparray
               sequence of points
        poly : list
               point numbers (index starting in 0)

        """
        fo = open(self.filename, 'w')
        npt = np.shape(p)[0]
        npoly = len(poly)
        fo.write("OFF\n")
        fo.write("%d 1 \n" % (npt + 1))
        fo.write("0.000 0.000 0.000 \n")
        for i in range(npt):
            fo.write("%6.3f %6.3f %6.3f \n" % (p[i, 0], p[i, 2], p[i, 1]))
        fo.write("%i " % (npoly - 1))
        for k in poly[:-1]:
            fo.write("%i  " % (k + 1))
        #fo.write(%6.3f %6.3f %6.3f 0.4\n" % (col[0],col[1],col[2]))
        fo.write("1.0 1.0 1.0 0.4\n")
        fo.close()

    def polygons(self, p, polys):
        """ create a gemoff file for a list of polygons

        Parameters
        ----------
        p    : nparray
               sequence of points
        poly : list
               point numbers (index starting in 0)

        Examples
        --------


        """

        fo = open(self.filename, 'w')
        npt = np.shape(p)[0]
        npoly = len(polys)
        fo.write("OFF\n")
        fo.write("%d %d \n" % (npt + 1, npoly))
        fo.write("0.000 0.000 0.000 \n")
        for i in range(npt):
            fo.write("%6.3f %6.3f %6.3f \n" % (p[i, 0], p[i, 2], p[i, 1]))

        for  poly in polys:
            nv = len(poly)
            fo.write("%i " % (nv))
            for k in poly:
                fo.write("%i  " % (k + 1))
            #fo.write(%6.3f %6.3f %6.3f 0.4\n" % (col[0],col[1],col[2]))
            fo.write("1.0 0.0 1.0 0.4\n")
        fo.close()

    def cylinder(self,r,l,nphi=20,nl=3,col=[1.,0.0,1.0],alpha=0.1):
        """ create a cylinder

        Parameters
        ----------

        r : radius
        l : length
        nphi : number of phi
        nl : number of l
        col : list [r,g,b]
        alpha : transparency

        """
        tphi = np.linspace(0,2*np.pi,nphi,endpoint=False)
        tz = np.linspace(-l/2.,l/2.,nl)
        npoly = nphi*(nl-1)
        nedges = nphi*(2*nl-1)
        fo = open(self.filename, 'w')
        #fo.write("OFF\n")
        fo.write("OFF %d %d %d\n" % (nphi*nl+1, npoly,nedges))
        fo.write("0.000 0.000 0.000 \n")
        for z in tz:
            for phi in tphi:
               x = r*np.cos(phi)
               y = r*np.sin(phi)
               fo.write("%6.3f %6.3f %6.3f \n" % (x, y, z))
        for k in range(npoly):
            il = k/nphi
            iphi = k % nphi
            a = il*nphi+iphi
            b = il*nphi+(iphi+1)%nphi
            c = (il+1)*nphi+iphi
            d = (il+1)*nphi+(iphi+1)%nphi

            fo.write("4 %i %i %i %i " %(a+1,b+1,d+1,c+1))
            str1 = str(col[0])+' '+str(col[1])+' '+str(col[2])+' '+ str(alpha)+'\n'
            fo.write(str1)

        fo.close()

    def box(self, extrem = np.array([-1,1,-1,1,-3,3])):
        """ create a box

        Parameters
        ----------

        extrem : ndarray
                 (1x6) [xmin,xmax,ymin,ymax,zmin,zmax]

        Examples
        --------

        >>> geo = Geomoff('test')
        >>> geo.box()

        """
        xmin = extrem[0]
        xmax = extrem[1]
        ymin = extrem[2]
        ymax = extrem[3]
        zmin = extrem[4]
        zmax = extrem[5]
        p = np.zeros((8, 3))
        p[0, :] = np.array([xmin, ymin, zmin])
        p[1, :] = np.array([xmax, ymin, zmin])
        p[2, :] = np.array([xmax, ymax, zmin])
        p[3, :] = np.array([xmin, ymax, zmin])
        p[4, :] = np.array([xmin, ymin, zmax])
        p[5, :] = np.array([xmax, ymin, zmax])
        p[6, :] = np.array([xmax, ymax, zmax])
        p[7, :] = np.array([xmin, ymax, zmax])
        fo = open(self.filename, 'w')
        fo.write("OFF\n")
        fo.write("8 6 12\n")
        for i in range(8):
            fo.write("%6.3f %6.3f %6.3f \n" % (p[i, 0], p[i, 2], p[i, 1]))

        fo.write("4 0 1 2 3 1 0 0 0.3\n")
        fo.write("4 7 4 0 3 1 0 0 0.3\n")
        fo.write("4 4 5 1 0 1 0 0 0.3\n")
        fo.write("4 5 6 2 1 1 0 0 0.3\n")
        fo.write("4 3 2 6 7 1 0 0 0.3\n")
        fo.write("4 6 5 4 7 1 0 0 0.3\n")
        fo.close()

    def pattern(self,theta,phi,E,**kwargs):
        """ export antenna pattern in a geomview format

        Parameters
        ----------

        theta : np.array (,Nt)
        phi   : np.array (,Np)
        E     : np.array complex  (Nt,Np)
        po    : origin (1x3)
        T     : rotation matrix (3x3)
        minr  : radius of minimum
        maxr  : radius of maximum
        ilog  : True (log) False (linear)

        Examples
        --------

           >>> from pylayers.util.geomutil import *
           >>> import numpy as np
           >>> th = np.arange(0,np.pi,0.05)
           >>> ph = np.arange(0,2*np.pi,0.05)
           >>> E = 1.5*np.sin(th[:,np.newaxis])*np.cos(0*ph[np.newaxis,:])
           >>> g = Geomoff('dipole')
           >>> g.pattern(th,ph,E)
           >>> g.show3()

        """

        defaults = { 'po': np.array([0,0,0]),
                     'T' : np.eye(3),
                     'minr' : 0.1,
                     'maxr' : 1 ,
                     'tag' : 'Pat',
                     'ilog' : False}

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        minr = kwargs['minr']
        maxr = kwargs['maxr']
        tag  = kwargs['tag']
        ilog = kwargs['ilog']
        po = kwargs['po']
        # T is an unitary matrix
        T  = kwargs['T']
        assert (abs(la.det(T))>0.99)
        # retrieving dimensions
        Nt = len(theta)#np.shape(theta)[0]
        Np = len(phi)#np.shape(phi)[1]

        theta = theta[:,np.newaxis]
        phi = phi[np.newaxis,:]

        if ilog:
            R = 10 * np.log10(abs(E))
        else:
            R = abs(E)

        #Th = np.outer(theta, np.ones(Np))
        #Ph = np.outer(np.ones(Nt), phi)

        if R.min()!=R.max():
            U = (R - R.min()) / (R.max() - R.min())
            Ry = minr + (maxr-minr) * U
        else:
            Ry = maxr
        # x (Nt,Np)
        # y (Nt,Np)
        # z (Nt,Np)

        x = Ry * np.sin(theta) * np.cos(phi)
        y = Ry * np.sin(theta) * np.sin(phi)
        z = Ry * np.cos(theta) * np.ones(phi.shape)

        # p : Nt x Np x 3
        p = np.concatenate((x[...,np.newaxis],y[...,np.newaxis],z[...,np.newaxis]),axis=2)
        #
        # antenna cs -> glogal cs
        # q : Nt x Np x 3
        q = np.einsum('ij,klj->kli',T,p)
        #
        # translation
        #
        q[...,0]=q[...,0]+po[0]
        q[...,1]=q[...,1]+po[1]
        q[...,2]=q[...,2]+po[2]

        Npoints = Nt * Np
        Nfaces = (Nt - 1) * Np
        Nedge = 0
        #
        # Colormap
        #
        colmap = plt.get_cmap()
        Ncol = colmap.N
        cmap = colmap(np.arange(Ncol))

        if R.min()!=R.max():
            g = np.round(U * (Ncol - 1)).astype(int)
        else:
            g = np.round(np.ones((Nt,Np))*(Ncol-1)).astype(int)

        fd = open(self.filename, 'w')
        fd.write('COFF\n')
        chaine = str(Npoints) + ' ' + str(Nfaces) + ' ' + str(Nedge) + '\n'
        fd.write(chaine)

        for ii in range(Nt):
            for jj in range(Np):
                cpos = str(q[ii, jj,0]) + ' ' + str(q[ii, jj,1]) + ' ' + str(q[ii, jj,2])
                cpos = cpos.replace(',', '.')
                ik = g[ii, jj]

                ccol = str(cmap[ik, 0]) + ' ' + str(cmap[ik, 1]) + \
                    ' ' + str(cmap[ik, 2])
                ccol = ccol.replace(',', '.')
                fd.write(cpos + ' ' + ccol + ' 0.2\n')

        for ii in range(Nt - 1):
            for jj in range(Np):
                p1 = ii * Np + jj
                p2 = ii * Np + np.mod(jj + 1, Np)
                p3 = (ii + 1) * Np + jj
                p4 = (ii + 1) * Np + np.mod(jj + 1, Np)
                chaine = '4 ' + str(p1) + ' ' + str(p2) + ' ' + \
                    str(p4) + ' ' + str(p3) + ' 0.5\n'
                fd.write(chaine)

        fd.close()

def angular(p1, p2):
    """ determine angle between p1 and p2 in [0 2pi]

    Parameters
    ----------
    p1
        point p1
    p2
        point p2 origin

    Notes
    -----

    weird the origin is p2

    Examples
    --------

    >>> import numpy as np
    >>> p1 = np.array([0,0])
    >>> p21 = np.array([1,0])
    >>> p22 = np.array([1,1])
    >>> p23 = np.array([0,1])
    >>> p24 = np.array([-1,1])
    >>> p25 = np.array([-1,0])
    >>> p26 = np.array([-1,-1])
    >>> p27 = np.array([0,-1])
    >>> p28 = np.array([1,-1])
    >>> a1  = angular(p21,p1)
    >>> a2  = angular(p22,p1)
    >>> a3  = angular(p23,p1)
    >>> a4  = angular(p24,p1)
    >>> a5  = angular(p25,p1)
    >>> a6  = angular(p26,p1)
    >>> a7  = angular(p27,p1)
    >>> a8  = angular(p28,p1)


    """
    #print DeprecationWarning('DEPRECATION WARNING : geomutil.angular going
    #                         deprecated  because wrong')
    if p1[0] < p2[0] and p1[1] < p2[1]:
        angle = np.arctan2((p2[1] - p1[1]), (p2[0] - p1[0])) + np.pi
    elif p1[0] > p2[0] and p1[1] < p2[1]:
        angle = np.arctan2((p2[1] - p1[1]), (p2[0] - p1[0])) + np.pi
    elif p1[0] > p2[0] and p1[1] > p2[1]:
        angle = np.arctan2((p2[1] - p1[1]), (p2[0] - p1[0])) + np.pi
    else:
        angle = np.arctan2((p2[1] - p1[1]), (p2[0] - p1[0])) + np.pi

    return(angle)



def vecang(v1,v2):
    """ angle between v1 and v2 , result in [0,2*pi]

    Parameters
    ----------
    v1 : np.array (3 x Np)
        vector
    v2 : np.array (3 x Np)
        vector
    Returns
    -------

    alpha : np.array (3 x Np)
        radians


    """
    if len(v1.shape) == 1:
        v1=v1.reshape(v1.shape[0],1)
    if len(v2.shape) == 1:
        v2=v2.reshape(v2.shape[0],1)


    ang =np.arctan2(v2[1,:], v2[0,:]) - np.arctan2(v1[1,:], v1[0,:])
    uneg = np.where(ang <0)[0]
    ang[uneg]= 2*np.pi+ ang[uneg]
    return ang
    # if ang <0 :
    #     return (2*np.pi+ang)
    # else :
    #     return ang


def SignedArea(p=np.array([[0, 10, 10, 0], [0, 0, -2, -2]])):

    """
        Calculate the signed area of a sequence of points in a  plane

        Parameters
        ----------
        p : array 2 x Np

        Examples
        --------

        >>> from pylayers.util.geomutil import *
        >>> p = np.array([[0,10,10,0],[0,0,-2,-2]])
        >>> A = SignedArea(p)
        >>> assert(A+20<1e-15)

    """
    return sum(np.hstack((p[0, 1::], p[0, 0:1])) * (np.hstack((p[1, 2::], p[1, 0:2])) - p[1, :])) / 2.


def Centroid(p=np.array([[0, 10, 10, 0], [0, 0, -2, -2]])):
    """ Determine the centroid of the polygon defined by a sequence of points in a plane

    References
    ----------
    http://en.wikipedia.org/wiki/Centroid

    Parameters
    ----------
        p : np array polygon (2xNp)

    Returns
    -------
        pc = Centroid()

    Examples
    --------

    >>> from pylayers.util.geomutil import *
    >>> p  = np.array([[0,10,10,0],[0,0,-2,-2]])
    >>> pc = Centroid(p)
    >>> d  = pc-np.array([5.,-1])
    >>> md = np.dot(d,d)
    >>> assert(md<1e-15)

    """
    A = SignedArea(p)
    assert(A<>0)
    T = p[0, :] * np.hstack((p[1, 1::], p[1, 0:1])) - \
        p[1, :] * np.hstack((p[0, 1::], p[0, 0:1]))
    Cx = sum(T * (p[0, :] + np.hstack((p[0, 1::], p[0, 0:1])))) / (6 * A)
    Cy = sum(T * (p[1, :] + np.hstack((p[1, 1::], p[1, 0:1])))) / (6 * A)
    pc = np.array([Cx, Cy])
    return(pc)


def Lr2n(p=np.array([[0, 10, 10, 0], [0, 0, -2, -2]]), closed=True):
    """  Linear Ring to normal

        Parameters
        ----------
        p      : np.array (2xN)
        closed : boolean
            default True

        Returns
        -------
        n   : np.array (2xN)
              normal

        Notes
        -----
            This function returns the internal normals to the LinearString of a Polygon

            The algoritm exploits the algebraic relation which exists
            between points coordinates and normal coordinates which involves
            the quasi toeplitz matrix M

            [-1  1  0 0 0  ...]
            [0  -1  1 0 0  ...]
            [
            [
            [           0 -1 1]  (truncate here if LineRing is open)
            -------------------
            [1          0 0 -1]  (add this line if LineRing is closed)

            p0                   p1
            x------------------x
            |        |         |
            |        v         |
            |                  |
            |->              <-|
            |        ^         |
            |        |         |
         p3 x------------------x p2

    Examples
    --------

    >>> import shapely.geometry as shg
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> points1 = shg.MultiPoint([(0, 0), (0, 1), (1, 1), (1,0 )])
    >>> points2 = shg.MultiPoint([(0, 0), (1, 0), (1, 1), (0,1 )])
    >>> poly1   = shg.Polygon(points1)
    >>> poly2   = shg.Polygon(points2)
    >>> lring1  = poly1.exterior
    >>> lring2  = poly2.exterior
    >>> x1,y1    = lring1.xy
    >>> x2,y2    = lring2.xy
    >>> p1 = np.array([x1[0:-1],y1[0:-1]])
    >>> p2 = np.array([x2[0:-1],y2[0:-1]])
    >>> n1 = Lr2n(p1)
    >>> n2 = Lr2n(p2)

    """
    Np = np.shape(p)[1]
    l = np.hstack((np.array([-1, 1]), np.zeros(Np - 2)))
    M = np.triu(toeplitz(l))
    if closed:
        M[Np - 1, 0] = 1
    else:
        M = M[0:Np - 1, :]

    n = np.dot(M, np.flipud(p).T)
    n[:, 1] = -n[:, 1]
    #
    # normalize normal
    #
    n = n.T
    modn = np.sqrt(np.sum(n * n, 0))
    assert(modn.all()>0)
    nn = n / modn
    #
    # enforce inwards normal whatever the linear ring orientation
    #
    sa = SignedArea(p)
    if sa > 0:
        nn = -nn
    return nn


def isBetween(p1, p2, p, epsilon=1e-5):
    """ test if p is between p1 and p2

    Parameters
    ----------

    p1 : np.array
    p2 : np.array
    p  : np.array
    epsilon : float
        tolerance default 1e-5

    Returns
    -------
    boolean

    Examples
    --------

    >>> p1 = np.array([0,0])
    >>> p2 = np.array([2,0])
    >>> p  = np.array([1,0])
    >>> assert(isBetween(p1,p2,p)),'error'

    """
    crossproduct = np.cross(p - p1, p2 - p1)
    if abs(crossproduct) > epsilon:
        return False
    dotproduct = np.dot(p - p1, p2 - p1)
    if dotproduct < 0:
        return False
    squaredlengthba = np.dot(p2 - p1, p2 - p1)
    if dotproduct > squaredlengthba:
        return False
    else:
        return True


def pvec(v1, v2):
    """ cross product between v1 and v2

     Parameters
     ----------
     v1 : numpy array
     v2 : numpy array
     Returns
     -------
     v3 = v1 x v2

     See Also
     --------
     np.cross

     Examples
     --------
     >>> v1 = np.array([1,0,0])
     >>> v2 = np.array([0,1,0])
     >>> v3 = pvec(v1,v2)

    """
    A = np.array(
        [[0., -v1[2], v1[1]], [v1[2], 0., -v1[0]], [-v1[1], v1[0], 0.]])
    v3 = np.dot(A, v2)
    return(v3)


def pvecn(v1, v2):
    """ cross product and normalization

    Parameters
    ----------

     v1 : numpy array
     v2 : numpy array

    Returns
    -------

     v3 = v1 x v2 / | v1 x v2 |

    Examples
    --------

    >>> v1 = np.array([2,0,0])
    >>> v2 = np.array([0,2,0])
    >>> v3 = pvecn(v1,v2)

    See Also
    --------

        numpy.cross

    """
    v3 = np.cross(v1, v2)
    try:
        v4 = v3 / np.sqrt(np.dot(v3, v3))
    except:
        print("error divide by zero in pvecn")
    return(v4)

def onb(A,B,v):
    """ orthonormal basis from 2 points defining an axe and a vector

    Parameters
    ----------

    A : np.array
        3 x n
    B : np.array
        3 x n
    v : np.array
        3 x n

    Returns
    -------


    T basis (un,vn,wn)
        3 x n x 3
    (un,vn) is a basis in the plane transverse to the axis vn
    wn is the unitary vector along vector AB

    Examples
    --------

    >>> A = np.array([[0,0,0,0],[1,2,3,4],[0,0,0,0]])
    >>> B = np.array([[0,0,0,0],[1,2,3,4],[10,10,10,10]])
    >>> v = np.array([[1,1,1,1],[0,0,0,0],[0,0,0,0]])
    >>> onb(A,B,v)
    array([[[ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.]]
    <BLANKLINE>
           [[ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.]],
    <BLANKLINE>
           [[ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.]],
    <BLANKLINE>
           [[ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.]]])


    see also
    --------

    pylayers.util.geomutil.Geomvect.geomBase
    pylayers.util.mobility.body

    """
    #np.random.seed(0)
    N = np.shape(A)[1]
    # modab 1xN
    modab = np.sqrt(np.sum((B-A)*(B-A),axis=0))
    # wn 3xN
    wn = (B - A) / modab
    #random_vector = np.random.rand(3,N)
    u = v - np.sum(v*wn,axis=0)*wn
    modu = np.sqrt(np.sum(u*u,axis=0))
    # un : 3xN
    un = u /modu
    # vn : 3xN
    vn = np.cross(wn,un,axis=0)
    #pdb.set_trace()
    T  = np.dstack((un,vn,wn))
    # reshape dimension for having index of cylinder axe first
    # N x 3 x 3
    T  = T.swapaxes(0,1)
    return T


def vec_sph(th, ph):
    """
    vec_sph(th,ph)

    return Spherical orthonormal frame

    [ [ eth]
      [ eph]    (theta,phi)
      [ er ] ]

    """

    e_th = np.array(
        (np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)))
    e_ph = np.array((-np.sin(ph), np.cos(ph), 0))
    e_r = np.array(
        (np.cos(ph) * np.sin(th), np.sin(ph) * np.sin(th), np.cos(th)))

    B = np.vstack((e_th, e_ph, e_r))
    return(B)


def ellipse(fd, p, vth, vph, Eth, Eph, N):
    """ build a geomview file of an ellipse

     Parameters
     ----------

     fd    : file descriptor
     p     : ellipse center
     vth   : unitary vector along theta
     vph   : unitary vector along phi
     Eth   : complex
     Eph   : complex
     N     : descretization step

    """

    pas = 2 * np.pi / N
    alpha = np.linspace(0, 2 * np.pi - pas, N)

    Rth = abs(Eth)
    Rph = abs(Eph)
    delta_th = np.arctan2(np.imag(Eth), np.real(Eth))
    delta_ph = np.arctan2(np.imag(Eph), np.real(Eph))

    pu1 = p + Rth * vth
    pu2 = p + Rph * vph

    u3 = np.ones(3)
    uN = np.ones(N)
    Al_th = np.outer(u3, alpha + delta_th)
    Al_ph = np.outer(u3, alpha + delta_ph)

    U1 = np.outer(vth, uN)
    U2 = np.outer(vph, uN)

    P = np.outer(p, uN)
    #
    # Un point de l'ellipse
    #
    pc = P + (Rth * U1 * np.cos(Al_th) + Rph * U2 * np.cos(Al_ph))

    vEre = p + (np.real(Eth) * vth + np.real(Eph) * vph)
    vEim = p + (np.imag(Eth) * vth + np.imag(Eph) * vph)

    fd.write("appearance { linewidth 3 }\n")
    fd.write("VECT\n")
    fd.write("%d %d %d \n" % (N, 2 * N, N))
    fd.write("\n")
    for i in range(N):
        fd.write("%d " % 2)
    fd.write("\n")
    for i in range(N):
        fd.write("%d " % 1)
    fd.write("\n")
    for i in range(N - 1):
        fd.write("%6.3f %6.3f %6.3f\n" % (pc[0, i], pc[1, i], pc[2, i]))
        fd.write("%6.3f %6.3f %6.3f\n" % (pc[0, i + 1], pc[1, i + 1], pc[2, i + 1]))
        fd.write("\n")
    fd.write("%6.3f %6.3f %6.3f\n" % (pc[0, N - 1], pc[1, N - 1], pc[2, N - 1]))
    fd.write("%6.3f %6.3f %6.3f\n" % (pc[0, 0], pc[1, 0], pc[2, 0]))
    fd.write("\n")
    for i in range(N):
        v = float(i - 1) / N
        fd.write("%g %g %g %g\n" % (v, v, v, 1))

def normalize(vec):
    """ normalize an array of N ndim  vectors

    Parameters
    ----------

    vec : ndarray (N x ndim)
        N ndim vectors

    Returns
    -------
    vecn : ndarray (N x ndim)
        N normalized ndim vectors

    Example
    -------

    >>> from pylayers.util.geomutil import *
    >>> vec = np.array([[1,1,0],[1,1,0],[1,0,1],[1,1,1]])
    >>> normalize(vec)
    array([[ 0.70710678,  0.70710678,  0.        ],
           [ 0.70710678,  0.70710678,  0.        ],
           [ 0.70710678,  0.        ,  0.70710678],
           [ 0.57735027,  0.57735027,  0.57735027]])

    Notes
    -----

    """
    N = np.shape(vec)[0]
    m = np.sqrt(np.sum(vec*vec,axis=1)).reshape(N,1)
    vecn = vec/m

    return(vecn)

def ptonseg(pta, phe, pt):
    """ return a point on the segment (pta,pte)

    Parameters
    ----------
    pta : ndarray
    phe : ndarray
    pt  : ndarray

    Returns
    -------
    p   : ndarray

    Example
    -------

    """
    v = phe - pta
    u = pt - pta
    Lv = np.sqrt(np.dot(v, v))
    Lu = np.sqrt(np.dot(u, u))
    assert(Lv<>0)
    assert(Lu<>0)
    vn = v / Lv
    un = u / Lu
    ctheta = np.dot(un, vn)
    alpha = ctheta * Lu
    if (alpha > 0) & (alpha < Lv):
        p = pta + alpha * vn
    else:
        p = np.array([])
    return p

def dptseg(p,pt,ph):
    """ distance between a set of points and a segment

    Parameters
    ----------
    ps  : ndim x Np
          array of Np points
    pt  : ndim x 1
          tail coordinates of segment
    ph  : ndim x 1
          head coordinates of segment

    Returns
    -------
    d1 : 1 x Np
        distance between pt and ortho projection of ps
    d2 : 1 x Np
        distance between ph and ortho projection of ps
    h  : distance between ps and ortho projection of ps

    Examples
    --------

    .. plot::
        :include-source:
        >>> import numpy as np
        >>> from pylayers.util.geomutil import *
        >>> pt = np.array([0,0])
        >>> ph = np.array([10,0])
        >>> p  = np.array([[-1,1 ,3,4,11],[8,1,2,3,3]])
        >>> d1,d2,h = dptseg(p,pt,ph)

    """
    ndim = len(pt)
    l = ph.reshape(ndim,1)-pt.reshape(ndim,1)
    norml = np.sqrt(np.dot(l.T,l))
    ln = l/norml

    ptp = p - pt.reshape(2,1)
    d1 = np.dot(ln.T,ptp)
    d2 = norml - d1
    ptpl = d1*ln
    ptpo = ptp - ptpl
    h = np.sqrt(np.sum(ptpo*ptpo,axis=0))

    return(d1,d2,h)


def linet(ax, p1, p2, al=0.9, color='blue', linewidth=1):
    """  draw a short line segment

    Parameters
    ----------

    ax : axes
    p1 : np.array
        start point
    p2 : np.array
        end point
    al : float
        0 < al < 1   percentage of drawing  default 0.9
    color :  string
        color default 'blue'
    linewidth : float
        line width default 1

    Returns
    -------

    ax  : Axes instance

    Examples
    --------

    .. plot::
        :include-source:

        >>> from pylayers.util.geomutil import *
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure()
        >>> ax  = fig.gca()
        >>> p1 = np.array([0,0])
        >>> p2 = np.array([1,0])
        >>> p3 = np.array([0,1])
        >>> p4 = np.array([1,1])
        >>> ax = linet(ax,p1,p2,al=0.7,color='red',linewidth=3)
        >>> ax = linet(ax,p2,p3,al=0.8,color='blue',linewidth=2)
        >>> ax = linet(ax,p3,p4,al=0.9,color='green',linewidth=1)
        >>> ax = linet(ax,p4,p1,al=1,color='cyan',linewidth=0.2)


    """
    v = p2 - p1
    L = np.sqrt(np.dot(v, v))
    vn = v / L
    pi = p1 + vn * (1 - al) * L
    pf = p2 - vn * (1 - al) * L
    ax.plot([pi[0], pf[0]], [pi[1], pf[1]], color=color, linewidth=linewidth)

    return(ax)


def ccw(a, b, c):
    """ counter clock wise order

    Parameters
    ----------

    a : ndarray (2,N)
    b : ndarray (2,N)
    c : ndarray (2,N)

    Returns
    -------

    array of booleans

    References
    ----------

    `Line Segment Intersection <http://www.bryceboe.com/2006/10/23/line-segment-intersection-algorithm/>`_


    Examples
    --------

    >>> import scipy as sp
    >>> a = sp.rand(2,100)
    >>> b = sp.rand(2,100)
    >>> c = sp.rand(2,100)
    >>> u = ccw(a,b,c)

    """
    assert a.shape[0]==2
    assert b.shape[0]==2
    assert c.shape[0]==2
    #return((c[1, :] - a[1, :]) * (b[0, :] - a[0, :]) > (b[1, :] - a[1, :]) * (c[0, :] - a[0, :]))
    return((c[1, ...] - a[1, ...]) * (b[0, ...] - a[0, ...]) >
           (b[1, ...] - a[1, ...]) * (c[0, ...] - a[0, ...]))

def intersect_line_seg(line,seg):
    """ intersect a line and a segment
    
    Parameters 
    ----------
    
    line : (point,vec)
    seg :  (pta,phe)
    
    Returns
    -------
    
    k : intersection parameter (0<k<1 if intersection)
    M : intersection point 
    
    """
    pt, v = line
    pta,phe = seg
    vseg = phe-pta
    xht = phe[0]-pta[0]
    yth = pta[1]-phe[1]
    num = -(v[1]*(pta[0]-pt[0])+v[0]*(pt[1]-pta[1]))
    den = (v[1]*xht+v[0]*yth)
    if (abs(den) > 0):
        k = num/den
    else:
        k = 2
    M = pta+k*vseg
    return(k,M)


def intersect(a, b, c, d):
    """ check if segment AB intersects segment CD

    Parameters
    ----------

    a : np.array (2xN)
    b : np.array (2xN)
    c : np.array (2xN)
    d : np.array (2xN)


    Examples
    --------

    .. plot::
        :include-source:

        >>> import scipy as sp
        >>> import numpy as np
        >>> from pylayers.util.geomutil import *
        >>> from pylayers.util.plotutil import *
        >>> import matplotlib.pylab as plt
        >>> N = 10
        >>> A = sp.rand(2,N)
        >>> B = sp.rand(2,N)
        >>> C = sp.rand(2,N)
        >>> D = sp.rand(2,N)
        >>> b1 = intersect(A,B,C,D)
        >>> pt1 = A[:,b1]
        >>> ph1 = B[:,b1]
        >>> pt2 = C[:,b1]
        >>> ph2 = D[:,b1]
        >>> f1,a1 = displot(pt1,ph1,'r')
        >>> f2,a2 = displot(pt2,ph2,'b')
        >>> ti = plt.title('test intersect')
        >>> A = np.array([[0],[0]])
        >>> B = np.array([[1],[1]])
        >>> C = np.array([[1],[0]])
        >>> D = np.array([[0],[1]])
        >>> intersect(A,B,C,D)
        array([ True], dtype=bool)
        >>> intersect(A,B,C,D)[0]
        True

    See Also
    --------

    ccw : counter clock wise detection

    """
    return ((ccw(a, c, d) != ccw(b, c, d)) & (ccw(a, b, c) != ccw(a, b, d)))



def isaligned(a,b,c):
    #return abs(((b[0,:]-a[0,:])*(c[1,:]-a[1,:]) - (b[1,:]-a[1,:])*(c[0,:]-a[0,:])))<1e-8
    val = abs(((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])))
    cond = val<1e-2
    #print val
    return cond

def isleft(a,b,c,tol=0.):
    """ Test point c is at left of the vector a-->b


    Parameters
    ----------

    a : np.array (2xN)
    b : np.array (2xN)
    c : np.array (2xN)
    tol : tolerance
    Returns
    -------

    boolean array (1xN)

    Examples
    --------

    .. plot::
        :include-source:

        >>> from pylayers.util.plotutil import *
        >>> import scipy as sp
        >>> import numpy as np
        >>> from pylayers.util.geomutil import *
        >>> from pylayers.util.plotutil import *
        >>> import matplotlib.pylab as plot
        >>> N = 20
        >>> A = sp.rand(2,N)
        >>> B = sp.rand(2,N)
        >>> C = np.array(([0.5,0.5])).reshape(2,1)
        >>> left=isleft(A,B,C)
        >>> il = np.where(left)[0]
        >>> inl = np.where(~left)[0]
        >>> plt.scatter(C[0],C[1],color='b',s=10)
        >>> displot(A[:,il],B[:,il],arrow=True,color='g')
        >>> displot(A[:,inl],B[:,inl],arrow=True,color='r')

    See Also
    --------

    pylayers.antprop.signature

    """
    return ((b[0,:]-a[0,:])*(c[1,:]-a[1,:])) - ((b[1,:]-a[1,:])*(c[0,:]-a[0,:]))>tol

def isleftorequal(a,b,c):
    return ((b[0,:]-a[0,:])*(c[1,:]-a[1,:])) - ((b[1,:]-a[1,:])*(c[0,:]-a[0,:]))>=0


def affine(X,Y):
    """ find affine transformation

    Parameters
    ----------

    X  : np.array
        3xN
    Y
        3xN

    Returns
    -------

    A : np.array
        3x3
    B : np.array
        3x1

    Notes
    -----

    Given X and Y find the affine transformation

    Y = A X + B

    """

    B = Y[:,0][:,np.newaxis]
    Yc = Y-B
    pX = la.pinv(X)
    A = np.dot(Yc,pX)
    return(A,B)

def cylmap(Y,r=0.0625,l=0.5):
    """ find affine transformation for a specific cylinder

    Parameters
    ----------

    Y
        3xN

    Returns
    -------

    A : np.array
        3x3
    B : np.array
        3x1

    Notes
    -----

    Y = A X + B

    """
    #X = np.array([[0,0,0],[0,0,-0.25],[0,0,0.25],[0.0625,0,0],[0,0.0625,0],[0.0625,0,0.25]]).T
    X = np.array([[0,0,0],[0,0,-l/2],[0,0,l/2],[r,0,0],[0,r,0],[r,0,l/2]]).T
    B = Y[:,0][:,np.newaxis]
    Yc = Y-B
    pX = la.pinv(X)
    A = np.dot(Yc,pX)
    return(A,B)

def mul3(A, B):
    """  matrix multiplication

    Parameters
    ----------

    A :
    B :

    Returns
    -------

    C :  A*B

    """
    sa = np.shape(A)
    sb = np.shape(B)
    la = len(sa)
    lb = len(sb)
    if ((la == 3) & (lb == 3)):
        if((sa[1] == sb[0]) & (sa[2] == sb[2])):
            C = np.zeros((sa[0], sb[1], sb[2]))
            for i in range(sa[2]):
                MA = A[:, :, i]
                MB = B[:, :, i]
                P = np.dot(MA, MB)
                C[:, :, i] = P
            return(C)
        else:
            print "wrong shape", sa, sb
    if ((la == 3) & (lb == 2)):
        if(sa[1] == sb[0]):
            C = np.zeros((sa[0], sb[1], sa[2]))
            for i in range(sa[2]):
                MA = A[:, :, i]
                P = np.dot(MA, B)
                C[:, :, i] = P
            return(C)
        else:
            print "wrong shape", sa, sb

    if ((la == 2) & (lb == 3)):
        if(sa[1] == sb[0]):
            C = np.zeros((sa[0], sb[1], sb[2]))
            for i in range(sb[2]):
                MB = B[:, :, i]
                P = np.dot(A, MB)
                C[:, :, i] = P
            return(C)
        else:
            print "wrong shape", sa, sb


def MRot3(a, axe):
    """
    Return a 3D rotation matrix along axe 0|1|2

    Parameters
    ----------

    a   :  angle (radians)
    axe :  0:x 1:y 2:z

    """
    M3 = np.eye(3)
    M2 = np.array(((np.cos(a), -np.sin(a)), (np.sin(a), np.cos(a))))
    if (axe == 0):
        M3[1:3, 1:3] = M2
    if (axe == 1):
        M3[0::2, 0::2] = M2
    if (axe == 2):
        M3[0:2, 0:2] = M2
    return(M3)


def MEulerAngle(alpha, beta, gamma):
    """ Calculate a rotation matrix from 3 Euler angles

    Parameters
    ----------

    alpha  : float
        rotation along axis z
    beta : float
        rotation along axis x
    gamma : float
        rotation along axis y

    Returns
    -------

    T    : np.array (3x3)
        rotation matrix

    Examples
    --------

    >>> import numpy as np
    >>> T=MEulerAngle(np.pi/2,np.pi/2,np.pi/2)

    Warnings
    --------

    Bizarre I was expected

    -1  0  0
     0  0  1
     0  1  0

    """
    Ra = MRot3(alpha, 2)
    Rb = MRot3(beta, 0)
    Rg = MRot3(gamma, 1)

    T = np.dot(np.dot(Ra, Rb), Rg)
    #T  = np.dot(np.dot(Rg,Rb),Ra)
    return(T)


def SphericalBasis(a):
    """
    SphericalBasis(a):

    a[:,0] : N x theta   theta angle
    a[:,1] : N x phi     phi angle

    M = N x [th,ph,s]  : 3 x 3 x N
    """

    tha = np.vstack((np.cos(a[:, 0]) * np.cos(
        a[:, 1]), np.cos(a[:, 0]) * np.sin(a[:, 1]), -np.sin(a[:, 0]))).T
    pha = np.vstack((-np.sin(a[:, 1]), np.cos(a[:, 1]), 0 * a[:, 0])).T
    sa = np.vstack((np.sin(a[:, 0]) * np.cos(
        a[:, 1]), np.sin(a[:, 0]) * np.sin(a[:, 1]), np.cos(a[:, 0]))).T

    M = np.dstack((tha, pha, sa)).T
    return M


def angledir(s):
    """ evaluate (theta,phi) from direction vector


    Parameters
    ----------

    s  : ndarray N x 3
         N direction vector

    Returns
    -------

      a  : ndarray 2xN
           N angle (theta,phi)

    Notes
    -----

    .. math::

        \\theta = \\arccos{(\\frac{\\mathbf{s}}{\\hat{\mathbf{z}})}}

    Example
    -------

    .. plot::
        :include-source:

        >>> import numpy as np
        >>> s = np.array([[2,0,0],[0,2,0],[0,0,1],[1,1,1]])
        >>> angledir(s)*180/np.pi
        array([[ 90.        ,   0.        ],
               [ 90.        ,  90.        ],
               [  0.        ,   0.        ],
               [ 54.73561032,  45.        ]])


    See Also
    --------

    BTB_Rx
    BTB_Tx

    """
    s = normalize(s)
    N = np.shape(s)[0]
    x = np.array((1, 0, 0)).reshape(1,3)
    y = np.array((0, 1, 0)).reshape(1,3)
    z = np.array((0, 0, 1)).reshape(1,3)
    u = np.dot(s,z.T)
    theta = np.arccos(u)
    v = s - z
    n = np.sqrt(np.sum(v * v, axis=1)).reshape(N,1)
    inull = np.where(n==0)[0]
    n[inull] = 1
    vn = v / n
    vnx = np.dot(vn, x.T)
    vny = np.dot(vn, y.T)
    phi = np.arctan2(vny, vnx)
    a_new = np.hstack((theta, phi))
    a_new[inull,0]=0
    a_new[inull,1]=0
    return(a_new)


def BTB_rx(a_g, T):
    """ Produce a set of rotation matrices for passage between global and
    local frames

    Parameters
    ----------

    a_g  :
        angle in global reference frame   2 x N  :  (theta,phi) x N
    T    :
        Rx rotation matrix     3 x 3

    Returns
    -------

    R  :  ndarray (3x3)
    al :  ndarray (r x 2)
        angle expressed in local basis

    See Also
    --------

    angledir
    SphericalBasis


    Notes
    -----

    N is the number or rays

    """
    G = SphericalBasis(a_g)
    th_g = G[0, :, :]
    ph_g = G[1, :, :]
    B_g = np.dstack((th_g, ph_g)).transpose((0, 2, 1))
    s_l = np.dot(T.T, G[2, :, :]).T
    al = angledir(s_l)
    L = SphericalBasis(al)
    th_l = L[0, :, :]
    ph_l = L[1, :, :]
    B_lT = np.dstack((th_l, ph_l)).transpose((2, 0, 1))
    R = mul3(B_lT, mul3(T.T, B_g))

    return R, al


def BTB_tx(a_g, T):
    """ Produce a set of rotation matrices for passage between global and local frame

    Parameters
    ----------

    a_g  : angle in global reference frame      2 x N  :  (theta,phi) x N
    T    : Tx rotation matrix     3 x 3

    """
    G = SphericalBasis(a_g)

    th_g = G[0, :, :]
    ph_g = G[1, :, :]

    B_gT = np.dstack((th_g, ph_g)).transpose((2, 0, 1))

    s_l = np.dot(T.T, G[2, :, :]).T

    al = angledir(s_l)

    L = SphericalBasis(al)
    th_l = L[0, :, :]
    ph_l = L[1, :, :]
    B_l = np.dstack((th_l, ph_l)).transpose((0, 2, 1))
    R = mul3(B_gT, mul3(T, B_l))

    return R, al

def plot_coords(ax, ob, color='#999999'):
    """ plotting coord of a `shapely` object

    Parameters
    ----------
    ax : matplotlib axes
    ob : shapely object
    color : string
        default '#999999'

    References
    ----------
    `Shapely <http://pypi.python.org/pypi/Shapely>`_

    """
    x, y = ob.xy
    ax.plot(x, y, 'o', color=color, zorder=2)  # color='#999999'


def plot_bounds(ax, ob, color='#000000'):
    """ plot bounds

    Parameters
    ----------

    ax : matplotlib axes
    ob : shapely object
    color : string
        default '#999999'


    References
    ----------
    `Shapely <http://pypi.python.org/pypi/Shapely>`_

    """
    x, y = zip(*list((p.x, p.y) for p in ob.boundary))
    ax.plot(x, y, color=color, zorder=0.1)  # '#000000'
    #ax.plot(x, y, 'o', color='#000000', zorder=0.1)   #'#000000'


def plot_line(ax, ob, color="#999999"):
    """ plot line

    Parameters
    ----------
    ax : matplotlib axes
    ob : shapely object
    color : string
        default '#999999'

    References
    ----------

    `Shapely <http://pypi.python.org/pypi/Shapely>`_

    Notes
    -----

    color = v_color(ob)

    Examples
    --------
    .. plot:
        :include-source:

    >>> from pylayers.util.geomutil import *
    >>> import matplotlib.pyplot as plt
    >>> seg = shg.LineString([(0,0),(1,1)])
    >>> fig = plt.figure()
    >>> ax  = fig.gca()
    >>> plot_line(ax,seg)
    >>> plt.show()

    """
    x, y = ob.xy
    ax.plot(x, y, color=color, alpha=0.7, linewidth=2,
            solid_capstyle='round', zorder=0.5)


def v_color(ob):
    """ return color

    Parameters
    ----------

    ob :

    References
    ----------

    http://pypi.python.org/pypi/Shapely
    """
    return COLOR[ob.is_simple]




#def createPolygons(
def plotPolygon(poly, color="#abcdef", alpha=0.8):
    """ plot a shapely Polygon

    Parameters
    ----------

    poly  : shapely polygon
    color : default "#abcdef"
    alpha : float
           transparency   (default 0.8)
    """
    fig = plt.gcf()
    gax = fig.get_axes()
    if len(gax) != 0:
        ax = gax[0]
    else:
        ax = fig.add_subplot(111)
    patch = PolygonPatch(poly, facecolor=color, alpha=alpha)
    ax.add_patch(patch)
    plt.show()


def shrinkPolygon(poly, d=0.1):
    """ shrink polygon

    Parameters
    ----------

    poly  : shapely polygon
    d     : float
        0.1

    Returns
    -------

    poly

    """
    poly1 = simplifyPolygon(poly)
    A1 = poly1.area
    p = np.array(poly1.exterior.xy)
    # enleve le dernier point
    q = p[:, 0:-1]
    n1 = Lr2n(q)
    N = np.shape(q)[1]
    for k in range(N):
        nrm = n1[:, k]
        p[:, k] = p[:, k] + d * nrm
        p[:, k + 1] = p[:, k + 1] + d * nrm
    q[:, 0] = q[:, 0] + d * nrm
    y = q.T.copy()
    ls = shg.asLineString(y)
    poly2 = shg.Polygon(ls)
    A2 = poly2.area
    if (A2 > A1):
        poly2 = shrinkPolygon(poly1, -d)

    return(poly2)


def shrinkPolygon2(poly1, d=0.1):
    """ shrink Polygon

    Parameters
    ----------

    poly1 Polygon

    """
    poly1 = simplifyPolygon(poly)
    p = np.array(poly1.exterior.xy)
    # enleve le dernier point
    q = p[:, 0:-1]
    n1 = Lr2n(q)
    N = np.shape(q)[1]
    for k in range(N):
        norm = n1[:, k]
        if k > 0:
            u = np.dot(norm, normold)
        else:
            u = 0
            norm1 = norm
        if u < 0.8:       # changement de direction
            p[:, k] = p[:, k] + d * norm
            if k != (N - 1):
                p[:, k + 1] = p[:, k + 1] + d * norm
            else:
                v = np.dot(norm, norm1)
                if v < 0.8:
                    p[:, k + 1] = p[:, k + 1] + d * norm
        else:           # meme direction
            p[:, k] = p[:, k] + d * normold
            if k != (N - 1):
                p[:, k + 1] = p[:, k + 1] + d * norm
            else:
                v = np.dot(norm, norm1)
                if v < 0.8:
                    p[:, k + 1] = p[:, k + 1] + d * norm

        normold = norm
    #n2  = np.hstack((n1[:,-1].reshape(2,1),n1[:,0:-1]))
    y = q.T.copy()
    ls = shg.asLineString(y)
    poly2 = shg.Polygon(ls)
    return(poly2)


def simplifyPolygon(poly1):
    """ Simplify polygon : suppress adjacent colinear segments

    Parameters
    ----------

    poly1

    """
    p = np.array(poly1.exterior.xy)
    N = np.shape(p)[1]
    q = p[:, 0].reshape(2, 1)
    for k in range(N - 2):
        v1 = p[:, k + 1] - p[:, k]
        v2 = p[:, k + 2] - p[:, k + 1]
        v1n = v1 / np.sqrt(np.dot(v1, v1))
        v2n = v2 / np.sqrt(np.dot(v2, v2))
        u = np.dot(v1n, v2n)
        if u < 0.98:
            q = np.hstack((q, p[:, k + 1].reshape(2, 1)))
    vini = q[:, 1] - q[:, 0]
    vin = vini / np.sqrt(np.dot(vini, vini))
    v = np.dot(v2n, vin)
    if v > 0.98:
        q = q[:, 1:]
    y = q.T.copy()
    ls = shg.asLineString(y)
    poly2 = shg.Polygon(ls)
    return(poly2)

#----------------------
# Taguhi
#----------------------


def wall_delta(x1, y1, x2, y2, delta=0.0001):
    """ Identification of new points

    After defining a tolerance length those points which are situated in
    the extremities of the walls at a distance equivalent to the
    tolerance length are identified.

    Parameters
    ----------
    x1 :  float
        The x component of the point of the first extremity
    y1 :  float
        The x component of the point of the first extremity
    x2 :  float
        The x component of the point of the second extremity
    y2 :  float
        The x component of the point of the second extremity.

    Returns
    -------
    bx  :  float
        The x component of the new point of the first extremity
    by  :  float
        The y component of the new point of the first extremity
    cx  :  float
        The x component of the new point of the second extremity
    cy  :  float
        The y component of the new point of the second extremity.

    Notes
    -----

    .. math:: bx=x1+(x2-x1) \\frac{\\delta}{mod(a)}.

    Examples
    --------

    >>> x1=-2.
    >>> y1=2.
    >>> x2=-1.
    >>> y2=1.
    >>> bx,by,cx,cy = wall_delta(x1,y1,x2,y2,delta=0.0001)
    >>> assert bx==-1.9999292893218814,'Mistake'
    >>> assert by==1.9999292893218814,'Mistake'
    >>> assert cx==-1.0000707106781186,'Mistake'
    >>> assert cy==1.0000707106781186,'Mistake'

    """

    ax = x2 - x1
    ay = y2 - y1
    a_mod = np.sqrt(ax ** 2 + ay ** 2)
    a_ch_x = ax / a_mod
    a_ch_y = ay / a_mod

    bx = x1 + ax * delta / a_mod
    by = y1 + ay * delta / a_mod
    cx = x2 - ax * delta / a_mod
    cy = y2 - ay * delta / a_mod

    return(bx, by, cx, cy)


def plot_coords2(ax, ob):
    """  plot point from coordinates

    References
    ----------

    http://pypi.python.org/pypi/Shapely

    """
    x, y = ob.xy

    ax.plot(x, y, 'o', color='#999999', zorder=2)


def plot_bounds2(ax, ob):
    """ plot bounds v2

    References
    ----------

    http://pypi.python.org/pypi/Shapely
    """
    x, y = zip(*list((p.x, p.y) for p in ob.boundary))
    ax.plot(x, y, color='#000000', zorder=0.1)


def plot_line2(ax, ob):
    """ plot line v2

    References
    ----------
    http://pypi.python.org/pypi/Shapely
    """
    x, y = ob.xy
    ax.plot(x, y, color=v_color(
        ob), alpha=0.7, linewidth=2, solid_capstyle='round', zorder=0.5)


def plot_coords3(ax, ob, color):
    """ plot coors v3

    References
    ----------

    http://pypi.python.org/pypi/Shapely
    """
    x, y = ob.xy

    ax.plot(x, y, 'o', color=color, zorder=2)


def plot_bounds3(ax, ob, color):
    """ plot bounds v3

    References
    ----------
    http://pypi.python.org/pypi/Shapely
    """
    x, y = zip(*list((p.x, p.y) for p in ob.boundary))
    ax.plot(x, y, color=color, zorder=1)


def plot_line3(ax, ob, color):
    """ plot lines v3

    References
    ----------

    http://pypi.python.org/pypi/Shapely

    """
    x, y = ob.xy
    ax.plot(x, y, color=color, alpha=0.7, linewidth=2,
            solid_capstyle='round', zorder=1)
#
## wedge functions
##


def valid_wedge(ps, pw, p1, p2, grazing):
    """ check set of N wedge sector validity for point ps

     Parameters
     ----------

     ps       : source point
     pw       : np.array (Nx2) wedge apex point
     p1       : np.array (Nx2) point 1 of wedge
     p2       : np.array (Nx2) point 2 of wedge
     grazing  : 0 (without grazing)
                1 (authorize grazing)

                         xps

            x pw
           / \
          /   \
         /     \
        x p1    x p2

    Returns
    -------

    valid : np.array (Nx1)
            valid = 1 if ps is in the convex sector
            valid = 0 if ps is in the concav sector

    Examples
    --------

    >>> p1 = np.array([-2,-2]).reshape(1,2)
    >>> p2 = np.array([2,-2]).reshape(1,2)
    >>> pw = np.array([0,0]).reshape(1,2)
    >>> ps1 = np.array([3,3]).reshape(1,2)
    >>> ps2 = np.array([0,-3]).reshape(1,2)
    >>> valid_wedge(ps1,pw,p1,p2,0)[0][0]
    1.0
    >>> valid_wedge(ps2,pw,p1,p2,0)[0][0]
    1.0


    Authors
    -------
    Bernard.uguen@univ-rennes1.fr
    """

    x1 = p1[:, 0] - pw[:, 0]
    y1 = p1[:, 1] - pw[:, 1]

    a1 = np.arctan2(y1, x1)

    x2 = p2[:, 0] - pw[:, 0]
    y2 = p2[:, 1] - pw[:, 1]

    a2 = np.arctan2(y2, x2)

    xs = ps[:, 0] - pw[:, 0]
    ys = ps[:, 1] - pw[:, 1]

    aas = np.arctan2(ys, xs)

    valid = np.zeros((1, len(x1)))

    b1_I2 = min(a1.all(), a2.all())
    b2_I2 = max(a1.all(), a2.all())

    mu_I2 = b2_I2 - b1_I2
    #
    # un >= ou <= permet de valider les points qui sont sur une tangente au diedre
    #
    if (grazing == 0):
        in_I2 = np.nonzero(
            (aas > b1_I2) & (aas < b2_I2) & (mu_I2 > np.pi))[0]
        valid[in_I2] = 1
        out_I2 = np.nonzero(
            ((aas < b1_I2) | (aas > b2_I2)) & (mu_I2 < np.pi))[0]
        valid[out_I2] = 1

    if (grazing == 1):
        in_I2 = np.nonzero(
            (aas >= b1_I2) & (aas <= b2_I2) & (mu_I2 > np.pi))[0]
        valid[in_I2] = 1
        out_I2 = np.nonzero(
            ((aas <= b1_I2) | (aas >= b2_I2)) & (mu_I2 < np.pi))[0]
        valid[out_I2] = 1

    return(valid)

def agwed_old(v,lwe):
    """

    Parameters:
    -----------
    lwe : np.array
        3x1 wedge vector
    v   : np.array(3x4)
        3x4 ( 4 stacked vectors)

        first vector of v is on face 0 perp to lwe
        second vector of v is on face n perp to lwe
        third vector is on the direction of incident ray  (-si)
        fourth vector is on the direction of diffracted ray (sd)

    all vectors of v are defined outgoing from the diffracting point

    Returns
    -------
        np.array([N*pi,phi0,phi])


    Example
    -------

    >>> import numpy as np
    >>> import pylayers.util.geomutil as geu
    >>> lwe = np.array([0,0,1])
    >>> u = np.array([1,0,0])
    >>> v1 = np.array([1,1,0])
    >>> si = np.array([-1,-1,0])
    >>> sd = np.array([-1,1,0])
    >>> v  = np.vstack([u,v1,si,sd]).T
    >>> M = geu.agwed(v,lwe)
    >>> print M*180/np.pi
    [ 315.  135.  225.]

    """
    print DeprecationWarning('Please use vectorized version : agwed')
    # lwe : (,3)
    lwe = lwe/np.sqrt(np.sum(lwe*lwe,axis=0))
    # v : (3,4)
    v  = v/np.sqrt(np.sum(v*v,axis=0))
    # ps (,4)
    ps  = np.dot(lwe,v)
    vp1  = v - v*ps
    vpn = vp1/np.sqrt(np.sum(vp1*vp1,axis=0))
    vpt  = vpn[0:2,:].T
    w = np.vstack((vpt[:,1],-vpt[:,0]))
    C = np.dot(vpt,w)
    D = np.dot(vpt,vpt.T)
    M = np.mod(2*np.pi-np.arctan2(np.dot(vpt,w),np.dot(vpt,vpt.T)),2*np.pi)[0,1:]
    return M

def agwed(v,lwe):
    """

    Parameters:
    -----------
    lwe : np.array
        3xNp wedge vector
    v   : np.array(3x4xNp)
        3x4xNp ( 4 stacked vectors)

        first vector of v is on face 0 perp to lwe
        second vector of v is on face n perp to lwe
        third vector is on the direction of incident ray  (-si)
        fourth vector is on the direction of diffracted ray (sd)

    all vectors of v are defined outgoing from the diffracting point

    Returns
    -------
        np.array([[N*pi,phi0,phi],...xNp])
                (3xNp)


    Example
    -------

    >>> import pylayers.util.geomutil as geu
    >>> import numpy as np
    >>> lwe = np.array([[0,0,1],[0,0,1]]).T
    >>> u = np.array([[1,0,0],[1,0,0]]).T
    >>> v1 = np.array([[1,1,0],[1,1,0]]).T
    >>> si = np.array([[-1,-1,0],[-1,1,0]]).T
    >>> sd = np.array([[-1,1,0],[1,-1,0]]).T
    >>> v  = np.hstack((u[:,None,:],v1[:,None,:],si[:,None,:],sd[:,None,:]))
    >>> M = geu.agwed(v,lwe)
    >>> print M*180/np.pi
    array([[ 315.,  315.],
       [ 135.,  225.],
       [ 225.,   45.]])

    """
    import ipdb
    ipdb.set_trace()
    # lwe : (3,N)
    lwe = lwe/np.sqrt(np.sum(lwe*lwe,axis=0))
    # v : (3,4,N)
    v  = v/np.sqrt(np.sum(v*v,axis=0))
    # ps (4,N)
    #ps  = np.dot(lwe,v)
    ps = np.einsum('ik,ijk->jk',lwe,v)
    vp1  = v - v*ps
    vpn = vp1/np.sqrt(np.sum(vp1*vp1,axis=0))
    # vpt = (N,4,2)
    vpt  = vpn[0:2,:,:]

    # w(4,N,2)
    w = np.dstack((vpt[1,:,:].T,-vpt[0,:,:].T)).T
    # C = np.dot(vpt,w)
    # D = np.dot(vpt,vpt.T)

    #vpt(2,4,N) x w(2,4,N) => C(4,4,N)
    C = np.einsum('kil,kjl->ijl',vpt,w)
    # D(4,4,N)
    D = np.einsum('kil,kjl->ijl',vpt,vpt)

    #M = np.mod(2*np.pi-np.arctan2(np.dot(vpt,w),np.dot(vpt,vpt.T)),2*np.pi)[0,1:]
    M = np.mod(2*np.pi-np.arctan2(C,D),2*np.pi)[0,1:,:]
    return M

def sector(p1, p2, pt):
    """ non signed angular sector  between
        (p1,pt) and (p2,pt)

    p1 x-----------x pt
               |  /
          alpha \/
                /
               x p2

    Parameters
    ----------
    p1 : np.array (3 x Np)
        point
    p2 : np.array (3 x Np)
        point
    pt : np.array (3 x Np)
        point

    Returns
    -------

    alpha : np.array (3 x Np)
        degree

    Notes
    -----

    Useful for AAS calculation


    """

    if len(p1.shape) == 1:
        p1=p1.reshape(p1.shape[0],1)
    if len(p2.shape) == 1:
        p2=p2.reshape(p2.shape[0],1)
    if len(pt.shape) == 1:
        pt=pt.reshape(pt.shape[0],1)
    p1pt = p1 - pt
    p2pt = p2 - pt
    u = p1pt / np.sqrt(np.sum((p1pt)*(p1pt),axis=0))
    v = p2pt / np.sqrt(np.sum((p2pt)*(p2pt),axis=0))
    # sum(a[i,j,:] * b[k,:,m])
    alpha = np.arctan2(u[1], u[0])
    beta = np.arctan2(v[1], v[0])
    v0 = abs(alpha - beta)
    v1 = 2 * np.pi - abs(alpha - beta)
    um0 = v0 < v1
    um1 = ~um0
    sector = np.empty(np.shape(u)[1])
    sector[um0]= v0[um0]
    sector[um1]= v1[um1]
    return sector*180/np.pi
    #if (abs(alpha + sector - sp.mod(beta, 2 * np.pi)) < 1e-3):
    #    return(np.array([alpha, beta]) * 180 / np.pi)
    #else:
    #    return(np.array([beta, alpha]) * 180 / np.pi)

def sectorold(p1, p2, pt):
    """ angular sector  p1 pt p2

    Parameters
    ----------
    p1 : np.array
        point
    p2 : np.array
        point
    pt : np.array
        point

    Returns
    -------

    alpha : np.array
        degree

    Notes
    -----

    Useful for AAS calculation


    """
    u = (p1 - pt) / np.sqrt(np.dot(p1 - pt, p1 - pt))
    v = (p2 - pt) / np.sqrt(np.dot(p2 - pt, p2 - pt))
    alpha = np.arctan2(u[1], u[0])
    beta = np.arctan2(v[1], v[0])
    sector = min(abs(alpha - beta), 2 * np.pi - abs(alpha - beta))
    return sector*180/np.pi


def dist(x,y,ax):
        """ calculates distance between two arrays along a given axis

        Parameters
        ----------
            x : numpy.ndarray
            y : numpy.ndarray
            ax : integer (0,1)

        Returns
        -------
            d : numpy.ndarray

        Examples
        --------
        .. plot::
            :include-source:

            >>> import numpy as np
            >>> x = np.array([[0., 0., 10., 10.],[0., 10., 10., 0.]])
            >>> y = np.array([[5.],[5.]])
            >>> ax = 0
            >>> d = dist(x,y,ax)
        """
        d = np.sqrt(np.sum((x-y)**2, axis=ax))
        return d

def line_intersection(l1,l2):
    """ intersection between two 2D lines using shapely

    Parameters
    ----------

    l1: numpy.ndarray
        coordinates of l1 points
    l2: numpy.ndarray
        coordinates of l2 points

    Returns
    -------

    p: numpy.ndarray
        coordinates of intersection point

    """
    shl1 = sh.LineString((l1[:,0],l1[:,1]))
    shl2 = sh.LineString((l2[:,0],l2[:,1]))
    if shl1.intersects(shl2):
        psh = shl1.intersection(shl2)
        return np.array([[psh.x],[psh.y]])
    else:
        return None

def linepoly_intersection(l,poly):
    """ intersection between a 2D line and a 2D polygon using shapely

    Parameters
    ----------

    l: numpy.ndarray
        coordinates of l points
    poly: numpy.ndarray
        coordinates of poly points

    Returns
    -------

    p: numpy.ndarray
        coordinates of intersection point

    """
    shl = sh.LineString((l[:,0],l[:,1]))
    shpoly = sh.polygon((poly[:,0],poly[:,1],poly[:,2]))
    psh = shl.intersection(shpoly)
    return np.array([[psh.x],[psh.y]])

def mirror(p,pa,pb):
    """ Compute the image of p wrt the segment pa pb

    Parameters
    ----------

    p : numpy.ndarray
        point to image
    pa : numpy.ndarray
        segment tail
    pb : numpy.ndarray
        segment head

    Returns
    -------

    M : numpy.ndarray

    Examples
    --------

    .. plot::
        :include-source:

    >>> from pylayers.util.geomutil import *
    >>> from pylayers.util.plotutil import *
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> p = np.random.randn(2,10)
    >>> pa  = np.array([-0.5,1])
    >>> pb  = np.array([0,0])
    >>> M = mirror(p,pa,pb)
    >>> plt.plot(p[0,:],p[1,:],'or',alpha=0.2)
    >>> plt.plot(M[0,:],M[1,:],'ob',alpha=0.2)
    >>> displot(p,M,alpha=0.2)
    >>> axis = np.vstack((pa,pb))
    >>> plt.plot(axis[:,0],axis[:,1])

    """

    if np.shape(pa)==(2,):
        pa = pa.reshape(2,1)
    if np.shape(pb)==(2,):
        pb = pb.reshape(2,1)

    pab = pb - pa
    alpha = np.sum(pab * pab, axis=0)
    zalpha = np.where(alpha == 0.)
    alpha[zalpha] = 1.

    a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
    b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
    c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
                        pa[1, :] * (pa[1, :] - pb[1, :]) *
                        (pb[0, :] - pa[0, :]))
    d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
                        pa[0, :] * (pa[1, :] - pb[1, :]) *
                        (pb[0, :] - pa[0, :]))

    N = 1
    S = np.zeros((2, 2))
    S[0, 0] = -a
    S[0, 1] = b
    S[1, 0] = b
    S[1, 1] = a
    A = np.eye(2)
    y = np.zeros(2)
    vc0 = np.array([c[0], d[0]]).reshape(2,1)
    v0 = np.dot(-S, p) + vc0
    x = la.solve(A, v0)
    return x

def axmat(pa,pb):
    """ Compute the image of p wrt the segment pa pb

    Parameters
    ----------


    pa : numpy.ndarray
        segment tail
    pb : numpy.ndarray
        segment head

    Returns
    -------

    S : numpy.ndarray
        symmetry matrix
    v : numpy.ndarray
        translatrion vector

    Notes
    -----

    fir x is corrdiante of the point to mirror, 
    the mirrored point x' from pa and pb  can be obtain with :

    x' = np.dot(x,S) + v



    Examples
    --------

    .. plot::
        :include-source:

    >>> from pylayers.util.geomutil import *
    >>> from pylayers.util.plotutil import *
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> p = np.random.randn(2,10)
    >>> pa  = np.array([-0.5,1])
    >>> pb  = np.array([0,0])
    >>> S,v = axmat(pa,pb)
    >>> M = np.dot(p,S) + v
    >>> plt.plot(p[0,:],p[1,:],'or',alpha=0.2)
    >>> plt.plot(M[0,:],M[1,:],'ob',alpha=0.2)
    >>> displot(p,M,alpha=0.2)
    >>> axis = np.vstack((pa,pb))
    >>> plt.plot(axis[:,0],axis[:,1])

    """

    if np.shape(pa)==(2,):
        pa = pa.reshape(2,1)
    if np.shape(pb)==(2,):
        pb = pb.reshape(2,1)

    pab = pb - pa
    alpha = np.sum(pab * pab, axis=0)
    zalpha = np.where(alpha == 0.)
    alpha[zalpha] = 1.

    dsal = (2. / alpha)
    pampby = pa[1, :] - pb[1, :]
    pbmpax = pb[0, :] - pa[0, :]
    prod = pbmpax * pampby

    a = 1 - dsal * (pampby ** 2)
    b = dsal * prod
    c = dsal * (pa[0, :] * (pampby ** 2) + pa[1, :] * prod)                   
    d = dsal * (pa[1, :] * (pbmpax ** 2) + pa[0, :] * prod) 
                        
    # a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
    # b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
    # c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
    #                     pa[1, :] * (pa[1, :] - pb[1, :]) *
    #                     (pb[0, :] - pa[0, :]))
    # d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
    #                     pa[0, :] * (pa[1, :] - pb[1, :]) *
    #                     (pb[0, :] - pa[0, :]))

    N = 1
    S = np.array([[a[0],-b[0]],[-b[0],-a[0]]])
    vc0 = np.array([c[0], d[0]])
    # v0 = np.dot(-S, p) + vc0
    return S,vc0


def distseg(a,b,c,d,alpha,beta):
    """ distance to segments

    Parameters
    ----------

    a : (3xN) initial point segment 1
    b : (3xN) end point segment 1
    c : (3xN) starting point segment 2
    d : (3xN) end point segment 2

    alpha :
    beta  :

    Returns
    -------

    f : square of the distance to the segment

    Examples
    --------

    >>> import numpy as np
    >>> np.random.seed(0)
    >>> a = np.random.rand(3,10)
    >>> b = np.random.rand(3,10)
    >>> c = np.random.rand(3,10)
    >>> d = np.random.rand(3,10)
    >>> alpha,beta,dmin = dmin3d(a,b,c,d)
    >>> alpha[alpha<0]=0
    >>> alpha[alpha>1]=1
    >>> beta[beta<0]=0
    >>> beta[beta>1]=1
    >>> f = distseg(a,b,c,d,alpha,beta)
    >>> p1 = a - alpha*(a-b)
    >>> p2 = c + beta*(d-c)
    >>> v = p1-p2
    >>> g = np.sum(v*v,axis=0)
    >>> diff = np.sum(f-g,axis=0)
    >>> np.testing.assert_almost_equal(diff,0)

    """

    ac = c-a
    cd = d-c
    ba = a-b

    u0 = np.sum(ac*ac,axis=0)
    u4 = np.sum(ba*ba,axis=0)
    u5 = np.sum(cd*cd,axis=0)
    u1 = np.sum(ba*ac,axis=0)
    u2 = np.sum(cd*ac,axis=0)
    u3 = np.sum(cd*ba,axis=0)

    f = u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5

    # m = a - alpha*ba
    # n = c + beta*cd
    # g = np.dot(m-n,m-n)

    return f

def dmin3d(a,b,c,d):
    """ evaluate the minimal distance between 2 set of segments

    Parameters
    ----------

    a : (3xN) initial point segment 1
    b : (3xN) end point segment 1
    c : (3xN) starting point segment 2
    d : (3xN) end point segment 2

    Returns
    -------

    alpha : segment parameterization
    beta  : segment parameterization
    dmin  : minimal distance between 2 segments

    Examples
    --------

    """

    ac = c-a
    cd = d-c
    ba = a-b

    u0 = np.sum(ac*ac,axis=0)
    u4 = np.sum(ba*ba,axis=0)
    u5 = np.sum(cd*cd,axis=0)
    u1 = np.sum(ba*ac,axis=0)
    u2 = np.sum(cd*ac,axis=0)
    u3 = np.sum(cd*ba,axis=0)

    den = u4*u5-u3*u3
    alpha = (u2*u3-u1*u5)/(1.*den)
    beta = (u1*u3-u2*u4)/(1.*den)
    dmin = np.sqrt(u0 + 2*(alpha*u1+beta*u2+alpha*beta*u3)+alpha*alpha*u4+ beta*beta*u5)

    return(alpha,beta,dmin)

# def gram_schmid(V): 
#     """ 
#     Gram-Schmid orthonormalization of a set of `M` vectors, in-place. 

#     Parameters 
#     ---------- 
#     V : array, shape (N, M) 

#     Notes
#     -----

#     from http://numpy-discussion.10968.n7.nabble.com/Efficient-orthogonalisation-with-scipy-numpy-td23635.html

#     """ 
#     # XXX: speed can be improved by using routines from scipy.lib.blas 
#     # XXX: maybe there's an orthonormalization routine in LAPACK, too, 
#     #      apart from QR. too lazy to check... 
#     n = V.shape[1] 
#     for k in xrange(n): 
#         V[:,k] /= np.linalg.norm(V[:,k]) 
#         for j in xrange(k+1, n): 
#             V[:,j] -= np.vdot(V[:,j], V[:,k]) * V[:,k] 
#     return V 


def gram_schmidt(Vini,force_direct=True):
    """
    Gram-Schmidt orthonormalization of a set of `M` vectors, in-place. 

    Parameters
    ----------

     Vini : array,
        shape (3,Nv,nf)  where number of vectors Nv = 3 and nf is an integer
    force_direct : boolean
        force basis to be direct (det>0)

    Example
    -------

    >>> import pylayers.util.geomutil as geu
    >>> import numpy as np
    >>> Nv = 3
    >>> Nframes = 10
    >>> V = np.random.rand(3,Nv,Nframes)
    >>> VG = geu.gram_schmid(V)
    """


    # check direct basis
    if force_direct:
        per=permutations((0,1,2),3)
        for p in per:
            P = np.vstack((Vini[:,p[0],0],Vini[:,p[1],0],Vini[:,p[2],0]))
            if np.linalg.det(P) >0 :
                Vini = Vini[:,p,:]
                break

    v0=Vini[:,0,:]
    v1=Vini[:,1,:]
    v2=Vini[:,2,:]

    n0 = np.linalg.norm(v0,axis=0)
    vn0=v0/n0

    pv10 = np.einsum('ij,ij->j',v1,vn0)
    v1p = v1 - pv10*vn0
    nv1 =  np.linalg.norm(v1p,axis=0)
    vn1 = v1p/nv1


    pv20 = np.einsum('ij,ij->j',v2,vn0)
    pv21 = np.einsum('ij,ij->j',v2,vn1)

    v2p = v2 - pv20*vn0 - pv21*vn1
    nv2 =  np.linalg.norm(v2p,axis=0)
    vn2 = v2p/nv2

    V= np.hstack((vn0[:,None,:],vn1[:,None,:],vn2[:,None,:]))

    if force_direct:
        # assert det >0
        assert len(np.where(np.linalg.det(np.rollaxis(V,2))<0)[0]) ==0
    # assert det != 0
    assert len(np.where(np.linalg.det(np.rollaxis(V,2))==0.)[0]) ==0
    return V



def qrdecomp(V):
    """
    Gram-Schmid orthonormalization of a set of `Nv` vectors, in-place.
    using qr decomp

    Parameters
    ----------

    V : array,
        shape (3,Nv,nf)  where number of vectors Nv = 3 and nf is an integer

    Returns
    -------

    V : array,


    References
    ----------

    from http://numpy-discussion.10968.n7.nabble.com/Efficient-orthogonalisation-with-scipy-numpy-td23635.html

    Example
    -------
    >>> import numpy as np
    >>> import pylayers.util.geomutil as geu
    >>> u=np.random.rand(3,1,10)
    >>> v=np.random.rand(3,1,10)
    >>> w=np.random.rand(3,1,10)
    >>> V = np.hstack((u,v,w))
    >>> W = geu.qrdecomp(V)
    >>> assert np.allclose(abs(np.linalg.det(W[:,:,0])),1.0)

    """
    #  speed can be improved by using routines from scipy.lib.blas
    #  maybe there's an orthonormalization routine in LAPACK, too,
    #      apart from QR. too lazy to check...

    import copy
    # nn = np.linalg.norm(V,axis=(0))

    # # for i in range(3):
    # #     V[i,:,:]=V[i,:,:]/nn

    # V=V/nn
    lv = np.shape(V)[2]
    V2 = copy.deepcopy(V)
    for k in xrange(lv):
        V[:,:,k],R = np.linalg.qr(V[:,:,k])
    # check where the vector along cylinder axis is colinear with the 1st basis axis
    # col = np.einsum('ij,ij->j',V[:,0,:],V2[:,0,:])
    # ucol = np.where(col < 0)
    # import ipdb
    # ipdb.set_trace()
    # V[:,:,ucol]=-V[:,:,ucol]
    import ipdb
    ipdb.set_trace()
    return V




def check_point_unicity(A):
    """ check if all rows of an array are unique

    Parameters
    ----------

        A : np.ndarray (Npt, 2|3)

    Return
    ------
        similar : list
        list of index of similar points
        if void list, all poitns are differents

    Example
    -------

    >>> import numpy as np
    >>> a = np.arange(10)
    >>> a = np.np.vstack((a,a))
    >>> check_point_unicity(a.T)
    []
    >>> b=np.array([4,4])
    >>> aa=np.concatenate((a,b[:,None]),axis=1)
    >>> check_point_unicity(aa.T)
    [4, 10]

    """
    similar = []
    for ua in xrange(len(A)):
        rA=np.roll(A,-ua,axis=0)
        # print rA
        if any((A[ua] == x).all() for x in rA[1:]):
            similar.append(ua)
    return similar


if __name__ == "__main__":
    plt.ion()
    doctest.testmod()


