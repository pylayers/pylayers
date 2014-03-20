#-*- coding:Utf-8 -*-
"""

Utility Functions
=================

.. autosummary::
    :toctree: generated/

    savestr2
    stretch
    segsplit
    extract
    inbracket
    incrochet
    geomLine
    geomFace
    ParseDirectionalLight
    ParseMaterial
    show
    parsevrml
    vrml2sha
    vrml2geom

VLayout Class
=============

.. autosummary::
    :toctree: generated/

    Vlayout.load
    Vlayout.show
    Vlayout.wallanalysis
    Vlayout.show3entity

"""
import os
import doctest
import glob
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylayers.util import geomutil as geo
from pylayers.util import pyutil as pyu
from descartes.patch import PolygonPatch
import shapely.geometry as shg
import shapely.ops as sho
import networkx as nx
from pylayers.gis.layout import Layout


def savestr2(dpt,dwall, _filename='struc.str2'):
    """ save walls in str2 format

    The default filename is struc.str2

    Parameters
    ----------

    dpt :
    dwall:

    """
    fd = open('struc.str2', 'w')
    s1 = str(len(dpt.keys()))+' '+str(len(dwall.keys()))+' 0\n'
    fd.write(s1)
    for ipt in dpt:
        p = dpt[ipt]
        s1 = str(p[0])+' '+str(p[1])+' 0 0 0 0\n'
        fd.write(s1)
    for iw in dwall:
        tail = dwall[iw]['tail']+1
        head = dwall[iw]['head']+1
        zmin = dwall[iw]['zmin']
        zmax = dwall[iw]['zmax']
        core = str(2)
        s1 = str(tail)+' '+str(
            head)+' 1 '+core+' '+'1 0 '+str(zmin)+' '+str(zmax)+'\n'
        fd.write(s1)
    fd.close()


def stretch(s, alphat=0.1, alphah=0.1):
    """ strech a LineString by a given perc on both terminations

        Parameters
        ----------
        s      : LineString shapely
        alphat : stretching coeff tail
        alphah : stretching coeff head

        Returns
        -------
        ss    : streched segment

        Examples
        --------

        >>> s1  = shg.LineString(((0,0),(1,0)))
        >>> ss1 = stretch(s1,0.1)
        >>> s2 = shg.LineString(((-0.1,0),(1.1,0)))
        >>> assert (ss1.equals(s2))
    """
    ls = s.length
    x, y = s.xy
    u = np.array((x[1]-x[0], y[1]-y[0]))
    un = u/ls
    pt = np.array((x[0], y[0]))
    ph = np.array((x[1], y[1]))
    # ppt = pt - un*ls*alpha
    # pph = ph + un*ls*alpha
    ppt = pt - un*alphat
    pph = ph + un*alphah
    ss = shg.LineString(((ppt[0], ppt[1]), (pph[0], pph[1])))
    return(ss)


def segsplit(s1, s2, tol=0.0001, alpha=0.1):
    """ split segment

    Parameters
    ----------

    s1    : shapely LineString
    s2    : shapely LineString
    tol   : tolerance for point equality test
    alpha : stretching factor

    Returns
    -------
    ts   : list of segment
    bks1 : boolean keep s1
    bks2 : boolean keep s2


    Examples
    --------

    >>> s1 = shg.LineString(((0,0),(1,0)))
    >>> s2 = shg.LineString(((1,0),(2,0)))
    >>> s3 = shg.LineString(((1,-10),(1,10)))
    >>> s4 = shg.LineString(((0.5,-10),(0.5,10)))
    >>> ts1 = segsplit(s1,s2)
    >>> ts2 = segsplit(s1,s3)
    >>> ts3 = segsplit(s1,s4)

    """
    ts1 = []
    ts2 = []
    bks1 = True
    bks2 = True
    beta = alpha/(2*alpha+1)
    if s1.intersects(s2):
        p1t, p1h = s1.boundary
        p2t, p2h = s2.boundary
        ls1 = s1.length
        ls2 = s2.length
        pi = s1.intersection(s2)
        if s1.touches(s2):  # touching segments
            if not (pi.equals_exact(p1t, tol) or pi.equals_exact(p1h, tol)):
                s11 = shg.LineString(((p1t.xy[0][0], p1t.xy[
                                     1][0]), (pi.xy[0][0], pi.xy[1][0])))
                s12 = shg.LineString(((pi.xy[0][0], pi.xy[
                                     1][0]), (p1h.xy[0][0], p1h.xy[1][0])))
                if s11.length > 0 and s11.length >= alpha:
                    ts1.append(s11)
                if s12.length > 0 and s12.length >= alpha:
                    ts1.append(s12)
                bks1 = False
            if not (pi.equals_exact(p2t, tol) or pi.equals_exact(p2t, tol)):
                s21 = shg.LineString(((p2t.xy[0][0], p2t.xy[
                                     1][0]), (pi.xy[0][0], pi.xy[1][0])))
                s22 = shg.LineString(((pi.xy[0][0], pi.xy[
                                     1][0]), (p2h.xy[0][0], p2h.xy[1][0])))
                if s21.length > 0 and s21.length > alpha:
                    ts2.append(s21)
                if s22.length > 0 and s22.length >= alpha:
                    ts2.append(s22)
                bks2 = False
        else:  # crossing segments
            s11 = shg.LineString(((p1t.xy[0][0], p1t.xy[
                                 1][0]), (pi.xy[0][0], pi.xy[1][0])))
            s12 = shg.LineString(((pi.xy[0][0], pi.xy[
                                 1][0]), (p1h.xy[0][0], p1h.xy[1][0])))
            s21 = shg.LineString(((p2t.xy[0][0], p2t.xy[
                                 1][0]), (pi.xy[0][0], pi.xy[1][0])))
            s22 = shg.LineString(((pi.xy[0][0], pi.xy[
                                 1][0]), (p2h.xy[0][0], p2h.xy[1][0])))
            ls11 = s11.length
            ls12 = s12.length
            ls21 = s21.length
            ls22 = s22.length

            if ls11 > ls12:
                ts1.append(s11)
            else:
                ts1.append(s12)

            if ls21 > ls22:
                ts2.append(s21)
            else:
                ts2.append(s22)
            # if s11.length>0 and s11.length>=alpha:
            #    ts1.append(s11)
            # if s12.length>0 and s12.length>=alpha:
            #    ts1.append(s12)
            # if s21.length>0 and s21.length>=alpha:
            #    ts2.append(s21)
            # if s22.length>0 and s21.length>=alpha:
            #    ts2.append(s22)
            bks1 = False
            bks2 = False
    return(ts1, ts2, bks1, bks2)


def extract(vrmlstrg, dico):
    """ converts recursively a vrml string into a dictionnary

     Parameters
     ----------

     vrmlstrg:
     dico :

     Returns
     ------
     dico : dictonnary associated with vrml string strg
    """
    val = vrmlstrg
    while len(val) != 0:
        key, val = inbracket(val)
        if val != '':
            dico[key] = val
            dico = extract(val, dico)

    return(dico)


def inbracket(strg):
    """ extraction of bracket content

    Parameters
    ----------
    strg : a string with a bracket

    Returns
    -------
    lbra : left part of the string
    inbr : string inside the bracket

    Examples
    --------

    >>> strg ='abcd{un texte}'
    >>> lbra,inbr = inbracket(strg)
    >>> assert(lbra=='abcd')
    >>> assert(inbr=='un texte')
    """
    strg = strg.replace('\r', '')
    ssp = strg.split('{')
    lbra = ssp[0]
    rbra = ''
    inbr = ''
    for k in ssp[1:]:
        rbra = rbra+k+'{'
    rsp = rbra.split('}')
    for k in rsp[:-1]:
        inbr = inbr+k+'}'
    inbr = inbr.rstrip('}')
    return(lbra, inbr)


def incrochet(strg):
    """ get content inside crochet
    Parameters
    ----------
    strg : string

    Returns
    -------
    lbra : left part of the string
    inbr : string inside the bracket

    Examples
    --------

    >>> strg ='abcd[un texte]'
    >>> lbra,inbr = incrochet(strg)
    >>> assert(lbra=='abcd')
    >>> assert(inbr=='un texte')

    """
    strg = strg.replace('\r', '')
    ssp = strg.split('[')
    lbra = ssp[0]
    rbra = ''
    inbr = ''
    for k in ssp[1:]:
        rbra = rbra+k+'['
    rsp = rbra.split(']')
    for k in rsp[:-1]:
        inbr = inbr+k+']'
    inbr = inbr.rstrip(']')
    return(lbra, inbr)


def geomLine(st):
    """ build a Line from string

    Parameters
    ----------
    st : string

    Returns
    -------
    tabindex   : array of indexes
    tcoord     : array of coordinates

    """
    st1 = st.split('coordIndex')
    index = st1[1]
    index = index.replace('[', '')
    index = index.replace(']', '')
    index = index.replace(' ', '')
    tindex = index.split(',')
    tabindex = []
    for st in tindex[:-1]:
        tabindex.append(int(st))
    tabindex = np.array(tabindex).reshape(len(tabindex)/3, 3)
    a, b = inbracket(st1[0])
    c = incrochet(b)[1]
    c = c.split(',')

    coord = []
    for st in c[:-1]:
        pt = st.split(' ')
        for ic in pt:
            coord.append(float(ic))
    tcoord = np.array(coord).reshape(len(coord)/3, 3)
    return(tabindex, tcoord)


def geomFace(st):
    """ build a Face from string

    Parameters
    ----------
    st : string

    Returns
    -------

    tabindex
    tcoord : ndarray

    """
    st1 = st.split('coordIndex')
    index = st1[1]
    index = index.replace('[', '')
    index = index.replace(']', '')
    index = index.replace(' ', '')
    tindex = index.split(',')
    tabindex = []
    for st in tindex[:-1]:
        tabindex.append(int(st))
    tabindex = np.array(tabindex)
    a, b = inbracket(st1[0])
    c = incrochet(b)[1]
    c = c.split(',')

    coord = []
    for st in c[:-1]:
        pt = st.split(' ')
        for ic in pt:
            coord.append(float(ic))
    tcoord = np.array(coord).reshape(len(coord)/3, 3)
    return(tabindex, tcoord)


def ParseDirectionalLight(st):
    """

    """
    d = {}
    s1 = st.split('intensity')
    t = s1[0]
    st1 = t.split(' ')
    d['on'] = bool(st1[1])
    t = s1[1]
    s2 = t.split('ambientIntensity')
    d['intensity'] = float(s2[0])
    t = s2[1]
    s2 = t.split('color')
    d['ambientIntensity'] = float(s2[0])
    t = s2[1]
    s2 = t.split('direction')
    d['color'] = s2[0]
    d['direction'] = s2[1]
    return(d)


def ParseMaterial(st):
    """

    """
    st = st.replace('material', '')
    st = st.replace('Material', '')
    st = st.replace('{', '')
    dst = st.split(',')
    d = {}
    for s in dst:
        # print s
        u = s.split(' ')
        # print u
    sp1 = st.split('diffuseColor')
    t = sp1[1]
    ts = t.split(',')
    print ts
    d['diffuseColor'] = ts[0]
    try:
        d['specularColor'] = ts[1].split('specularColor')[1]
        d['shininess'] = ts[2].split('shininess')[1]
        d['transparency'] = ts[3].split('transparency')[1]
    except:
        d['specularColor'] = ts[2].split('specularColor')[1]
        d['shininess'] = ts[3].split('shininess')[1]
        d['transparency'] = ts[4].split('transparency')[1]
    return(d)


def show(dico):
    """ show dico
    """
    for key in dico.keys():
        if key != 'name':
            plt.plot(dico[key]['coord'][:, 0], dico[key]['coord'][:, 2])
        else:
            print dico[key]
    plt.show()


def parsevrml(filename):
    """ parse a vrml file

    Parameters
    ----------

    filename : vrml filename

    Returns
    -------

    dg : dictionnaries of group

    """
    fd = open(filename, 'r')
    lignes = fd.readlines()
    fd.close()
    tobj = {}
    tnum = {}
    k = -1
    for li in lignes:
        li = li.replace('\n', '')
        li = li.replace('\r,', '')  # line finishing with a comma
        if li.find('DEF') != -1:      # if line contains DEF keyword new dictionnary entry
            k = k+1
            tmp = li.split(
                ' ')    # split with space the index 1 is the dictionnary entry name
            obid = tmp[1]
            tobj[k] = li
            tnum[k] = obid            # name of the entry
        else:
            try:
                tobj[k] = tobj[
                    k]+li  # concatenation of the line in the current dictionnary entry k
            except:
                pass
    td = {}
    dg = {}                         # list of groups
    targetkey = ' geometry IndexedLineSet '
    targetkey = ' geometry IndexedFaceSet '
    for on in tnum.keys():
        dico = {}
        ob = tnum[on]
        val = tobj[on]
        dico = extract(val, dico)
        if ob.find('LIGHT') != -1:
            td[ob] = ParseDirectionalLight(dico.values()[0])
        if ob.find('APPEARANCE') != -1:
            td[ob] = ParseMaterial(dico.values()[0])
        if val.find('Group') != -1:
            try:
                # tg.append(g)    # save previous group
                dg[name] = g
            except:
                pass
            g = {}
            name = val.split('Group')[0].split(' ')[1]  # name of the group
        if ob.find('ID_') != -1:
            st = dico[targetkey]
            i, c = geomFace(st)
            g[tnum[on]] = {'index': i, 'coord': c}

    dg[name] = g
    return(dg)


def vrml2sha(tg):
    """ convert vrml object into shapely polygons

        Parameters
        ----------
        tg   : list of objects
    """
    for l in tg:
        for k in l.keys():
            if k != 'name':
                c = l[k]['coord']
                i = l[k]['index']
                tt = []
                ltt = []
                for u in i:
                    if u == -1:  # -1 closure indicator
                        ltt.append(tt)
                        tt = []
                    else:
                        tt.append(u)
                P = geo.Polygon(c)


def vrml2geom(tg, rac):
    """ convert vrml object into geomview files

    Parameters
    ----------
    tg   : list of objects
    rac  : filename prefix
    """
    for l in tg:
        # prepare geomview file
        fi1 = rac+l['name']
        _fi1 = fi1+'.list'
        fina = pyu.getlong(_fi1, 'geom')
        fd = open(fina, 'w')
        fd.write('LIST')
        for k in l.keys():
            if k != 'name':
                filename = fi1+'-'+str(k)
                fd.write('{<'+filename+'.off}\n')
                G = geo.Geomoff(filename)
                c = l[k]['coord']
                i = l[k]['index']
                tt = []
                ltt = []
                for u in i:
                    if u == -1:
                        ltt.append(tt)
                        tt = []
                    else:
                        tt.append(u)
                # build a geomview list of polygons
                G.polygons(c, ltt)
        fd.close()


class VLayout(object):
    def load(self, filename):
        """
        Parameters
        ----------

        filename : str

        """
        dg = parsevrml(filename)
        self.entity = {}
        for t in dg:     # WALL, COLUMN, DOOR , STAIR , SPACE
            self.entity[t] = {}
            k = 0
            for ID in dg[t]:
                c = dg[t][ID]['coord']
                self.entity[t][k] = {}
                self.entity[t][k]['ID'] = ID
                self.entity[t][k]['coord'] = c
                l = dg[t][ID]['index']
                dp = {}
                p = []
                kk = 0
                for il in l:
                    if il == -1:
                        dp[kk] = p
                        p = []
                        kk = kk + 1
                    else:
                        p.append(il)

                self.entity[t][k]['index'] = dp
                k = k + 1
        #
        # Simplify Coord (Projection in 0xy plane)
        #
        for g in self.entity:
            for ID in self.entity[g]:
                x = self.entity[g][ID]['coord'][:, 0]
                y = -self.entity[g][ID]['coord'][:, 2]
                z = self.entity[g][ID]['coord'][:, 1]
                tp = np.vstack((x, y))
                Np = np.shape(tp)[1]
                tp2 = {}
                already = False
                iop = 0
                for ip in range(Np):
                    p = tp[:, ip]
                    for iold in tp2.keys():
                        pold = tp2[iold]['coord']
                        if np.shape(pold) == (2,):
                            dist = np.dot(p-pold, p-pold)
                            if dist < 1e-15:
                                already = True
                                tp2[iold]['nump'].append(ip)
                    if not already:
                        tp2[iop] = {}
                        tp2[iop]['coord'] = p
                        tp2[iop]['nump'] = [ip]
                        iop = iop + 1
                    already = False
                self.entity[g][ID]['c2d'] = tp2

        #
        # create a transcode between 3d index point and 2d index point
        #
        for g in self.entity:
            for ID in self.entity[g]:
                tr = {}
                dr2 = self.entity[g][ID]['c2d']
                for k in dr2.keys():
                    for u in dr2[k]['nump']:
                        tr[u] = k
                self.entity[g][ID]['tr'] = tr
        #
        # Create a new index for 2d points
        #
        for g in self.entity:
            for ID in self.entity[g]:
                di = self.entity[g][ID]['index']
                di2 = {}
                for k in di.keys():  # for all polygons
                    ti2 = []      # reserve a list= for 2d indexes
                    lpoly = di[k]  # get sequence of 3d points
                    for ip in lpoly:
                        ti2.append(self.entity[g][ID]['tr'][ip])
                    if len(np.unique(ti2)) == len(ti2):
                        di2[k] = ti2
                self.entity[g][ID]['index2'] = di2

        #
        # Create Polygon2D
        #

        for g in self.entity:
            for ID in self.entity[g]:
                dp = {}
                for ip in self.entity[g][ID]['index2']:
                    lp = self.entity[g][ID]['index2'][ip]
                    tp = []
                    for ip in lp:
                        tp.append(self.entity[g][ID]['c2d'][ip]['coord'])
                    poly = geo.Polygon(tp)
                    dp[ip] = poly
                self.entity[g][ID]['poly'] = dp
    #
    #
    #

    def show(self, num=100):
        """ show entities

        Parameters
        ----------
        num : int

        """
        if num > len(self.entity.keys()):
            group = self.entity.keys()
        else:
            group = [self.entity.keys()[num]]
        for g in group:
            for ID in self.entity[g]:
                for p in self.entity[g][ID]['poly']:
                    colrand = hex(int((2**23)*sp.rand(
                        1))+2**20).replace('0x', '#')
                    self.entity[g][ID]['poly'][
                        p].plot(color=colrand, alpha=0.3)
        plt.axis('scaled')

    def wallanalysis(self):
        """ walls analysis
        """
        w = self.entity['WALL']
        dwall = {}
        for k in w.keys():
            dwall[k] = {}
            dwall[k]['ID'] = w[k]['ID']
            #
            # height retrieval
            #
            z = np.sort(np.unique(w[k]['coord'][:, 1]))
            dwall[k]['zmin'] = min(z)
            dwall[k]['zmax'] = max(z)
            if len(z) == 2:
                dwall[k]['door'] = False
            else:
                dwall[k]['door'] = True
                dwall[k]['zdoor'] = z[1]
            #
            # Length,width retrieval
            #
            dp = w[k]['poly']
            gxmin = 1e15
            gymin = 1e15
            gxmax = -1e15
            gymax = -1e15
            for ik in dp.keys():
                pol = dp[ik]
                bounds = pol.bounds
                xmin = bounds[0]
                ymin = bounds[1]
                xmax = bounds[2]
                ymax = bounds[3]
                gxmin = min(gxmin, xmin)
                gymin = min(gymin, ymin)
                gxmax = max(gxmax, xmax)
                gymax = max(gymax, ymax)
            dx = gxmax-gxmin
            dy = gymax-gymin
            length = max(dx, dy)
            width = min(dx, dy)
            xmoy = (gxmin+gxmax)/2.
            ymoy = (gymin+gymax)/2.
            if (dx == width):
                seg = shg.LineString(((xmoy, gymin), (xmoy, gymax)))
            elif (dy == width):
                seg = shg.LineString(((gxmin, ymoy), (gxmax, ymoy)))
            # seg = stretch(seg,0.1)
            dwall[k]['bounds'] = bounds
            dwall[k]['thickness'] = width
            dwall[k]['length'] = length
            dwall[k]['seg'] = seg

        return(dwall)

    def show3entity(self, group, IDs):
        """ geomview vizualisation of entity

        Parameters
        ----------

        group :
        IDs :

        """
        te = self.entity[group]
        fi1 = 'entity'
        GL = geo.Geomlist(fi1)
        for k in IDs:
            ID = te.keys()[k]
            filename = fi1+'-'+str(k)
            GL.append('{<'+filename+'.off}\n')
            G = geo.Geomoff(filename)
            c = te[ID]['coord']
            i = te[ID]['index']
            tt = []
            ltt = []
            print i
            for u in i:
                ltt.append(i[u])
            # build a geomview list of polygons
            print ltt
            G.polygons(c, ltt)
        GL.show3()

if __name__ == "__main__":
    doctest.testmod()
