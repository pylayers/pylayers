# -*- coding:Utf-8 -*-
""" Handle Furnitures


Summary
-------

This module handle the description of furnitures
For now furnitures are rectangular objects.
In the future it should be able to accept any polygon as a furniture element

Furniture are stored in an ini file of dirstruc.

    from pylayers.gis.layout import *
    import matplotlib.pylab as plt
    L = Layout()
    L.loadstr('where1.ini')
    L.showGs()
    L.loadfur('Furw1.ini')
    fig = plt.gcf()
    ax  = fig.get_axes()[0]
    for fur in L.lfur:
        fur.show(fig,ax)


"""

import numpy as np
import re
import doctest
import ConfigParser
#import Graph
from   matplotlib.path import Path
from   matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
import pylayers.util.easygui
from   pylayers.util import pyutil as pyu
from   pylayers.util.project import *


class Furniture(object):
    """ Class Furniture

    Attributes
    ----------

    name
    desc
    origin
    height
    length
    width
    thichness
    Matname

    Notes
    -----

    A piece of furniture is a rectangular object with the following parameters

    """
    def __init__(self, name='',
                 desc='',
                 origin=np.array([0, 0]),
                 angle=0,
                 H=1.0,
                 L=1.0,
                 W=1,
                 T=0,
                 Matname='METAL'):
        """
            Parameters
            ----------
            desc    : furniture reference name
            origin
            angle
            H
            L
            W
            T
            Matname
        """
        self.name = name
        self.desc = desc
        self.origin = origin
        self.angle = angle
        self.height = H
        self.length = L
        self.width = W
        self.thickness = T
        self.Matname = Matname

    def info(self):
        """ display furniture information

        Provides info about Furniture object

        """
        for key in self.__dict__:
            print key, '\t\t:\t', self.__dict__[key]

    def set_position(self, p=[], angle=0):
        """ set position

        Parameters
        ----------

        p : ndarray (1x2)

        """
        if p == []:
            data = multenterbox('Choose position and angle ', 'Enter',
                                ('x', 'y', 'angle(deg)'),
                                (self.origin[0], self.origin[1], self.angle))
            self.origin[0] = eval(data[0])
            self.origin[1] = eval(data[1])
            self.angle = eval(data[2])
        else:
            self.origin[0] = p[0]
            self.origin[1] = p[1]
            self.angle = angle

    def load(self, _filename, name):
        """ load furniture from file

        Parameters
        ----------
        _filename
            file of furnitures
        name
            furniture id

        Examples
        --------
        .. plot::
            :include-source:

            >>> from pylayers.gis.furniture import *
            >>> import matplotlib.pylab as plt
            >>> F = Furniture()
            >>> F.load('Furw1.ini','R1_A')
            >>> fig,ax = F.show()
            >>> axis = plt.axis('scaled')
            >>> plt.show()
        """
        filefurn = pyu.getlong(_filename, "struc/furnitures")
        config = ConfigParser.ConfigParser()
        config.read(filefurn)

        self.name = config.get(name, "name")
        self.desc = config.get(name, "desc")

        ch = config.get(name, "origin")
        ch = ch.replace('[', '')
        ch = ch.replace(']', '')
        ch = ch.split(',')
        self.origin = np.array([eval(ch[0]), eval(ch[1])])

        self.angle = config.getfloat(name, "angle")
        self.height = config.getfloat(name, "height")
        self.width = config.getfloat(name, "width")
        self.length = config.getfloat(name, "length")
        self.thickness = config.getfloat(name, "thickness")
        self.Matname = config.get(name, "Matname")

    def save(self, _filename):
        """ save furniture file

        Parameters
        ----------

        _filename : string

        Furniture files are stored in  the struc directory of the current
        project

        .. todo::

            save a furniture file (struc)

        """
        filefurn = pyu.getlong(_filename, "struc")
        fd = open(filefurn, "a")
        config = ConfigParser.ConfigParser()

        secname = self.name
        config.add_section(secname)

        config.set(secname, "name", self.name)
        config.set(secname, "desc", self.desc)
        sorigin = '[ ' + str(self.origin[0]) + ' , ' + str(
            self.origin[1]) + ' ]'
        config.set(secname, "origin", sorigin)
        config.set(secname, "angle", self.angle)
        config.set(secname, "height", self.height)
        config.set(secname, "length", self.length)
        config.set(secname, "width", self.width)
        config.set(secname, "thickness", self.thickness)
        config.set(secname, "Matname", self.Matname)

        config.write(fd)

    def position(self):
        """ Calculate the coordinate of the furniture

        Returns
        -------
           position : a list of points

        """
        position = []
        p0 = self.origin
        u = np.array([np.cos(self.angle * np.pi / 180),
                      np.sin( self.angle * np.pi / 180)])
        v = np.array([-np.sin(self.angle * np.pi / 180),
                      np.cos(self.angle * np.pi / 180)])

        p1 = p0 + u * self.length
        p2 = p1 + v * self.width
        p3 = p2 - u * self.length
        position = [p0, p1, p2, p3]

        return position

    def show(self, fig=[], ax=[], offx=0, offy=10):
        """ show furnitures using shapely 

        Parameters
        ----------
            fig
                figure
            ax
                axes
            offx
                offset x (0)
            offy
                offset y  (10)

        """

        vertices = []
        codes = []
        codes = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]

        p0 = self.origin
        u = np.array([np.cos(self.angle * np.pi / 180), np.sin(
            self.angle * np.pi / 180)])
        v = np.array([-np.sin(self.angle * np.pi / 180), np.cos(
            self.angle * np.pi / 180)])
        p1 = p0 + u * self.length
        p2 = p1 + v * self.width
        p3 = p2 - u * self.length

        vertices = [p0, p1, p2, p3, p0]

        vertices = np.array(vertices, float)
        path = Path(vertices, codes)
        if self.Matname == 'METAL':
            pathpatch = PathPatch(path, facecolor='grey',
                                  edgecolor='blue', alpha=0.5)
        else:
            pathpatch = PathPatch(path, facecolor='green',
                                  edgecolor='black', alpha=0.2)

# create figure  do not exist create it
        if not isinstance(fig, plt.Figure):
            fig = plt.figure()
            ax = fig.add_subplot(111)
        #ax = fig.get_axes()[0]
        ax = fig.gca()
        ax.add_patch(pathpatch)
        #ax.dataLim.update_from_data_xy(vertices)
        ax.autoscale_view()
        return(fig,ax)


if __name__ == "__main__":
    plt.ion()
    doctest.testmod()

#   fig.savefig(figuredir+filename+ext1,orientation='portrait')
#   fig.savefig(figuredir+filename+ext2,orientation='portrait')
#   fig.savefig(figuredir+filename+ext3,orientation='portrait')
