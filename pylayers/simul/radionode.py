#/usr/bin/python
# -*- coding: utf-8 -*-
"""


RadioNode Class
===============

A radio Node is a data structure which store positions of 
a radio node at several time step. It contains also the information 
about the antenna file.

This data structure is used in DLink.  

Members of a radionode are 

position    : 3xNt
time        : 1xNt 
orientation : 3x3xNt
typ         : 'tx' or 'rx'
name        : string 
fileant     : antenna file 


.. autosummary::
    :toctree: generated/


    RadioNode.__init__
    RadioNode.__repr__
    RadioNode.pos2pt
    RadioNode.transform
    RadioNode.info
    RadioNode.clear
    RadioNode.points
    RadioNode.point
    RadioNode.linevect
    RadioNode.line
    RadioNode.surface
    RadioNode.volume
    RadioNode.loadini
    RadioNode.loadspa
    RadioNode.save
    RadioNode.gpoint
    RadioNode.gline
    RadioNode.gsurface
    RadioNode.gvolume
    RadioNode.show
    RadioNode.show3
    RadioNode.move
    RadioNode.extract
    RadioNode.loadvsh
    RadioNode.gantenna

"""

import os
import doctest
import numpy as np
import ConfigParser
import pylayers.util.easygui as eg
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geo
from pylayers.antprop.antenna import *
from pylayers.mobility.trajectory import *
from pylayers.util.project import *
import numpy as np
import scipy as sp


class RadioNode(object):
    """ container for a Radio Node

     This class manages the spatial and temporal behavior of a radio node

     Attributes
     ----------

     position
        position of the RadioNode np.array([],dtype=float)
     time
        time tag of the RadioNode np.array([],dtype=float)
     orientation
        orientation 3x3xn (rotation matrix for each position)
     points
        dictionnary of points (redundant information)
     antenneid
        id of the antenna
     type
        0: undefined 1: Tx 2 : Rx

     Methods
     -------

     info     : display information about a RadioNode
     loadspa  : load a spa file in PulsRay data format
     save     : save a RadioNode file in .spa, .ini, .vect data format
     point    : set a RadioNode point position
     points   : set a RadioNode set of points position
     line     : set a RadioNode route
     surface  : set a RadioNode area
     volume   : set a RadioNode volume
     gpoint   : set a RadioNode point position (gui)
     gline    : set a RadioNode route          (gui)
     gsurface : set a RadioNode area           (gui)
     gvolume  : set a RadioNode volume         (gui)
     show3    : display the RadioNode in the associated structure

    """

    def __init__(self, name = '',typ='undefined',
                 _fileini='radionode.ini',
                 _filetraj='pos.csv',
                 _fileant='defant.vsh3',
                 _filestr='defstr.str2'):
        """

        the _fileini file must be placed in the ini directory

        Parameters
        ----------

        typ : int
            0 : undefined
            1 : tx
            2 : rx
        _fileini : string
            file of RadioNode coordinates
        _fileant : string
            file of antenna VSH
        _filestr : string
            file of layout structure

        Notes
        -----

        The point [0,0,0] is defined as the first point (index 0)

        """

        self.position = np.array([], dtype=float)
        self.position = np.array([0, 0, 0]).reshape(3, 1)
        self.time = np.array([], dtype=float)
        self.orientation = np.eye(3).reshape(3, 3, 1)

        self.typ = typ
        self.N = 1
        self.name=name
        #
        # clean existing .ini file
        #
        if _fileini == 'radionode.ini':
            if typ == 'tx':
                _fileini = _fileini.replace('node', 'tx')
            if typ == 'rx':
                _fileini = _fileini.replace('node', 'rx')
            fileini = pyu.getlong(_fileini, 'ini')
            # delete radionode.ini if it exists
            try:
                os.remove(fileini)
            except:
                pass

        prefix = _fileini.replace('.ini','')
        prefix = prefix.replace('.spa','')
        self.fileini = prefix + '.ini'
        self.filespa = prefix + '.spa'
        self.filegeom = prefix + '.vect'
        self.filestr = _filestr
        fileini = pyu.getlong(self.fileini,'ini')
        # if file _fileini exists it is loaded
        try:
            fd = open(fileini,'r')
            fd.close()
            self.loadini(self.fileini, 'ini')
        except:
            pass

        #
        # print _fileant
        #
        self.fileant = _fileant
        try:
            self.loadvsh() # is it still necessary ? 
        except:
            raise NameError('antenna file does not exist')
        self.save()

    def __repr__(self):
        """ representation of radio node
        Only position if shown 
        """
        st = ''
        for k in range(self.N):
            st = st + str(k)+ ' : ' \
                    + str(round(self.position[0,k]*1000)/1000.) + ' ' + \
                      str(round(self.position[1,k]*1000)/1000.) + ' ' + \
                      str(round(self.position[2,k]*1000)/1000.) + '\n'
        return(st)

    def pos2pt(self):
        """ convert position to points dict

        """
        npt = np.shape(self.position)[1]
        self.points = {}
        for k in range(npt):
            self.points[k + 1] = self.position[:, k]

    def transform(self,alpha,trans):
        """ tranform position rotation + translation

        Parameters
        ----------

        alpha : float 
            angle (rad)
        trans : np.array()  (,2)

        """
        d2r = np.pi/180
        Rot = np.array([[np.cos(d2r*alpha),-np.sin(d2r*alpha)],
                        [np.sin(d2r*alpha),np.cos(d2r*alpha)]])
        self.position[0:2,:] = np.dot(Rot,self.position[0:2,:])
        self.position[0:2,:] = self.position[0:2,:]+trans[:,np.newaxis]

    def info(self):
        """ display RadioNodes informations

        """
        print "npos       : ", self.N
        print "position   : ", self.position
        print "name       : ", self.name
        #print "orientation : ", self.orientation
        print "type       : ", self.typ
        print "fileini    : ", self.fileini
        print "filespa    : ", self.filespa
        print "filegeom   : ", self.filegeom
        print "fileant    : ", self.fileant
        try:
            print "filestr    : ", self.filestr
        except:
            pass

    def clear(self):
        """ clear positions

        The origin [0,0,0] is always defined as the first point

        """
        self.position = np.array([], dtype=float)
        self.position = np.array([0., 0., 0.]).reshape(3, 1)
        self.N = 1

    def points(self, pt=np.array([[0], [0], [0]])):
        """ add a set of points to RadioNode

        Parameters
        ----------

        pt : ndarray
             point position (3 x Npt)

        """
        if type(pt) == list:
            pt = np.array(pt)

        self.position = pt
        self.N = np.shape(self.position)[1]
        self.save()

    def point(self, pt=[0, 0, 0], time=[1], orientation=[], mode='subst'):
        """ add a position to RadioNode

        The new RadioNode is saved in .spa

        Parameters
        ----------

        pt : ndarray
            point position (1 x 3)
        time : ndarray
            1x1
        orientation : ndarray
            3x3 matrix
        mode: string
            'subst' for replacement (default)
            'append' for appending

         Examples
         --------

         >>> from pylayers.simul.radionode import *
         >>> import numpy as np
         >>> tx = RadioNode()
         >>> tx.point([1,1,1],[1],np.eye(3),'subst')
         >>> tx.position
         array([[1],
                [1],
                [1]])

        """
        if isinstance(pt, list):
            pt = np.array(pt)
        if isinstance(time, list):
            time = np.array(time)
        
        orientation = np.reshape(np.eye(3), (3, 3, 1))
        pt = np.array(pt)
        time = np.array(time)
        pt = np.reshape(pt, (3, 1))

        if mode == 'subst':
            self.time = time
            self.position = pt
            self.orientation = orientation
        else:
            try:
                self.time = np.append(self.time, time, axis=0)
                self.position = np.append(self.position, pt, axis=1)
                self.orientation = np.append(self.orientation, orientation, axis=2)
            except:
                self.time = time
                self.position =  pt
                self.orientation = orientation

        self.pos2pt()
        self.save()

    def linevect(self,npt=1, step=1.0 , ptt=[0, 0, 0], vec=[1, 0, 0], mode='subst'):
        """ create a line along a direction 

        Parameters
        ----------

        npt : int 
            number of points
        step : float
            incremental distance in meters
        ptt  : list or array
            1x3 point tail (starting point)
        vec : list or arry 
            1x3 unitary vector 
        mode : string 
            'subst'
            'append'

        Examples
        --------
        
        >>> from pylayers.simul.radionode import *
        >>> r = RadioNode()
        >>> r.linevect(npt=3)
        >>> r
        0 : 0.0 0.0 0.0
        1 : 1.0 0.0 0.0
        2 : 2.0 0.0 0.0
        <BLANKLINE>

        """
        if isinstance(ptt, list):
            ptt = np.array(ptt)
        if isinstance(vec, list):
            vec = np.array(vec)
        if (npt <= 1):
            raise ValueError('npt should be greater than 1')
        ptt = np.reshape(ptt, (3, 1))
        vec = np.reshape(vec, (3, 1))
        k = np.arange(npt)
        pt = ptt + k*step*vec
        if mode == 'subst':
            self.position = pt
        else:
            self.position = np.append(self.position, pt, axis=1)
        self.pos2pt()
        self.N = np.shape(self.position)[1]
        self.save()
    
    def line(self, npt, ptt=[0, 0, 0], pth=[1, 0, 0], mode='subst'):
        """ build a line trajectory for a RadioNode

        Parameters
        ----------

        npt : integer
            number of points
        ptt : list or ndarray
            starting point coordinates  (default [0,0,0])
        ptf : list or ndarray
            ending point coordinates
        mode : string
            'subst' for replacement (default)
            'append' for appending

        Examples
        --------

        >>> from pylayers.simul.radionode import *
        >>> r = RadioNode()
        >>> r.line(3,[0,0,0],[1,0,0])
        >>> r.position
        array([[ 0. ,  0.5,  1. ],
               [ 0. ,  0. ,  0. ],
               [ 0. ,  0. ,  0. ]])


        """
        if isinstance(ptt, list):
            ptt = np.array(ptt)
        if isinstance(pth, list):
            pth = np.array(pth)
        if (npt <= 1):
            raise ValueError('npt should be greater than 1')
        ptt = np.reshape(ptt, (3, 1))
        pth = np.reshape(pth, (3, 1))
        pas = 1.0 / (npt - 1)
        k = np.arange(0.0, 1.0 + pas, pas)
        pt = ptt + k * (pth-ptt)
        if mode == 'subst':
            self.position = pt
        else:
            self.position = np.append(self.position, pt, axis=1)
        self.pos2pt()
        self.N = np.shape(self.position)[1]
        self.save()

    def surface(self, N1=2, N2=2, p0=[0, 0, 0], p1=[1, 0, 0], p2=[0, 1, 0], mode='subst'):
        """ add a surface to RadioNode

        add a surface with basis (p0p1,p0p2)

        Parameters
        ----------

        N1 : int
            default 2
        N2 : int
            default 2
        p0 : array or list
            first point
        p1 : array or list
            second point
        p2 : array or list
            third point
        mode : string
            'subst'
            'append'


        Examples
        --------

        >>> from pylayers.simul.radionode import *
        >>> tx= RadioNode()
        >>> tx.surface(10,10,[0,0,1.5],[3.0,0,1.5],[0.0,3.0,1.5],'subst')



        """
        p0 = np.array(p0)
        p1 = np.array(p1)
        p2 = np.array(p2)
        p0 = np.reshape(p0, (3, 1))
        p1 = np.reshape(p1, (3, 1))
        p2 = np.reshape(p2, (3, 1))
        pas1 = 1.0 / (N1 - 1)
        k1 = np.arange(0.0, 1.0 + pas1, pas1)
        pas2 = 1.0 / (N2 - 1)
        k2 = np.arange(0.0, 1.0 + pas2, pas2)
        n1 = len(k1)
        n2 = len(k2)
        kk1 = np.kron(np.ones(n2), k1)
        kk2 = np.kron(k2, np.ones(n1))
        pt = p0 + kk1 * (p1 - p0) + kk2 * (p2 - p0)
        if mode == 'subst':
            self.position = pt
        else:
            self.position = np.append(self.position, pt, axis=1)
        self.pos2pt()
        self.N = np.shape(self.position)[1]
        self.save()

    def volume(self,N1=2,N2=2,N3=2,p0=[0, 0, 0],p1=[1, 0, 0], p2=[0, 1, 0], p3=[0, 0, 1], mode='subst'):
        """ add a volume to RadioNode

        build a volume with edges : p0p1, p0p2, p0p3

        Parameters
        ----------

        N1 : int
            number of points on axis 1
        N2 : int
            number of points on axis 2
        N3 : int
            number of points on axis 3
        p0 : list or ndarray
            first point 
        p1 : list or ndarray
            second point 
        p2 : list or ndarray
            third point 
        p3 : list or ndarray
            fourth point 

        Examples
        --------

        >>> from pylayers.simul.radionode import *
        >>> tx = RadioNode()
        >>> tx.volume(10,10,10,[0,0,1.0],[3.0,0,1.1],[0.0,3.0,1.1],[0.0,0.0,2.0])

        """
        if isinstance(p0, list):
            p0 = np.array(p0)
        if isinstance(p1, list):
            p1 = np.array(p1)
        if isinstance(p2, list):
            p2 = np.array(p2)
        if isinstance(p3, list):
            p3 = np.array(p3)

        p0 = np.reshape(p0, (3, 1))
        p1 = np.reshape(p1, (3, 1))
        p2 = np.reshape(p2, (3, 1))
        p3 = np.reshape(p3, (3, 1))

        pas1 = 1.0 / (N1 - 1)
        k1 = np.arange(0.0, 1.0 + pas1, pas1)
        pas2 = 1.0 / (N2 - 1)
        k2 = np.arange(0.0, 1.0 + pas2, pas2)
        pas3 = 1.0 / (N3 - 1)
        k3 = np.arange(0.0, 1.0 + pas3, pas3)
        n1 = len(k1)
        n2 = len(k2)
        n3 = len(k3)
        kk1 = np.kron(np.ones(n2 * n3), k1)
        kk2 = np.kron(np.kron(np.ones(n1), k2), np.ones(n3))
        kk3 = np.kron(k3, np.ones(n1 * n2))
        pt = p0 + kk1 * (p1 - p0) + kk2 * (p2 - p0) + kk3 * (p3 - p0)
        if mode == 'subst':
            self.position = pt
        else:
            self.position = np.append(self.position, pt, axis=1)
        self.pos2pt()
        self.N = np.shape(self.position)[1]
        self.save()

    def loadini(self, _filespa, rep='ini'):
        """ load an .ini file

        Parameters
        ----------

        _filespa : string
            short filename
        rep : string
            directory name

        """
        filespa = pyu.getlong(_filespa, rep)
        #print filespa+  "   loadini"
        space = ConfigParser.ConfigParser()
        space.read(filespa)

        points = space.items("coordinates")
        self.points = pyu.lt2idic(points)
        self.N = len(self.points.keys())
        del self.position
        for k in self.points.keys():
            try:
                self.position = np.hstack((self.position,
                                           self.points[k].reshape(3,1)))
            except:
                self.position = self.points[k].reshape(3,1)

        self.traj=Trajectory(pt=self.position.T)

    def loadspa(self, _filespa, rep=pstruc['DIRLCH']):
        """ load a spa file

        Parameters
        ----------

           _filespa : string
               short filename
           rep : string 
               directory name  


        """

        self.filespa = _filespa
        filespa = pyu.getlong(_filespa, rep)

        try:
            fid = open(filespa)
        except:
            print "filespa does not exist"
            return()

        lig = fid.readlines()
        typ = int(lig[0])
        if typ == 0:
            nnpt = int(lig[1])
            coord = lig[2:]
            for index in range(len(coord)):
                point = map(float, coord[index].split())
                ndpoint = np.array([[point[0]], [point[1]], [point[2]]])
                self.position = np.append(self.position, ndpoint, axis=1)
        self.time = np.arange(nnpt)
        ident = np.eye(3)
        tmp = np.zeros(9 * nnpt)
        self.orientation = np.reshape(tmp, (3, 3, nnpt))
        self.N = nnpt
        for i in range(nnpt):
            self.orientation[:, :, i] = ident
        fid.close()

    def save(self):
        """ save RadioNode in  .ini, .spa, .vect file

        This function save the RadioNode in different files format

        .spa  : pulsray format
        .vect : geomview format

        """
        
        _filespa = self.filespa
        _fileini = self.fileini
        fileini = pyu.getlong(_fileini, 'ini')
        try:
            fd = open(fileini, "w")
        except:
            print fileini + ' does not exist'

        space = ConfigParser.ConfigParser()
        space.add_section("coordinates")
        npt = np.shape(self.position)[1]

        for k in range(npt):
            x = self.position[0, k]
            y = self.position[1, k]
            z = self.position[2, k]
            space.set("coordinates", str(k + 1), str(x) + ' ' +
                      str(y) + ' ' + str(z))
        space.write(fd)
        fd.close()

        points = space.items("coordinates")

        if self.typ == 'undefined':
            filespa = pyu.getlong(_filespa, 'ini')
            colorname = 'green'
        elif self.typ == 'tx':
            filespa = pyu.getlong(_filespa, pstruc['DIRLCH'])
            colorname = 'red'
        elif self.typ == 'rx':
            filespa = pyu.getlong(_filespa, pstruc['DIRTRA'])
            colorname = 'blue'

        # save points in GeomVect container
        filename = self.filegeom.replace('.vect', '')
        filename = filename.replace('.off', '')
        try: 

            gv = geo.Geomoff(filename)
            ant = self.A
            # if not hasattr(ant,'theta'):
            if ant.fromfile:
                ant.Fsynth3(pattern=True)
            else: 
                ant.Fpatt(pattern=True)
            V = ant.SqG[ant.Nf/2,:,:]
 
            if not hasattr(self,'position'):
                print "no position available"
            T=self.orientation.reshape(3,3)
            gv.pattern(ant.theta,ant.phi,V,po=self.position,T=T,ilog=False,minr=0.01,maxr=1.)
            self.filegeom=filename + '.off'
        except:
            #pdb.set_trace()
            if hasattr(self,'position'):
                pass
            #    problem gv is not a geomvect    
            #    gv.points(self.position, colorname)
            else :
                print " no position available "

        if _filespa.split('.')[1] == 'spa':
            fi_spa = open(filespa, 'w')
            npt = np.shape(self.position)[1]
            snpt = str(npt) + "\n"
            snpt2 = str(npt) + " " + str(npt) + " " + str(npt) + "\n"
            fi_spa.write("0\n")
            fi_spa.write(snpt)
            for i in range(npt):
                x = str(self.position[0, i]).replace(',', '.')
                y = str(self.position[1, i]).replace(',', '.')
                z = str(self.position[2, i]).replace(',', '.')
                chaine = x + " " + y + " " + z + "\n"
                chaine2 = chaine.replace(',', '.')
                fi_spa.write(chaine)

            fi_spa.close()

    def gpoint(self, mode='subst', display=False):
        """ gui point 

        Parameters
        ----------

        mode : string 
            subst
        display : boolean 
            False

        """
        p0 = self.position[:, 0]
        (p0, n1) = eg.pointbox(p0, 1)
        self.point(p0, [1], np.eye(3), mode)
        self.save()
        if display:
            self.show3()

    def gline(self, mode='subst', display=False):
        """ gui line

        A line is built between the first point and the gui point

        Parameters
        ----------

        mode : string 
            subst
        display : boolean 
            False

        """
        p0 = self.position[:, 0]
        (p1, N1) = eg.pointbox(p0, 10)
        self.line(N1, p0, p1, mode)
        self.save()
        if display:
            self.show3()

    def gsurface(self, mode='subst', display=False):
        """ gui surface

        Parameters
        ----------

        mode : string 
            subst
        display : boolean 
            False

        """
        p0 = self.position[:, 0]
        (p1, N1) = eg.pointbox(p0, 10, 'Enter Surface second point')

        (p2, N2) = eg.pointbox(p1, 10, 'Enter Surface third point')
        self.surface(N1, N2, p0, p1, p2, mode)
        self.save()
        if display:
            self.show3()

    def gvolume(self, mode='subst', display=False):
        """ gui volume 

        Parameters
        ----------

        mode : string 
            subst
        display : boolean 
            False

        """
        p0 = self.position[:, 0]
        (p1, N1) = eg.pointbox(p0, 10, 'Enter Volume second point')
        (p2, N2) = eg.pointbox(p1, 10, 'Enter Volume third point')
        (p3, N3) = eg.pointbox(p2, 10, 'Enter Volume fourth point')
        self.volume(N1, N2, N3, p0, p1, p2, p3, mode)
        self.save()
        if display:
            self.show3()
#        def savevect(self):
#                """
#                Create a .vect file
#                Le type de format est 0 . Coordonnees explicites de tous les points.
#
#       save(_filespa)
#
#       _filespa : file short name
#
#               """
#
#                if self.typ==0:
#                        self.filegeom="RadioNode.vect"
#                        filegeom   = pyu.getlong("RadioNode.vect","geom")
#                elif self.typ==1:
#                        self.filegeom = "RadioTx.vect"
#                        filegeom   = pyu.getlong("RadioTx.vect","geom")
#                elif self.typ==2:
#                        self.filegeom = "RadioRx.vect"
#                        filegeom   = pyu.getlong("RadioRx.vect","geom")
#
#                fi_geom = open(filegeom,'w')
#
#                npt = shape(self.position)[1]
#                snpt2 = str(npt)+" "+str(npt)+" "+str(npt)+"\n"
#                if npt>1:
#                        fi_geom.write("appearance{\n")
#                        fi_geom.write("linewidth 8}\n")
#                        fi_geom.write("VECT\n")
#                        fi_geom.write(snpt2)
#                        fi_geom.write("1 "*npt+"\n")
#                        fi_geom.write("1 "*npt+"\n")
#                else:
#                        fi_geom.write("ESPHERE\n")
#                        fi_geom.write("0.2\n")
#
#                for i in range(npt):
#                        x = str(self.position[0,i])
#                        y = str(self.position[1,i])
#                        z = str(self.position[2,i])
#                        chaine = x+" "+y+" "+z+"\n"
#                        fi_geom.write(chaine)
#                if npt>1:
#                        if self.typ==0:
#                                fi_geom.write(npt*"0 0 1 1\n")
#                        elif self.typ==1:
#                                fi_geom.write(npt*"1 0 0 1\n")
#                        elif self.typ==2:
#                                fi_geom.write(npt*"0 1 0 1\n")
#                fi_geom.close()

#        def savespa(self):
#                """
#                Create a  .spa file
#                Le type de format est 0 . Coordonnees explicites de tous les points.
#
#       savespa(_filespa)
#
#       _filespa : file short name
#
#                """
#
#               _filespa = self.filespa
#
#                elif self.typ==1:
#                        filespa    = pyu.getlong(_filespa,'launch')
#                elif self.typ==2:
#                        filespa    = pyu.getlong(_filespa,'trace')
#
#                fi_spa  = open(filespa,'w')
#
#                npt = shape(self.position)[1]
#                snpt = str(npt)+"\n"
#
#                fi_spa.write("0\n")
#                fi_spa.write(snpt)
#                for i in range(npt):
#                        x = str(self.position[0,i])
#                        y = str(self.position[1,i])
#                        z = str(self.position[2,i])
#                        chaine = x+" "+y+" "+z+"\n"
#                        fi_spa.write(chaine)
#                fi_spa.close()

    def show(self, num = [], fig=[], ax =[],size=5,marker='o',color='b'):
        """ Display RadioNode position in the 2D strucure

        Parameters
        ----------


        """
        if num ==[]:
            num = np.arange(np.shape(self.position)[1])
        
        x = self.position[0, num]
        y = self.position[1, num]

        if ax == []:
            fig = plt.gcf()
            ax = fig.gca()
        ax.scatter(x, y,s=size,c=color,linewidth=0)
        return fig,ax

    def show3(self, _filestr='DLR.off'):
        """ display RadioNode position in geomview

        Parameters
        ----------

        _filestr : string
            structure file
        """
        filename = pyu.getlong("strucRN.off", pstruc['DIRGEOM'])
        fo = open(filename, "w")
        filegeom = pyu.getlong(self.filegeom, pstruc['DIRGEOM'])

        # get .off filename from .str or .str2 filename
        _fileoff, ext = os.path.splitext(self.filestr)
        _fileoff = _fileoff + '.off'
        fo.write("LIST\n")
        if os.path.isfile(pyu.getlong(_fileoff,'geom')):
            fo.write("{<" + _fileoff + "}\n")
        fo.write("{<" + self.filegeom + "}\n")
        fo.write("{</usr/share/geomview/geom/xyz.vect}\n")
        fo.close()
        chaine = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
        os.system(chaine)

    def move(self, dx, dy, dz):
        """ move RadioNode with a specified offset over each cartesian axis

        Parameters
        ----------
        dx : float
        dy : float
        dz : float

        """
        self.position[0, :] += dx
        self.position[1, :] += dy
        self.position[2, :] += dz
        self.save()

    def extract(self, i):
        """ extract the i-th radionode component 

        Parameters
        ----------

        i : integer

        Returns
        -------
        
        u : RadioNode


        """

        if self.typ == 'undefined':
            u = RadioNode(self.filestr)
        elif self.typ == 'tx':
            u = RadioTx(self.filestr, self.signal)
        elif self.typ == 'rx':
            u = RadioRx(self.filestr, self.fc, self.bandwidth, self.NF)

        u.position = self.position[:, i]
#
#
#       u.time        = self.time[i]
#       u.orientation = self.orientation[:,:,i]
        u.filespa = "filespa.spa"

#
# Write the RadioNode Coordinate in filespa
#
        if self.typ != 'undefined':
            if self.typ == 'tx':
                filespa = pyu.getlong("filespa.spa", pstruc['DIRLCH'])
            elif self.typ == 'rx':
                filespa = pyu.getlong("filespa.spa", pstruc['DIRTRA'])
            fi = open(filespa, 'w')
            fi.write("0\n")
            fi.write("1\n")
            x = str(self.position[0, i])
            y = str(self.position[1, i])
            z = str(self.position[2, i])
            chaine = x + " " + y + " " + z + "\n"
            fi.write(chaine)
            fi.close()

        return u

    def loadvsh(self):
        """ load an antenna 


        """
        #print self.fileant
        A = Antenna(self.fileant)
        self.A = A

    def gantenna(self, mode='subst'):
        """ get antenna file
        """
        import tkFileDialog
        FD = tkFileDialog

        fileant = FD.askopenfilename(filetypes=[("Fichiers vsh3", "*.vsh3"),
                                                ("All", "*")],
                                     title="Please choose an antenna file",
                                     initialdir=antdir)

        _fileant = os.path.split(fileant)[1]
        self.fileant = _fileant
        self.loadvsh()

if (__name__ == "__main__"):
    tx = RadioNode(_fileini='w2m1rx.ini')
    doctest.testmod()
