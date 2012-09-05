#/usr/bin/python
# -*- coding: utf-8 -*-
import os
import doctest
import numpy as np
import ConfigParser
import pylayers.util.geomutil as geo
import pylayers.util.easygui
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu 
from pylayers.antprop.antenna import *
from numpy import *
from scipy import *


class RadioNode(object):
    """ container for a Radio Node

     a RadioNode is either a transmitter or a receiver
     This class manages the spatial and temporal behavior of a radio node

     Attributes
     ----------

     position
        position of the RadioNodenp.array([],dtype=float)
     time
        time tag of the RadioNodenp.array([],dtype=float)
     orientation
        orientation 3x3xn (rotation matrix for each position)
     antenneid
        id of the antenna
     type
        0: undefined 1: Tx 2 : Rx

     Methods
     -------
     info     : display information about a RadioNode
     loadspa  : load a spa file in PulsRay data format
     savespa  : save a spa file in PulsRay data format
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

    def __init__(self, _filespa, _fileant, type=0):
        """
        Parameters
        ----------
        _filespa : string
            short filename of the

        type : integer
             1=Transmitter 2=Receiver

        This function can read either a spa or an ini file 
        For compatibility with pulsray we have to keep the old
        .spa format 

        """

        self.position = np.array([], dtype=float)
        self.position.shape = (3, 0)
        self.time = np.array([], dtype=float)
        self.orientation = np.array([], dtype=float)
        self.orientation.shape = (3, 3, 0)
        self.antenneid = 0
        self.type = type
        if type == 1:
            self.filespa = _filespa
            self.fileant = _fileant
            if self.filespa.split('.')[1] == 'spa':
                self.loadspa(self.filespa, 'launch')
            else:
                self.loadini(self.filespa, 'launch')
            self.loadvsh()
            self.savespa()
        if type == 2:
            self.filespa = _filespa
            self.fileant = _fileant
            if self.filespa.split('.')[1] == 'spa':
                self.loadspa(self.filespa, 'trace')
            else:
                self.loadini(self.filespa, 'trace')
                self.N = len(self.points.keys())
            self.loadvsh()
            self.savespa()
            
        self.N = len(self.points.keys())
        for k in self.points.keys():
            try:
                self.position = np.vstack((self.position, self.points[k]))
            except:
                self.position = self.points[k]

    def info(self):
        """ display RadioNodes informations
        """
        print "Npos       : ", self.N
        print "Position   : ", self.position
        print "type       : ", self.type
        print "filespa    : ", self.filespa
        print "filegeom   : ", self.filegeom
        print "fileant    : ", self.fileant

    def points(self, pt=array([0, 0, 0])):
        """ Add a position to RadioNode

        Parameters
        ----------
        pt : ndarray
             point position (3 x Npt)
        """
        pt = np.array(pt)
        self.position = pt
        self.savespa()

    def point(self, pt=array([0, 0, 0]), time=array([1]), orientation=eye(3), mode='subst'):
        """
        Add a position to RadioNode

        Parameters
        ----------
            pt
                point position (1 x 3)
            time
                1x1
            orientation
                3x3 matrix
            mode
                'subst' for deplacement (default)

         >>> from pylayers.simul.radionode import *
         >>> tx = RadioNode()
         >>> tx.point([1,1,1],[1],eye(3),'subst')
         >>> tx.info()

        """
        orientation = np.reshape(orientation, (3, 3, 1))
        pt = np.array(pt)
        time = np.array(time)
        pt = np.reshape(pt, (3, 1))

        if mode == 'subst':
            self.time = time
            self.position = pt
            self.orientation = orientation
        else:
            self.time = append(self.time, time, axis=0)
            self.position = append(self.position, pt, axis=1)
            self.orientation = append(self.orientation, orientation, axis=2)

        self.savespa()

    def line(self, Npt, pti=[0, 0, 0], ptf=[1, 0, 0], mode='subst'):
        """
        Add a line to RadioNode

        Usage :  tx.route(Npt,[0,0,0],[1,1,1])

        """
        pti = np.array(pti)
        pff = np.array(ptf)
        pti = np.reshape(pti, (3, 1))
        ptf = np.reshape(ptf, (3, 1))
        pas = 1.0 / (Npt - 1)
        k = arange(0.0, 1.0 + pas, pas)
        pt = ptf + (1.0 - k) * (pti - ptf)
        if mode == 'subst':
            self.position = pt
        else:
            self.position = append(self.position, pt, axis=1)

    def surface(self, N1=2, N2=2, p0=[0, 0, 0], p1=[1, 0, 0], p2=[0, 1, 0], mode='subst'):
        """
        Add a surface to RadioNode

        >>> tx.surface(10,10,[0,0,1.5],[3.0,0,1.5],[0.0,3.0,1.5],'subst')

        mode = {'subst','append' }

        """
        p0 = np.array(p0)
        p1 = np.array(p1)
        p2 = np.array(p2)
        p0 = np.reshape(p0, (3, 1))
        p1 = np.reshape(p1, (3, 1))
        p2 = np.reshape(p2, (3, 1))
        pas1 = 1.0 / (N1 - 1)
        k1 = arange(0.0, 1.0 + pas1, pas1)
        pas2 = 1.0 / (N2 - 1)
        k2 = arange(0.0, 1.0 + pas2, pas2)
        n1 = len(k1)
        n2 = len(k2)
        kk1 = kron(ones(n2), k1)
        kk2 = kron(k2, ones(n1))
        pt = p0 + kk1 * (p1 - p0) + kk2 * (p2 - p0)
        if mode == 'subst':
            self.position = pt
        else:
            self.position = append(self.position, pt, axis=1)

    def volume(self, N1=2, N2=2, N3=2, p0=[0, 0, 0], p1=[1, 0, 0], p2=[0, 1, 0], p3=[0, 0, 1], mode='subst'):
        """
        Add a volume to RadioNode

        Usage :  tx.volume(N1,N2,N3,p0,p1,p2,p3,mode)

        >>> tx.volume(10,10,10,[0,0,1.0],[3.0,0,1.1],[0.0,3.0,1.1],[0.0,0.0,2.0])

        """
        p0 = np.array(p0)
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
        p0 = np.reshape(p0, (3, 1))
        p1 = np.reshape(p1, (3, 1))
        p2 = np.reshape(p2, (3, 1))
        p3 = np.reshape(p3, (3, 1))
        pas1 = 1.0 / (N1 - 1)
        k1 = arange(0.0, 1.0 + pas1, pas1)
        pas2 = 1.0 / (N2 - 1)
        k2 = arange(0.0, 1.0 + pas2, pas2)
        pas3 = 1.0 / (N3 - 1)
        k3 = arange(0.0, 1.0 + pas3, pas3)
        n1 = len(k1)
        n2 = len(k2)
        n3 = len(k3)
        kk1 = kron(ones(n2 * n3), k1)
        kk2 = kron(kron(ones(n1), k2), ones(n3))
        kk3 = kron(k3, ones(n1 * n2))
        pt = p0 + kk1 * (p1 - p0) + kk2 * (p2 - p0) + kk3 * (p3 - p0)
        if mode == 'subst':
            self.position = pt
        else:
            self.position = append(self.position, pt, axis=1)

    def loadini(self, _filespa, rep='launch'):
        """ load an ini file

        Parameters
        ----------
        _filespa : string
            short filename
        rep : string
            directory name

        """
        filespa = getlong(_filespa, rep)
        self.space = ConfigParser.ConfigParser()
        self.space.read(filespa)

        points = self.space.items("coordinates")
        self.points = lt2idic(points)

    def loadspa(self, _filespa, rep='launch'):
        """
        load a spa file

        Parameters
        ----------
           _filespa : short filename

          .. todo::
              gérer le problème des lignes blanches !!
        """

        self.filespa = _filespa
        filespa = getlong(_filespa, rep)

        try:
            fid = open(filespa)
        except:
            print "filespa does not exist"
            return()

        lig = fid.readlines()
        type = int(lig[0])
        if type == 0:
            nnpt = int(lig[1])
            coord = lig[2:]
            for index in range(len(coord)):
                point = map(float, coord[index].split())
                ndpoint = np.array([[point[0]], [point[1]], [point[2]]])
                self.position = append(self.position, ndpoint, axis=1)
        self.time = arange(nnpt)
        ident = eye(3)
        tmp = zeros(9 * nnpt)
        self.orientation = np.reshape(tmp, (3, 3, nnpt))
        self.N = nnpt
        for i in range(nnpt):
            self.orientation[:, :, i] = ident
        fid.close()

    def savespa(self, k=-1):
        """ save RadioNode in a .spa file and a .vect file

        Parameters
        ----------
        k : integer 
            
        .. todo: save as an ini file 

        """
        _filespa = self.filespa
        # appent rx number to filename
        if k != -1:
            _filespa = _filespa.replace('.spa', str(k))
            _filespa = _filespa + '.spa'
        if self.type == 0:
            self.filegeom = "RadioNode"
            colorname = 'green'
        elif self.type == 1:
            self.filegeom = "RadioTx"
            filespa = getlong(_filespa, 'launch')
            colorname = 'red'
        elif self.type == 2:
            self.filegeom = "RadioRx"
            filespa = getlong(_filespa, 'trace')
            colorname = 'blue'

        gv = geo.GeomVect(self.filegeom)
        try:
            gv.points(self.points, colorname)
        except:
            gv.points(self.position, colorname)
        ##
        ## ..todo::
        ##    This is a fix : temporary
        ##
        if _filespa.split('.')[1] == 'spa':
            fi_spa = open(filespa, 'w')
            npt = shape(self.position)[1]
            snpt = str(npt) + "\n"
            snpt2 = str(npt) + " " + str(npt) + " " + str(npt) + "\n"
            fi_spa.write("0\n")
            if k == -1:
                fi_spa.write(snpt)
                for i in range(npt):
                    x = str(self.position[0, i]).replace(',', '.')
                    y = str(self.position[1, i]).replace(',', '.')
                    z = str(self.position[2, i]).replace(',', '.')
                    chaine = x + " " + y + " " + z + "\n"
                    chaine2 = chaine.replace(',', '.')
                    fi_spa.write(chaine)
            else:
                fi_spa.write("1\n")
                x = str(self.position[0, k]).replace(',', '.')
                y = str(self.position[1, k]).replace(',', '.')
                z = str(self.position[2, k]).replace(',', '.')
                chaine = x + " " + y + " " + z + "\n"
                chaine2 = chaine.replace(',', '.')
                fi_spa.write(chaine)

            fi_spa.close()
#        if npt>1:
#            if self.type==0:
#                fi_geom.write(npt*"0 0 1 1\n")
#            elif self.type==1:
#                fi_geom.write(npt*"1 0 0 1\n")
#            elif self.type==2:
#                fi_geom.write(npt*"0 1 0 1\n")
#        fi_geom.close()

    def gpoint(self, mode='subst', dispaly=False):
        """
        """
        p0 = self.position[:, 0]
        (p0, n1) = easygui.pointbox(p0, 1)
        self.point(p0, [1], eye(3), mode)
        self.savespa()
        if display:
            self.show3()

    def gline(self, mode='subst', display=False):
        """
        gline(mode='subst',display=False)
        """
        p0 = self.position[:, 0]
        (p1, N1) = easygui.pointbox(p0, 10)
        self.line(N1, p0, p1, mode)
        self.savespa()
        if display:
            self.show3()

    def gsurface(self, mode='subst', display=False):
        """
        gline(mode='subst',display=False)
        """
        p0 = self.position[:, 0]
        (p1, N1) = easygui.pointbox(p0, 10, 'Enter Surface second point')

        (p2, N2) = easygui.pointbox(p1, 10, 'Enter Surface third point')
        self.surface(N1, N2, p0, p1, p2, mode)
        self.savespa()
        if display:
            self.show3()

    def gvolume(self, mode='subst', display=False):
        """
        gline(mode='subst',view=0)
        """
        p0 = self.position[:, 0]
        (p1, N1) = easygui.pointbox(p0, 10, 'Enter Volume second point')
        (p2, N2) = easygui.pointbox(p1, 10, 'Enter Volume third point')
        (p3, N3) = easygui.pointbox(p2, 10, 'Enter Volume fourth point')
        self.volume(N1, N2, N3, p0, p1, p2, p3, mode)
        self.savespa()
        if display:
            self.show3()
#        def savevect(self):
#                """
#                Create a .vect file
#                Le type de format est 0 . Coordonnées explicites de tous les points.
#
#       savespa(_filespa)
#
#       _filespa : file short name
#
#               """
#
#                if self.type==0:
#                        self.filegeom="RadioNode.vect"
#                        filegeom   = getlong("RadioNode.vect","geom")
#                elif self.type==1:
#                        self.filegeom = "RadioTx.vect"
#                        filegeom   = getlong("RadioTx.vect","geom")
#                elif self.type==2:
#                        self.filegeom = "RadioRx.vect"
#                        filegeom   = getlong("RadioRx.vect","geom")
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
#                        if self.type==0:
#                                fi_geom.write(npt*"0 0 1 1\n")
#                        elif self.type==1:
#                                fi_geom.write(npt*"1 0 0 1\n")
#                        elif self.type==2:
#                                fi_geom.write(npt*"0 1 0 1\n")
#                fi_geom.close()

#        def savespa(self):
#                """
#                Create a  .spa file
#                Le type de format est 0 . Coordonnées explicites de tous les points.
#
#       savespa(_filespa)
#
#       _filespa : file short name
#
#                """
#
#               _filespa = self.filespa
#
#                elif self.type==1:
#                        filespa    = getlong(_filespa,'launch')
#                elif self.type==2:
#                        filespa    = getlong(_filespa,'trace')
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
    def show(self, num, sp):
        """
         Display RadioNode position in the 2D strucure
        """
        x = self.position[0, num]
        y = self.position[1, num]
        sp.plot(x, y, 'ob')
        #indoor = IndoorStr()
        #filename = getlong(self.simul.filestr,'struc')
        #indoor.load(filename)
        #pt =self.position
        #indoor.show([0],[0],pt)

    def show3(self):
        """
         Display RadioNode position in the 3D strucure
         ..todo::
             Dans GeomUtil classe liste construire un fichier .list à partir d'une liste
             de fichiers
        """
        filename = getlong("strucRN.off", "geom")
        filegeom = getlong(self.filegeom, "geom")
        fo = open(filename, "w")
        fo.write("LIST\n")
    #       fo.write("{<RadioNode.vect}\n")
        fo1 = open(filegeom)
        fo.write(fo1.read())
        fo.write("\n")
        fo1.close()

        fo.write("{<struc.off}\n")
        fo.write("{</usr/share/geomview/geom/xyz.vect}\n")
        fo.close()
        chaine = "geomview -nopanel -b 1 1 1 " + filename + " 2>/dev/null &"
        os.system(chaine)

    def move(self, dx, dy, dz):
        """
          Move RadioNode with a specified offset over each cartesian axis
        """
        self.position[0, :] += dx
        self.position[1, :] += dy
        self.position[2, :] += dz
        self.savespa()

    def extract(self, i):
        """
        Extract the i-th radionode component (i=0 first position)
        .. todo::
            Vérifier si deepcopy peut faire la même chose en + court
        """

        if self.type == 0:
            u = RadioNode(self.filestr)
        elif self.type == 1:
            u = RadioTx(self.filestr, self.signal)
        elif self.type == 2:
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
        if self.type != 0:
            if self.type == 1:
                filespa = getlong("filespa.spa", "launch")
            elif self.type == 2:
                filespa = getlong("filespa.spa", "trace")
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
        """ load an antenna .vsh3 file

        """

        A = Antenna('vsh3',self.fileant)
        self.A = A

    def gantenna(self, mode='subst'):
        """
        gantenna(self,mode='subst'):
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

        if self.type == 1:
            self.fileantTx = _fileant
        if self.type == 2:
            self.fileantRx = _fileant


if (__name__ == "__main__"):
    import os
    from pylayers.simul.simulem import *
    doctest.testmod()

    _filestr = 'ceamimo2.str'
    _fileslab = "def.slab"
    _filemat = "def.mat"
    _filefreq = "def.freq"
    _filepalch = "def.palch"
    _filepatra = "def.patra"
    _filespaTx = "Tx.spa"
    _filespaRx = "Rx.spa"
