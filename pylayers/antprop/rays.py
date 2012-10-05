#!/usr/bin/python
# -*- coding: latin1 -*-
import pdb
import os
import pdb
import glob
import doctest
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import struct as stru
#from math import *
#from Indoor import *
import pylayers.util.geomutil as geu
import pylayers.util.pyutil as pyu
from   pylayers.util.project import *
import pylayers.antprop.slab as slab
#
#  This file contains
#
#       class Interaction
#           info
#           class IntB(Interaction)
#              eval(f)
#           class IntL(Interaction)
#              eval(f)
#           class IntR(Interaction)
#              eval(f)    TBF
#           class IntT(Interaction)
#           class IntD(Interaction)
#
#       class Ray3D
#           info()
#           locbas(L)
#           delay()
#           show()
#           show3(bdis=True,bbas=False,bstruc=True,col=np.array([1,0,1]),id=0)
#
#       class RayTud
#           info()
#           delay()  compliqué !!
#           eval(f)
#
#       class GrRayTud
#           delay()
#           choose()
#           load(_filename,sl)
#           info(n=-1)
#
#       class GrRay3D
#           info(level=0)
#           choose()
#           load(_filename,struc)
#           delay()
#           save(_filename)
#           reverse(struc)
#           show()
#           show3()
#
# A Faire :
#
#    Finir de coder les évaluations d'interaction
#    Construire le pont entre Ray3D et RayTud  ( Validation )
#        equivalent de tra2tud
#               CreateRays_TudFile2 dans Tud/Tud_rayon.c
#               GetRayonsParam()
#
#    1 ) Lire simultanément .tra .tud .field et visualiser le tout
#    2 ) Passer de .tra --> .tud
#    3 ) Passer de .tud --> .field
#
#class Interaction(object):
#    """
#       Interaction coefficient
#       """
#         def __init__(self,intpara,freq):
#        self.intpara=intpara
#        Nf = len(freq)
#        self.Ci = np.array(np.zeros([Nf,2,2]))
#    def eval(self,freq):
#        """
#        Evaluate Interaction Matrix
#        """
#        typ=self.inpara[0]
#        if (typ==-1):
#            self.Ci=intpara[1]
#        elif (typ==0):
#            dist = self.intpara[1]
#            self.Ci=eye(2)/dist
#        elif (typ==1):
#            theta = self.intpara[1]
#            Si    = self.intpara[2]
#            Sr    = self.intpara[3]
#            if (rho1 == 0.0):
#                rho1  = 1/Si
#                rho2  = 1/Si
#
#        #elif (typ==2):
#        #elif (typ==3):


class Interaction(object):
    """ Interaction parameters

    Notes
    -----
    This class contains all the informations about a given Interaction

    type :  Interaction type
            -1  -  Ray's local basis
                   m11 m12
                   m21 m22
             0  -  Direct (LOS)
                   LOS distance
             1  -  Reflexion
                   theta
                   Si
                   Sr
             2  -  Transmission
                   theta
                   Si
                   St
             3  -  Difraction
                   theta
                   theta_d
                   Si
                   Sd
                   beta0
                   N
                   typed

    Mat  : ( Mav1 , Mslab1 , Map1 , sense1 , Mav2 , Mslab2 , Map2 , sense2 )

    Methods
    -------

    info()

    """
    def __init__(self, typ=0):
        self.typ = typ

    def info(self):
        if (self.typ == -1):
            print "local basis"
            print "---"
        if (self.typ == 0):
            print "LOS"
            print "---"
        if (self.typ == 1):
            print "Reflexion"
            print "---"
        if (self.typ == 2):
            print "Transmission"
            print "---"
        if (self.typ == 3):
            print "Diffraction"
            print "---"

        if (self.typ == -1):
            print "M : ", self.M
        if (self.typ == 0):
            print "dist : ", self.dist
        if (self.typ == 1):
            print "theta : ", self.theta
            print "si : ", self.si
            print "sr : ", self.sr
            for i in range(len(self.Mat1)):
                print self.Mat1[i]['name']
        if (self.typ == 2):
            print "theta : ", self.theta
            print "si : ", self.si
            print "st : ", self.st
        if (self.typ == 3):
            print "theta ", self.theta
            print "thetad ", self.thetad
            print "self.si ", self.si
            print "self.sd ", self.sd
            print "self.beta ", self.beta
            print "self.N ", self.N
            print "self.typed ", self.typed


#        print "Mat1 : left  <=> Mav1    : ",self.Mat[0]
#        print "Mat1 : core  <=> Mslab1  : ",self.Mat[1]
#        print "Mat1 : right <=> Map1    : ",self.Mat[2]
#        print "Mat1 : sense <==> sense1 : ",self.Mat[3]
#        print "Mat2 : left  <=> Mav2    : ",self.Mat[4]
#        print "Mat2 : core  <=> Mslab2  : ",self.Mat[5]
#        print "Mat2 : right <=> Map2    : ",self.Mat[6]
#        print "Mat2 : sense <==> sense2 : ",self.Mat[7]

class IntB(Interaction):
    """ Local Basis interaction class

    Notes
    ------

    """
    def __init__(self, typ, M):
        Interaction.__init__(self, typ)
        self.M = M

    def eval(self, fGHz):
        """ evaluate interaction

        Parameters
        ----------
        fGHz :  float

        Notes
        -----
        repeat the M matrix along the frequency axis (axis = 0 )

        .. todo::  consider using broadcasting
        """

        nf = len(fGHz)
        Co = np.array([self.M])
        return(Co.repeat(nf, axis=0))


class IntL(Interaction):
    """ LOS interaction class


    """
    def __init__(self, typ, dist):
        Interaction.__init__(self, typ)
        self.dist = dist

    def eval(self, f):
        nf = len(f)
        Co = np.array(np.zeros([nf, 2, 2]), dtype=complex)
        div = 1.0 / self.dist
        Co[:, 0, 0] = div
        Co[:, 1, 1] = div
        return(Co)


class IntR(Interaction):
    """ Reflexion interaction class
    """
    def __init__(self, typ, data):
        Interaction.__init__(self, typ)
        self.theta = data[0]
        self.si = data[1]
        self.sr = data[2]

    def eval(self, fGHz):
        """
        .. todo:: Reflexion interaction is not implemented yet

        Parameters
        ----------
        fGHz : float
            frequency in GHz

        """
        div = np.sqrt((ro1 * ro2) / ((dr + ro1) * (dr + ro2)))
        si = self.si
        sr = self.sr
        theta = self.theta


class IntT(Interaction):
    """ Transmission  interaction class
    """
    def __init__(self, typ, data):
        Interaction.__init__(self, typ)
        self.theta = data[0]
        self.si = data[1]
        self.st = data[2]


class IntD(Interaction):
    """ Diffraction interaction class
    """
    def __init__(self, typ, data1, data2):
        Interaction.__init__(self, typ)
        self.theta = data1[0]
        self.thetad = data1[1]
        self.si = data1[2]
        self.sd = data1[3]
        self.beta = data1[4]
        self.N = data1[5]
        self.typed = data2[0]

#    def eval(self,freq):
#        """
#        Evaluate Interaction Matrix
#        """
#        typ=self.inpara[0]
#        if (typ==-1):
#            self.Ci=intpara[1]
#
#        elif (typ==0):
#            dist = self.intpara[1]
#            self.Ci=eye(2)/dist
#
#        elif (typ==1):
#            theta = self.intpara[1]
#            Si    = self.intpara[2]
#            Sr    = self.intpara[3]
#
#        elif (typ==2):
#            theta = self.intpara[1]
#            Si    = self.intpara[2]
#            St    = self.intpara[3]
#
#        elif (typ==3):
#            theta = self.intpara[1]
#            thetad = self.intpara[2]
#            Si    = self.intpara[3]
#            Sd    = self.intpara[4]
#            beta0 = self.intpara[5]
#            n = self.intpara[6]
#            typd = self.intpara[7]


class Ray2D(object):
    """ 2D Ray class

    Attributes
    ----------
    nn      :
        number of nodes (number of interactions + 2)
    pt      :
        point coordinates
    signature   :
        (new name for nstr)
    """
    def __init__(self):
        pass

    def eval(self, L, signature):
        pass


class Ray3D(object):
    """ Ray3D class

    This class handles a 3D ray

    Attributes
    ----------
    nn      : int
        number of nodes (number of interaction + 2)
    pt      : np.array
        point coordinates
    nstr    : element of a signature
        structure number (>0 edge) (<0 node) (0 Tx or Rx)
    deep    : int
    beta0   : float
    phii    : float
    phid    : float
    length  : float
    Bi   : interaction input basis    shape=( nn,3,3)
    Bo   : intercation output basis   shape=( nn,3,3)

    Bi(l,:,0) = s_in
    Bi(l,:,1) = etheta_in
    Bi(l,:,2) = ephi_in

    Bo(l,:,0) = s_out
    Bo(l,:,1) = etheta_out
    Bo(l,:,2) = ephi_out

    Methods
    -------
    show3(bdis=True,col=np.array([1,0,1]),id=0) : geomview visualization
    info()  : Ray information

    """
    def info(self):
        """ info
        """
        print "Number of nodes : " + str(self.nn)
        print "points : "
        print self.pt.transpose()
        print "nstr   : "
        print self.nstr
        print "B in     : "
        print self.Bi
        print "B out    : "
        print self.Bo
        print "AOD : "
        print "thetat (deg)   : " + str(self.tht * 180 / np.pi)
        print "phit   (deg)   : " + str(self.pht * 180 / np.pi)
        print "AOA : "
        print "thetar (deg)   : " + str(self.thr * 180 / np.pi)
        print "phir   (deg)   : " + str(self.phr * 180 / np.pi)

    def locbas(self, L):
        """ evaluation of local basis over the ray

        Parameters
        -----------
        L : Layout Structure Graph

        Notes
        -----
            evaluate Bi and Bo
        """

        # initialize RayTud
        self.gt = RayTud()
        self.gt.inter = []

        nn = self.nn
        #
        # si : (N-1,3)
        #
        #si  = np.array(np.zeros([nn-1,3],dtype=np.float64))
        #lsi = np.array(np.zeros(nn-1,dtype=np.float64))
        #
        # Calcul des vecteurs unitaires si
        # .. todo:: remove this for loop
        #   not tested

        v = self.pt[1:, :] - self.pt[0:-1, :]
        lsi = np.sqrt(np.sum(v * v, axis=1))
        # reshape is required for broadcasting
        lsir = lsi.reshape(nn - 1, 1)
        si = v / lsir
        #for k in range(nn-1):
        #    v  = self.pt[k+1,:]-self.pt[k,:]
        #    nv = np.sqrt(np.dot(v,v))
        #    try:
        #        si[k,:]  = v/nv
        #        lsi[k,:] = nv
        #    except:
        #        print("error divide by 0 in Ray3D.geom ")

        self.Bi = np.array(np.zeros([nn - 1, 3, 3], dtype=np.float64))
        self.Bo = np.array(np.zeros([nn - 1, 3, 3], dtype=np.float64))

        nint = nn - 2

        #
        # Repere de sortie du Tx
        #

        BoO = np.array(np.zeros([3, 3]))
        th = np.arccos(si[0, 2])
        ph = np.arctan2(si[0, 1], si[0, 0])
        eth = np.array([np.cos(
            th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)])
        eph = np.array([-np.sin(ph), np.cos(ph), 0.0])
        Bo0 = np.array([si[0, :], eth, eph]).transpose()

        self.Bo[0, :, :] = Bo0
        #
        # Ray AOD
        #
        self.tht = th
        self.pht = ph
        #
        # Repere d'entree du Rx
        #
        # On impose que le vecteur unitaire d'entree sur
        # le recepteur soit confondu avec le vecteur unitaire
        # de sortie de la derniere interaction
        #
        th = np.arccos(si[nn - 2, 2])
        ph = np.arctan2(si[nn - 2, 1], si[nn - 2, 0])
        eth = np.array([np.cos(
            th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)])
        eph = np.array([-np.sin(ph), np.cos(ph), 0.0])
        Bini = np.array([si[nn - 2, :], eth, eph]).transpose()
        self.Bi[nint, :, :] = Bini
        #
        # Ray AOA
        #
        self.thr = np.pi - th
        self.phr = np.pi + ph
        #
        # loop over the nn-2 interactions
        #

        for l in range(nn - 2):
            typ = self.etype[l + 1]
            nstr = self.nstr[l + 1]

            #
            # Retrieve the interaction vector
            #   normal (R or T) or edge direction (D)
            #
            if nstr > 0:
                # if wall
                if nstr <= L.Ne:
                    vn = L.Gs.node[nstr]['norm']
                else:  # if ceil or floor
                    if nstr == L.Ne + 1:
                        # ceil
                        vn = np.array([0., 0., -1.0])
                    if nstr == L.Ne + 2:
                        # floor
                        vn = np.array([0., 0., 1.0])
            else:  # diffaction
                vn = np.array([0., 0., -1.0])

            # .. todo:: regler le proble de nstr sur rayon 358 : sircut.str
            #try:
            #    ps = np.dot(vn,si[l,:])
            #except:
            #    pdb.set_trace()
            #
            ps = np.dot(vn, si[l, :])
            # Si Reflexion ou Transmission inversion signe normale si.n <0
            #
            if ((typ == 1) | (typ == 2)):
                if (ps < 0):
                    vn = -1.0 * vn
            #
            # Si si et n sont colineaires, le plan d'incidence n'est pas defini
            #
            if (abs(ps) > 1 - 1e-7):
                self.Bi[l + 1, :] = self.Bo[l, :]
            else:
                s_in = si[l, :]
                w = geu.pvecn(s_in, vn)
                v = np.cross(w, s_in)
                M = np.array([s_in, v, w]).transpose()
                self.Bi[l, :, :] = M
                s_out = si[l + 1, :]
                w = geu.pvecn(s_out, vn)
                v = np.cross(w, s_out)
                M = np.array([s_out, v, w]).transpose()
                self.Bo[l + 1, :, :] = M
            #
            # Create interaction
            # see GetRayonPram de tratotud.c
            #
            if typ == 1:  # Reflexion
                #theta  = np.arccos(abs(ps))  'validated'
                theta = self.phii[l + 1]
                siR = lsi[l]
                srR = lsi[l + 1]
                I = IntR(typ, [theta, siR, srR])
            if typ == 2:  # Transmission
                #theta  = np.arccos(abs(ps))  'validated'
                theta = self.phii[l + 1]
                siT = lsi[l]
                srT = lsi[l + 1]
                I = IntT(typ, [theta, siT, srT])
            if typ == 3:  # Diffraction
                theta = self.phii[l + 1]
                thetad = self.phid[l + 1]
                siD = lsi[l]
                sdD = lsi[l + 1]
                beta = self.beta0[l + 1]
                #
                # Attributes for nst <0 of L.Gs
                #
                N = 0
                typed = 0
                I = IntD(typ, [theta, thetad, siD, sdD, beta, N], [typed])
            #
            # add a new interaction in raytud object
            #
            self.gt.inter.append(I)
        #
        # Insert rotation matrices between interaction
        #
        # M = np.zeros((nn-1,2,2))
        #gt.inter.append(I=[)
        for k in range(nn - 1):
            M = np.dot(self.Bo[k, :, 1::].T, self.Bi[k, :, 1::])
            I = IntB(-1, M)
            self.gt.inter.insert(2 * k, I)
            #M[k,:,:] = np.dot(self.Bo[k,:,1::].T,self.Bi[k,:,1::])

        self.gt.ni = 2 * nn - 3
    #def show3(self,bdis=True,bbas=False,col=np.array([1,0,1]),id=0):

    def delay(self):
        """ delay

        Returns
        -------

        delay : float
            delay in ns for each segment of the ray

        """
        pt = self.pt[0:-1, :].T
        ph = self.pt[1::, :].T
        d = pt - ph
        d2 = d * d
        delay = sum(np.sqrt(sum(d * d, axis=0)) / 0.3)
        return(delay)

    def show(self, fig=[], ax=[], col='b', node=False):
        """
        show(ax,vol='b')
        show a Ray projection in 2D

        """
        if fig ==[]:
            fig = plt.gcf()
        if ax==[]:
            ax = fig.gca()


        Nseg = self.nn - 1
        pt = self.pt[0:-1, 0:2].T
        ph = self.pt[1::, 0:2].T
        pz = np.empty((2,))
        pn = np.zeros((2,))
        for i in range(Nseg):
            pz = np.vstack((pz, pt[:, i], ph[:, i], pn))
        m1 = np.array([0, 0, 1])
        mask = np.kron(np.ones((2, Nseg)), m1)
        pzz = pz[1:, :].T
        vertices = np.ma.masked_array(pzz, mask)
        ax.plot(vertices[0, :], vertices[1, :], color=col)
        if node:
            ax.plot(self.pt[:, 0], self.pt[:, 1], 'ok')
        return fig,ax

    def show3(self, bdis=True, bbas=False, bstruc=True, col=np.array([1, 0, 1]), id=0, linewidth=1):
        """ show3(bdis=True,bbas=False,bstruc=True,col=np.array([1,0,1]),id=0)

        Parameters
        ----------

        bdis :
            display boolean - if False return .vect filename
        bbas :
            display local basis
        bstruc :
            display structure
        col  :
            color of the ray
        id   :
            id of the ray
        linewidth :
        """
        filerac = pyu.getlong("ray" + str(id), pstruc['DIRGEOM'])
        _filerac = pyu.getshort(filerac)
        filename_list = filerac + '.list'
        filename_vect = filerac + '.vect'
        try:
            fo = open(filename_vect, "w")
        except:
            raise NameError(filename)

        fo.write("appearance { linewidth %d }\n" % linewidth)

        fo.write("VECT\n")

        fo.write("1 %d 1\n\n" % self.nn)
        fo.write("%d\n" % self.nn)
        fo.write("1\n")
        for i in range(self.nn):
            fo.write("%g %g %g\n" % (self.pt[i, 0], self.pt[i,
                                                            1], self.pt[i, 2]))
        #fo.write("%d %d %d 0\n" % (col[0],col[1],col[2]))
        fo.write("%g %g %g 0\n" % (col[0], col[1], col[2]))
        fo.close()

        #
        # Ajout des bases locales
        #

        fo = open(filename_list, "w")
        fo.write("LIST\n")
        fo.write("{<" + filename_vect + "}\n")
        if (bstruc):
            #fo.write("{<strucTxRx.off}\n")
            fo.write("{<struc.off}\n")
        if (bbas):
            for i in range(self.nn - 1):
                ptb = (self.pt[i + 1, :] + self.pt[i, :]) / 2
                fibi = _filerac + "Bi" + str(i)
                vfibi = GeomVect(fibi)
                colbi = np.array([[1, 0, 0], [1, 0.25, 0.25], [1, 0.5, 0.5]])
                vfibi.geomBase(self.Bi[i, :, :], ptb, colbi, 2)
                fibo = _filerac + "Bo" + str(i)
                vfibo = GeomVect(fibo)
                colbo = np.array([[0, 0, 1], [0.25, 0.25, 1], [0.5, 0.5, 1]])
                vfibo.geomBase(self.Bo[i, :, :], ptb, colbo, 2)
                fo.write("{<" + fibi + ".vect" + "}\n")
                fo.write("{<" + fibo + ".vect" + "}\n")

        filename = filename_list
        fo.close()

        if (bdis):
        #
        # Geomview Visualisation
        #
            chaine = "geomview -nopanel -b 1 1 1 " + filename + \
                " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)


class RayTud(object):
    """ Ray Tud

    Attributes
    ----------

    ni      :
        number of interactions
    inter   :
        Interaction object
    C       :
        Ray Transfer function matrix :  np.array(np.zeros([nf,2,2]),dtype=complex)
    nf      :
        number of frequencies

    Methods
    -------
    info         :
        Ray informations
    eval(f)      :
        Evaluate C tilde matrix for a set of frequencies
    comp(raytud) :
        comparaison de deux rayons Tud

    """
    def __init__(self):
        self.ni = 0

    def info(self):
        print "Number of interactions : ", self.ni
        nbi = self.ni
        for i in range(nbi):
            print "Interaction : ", i
            self.inter[i].info()

    def delay(self):
        """
            calculate delay of the ray
        """
        nbi = self.ni
        d = 0
        for i in range(nbi):
            # LOS
            if self.inter[i].typ == 0:
                d = self.inter[i].dist
            # Reflexion
            if self.inter[i].typ == 1:
                if d == 0:
                    l1 = self.inter[i].si
                    l2 = self.inter[i].sr
                    d = d + l1 + l2
                else:
                    l = self.inter[i].sr
                    d = d + l
            # Transmission
            if self.inter[i].typ == 2:
                if d == 0:
                    l1 = self.inter[i].si
                    l2 = self.inter[i].st
                    d = d + l1 + l2
                else:
                    l = self.inter[i].st
                    d = d + l
            # Diffraction
            if self.inter[i].typ == 3:
                if d == 0:
                    l1 = self.inter[i].si
                    l2 = self.inter[i].sd
                    d = d + l1 + l2
                else:
                    l = self.inter[i].sd
                    d = d + l
        return(d / 0.3)

    def signature(self):
        """
            return ray signature
        """
        Signa = []
        for k in range(self.ni):
            Signa.append(self.inter[k].typ)
        return(Signa)

    def eval(self, fGHz=[2.4]):
        """
        evaluate the field over the ray

        Parameters
        ----------

        Evaluate C tilde matrix for a set of frequencies

        C  :  np.array(np.zeros([nf,2,2]),dtype=complex)

        """
        nf = len(fGHz)
        self.C = np.array(np.zeros([nf, 2, 2]), dtype=complex)
        Co = np.array(np.zeros([nf, 2, 2]), dtype=complex)
        Co[:, 0, 0] = 1
        Co[:, 1, 1] = 1
        #
        # Loop over all the ray interactions
        # ..
        for i in range(self.ni):
            I = self.inter[i]
            CI = I.eval(fGHz)
            for k in range(nf):
                U = np.dot(Co[k, :, :], CI[k, :, :])
                Co[k, :, :] = U
        self.nf = nf
        self.C = Co


class GrRayTud(object):
    """  a cluster of Rays in Tud format 

    Attributes
    ----------
    nray    : number of rayTud
    rayTud  : list of RayTud object

    Methods
    -------

    load(_filetud,sl)  : load from a .tud file
    info(number)

    """
    def __init__(self):
        self.nray = 0

    def dir(self):
        """ list the available file in tuddir


        Returns
        -------
        lfile_s : list
            sorted list of all the .tud file of tuddir

        Notes
        -----
        tuddir is defined in the Project module

        Example
        -------

        >>> from pylayers.antprop.rays import *
        >>> g = GrRay3D()
        >>> l1,l2,l3,l4  = g.dir()

        """

        pathname = tuddir + '/*.tud'
        lfile_l = glob.glob(pathname)
        lfile_tud = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_tud.append(fis)
        lfile_tud.sort()

        pathname = tuddir + '/*.tang'
        lfile_l = glob.glob(pathname)
        lfile_tang = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_tang.append(fis)
        lfile_tang.sort()

        pathname = tuddir + '/*.rang'
        lfile_l = glob.glob(pathname)
        lfile_rang = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_rang.append(fis)
        lfile_rang.sort()

        pathname = tuddir + '/*.tauk'
        lfile_l = glob.glob(pathname)
        lfile_tauk = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_tauk.append(fis)
        lfile_tauk.sort()

        return lfile_tud, lfile_tang, lfile_rang, lfile_tauk

    def delay(self):
        """ return ray delays in an array
        """
        dt = np.array([])
        for nr in range(self.nray):
            d = self.rayTud[nr].delay()
            dt = np.hstack((dt, d))
        return(dt)

    def choose(self):
        """ Choose a tud  file in tuddir
        """
        import tkFileDialog
        FD = tkFileDialog
        filetud = FD.askopenfilename(filetypes=[("Fichiers  ", "*.tud"),
                                     ("All", "*")],
                                     title="Please choose a Tracing file",
                                     initialdir=tuddir)
        _filetud = pyu.getshort(filetud)
        tabc = _filetud.split('_')
        _filestr = tabc[0] + '.str'
        _fileslab = tabc[1] + '.slab'
        _filemat = tabc[2] + '.mat'
        sl = Slab.SlabDB()
        sl.mat = Slab.MatDB()
        sl.mat.load(_filemat)
        sl.load(_fileslab)
        #indoor = Indoor(sl,_filestr)
        self.load(_filetud, sl)

    def save(self, _filetud='ftud', _filetang='ftang', _filerang='frang'):
        """
            save .tud .tang .rang

            Parameters
            ----------

            _filetud  :
            _filetang :
            _filerang :

        """

        filetud = pyu.getlong(_filetud, pstruc['DIRTUD'])
        filetang = pyu.getlong(_filetang,pstruc['DIRTUD'] )
        filerang = pyu.getlong(_filerang,pstruc['DIRTUD'])

        fo = open(filetud, "wb")
        data = stru.pack('i', self.nray)

        for k in range(self.nray):
            rk = self.rayTud[k]
            ni = rk.ni
            dt = stru.pack('i', ni)
            data = data + dt
            for it in rk.inter:
                if (it.typ == -1):
                    dt = stru.pack('i', it.typ)
                    data = data + dt
                    dt = stru.pack('4d', it.M[0, 0], it.M[
                        0, 1], it.M[1, 0], it.M[1, 1])
                    data = data + dt
                elif (it.typ == 0):
                    dt = stru.pack('i', it.typ)
                    data = data + dt
                    dt = stru.pack('d', it.dist)
                    data = data + dt
                elif (it.typ == 1):
                    dt = stru.pack('i', it.typ)
                    data = data + dt
                    dt = stru.pack('3d', it.theta, it.si, it.sr)
                    data = data + dt
                elif (it.typ == 2):
                    dt = stru.pack('i', it.typ)
                    data = data + dt
                    dt = stru.pack('3d', it.theta, it.si, it.st)
                    data = data + dt
                elif (it.typ == 3):
                    dt = stru.pack('i', it.typ)
                    data = data + dt
                    dt = stru.pack('6d', it.theta, it.thetad,
                                   it.si, it.sd, it.beta, it.N)
                    data = data + dt
                    dt = stru.pack('i', it.typed)
                    data = data + dt
                dt = stru.pack('8i', it.datMat[0], it.datMat[1], it.datMat[2], it.datMat[3],
                               it.datMat[4], it.datMat[5], it.datMat[6], it.datMat[7])
                data = data + dt
        fo.write(data)
        fo.close()

        fo = open(filetang, "wb")
        fo.write(stru.pack('i', self.nray))
        for ag in self.tang:
            fo.write(stru.pack('2d', ag[0], ag[1]))
        fo.close()

        fo = open(filerang, "wb")
        fo.write(stru.pack('i', self.nray))
        for ag in self.rang:
            fo.write(stru.pack('2d', ag[0], ag[1]))
        fo.close()

    def load(self, _filetud, _filetang, _filerang, sl):
        """ Load a set of Ray from the PulsRay .tud file

        Parameters
        ----------
        _filename  : tud filename
        _filetang  : tang filename
        _filerang  : rang filename
        sl         : slab database object

        Notes
        -----
        a filename beginning with _ is a short filename

        Examples
        --------

        >>> from pylayers.gis.layout import *
        >>> from pylayers.antprop.rays import *
        >>> L = Layout()
        >>> g3 = GrRay3D()
        >>> l = g3.dir()
        >>> nr = 10
        >>> file0 = l[nr]
        >>> s1 = file0.split('_')
        >>> _filestr = s1[0]+'.str'
        >>> L.loadstr(_filestr)
        >>> g3.load(l[nr],L)
        >>> gt = GrRayTud()
        >>> l1,l2,l3,l4 = gt.dir()
        >>> gt.load(l1[nr],l2[nr],l3[nr],L.sl)
        >>> r30 = g3.ray3d[0]
        >>> rt0 = gt.rayTud[0]


        """

        valerr = False

        filetud = pyu.getlong(_filetud, pstruc['DIRTUD'])
        filetang = pyu.getlong(_filetang,pstruc['DIRTUD'])
        filerang = pyu.getlong(_filerang,pstruc['DIRTUD'])

        fo = open(filetud, "rb")
        data = fo.read()
        fo.close()
        print sl.mat.di

        start = 0
        stop = start + 4
        dt = data[start:stop]
        self.nray = stru.unpack('i', dt)[0]

        self.rayTud = []

        for k in range(self.nray):
            raytud = RayTud()
            start = stop
            stop = start + 4
            dt = data[start:stop]
            #
            # Interaction number over the ray
            #
            # if ni==0  : LOS case no interaction
            #
            nbint = stru.unpack('i', dt)[0]
            raytud.ni = nbint
            Inter = []
            #print "Nombre interactions : ",nbint
            for i in range(nbint):
                start = stop
                stop = start + 4
                dt = data[start:stop]
                caract = stru.unpack('i', dt)[0]
                if (caract == -1):
                    start = stop
                    stop = start + 32
                    dt = data[start:stop]
                    m = stru.unpack('4d', dt)
                    M = np.array([[m[0], m[1]], [m[2], m[3]]])
                    inter = IntB(-1, M)
            #        inter.data = M
                elif (caract == 0):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    dist = stru.unpack('d', dt)
                    inter = IntL(0, dist[0])
            #        inter.data = dist
                elif (caract == 1):
                    start = stop
                    stop = start + 24
                    dt = data[start:stop]
                    datR = stru.unpack('3d', dt)
                    inter = IntR(1, datR)
            #        inter.data = datR
                elif (caract == 2):
                    start = stop
                    stop = start + 24
                    dt = data[start:stop]
                    datT = stru.unpack('3d', dt)
                    inter = IntT(2, datT)
            #        inter.data = datT
                elif (caract == 3):
            #        inter.data = []
                    start = stop
                    stop = start + 48
                    dt = data[start:stop]
                    datD = stru.unpack('6d', dt)
            #        (inter.data).append(datD)
                    start = stop
                    stop = start + 4
                    dt = data[start:stop]
                    typD = stru.unpack('i', dt)
                    inter = IntD(3, datD, typD)
            #        (inter.data).append(typD)

                start = stop
                stop = start + 32
                dt = data[start:stop]
                datMat = stru.unpack('8i', dt)
                inter.datMat = datMat
                l1 = datMat[0]
                c1 = datMat[1]
                r1 = datMat[2]
                s1 = datMat[3]

                l2 = datMat[4]
                c2 = datMat[5]
                r2 = datMat[6]
                s2 = datMat[7]

                # evalfield bug fix
                # slab material inconsistency
                #
                if l1 not in sl.mat.di:
                    valerr = True
                    break

                if ((caract == 1) | (caract == 2) | (caract == 3)):
                    inter.Mat1 = []
                    dim = sl.mat.di
                    #print dim.keys()
                    #print l1
                    #print r1
                    dis = sl.di
                    if s1 == 0:
                        matl = sl.mat[dim[l1]]
                        matr = sl.mat[dim[r1]]
                    else:
                        matl = sl.mat[dim[r1]]
                        matr = sl.mat[dim[l1]]

                    inter.Mat1.append(matl)
                    slab1 = sl[dis[c1]]
                    for i in range(slab1['nbmat']):
                        im = slab1['imat'][i]
                        th = slab1['thickness'][i]
                        matc = sl.mat[dim[im]]
                        #matc.thick = th       # !!! thick existe plus
                        inter.Mat1.append(matc)

                    inter.Mat1.append(matr)
#
#  Mat 2 is used only for diffraction
#
                if (caract == 3):
                    inter.Mat2 = []
                    if s2 == 0:
                        matl = sl.mat[dim[l2]]
                        matr = sl.mat[dim[r2]]
                    else:
                        matl = sl.mat[dim[r2]]
                        matr = sl.mat[dim[l2]]

                    inter.Mat2.append(matl)
                    slab2 = sl[dis[c2]]
                    for i in range(slab2.nbmat):
                        im = slab2.imat[i]
                        th = slab2.thickness[i]
                        matc = sl.mat[dim[im]]
                        matc.thick = th
                        inter.Mat2.append(matc)

                    inter.Mat2.append(matr)

                Inter.append(inter)
            if valerr:
                break
            raytud.inter = Inter
            delay = raytud.delay()
            #
            # this is fix of a bug in pulsray delay discontinuities²
            #
            # impose que les délais soient croissants
            #
            if k == 0:
                self.rayTud.append(raytud)
                delayold = delay
            if (k > 0) & (delay > delayold):
                self.rayTud.append(raytud)
                delayold = delay
        nray = len(self.rayTud)
        self.nray = nray
        # decode the angular files (.tang and .rang)

        self.fail = False
        try:
            fo = open(filetang, "rb")
        except:
            self.fail = True
            print "file ", filetang, " is unreachable"
        if not self.fail:
            nray_tang = stru.unpack('i', fo.read(4))[0]
            buf = fo.read()
            fo.close()
            # coorectif Bug evalfield
            tmp = np.ndarray(shape=(nray_tang, 2), buffer=buf)
            self.tang = tmp[0:nray, :]
        try:
            fo = open(filerang, "rb")
        except:
            self.fail = True
            print "file ", filerang, " is unreachable"

        if not self.fail:
            nray_rang = stru.unpack('i', fo.read(4))[0]
            buf = fo.read()
            fo.close()
            # correctif Bug evalfield
            tmp = np.ndarray(shape=(nray_rang, 2), buffer=buf)
            self.rang = tmp[0:nray, :]

#    def info(self,num=0):
#            r = self.ray[num]
#        nbint = r[0]
#        for l in range(nbint):
#            interaction = r[l+1]
#            print interaction[0]

    def info(self, n=-1):
        """ info

        Parameters
        ----------

        n : int
            ray index (default = -1 all rays)

        """
        print "Nb rayons : ", self.nray
        if n == -1:
            for i in range(self.nray):
                print "Rayon : ", i
                self.rayTud[i].info()
                print "\n"
        else:
            print "rayon no : ", n
            self.rayTud[n].info()


class GrRay3D(object):
    """ A set of Ray3D with the same Tx and Rx

    Attributes
    ----------

    A GrRay is a group of rays sharing the same end points (Tx and Rx)
    n     : number of ray3D
    Tx    : Tx point (3x1)
    Rx    : Rx point (3x1)
    ray3d : list of Ray3D object

    Methods
    -------

    choose()           : choose tra file
    load(_filename,L)  : load a grRay3d from .tra file
    save(_filename)    : load a grRay3d from .tra file
    info(level)        : display information about GrRay3D
    show3(self,bdis=True,bstruc=True,id=0) : geomview visualization
    reciprocal(Gs)     :

    """
    def __init__(self):
        self.n = 0
        self.Tx = np.array([0.0, 0.0, 0.0])
        self.Rx = np.array([0.0, 0.0, 0.0])

    def dir(self):
        """ list the available file in tradir


        Returns
        -------
        lfile_s : list
            sorted list of all the .str file of tradir

        Notes
        -----
        tradir is defined in the Project module

        Example
        -------
        >>> from pylayers.antprop.rays import *
        >>> g = GrRay3D()
        >>> lfile = g.dir()

        .. todo:
            limit to a given filestruc
        """

        pathname = tradir + '/*.tra'
        lfile_l = glob.glob(pathname)
        lfile_s = []
        for fi in lfile_l:
            fis = pyu.getshort(fi)
            lfile_s.append(fis)
        lfile_s.sort()
        return lfile_s

    def info(self, level=1):
        """
        Parameters
        ----------
        level : int
            level of information

        """
        print "Number of Ray : ", self.n
        if level == 1:
            for i in range(self.n):
                self.ray3d[i].info()

    def choose(self):
        """ Choose a Tracing  file in tradir

        """
        import tkFileDialog
        FD = tkFileDialog
        filetra = FD.askopenfilename(filetypes=[("Fichiers Launching ", "*.tra"), ("All", "*")], title="Please choose a Tracing file", initialdir=tradir)
        _filetra = pyu.getshort(filetra)
        tabc = _filetra.split('_')
        _filestr = tabc[0] + '.str'
        _fileslab = tabc[1] + '.slab'
        _filemat = tabc[2] + '.mat'
        sl = Slab.SlabDB()
        sl.mat = Slab.MatDB()
        sl.mat.load(_filemat)
        sl.load(_fileslab)
        self.L = Layout()
        self.L.loadstr(_filestr)

    def load(self, _filename, L):
        """ load a .tra de PulsRay

        Parameters
        ----------
        _filename : PulsRay .tra data format  filename
        L         : Layout object

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.antprop.rays import *
            >>> from pylayers.gis.layout import *
            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> import scipy as sp
            >>> g = GrRay3D()
            >>> lfile = g.dir()
            >>> n = len(lfile)
            >>> k = np.ceil(n*sp.rand()).astype('int')
            >>> file0  = lfile[k]
            >>> s1 = file0.split('_')
            >>> _filestr = s1[0]+'.str'
            >>> L = Layout()
            >>> L.load(_filestr)
            >>> f,a = L.showGs()
            >>> g.load(file0,L)
            >>> g.show(a,np.arange(10))
            >>> plt.show()

        """
        filename = pyu.getlong(_filename, pstruc['DIRTRA'])
        try:
            fo = open(filename, "rb")
        except:
            raise NameError(filename)

        data = fo.read()
        fo.close()
        #
        # decode read data
        #
        start = 0
        stop = start + 1024
        dt = data[start:stop]
        self.flch = dt.replace("\x00", "")

        start = stop
        stop = start + 1024
        dt = data[start:stop]
        self.fpatra = dt.replace("\x00", "")

        start = stop
        stop = start + 1024
        dt = data[start:stop]
        self.fspa = dt.replace("\x00", "")

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Tx[0] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Tx[1] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Tx[2] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Rx[0] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Rx[1] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 8
        dt = data[start:stop]
        self.Rx[2] = stru.unpack('d', dt)[0]

        start = stop
        stop = start + 4
        dt = data[start:stop]
        self.Tracing_exist = stru.unpack('i', dt)[0]

        if (self.Tracing_exist != 0):
            start = stop
            stop = start + 4
            dt = data[start:stop]
            self.n = stru.unpack('i', dt)[0]

            self.ray3d = []

            for i in range(self.n):
                #print "Rayon N° : ",i
                ray3D = Ray3D()
                start = stop
                stop = start + 4
                dt = data[start:stop]
                ray3D.nn = stru.unpack('i', dt)[0]

                ray3D.pt = np.array(np.zeros([ray3D.nn, 3]), dtype=np.float64)

                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    ray3D.pt[k, 0] = stru.unpack('d', dt)[0]

                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    ray3D.pt[k, 1] = stru.unpack('d', dt)[0]

                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    ray3D.pt[k, 2] = stru.unpack('d', dt)[0]

                ray3D.nstr = np.array(np.zeros(ray3D.nn), dtype=int)
                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 4
                    dt = data[start:stop]
                    ray3D.nstr[k] = stru.unpack('i', dt)[0]

                ray3D.deep = np.array(np.zeros(ray3D.nn), dtype=int)
                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 4
                    dt = data[start:stop]
                    ray3D.deep[k] = stru.unpack('i', dt)[0]

                ray3D.beta0 = np.array(np.zeros(ray3D.nn), dtype=np.float64)
                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    ray3D.beta0[k] = stru.unpack('d', dt)[0]

                ray3D.phii = np.array(np.zeros(ray3D.nn), dtype=np.float64)
                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    ray3D.phii[k] = stru.unpack('d', dt)[0]

                ray3D.phid = np.array(np.zeros(ray3D.nn), dtype=np.float64)
                for k in range(ray3D.nn):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    ray3D.phid[k] = stru.unpack('d', dt)[0]

                ray3D.elength = np.array(
                    np.zeros(ray3D.nn - 1), dtype=np.float64)
                for k in range(ray3D.nn - 1):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    ray3D.elength[k] = stru.unpack('d', dt)[0]

                ray3D.etype = np.array(np.zeros(ray3D.nn - 1), dtype=int)
                for k in range(ray3D.nn - 1):
                    start = stop
                    stop = start + 4
                    dt = data[start:stop]
                    ray3D.etype[k] = stru.unpack('i', dt)[0]
                ray3D.etype = np.append(ray3D.etype, [0])
                #
                # local basis creation
                #
                # DEBUG
                #if i < 410:
                #    print ray3D.nstr
                #    n1 = []
                #    for kk in range(len(ray3D.nstr)-2):
                #        p1 = ray3D.pt[kk+1,0:2]
                #        nn = struc.onseg(p1,0.2)
                #        try:
                #            num = nn[0]+1
                #            n1.append(num)
                #        except:
                #            n1.append(-7)
                #    print n1
                iz = np.nonzero(ray3D.nstr == 0)
                niz = len(iz[0])
                self.ray3d.append(ray3D)
                #if niz==2:
                #    ray3D.locbas(L)
                #    self.ray3d.append(ray3D)
                #else:
                #    self.n = self.n-1
                    #print "Problem with ray : ",i
                    #print ray3D.nstr
                    #print ray3D.pt
        fo.close()
#           except:
#               print "Le fichier", filename, "est introuvable"

    def delay(self):
        """ delay

        Returns
        -------
        td  : np.array
            delays
        """
        td = np.array([])
        for n in range(self.n):
            td = np.hstack((td, self.ray3d[n].delay()))
        return(td)

    def save(self, _filename):
        """
        save

        Parameters
        ---------
        _filename : str

        Save a  GrRay3d object in a .tra de PulsRay
        filename : PulsRay .tra data format  filename
        """

        filename = pyu.getlong(_filename, pstruc['DIRTRA'])
        try:
            fo = open(filename, "wb")
        except:
            raise NameError(filename)

        dt_lch = self.flch
        L = len(dt_lch)
        if L < 1024:
            for i in range(1024 - L):
                dt_lch = dt_lch + "\x00"

        dt_patra = self.fpatra
        L = len(dt_patra)
        if L < 1024:
            for i in range(1024 - L):
                dt_patra = dt_patra + "\x00"

        dt_spa = self.fspa
        L = len(dt_spa)
        if L < 1024:
            for i in range(1024 - L):
                dt_spa = dt_spa + "\x00"

        data = dt_lch + dt_patra + dt_spa

        dt = stru.pack('d', self.Tx[0])
        data = data + dt

        dt = stru.pack('d', self.Tx[1])
        data = data + dt

        dt = stru.ack('d', self.Tx[2])
        data = data + dt

        dt = stru.pack('d', self.Rx[0])
        data = data + dt

        dt = stru.pack('d', self.Rx[1])
        data = data + dt

        dt = stru.pack('d', self.Rx[2])
        data = data + dt

        dt = stru.pack('i', self.Tracing_exist)
        data = data + dt

        if (self.Tracing_exist != 0):

            dt = stru.pack('i', self.n)
            data = data + dt

            for i in range(self.n):
                ray3D = self.ray3d[i]

                dt = stru.pack('i', ray3D.nn)
                data_ray = dt

                for k in range(ray3D.nn):
                    dt = stru.pack('d', ray3D.pt[k, 0])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn):
                    dt = stru.pack('d', ray3D.pt[k, 1])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn):
                    dt = stru.pack('d', ray3D.pt[k, 2])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn):
                    dt = stru.pack('i', ray3D.nstr[k])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn):
                    dt = stru.pack('i', ray3D.deep[k])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn):
                    dt = stru.pack('d', ray3D.beta0[k])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn):
                    dt = stru.pack('d', ray3D.phii[k])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn):
                    dt = stru.pack('d', ray3D.phid[k])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn - 1):
                    dt = stru.pack('d', ray3D.elength[k])
                    data_ray = data_ray + dt

                for k in range(ray3D.nn - 1):
                    dt = stru.pack('i', ray3D.etype[k])
                    data_ray = data_ray + dt

                data = data + data_ray
        fo.write(data)
        fo.close()

    def reciprocal(self, L):
        """  Reciprocity channel

        Parameters
        ----------
        L   : Layout object
        return a reciprocal GrRay3D
        """
        G = GrRay3D()

        G.Tx = self.Rx
        print "G.Tx = ", G.Tx
        G.Rx = self.Tx
        print "G.Rx = ", G.Rx
        G.n = self.n

        print "G.n = ", G.n
        print "self.n = ", self.n
        ray3D = []
        for i in range(G.n):
            ray3d = Ray3D()
            ray3d.nn = self.ray3d[i].nn
            nn = ray3d.nn

            ray3d.pht = self.ray3d[i].phr
            ray3d.phr = self.ray3d[i].pht
            ray3d.tht = self.ray3d[i].thr
            ray3d.thr = self.ray3d[i].tht
            ray3d.deep = self.ray3d[i].deep

            ray3d.elength = np.array(np.zeros(ray3d.nn - 1))
            ray3d.etype = np.array(np.zeros(ray3d.nn - 1))

            for j in range(ray3d.nn - 1):
                ray3d.elength[j] = self.ray3d[i].elength[nn - j - 2]
                ray3d.etype[j] = self.ray3d[i].etype[nn - j - 2]

            ray3d.pt = np.array(np.zeros([ray3d.nn, 3]))
            ray3d.phii = np.array(np.zeros(ray3d.nn))
            ray3d.phid = np.array(np.zeros(ray3d.nn))
            ray3d.beta0 = np.array(np.zeros(ray3d.nn))
            ray3d.nstr = np.array(np.zeros(ray3d.nn))

            for j in range(nn):
                ray3d.pt[j] = self.ray3d[i].pt[nn - j - 1]
                ray3d.phid[j] = self.ray3d[i].phid[nn - j - 1]
                ray3d.phii[j] = self.ray3d[i].phii[nn - j - 1]
                ray3d.beta0[j] = self.ray3d[i].beta0[nn - j - 1]

                ray3d.nstr[j] = self.ray3d[i].nstr[nn - j - 1]

            ray3d.locbas(L.Gs)
            ray3D.append(ray3d)
        G.ray3d = ray3D
        return(G)

    def show(self,fig=[], ax=[], rayset=np.array([]), col='b', node=False):
        """

        Parameters
        ----------
        ax     :
            axes object
        rayset :
            set of rays np.array([])
        col  : string
            default  {'b'}
        node : boolean

        """
        if fig ==[]:
            fig = plt.gcf()
        if ax==[]:
            ax = fig.gca()


        for i in rayset:
            r = self.ray3d[i]
            fig,ax=r.show(fig=fig,ax=ax, col=col, node=node)
        return fig, ax


    def show3(self, rayset=np.array([]), bdis=True, bstruc=True, id=0):
        """ 3D show using geomview

        Parameters
        ----------
        rayset : set of index of rays to be displayed
        bdis : display boolean - if False return .vect filename
        id   : id of the grRray

        """
        if (len(rayset) == 0):
            rayset = range(self.n)

        filename = pyu.getlong("grRay" + str(id) + ".list", pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        if bstruc:
            fo.write("{<defstr.off}\n")
            #fo.write("{<strucTxRx.off}\n")
            for i in rayset:
                r = self.ray3d[i]
                col = np.array([0, 0, 1])
                fileray = r.show3(False, False, False, col, i)
                fo.write("{< " + fileray + " }\n")
        fo.close()
        if (bdis):
            chaine = "geomview " + filename + " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)


if __name__ == "__main__":
    doctest.testmod()
