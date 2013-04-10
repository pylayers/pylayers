#!/usr/bin/python
# -*- coding: latin1 -*-
import pdb
import os
import pdb
import glob
import doctest
import ConfigParser
import networkx as nx
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import struct as stru
import pylayers.util.geomutil as geu
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.antprop.slab import *
from pylayers.antprop.interactions import *

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
        # si  = np.array(np.zeros([nn-1,3],dtype=np.float64))
        # lsi = np.array(np.zeros(nn-1,dtype=np.float64))
        #
        # Calcul des vecteurs unitaires si
        # .. todo:: remove this for loop
        #   not tested

        v = self.pt[1:, :] - self.pt[0:-1, :]
        lsi = np.sqrt(np.sum(v * v, axis=1))
        # reshape is required for broadcasting
        lsir = lsi.reshape(nn - 1, 1)
        si = v / lsir
        # for k in range(nn-1):
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
        eth = np.array([np.cos(th) * np.cos(ph),
                        np.cos(th) * np.sin(ph),
                        -np.sin(th)])
        eph = np.array([-np.sin(ph),
                        np.cos(ph),
                        0.0])
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
        # input unitary vector  impose que le vecteur unitaire d'entree sur
        # le recepteur soit confondu avec le vecteur unitaire
        # de sortie de la derniere interaction
        #
        th = np.arccos(si[nn - 2, 2])
        ph = np.arctan2(si[nn - 2, 1], si[nn - 2, 0])
        eth = np.array([np.cos(th) * np.cos(ph),
                        np.cos(th) * np.sin(ph),
                        -np.sin(th)])
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
            # try:
            #    ps = np.dot(vn,si[l,:])
            # except:
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
                # theta  = np.arccos(abs(ps))  'validated'
                theta = self.phii[l + 1]
                siR = lsi[l]
                srR = lsi[l + 1]
                I = IntR([theta, siR, srR])
            if typ == 2:  # Transmission
                # theta  = np.arccos(abs(ps))  'validated'
                theta = self.phii[l + 1]
                siT = lsi[l]
                srT = lsi[l + 1]
                I = IntT([theta, siT, srT])
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
                I = IntD([theta, thetad, siD, sdD, beta, N], [typed])
            #
            # add a new interaction in raytud object
            #
            self.gt.inter.append(I)
        #
        # Insert rotation matrices between interaction
        #
        # M = np.zeros((nn-1,2,2))
        # gt.inter.append(I=[)
        for k in range(nn - 1):
            M = np.dot(self.Bo[k, :, 1::].T, self.Bi[k, :, 1::])
            I = IntB(M)
            self.gt.inter.insert(2 * k, I)
            # M[k,:,:] = np.dot(self.Bo[k,:,1::].T,self.Bi[k,:,1::])

        self.gt.ni = 2 * nn - 3
    # def show3(self,bdis=True,bbas=False,col=np.array([1,0,1]),id=0):

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
        """ show a Ray projection in 2D

        """
        if fig == []:
            fig = plt.gcf()
        if ax == []:
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
        return fig, ax

    def show3(self, _filestr='defstr', bdis=True, bbas=False, bstruc=True, col=np.array([1, 0, 1]), id=0, linewidth=1):
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
        # fo.write("%d %d %d 0\n" % (col[0],col[1],col[2]))
        fo.write("%g %g %g 0\n" % (col[0], col[1], col[2]))
        fo.close()

        #
        # Ajout des bases locales
        #

        fo = open(filename_list, "w")
        fo.write("LIST\n")
        fo.write("{<" + filename_vect + "}\n")
        if (bstruc):
            # fo.write("{<strucTxRx.off}\n")
            fo.write("{<" + _filestr + ".off}\n")
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
        """ returns ray signature
        """
        Signature = []
        for k in range(self.ni):
            Signature.append(self.inter[k].typ)
        return(Signature)

    def eval(self, fGHz=[2.4]):
        """ evaluate the field over the ray

        Parameters
        ----------
        fGHz : ndarray
            frequency axis


        """
        nf = len(fGHz)
        self.C = np.zeros((2, 2), dtype=complex).reshape(1, 2, 2)
        Co = np.eye(2, 2, dtype=complex).reshape(1, 2, 2)
        #
        # Loop over all the ray interactions
        #   + left matrix multiplication
        #   + broadcasting along frequency axis
        #
        for i in range(self.ni):
            I = self.inter[i]
            Ci = I.eval(fGHz)
            Co = np.einsum('kpq,kqr->kpr', Ci, Co)
        self.nf = nf
        self.C = Co

class GrRayTud(dict):
    """  a cluster of Rays in Tud format

    Attributes
    ----------
    nray    : number of rayTud

    rayTud  : list of RayTud object DEPRECATED

    I : class interaction which contrinats all the interactionn of The
        GrRayTud ( see class Interactions help)

    Methods
    -------

    load(_filetud,sl)  : load from a .tud file
    eval()

    """
    def __init__(self):
        self.nray = 0
        # Interactions instance
        self.I = Interactions()
        # dictionnay of interaction legth

    def __repr__(self):
        s = ''
        ni = 0
        nl = 0
        for k in self:
            r = self[k]['rayidx']
            nr = len(r)
            s = s + str(k)+' / '+str(nr)+ ' : '+str(r)+'\n'
            ni = ni + nr*k
            nl = nl + nr*(2*k+1)
        s = s + '-----'+'\n'
        s = s+'ni : '+str(ni)+'\n'
        s = s+'nl : '+str(nl)+'\n'
        return(s)

    def dir(self):
        """ list the available file in tuddir


        Returns
        -------
        lfile_tud : list
            sorted list of all the .tud file of tuddir
        lfile_tang : list
            sorted list of all the .tang file of tuddir
        lfile_rang : list
            sorted list of all the .rang file of tuddir
        lfile_tauk : list
            sorted list of all the .tauk file of tuddir

        Notes
        -----
        tuddir is defined in the Project module

        Example
        -------

        >>> from pylayers.antprop.raysc import *
        >>> g = GrRay3D()

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
        # indoor = Indoor(sl,_filestr)
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
        filetang = pyu.getlong(_filetang, pstruc['DIRTUD'])
        filerang = pyu.getlong(_filerang, pstruc['DIRTUD'])

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
                dt = stru.pack(
                    '8i', it.datMat[0], it.datMat[
                        1], it.datMat[2], it.datMat[3],
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

        .. plot::
            :include-source:

            >>> from pylayers.gis.layout import *
            >>> from pylayers.antprop.raysc import *



        """

        valerr = False

        filetud = pyu.getlong(_filetud, pstruc['DIRTUD'])
        filetang = pyu.getlong(_filetang, pstruc['DIRTUD'])
        filerang = pyu.getlong(_filerang, pstruc['DIRTUD'])

        fo = open(filetud, "rb")
        data = fo.read()
        fo.close()

        start = 0
        stop = start + 4
        dt = data[start:stop]
        self.nray = stru.unpack('i', dt)[0]
        self.rayTud = []

        nimax = 0
        index = 0
        # in order to inialize all type of interaction with the correct
        # frequencies

        B = IntB()
        L = IntL()
        R = IntR()
        T = IntT()
        D = IntD()

        self.mapp = np.zeros((3))
        for k in range(self.nray):
            nir = 0
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
            nbi = nbint
#            raytud.ni = nbint
#            Inter = []

            for i in range(nbint):
                start = stop
                stop = start + 4
                dt = data[start:stop]
                caract = stru.unpack('i', dt)[0]
                ii = i
                if i == 0:
                    decal = False
                    if (caract != -1):
                        # check if first interaction is a IntB.
                        # if not, itthis first interaction is forced to be
                        # an identity matrix
                        M = np.array((1, 0, 0, 1))
                        try:
                            B.stack(M, index)  # inter = IntB(M)
                        except:
                            B = IntB(data=M, idx=index)
                        # because this id matrix, the ray has 1 extra interaction
                        # so we need to
                        # 1.create the correct key dictionnary
                        # 2.remap the index
                        # 3. use another variable loop (ii) incremented

                        nbi = nbi+1

                    try:
                        self[nbi]['rays'] = np.vstack((self[nbi]['rays'],
                                                       np.zeros((1, nbi), dtype=int)))
                        self[nbi]['nbrays'] = self[nbi]['nbrays']+1
                        self[nbi]['rayidx'] = np.hstack((
                            self[nbi]['rayidx'], np.array(([k]))))

                    except:
                        self[nbi] = {}
                        self[nbi]['rays'] = np.zeros((1, nbi), dtype=int)
                        self[nbi]['nbrays'] = 1
                        self[nbi]['rayidx'] = np.array(([k]))

                    try:
                        self.rayidx = np.hstack((self.rayidx, np.array((nbi))))
                    except:
                        self.rayidx = np.array((nbi))

                    if nbi != nbint:
                        self[nbi]['rays'][-1][ii] = index
                        self.mapp = np.vstack((
                            self.mapp, np.array([index, k, -1])))
                        index = index+1
                        decal = True

                # decal = True if an extra B interaction is added at the begining
                # of the ray
                if decal:
                    ii = i + 1
                self[nbi]['rays'][-1][ii] = index
                # B interaction
                if (caract == -1):
                    start = stop
                    stop = start + 32
                    dt = data[start:stop]
                    m = stru.unpack('4d', dt)
                    M = np.array((m[0], m[1], m[2], m[3]))
                    try:
                        B.stack(data=M, idx=index)  # inter = IntB(M)
                    except:
                        print 'except B'
                        B = IntB(data=M, idx=index)
                    index = index+1
#                   M = np.array([[m[0], m[1]], [m[2], m[3]]])
#                    inter = IntB(M)
            #        inter.data = M

                # Los Interaction
                elif (caract == 0):
                    start = stop
                    stop = start + 8
                    dt = data[start:stop]
                    dist = stru.unpack('d', dt)
                    dist = np.array((dist))
                    try:
                        L.stack(data=dist, idx=index)  # inter = IntB(M)
                    except:
                        print 'except L'
                        L = IntL(dist, index)
                    index = index+1
#                    inter = IntL(dist[0])
            #        inter.data = dist

                # Reflexion Interaction
                elif (caract == 1):
                    start = stop
                    stop = start + 24
                    dt = data[start:stop]
                    datR = stru.unpack('3d', dt)
                    dat = np.array((datR[0], datR[1], datR[2]))
                    try:
                        R.stack(data=dat, idx=index)  # inter = IntB(M)
                    except:
                        print 'except R'
                        R = IntR(data=dat, idx=index)
#                    inter = IntR(datR)
#            #        inter.data = datR
                    index = index+1

                # Transmission interaction
                elif (caract == 2):
                    start = stop
                    stop = start + 24
                    dt = data[start:stop]
                    datT = stru.unpack('3d', dt)
                    dat = np.array((datT[0], datT[1], datT[2]))
                    try:
                        T.stack(data=dat, idx=index)
                    except:
                        print 'except T'
                        T = IntT(data=dat, idx=index)
                    index = index+1
#                    inter = IntT(datT)
#            #        inter.data = datT

                # Diffraction interaction
######################## TODO AFTER THIS POINT !!!!!!!!!!!!!!
                elif (caract == 3):
                    print 'catratct 3'
            #        inter.data = []
                    start = stop
                    stop = start + 48
                    dt = data[start:stop]
                    datD = stru.unpack('6d', dt)
#            #        (inter.data).append(datD)
                    start = stop
                    stop = start + 4
                    dt = data[start:stop]
                    typD = stru.unpack('i', dt)
                    index = index+1
#                    inter = IntD(datD, typD)
#            #        (inter.data).append(typD)

                nir = nir+1
                start = stop
                stop = start + 32
                dt = data[start:stop]
                datMat = stru.unpack('8i', dt)
                l1 = datMat[0]
                c1 = datMat[1]
                r1 = datMat[2]
                s1 = datMat[3]

                l2 = datMat[4]
                c2 = datMat[5]
                r2 = datMat[6]
                s2 = datMat[7]

#  if interaction is reflexion or transmission
#                if caract == 1 or caract == 2 or caract == 3:

                    ## find corresponding name
#                    slname = self.I.slab.di[slidx]
#                    # if the material hasn't be evaluated before
#                    if not self.I.slab[slname]['evaluated']:
#                        # evaluate it
#                        print "evaluate",slname
#                        self.I.slab[slname].ev(fGHz=self.I.f, theta=self.I.t)
                ### fill slab index dictionnary for each type of interactions
                if caract == 1:
                    ## read material index

                    slidx = c1
                    ## find corresponding name
                    slname = R.slab.di[slidx]
                    try:
                        R.dusl[slname].append(R.uslidx)
                    except:
                        R.dusl[slname] = [R.uslidx]
                    R.uslidx = R.uslidx + 1
                if caract == 2:
                    slidx = c1
                    slname = T.slab.di[slidx]
                    try:
                        T.dusl[slname].append(T.uslidx)
                    except:
                        T.dusl[slname] = [T.uslidx]
                    T.uslidx = T.uslidx + 1
            # in case of diffraction
                if caract == 3:
                    pass

        ### the total number of interaction( for all rays)
#        self.I.nimax=index
        ### Add all type of interactions into the Interaction class
        self.I.add([B, L, R, T, D])

    def ray(self, r):
        """
            Give the ray number and it returns the index of its interactions
        """
        raypos = np.nonzero(self[self.rayidx[r]]['rayidx'] == r)
        return(self[self.rayidx[r]]['rays'][raypos][0])

    def typ(self, r):
        """
            return the list of interaction type of a given ray
        """

        a = self.ray(r)
        return(self.I.typ[a])

    def eval(self):
        """
        evaluation of Ctilde
        """
        print 'GrRayTUD evaluation'
        if not self.I.evaluated:
            self.I.eval()

        self.Ctilde = np.zeros((self.I.nf, self.nray, 2, 2), dtype=complex)
        self.delays = np.zeros((self.nray))
        self.dis = np.zeros((self.nray))
        nf = self.I.nf  # number of frequence

        for l in self.keys():
            # l stands for the number of interactions
            r = self[l]['nbrays']
            # reshape in order to have a 1D list of insde
            # reshape ray index
            rrl = self[l]['rays'].reshape(r*l)
            # get the corresponding evaluated interactions
            A = self.I.I[:, rrl, :, :].reshape(self.I.nf, r, l, 2, 2)
            alpha = self.I.alpha[rrl].reshape(r, l)
            gamma = self.I.gamma[rrl].reshape(r, l)
            si0 = self.I.si0[rrl].reshape(r, l)
            sout = self.I.sout[rrl].reshape(r, l)

            try:
                del Z
            except:
                pass

            ## loop on the all the interactions of ray with l interactions
            for i in range(1, l-1, 2):

###########################################
#                # Divergence factor D
##                 not yet implementented
###########################################
#                if i == 1:
#                    D0=1./si0[:,1]
#                    rho1=si0[:,1]*alpha[:,i]
#                    rho2=si0[:,1]*alpha[:,i]*gamma[:,i]
#                    D=np.sqrt(
#                     ( (rho1 ) / (rho1 + sout[:,i]) )
#                     *( (rho2) / (rho2 + sout[:,i])))
#                    D=D*D0
#                    rho1=rho1+(sout[:,i]*alpha[:,i])
#                    rho2=rho2+(sout[:,i]*alpha[:,i]*gamma[:,i])
#                    pdb.set_trace()
##                     gerer le loss
#                    if np.isnan(D).any():
#                        p=np.nonzero(np.isnan(D))[0]
#                        D[p]=1./sout[p,1]
#                else :
#                    D=np.sqrt(
#                     ( (rho1 ) / (rho1 + sout[:,i]) )
#                     *( (rho2) / (rho2 + sout[:,i])))

#                    rho1=rho1+(sout[:,i]*alpha[:,i])
#                    rho2=rho2+(sout[:,i]*alpha[:,i]*gamma[:,i])
###########################################

                #  A0  (X dot Y)
                #  |    |     |
                #  v    v     v
                ##########################
                ## B  # I  # B  # I  # B #
                ##########################
                #      \_____/   \______/
                #         |         |
                #       Atmp(i)   Atmp(i+1)
                #
                # Z=Atmp(i) dot Atmp(i+1)

                X = A[:, :, i, :, :]
                Y = A[:, :, i+1, :, :]
                ## Dot product interaction X Basis
                Atmp = np.sum(X[..., :, :, np.newaxis]*Y[
                              ..., np.newaxis, :, :], axis=-2)  # *D[np.newaxis,:,np.newaxis,np.newaxis]
                if i == 1:
                ## First Baspdis added
                    A0 = A[:, :, i-1, :, :]
                    Z = np.sum(A0[..., :, :, np.newaxis]*Atmp[
                               ..., np.newaxis, :, :], axis=-2)
                else:
                    # dot product previous interaction with latest
                    Z = np.sum(Z[..., :, :, np.newaxis]*Atmp[
                               ..., np.newaxis, :, :], axis=-2)

            # fill the C tilde
            self.Ctilde[:, self[l]['rayidx'], :, :] = Z[:, :, :, :]

            # delay computation:
            self[l]['dis'] = self.I.si0[self[l]['rays'][
                :, 1]] + np.sum(self.I.sout[self[l]['rays']], axis=1)

            # Power losses due to distances
            # will be removed once the divergence factor will be implemented
            self.Ctilde[:, self[l]['rayidx'], :, :] = self.Ctilde[:, self[l][
                'rayidx'], :, :]*1./(self[l]['dis'][np.newaxis, :, np.newaxis, np.newaxis])
            self.delays[self[l]['rayidx']] = self[l]['dis']/0.3
            self.dis[self[l]['rayidx']] = self[l]['dis']

    def info(self, r):
        '''
            information for a given ray r

        Attributes
        ----------
        r : a ray number

        '''
        print '-------------------------'
        print 'Informations of ray #', r
        print '-------------------------\n'

        ray = self.ray(r)
        typ = self.typ(r)
        print '{0:5} , {1:4}, {2:10}, {3:7}, {4:10}, {5:10}'.format('Index', 'type', 'material', 'th(rad)', 'alpha', 'gamma2')
        for iidx, i in enumerate(typ):
            if i == 'T' or i == 'R':
                I = getattr(self.I, i)
                for m in I.dusl.keys():
                    midx = I.dusl[m]
                    Iidx = np.array((I.idx))[midx]
                    th = I.data[I.dusl[m], 0]
                    gamma = I.gamma[midx]
                    alpha = I.alpha[midx]
                    for ii, Ii in enumerate(Iidx):
                        if Ii in ray:
                            print '{0:5} , {1:4}, {2:10}, {3:7.2}, {4:10.2}, {5:10.2}'.format(Ii, i, m, th[ii], alpha[ii], gamma[ii])
            else:
                print '{0:5} , {1:4}, {2:10}, {3:7}, {4:10}, {5:10}'.format(ray[iidx], i, '-', '-', '-', '-')

        print '\n----------------------------------------'
        print ' Matrix of ray #', r, 'at f=', self.I.f[0]
        print '----------------------------------------'

        for iidx, i in enumerate(typ):
            print 'interaction #', ray[iidx], 'type:', i
            print self.I.I[0, ray[iidx], :, :]

    def imshowinter(self, r, f, evaluated=False, show=True):
        """
            im show interactions for
            a given ray r and
            a given freq index f
        """
        plt.ion()
        nbinter = self.rayidx[r]
        # find the position of the ray into the dictionnary of legnth of
        # interactions
        pray = np.nonzero(self[nbinter]['rayidx'] == r)
        # get the interaciton idx
        inter = self[nbinter]['rays'][pray][0]
        fig, axs = plt.subplots(
            nrows=2, ncols=nbinter, sharex=True, sharey=True)
        for i, data in enumerate(self.I.I[f, inter]):
            axs[0, i].imshow(np.imag(data), interpolation='none')
            axs[1, i].imshow(np.real(data), interpolation='none')
            axs[1, i].set_xlabel('inter'+self.I.typ[inter[i]])
            if i == 0:
                axs[0, i].set_ylabel('imag')
                axs[1, i].set_ylabel('real')
            if i == nbinter-1:
                axs[0, i].colorbar()

        if self.I.evaluated and evaluated:
            fig2, axs2 = plt.subplots(
                nrows=2, ncols=1, sharex=True, sharey=True)
            data = self[nbinter]['Ctilde'][f, pray, :, :]
            data = data.reshape(2, 2)
            axs2[0].imshow(np.imag(data), interpolation='none')
            axs2[1].imshow(np.real(data), interpolation='none')
            axs2[1].set_xlabel('Ctilde ray number :' + str(r))
            axs2[0].set_ylabel('imag')
            axs2[1].set_ylabel('real')
            if i == nbinter-1:
                axs[0, i].colorbar()

        if show:
            plt.show()

    def get_thetas(self):
        """
        Get all thetas ( incidence angles ) from all computed RayTud

        Returns
        -------
            np.array containg all incidence angles theta from all rays

        """
        th = []
        for r in self.rayTud:
            for idx, i in enumerate(r.inter):
                if not (isinstance(i, IntB) or isinstance(i, IntL)):
                    th.append(r.inter[idx].theta)
        return(np.unique(th))

    def get_mat(self):
        """
        Get all used materials from all computed Raytud

        Returns
        -------
            np.array containg all materials involved in raytracing

        """
        mat = []
        for r in self.rayTud:
            for idx, i in enumerate(r.inter):
                if not (isinstance(i, IntB) or isinstance(i, IntL)):
                    for midx, m in enumerate(r.inter[idx].Mat1):
                        mat.append(r.inter[idx].Mat1[midx]['name'])
        return(np.unique(mat))


#    def info(self, n=-1):
#        """ info
#        Parameters
#        ----------
#        n : int
#            ray index (default = -1 all rays)
#        """
#        print "Nb rayons : ", self.nray
#        if n == -1:
#            for i in range(self.nray):
#                print "Rayon : ", i
#                self.rayTud[i].info()
#                print "\n"
#        else:
#            print "rayon no : ", n
#            self.rayTud[n].info()

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
        >>> from pylayers.antprop.raysc import *
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

            >>> from pylayers.antprop.raysc import *
            >>> from pylayers.util.project import *
            >>> from pylayers.gis.layout import *
            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> import scipy as sp
            >>> g = GrRay3D()
            >>> lfile = g.dir()
            >>> n = len(lfile)
            >>> file0  = lfile[0]
            >>> s1 = file0.split('_')
            >>> _filestr = s1[0]+'.str'
            >>> L = Layout()
            >>> L.load(_filestr)
            >>> f,a = L.showGs()
            >>> g.load(file0,L)
            >>> f1 = g.show(f,a,np.arange(10))
            >>> plt.show()
            >>> f,a = L.showGs()
            >>> f2 = g.show(f,a,np.arange(100))
            >>> plt.show()
            >>> f,a = L.showGs()
            >>> f3 = g.show(f,a,np.arange(300))
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
                # print "Rayon NÂ° : ",i
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
                # if i < 410:
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
                # if niz==2:
                #    ray3D.locbas(L)
                #    self.ray3d.append(ray3D)
                # else:
                #    self.n = self.n-1
                    # print "Problem with ray : ",i
                    # print ray3D.nstr
                    # print ray3D.pt
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
        """ show a cluster of rays 

        Parameters
        ----------
        fig    : figure instance 
        ax     : axes instance  
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
            # fo.write("{<strucTxRx.off}\n")
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
