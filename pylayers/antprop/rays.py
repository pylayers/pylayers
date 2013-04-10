#!/usr/bin/python
# -*- coding: latin1 -*-
import pdb
import os
import pdb
import ConfigParser
import glob
import doctest
import networkx as nx
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import struct as stru
import pylayers.util.geomutil as geu
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.antprop.interactions import *
from pylayers.antprop.slab import *


class Rays(dict):
    def __init__(self, pTx, pRx):
        self.pTx = pTx
        self.pRx = pRx
        self.nray = 0

    def show(self, L):
        """  plot 2D rays within the simulated environment

        Parameters
        ----------

            L : Layout

        """

        fig = plt.figure()
        ax = fig.add_subplot(111)
        L.showGs(fig, ax)
        ax.plot(self.pTx[0], self.pTx[1], 'or')
        ax.plot(self.pRx[0], self.pRx[1], 'og')
        for i in self.keys():
            for j in range(len(self[i]['pt'][0, 0, :])):
                ray = np.hstack((self.pTx[0:2].reshape((2, 1)),
                                 np.hstack((self[i]['pt'][0:2, :, j],
                                            self.pRx[0:2].reshape((2, 1))))
                                 ))
                ax.plot(ray[0, :], ray[1, :], alpha=0.6, linewidth=1.)

    def mirror(self, H=3, N=1):
        """ mirror

        Parameters
        ----------

        H : float
            ceil height (default 3m)

        N : int
            handle the number of mirror reflexions

        """
        km = np.arange(-N+1, N+1, 1)
        kp = np.arange(-N, N+1, 1)
        ht = self.pTx[2]
        hr = self.pRx[2]
        zkp = 2*kp*H + ht
        zkm = 2*km*H - ht

        d = {}

        for zm in zkm:
            if zm < 0:
                bup = H
                pas = H
                km = int(np.ceil(zm/H))
            else:
                bup = 0
                pas = -H
                km = int(np.floor(zm/H))
            thrm = np.arange(km*H, bup, pas)
            d[zm] = abs(thrm-zm)/abs(hr-zm)

        for zp in zkp:
            if zp < 0:
                bup = H
                pas = H
                kp = int(np.ceil(zp/H))
            else:
                bup = 0
                pas = -H
                kp = int(np.floor(zp/H))
            thrp = np.arange(kp*H, bup, pas)
            d[zp] = alphap = abs(thrp-zp)/abs(hr-zp)
            # print "zp",zp
            # print "kp",kp
            # print "thrp",thrp
            # print "alphap",d[zp]

        return(d)

    def to3D(self, H=3, N=1):
        """ transform 2D ray to 3D ray

        Parameters
        ----------

        H : float
            ceil height (default 3m)

        N : int
            handle the number of mirror reflexions

        Notes
        -----


        """

        tx = self.pTx
        rx = self.pRx
        #
        # Phase 1 : calculate Tx images height and vertical parameterization
        #

        d = self.mirror(H=H, N=N)

        #
        # Phase 2 : calculate 2D parameterization
        #

        for i in self:
            pts = self[i]['pt'][0:2, :, :]
            sig = self[i]['sig']
            t = self.pTx[0:2].reshape((
                2, 1, 1)) * np.ones((1, 1, len(pts[0, 0, :])))
            r = self.pRx[0:2].reshape((
                2, 1, 1)) * np.ones((1, 1, len(pts[0, 0, :])))
            pts1 = np.hstack((t, np.hstack((pts, r))))
            si1 = pts1[:, 1:, :] - pts1[:, :-1, :]
            si = np.sqrt(np.sum(si1 * si1, axis=0))
            al1 = np.cumsum(si, axis=0)
            self[i]['alpha'] = np.zeros(np.shape(si[:-1, :]))

            for j in range(len(self[i]['alpha'][:, 0])):
                self[i]['alpha'][j, :] = np.sum(si[
                                                0:j+1, :], axis=0)/np.sum(si, axis=0)
                self[i]['pt'][2, j, :] = tx[2] + self[
                    i]['alpha'][j, :] * (rx[2] - tx[2])

        #
        #  Phase 3 : Initialize 3D rays dictionnary
        #
        r3d = Rays(tx, rx)
        r3d.is3D = True

        #
        # Phase 4 : Fill 3D rays information
        #
        # Two nested loops
        #
        #      for all interaction group
        #          for all type of 3D rays
        #             extension
        #             sort
        #             coordinates as a function of parameter
        #
        for k in self:   # for all interaction group k
            # k = int(k)
            Nrayk = np.shape(self[k]['alpha'])[
                1]  # Number of rays in interaction group k
            a1 = self[k]['alpha']           # get  2D parameterization
            sig = self[k]['sig']            # get  2D signature
            a1 = np.concatenate((np.zeros((1, Nrayk)), a1, np.ones((
                1, Nrayk))))    # add parameterization of tx and rx (0,1)
            sig = np.hstack((np.zeros((2, 1, Nrayk), dtype=int),
                             sig,
                             np.zeros((2, 1, Nrayk), dtype=int)))  # add signature of Tx and Rx (0,0))
            Tx = tx.reshape(3, 1, 1)*np.ones((1, 1, Nrayk))
            Rx = rx.reshape(3, 1, 1)*np.ones((1, 1, Nrayk))
            pte = self[k]['pt']             # ndim x k x Nrayk
            pte = np.hstack((Tx, pte, Rx))       # ndim x k+2 x Nrayk
            for l in d:                          # for each vertical pattern (C,F,CF,FC,....)
                Nint = len(d[
                           l])                 # number of additional interaction
                if Nint > 0:                     # if new interaction ==> need extension
                    a1e = np.concatenate((a1, d[l].reshape(len(d[l]), 1)*np.ones(
                        (1, Nrayk))))  # extended old parameterization
                    ks = np.argsort(a1e, axis=0)  # get sorted indices
                    a1es = np.sort(
                        a1e, axis=0)  # sorted extended parameterization
                    ptee = np.hstack((pte, np.zeros((
                        3, Nint, Nrayk))))     # ndim x (Nint+k+2) x Nrayk
                    if l < 0:
                        u = np.mod(range(Nint), 2)
                    else:
                        u = 1 - np.mod(range(Nint), 2)

                    esigs = np.zeros((1, Nint, Nrayk), dtype=int)
                    esigi = (u+4).reshape(1, Nint, 1)*np.ones(
                        (1, 1, Nrayk), dtype=int)
                    esig = np.vstack((esigs, esigi))
                    # sige   = np.hstack((sig,np.zeros((2,Nint,Nrayk))))
                    # # 2 x (Nint+k+2) x Nrayk
                    sige = np.hstack((sig, esig))
                                       # 2 x (Nint+k+2) x Nrayk
                    ptees = ptee[:, ks, range(Nrayk)]     # sorted points
                    siges = sige[:, ks, range(Nrayk)]     # sorted signature
                    iint_f, iray_f = np.where(siges[
                                              1, :] == 4)  # floor interaction
                    iint_c, iray_c = np.where(siges[
                                              1, :] == 5)  # ceil interaction

                    coeff_f = (a1es[iint_f, iray_f]-a1es[iint_f-1, iray_f])/(
                        a1es[iint_f+1, iray_f]-a1es[iint_f-1, iray_f])
                    coeff_c = (a1es[iint_c, iray_c]-a1es[iint_c-1, iray_c])/(
                        a1es[iint_c+1, iray_c]-a1es[iint_c-1, iray_c])
                    ptees[0:2, iint_f, iray_f] = ptees[0:2, iint_f-1, iray_f] + coeff_f*(
                        ptees[0:2, iint_f+1, iray_f]-ptees[0:2, iint_f-1, iray_f])
                    # ptees[2,iint_f,iray_f]   = 0
                    ptees[0:2, iint_c, iray_c] = ptees[0:2, iint_c-1, iray_c] + coeff_c*(
                        ptees[0:2, iint_c+1, iray_c]-ptees[0:2, iint_c-1, iray_c])
                    # ptees[2,iint_c,iray_c]   = H
                    z = np.mod(l+a1es*(rx[2]-l), 2*H)
                    pz = np.where(z > H)
                    z[pz] = 2*H-z[pz]
                    ptees[2, :] = z
                else:
                    a1es = a1                        # recopy old 2D parameterization (no extension)
                    ks = np.argsort(a1es, axis=0)
                    ptees = pte
                    siges = sig
                try:
                    # r3d[k+Nint]['alpha'] = np.hstack((r3d[k+Nint]['alpha'],a1es))
                    # r3d[k+Nint]['ks'] = np.hstack((r3d[k+Nint]['ks'],ks))
                    r3d[k+Nint]['pt'] = np.dstack((r3d[k+Nint]['pt'], ptees))
                    r3d[k+Nint]['sig'] = np.dstack((r3d[k+Nint]['sig'], siges))
                except:
                    r3d[k+Nint] = {}
                    # r3d[k+Nint]['alpha'] = a1es
                    # r3d[k+Nint]['ks'] = ks
                    r3d[k+Nint]['pt'] = ptees
                    r3d[k+Nint]['sig'] = siges
        return(r3d)

    def locbas(self, L):
        """
        Parameters
        ----------

        L : Layout

        """

        #
        # extract normal in np.array
        #

        # nsegment x 3
        norm = np.array(nx.get_node_attributes(
            L.Gs, 'norm').values())

        # nsegment x k
        key = np.array(nx.get_node_attributes(
            L.Gs, 'norm').keys())

        nmax = max(L.Gs.node.keys())
        mapping = np.zeros(nmax+1, dtype=int)
        mapping[key] = np.arange(len(key), dtype=int)

        for k in self:

            nstr = self[k]['sig'][0, 1:-1, :]      # nint x nray
            ityp = self[k]['sig'][1, 1:-1, :]      # nint x nray

            nray = np.shape(nstr)[1]

            uwall = np.where((ityp == 1) | (ityp == 2))
            udiff = np.where((ityp == 3))
            ufloor = np.where((ityp == 4))
            uceil = np.where((ityp == 5))

            nstrwall = nstr[uwall[0], uwall[1]]   # nstr of walls
            self[k]['nstrwall'] = nstrwall       # store

            self[k]['norm'] = np.zeros((3, k, nray))   # 3 x int x nray

            #
            # Warning : The following commented line assumes that all the segment number are contiguous
            # self[k]['norm'][:,uwall[0],uwall[1]] = norm[nstrwall-1,:].T
            #

            self[k]['norm'][:, uwall[0], uwall[
                1]] = norm[mapping[nstrwall], :].T
            self[k]['norm'][2, ufloor[0], ufloor[1]] = np.ones(len(ufloor[0]))
            self[k]['norm'][2, uceil[0], uceil[1]] = -np.ones(len(uceil[0]))

            v = self[k]['pt'][:, 1:, :]-self[k]['pt'][:, 0:-1, :]
            lsi = np.sqrt(np.sum(v*v, axis=0))
            si = v/lsi             # ndim , nint - 1 , nray
            self[k]['si'] = si

            vn = self[k]['norm']
            s_in = si[:, 0:-1, :]
            s_out = si[:, 1:, :]

            #
            # AOD (rad)
            #

            # th : ,r
            th = np.arccos(si[2, 0, :])

            # ph : ,r
            ph = np.arctan2(si[1, 0, :], si[0, 0, :])

            # aod : 2 x r  (radians)
            self[k]['aod'] = np.vstack((th, ph))

            # eth : 3 x r
            eth = np.array([np.cos(th) * np.cos(ph),
                           np.cos(th) * np.sin(ph),
                            -np.sin(th)])
            # eph : 3 x r
            eph = np.array([-np.sin(ph),
                            np.cos(ph),
                            np.zeros(len(ph))])

            # Bo0 : 3 x 2 x r
            Bo0 = np.concatenate((eth[:, np.newaxis, :],
                                  eph[:, np.newaxis, :]), axis=1)

            self[
                k]['Bo0'] = np.concatenate((si[:, 0, np.newaxis, :], eth[:, np.newaxis, :],
                                            eph[:, np.newaxis, :]), axis=1)

            #
            # scalar product si . norm
            #

            scpr = np.sum(vn*s_in, axis=0)
            self[k]['scpr'] = scpr
            self[k]['theta'] = np.arccos(abs(scpr))  # *180/np.pi

            #
            # Warning need to handle singular case when s_in//vn
            #
            w = np.cross(s_in, vn, axisa=0, axisb=0, axisc=0)
            wn = w/np.sqrt(np.sum(w*w, axis=0))
            v = np.cross(wn, s_in, axisa=0, axisb=0, axisc=0)

            es_in = np.expand_dims(s_in, axis=1)
            ew = np.expand_dims(wn, axis=1)
            ev = np.expand_dims(v, axis=1)

            #  Bi 3 x 2 x i x r
            Bi = np.concatenate((ew, ev), axis=1)
            # self[k]['Bi'] = np.concatenate((es_in,ew,ev),axis=1)

            w = np.cross(s_out, vn, axisa=0, axisb=0, axisc=0)
            wn = w/np.sqrt(np.sum(w*w, axis=0))
            v = np.cross(wn, s_out, axisa=0, axisb=0, axisc=0)

            es_out = np.expand_dims(s_out, axis=1)
            ew = np.expand_dims(wn, axis=1)
            ev = np.expand_dims(v, axis=1)

            #  Bi 3 x 2 x i x r
            Bo = np.concatenate((ew, ev), axis=1)
            # self[k]['Bo'] = np.concatenate((es_out,ew,ev),axis=1)

            #
            # AOA (rad)
            #

            # th : ,r
            th = np.arccos(si[2, -1, :])

            # th : ,r
            ph = np.arctan2(si[1, -1, :], si[0, -1, :])

            # aoa : 2 x r  (radians)
            self[k]['aoa'] = np.vstack((th, ph))
            eth = np.array([np.cos(th) * np.cos(ph),
                           np.cos(th) * np.sin(ph),
                            -np.sin(th)])
            eph = np.array([-np.sin(ph),
                            np.cos(ph),
                            np.zeros(len(ph))])
            # Bo0 : 3 x 2 x r
            BiN = np.concatenate((eth[:, np.newaxis, :],
                                  eph[:, np.newaxis, :]), axis=1)

            # self[k]['BiN'] = np.concatenate((si[:,-1,np.newaxis,:],eth[:,np.newaxis,:],
            #                                    eph[:,np.newaxis,:]),axis=1)

            #
            # pasting (Bo0,B,BiN)
            #
            # B : 3 x 2 x i x r

            Bo = np.concatenate((Bo0[:, :, np.newaxis, :], Bo), axis=2)
            Bi = np.concatenate((Bi, BiN[:, :, np.newaxis, :]), axis=2)

            # M : 2 x 2 x i x r

            self[k]['M'] = np.einsum('xv...,xw...->vw...', Bo, Bi)

            # BiN = np.array([si[:,-1,:], eth, eph])    # ndim x 3 x Nray
            # self[k]['BiN']=BiN
            # self[k]['B']=np.sum(self[k]['Bi'][:2,:2,np.newaxis]*self[k]['Bo'][np.newaxis,:2,:2],axis=1)

    def fillinter(self, L):
        """  docstring for fillinter

        Parameters
        ----------

        L : Layout
        """
        I = Interactions()
        B = IntB()
        Los = IntL()
        R = IntR()
        T = IntT()
        D = IntD()
        idx = np.array(())
        idxts = 0
        nbrayt = 0

        # Transform dictionnary of slab name to array
        slv = nx.get_node_attributes(L.Gs, "name").values()
        slk = nx.get_node_attributes(L.Gs, "name").keys()

        # find all material used in simulation
        uslv = np.unique(slv)
        uslv = np.hstack((uslv, np.array(('CEIL', 'FLOOR'))))

#       create reverse dictionnary with all material as a key
#       and associated point/segment as a value

        dsla = {}
        for s in uslv:
            dsla[s] = np.where(s == np.array(slv))[0]

        nmax = max(L.Gs.node.keys())
        sla = np.zeros((nmax+1), dtype='S20')

        # array type str with more than 1 character
        # warning use zeros instead of empty because slab zero
        # is virtually used before assigning correct slab to ceil and floor

        # slab is now an array of string.
        # each value of Gs node is the index of the coresponding slab
        sla[slk] = np.array(slv)

        R.dusl = dict.fromkeys(uslv, np.array((), dtype=int))
        T.dusl = dict.fromkeys(uslv, np.array((), dtype=int))

        tsl = np.array(())
        rsl = np.array(())
        for k in self:

            uR = uT = uD = uRf = uRc = 0.
            nstr = self[k]['sig'][0, 1:-1, :]      # nint x nray
            ityp = self[k]['sig'][1, 1:-1, :]      # nint x nray
            theta = self[k]['theta']
            si = self[k]['si']

            ## flatten information
            ######################

            # reshape nstr in order to be flat (1 dimension)
            nstrf = np.reshape(nstr, nstr.size)
            # flatten ityp
            itypf = ityp.reshape(ityp.size)
            thetaf = theta.reshape(theta.size)
            sif = si[0, :, :].reshape(si[0, :, :].size)
            b = self[k]['M'].reshape(2, 2, ityp.size)

            ## index creation
            ##################
            # create index for retrieve interactions

            idxts = idxts + idx.size  # total size idx

            # idx is an abolute index of the interaction position
            idx = idxts + np.arange(ityp.size).reshape(np.shape(ityp)).T

            nbray = np.shape(idx)[0]

            self[k]['rays'] = idx
            self[k]['nbrays'] = nbray
            self[k]['rayidx'] = nbrayt + np.arange(nbray)
            nbrayt = nbrayt + nbray

            self.nray = self.nray + self[k]['nbrays']
            idxf = idx.reshape(idx.size)

            ## find used slab
            ##################

            # find slab type for the rnstr

            sl = sla[nstrf]

            # seek interactions position
            ############################
            uR = np.where((itypf == 1))[0]
            uT = np.where((itypf == 2))[0]
            uD = np.where((itypf == 3))[0]
            uRf = np.where((itypf == 4))[0]
            uRc = np.where((itypf == 5))[0]

            # assign floor and ceil slab
            ############################

            # WARNING
            # in future version floor and ceil could be different for each cycle.
            # this information would be directly obtained from L.Gs
            # then the two following lines would have to be  modified

            sl[uRf] = 'FLOOR'
            sl[uRc] = 'CEIL'

#            # Fill the used slab
#            ####################

            tsl = np.hstack((tsl, sl[uT]))
            rsl = np.hstack((rsl, sl[uR], sl[uRf], sl[uRc]))
##            for s in uslv:
##
##                T.dusl[s]=np.hstack((T.dusl[s],len(T.idx) + np.where(sl[uT]==s)[0]))
##                R.dusl[s]=np.hstack((R.dusl[s],len(R.idx) + np.where(sl[uR]==s)[0]))
##            R.dusl['FLOOR']=np.hstack((R.dusl['FLOOR'],len(R.idx)+len(uR) + np.where(sl[uRf]=='FLOOR')[0]))
# R.dusl['CEIL']=np.hstack((R.dusl['CEIL'],len(R.idx)+len(uR)+len(uRf) +
# np.where(sl[uRc]=='CEIL')[0]))

            # Basis
            # Hugr issue with B index
            # Friedman version Bs was entering in the index
            # maybe B can have the same index that interactions
            # but this must be manage when evaluation of CIR is made
            B.stack(data=b.T, idx=idxf)

##            #Reflexion
##            ##########
##          # wall reflexion
            R.stack(data=np.array((thetaf[uR], sif[uR], sif[uR+1])).T,
                    idx=idxf[uR])
            # floor reflexion
            R.stack(data=np.array((thetaf[uRf], sif[uRf], sif[uRf+1])).T,
                    idx=idxf[uRf])
            # ceil reflexion
            R.stack(data=np.array((thetaf[uRc], sif[uRc], sif[uRc+1])).T,
                    idx=idxf[uRc])

###         sl[idxf[uT]]

            # Transmision
            ############
            T.stack(data=np.array((thetaf[uT], sif[
                    uT], sif[uT+1])).T, idx=idxf[uT])

        T.create_dusl(tsl)
        R.create_dusl(rsl)
        self.I = I
        self.I.add([T, R])
        self.B = B

    def eval(self):
        """docstring for eval"""

        print 'GrRayTUD evaluation'
        if not self.I.evaluated:
            self.I.eval()
            self.B.eval()

        self.Ctilde = np.zeros((self.I.nf, self.nray, 2, 2), dtype=complex)
        self.delays = np.zeros((self.nray))
        self.dis = np.zeros((self.nray))
        nf = self.I.nf  # number of frequence
        # loop on interations
        for l in self:
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

        # To be corrected in a futeur version
        self.Ctilde = np.swapaxes(self.Ctilde, 1, 0)

    def signature(self, L):
        """
        """
        sig = Signatures(L, self.pTx, self.pRx)
        for k in self:
            sig[k] = self[k]['sig']
        return(sig)

    def show3d(self,
               ray,
               bdis=True,
               bbas=False,
               bstruc=True,
               col=np.array([1, 0, 1]),
               id=0,
               linewidth=1):
        """
        plot a 3D ray
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

        fo.write("1 %d 1\n\n" % len(ray[0, :]))
        fo.write("%d\n" % len(ray[0, :]))
        fo.write("1\n")
        for i in range(len(ray[0, :])):
            fo.write("%g %g %g\n" % (ray[0, i], ray[1, i],
                                     ray[2, i]))
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

    def show3(self, bdis=True, bstruc=True, id=0,
              strucname='defstr', ilist=[], raylist=[]):
        """ plot 3D rays within the simulated environment

        Parameters
        ----------

            raysarr: numpy.ndarray

        """
        if ilist == []:
            ilist = self.keys()
        pTx = self.pTx.reshape((3, 1))
        pRx = self.pRx.reshape((3, 1))
        filename = pyu.getlong("grRay" + str(id) + ".list", pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        if bstruc:
            fo.write("{<"+strucname+".off}\n")
            # fo.write("{<strucTxRx.off}\n")
            k = 0
            for i in ilist:
                if raylist == []:
                    rlist = range(np.shape(self[i]['pt'])[2])
                else:
                    rlist = raylist
                for j in rlist:
                    ray = np.hstack((pTx,
                                     np.hstack((self[i]['pt'][:, :, j], pRx))))
                    # ray = rays[i]['pt'][:,:,j]
                    col = np.array([2, 0, 1])
                    # print ray
                    fileray = self.show3d(ray=ray, bdis=False,
                                          bstruc=False, col=col, id=k)
                    k += 1
                    fo.write("{< " + fileray + " }\n")
        fo.close()
        if (bdis):
            chaine = "geomview " + filename + " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)


if __name__ == "__main__":
    doctest.testmod()
