# -*- coding:Utf-8 -*-
#####################################################################
#This file is part of RGPA.

#Foobar is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#Foobar is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------------------
#authors :
#Nicolas AMIOT          : nicolas.amiot@univ-rennes1.fr
#Bernard UGUEN          : bernard.uguen@univ-rennes1.fr
#Mohamed LAARAIEDH      : mohamed.laaraiedh@univ-rennes1.fr
#####################################################################
import numpy as np
import scipy as sp
from copy import copy
from pylayers.location.geometric.util.boxn import *
from pylayers.location.geometric.constraints.constraint import *


class TDOA(Constraint):
    def __init__(self, value=45, std=4, vcw=3, p=np.array([[0, 0, 0], [10, 10, 10]])):
        """



        Attributes
        ----------

        tdoa = (|F1M| - |F2M|)/0.3 = (|p1M| - |p2M|)/0.3    (ns)


        value : Time difference of Arrival  (ns)    30 ns  = 9 m
        std   : Standard deviation on tdoa  (ns)    1 ns   = 30cm
        vcw   : constraint with factor


        Methods
        -------

        vrai(deltar)  : check constraint validity

        """
        Constraint.__init__(self, 'TDOA', p)

        #
        # (vn,wn,tn) triedre orthonormal
        #
    #       self.Dmax = 25 # limit of tdoa box in meter
        ##
        self.tdoa_axes(p)

        self.f = self.nv / 2
        self.Dmax = self.nv
        self.value = min(value, 2 * self.f / 0.3)
        self.value = max(self.value, -2 * self.f / 0.3)
        self.std = std
        self.vcw = vcw
        self.drange = self.value * 0.3
        self.sstd = self.std * 0.3
        self.tdoa_box(vcw)
#               if self.ndim == 3:
#                       BOUND1 = np.array([0.0,0.0,-2.0])
#                       BOUND2 = np.array([20.0,20.0,2.0])
#                       box         = BoxN(np.vstack((BOUND1,BOUND2)),ndim=np.shape(self.p)[1])
#               else:
#                       BOUND1 = np.array([0.0,0.0])
#                       BOUND2 = np.array([20.0,20.0])
#                       box         = BoxN(np.vstack((BOUND1,BOUND2)),ndim=np.shape(self.p)[1])

#               self.lbox    = LBoxN([box],ndim=np.shape(self.p)[1])

        self.annulus_bound()
        self.Id = copy.copy(self.C_Id)
        Constraint.C_Id = Constraint.C_Id + \
            1   # constraint counter is incremented

    def tdoa_axes(self, p):
        """
        tdoa_axes(self,p) : triedre [vn,wn,tn], support of the contraint
        """
        #
        # Dmax
        #

        #
        # La boite du tdoa est tronquee
        #
        self.F1 = p[0, :]
        self.F2 = p[1, :]
        #
        #
        #
        v = self.F2 - self.F1
        self.nv = np.sqrt(np.dot(v, v))
        vn = v / (self.nv * 1.0)

        if self.ndim > 2:
            if np.abs(v[2]) < 0.9:
                w = np.array([v[1], -v[0], 0])
            else:
                w = np.array([v[2], 0, -v[0]])
            nw = np.sqrt(np.dot(w, w))
            wn = w / (nw * 1.0)
            tn = np.cross(vn, wn)
            self.triedre = [wn, tn, vn]  # [z,x,y]
        else:
            w = np.array([v[1], -v[0]])
            nw = np.sqrt(np.dot(w, w))
            wn = w / (nw * 1.0)
            self.triedre = [wn, vn]

    def tdoa_box(self, vcw):
        """
        tdoa_box(self,vcw) : create the inclusive box for a given vcw
        """

        if self.ndim == 3:
            wn = self.triedre[0]
            tn = self.triedre[1]
            vn = self.triedre[2]

        if self.ndim == 2:
            wn = self.triedre[0]
            vn = self.triedre[1]

        eps = vcw * self.sstd
        delta = self.drange
        deltap = min(delta + eps, self.nv)
        deltam = max(delta - eps, -self.nv)
        c = delta / 2
        cp = deltap / 2
        cm = deltam / 2
        arge = self.f ** 2 - c ** 2
        argep = self.f ** 2 - cp ** 2
        argem = self.f ** 2 - cm ** 2
        e = np.sqrt(arge)
        ep = np.sqrt(argep)
        em = np.sqrt(argem)

        if cp < 0:
            pp = self.F1 + (self.f + cp) * vn
        else:
            if ep > 0:
                offset = (cp * np.sqrt((self.Dmax / ep) ** 2 + 1) - self.f)
                #print "ep >0 : offset ",offset
            else:
                offset = -self.Dmax
            pp = self.F2 + offset * vn

        if cm < 0:
            if em > 0:
                offset = (cm * np.sqrt((self.Dmax / em) ** 2 + 1) + self.f)
                #print "em >0 : offset ",offset
            else:
                offset = self.Dmax
            pm = self.F1 + offset * vn
        else:
            pm = self.F2 - (self.f - cm) * vn

        if self.ndim == 3:
            p1 = pp + self.Dmax * wn - self.Dmax * tn
            p2 = pp + self.Dmax * wn + self.Dmax * tn
            p3 = pp - self.Dmax * wn + self.Dmax * tn
            p4 = pp - self.Dmax * wn - self.Dmax * tn
            p5 = pm + self.Dmax * wn - self.Dmax * tn
            p6 = pm + self.Dmax * wn + self.Dmax * tn
            p7 = pm - self.Dmax * wn + self.Dmax * tn
            p8 = pm - self.Dmax * wn - self.Dmax * tn
            pquad = np.vstack((p1, p2, p3, p4, p5, p6, p7, p8))

        if self.ndim == 2:
            p1 = pp + self.Dmax * wn
            p2 = pp - self.Dmax * wn
            p3 = pm + self.Dmax * wn
            p4 = pm - self.Dmax * wn
            pquad = np.vstack((p1, p2, p3, p4))

        imin = np.min(pquad, axis=0)
        imax = np.max(pquad, axis=0)

        self.ep = ep
        self.em = em
        self.cp = cp
        self.cm = cm

        self.lbox = LBoxN(
            [BoxN(np.vstack((imin, imax)), ndim=np.shape(self.p)[1])])

    def annulus_bound(self):
        """
        annulus_bound():
        Compute the minimum and maximum distance of the enclosing annulus of the constraint for a given self.vcw
        """
        if self.value > 0:
            self.cmin = self.drange - self.vcw * self.sstd
            self.cmax = self.drange + self.vcw * self.sstd
        else:
            self.cmin = self.drange + self.vcw * self.sstd
            self.cmax = self.drange - self.vcw * self.sstd

        self.mean = (self.cmin + self.cmax) / 2

    def repart(self, DD):
        return(1. / (self.sstd * np.sqrt(2 * np.pi)) * np.exp(-(DD - self.mean) ** 2 / (2 * self.sstd ** 2)))

    def rescale(self, vcw):
        """
        rescale(vcw) :  rescale constraint with vcw factor
        """
        self.vcw = vcw
        #print self.vcw
#               pdb.set_trace()
        self.tdoa_box(self.vcw)
        print 'TDOA', self.vcw
        #self.estvol() <= TO BE DONE IN TDOA

    def inclusive(self, b):
        """
        inclusive(b) : A box b is inclusive for the constraint if self.p is included in the box

        """
        if b.inbox(self.p):
            return True
        else:
            return False

    def valid(self, b):
        """
        valid(b) : check if box b is valid for the given constraint

        A box is valid if it not not valid

        A box is not valid if all distances are greater than rangemax
                           or all distances are less than rangemin
        """

        v = b.bd2coord()
        P0 = np.outer(np.ones(len(v)), self.p[0, :])
        P1 = np.outer(np.ones(len(v)), self.p[1, :])
        F1v = np.sqrt(np.sum((P0 - v) * (P0 - v), axis=1))
        F2v = np.sqrt(np.sum((P1 - v) * (P1 - v), axis=1))
        D = (F1v - F2v)

        if self.value > 0:
            DDcmin = sum(D >= self.cmin)
            DDcmax = sum(D <= self.cmax)
        else:
            DDcmin = sum(D >= self.cmax)
            DDcmax = sum(D <= self.cmin)

        if DDcmin + DDcmax > 15:
            return(True)
        elif (DDcmin < 1) | (DDcmax < 1):  # si toute points sont inf a cmin ou sup a cmax
            return('out')
        else:
            return(False)

    def valid_v(self, v):
        """
        valid_v(v) : check if vertex are valid for the given constraint

        A box is valid if it not not valid

        A box is not valid if all distances are greater than rangemax
                           or all distances are less than rangemin
        """
        ppb = pow(2, len(self.p[0, :]))
        nbbox = len(v) / ppb
        DDbound = np.zeros((4, len(v)), dtype='bool')
        TB = np.zeros((4, nbbox), dtype='bool')

        P0 = np.outer(np.ones(len(v)), self.p[0, :])
        P1 = np.outer(np.ones(len(v)), self.p[1, :])
        F1v = np.sqrt(np.sum((P0 - v) * (P0 - v), axis=1))
        F2v = np.sqrt(np.sum((P1 - v) * (P1 - v), axis=1))
        DD = (F1v - F2v)

#               if self.value > 0:
#                       DD2 = (D>=self.cmin) & (D<=self.cmax)
#               else :
#                       DD2 = (D>=self.cmax) & (D<=self.cmin)

        if self.value > 0:
            # calculate all distance from constraint origin to all vertexes
            #DD = np.sqrt(np.sum(D*D,axis=1))
            # for each box , find the vertex closest to the constraint origin and the farest.
            T = np.array((np.min(DD.reshape(nbbox, ppb),
                                 axis=1), np.max(DD.reshape(nbbox, ppb), axis=1)))

            TB[0, :] = (T[0, :] <= self.cmin)
            TB[1, :] = (T[0, :] <= self.cmax)
            TB[2, :] = (T[1, :] <= self.cmin)
            TB[3, :] = (T[1, :] <= self.cmax)
            DDbound[0, :] = (DD >= self.cmin)
            DDbound[1, :] = (DD <= self.cmax)

        else:
            # calculate all distance from constraint origin to all vertexes
            #DD = np.sqrt(np.sum(D*D,axis=1))
            # for each box , find the vertex closest to the constraint origin and the farest.
            T = np.array((np.min(DD.reshape(nbbox, ppb),
                                 axis=1), np.max(DD.reshape(nbbox, ppb), axis=1)))

            TB[0, :] = (T[0, :] <= self.cmax)
            TB[1, :] = (T[0, :] <= self.cmin)
            TB[2, :] = (T[1, :] <= self.cmax)
            TB[3, :] = (T[1, :] <= self.cmin)
            DDbound[0, :] = (DD >= self.cmax)
            DDbound[1, :] = (DD <= self.cmin)

        return DDbound, TB

#
#               return DD2
    def inclusive(self, b):
        """
        inclusive(b) : A box b is inclusive for the constraint if self.p is included in the box

        """
        if b.inbox(self.p[0]) | b.inbox(self.p[1]):
            return True
        else:
            return False
