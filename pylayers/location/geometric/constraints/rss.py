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
import ConfigParser
from pylayers.util.project import *
import pylayers.util.pyutil as pyu

from pylayers.location.geometric.util.boxn import *
from pylayers.location.geometric.constraints.constraint import *
#from pylayers.location.algebraic.RSSLocation import *
from pylayers.network.model import *


class RSS(Constraint):
    """     RSS Constraint

    Description and evaluation of TOA constraints

    :Parameters:
            value   : float
                    Constraint value in dB. Default = 40
            std     : float
                    Value standard deviation in dB. default = 1.0
            vcw     : float
                    scale factor. Default = 1.0
            model   : dictionnary
                    keys:   self.model['PL0'] =-34.7
                            self.model['d0']  = 1.0
                            self.model['RSSnp'] = 2.64
                            self.model['RSSStd'] = 4.34
                            self.model['Rest'] = 'mode'
            p       : np.array 1 x ndim
                    constraint center

    :Attributes:

            range   : distance conversion from time self.value.
            sstd    : distance conversion from time self.std
            runable : True NOT USED
            evaluated :False NOT USED
            self.id : Constraint ID

            from annulus bound:
            min     : minimum value of observable
            max     : maximum value of observable
            mean    : mean value of observable


    :Methods:
            annulus_bound(self)     : Compute the minimum and maximum distance of the enclosing annulus of the constraint
            rescale(self,vcw)       : rescale contraint boundary with a given scale factor 'vcw'
            inclusive(self,b)       : Is constraint center is inside a box ?
            valid(self,b)           : Test if Lbox is compatible with the constraint
            valid_v(self,lv)        : Test if a liste of a vertexes from a box is compatible with the constraint. vertexes are obtained thanks to LBoxN.bd2coordinates()
            estvol(self)            : Constraint Volume estimation

    """

    def __init__(self, id='0', value=0, std=0, vcw=3, p=np.array([]), model={}, origin={}):
        Constraint.__init__(self, type='RSS', id=id, p=p, origin=origin)
#               Constraint.C_id = Constraint.C_id+1   # constraint counter is incremented
        self.value = value  # attennation (dB)
        self.std = std
        self.vcw = vcw
        if model == {}:
            self.config = ConfigParser.ConfigParser()
            self.config.read(pyu.getlong('EMSolver.ini', 'ini'))
            param = dict(self.config.items('rat1_PLM'))
            self.model = PLSmodel(f=eval(param['f']), 
                                  rssnp=eval(param['rssnp']), 
                                  d0=eval(param['d0']), 
                                  method=param['method'])
#                       self.model={}
#                       self.model['PL0'] =-34.7
#                       self.model['d0']  = 1.0
#                       self.model['RSSnp'] = 2.64
#                       self.model['RSSStd'] = 4.34
#                       self.model['Rest'] = 'mode'
            self.param = self.model.param
        else:
            self.model = model
            self.param = self.model.param
        self.usable = True
        self.visible = True
        self.obsolete = False
        self.update()

    def update(self):
        """
        update constraint inforamtion
        """
        if self.p.any():
            self.runable = True
        else:
            self.runable = False
#               self.LOC = RSSLocation(self.p)
        self.sstd = self.model.getRangeStd(self.value, self.std)  # (self.LOC.getRangeStd(self.p, self.model['PL0'], self.model['d0'], self.value,self.model['RSSnp'], self.model['RSSStd'], self.model['Rest']))/0.3
        self.range = self.model.getRange(self.value, self.std)  # self.LOC.getRange(self.p, self.model['PL0'], self.model['d0'], self.value, self.model['RSSnp'], self.model['RSSStd'], self.model['Rest'])
#               self.sstd   = self.std*0.3
        self.rescale(self.vcw)
        self.evaluated = False
        self.annulus_bound()

    def annulus_bound(self):
        """
        annulus_bound():
        Compute the minimum and maximum distance of the enclosing annulus of the constraint

        :Returns:
                Nothing but update cmin, cmax
        """
        self.cmin = max(0, self.range - self.vcw * self.sstd)
        self.cmax = self.range + self.vcw * self.sstd

    def rescale(self, vcw):
        """rescale bounds of the constraint with vcw factor

        :Parameters:
                vcw : float
                        scale factor

        """
        self.vcw = vcw
        dx = self.range + (self.vcw * self.std) * 0.3 * \
            np.ones(len(self.p))
        box = BoxN(
            np.array([self.p - dx, self.p + dx]), ndim=len(self.p))

        self.lbox = LBoxN([box])
        self.estvol()

    def inclusive(self, b):
        """ Is constraint center is inside a box ?
        :Parameters:
                b       : LBoxN

        :Returns:
                Boolean

        """
        if b.inbox(self.p):
            return True
        else:
            return False



    def valid_v(self, v):
        """ Test if a liste of a vertexes from a box is compatible with the constraint. Vertexes are obtained thanks to LBoxN.bd2coordinates()

        valid_v(v) : check if a set of vertex are valid for the given constraint

        Returns

        - DDbound ( boxes validity)
        - TB    for error checker

        DDbound = list[[tested vertex DD>self.cmin],[tested vertex DD<self.cmax]]


                                    nb box
                                <-------------->
                pmin<cmin       |
                pmin<cmax       |
        TB =    pmax<cmin       |
                pmin<cmax       |



                                            nb vertexes
                                        <------------------->
                        pmin<Dmin       |
                        pmin<Dmax       |
        DDbound =       pmax<Dmin       |
                        pmin<Dmax       |


        Parameters
        ----------
                v       : np.arrays
                        a vertexes arrays

        Returns
        -------

                DDbound : np.array 2 x vertexes containing boolean
                        Test Lboxn boundaries  with contraint
                TB      : np.array 4 v vertexes
                        Multiple test for error checker in Lboxn boundaries  with contraint



        """

#               DDbound = []
        p0 = self.p
        ppb = pow(2, len(p0))
        nbbox = len(v) / ppb

        DDbound = np.zeros((4, len(v)), dtype='bool')
        TB = np.zeros((4, nbbox), dtype='bool')

        P = np.outer(np.ones(len(v)), p0)
        D = P - v

        # calculate all distance from constraint origin to all vertexes
        DD = np.sqrt(np.sum(D * D, axis=1))
        # for each box , find the vertex closest to the constraint origin and the farest.
        T = np.array((np.min(DD.reshape(nbbox, ppb), axis=1),
                      np.max(DD.reshape(nbbox, ppb), axis=1)))

        TB[0, :] = (T[0, :] <= self.cmin)
        TB[1, :] = (T[0, :] <= self.cmax)
        TB[2, :] = (T[1, :] <= self.cmin)
        TB[3, :] = (T[1, :] <= self.cmax)
        DDbound[0, :] = (DD >= self.cmin)
        DDbound[1, :] = (DD <= self.cmax)
#               DDbound[2,:]=(DD<=self.cmin)
#               DDbound[3,:]=(DD<=self.cmax)

        return DDbound, TB









    def valid(self, b):
        """check if a box is valid for the given constraint

        :Parameters:
                b       : LBoxN

        :Return:
                True    : if box is enclosed
                False   : if the box is ambiguous
                'Out'   : if the box is out

        """
        p0 = self.p
        v = self.bd2coord(b)
        P = np.outer(np.ones(len(v)), p0)
        D = P - v
        DD = np.sqrt(np.sum(D * D, axis=1))

        DDcmin = sum(DD >= cmin)
        DDcmax = sum(DD <= cmax)

        if DDcmin + DDcmax > 15:
            return(True)
        elif (DDcmin < 1) | (DDcmax < 1):  # (bmax<cmin)|(bmin>cmax):

            return('out')
        else:
            return(False)

    def valid_vold(self, v):
        """
        .. todo:
                update as TOA valid_v !!

        valid_v(v) : check if vertex are valid for the given constraint

        A box is valid if it not not valid

        A box is not valid if all distances are greater than rangemax
                           or all distances are less than rangemin
        """
        p0 = self.p
        P = np.outer(np.ones(len(v)), p0)
        D = P - v
        DD = np.sqrt(np.sum(D * D, axis=1))
        DD2 = (DD >= self.cmin) & (DD <= self.cmax)

        return DD2







    def estvol(self):
        """ Constraint Volume estimation

        """
        R2 = self.range + self.vcw * self.sstd
        R1 = self.range - self.vcw * self.sstd
        V2 = (4. * np.pi * R2 ** 3) / 3
        V1 = (4. * np.pi * R1 ** 3) / 3
        self.estvlm = V2 - V1
