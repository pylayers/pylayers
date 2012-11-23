#-*- coding:Utf-8 -*-
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
from pylayers.util.project import *
import numpy as np
import scipy as sp
import time
from pylayers.location.geometric.util.boxn import *
from pylayers.location.geometric.util import geomview as g
from pylayers.location.geometric.util.scene import *
import os
import sys
import numpydoc

#__docformat__ = 'reStructuredText'


class CLA(object):
    """
    Constraint Layer Array class
    The Constraint Layer Array gather all constraints and process them.

    :Parameters:
    ============
            c  : list
                    contraints contained in CLA
            type : list
                    types of contraints contained in CLA
            std  : list
                    standard deviation of constraints
            vcw   : list
                    scale factor of constraints
            Nc : integer
                    Layer number of current processing
            pe :  np.array
                    Position estimated
            dlayer  : dictionnary
                    key : Layer number
                    value : list of Enclosed (0) and ambiguous (1) boxes.
            iter : integer
                    current iteration of refine process
            erronous : list
                    fills with number constraint which are not compatibleselfselfselfself.

    :Methods:
    =========
            info(self)                              : Give info

            rescale(self,f_vcw,cid=None)            : rescale Constraint Box

            annulus_bound(self,cid=None)            : rescale Constraint

            append(self,c)                          : Append a Constraint to CLA

            setvcw(self,vcw):                       : Set vcw for all constraint

            merge2(self,vcw_init=1.0)               : Merge all constraint from the CLA

            valid_v(self,lv,N)                      : Test vertexes with all constraints

            refine(self,l,NBOXMAX=100,VOLMIN=0.1)   : reduce the validity zone

            show3(self,l=-1,amb=False,sc='all')     : show3

            prob(self,c,d)                          : Compute DDP for the given vertexes

            gapdetect(self,l,dlindx)                : Gap detection for bimodal solution

            min_dist(self,a,b)                      : OBSOLETE

            estpos2(self,l=-1,amb=False)            : Position estimation
    """
#       MEMBERS
#               Nc        : number of constraints
#               c         : list of constraints                                  1 x Nc
#               std       : list of standard deviation of constraints            1 x Nc
#                          if std = 0 it means the constraint is hard and it force the
#                          the estimated point to belong to the bounding box of this
#                          constraint
#               w         : list of weight of constraints                        1 x Nc
#                          if w = 0 it means the constraint is hard and it force the
#                          the estimated point to belong to the bounding box of this
#                          constraint
#
#               validity  : validity array (N x Nc)
#               dlayer    : dictionnary containing a list of 2 elements :
#                               - the list of boxes that are inside the validity area (VA)
#                               - the list of boxes which at least an edge is inside the validity area(VA)
#               dpe       : dictionnary containing the estimated points
#       :Methods:
#               info()
#               append(c,std)
#               remove(c,k)
#               merge2()
#               layer(lbox,l=-1)
#               grid(l=-1,Msz=1000)
#               eval(Msz=1000)
#               show3()
#               estpos(amb=False)

#       List of elementary Constraints

    def __init__(self, parmsh={}):
        self.c = []
        self.type = []
        self.std = []
        self.w = []
        self.vcw = []
        self.Nc = 0
        self.pe = np.array([])
        self.dlayer = {}
        self.iter = 0
        self.erronous = []
        self.id = []
        self.origin = []
        self.runable = []

        if len(parmsh) == 0:
            self.parmsh = parmsh
            self.parmsh['display'] = False     # launch geomview interactively
            self.parmsh['scene'] = False      # display whole scene
            self.parmsh['boxes'] = True       # display constraint box
            self.parmsh['constr_boxes'] = True       # display constraint box
            self.parmsh['estimated'] = True  # display estimated point

        else:
            self.parmsh = parmsh

    def info(self):
        """gives info on the CLA

        .. todo::

                improve info

        """
        N = len(self.c)
        print "CLA : ", N, " constraints"
        print "-----------------------"
        for k in range(N):
            print "Constraint N° ", k
            self.c[k].info()

    def update(self):
        """update
                update all constraints of the CLA
        """
        [c.update() for c in self.c if c.runable]


    def compute(self,pe=True):
        """
        Compute the cla to estimate the postion
    
        Parameters
        ----------
            pe : boolean 
               set to True to compute the position estimation store into self.pe

        Returns
        -------
            boolean
                True if the position estimation has been performed.

        """
        self.merge2()
        self.refine(self.Nc)
        if (sum(self.runable) >= 3) and (pe == True):
            self.estpos2()
            return True
        else:
            return False
        self.Nc=len(self.c)


    def rescale(self, f_vcw, cid=None):
        """idem setvcw but update current vcw with a multiplier factor

        change vcw for all constraints of the CLA

        ...

        .. todo::

                Thing to PROPERLY merge with self.setvcw

        Parameters
        ----------
                f_vcw : a scale factor of the current vcw of the constraint.
                cid : a list of constraints for which the self.vcw will be applied. If cid=None, all constraints are updates. default=None

        Returns
        -------
                Nothing but update vcw either for each constraints from cid list either for all contraints in the CLA list self.c.
        """
        #print "rescale",vcw

        if cid is None:
            [c.rescale(f_vcw * c.vcw) for c in self.c]
        else:
            [c.rescale(f_vcw * c.vcw) for c in self.c if c.Id in cid]

    def annulus_bound(self, cid=None):
        """ adapt cmin and cmax of constraints

        Update cmin and cmax of constraints for a given self.vcw

        ...


        :Parameters:
                cid : a list of constraints for which the self.vcw will be applied. If cid=None, all constraints are updates. default=None

        :Returns:
                Nothing but update boxe size either for each constraints from cid list either for all contraints in the CLA list self.c.
        """
        #print "rescale",vcw
        if cid is None:
            [c.annulus_bound() for c in self.c]
        else:
            [c.annulus_bound() for c in self.c if c.Id in cid]

    def append(self, c):
        """add a constraint into the CLA

        add a constraint into the CLA

        ...


        :Parameters:
                c       : any constraint wichi heritates from Constraint object

        :Returns:
                Nothing but fills self.c list of constraints

        """
        self.c.append(c)
        self.id.append(c.id)
        self.origin.append(c.origin)
        self.type.append(c.type)
        self.runable.append(c.runable)
        self.std.append(c.std)
        self.Nc = self.Nc + 1
        self.vcw.append(c.vcw)
        #
        # Reevaluate weights
        #
        u = np.nonzero(np.array(self.std) > 0)   # std >0
        sumstd = np.sum(np.array(self.std)[u], axis=0).astype('float')
        self.w = np.array(self.std) / sumstd
        self.ndim = c.lbox.ndim

    def remove(self, k):
        """OBSOLETE/ TO BE DEVELOPPED

        ...

        remove(k) : remove a constraint to cla
        """
        self.c.remove(self.c[k])
        self.std.remove(c.std[k])
        sumstd = np.sum(np.array(self.std)[u], axis=0).astype('float')
        self.Nc = self.Nc - 1
        #
        # Reevaluate weights
        #
        u = np.nonzero(np.array(self.std) > 0)   # std >0
        sumstd = np.sum(np.array(self.std)[u], axis=0)
        self.w = np.array(self.std) / sumstd

    def setvcw(self, vcw):
        """update scale factor of all constraint

        rescale all the constraints's boxes according to the given vcw

        ...


        .. todo::

                rename 'setvcw' or trying to merge with 'self.rescale' in order to be homogene with all constraints (TOA,TDOA,RSS)

        :Parameters:
                vcw     : a vcw value

        :Returns:
                Nothing but update all constraint from the CLA
        """
        for c in self.c:
            c.rescale(vcw)

    def merge2(self, vcw_init=1.0):
        """Merge all constraints from the CLA2_reduc2

        Inteligent merging of  constraints in the CLA and look for the smallest intersection box of all the constraints through a dichotomous process.

        - if the result of this merging is empty (no common intersections between all the boxes), all the constraints's vcw are increased (x2) and this processing is operated until an intersection exists (physically intersection MUST exist)
        - if the result of this merging is not empty (intersection exists between all the boxes), all the constraints's vcw are decreased and this processing is operated until no intersection exists. the previous value of vcw is thus used for all constraints.

        This method ensure to find the smallest instersection box satisfaying all the constraints


        Also here is initialized self.dlayer.
        use of dlayer dictionnary:

        self.dlayer[Layer][type of boxes]

        Layer = number of intersecting constraints
        type of boxes : 0 = enclose boxes (EB)
                        1 = ambiguous boxes (AB)

        After the merging, all constraints boxes are store as AB list. EB list is void.


        ...


        :Parameters:
                vcw_init        : intial value of scale factor vcw. This value is updated during the process and affect all constraints ! default =1.0

        :Returns:
                Nothing but fills self.dlayer[Nc][0] (with a void list)  and self.dlayer[Nc][1] (with the intial restricted box). Nc is the number of intersecting constraints
        """

        Nc = self.Nc - len(np.nonzero(np.array(self.type) == 'RSS')[0]) - len(np.nonzero(np.array(self.runable) == False)[0]) 
        self.Nc = Nc
        vcwmin = 1.0  # max(self.vcw)
        step = 1.0
        vcw1 = vcwmin + step

        onlyRSS = False
        if 'RSS' in self.type:
            if 'TOA' not in self.type:
                if 'TDOA' not in self.type:
                    onlyRSS = True
            elif 'TDOA' not in self.type:
                if 'TOA' not in self.type:
                    onlyRSS = True

        while (step > 0.05) | (vcw1 == vcwmin):

            self.setvcw(vcw1)
                #constraints vcw set to current value

            try:
                del tlb
            except:
                pass

            for c in self.c:                # find intersection between all constraints for the current vcw
                if (c.type != 'Exclude'):
                    if (c.type != 'RSS') or onlyRSS:
                        if c.runable:
                            lb = c.lbox
                            try:
                                tlb = tlb.intersect(lb)
                            except:
                                tlb = lb
                    else:
                        pass
                else:
                    ex = c

            if len(tlb.box) == 0:             # if the list is empty (no intersection ) vcw1 is increased
                vcw1 = vcw1 + step
                step = step * 1.2
#                print step, vcw1
            else:                           # if the list is not empty (intersection exist) vcw1 is decreased
                vcw1 = max(vcw1 - step / 2., vcwmin)  # vcw > vcwmin
                step = step / 4.
#                print step, vcw1

        if (np.diff(tlb.box[0].bd, axis=0)[0][0] == 0) | (np.diff(tlb.box[0].bd, axis=0)[0][1] == 0):
            self.setvcw(vcw1 + 1.0)

        try:
            tlb = tlb.intersect(ex.lbox)

        except:
            pass
        self.vcw_init = vcw_init
        self.dlayer[Nc] = [LBoxN([]), tlb]
        self.dlayer[Nc][1].volume()

    def valid_v(self, lv, N):
        """test a vertex list with constraints

        Each vertexes from boxes pass into the list are tested to determine if the box is out (OB), ambiguous (AB) or enclosed (EB)

        ...


        :Parameters:
                self
                lv : a vertex list from BOXN.octants
                N  : number of constraints aka layer number

        :Returns:
                AB : a list with the numerous of Ambiguous Boxes
                EB : a list with the numerous of Enclosed Boxes
        """
        assert N <= self.Nc, " N > Number of Constraints "

        Nmiss = self.Nc - N
        miss_cpt = 0
        f_flag = 0
        o_flag = 0
        pndim = pow(2, self.ndim)
        sDDB = np.ones((4, len(lv)), dtype='bool')
        sT = np.ones((4, len(lv) / pndim), dtype='bool')
        sTAB = np.ones((len(lv) / pndim), dtype='bool')

        TT = []
        Ds = []

        for c in self.c:                # for each constraints

            if (c.type != 'RSS') & (c.type != 'Exclude') & (c.runable):

                DDB, TB = c.valid_v(
                    lv)  # .reshape(2,len(lv)/4,pow(2,self.ndim))
                TT.append(TB)

                if not (DDB[0].any()) | (DDB[1].any()):         # if all  boxes  are out
                    self.erro[c.Id] = self.erro[c.Id] + 1

                sDDB = DDB * sDDB
                # ERROR CHECKER
                AA = TB[0, :]
                BB = TB[1, :]
                CC = TB[2, :]
                DD = TB[3, :]
                TAB = (((~AA) * (~BB) * (DD)) + (BB * (~CC)
                                                 * (~DD)) + (AA * BB * (~CC)))

                sTAB = (sTAB * TAB)

        if self.ndim == 3:
            B = (sDDB[0] * sDDB[1]).reshape(len(lv) / 8, 8)
            sB = np.sum(B, axis=1)
            EB = np.nonzero((sB) > 7)[0]
            AB = np.nonzero((sB > 0) & (sB < 8))[0]
            # error checker
            ABt = np.nonzero(sTAB)[0]
            AB = np.unique(np.hstack((AB, ABt)))
            return (EB, AB)
        if self.ndim == 2:

            B = (sDDB[0] * sDDB[1]).reshape(len(lv) / 4, 4)
            sB = np.sum(B, axis=1)
            EB = np.nonzero((sB) > 3)[0]
            AB = np.nonzero((sB > 0) & (sB < 4))[0]
            # error checker
#                       ABt = np.nonzero(sTAB)[0]
#                       AB  = np.unique(np.hstack((AB,ABt)))
            return (EB, AB)

    def refine(self, l, NBOXMAX=32, VOLMIN=0.001):

        """refine the l layer of the CLA

        Refine the l layer of the CLA  until the maximum number of boxes (NBOXMAX) or the minimal volume of boxes (VOLMIN) has been reached.

        Once the CLA has been merged, this method aims to enclose the solution thanks to an octree/quadtreee process

        self.dlayer[l][0] : LBox which contains boxes inside of the validity area (VA)
        self.dlayer[l][1] : LBox which contains ambiguous boxes (partially inside of the VA == at least 1 edge inside the VA)
        All boxes partially inside of the VA are divided into octants. Each octants are tested into the self.valid.


        ...

        :Parameters:
                self
                l : the layer number
                NBOXMAX : the maximum number of obtained boxes
                VOLMIN :  the minimum volume achievable by the obtained boxes


        :Returns:
                Nothing, but fills self.dlayer[l][0] and self.dlayer[l][1] respectively with enclosed boxes and ambiguous boxes
        """

        self.iter = self.iter + 1

        Nc = self.Nc
        if self.iter == 1:
            #print NBOXMAX
            self.FINISHED = 0

        self.erro = np.zeros(self.Nc)

        a = []
#        print 'iter', self.iter
        B = self.dlayer[l][1].octant()

        lv = B.bd2coord()

        EB, AB = self.valid_v(lv, l)
        del lv
        self.erronous.append(self.erro)

        nbox = len(EB)
        nboxamb = len(AB)

#        print nbox
#        print nboxamb
        # if all boxes are out of the VA...
#               if  ((nboxamb==0)&(nbox==0)) and len(self.dlayer[l][0].box) == 0:
        if  ((nboxamb == 0) & (nbox == 0)) and len(self.dlayer[l][0].box) == 0:
            if self.iter < 25:

                pb = np.nonzero(self.erro != 0)[0]
                if len(pb) != 0:
#                    print "specific size up", pb
                    self.rescale(1.2, pb)
                    self.annulus_bound(pb)
                else:

#                    print 'all contraints size up '
                    self.rescale(1.2)
                    self.annulus_bound()

                self.refine(l)

            else:

                self.iter = 0
                self.dlayer[l - 1] = self.dlayer[l]
                    # unstack to a lower the layer
                l = l - 1

                assert l >= 0, pdb.set_trace()

                self.refine(l)

        # if it exists at least a box ambiguous or not in the VA
        else:

            if (nbox != 0 and nboxamb == 0):
                self.FINISHED = 1

            # Update EB
            if len(EB) != 0:
                self.dlayer[l][0].append_l(LBoxN(B.box[EB]))

            # Update AB
            self.dlayer[l][1] = LBoxN(B.box[AB], ndim=self.ndim)

            # check if it remains is more AB to refine
            if nboxamb != 0:
                lv = 1
            else:
                lv = 0

            # while the max number of boxes (NBOXMAX) is not reached or the elementary volume of boxes (VOLMIN) is not reached
            # self.refine is executed.
            # else  self.refine is over.

            if (((nboxamb + nbox) < NBOXMAX) and (self.dlayer[l][lv].box[-1].vol > VOLMIN)) and self.FINISHED == 0:
                self.refine(l)
            else:
                self.iter = 0
                self.Nc = l

    def show3(self, l=-1, amb=False, sc='all'):
        """ Display constraints and theirs boxes through geomview.


        geomview parameters are the following

        self.parmsh['display']=False            # launch geomview interactively
        self.parmsh['scene']=True               # display whole scene
        self.parmsh['boxes']=True               # display constraint box
        self.parmsh['constr_boxes']=False       # display constraint box
        self.parmsh['estimated']=True           # display estimated point

        ...


        .. todo:

                create a .ini file for geomview parameters

        :Parameters:
                l       : layer number to observe. If -1 estimation is made on the highest available layer. default = -1
                amb     : display ambiguous boxes. default = false
                sc      : display all constraint or give a list with the constrinat number to observe ex: [0,1,3]. default 'all'

        :Returns:
                Nothing but calls a geomview instance


        """
        Nc = self.Nc
        filename = basename + "/geom/cla.list"
        fd = open(filename, "w")
        fd.write("LIST\n")
        par = self.parmsh
        if par['constr_boxes']:

            if l == -1:
                if sc == 'all':
                    for c in self.c:
                        if c.runable:
                            c.parmsh['display'] = False
                            c.parmsh['scene'] = False
                            fname = c.show3()
                            fd.write("{<" + fname + ".list}\n")

                else:
                    try:
                        for vsc in sc:
                            if vsc.runable:
                                self.c[vsc].parmsh['display'] = False
                                self.c[vsc].parmsh['scene'] = False
                                fname = self.c[vsc].show3()
                                fd.write("{<" + fname + ".list}\n")
                    except:
                        if sc.runable:
                            self.c[sc].parmsh['display'] = False
                            self.c[sc].parmsh['scene'] = False
                            fname = self.c[sc].show3()
                            fd.write("{<" + fname + ".list}\n")

            else:
                if c[l].runable:
                    self.c[l].parmsh['dispay'] = False
                    self.c[l].parmsh['scene'] = False
                    fname = self.c[l].show3()
                    fd.write("{<" + fname + ".list}\n")

        col = ['r', 'b', 'g', 'm', 'y', 'b', 'r']

        if par['scene']:
            an = np.zeros(len(self.bn))
            for c in self.c:
                if c.runable:
                    an = np.vstack((an, c.p))

            S = Scene(an=an, bn=self.bn)
            sce = S.generate()

        if par['estimated']:
            try:
                sce = g.cloud(self.pe, display=False, name='scene',
                              color='k', dice=6, access='append')
                fd.write("{<" + sce + "}\n")
            except:
                pass

        if par['boxes']:

            for l in self.dlayer.keys():
                self.dlayer[l][0].parmsh['display'] = False
                self.dlayer[l][1].parmsh['display'] = False

                try:
                    fname = self.dlayer[l][0].show3(col=col[Nc - l + 1], Id=l)
                    fd.write("{<" + fname + "}\n")
                except:
                    pass

                if amb:
                    coco = ['r', 'v', 'b', 'y']
                    fname = self.dlayer[l][1].show3(col=col[Nc - l], Id=l + 1)
#                                       fname = self.dlayer[l][1].show3(col=coco,Id=l+1)
                    fd.write("{<" + fname + "}\n")

        fd.close()
        chaine = "geomview  -nopanel  -b 1 1 1 " + filename + " 2>/dev/null &"

        os.system(chaine)

    def prob(self, c, d):
        """ determine probability of list of vertex

        Return the probability of each vertex from an array in regard of the constraint origin, standard deviation and vcw

        ...


        :Parameters:
                c       : contraint number in the self.c list
                d       : an array of vertex

        Retun:
                v       : probability of each vertex

        """

        if self.c[c].type == 'TDOA':
            v = (1 / ((self.c[c].sstd * self.c[c].vcw) * np.sqrt(2 * np.pi))) * np.exp(-(d - self.c[c].value * 0.3) ** 2 / (2 * (self.c[c].sstd) * self.c[c].vcw) ** 2)

        elif self.c[c].type == 'TOA':

            v = (1 / (((self.c[c].sstd) * self.c[c].vcw) * np.sqrt(2 * np.pi))) * np.exp(-(d - self.c[c].value * 0.3) ** 2 / (2 * (self.c[c].sstd) * self.c[c].vcw) ** 2)

        elif self.c[c].type == 'RSS':
#
#                       v = (1/(((self.c[c].sstd)*self.c[c].vcw)*np.sqrt(2*np.pi)))*np.exp(-(d-self.c[c].value*0.3)**2/(2*(self.c[c].sstd)*self.c[c].vcw)**2)
#                       v=v[0]
            S = (-self.c[c].sstd * np.log(10)) / (-10 * self.c[c].model.RSSnp)
            M = ((self.c[c].model.PL0 - self.c[c].value) *
                 np.log(10)) / (10 * self.c[c].model.RSSnp)
            v = 1 / (d * S * np.sqrt(2 * np.pi)) * np.exp(
                -(((np.log(d) - M) ** 2) / (2. * (S ** 2))))

#                       std = self.c[c].sstd#10**(self.c[c].model['RSSnp']/20.)
##                      mean = self.c[c].range
#                       mean = np.log(self.c[c].range)+std**2
#                       v = 1/(d*np.sqrt(2*np.pi))*np.exp(-(np.log(d)-mean)**2/(2*std**2))
        return(v)

    def gapdetect(self, l, dlindx):
        """basic gap detection

        Detects if separated clusters of boxes are observables. his situation is usual in under determined estimation.
        This only test on each axis if all boxes are contiguous. If not, a gap is declared and clusters are created.


        ...


        :Parameters:
                l       : layer number
                dlindx  : select the boxes type ( from self.dlayer) for gap detection 0=enclose or 1=ambigous boxes

        :Returns:
                clust   : a list of array. each array contains boxes from the same cluster
                axis    : axis/axes where gap has/have been detectes

        """

        gcoord = []
        axis = np.zeros(self.ndim, dtype='int8')
        clust = []
        for i in range(self.ndim):
            uni, inv, idd = np.unique(self.dlayer[l][dlindx]
                                      .bd[:, i], return_inverse=True, return_index=True)
#                       uni,inv,idd =np.unique(self.dlayer[l][dlindx].ctr[:,i],return_inverse=True,return_index=True)

            slope = np.diff(np.diff(uni))
            if len(slope) != 0:
                if abs(np.min(slope)) > 1e-9:

                    gidx = np.nonzero(np.min(slope) == slope)[0]
#                                       print 'GAP DETECTED in AXIS',i
                    axis[i] = 1

                    try:
                        # divisé par 2 pour pouvoir aveir les index  de cluster comme les centre des box
                        clust.append(np.nonzero(uni[gidx[0]] < self.dlayer[l]
                                                [dlindx].bd[:, i])[0] / 2)
                        clust.append(np.nonzero(uni[gidx[0]] > self.dlayer[l]
                                                [dlindx].bd[:, i])[0] / 2)

                    except:
                        pdb.set_trace()
            else:
                clust = []

        return clust, axis

    def min_dist(self, a, b):
        """OBSOLETE
        """
        print 'min dist'
        pdb.set_trace()
        # recherche distance entre barycentre et les centres des boites distance2barycentre(d2b)
        d2b = np.sqrt(np.sum((a - b) ** 2, axis=1))
        # on retourne pe comme etant le centre de la boite ayant le plus faible distrance avec barycentre
        indx = np.nonzero(d2b == min(d2b))[0]
        return(indx[0])

    def estpos(self, l=-1, amb=False, test=False):
        """
        estpos(l,amb=True) : estimate position

        l : layer number
        amb : if True include ambigous boxes of VA in gravity center computation
        """

        if l == -1:
            l = np.max(self.dlayer.keys())

        PP = []
        self.saveP = []

        for dlindx in range(2):
            for i in range(len(self.dlayer[l][dlindx].box)):
                poids = []
                for j in range(len(self.c)):
                    if self.c[j].type != 'Exclude':

                        d = np.sqrt(np.sum((self.dlayer[l][dlindx].box[i].ctr - self.c[j].p) * (self.dlayer[l][dlindx].box[i].ctr - self.c[j].p)))

                        poids.append(self.prob(j, d))

                P = sum(np.array(poids) * np.array(poids)) / (len(poids))
                self.saveP.append(P)

                PP.append(P * self.dlayer[l][dlindx].box[i].ctr)

        self.pe = np.sum(PP, axis=0) / np.sum(self.saveP)

    def estpos2(self, l=-1, amb=False):
        """ Position estimation

        estimate position from the enclosed or/and ambibuous boxes

        ...


        :Parameters:

                l       : Layer of the estimation. If -1 estimation is made on the highest available layer
                amb     : Use ambiguous boxes (if available) to perform the position estimation. default = False

        :Returns :
                Nothing but fills self.pe with an array

        """

        if l == -1:
            l = np.max(self.dlayer.keys())

        PP = []
        poids = []

        if len(self.dlayer[l][0].box) != 0:  # si enclosed box exists
            dlindx = 0
#            print 'Enclosed pos estim'
        else:
            dlindx = 1
#            print 'Amiguous pos estim'
        self.saveP = np.zeros((len(self.dlayer[l][dlindx].box)))
        clust, axis = self.gapdetect(l, dlindx)

        box_center = self.dlayer[l][dlindx].ctr

        for j in range(len(self.c)):
            #if self.c[j].type != 'Exclude':
            if (self.c[j].type != 'Exclude') & (self.c[j].runable):

                # compute distance between contraint center and all vertexes
                if self.c[j].type == 'TOA' or self.c[j].type == 'RSS':
                    d = np.sqrt(np.sum((box_center - self.c[j].p * np.ones((len(box_center), 1))) ** 2, axis=1))
                elif self.c[j].type == 'TDOA':
                    F1v = np.sqrt(np.sum((self.c[j].p[0] - box_center) * (self.c[j].p[0] - box_center), axis=1))
                    F2v = np.sqrt(np.sum((self.c[j].p[1] - box_center) * (self.c[j].p[1] - box_center), axis=1))
                    d = (F1v - F2v)

                try:
                    poids = (poids * (self.prob(j, d)))
                    poids = (poids * poids.T) / len(poids)

                except:
                    poids = (self.prob(j, d))
                    poids = (poids * poids.T) / len(poids)
#                                       poids.append(self.prob(j,d))

#                       pdb.set_trace()
#                       P=sum(np.array(poids)*np.array(poids))/(len(poids))
#                       self.saveP[i]=P
        self.saveP = poids
#                       PP.append(P*self.dlayer[l][dlindx].box[i].ctr)
##########################################

        if clust != []:
            self.pecluster=[]
            print 'cluster'
            lclust = []
            dd = []
            mps = -1.0
            saxis = sum(axis)

            p = 1

            for i in range(len(axis)):
                if axis[i] != 0:
                    try:
                        count = np.vstack((count, np.repeat(range(2 * (p - 1), (2 * (p - 1)) + 2) * (pow(2, saxis - p)), p)))
                    except:
                        count = np.repeat(range(2 * (p - 1), (2 * (p - 1)) + 2)
                                          * (pow(2, saxis - p)), p)
                    p = p + 1
            count = count.T

            for i in range(len(clust)):
                if len(clust) < 3:
                    clusters = clust[i]
                else:

                    if len(np.shape(count)) > 1:
                        clusters = np.intersect1d(clust[count[i,
                                                              0]], clust[count[i, 1]])
                    else:
                        clusters = np.intersect1d(clust[count[
                            0]], clust[count[1]])

                clust_vol = np.sum(np.array(self.dlayer[l][
                    dlindx].vol)[np.unique(clusters)])

                if len(clusters) != 0:
                    mp = np.max(self.saveP[clusters])

                    if mps < mp:
                        mps = mp
                        estclu = clusters

                if clust_vol != 0:
                    lclust.append(clusters)
                    pc = np.sum(np.array(self.dlayer[l][dlindx].ctr)[np.unique(clusters)], axis=0) / len(np.unique(clusters))
                    dd.append(np.sqrt(np.sum((pc - self.c[0].p) ** 2)))

#                       try:
#                               vmax=[]
#                               for i in range(len(lclust)):
#                                       vmax.append(np.max(poids[np.unique(lclust[i])]))
#                               peindx = np.nonzero(poids==max(vmax))[0][0]
#                               self.pe = self.dlayer[l][dlindx].ctr[peindx]

            try:
                M = (((-self.c[0].model['PL0'] - self.c[0].value) * np.log(10)
                      ) / (10. * self.c[0].model['RSSnp']))[0]
                LL = np.log(dd[1] / dd[0]) * (1 + np.log(
                    dd[0] * dd[1]) - 2 * M)

                if LL > 0:
#                                       vmax = np.max(poids[np.unique(lclust[0])])
#                                       peindx=np.nonzero(poids[vmax]==poids)[0][0]
#                                       self.pe = self.dlayer[l][dlindx].ctr[np.unique(lclust[0])[peindx]]

                    self.pe = np.mean(self.dlayer[l][dlindx].ctr[
                        np.unique(lclust[0])], axis=0)
                    pestdmax = np.max(self.dlayer[l][
                        dlindx].ctr[np.unique(lclust[0])])
                    pestdmin = np.min(self.dlayer[l][
                        dlindx].ctr[np.unique(lclust[0])])
                    self.pestd = pestdmax - pestdmin
                else:
#                                       vmax = np.max(poids[np.unique(lclust[1])])
#                                       peindx=np.nonzero(poids[vmax]==poids)[0][0]
#                                       self.pe = self.dlayer[l][dlindx].ctr[np.unique(lclust[1])[peindx]]
                    

                    self.pe = np.mean(self.dlayer[l][dlindx].ctr[
                        np.unique(lclust[1])], axis=0)
                    pestdmax = np.max(self.dlayer[l][
                        dlindx].ctr[np.unique(lclust[1])])
                    pestdmin = np.min(self.dlayer[l][
                        dlindx].ctr[np.unique(lclust[1])])
                    self.pestd = pestdmax - pestdmin

                    for cl in lclust:
                        self.pecluster.append(np.mean(self.dlayer[l][dlindx].ctr[
                        np.unique(cl)], axis=0))

            except:

                if np.sum(poids) > 0.:
                    self.pe = np.sum(poids * self.dlayer[l][dlindx]
                        .ctr.T, axis=1) / np.sum(poids)
                else:
                    self.pe = np.sum(self.dlayer[l][dlindx].ctr, axis=0) / \
                        len(self.dlayer[l][dlindx].ctr)
                pestdmax = np.max(self.dlayer[l][dlindx].bd, axis=0)
                pestdmin = np.min(self.dlayer[l][dlindx].bd, axis=0)
                self.pestd = pestdmax - pestdmin
                for cl in lclust:	
                    self.pecluster.append(np.mean(self.dlayer[l][dlindx].ctr[
                    np.unique(cl)], axis=0))
                #self.pe = np.sum(self.dlayer[l][dlindx].ctr,axis=0)/len(self.dlayer[l][dlindx].ctr)
#                       try:
#                               print 'clust vol',clust_vol
#                               self.pe = np.mean(self.dlayer[l][dlindx].ctr[estclu],axis=0)
#
#                       except:
#                               pdb.set_trace()
#
#                               PP=self.saveP[clust[estclu]/2]*self.dlayer[l][dlindx].ctr[clust[estclu]/2].T
#                               self.pe = np.sum(PP.T,axis=0)/np.sum(self.saveP[clust[estclu]/2])
#                               pdb.set_trace()
        else:
            if np.sum(poids) > 0.:
                self.pe = np.sum(poids * self.dlayer[l][
                    dlindx].ctr.T, axis=1) / np.sum(poids)
            else:
                self.pe = np.sum(self.dlayer[l][dlindx].ctr,
                                 axis=0) / len(self.dlayer[l][dlindx].ctr)
            pestdmax = np.max(self.dlayer[l][dlindx].bd, axis=0)
            pestdmin = np.min(self.dlayer[l][dlindx].bd, axis=0)
            self.pestd = pestdmax - pestdmin
