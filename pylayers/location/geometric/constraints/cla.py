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
try:
    from interval import interval,inf
    pyinterval_installed=True
except:
    pyinterval_installed=False
import os
import sys

#__docformat__ = 'reStructuredText'


class CLA(object):
    """  Constraint Layer Array class
    The Constraint Layer Array gather all constraints and process them.

    Attributes
    ----------

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

    Methods
    -------

    info(self)                              : Give info

    compute(pe=True,mergeRSS=False,refineRSS=True, NBOXMAX=50, VOLMIN=0.001,HT=True,forceamb=False):
                                              compute the CLA to estimate the positon.
    
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
        self.runable = [] # does pe is known ?
        self.visible = [] # does link physically exist ? aka 2 nodes are in visiblity ?
        self.obsolete = [] # is the ldp has been obtain a long time ago
        self.usable=[] # constraints are usable = runable + visible



        if len(parmsh) == 0:
            self.parmsh = parmsh
            self.parmsh['display'] = False     # launch geomview interactively
            self.parmsh['scene'] = False      # display whole scene
            self.parmsh['boxes'] = True       # display constraint box
            self.parmsh['constr_boxes'] = True       # display constraint box
            self.parmsh['estimated'] = True  # display estimated point

        else:
            self.parmsh = parmsh

    def __repr__(self):
        np.set_printoptions(precision=3)
        s = '{0:4} | {1:6} |{2:4} | {3:4} | {4:15}| {5:9}| {6:5}| {7:7}| {8:6}|'.format('node','peer','type', 'rat', 'p', 'value', 'std', 'runable' , 'usable' )
        for c in self.c:
            node = c.origin['id']
            peer = c.origin['link']
            rat  = c.origin['rat']
            if c.type != 'TDOA':
                s = s + '\n' + '{0:4} | {1:6} |{2:4} | {3:4} | {4:15}| {5:9}| {6:5}| {7:7}| {8:6}|'.format(node,peer,c.type,rat, c.p, c.value, c.std, c.runable, c.usable)
            else:
                s = s + '\n' + '{0:4} | {1:6} |{2:4} | {3:4} | {4:15}| {5:9}| {6:5}| {7:7}| {8:6}|'.format(node,peer,c.type,rat, c.p[0], c.value, c.std, c.runable, c.usable)
                s = s + '\n' + '                            '+str(c.p[1])
                # s = s + '\n' + '{0:4} | {1:15}| {2:9}| {3:5}| {4:7}| {5:6}| {6:8}| {7:9}'.format(c.type, c.p[0], c.value, c.std, c.runable, c.usable , c.obsolete , c.evaluated)
                # s = s + '\n' + '       '+str(c.p[1])

        s = s + '\n\n' + 'position evaluated by the CLA\n' + str(self.pe)
        return s

    def info(self):
        for c in self.c:
            c.info()


    def update(self):
        """update
                update all constraints of the CLA
        """
        [c.update() for c in self.c if c.runable]
        self.runable=[c.runable for c in self.c]
        self.obsolete=[c.obsolete for c in self.c]
        self.visible=[c.visible for c in self.c]
        self.usable=[c.usable for c in self.c]

    def compute(self,pe=True,mergeRSS=False,refineRSS=True, NBOXMAX=50, VOLMIN=0.001,HT=True,forceamb=False):
        """
        Compute the cla to estimate the postion
    
        Parameters
        ----------


        pe : boolean 
           set to True to compute the position estimation store into self.pe

        mergeRSS : boolean
            True : if there is RSS in cla, they are used to find the smallest merge
            False (default): even if there is RSS in cla, they are neglected during the merge process

        refineRSS :boolean
            True (default): if there is RSS in cla, they are used to decide if boxes are enclosed of ambiguous
                            during the refine process
            False: if there is RSS in cla, they are ignore during the refine process 

        NBOXMAX : integer 
            Choose the maximum boxes generated during the refine process (escape value of the while and recursive function)

        NVOLMIN : float 
            Choose the minimum volume of the boxes obtained  during the refine process (escape value of the while and recursive function)

        HT : boolean
            True : if a cluster ppears (2 sets of distinct boxes ) an hypthesis testuing method is applied
                    in estpos2 method 
            False : no HT methos is applied 

            Description of the hypothesis testing (HT) method in:
            Hybrid positioning based on hypothesis thesting
            N. Amiot, T. Pedersen, M. Laaraiedh, B. Uguen. 
            A Hybrid Positioning Method Based on Hypothesis Testing
            ,Wireless Communications Letters, IEEE, vol.1, no.4, pp.348-351, August 2012
            http://ieeexplore.ieee.org.passerelle.univ-rennes1.fr/stamp/stamp.jsp?tp=&arnumber=6205594


        Returns
        -------

        return : boolean
            True if the position estimation has been performed.

        update a self.pe which contain the estimated position   


        """
        self.merge2(RSS=mergeRSS)
        self.refine(l=self.Nc,NBOXMAX=NBOXMAX, VOLMIN=VOLMIN,RSS=refineRSS)
        self.update()
        if (sum(self.usable) >= 3) and (pe == True):
            self.estpos2(HT=HT)
            self.Nc=len(np.where(self.usable)[0])
            return True
        elif forceamb:
            self.estpos2(HT=HT)
            return False

        else:
            self.Nc=len(np.where(self.usable)[0])
            return False



#    def compute_amb(self,pe=True,HT=True):

#        self.merge2(RSS=False)
#        self.refine(self.Nc,RSS=False)
#        self.estpos2(HT=HT)
#        self.Nc=len(np.where(self.usable)[0])
#        return True

    def rescale(self, f_vcw, cid=None):
        """idem setvcw but update current vcw with a multiplier factor

        change vcw for all constraints of the CLA


        Parameters
        ----------

        f_vcw : a scale factor of the current vcw of the constraint.
        cid : a list of constraints for which the self.vcw will be applied. If cid=None, all constraints are updates. default=None

        Returns
        -------
        
        Nothing but update vcw either for each constraints from cid list either for all contraints in the CLA list self.c

        """
        #print "rescale",vcw

        if cid is None:
            [c.rescale(f_vcw * c.vcw) for c in self.c]
        else:
            [c.rescale(f_vcw * c.vcw) for c in self.c if c.Id in cid]

    def annulus_bound(self, cid=None):
        """ adapt cmin and cmax of constraints

        Update cmin and cmax of constraints for a given self.vcw

        


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

        


        Parameters
        ----------
                c       : any constraint wichi heritates from Constraint object

        Returns
        -------
                Nothing but fills self.c list of constraints

        """
        self.c.append(c)
        self.id.append(c.id)
        self.origin.append(c.origin)
        self.type.append(c.type)
        self.runable.append(c.runable)
        self.visible.append(c.runable)
        self.obsolete.append(c.obsolete)
        # by default, if a constraint is runable, it will be used
        self.usable.append(c.runable and c.visible and not c.obsolete)
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


        Parameters
        -----------

        vcw     : a vcw value
        RSS : boolean
            True : RSS are considered in merging
            False : RSS are excluded from merging
        Returns
        -------

        Nothing but update all constraint from the CLA

        """
        for c in self.c:
            c.rescale(vcw)


    def merge2(self, vcw_init=1.0, RSS=False):
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


        Parameters
        ----------

        vcw_init : float
            intial value of scale factor vcw. This value is updated during the process and affect all constraints ! default =1.0

        Returns
        -------
        Nothing but fills self.dlayer[Nc][0] (with a void list)  and self.dlayer[Nc][1] (with the intial restricted box). Nc is the number of intersecting constraints
        """

#        Nc = self.Nc - len(np.nonzero(np.array(self.type) == 'RSS')[0]) - len(np.nonzero(np.array(self.runable) == False)[0]) 
#        Nc = self.Nc - len(np.nonzero(np.array(self.runable) == False)[0]) 

        Nc = len(np.where(self.usable)[0])#self.Nc - len(np.nonzero(np.array(self.usable) == False)[0]) 
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
                    if (c.type != 'RSS') or onlyRSS or RSS:
                        if c.usable:
                            lb = c.lbox
                            try:
                                tlb = tlb.intersect(lb)
                            except:
                                tlb = lb
                        else:
                            pass
                else:
                    ex = c
            try:
                tlb = tlb.intersect(ex.lbox)
            except:
                pass

            if len(tlb.box) == 0:             # if the list is empty (no intersection ) vcw1 is increased
                vcw1 = vcw1 + step
                step = step * 1.2
                #print step, vcw1
            else:                           # if the list is not empty (intersection exist) vcw1 is decreased
                vcw1 = max(vcw1 - step / 2., vcwmin)  # vcw > vcwmin
                step = step / 4.
                #print step, vcw1
        try:
            if (np.diff(tlb.box[0].bd, axis=0)[0][0] == 0) | (np.diff(tlb.box[0].bd, axis=0)[0][1] == 0):
                self.setvcw(vcw1 + 1.0)
        except:
            pass
#        try:
#            tlb = tlb.intersect(ex.lbox)

#        except:
#            pass
#        pdb.set_trace()
        self.vcw_init = vcw_init
        self.dlayer[Nc] = [LBoxN([]), tlb]
        self.dlayer[Nc][1].volume()

    def valid_v(self, lv, N, RSS=True):
        """test a vertex list with constraints

        Each vertexes from boxes pass into the list are tested to determine if the box is out (OB), ambiguous (AB) or enclosed (EB)


        Parameters
        ----------

        lv : a vertex list from BOXN.octants
        N  : number of constraints aka layer number
        RSS : boolean
            True : RSS constraints are kept as any other constraints for boxes evaluation (ambigous /enclosed)
            False : RSS constraints are ignored in boxes evaluation (ambigous /enclosed)

        Returns
        -------

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

        if RSS:
            loop_condition="(c.type != 'Exclude') & (c.usable)"
        else :
            loop_condition="(c.type != 'RSS') & (c.type != 'Exclude') & (c.usable)"
        for c in self.c:                # for each constraints
            if eval(loop_condition):

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
            ABt = np.nonzero(sTAB)[0]
            AB  = np.unique(np.hstack((AB,ABt)))
            return (EB, AB)

    def refine(self, l, NBOXMAX=50, VOLMIN=0.001,RSS=True):

        """refine the l layer of the CLA

        Refine the l layer of the CLA  until the maximum number of boxes (NBOXMAX) or the minimal volume of boxes (VOLMIN) has been reached.

        Once the CLA has been merged, this method aims to enclose the solution thanks to an octree/quadtreee process

        self.dlayer[l][0] : LBox which contains boxes inside of the validity area (VA)
        self.dlayer[l][1] : LBox which contains ambiguous boxes (partially inside of the VA == at least 1 edge inside the VA)
        All boxes partially inside of the VA are divided into octants. Each octants are tested into the self.valid.


        
        Parameters
        ----------

        l : the layer number
        NBOXMAX : the maximum number of obtained boxes
        VOLMIN :  the minimum volume achievable by the obtained boxes


        Returns
        -------

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

        EB, AB = self.valid_v(lv, l,RSS=RSS)
        del lv
        self.erronous.append(self.erro)

        nbox = len(EB)
        nboxamb = len(AB)

#        print nbox
#        print nboxamb
        # if all boxes are out of the VA
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

                self.refine(l,NBOXMAX, VOLMIN,RSS)

            else:

                self.iter = 0
                self.dlayer[l - 1] = self.dlayer[l]
                    # unstack to a lower the layer
                l = l - 1

                assert l >= 0, pdb.set_trace()

                self.refine(l,NBOXMAX, VOLMIN,RSS)

        # if it exists at least a box ambiguous or not in the VA
        else:

            if (nbox != 0 and nboxamb == 0):
                self.FINISHED = 1

            # Update EB
            if len(EB) != 0:
                self.dlayer[l][0].append_l(LBoxN(B.box[EB], ndim=self.ndim))

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
                self.refine(l,NBOXMAX, VOLMIN,RSS)
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

        


        Parameters
        ----------

        l       : layer number to observe. If -1 estimation is made on the highest available layer. default = -1
        amb     : display ambiguous boxes. default = false
        sc      : display all constraint or give a list with the constrinat number to observe ex: [0,1,3]. default 'all'

        Returns
        -------        

        Nothing but calls a geomview instance


        """
        Nc = self.Nc
        filename = basename + "/geom/cla.list"
        fd = open(filename, "w")
        fd.write("LIST\n")
        par = self.parmsh
        

        if l == -1:
            if sc == 'all':
                for c in self.c:
                    if c.runable:
                        c.parmsh['display'] = False
                        c.parmsh['scene'] = False
                        # if constrinat boxes has to be displayed 
                        if par['constr_boxes']:
                            c.parmsh['boxes'] = False
                        else :
                            c.parmsh['boxes'] = True
                        fname = c.show3()
                        fd.write("{<" + fname + ".list}\n")

            else:
                try:
                    for vsc in sc:
                        if self.c[vsc].runable:
                            self.c[vsc].parmsh['display'] = False
                            self.c[vsc].parmsh['scene'] = False
                        if par['constr_boxes']:
                            self.c[vsc].parmsh['boxes'] = False
                        else :
                            self.c[vsc].parmsh['boxes'] = True
                            fname = self.c[vsc].show3()
                            fd.write("{<" + fname + ".list}\n")
                except:
                    if self.c[sc].runable:
                        self.c[sc].parmsh['display'] = False
                        self.c[sc].parmsh['scene'] = False
                        if par['constr_boxes']:
                            self.c[sc].parmsh['boxes'] = False
                        else :
                            self.c[sc].parmsh['boxes'] = True
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

        


        Parameters
        ----------

        c       : contraint number in the self.c list
        d       : an array of vertex

        Returns
        -------

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
            S = (-self.c[c].sstd * np.log(10)) / (-10 * self.c[c].model.rssnp)
            M = ((self.c[c].model.PL0 - self.c[c].value) *
                 np.log(10)) / (10 * self.c[c].model.rssnp)
            v = 1 / (d * S * np.sqrt(2 * np.pi)) * np.exp(
                -(((np.log(d) - M) ** 2) / (2. * (S ** 2))))

#                       std = self.c[c].sstd#10**(self.c[c].model['RSSnp']/20.)
##                      mean = self.c[c].range
#                       mean = np.log(self.c[c].range)+std**2
#                       v = 1/(d*np.sqrt(2*np.pi))*np.exp(-(np.log(d)-mean)**2/(2*std**2))
        return(v)

#    def gapdetect(self, l, dlindx):
#        """basic gap detection

#        Detects if separated clusters of boxes are observables. his situation is usual in under determined estimation.
#        This only test on each axis if all boxes are contiguous. If not, a gap is declared and clusters are created.


#        


#        Parameters
#        ----------
#                l       : layer numbero
#                dlindx  : select the boxes type ( from self.dlayer) for gap detection 0=enclose or 1=ambigous boxes

#        Return
#        ------
#                clust   : a list of array. each array contains boxes from the same cluster
#                axis    : axis/axes where gap has/have been detectes

#        """

#        gcoord = []
#        axis = np.zeros(self.ndim, dtype='int8')
#        clust = []
##        c2={}
#        for i in range(self.ndim):
#            uni, inv, idd = np.unique(self.dlayer[l][dlindx]
#                                      .bd[:, i], return_inverse=True, return_index=True)
##                       uni,inv,idd =np.unique(self.dlayer[l][dlindx].ctr[:,i],return_inverse=True,return_index=True)

#            slope = np.diff(np.diff(uni))

##            if len(slope) != 0:
#            if len(slope) >1:
#                if abs(np.min(slope)) > 1e-9:
##                    c2[i]=[]
#                    gidx = np.nonzero(np.min(slope) == slope)[0]
##                                       print 'GAP DETECTED in AXIS',i
#                    axis[i] = 1

#                    try:
#                        # divisé par 2 pour pouvoir aveir les index  de cluster comme les centre des box
#                        clust.append(np.nonzero(uni[gidx[0]] < self.dlayer[l]
#                                                [dlindx].bd[:, i])[0] / 2)
#                        clust.append(np.nonzero(uni[gidx[0]] > self.dlayer[l]
#                                                [dlindx].bd[:, i])[0] / 2)
##                        c2[i].append(np.nonzero(uni[gidx[0]] < self.dlayer[l]
##                                                [dlindx].bd[:, i])[0] / 2)
##                        c2[i].append(np.nonzero(uni[gidx[0]] < self.dlayer[l]
##                                                [dlindx].bd[:, i])[0] / 2)
#                    except:
#                        pdb.set_trace()
#            else:
#                clust = []
#        if clust !=[]:
#            pdb.set_trace()
#        return clust, axis


    def gapdetect(self, l, dlindx):
        """basic gap detection

        Detects if separated clusters of boxes are observables. his situation is usual in under determined estimation.
        This only test on each axis if all boxes are contiguous. If not, a gap is declared and clusters are created.


        

        Parameters
        ----------

        l       : layer number
        dlindx  : select the boxes type ( from self.dlayer) for gap detection 0=enclose or 1=ambigous boxes

        Returns
        ------

        clust   : a list of array. each array contains boxes from the same cluster
        axis    : axis/axes where gap has/have been detectes

        Example
        -------

        >>> from pylayers.location.geometric.constraints.cla import *
        >>> from pylayers.location.geometric.constraints.toa import *
        >>> from pylayers.location.geometric.constraints.exclude import *
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np

        >>> a=np.array(([1,0,0]))
        >>> b=np.array(([10,0,0]))
        >>> nodes=np.array(([-10,10],[-10,10],[-1,1]))
        >>> n= np.array((5,5,0))
        >>> d1=np.sqrt(np.sum((a-n)**2))
        >>> d2=np.sqrt(np.sum((b-n)**2))
        >>> T1=TOA(id=1,value=d1/0.3,std=0.5,p=a)
        >>> T2=TOA(id=2,value=d2/0.3,std=0.5,p=b)
        >>> E=Exclude(nodes.T)
        >>> T1.runable=True
        >>> T2.runable=True
        >>> C=CLA()
        >>> C.append(T1)
        >>> C.append(T2)
        >>> C.append(E)
        >>> C.merge2()
        >>> C.refine(C.Nc)
        >>> C.gapdetect(C.Nc,1)
        

        """
        gcoord = []
        axis = np.zeros(self.ndim, dtype='int8')
        clust = []
        c2={}
        axis=np.zeros(self.ndim, dtype='int8')

        for i in range(self.ndim):
            # find all begining point on axis i
            uA,iuA=np.unique(self.dlayer[l][dlindx].bd[::2,i],return_index=True)
            # find all ending point on axis i
            uB,iuB=np.unique(self.dlayer[l][dlindx].bd[1::2,i],return_index=True)
            # remove 1st point in uA
            uAA = uA[1:]
            iuAA = iuA[1:]
            # remove last point in uA
            uBB = uB[:-1]
            iuBB = iuB[:-1]

#            u=[]
#            # find center of all these segment  
#            [u.append((uA[k]+uA[k+1])/2) for k in range(len(uA)-1) ]

#            # get all center of the boxes
#            C=self.dlayer[l][dlindx].ctr[:,i]
#            v=np.unique(C)


            # if no gap, all begining point must also be ending point, otherwise,
            # a gap exists
            igap=[]
#            [igap.append(ik) for ik,k in enumerate(u) if k not in v]
            [igap.append(ik) for ik,k in enumerate(uAA) if k not in uBB]
            if len(igap) > 1:
                igap=[igap[0]]
            # if a segment has a center which is not a box center , there is a gap
            # indexes are split into 2 set
            if not len(igap) ==0:

                # in a futur version it will be more convenient to stock each 
                # detected cluster in a given axis with a dictionary as the given
                # axis as a key.
#               c2[i].append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]<=cm[igap]))
#               c2[i].append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]>cm[igap]))
#                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]<=gap)[0]/2)
#                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]>gap)[0]/2)

                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]<=uA[igap])[0]/2)
                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]>uA[igap])[0]/2)
                axis[i]=1
#            else :
#                clust = []
        return clust,axis


    def gapdetect2(self, l, dlindx):
        """basic gap detection

        Detects if separated clusters of boxes are observables. his situation is usual in under determined estimation.
        This only test on each axis if all boxes are contiguous. If not, a gap is declared and clusters are created.
        requiere  pyinterval class 

        

        Parameters
        ----------

        l       : layer number
        dlindx  : select the boxes type ( from self.dlayer) for gap detection 0=enclose or 1=ambigous boxes

        Return
        ------

        clust   : a list of array. each array contains boxes from the same cluster
        axis    : axis/axes where gap has/have been detectes

        Example
        -------

        >>> from pylayers.location.geometric.constraints.cla import *
        >>> from pylayers.location.geometric.constraints.toa import *
        >>> from pylayers.location.geometric.constraints.exclude import *
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np

        >>> a=np.array(([1,0,0]))
        >>> b=np.array(([10,0,0]))
        >>> nodes=np.array(([-10,10],[-10,10],[-1,1]))
        >>> n= np.array((5,5,0))
        >>> d1=np.sqrt(np.sum((a-n)**2))
        >>> d2=np.sqrt(np.sum((b-n)**2))
        >>> T1=TOA(id=1,value=d1/0.3,std=np.array((0.5)),p=a)
        >>> T2=TOA(id=2,value=d2/0.3,std=np.array((0.5)),p=b)
        >>> E=Exclude(nodes.T)
        >>> T1.runable=True
        >>> T2.runable=True
        >>> C=CLA()
        >>> C.append(T1)
        >>> C.append(T2)
        >>> C.append(E)
        >>> C.merge2()
        >>> C.refine(C.Nc)
        >>> C.gapdetect2(C.Nc,1)
        

        """
        gcoord = []
        axis = np.zeros(self.ndim, dtype='int8')
        clust = []
        c2={}
        axis=np.zeros(self.ndim, dtype='int8')

        for i in range(self.ndim):
            # reshape boxes to be compliant with interval
            Z=self.dlayer[l][dlindx].bd[:,i]
            Zr=Z.reshape(len(Z)/2,2)
            # create intervals
            I=[interval(Zr[k]) for k in range(len(Zr))]
            ii=interval()

            # gather interval
            for j in I:
                ii=ii|j
            # if a gap appears (more than a unique interval) 
            if len(ii)>1:

                # in a futur version it will be more convenient to stock each 
                # detected cluster in a given axis with a dictionary as the given
                # axis as a key.
#               c2[i].append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]<=cm[igap]))
#               c2[i].append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]>cm[igap]))
#                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]<=gap)[0]/2)
#                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]>gap)[0]/2)

                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]<=ii[0][1])[0]/2)
                clust.append(np.nonzero(self.dlayer[l][dlindx].bd[:,i]>=ii[1][0])[0]/2)
                axis[i]=1


        return clust,axis


    def min_dist(self, a, b):
        """
        OBSOLETE

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

        DEPRECATED !
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

    def estpos2(self, l=-1, amb=False,HT=False):
        """ Position estimation

        estimate position from the enclosed or/and ambibuous boxes

        


        Parameters
        ----------

        l       : Layer of the estimation. If -1 estimation is made on the highest available layer
        amb     : Use ambiguous boxes (if available) to perform the position estimation. default = False
        HT      : boolean
                True : if a cluster ppears (2 sets of distinct boxes ) an hypthesis testuing method is applied
                        in estpos2 method 
                False : no HT methos is applied 

                Hybrid positioning based on hypothesis thesting
                N. Amiot, T. Pedersen, M. Laaraiedh, B. Uguen. 
                A Hybrid Positioning Method Based on Hypothesis Testing
                ,Wireless Communications Letters, IEEE, vol.1, no.4, pp.348-351, August 2012

        Returns 
        -------

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

        if pyinterval_installed:
            clust, axis = self.gapdetect2(l, dlindx)
        else:
            clust, axis = self.gapdetect(l, dlindx)

        box_center = self.dlayer[l][dlindx].ctr
        uc = np.where(self.usable)[0]

        
        # proba computation for all center of each boxes
        for j in uc:#range(len(self.c)):
            #if self.c[j].type != 'Exclude':
            if (self.c[j].type != 'Exclude') & (self.c[j].usable):
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
        self.pecluster=[]
        if clust != []:



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
            lpc=[]
            for i in range(len(clust)):
                if len(clust) < 3:
                    clusters = clust[i]
                else:

                    if len(np.shape(count)) > 1:
                        clusters = np.intersect1d(clust[count[i,0]], clust[count[i, 1]])
                    else:
                        clusters = np.intersect1d(clust[count[0]], clust[count[1]])

                clust_vol = np.sum(np.array(self.dlayer[l][
                    dlindx].vol)[np.unique(clusters)])

                if len(clusters) != 0:
                    mp = np.max(self.saveP[clusters])

                    if mps < mp:
                        mps = mp
                        estclu = clusters


                itoas=np.where(np.array(self.type)=='TOA')[0]
                if clust_vol != 0 and len(itoas) == 2:
                    lclust.append(clusters)
                    pc = np.sum(np.array(self.dlayer[l][dlindx].ctr)[np.unique(clusters)], axis=0) / len(np.unique(clusters))
                    lpc.append(pc)
                    # verifier que les contraintes utilisées sont les bonne ( ce n'est pas le cas)
                    # ne marche que si 2 constriantes genere le cluster ( a robustifier)   
                    pu = np.where(self.usable)[0]
                    # try:
                    #     dd.append(np.sqrt(np.sum((pc - self.c[itoas[0]].p) ** 2)))
                    # except:
                    #     dd.append(np.sqrt(np.sum((pc - self.c[itoas[1]].p) ** 2)))
                    # print pc
                    try:
                        dd.append(np.sqrt(np.sum((pc - self.c[itoas[0]].p) ** 2)))
                    except:
                        dd.append(np.sqrt(np.sum((pc - self.c[itoas[1]].p) ** 2)))
                    print pc
            #                       try:
#                               vmax=[]
#                               for i in range(len(lclust)):
#                                       vmax.append(np.max(poids[np.unique(lclust[i])]))
#                               peindx = np.nonzero(poids==max(vmax))[0][0]
#                               self.pe = self.dlayer[l][dlindx].ctr[peindx]

            if HT:

                print "enter in HT processing"
                try:

                    # for now, it is supposed that all RSS share the same model
                    rssvalues=[]
                    icr=np.where(np.array(self.type)=='RSS')[0]
                    for irss in range(len(icr)):
                        d0=np.sqrt(np.sum((self.c[icr[irss]].p-lpc[0])**2))
                        d1=np.sqrt(np.sum((self.c[icr[irss]].p-lpc[1])**2))
                        rssvalues.append(self.c[icr[irss]].value)
                        try:
                            drss= np.vstack((drss,np.array((d0,d1))))
                        except:
                            drss= np.array((d0,d1))
                    if len(np.shape(drss))==1:
                        drss=drss.reshape(1,2)    

                    M = (((-self.c[icr[0]].model.PL0 - self.c[icr[0]].value) * np.log(10) ) / (10. * self.c[icr[0]].model.rssnp))
                    PL0= -self.c[icr[0]].model.PL0 
                    NP = self.c[icr[0]].model.rssnp
                    
                    mu1=PL0-10*NP*np.log10(drss[:,0])
                    mu2=PL0-10*NP*np.log10(drss[:,1])
                    sig=self.c[icr[0]].model.sigrss
                    values=np.array((rssvalues))
                    LT=np.sum(1/(2.*sig**2)*(mu2**2-mu1**2))
                    RT=np.sum((1/(1.*sig))*values*(mu1-mu2))


                    # LL = np.log(dd[1] / dd[0]) * (1 + np.log(dd[0] * dd[1]) - 2 * M)


                    # if LL > 0:
                    if LT>RT:
    #                                       vmax = np.max(poids[np.unique(lclust[0])])
    #                                       peindx=np.nonzero(poids[vmax]==poids)[0][0]
    #                                       self.pe = self.dlayer[l][dlindx].ctr[np.unique(lclust[0])[peindx]]

                        #if LL>0  cluster 0 is selctionned and tits centroids is chosen as position estimation

                        self.pe = np.mean(self.dlayer[l][dlindx].ctr[
                            np.unique(lclust[0])], axis=0)
                        print "HT processing done"
                        pestdmax = np.max(self.dlayer[l][
                            dlindx].ctr[np.unique(lclust[0])])
                        pestdmin = np.min(self.dlayer[l][
                            dlindx].ctr[np.unique(lclust[0])])
                        self.pestd = pestdmax - pestdmin


                    else:


                        #if LL<0  cluster 1 is selctionned and tits centroids is chosen as position estimation

                        self.pe = np.mean(self.dlayer[l][dlindx].ctr[
                            np.unique(lclust[1])], axis=0)
                        pestdmax = np.max(self.dlayer[l][
                            dlindx].ctr[np.unique(lclust[1])])
                        pestdmin = np.min(self.dlayer[l][
                            dlindx].ctr[np.unique(lclust[1])])
                        self.pestd = pestdmax - pestdmin


                # if HT fail for some reasons , a classical position estimation  is performed 
                except:
                    print "!!!!! HT FAIL !!!!!!!"
                    print "2 first constraint of CLA have to be TOA and others RSS in order to use HT"

                    if np.sum(poids) > 0.:
                        self.pe = np.sum(poids * self.dlayer[l][dlindx]
                            .ctr.T, axis=1) / np.sum(poids)
                    else:
                        self.pe = np.sum(self.dlayer[l][dlindx].ctr, axis=0) / \
                            len(self.dlayer[l][dlindx].ctr)
                    pestdmax = np.max(self.dlayer[l][dlindx].bd, axis=0)
                    pestdmin = np.min(self.dlayer[l][dlindx].bd, axis=0)
                    self.pestd = pestdmax - pestdmin

 



            # if no HT
            else:
                if np.sum(poids) > 0.:
                    self.pe = np.sum(poids * self.dlayer[l][dlindx]
                        .ctr.T, axis=1) / np.sum(poids)
                else:
                    self.pe = np.sum(self.dlayer[l][dlindx].ctr, axis=0) / \
                        len(self.dlayer[l][dlindx].ctr)
                pestdmax = np.max(self.dlayer[l][dlindx].bd, axis=0)
                pestdmin = np.min(self.dlayer[l][dlindx].bd, axis=0)
                self.pestd = pestdmax - pestdmin


            # store the centroid of clusters into self.peclsuter
            for cl in lclust:	
                self.pecluster.append(np.mean(self.dlayer[l][dlindx].ctr[
                np.unique(cl)], axis=0))


        # if not cluster
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
            self.pecluster=[self.pe]

