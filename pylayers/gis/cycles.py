# -*- coding:Utf-8 -*-
#
# Class Cycles
# Class Cycle
# 
# Author : B.Uguen 
# 
import doctest 
import networkx as nx
import numpy as np 
from  pylayers.util import geomutil as geu
from  pylayers.util import pyutil as pyu
import matplotlib.pyplot as plt
import pdb

class Cycles(object):
    """ List of cycles

    Attributes
    ----------
    LC    : List of Cycles
    Nc    : number of Cycles
    Area1 : list of cycles area

    Methods
    -------
        inclusion 
        simplify
    """
    def __init__(self,LC,Gs):
        """
        LC : list of Cycle
        Gs : Stucture Graph of a Layout
        """
        self.Gs = Gs
        self.LC = LC
        self.Nc = len(LC)
        self.Gt = {}
        self.inclusion()
        self.Area1 = []
        for k in range(len(LC)):
            self.Area1.append(LC[k].area)
        self.simplify()

    def info(self):
        """ get info from Cycles object
        """
        print "Number of cycles : ",self.Nc
        for k in self.Gi:
            print k , self.Gi.neighbors(k)

    def inclusion(self):
        """

        Notes
        -----

        update cycle inclusion graph
        inclusion is not a reciprocal relation that is the reason why 
        we loop twice over the entire cycles

        Algorithme :
            Pour tous les cycles
                Ck
                Pour tous les cycles
                      Cl
                      Si c'est le meme cycle on passe
                      Sinon
                    check inclusion of Ck in Cl
                    If Inclusion detected
                        update graphe
                        d'inclusion Cl in Ck

        Attention :
            Le test d'inclusion de cycle de fonctionne que si les
            cycles ont un brin en commun. Si deux cycles ne
            partagent aucun brin le test d'inclusion renvoie False
            meme si ce n'est pas le cas

        Les relations d'inclusion reelles are recovered par
        transitivity de la relation d'inclusion

        Cette fonction met a jour un Digraph
        """
        self.Gi     = nx.DiGraph()
        self.contained_in = {}
        self.contained_by = {}
        self.l1     = {}
        self.l2     = {}
        self.Area2  = []
        for k in range(self.Nc):
            ck = self.LC[k]
            self.Gi.add_node(k)
            Area     = 0
            for l in range(self.Nc):
                if l != self.Nc:
                    cl = self.LC[l]
                    incl,path = ck.inclusion(cl)
                    if incl:
                        try:
                            self.contained_by[k].append(l)
                        except:
                            self.contained_by[k]=[l]
                        try:
                            self.contained_in[l].append(k)
                        except:
                            self.contained_in[l]=[k]
                        Area = abs(cl.area) + Area
            self.Area2.append(Area)
            try:
                self.l1[k]       = len(self.contained_by[k])
            except:
                self.l1[k]       = 0

        #
        # l2 : dictionnaire  longueur d'inclusion <--> cycles inclus
        #
        for k in self.l1.keys():
            l = self.l1[k]
            try:
                self.l2[l].append(k)
            except:
                self.l2[l]=[k]
        # pour toutes les longueurs d'inclusion en partant de 0
        for l in self.l2.keys():
            # liste des no de cycle de longueur d'inclusion l
            lnode = self.l2[l]
            for n in lnode:
                if self.contained_in.has_key(n):
                    v = self.contained_in[n]
                    #print "cycle :",n," is contained in cycles ",v
                    cmin = self.Nc
                    for cy in v:
                        # cherche le cycle incluant de
                        # longueur minimale
                        if self.l1[cy] < cmin:
                            c    = cy
                            cmin = self.l1[cy]
                    self.Gi.add_edge(c,n)

    def simplify(self):
        """ simplify Cycles

        Notes 
        ------

            Pour tous les cycles
                Ck
                Recherche les voisins du cycle k

                Si pas de voisins
                    Le cycle est simple ( non incluant)
                    Mise a jour du graphe de topologie Gt
                Sinon
                    Le cycle est incluant
                    Construit la liste M de tous les cycles inclus dans Ck

                    Simplification de Ck avec M
                    Mise a jour du graphe de topologie Gt

        """
        self.Gt     = nx.Graph()
        self.Gt.pos = {}
        #
        # Simplification des cycles incluants
        #
        # Ce n'est pas le bon algo
        #
        NC = []
        for k in self.Gi.node.keys():
            ck = self.LC[k]
            l  = nx.neighbors(self.Gi,k)
            if len(l)!=0: # cycle incluant
                cl = self.LC[l[0]]
                C  = ck.untie(self.Gs,cl)
                #AC = abs(C.area)
                #Al = abs(cl.area)
                #Ak = abs(ck.area)
                #if (AC+Al!=Ak):
                #    print "Ak    :",Ak
                #    print "AC+Al :",AC + Al
                #    print "diff :",AC + Al - Ak
                NC.append(C)
                self.Gt.add_node(k,vnodes=C.cycle)
                self.Gt.pos[k]=C.g
                #C.show()
            else:  # cycle simple
                NC.append(ck)
                self.Gt.add_node(k,vnodes=ck.cycle)
                self.Gt.pos[k]=ck.g
                #ck.show('red')
        #self.LC = NC
        MC = []
        for j in self.l2:
            for k in self.l2[j]:
                ck = self.LC[k]
                if j>0:
                    l  = nx.neighbors(self.Gi,k)
                    C  = ck
                    # simplification cycles inclus multiples
                    for nl in l:
                        cl = self.LC[nl]
#                        if (j == 16) and (k ==48) and (nl == 52):
#                            pdb.set_trace()
                        C  = C.untie(self.Gs,cl)


                    MC.append(C)
                    self.Gt.add_node(k,vnodes=C.cycle)
                    self.Gt.pos[k]=C.g
                else:
                    MC.append(ck)
                    self.Gt.add_node(k,vnodes=ck.cycle)
                    self.Gt.pos[k]=ck.g

        self.LC = MC
class Cycle(object):
    """ Graph cycle class

    Attributes
    ----------
    vertices : np.array
    edges   : np.array
    cycle   : np.array
    """
    def __init__(self,G,cycle):
        # This is to obtained an ordered cycle
        self.G  = G.subgraph(cycle)
        cycle   = nx.algorithms.cycles.cycle_basis(self.G)[0]
        self.G  = G.subgraph(cycle)
        self.G.pos = {}
        for c in cycle:
            self.G.pos[c] = G.pos[c]
        self.cycle     = np.array(cycle)
        self.size      = len(self.cycle)
        self.vertices  = self.cycle[np.nonzero(np.array(self.cycle)<0)[0]]
        self.edges     = self.cycle[np.nonzero(np.array(self.cycle)>0)[0]]
        self.Nv    = len(self.vertices)
        self.Ne    = len(self.edges)
        self.p     = np.array([])
        for v in  self.vertices:
            try:
                self.p  = np.hstack((self.p,np.array(self.G.pos[v]).reshape(2,1)))
            except:
                self.p  = np.array(self.G.pos[v]).reshape(2,1)

        self.area = geu.SignedArea(self.p)
        self.g    = geu.Centroid(self.p)

    def info(self):
        """
        """
        print self.cycle
        print self.vertices
        print self.edges
        print self.size
        print self.area

    def flip(self):
        """ flip the Cycle
        """
        self.cycle    = self.cycle[::-1]
        self.vertices = self.vertices[::-1]
        self.edges    = self.edges[::-1]
        self.p        = self.p[:,::-1]
        self.area     = geu.SignedArea(self.p)

    def inclusion(self,Cy):
        """  Return True if self includes Cy

        Parameters
        ----------
        Cy  : Cycle 

        Returns
        -------
        boolean  : answer to the inclusion question 
        path     : path 

        Notes 
        -----

        Inclusion implies :

            1. The area of the including cycle is larger than the area of included cycle
            2. The 2 cycles have at least  a common segment 
            3. When the  common segment is travelled with the same orientation both cycle SignedArea have the same sign

        """
        if abs(self.area)>abs(Cy.area):
            flip,path = self.intersect(Cy)
            if len(path)>0:
                s1 = Cy.area*self.area
                if flip:  # sens oppose
                    if s1>0:
                        return False,path
                    else:
                        return True,path
                else:     # meme sens
                    if s1>0:
                        return True,path
                    else:
                        return False,path
            else:
                return False,np.array([])
        else:
            return False,np.array([])

    def simplify(self,G,L):
        """ 
        Parameters
        ----------

        G : 
        L : list of cycles

        """
        T = self
        for Cy in L:
            b,p = T.inclusion(Cy)
            if b:
                T = T.untie(G,Cy)
        return(T)

    def intersect(self,Cy):
        """ Find intersection with an other cycle

        Parameters
        ----------
        Cy : Cycle 

        Returns
        -------

        flip  :
        brin  :

        Notes 
        ------

        flip,path = C1.intersect(C2)

        if flip == True --> The 2 cycles have a reverse travel
        direction for the common path

        path is the common path along the cycle

        """
        e1  = self.edges
        e2  = Cy.edges
        v1  = self.vertices
        v2  = Cy.vertices
        u   = np.intersect1d(e1,e2)
        Ni  = len(u)
        #print Ni," segments en commun"
        if Ni>0:
            if Ni>1:
                brin1 = self.cycle
                brin2 = Cy.cycle
                ####
                #### Faux
                ####
                #tk1 = np.array([],dtype=int)
                #tk2 = np.array([],dtype=int)
                #for ku in u:
                #    ik1 = np.nonzero(e1==ku)[0]
                #    ik2 = np.nonzero(e2==ku)[0]
                #    tk1 = np.hstack((tk1,ik1))
                #    tk2 = np.hstack((tk2,ik2))
                #tk1.sort()
                #tk2.sort()
                #brin1 = e1[tk1]
                #brin2 = e2[tk2]
                #if Ni==2:
                #    if prod(brin1==brin2):
                #        flip = False
                #    else:
                #        flip = True
                #else:
                if max(pyu.corrcy(brin1,brin2))>max(pyu.corrcy(brin1,brin2[::-1])):
                    flip = False
                else:
                    flip = True
                brin = u
            else:
                v   = np.intersect1d(v1,v2)
                tk1 = np.array([],dtype=int)
                tk2 = np.array([],dtype=int)
                for kv in v:
                    ik1 = np.nonzero(v1==kv)[0]
                    ik2 = np.nonzero(v2==kv)[0]
                    tk1 = np.hstack((tk1,ik1))
                    tk2 = np.hstack((tk2,ik2))
                dtk1 = tk1[1]-tk1[0]
                dtk2 = tk2[1]-tk2[0]
                if abs(dtk1)>1:
                    dtk1 = -dtk1
                if abs(dtk2)>1:
                    dtk2 = -dtk2
                if dtk1*dtk2>0:
                #    print "meme sens"
                    flip = False
                else:
                #    print "sens oppose"
                    flip = True

                #tk1.sort()
                #tk2.sort()
                #brin1 = v1[tk1]
                #brin2 = v2[tk2]
                #print brin1,brin2
                #if max(corrcy(brin1,brin2))>max(corrcy(brin1,brin2[::-1])):
                #    flip = False
                #else:
                #    flip = True
                brin = u
        else:
            flip  = False
            brin = np.array([])

        return(flip,brin)


    def untie(self,G,Cy):
        """
          untie(G,Cy)

          Parameters
          ----------
          G  : 
          Cy :  Cycle 

          Notes 
          -----

          Simplify 2 cycles which are in a mutual relation of inclusion 

        """
        e1 = self.edges
        e2 = Cy.edges
        incl,path = self.inclusion(Cy)
        if incl:
            ee1   = e1[~np.in1d(e1,path)]
            ee2   = e2[~np.in1d(e2,path)]
            r     = np.hstack((ee1,ee2))
            cycle = np.array([],dtype=int)
            for kt in r:
                nb = G.neighbors(kt)
                cycle = np.hstack((cycle,nb[0],kt,nb[1]))
            # numpy old version : cycle = np.unique1d(cycle)
            cycle = np.unique(cycle)
            NCy   = Cycle(G,cycle)
            return(NCy)


    def corrcy(self,Cy,flip=False):
        """ cyclic matching correlation

        Parameters
        ----------

        Cy : Cycle
        flip : boolean 
            default = False 

        """
        e1 = self.edges
        e2 = Cy.edges
        u  = np.intersect1d(e1,e2)
        Ni = len(u)
        print Ni," intersections"
        v1 = self.vertices
        v2 = Cy.vertices
        n1 = len(e1)
        n2 = len(e2)
        tk1  = np.array([])
        tk2  = np.array([])
        brin = np.array([])
        if n1>n2:
            for k in range(n1):
                c1k   = np.hstack((e1[k::],e1[0:k]))
                diff1 = c1k[0:n2] - e2[::-1]
                diff2 = c1k[0:n2] - e2

                u1    = np.nonzero(diff1==0)[0]
                u2    = np.nonzero(diff2==0)[0]
                l1    = len(u1)
                l2    = len(u2)
                if l1==Ni:
                    brin = e1[::-1][u1]
                    flip = True
                if l2==Ni:
                    brin = e2[u2]
                    flip = False
                tk1   = np.hstack((tk1,l1))
                tk2   = np.hstack((tk2,l2))

            if (Ni>0) & len(brin)==0:
                mt1   = max(tk1)
                ut1   = np.nonzero(tk1==mt1)[0][0]
                mt2   = max(tk2)
                ut2   = np.nonzero(tk2==mt2)[0][0]

                if (mt1>mt2):
                    brin = e1[ut1:ut1+mt1]
                    flip = True
                else:
                    brin = e1[ut2:ut2+mt2]
                    flip = False


        else:
            for k in range(n2):
                c2k   = np.hstack((e2[k::],e2[0:k]))
                diff1 = c2k[0:n1]- e1[::-1]
                diff2 = c2k[0:n1]- e1
                u1    = np.nonzero(diff1==0)[0]
                u2    = np.nonzero(diff2==0)[0]
                l1    = len(u1)
                l2    = len(u2)
                if l1==Ni:
                    brin = e2[::-1][u1]
                    flip = True
                if l2==Ni:
                    brin = e2[u2]
                    flip = False
                tk1   = np.hstack((tk1,l1))
                tk2   = np.hstack((tk2,l2))

            if (Ni>0) & len(brin)==0:
                mt1   = max(tk1)
                ut1   = np.nonzero(tk1==mt1)[0][0]
                mt2   = max(tk2)
                ut2   = np.nonzero(tk2==mt2)[0][0]
                if (mt1>mt2):
                    brin = e2[ut1:ut1+mt1]
                    flip = True
                else:
                    brin = e2[ut2:ut2+mt2]
                    flip = False


        return(flip,brin,tk1,tk2)


    def show(self,color='b'):
        nx.draw_networkx_edges(self.G,self.G.pos,width=2,edge_color=color)
        plt.draw()

if __name__ == "__main__":
    doctest.testmod()
