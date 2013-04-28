# -*- coding:Utf-8 -*-
#
# Class Cycles(list)
# Class Cycle(Graph)
#
# Author : B.Uguen
#  2nd implementation : April 2013
#
import doctest
import networkx as nx
import numpy as np
from  pylayers.util import geomutil as geu
from  pylayers.util import pyutil as pyu
import matplotlib.pyplot as plt
import logging
import pdb

class Cycles(nx.DiGraph):
    """ Graph of Cycles

    Methods
    -------
        inclusion
        simplify
        decompose

    """
    def __init__(self):
        """
            graph of cycle
        """
        nx.DiGraph.__init__(self)

    def __repr__(self):
        s = ''
        s = s + 'Number of cycles :' +str(len(self.node)) + '\n'
        for k in self:
            s = s + str(self.node[k]['cycle']) + '\n'
        return(s)

    def info(self):
        """ get info from Cycles object
        """
        print "Number of cycles : ",len(self)
        for k in self.Gi:
            print k , self.Gi.neighbors(k)

    def inclusion(self,full=True):
        """  Inform the inclusion relations between cycles

        Notes
        -----

        update cycle inclusion graph

        inclusion is not a reciprocal relation that is the reason why
        we loop twice over the entire cycles list

        Algorithm :
            for all cycles ck
                for all cycles cl
                    if ck <> cl
                        check inclusion of Ck in Cl
                        If Inclusion detected
                            update graphe

        Warning
        -------
            The inclusion test works only for cycles which share
            at least one portion of their cycle.
            If one cycle is completely included in the other the return
            will be false even if there is inclusion

        actual inclusion relationships are recovered by transitivity of inclusion

        This function updates a Digraph of inclusion

        """
        self.contained_in = {}
        self.contained_by = {}
        self.l1     = {}
        self.l2     = {}
        self.Area2  = []
        #
        # the list of cycles is either the whole set of nodes
        # or only connected nodes (nodes with degree > 0)
        #
        if full:
            lcyt= [ k for k in self]
        else:
            lcyt = [ k for k in self if nx.degree(self)[k]>0]
        #print lcyt
        #
        # loop on all cycle in list of cycles
        for k in lcyt:
            #Areak = 0
            ck = self.node[k]['cycle']
            for l in self:
                if l!=k:
                    cl = self.node[l]['cycle']
                    # test inclusion of cl(small) by ck(big)
                    # non reciprocal
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
                        # update Area iteration k
                        #Areak = abs(cl.area) + Areak
            #self.Area2.append(Areak)
            try:
                self.l1[k]       = len(self.contained_by[k])
            except:
                self.l1[k]       = 0

        #
        # l1 : dictionnary
        #      key : cycle index   value : number of included cycles
        #
        # l2 : dictionnary
        #      key : number of included cycles  value : cycle index
        #

        for k in self.l1.keys():
            l = self.l1[k]
            try:
                self.l2[l].append(k)
            except:
                self.l2[l]=[k]
        # for all inclusion length (starting in 0) 
        for l in self.l2.keys():
            # retrieve the list of included cycles of length l
            lnode = self.l2[l]
            for cy1 in lnode:
                if self.contained_in.has_key(cy1):
                    v = self.contained_in[cy1]
                    logging.info("cycle %d is included in cycles %s",cy1,str(v)) 
                    cmin = len(self) # ???
                    for cy2 in v:
                        # seek the including cycle of minimal length
                        if self.l1[cy2] < cmin:
                            c    = cy2
                            cmin = self.l1[cy2]
                    # c includes cy1
                    self.add_edge(c,cy1)


    def decompose2(self):
        """ recursive function to transform the cycle basis

        """
        #
        # root if degree == number of successors
        #
        plt.ion()
        self.lcyroot =[k for k in nx.degree(self)
                  if (((nx.degree(self)[k]==len(self.neighbors(k)))) and
                  (nx.degree(self)[k]>0))]

        #for ncyroot in lcyroot:
        while len(self.lcyroot)>0:
            # get big cycle
            ncyroot = self.lcyroot.pop()
            cybig = self.node[ncyroot]['cycle']
            #cybig.show('b')
            # get the list of the successor cycle indices
            lninc = nx.neighbors(self,ncyroot)
            for ninc in lninc:
                cyinc = self.node[ninc]['cycle'] # cycle included
                cysmall = cybig.split(cyinc) # small cycle
                if cysmall !=None:
                    cybig = cysmall
            #cyinc.show('g')
            #
            # split cycle :  cybig = cyinc \cup cysmall
            #
            if cysmall != None:
                self.node[ncyroot]['cycle']=cysmall
                self.pos[ncyroot]=tuple(cysmall.g)
                cysmall.show('r')
                plt.show()
                for ninc in lninc:
                    print ncyroot,ninc
                    self.remove_edge(ncyroot,ninc)

            #
            # recursive call of decompose
            #
            self = self.decompose2()
            #
            # Two lines to kill the monster
            #
            if len(self.lcyroot)==0:
                return(self)
        return(self)

    def decompose(self):
        """ recursive function to transform the cycle basis

        """

        self.inclusion(True)
        #
        # root if degree == number of successors
        #
        self.lcyroot =[k for k in nx.degree(self)
                  if (((nx.degree(self)[k]==len(self.neighbors(k)))) and
                  (nx.degree(self)[k]>0))]
        print self.lcyroot

        #for ncyroot in lcyroot:
        while len(self.lcyroot)>0:
            # get big cycle
            ncyroot = self.lcyroot.pop()
            cybig = self.node[ncyroot]['cycle']
            #cybig.show('b')
            # buid dictionnary dsuc
            #   key : cycle number
            #   value : list of successors
            lninc = nx.neighbors(self,ncyroot)
            # get the list of the successor cycle indices
            for ninc in lninc:
                cyinc = self.node[ninc]['cycle']     # cycle included
                #cyinc.show('g')
                #
                # split cycle :  cybig = cyinc \cup cysmall
                #
                cysmall = cybig.split(cyinc) # small cycle
                #self.add_node(ninc,cycle=cyinc)
                if cysmall != None:
                    self.node[ncyroot]['cycle']=cysmall
                    self.pos[ncyroot]=tuple(cysmall.g)
                    self.remove_edge(ncyroot,ninc)
                    #cysmall.show('r')
                    #
                    # recursive call of decompose
                    #
                self = self.decompose()
                #
                # Two lines to kill the monster
                #
                if len(self.lcyroot)==0:
                    return(self)
        return(self)

class Cycle(object):
    """ Graph cycle class

    Attributes
    ----------
    vertices : np.array
    edges   : np.array
    cycle   : np.array
    """
    def __init__(self,G):
        # This call to cycle_basis is to obtained an ordered cycle
        self.G  = G
        cycle = nx.algorithms.cycles.cycle_basis(self.G)[0]
        self.G  = G.subgraph(cycle)
        self.G.pos = {}
        self.G.pos.update({ node : G.pos[node] for node in cycle})
        #for node in cycle:
        #    self.G.pos[node] = G.pos[node]
        self.cycle = np.array(cycle)
        self.size      = len(self.cycle)
        self.vertices  = self.cycle[np.nonzero(np.array(self.cycle)<0)[0]]
        self.edges     = self.cycle[np.nonzero(np.array(self.cycle)>0)[0]]
        self.Nv    = len(self.vertices)
        self.Ne    = len(self.edges)
        for v in  self.vertices:
            try:
                pt =  np.array(self.G.pos[v]).reshape(2,1)
            except:
                pdb.set_trace()
            try:
                self.p  = np.hstack((self.p,pt))
            except:
                self.p  = pt

        self.area = geu.SignedArea(self.p)
        self.g    = geu.Centroid(self.p)


    def __repr__(self):
        s = ''
        s = s + 'cycle nstr'+str(self.cycle)+'\n'
        s = s + 'point number '+str(self.Nv)+'\n'
        s = s + 'segment number '+str(self.Ne)+'\n'
        s = s + 'area : '+ str(self.area)+'\n'
        s = s + 'centroid : '+ str(self.g)+'\n'
        #s = s + 'G :'+str(self.G.nodes())
        return(s)



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

    def simplify(self,L):
        """
        Parameters
        ----------

        L : list of cycles

        """
        T = self
        for Cy in L:
            b,p = T.inclusion(Cy)
            if b:
                T = T.split(Cy)
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


    def split(self,cyin):
        """ split(Cy)

          Parameters
          ----------
          cyin :  input cycle

          Returns
          -------
          cyout : list of output cycles

          Notes
          -----

          Simplify 2 cycles which are in a mutual relation of inclusion

        """
        e1 = self.edges
        e2 = cyin.edges
        Gunion = nx.compose(self.G,cyin.G)
        # inform nodes position
        Gunion.pos = {}
        Gunion.pos.update({node:self.G.pos[node] for node in self.G})
        Gunion.pos.update({node:cyin.G.pos[node] for node in cyin.G})
        incl,path = self.inclusion(cyin)
        if incl:
            ee1   = e1[~np.in1d(e1,path)]
            ee2   = e2[~np.in1d(e2,path)]
            r     = np.hstack((ee1,ee2))
            cycle = np.array([],dtype=int)
            for kt in r:
                nb = Gunion.neighbors(kt)
                cycle = np.hstack((cycle,nb[0],kt,nb[1]))
            cycle = np.unique(cycle)
            # extract subgraph from cycle
            G     = nx.subgraph(Gunion,cycle)
            G.pos = {}
            G.pos.update({node: Gunion.pos[node] for node in Gunion})
            cyout = Cycle(G)
            diffarea = abs(self.area)-(abs(cyin.area)+abs(cyout.area))
            if diffarea>1e-10:
                print "area error",diffarea
                pdb.set_trace()

            return(cyout)


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
        nx.draw_networkx_edges(self.G,self.G.pos,width=2,edge_color=color,alpha=0.4)
        plt.draw()

if __name__ == "__main__":
    doctest.testmod()
