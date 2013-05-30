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
import shapely.geometry as sh
import numpy as np
from  pylayers.util import geomutil as geu
from  pylayers.util import pyutil as pyu
import matplotlib.pyplot as plt
import logging
import copy
import pdb

class Cycles(nx.DiGraph):
    """ Graph of Cycles

    Methods
    -------
        inclusion
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
            lcyt = [ k for k in self]
        else:
            lcyt = [ k for k in self if nx.degree(self)[k]>0]
        #print lcyt
        #
        # loop on all cycle in list of cycles
        for k in lcyt:
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

        # check transitivity
        for ncy1 in self.contained_by:
            lcy1 = self.contained_by[ncy1]
            for ncy2 in lcy1:
                if ncy2 in self.contained_by:
                    lcy2 = self.contained_by[ncy2]
                    for ncy2 in lcy2:
                        if ncy2 not in lcy1:
                            ##print "adding ",ncy2," in ",ncy1
                            self.contained_by[ncy1].append(ncy2)

        #
        # l1 : dictionnary
        #      key : cycle index   value : number of included cycles
        #
        for k in lcyt:
            try:
                self.l1[k]       = len(self.contained_by[k])
            except:
                self.l1[k]       = 0
        # l2 : dictionnary
        #      key : number of included cycles  value : cycle index
        #

        for k in self.l1.keys():
            l = self.l1[k]
            try:
                self.l2[l].append(k)
            except:
                self.l2[l]=[k]
        #pdb.set_trace()
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


    def decompose(self):
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
            # get the list of the successor cycle indices
            lninc = nx.neighbors(self,ncyroot)
            #
            # reduce to the smallest cycle
            #
            # when a punctual contact is detected, the simplification order is important
            #
            reloop = False
            #
            # first attempt
            #
            for ninc in lninc:
                cyinc1 = self.node[ninc]['cycle'] # cycle included
                punctual,cyinc2 = cybig.split(cyinc1) # small cycle
                if punctual:
                    reloop = True
                    logging.warning("punctual contact detected proceed in reverse order")
                if cyinc2 !=None: # divide again with updated big cycle
                    cybig = cyinc2
                else:
                    pass
                    #print "None"
            #
            # second attempt if a punctual contact has been detected
            #
            if reloop:
                lninc.reverse()
                reloop=False
                for ninc in lninc:
                    cyinc1 = self.node[ninc]['cycle'] # cycle included
                    punctual,cyinc2 = cybig.split(cyinc1) # small cycle
                    if punctual:
                        logging.critical("contact detected reconsider the layout description ")
                    if cyinc2 !=None: # divide again with updated big cycle
                        cybig = cyinc2
                    else:
                        pass
                        #print "None"

            if cyinc2==None:
                cyinc2 = cybig
            #
            # split cycle :  cybig = cyinc \cup cysmall
            #
            #a1 = abs(cyinc1.area)
            #a2 = abs(cyinc2.area)
            #if a1<a2:
            #    cysmall = cyinc1
            #else:
            cysmall = cyinc2
            self.node[ncyroot]['cycle']=cysmall
            self.pos[ncyroot]=tuple(cysmall.g)
            #plt.figure()
            #plt.title(str(ncyroot))
            #cysmall.show('r')
            #plt.show()
            for ninc in lninc:
                self.remove_edge(ncyroot,ninc)

            #
            # recursive call
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
        self.update()

    def __add__(self,cy):
        """
        addition of 2 disjoint cycles is not a cycle
        """
        flip,path = self.intersect(cy)
        if len(path)>0:
            if path[0]<0:
                g1  = copy.deepcopy(self.G)
                g2  = copy.deepcopy(cy.G)
                #
                pt  = g1.pos[path[0]]
                #
                n1  = nx.neighbors(g1,path[0])
                n2  = nx.neighbors(g2,path[0])
                #
                g1.remove_node(path[0])
                g2.remove_node(path[0])
                #
                g1.add_node(-1000000)
                g1.add_node(-1000001)
                g2.add_node(-1000002)
                g2.add_node(-1000003)
                #
                g1.add_edge(n1[0],-1000000)
                g1.add_edge(n1[1],-1000001)
                g2.add_edge(n2[0],-1000002)
                g2.add_edge(n2[1],-1000003)
                #
                # Need to check an intersection for determining  proper connection
                #
                p10 = np.array(g1.pos[n1[0]]).reshape(2,1)
                p11 = np.array(g1.pos[n1[1]]).reshape(2,1)
                p20 = np.array(g2.pos[n2[0]]).reshape(2,1)
                p21 = np.array(g2.pos[n2[1]]).reshape(2,1)
                #
                boolinter = geu.intersect(p10,p20,p11,p21)[0]
                #
                Gc = nx.Graph()
                Gc = nx.compose(g1,g2)
                if boolinter:
                    Gc.add_edge(-1000000,-1000003)
                    Gc.add_edge(-1000001,-1000002)
                else:
                    Gc.add_edge(-1000000,-1000002)
                    Gc.add_edge(-1000001,-1000003)

                Gc.pos = {}
                Gc.pos[-1000000] = pt
                Gc.pos[-1000001] = pt
                Gc.pos[-1000002] = pt
                Gc.pos[-1000003] = pt
                Gc.pos.update(self.G.pos)
                Gc.pos.update(cy.G.pos)
                newcy = Cycle(Gc)
                return(newcy)
            else:
                g1  = copy.deepcopy(self.G)
                g2  = copy.deepcopy(cy.G)
                g1.remove_node(path[0])
                g2.remove_node(path[0])
                Gc = nx.Graph()
                Gc  = nx.compose(g1,g2)
                Gc.pos = {}
                Gc.pos.update(self.G.pos)
                Gc.pos.update(cy.G.pos)
                newcy = Cycle(Gc)
                newcy.update()
                return(newcy)
        else:
            return


    def update(self):
        """ update

        """
        self.size  = len(self.cycle)
        self.vertices = self.cycle[np.nonzero(np.array(self.cycle)<0)[0]]
        self.edges = self.cycle[np.nonzero(np.array(self.cycle)>0)[0]]
        self.Nv = len(self.vertices)
        self.Ne = len(self.edges)

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
        boolean  : True if Cy \subset self
        path     : path

        Notes
        -----

        Inclusion implies :

            1. The area of the including cycle is larger than the area of included cycle
            2. The 2 cycles have at least  a common segment
            3. When the  common segment is travelled with the same orientation
               both cycle SignedArea have the same sign

        """
        if abs(self.area)>abs(Cy.area):
            flip,path = self.intersect(Cy)
            if len(path)>0:
                if path[0]>0: # cycles share a segments
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
                else: #cycles share a point 
                    areabig   = abs(self.area)
                    #areasmall = abs(Cy.area)
                    cycomb    = self + Cy
                    areacomb  = abs(cycomb.area)
                    if areacomb < areabig:
                        return True,path
                    else:
                        return False,np.array([])
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

        flip  : boolean
                indicates the orientation of the common brin
                True : flipped direction
                False : same direction
        brin  : common portion of th

        Notes
        -----

        flip,path = C1.intersect(C2)

        if flip == True --> The 2 cycles have a reverse travel
        direction for the common path

        path is the common path along the cycle

        """
        #
        # segment number
        #
        e1  = self.edges
        e2  = Cy.edges
        #
        # point number
        #
        v1  = self.vertices
        v2  = Cy.vertices
        #
        # intersection of segments and vertices
        #
        u   = np.intersect1d(e1,e2)
        v   = np.intersect1d(v1,v2)
        #
        Nis  = len(u)
        Nip  = len(v)
        #print Nis," segments en commun"
        if Nis>0:
            if Nis>1:
                brin1 = self.cycle
                brin2 = Cy.cycle
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

                brin = u
        else:
            flip  = False
            if Nip > 1:
                brin  = np.array([])
            else:
                brin  = v

        return(flip,brin)


    def split(self,cyin):
        """ split(Cy)

          Parameters
          ----------
          cyin :  input cycle

          Returns
          -------
          cyout : list of output cycles
          punctual : boolean
                True if punctual connection

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
        u = np.where(path>0)[0] # test share segment condition
        v = np.where(path<0)[0] # test share point condition
        if (incl) & (len(u)>0):
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

            return(False,cyout)
        if (incl) & (len(v)>0):
            return(True,None)
        else:
            return(False,None)


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


    def show(self,**kwargs):
        """
        show cycle

        Acceleration can be obtained if the polygon is calculated once
        """
        #nx.draw_networkx_edges(self.G,self.G.pos,width=2,edge_color=color,alpha=0.4)
        npoints = filter(lambda x : x <0 ,self.cycle)
        coords  = map(lambda x : self.G.pos[x],npoints)
        poly = geu.Polygon(sh.MultiPoint(tuple(coords)),self.cycle)
        fig,ax = poly.plot(**kwargs)
        return fig,ax
        #plt.draw ()

if __name__ == "__main__":
    doctest.testmod()
