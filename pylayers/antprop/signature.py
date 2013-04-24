#-*- coding:Utf-8 -*-
import doctest
import numpy as np
#import scipy as sp
import scipy.linalg as la
import pdb
import networkx as nx
import pylayers.gis.layout as layout
import pylayers.util.geomutil as geu
#import pylayers.util.graphutil as gph
import pylayers.util.pyutil as pyu
import matplotlib.pyplot as plt
from pylayers.util.project import *
from mpl_toolkits.mplot3d import Axes3D
from pylayers.antprop.rays import Rays
import copy
#from numba import autojit

def showsig(L,s,tx,rx):
    """
    """
    L.display['thin']=True
    fig,ax = L.showGs()
    L.display['thin']=False
    L.display['edlabel']=True
    L.showGs(fig=fig,ax=ax,edlist=s,width=4)
    plt.plot(tx[0],tx[1],'x')
    plt.plot(rx[0],rx[1],'+')
    plt.title(str(s))
    plt.show()
    L.display['edlabel']=False


def gidl(g):
    """ gi without diffraction

    Return
    ------
   """

    edlist=[]
    for n in g.nodes():
        en = eval(n)
        if type(en) == tuple :
            edlist.append(n)
    return(g.subgraph(edlist))


def edgeout(L,g):
    """

    Parameters
    ----------

    L : Layout
    g : Digraph Gi
    """

    for e in g.edges():
        # extract  both interactions
        i0 = eval(e[0])
        i1 = eval(e[1])
        try:
            nstr0 = i0[0]
        except:
            nstr0 = i0


        try:
            nstr1 = i1[0]
            if len(i1)>2:
                typ=2
            else :
                typ=1

        except:
            nstr1 = i1
            typ = 3


        output = []
        if nstr1>0:
            # segment unitary vector

            l1 = L.seguv(np.array([nstr1]))
            p0 = np.array(L.Gs.pos[nstr0])
            p1 = np.array(L.Gs.pos[nstr1])
            v01  = p1-p0
            v01m = np.sqrt(np.dot(v01,v01))
            v01n = v01/v01m
            v10n = -v01n
            # next interaction
            for i2 in nx.neighbors(g,str(i1)):
                i2 = eval(i2)
                if type(i2)==int:
                    nstr2 = i2
                else:
                    nstr2 = i2[0]
                p2 = np.array(L.Gs.pos[nstr2])
                v12 = p2-p1
                v12m = np.sqrt(np.dot(v12,v12))
                v12n = v12/v12m
                d1 = np.dot(v01n,l1)
                d2 = np.dot(l1,v12n)
#                if nstr0==32 and nstr1 == 42  and nstr2 ==50:
#                    pdb.set_trace()
                if d1*d2>=0 and typ == 1:
                    output.append(str(i2))
#                elif d1*d2>=-0.2 and typ ==2:
                elif typ == 2 :
                    if abs(d1) <0.9 and abs(d2) <0.9 :
                        if d1*d2 >= -0.2:
                            output.append(str(i2))
                else:
                    pass
        g.add_edge(str(i0),str(i1),output=output)

    return(g)

class Signatures(dict):
    """
    gathers all signatures from a layout given tx and rx

    Attributes
    ----------
        L : gis.Layout
        pTx : numpy.ndarray
            position of Tx
        pRx : numpy.ndarray
            position of Rx
    """

    def __init__(self,L,source,target):
        """
        Parameters
        ----------
        L : Layout
        source :
        target :
        """
        self.L = L
        self.source = source
        self.target = target

    def __repr__(self):
        size = {}
        s = self.__class__.__name__ + ' : '  + '\n' + '------------------'+'\n'
        s = s + str(self.__sizeof__())+'\n'
        for k in self:
            size[k] = len(self[k])
        s = s + 'from : '+ str(self.source) + ' to ' + str(self.target)+'\n'
        for k in self:
            s = s + str(k) + ' : ' + str(len(self[k])) + '\n'

        return(s)

    def info(self):
        """
        """
        print "Signatures for scenario defined by :"
        print "Layout"
        print "======"
        L = self.L.info()
        print "================================"
        print "source : ", self.source
        print "target : ", self.target



    def propaths(self,G, source, target, cutoff=None):
        """ all_simple_paths

        Parameters
        ----------

        G : networkx Graph Gi
        source : int
        target : int 
        cutoff : int

        Notes
        -----

        adapted from all_simple_path of networkx 


        """
        #print "source :",source
        #print "target :",target

        if cutoff < 1:
            return

        visited = [source]
        # stack is a list of iterators
        stack = [iter(G[source])]
        # while the list of iterators is not void


        while stack: #
            # children is the last iterator of stack

            children = stack[-1]
            # next child
            child = next(children, None)
            #print "child : ",child
            #print "visited :",visited
            if child is None  : # if no more child
                stack.pop()   # remove last iterator
                visited.pop() # remove from visited list
            elif len(visited) < cutoff: # if visited list is not too long
                if child == target:  # if child is the target point
                    #print visited + [target]
                    yield visited + [target] # output signature
                elif child not in visited: # else visit other node
                    stack.append(iter(G[visited[-1]][child]['output']))
                    visited.append(child)

            else: #len(visited) == cutoff (visited list is too long)
                if child == target or target in children:
                    #print visited + [target]
                    yield visited + [target]
                stack.pop()
                visited.pop()


    def calsig(self,G,dia={},cutoff=None):
        """

        G   : Gf graph
        dia : dictionnary of interactions
        cutoff : integer

        """
        if cutoff < 1:
            return

        di=copy.deepcopy(dia)
        source = 'Tx'
        target = 'Rx'
        d={}

        visited = [source]
        stack = [iter(G[source])]

        out=[]

        while stack:
#            pdb.set_trace()
            children = stack[-1]
            child = next(children, None)
            if child is None:
                stack.pop()
                visited.pop()
                if len(out) !=0:
                    out.pop()
                    out.pop()
            elif len(visited) < cutoff:
                if child == target:
                    lot = len(out)
                    try:
                        d.update({lot:d[lot]+(out)})
                    except:
                        d[lot]=[]
                        d.update({lot:d[lot]+(out)})
#                    yield visited + [target]
                elif child not in visited:
                    visited.append(child)
                    out.extend(di[child])
                    stack.append(iter(G[child]))
            else: #len(visited) == cutoff:
                if child == target or target in children:
#                    yield visited + [target]
                    lot = len(out)
                    try:
                        d.update({lot:d[lot]+(out)})
                    except:
                        d[lot]=[]
                        d.update({lot:d[lot]+(out)})
                stack.pop()
                visited.pop()
                if len(out) !=0:
                    out.pop()
                    out.pop()
        return d

#    def all_simple_paths(self,G, source, target, cutoff=None):
#        """ all_simple_paths

#        Parameters
#        ----------

#        G : networkx Graph Gi
#        source : int
#        target : int 
#        cutoff : int

#        Notes
#        -----

#        adapted from all_simple_path of networkx 


#        """
#        if cutoff < 1:
#            return

#        visited = [source]
#        # stack is a list of iterators
#        stack = [iter(G[source])]
#        # while the list of iterators is not void
#        while stack: #
#            # children is the last iterator of stack
#            children = stack[-1]
#            # next child
#            child = next(children, None)
#            if child is None: # if no more child
#                stack.pop()   # remove last iterator
#                visited.pop() # remove from visited list
#            elif len(visited) < cutoff: # if visited list is not too long 
#                if child == target:  # if child is the target point 
#                    yield visited + [target] # output signature
#                elif child not in visited: # else visit other node
#                    visited.append(child)
#                    # explore all child connexion
#                    # TODO : limit the explorable childs
#                    lc = [k for k in iter(G[child])]
#                    n  = len(lc)
#                    if n >12:
#                        explore=iter([lc[k] for k in
#                                      np.unique(np.random.randint(0,n,12))])
#                    else:
#                        explore=iter(G[child])
#                    stack.append(explore)
#                    #stack.append(iter(G[child]))
#            else: #len(visited) == cutoff (visited list is too long)
#                if child == target or target in children:
#                    yield visited + [target]
#                stack.pop()
#                visited.pop()

#    def run(self,metasig,cutoff=1):
#        """ get signatures (in one list of arrays) between tx and rx

#        Parameters
#        ----------

#            cutoff : limit the exploration of all_simple_path

#        Returns
#        -------

#            sigslist = numpy.ndarray

#        """
#        try:
#            self.L.dGi
#        except:
#            self.L.buildGi2()
#        # all the vnodes >0  from the room
#        #
#        #NroomTx = self.L.pt2ro(tx)
#        #NroomRx = self.L.pt2ro(rx)
#        #print NroomTx,NroomRx

#        #if not self.L.Gr.has_node(NroomTx) or not self.L.Gr.has_node(NroomRx):
#        #    raise AttributeError('Tx or Rx is not in Gr')

#        # list of interaction in roomTx 
#        # list of interaction in roomRx
#        #ndt = self.L.Gt.node[self.L.Gr.node[NroomTx]['cycle']]['inter']
#        #ndr = self.L.Gt.node[self.L.Gr.node[NroomRx]['cycle']]['inter']

#        lis = self.L.Gt.node[self.source]['inter']
#        lit = self.L.Gt.node[self.target]['inter']

#        # source
#        #ndt1 = filter(lambda l: len(eval(l))>2,ndt) # Transmission
#        #ndt2 = filter(lambda l: len(eval(l))<3,ndt) # Reflexion

#        lisT = filter(lambda l: len(eval(l))>2,lis) # Transmission
#        lisR = filter(lambda l: len(eval(l))<3,lis) # Reflexion

#        # target
#        # ndr1 = filter(lambda l: len(eval(l))>2,ndr) # Transmission
#        # ndr2 = filter(lambda l: len(eval(l))<3,ndr) # Reflexion

#        litT = filter(lambda l: len(eval(l))>2,lit) # Transmission
#        litR = filter(lambda l: len(eval(l))<3,lit) # Reflexion

#        # tx,rx : attaching rule
#        #
#        # tx attachs to out transmisision point
#        # rx attachs to in transmission point

#        #
#        # WARNING : room number <> cycle number
#        #

#        #ncytx = self.L.Gr.node[NroomTx]['cycle']
#        #ncyrx = self.L.Gr.node[NroomRx]['cycle']

#        #ndt1 = filter(lambda l: eval(l)[2]<>ncytx,ndt1)
#        #ndr1 = filter(lambda l: eval(l)[1]<>ncyrx,ndr1)

#        lisT = filter(lambda l: eval(l)[2]<>self.source,lisT)
#        litT = filter(lambda l: eval(l)[1]<>self.target,litT)

#        #ndt = ndt1 + ndt2
#        #ndr = ndr1 + ndr2
#        lis  = lisT + lisR
#        lit  = litT + litR

#        #ntr = np.intersect1d(ndt, ndr)
#        li = np.intersect1d(lis, lit)

##        for meta in metasig:
#        Gi = nx.DiGraph()
#        for cycle in metasig:
#            Gi = nx.compose(Gi,self.L.dGi[cycle])

#        # facultative update positions
#        Gi.pos = {}
#        for cycle in metasig:
#            Gi.pos.update(self.L.dGi[cycle].pos)
#        pdb.set_trace()
#        #
#        # 
#        #
#        # remove diffractions from Gi
#        Gi = gidl(Gi)
#        # add 2nd order output to edges
#        Gi = edgeout(self.L,Gi)
#        #for interaction source  in list of source interaction 
#        for s in lis:
#            #for target interaction in list of target interaction
#            for t in lit:

#                if (s != t):
#                    #paths = list(nx.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
#                    #paths = list(self.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
#                    paths = list(self.propaths(Gi,source=s,target=t,cutoff=cutoff))
#                    #paths = [nx.shortest_path(Gi,source=s,target=t)]
#                else:
#                    #paths = [[nt]]
#                    paths = [[s]]
#                ### supress the followinfg loops .
#                for path in paths:
#                    sigarr = np.array([],dtype=int).reshape(2, 0)
#                    for interaction in path:

#                        it = eval(interaction)
#                        if type(it) == tuple:
#                            if len(it)==2: #reflexion
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[1]],dtype=int)))
#                            if len(it)==3: #transmission
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[2]],dtype=int)))
#                        elif it < 0: #diffraction
#                            sigarr = np.hstack((sigarr,
#                                                np.array([[it],[3]],dtype=int)))
#                    #print sigarr
#                    try:
#                        self[len(path)] = np.vstack((self[len(path)],sigarr))
#                    except:
#                        self[len(path)] = sigarr


    def run2(self,cutoff=1,dcut=2):
        """ get signatures (in one list of arrays) between tx and rx

        Parameters
        ----------

            cutoff : limit the exploration of all_simple_path

        Returns
        -------

            sigslist = numpy.ndarray

        """
#        try:
#            self.L.dGi
#        except:
#            self.L.buildGi2()

        # all the vnodes >0  from the room
        #
        #NroomTx = self.L.pt2ro(tx)
        #NroomRx = self.L.pt2ro(rx)
        #print NroomTx,NroomRx

        #if not self.L.Gr.has_node(NroomTx) or not self.L.Gr.has_node(NroomRx):
        #    raise AttributeError('Tx or Rx is not in Gr')

        # list of interaction in roomTx 
        # list of interaction in roomRx
        #ndt = self.L.Gt.node[self.L.Gr.node[NroomTx]['cycle']]['inter']
        #ndr = self.L.Gt.node[self.L.Gr.node[NroomRx]['cycle']]['inter']


#        lisT = filter(lambda l: len(eval(l))>2,lis) # Transmission
#        lisR = filter(lambda l: len(eval(l))<3,lis) # Reflexion

#        # target
#        # ndr1 = filter(lambda l: len(eval(l))>2,ndr) # Transmission
#        # ndr2 = filter(lambda l: len(eval(l))<3,ndr) # Reflexion

#        litT = filter(lambda l: len(eval(l))>2,lit) # Transmission
#        litR = filter(lambda l: len(eval(l))<3,lit) # Reflexion

#        # tx,rx : attaching rule
#        #
#        # tx attachs to out transmisision point
#        # rx attachs to in transmission point

#        #
#        # WARNING : room number <> cycle number
#        #

#        #ncytx = self.L.Gr.node[NroomTx]['cycle']
#        #ncyrx = self.L.Gr.node[NroomRx]['cycle']

#        #ndt1 = filter(lambda l: eval(l)[2]<>ncytx,ndt1)
#        #ndr1 = filter(lambda l: eval(l)[1]<>ncyrx,ndr1)

#        lisT = filter(lambda l: eval(l)[2]<>self.source,lisT)
#        litT = filter(lambda l: eval(l)[1]<>self.target,litT)

#        #ndt = ndt1 + ndt2
#        #ndr = ndr1 + ndr2
#        lis  = lisT + lisR
#        lit  = litT + litR

#        #ntr = np.intersect1d(ndt, ndr)
#        li = np.intersect1d(lis, lit)


#        for meta in metasig:
#            Gi = nx.DiGraph()
#            for cycle in meta:
#                Gi = nx.compose(Gi,self.L.dGi[cycle])
#            # facultative update positions
#            Gi.pos = {}
#            for cycle in meta:
#                Gi.pos.update(self.L.dGi[cycle].pos)
#            #
#            #
#            #
#        metasig=self.lineofcycle(cs,ct)

############################################################
##       obtain the list of cycle in line

        cs = self.source
        ct = self.target
        lcil=self.L.cycleinline(cs,ct)
        lca = [] # list of cycle around
        for cy in lcil:
            ncy = nx.neighbors(self.L.Gt,cy)
            lca = lca+ncy
        lca = list(np.unique(np.array(lca)))
        lca = lcil

############################################################
##       Compose graph of interactions with the list lca of
##       cycles around line of cycles

        Gi = nx.DiGraph()
        for cycle in lca:
            Gi = nx.compose(Gi,self.L.dGi[cycle])

        # facultative update positions
        Gi.pos = {}
        for cycle in lca:
            Gi.pos.update(self.L.dGi[cycle].pos)

#        #
#        #
#        #
        Gf = nx.DiGraph()
        Gf.pos = {}
        # remove diffractions from Gi
        Gi = gidl(Gi)
        # add 2nd order output to edges
        Gi = edgeout(self.L,Gi)
        #for interaction source  in list of source interaction

####################################################
#        filter list of interactions in termination cycles

        lis = self.L.Gt.node[lcil[0]]['inter']
        lit = self.L.Gt.node[lcil[-1]]['inter']
        # filter lis remove transmission coming from outside
        lli = []
        for li in lis:
            ei = eval(li)
            if len(ei)==2:
                lli.append(li)
            if len(ei)==3:
                if ei[2]<>cs:
                   lli.append(li)
        # filter lit remove transmission going outside
        llt = []
        for li in lit:
            ei = eval(li)
            if len(ei)==2:
                llt.append(li)
            if len(ei)==3:
                if ei[2]==ct:
                   llt.append(li)
        lis = lli
        lit = llt


#################################################
#       propaths (a.k.a. all simple path) per adjacent cycles along cycles in line
#       Obtaining Gf: filtred graph of Gi with Gc ( rename Gt in Gc )

        for ic in np.arange(len(lcil)-2):
            lsource = []
            ltarget = []
            linter = self.L.Gt.node[lcil[ic]]['inter']
            # determine list of sources
            if ic>0:
                ls = self.L.Gt[lcil[ic]][lcil[ic+1]]['segment']
                for source in ls:
                    lsource.append(str((source, lcil[ic], lcil[ic+1])))
            else:
                lsource = lis

            # determine list of targets
            if ic+2 < len(lcil)-1:
            #if ic+3 < len(lcil)-1:
                lt = self.L.Gt[lcil[ic+1]][lcil[ic+2]]['segment']
                #lt = self.L.Gt[lcil[ic+2]][lcil[ic+3]]['segment']
                for target in lt:
                    ltarget.append(str((target , lcil[ic+1], lcil[ic+2])))
                    #ltarget.append(str((target , lcil[ic+2], lcil[ic+3])))
            else:
                ltarget = lit

            lt   = filter(lambda l: len(eval(l))==3,linter)
            #lti = filter(lambda l: eval(l)[2]==lcil[ic+1],lt)
            lto = filter(lambda l: eval(l)[2]<>lcil[ic],lt)
            ltom = filter(lambda l: eval(l)[2]!=lcil[ic-1],lto)
            ltomp = filter(lambda l: eval(l)[2]!=lcil[ic+1],ltom)

            lsource = lsource + ltomp
            #pdb.set_trace()
            for s in lsource :
                #print s
                for t in ltarget:
                    #print t
                    paths = list(self.propaths(Gi,source=s,target=t,cutoff=cutoff))

                    for path in paths:
                        itm1 = path[0]
                        if itm1 not in Gf.node.keys():
                            Gf.add_node(itm1)
                            Gf.pos[itm1]=self.L.Gi.pos[itm1]
                        for it in path[1:]:
                            if it not in Gf.node.keys():
                                Gf.add_node(it)
                                Gf.pos[it]=self.L.Gi.pos[it]
                            Gf.add_edge(itm1,it)
                            itm1 = it
#                        else:
#                            #paths = [[nt]]
#                            paths = [[s]]


################################################################
#       Obtain position of centroid of cycles source and target


        poly1 = self.L.Gt.node[cs]['polyg']
        cp1 = poly1.centroid.xy

        poly2 = self.L.Gt.node[ct]['polyg']
        cp2 = poly2.centroid.xy
        pcs = np.array([cp1[0][0],cp1[1][0]])
        pct = np.array([cp2[0][0],cp2[1][0]])

        Gf.add_node('Tx')
        Gf.pos['Tx']=tuple(pcs[:2])

        for i in self.L.Gt.node[cs]['inter']:
            if i in  Gf.nodes():
                Gf.add_edge('Tx',i)

        Gf.add_node('Rx')
        Gf.pos['Rx']=tuple(pct[:2])

        for i in self.L.Gt.node[ct]['inter']:
            if i in  Gf.nodes():
                Gf.add_edge(i,'Rx')
        # a =[ 0,  1,  2,  1,  4,  1,  6,  1,  8,  1, 10, 1]
        # aa = np.array(a)
        # X=aa.reshape((2,3,2)) # r x i x 2
        # Y=X.swapaxes(0,2) # 2 x i x r



        self.Gf = Gf
        print 'signatures'
        co = nx.dijkstra_path_length(Gf,'Tx','Rx')
        sig = self.calsig(Gf,dia=self.L.di,cutoff=co+dcut)


        for k in sig:
            ns = len(sig[k])
            nbi = k/2
            nr = ns/k
            self[nbi]=(np.array(sig[k]).reshape(nr,nbi,2)).swapaxes(0,2)


        d={}

        for k in self :
            a= self[k]
            nbr = np.shape((a[0]))[1]
            d[k]=np.zeros((2*nbr,k),dtype=int)
            for r in range(nbr):
                for i in range(k):
                    d[k][2*r,i]=a[0,i,r]
                    d[k][2*r+1,i]=a[1,i,r]
        self.update(d)

        # self [nbi] = 2 x i x r


                ### supress the following loops .
#                for path in paths:
#                    sigarr = np.array([],dtype=int).reshape(2, 0)
#                    for interaction in path:
#
#                        it = eval(interaction)
#                        if type(it) == tuple:
#                            if len(it)==2: #reflexion
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[1]],dtype=int)))
#                            if len(it)==3: #transmission
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[2]],dtype=int)))
#                        elif it < 0: #diffraction
#                            sigarr = np.hstack((sigarr,
#                                                np.array([[it],[3]],dtype=int)))
#                    #print sigarr
#                    try:
#                        self[len(path)] = np.vstack((self[len(path)],sigarr))
#                    except:
#                        self[len(path)] = sigarr


#        for s in lis:
#            #for target interaction in list of target interaction
#            for t in lit:

#                if (s != t):
#                    #paths = list(nx.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
#                    #paths = list(self.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
#                    paths = list(self.propaths(Gi,source=s,target=t,cutoff=cutoff))
#                    #paths = [nx.shortest_path(Gi,source=s,target=t)]
#                else:
#                    #paths = [[nt]]
#                    paths = [[s]]
#                ### supress the followinfg loops .
#                for path in paths:
#                    sigarr = np.array([],dtype=int).reshape(2, 0)
#                    for interaction in path:

#                        it = eval(interaction)
#                        if type(it) == tuple:
#                            if len(it)==2: #reflexion
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[1]],dtype=int)))
#                            if len(it)==3: #transmission
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[2]],dtype=int)))
#                        elif it < 0: #diffraction
#                            sigarr = np.hstack((sigarr,
#                                                np.array([[it],[3]],dtype=int)))
#                    #print sigarr
#                    try:
#                        self[len(path)] = np.vstack((self[len(path)],sigarr))
#                    except:
#                        self[len(path)] = sigarr

    def run3(self,cs,ct,cutoff=1):

        ns = nx.neighbors(self.L.Gt,cs)
        nt = nx.neighbors(self.L.Gt,ct)
        print ns
        print nt
        path=[]
        for s in ns:
            for t in nt:
                p=nx.dijkstra_path(self.L.Gt,s,t)
                if not cs in p and not ct in p:
                    path.append(p)
        return path


    def meta(self):
        G = self.L.Gt
        # metasig = list(nx.all_simple_paths(self.L.Gt,source=self.source,target=self.target,cutoff=cutoff))
        #for cutoff in range(1,5):
        metasig = nx.shortest_path(G,source=self.source,target=self.target)
        for c in metasig:
            try :
                n = np.hstack((n,np.array(G.neighbors(c))))
            except:
                n = np.array(G.neighbors(c))
        n = np.unique(n)
        for d in n:
            try :
                n = np.hstack((n,np.array(G.neighbors(d))))
            except:
                n = np.array(G.neighbors(d))

        return np.unique(n)


    def lineofcycle(self,cs=[],ct=[]):
        """ shortest path between 2 cycle
        """
        if cs == []:
            cs = self.source
        if ct == []:
            ct = self.target
        return nx.shortest_path(self.L.Gt,source=cs,target=ct)



    def showi(self,uni=0,us=0):
        """ interactive show

        press n to visit signatures sequentially

        Parameters
        ----------
            uni : index of interaction dictionnary keys
            us : signature index

        """
        plt.ion()
        fig=plt.figure()
        ax=fig.add_subplot(111)
#        fig,ax=self.L.showG(fig=fig,ax=ax,graph='s')
#        plt.draw()

        nit = self.keys()
        ni = nit[uni]
        ust = len(self[ni])/2

        poly1 = self.L.Gt.node[self.source]['polyg']
        cp1 = poly1.centroid.xy

        poly2 = self.L.Gt.node[self.target]['polyg']
        cp2 = poly2.centroid.xy

        ptx = np.array([cp1[0][0],cp1[1][0]])
        prx = np.array([cp2[0][0],cp2[1][0]])

        st='a'

        while st != 'q':
            inter=[]
            ax=fig.add_subplot(111)
            fig,ax=self.L.showG(fig=fig,ax=ax,graph='s')
            title = '# interaction :', ni, 'signature #',us,'/',ust
            ax.set_title(title)

            line = ptx
            # draw terminal points (centroid of source and target cycle)

            ax.plot(ptx[0],prx[1],'xr')
            ax.plot(prx[0],prx[1],'xb')

            if ni not in self.keys():
                print "incorrect number of interactions"
            pos={}

            try:
                for u in self[ni][us*2]:
                    pos.update({u:self.L.Gs.pos[u]})
                    line = np.vstack((line,np.array((self.L.Gs.pos[u]))))
                nx.draw_networkx_nodes(self.L.Gs,pos=pos,nodelist=pos.keys(),node_color='r',ax=ax)

                for ii in self[ni][(us*2)+1]:
                    if ii == 1:
                        inter.append('R')
                    if ii == 2:
                        inter.append('T')
            except:
                print "signature index out of bounds of signature"

            line = np.vstack((line,prx))
            ax.plot(line[:,0],line[:,1])
            plt.draw()
            print inter
            st = raw_input()
            ax.cla()
            if st == 'n':
                if us+2 <= ust:
                    us=us+2

                else:
                    uni = uni+1
                    try:
                        ni = nit[uni]
                        ust = len(self[ni])/2
                        us=0
                    except:
                        uni=0
                        ni=nit[uni]
                        us = 0


            else:
                print 'press n for next signature'


    def rays(self,ptx,prx):
        """ from signatures dict to 2D rays

        Parameters
        ----------

            dsig : dict

        Returns
        -------

            rays : dict

        """
        rays = Rays(ptx,prx)
        for k in self:
            tsig = self[k]
            shsig = np.shape(tsig)
            for l in range(shsig[0]/2):
                sig = tsig[2*l:2*l+2,:]
                s   = Signature(sig)
                Yi  = s.sig2ray(self.L, ptx[:2], prx[:2])
                if Yi is not None:
                    Yi = np.fliplr(Yi)
                    nint = len(sig[0, :])
                    if nint in rays.keys():
                        Yi3d = np.vstack((Yi[:, 1:-1], np.zeros((1, nint))))
                        Yi3d = Yi3d.reshape(3, nint, 1)
                        rays[nint]['pt'] = np.dstack(( rays[nint]['pt'], Yi3d))
                        rays[nint]['sig'] = np.dstack(( rays[nint]['sig'], sig.reshape(2, nint, 1)))
                    else:
                        rays[nint] = {'pt': np.zeros((3, nint, 1)),
                                      'sig': np.zeros((2, nint,
                                                            1),dtype=int)}
                        rays[nint]['pt'][0:2, :, 0] = Yi[:, 1:-1]
                        rays[nint]['sig'][:, :, 0] = sig
        return rays

class Signature(object):
    """ class Signature

    A signature contains two lists

    seq : list of interaction numbers
    typ : list of interaction type
    """
    def __init__(self, sig):
        """
        pa  : tail point of interaction segment
        pb  : head point of interaction segment
        pc  : center point of interaction segment
        typ : type of interaction 1-R 2-T 3-D
        seq : sequence of interaction point (edges (>0)  or vertices (<0)
        """
        self.seq = sig[0, :]
        self.typ = sig[1, :]

#    def __repr__(self):
#        s = self.__class__ + ':' + str(self.__sizeof__())+'\n'
#        s = s + self.seq + '\n' + self.typ
#        return s

    def info(self):
        """
        """
        for k in self.__dict__.keys():
            print k, ':', self.__dict__[k]

    def split(self):
        """
        split signature
        """
        pass

    def ev(self, L):
        """  evaluation of Signature

        Parameters
        ----------
            L : Layout

        Notes
        -----
        Le type des interactions d extremite reste indetermine a ce stade
        """
        N = len(self.seq)
        self.pa = np.zeros((2, N))  # tail
        self.pb = np.zeros((2, N))  # head
        self.pc = np.zeros((2, N))  # center
        #self.typ = np.zeros(N)
        self.norm = np.zeros((2, N))

        for n in range(N):
            k = self.seq[n]
            if k > 0:  # segment
                ta, he = L.Gs.neighbors(k)
                norm1 = L.Gs.node[k]['norm']
                norm = np.array([norm1[0], norm1[1]])
                self.pa[:, n] = np.array(L.Gs.pos[ta])
                self.pb[:, n] = np.array(L.Gs.pos[he])
                self.pc[:, n] = np.array(L.Gs.pos[k])
                self.norm[:, n] = norm
                #self.typ[n] = 1
            else:      # node
                pa = np.array(L.Gs.pos[k])
                norm = np.array([0, 0])
                self.pa[:, n] = pa
                self.pb[:, n] = pa
                self.pc[:, n] = pa
                self.norm[:, n] = norm
                #self.typ[n] = 3
        #
        #  vecteurs entre deux points adjascents de la signature
        #
        #self.v   = s.pc[:,1:]-s.pc[:,:-1]
        #self.vn  = self.v /np.sqrt(np.sum(self.v*self.v,axis=0))
        #u1       = np.sum(self.norm*self.vn[:,0:-1],axis=0)
        #u2       = np.sum(self.norm*self.vn[:,1:],axis=0)
        #self.typ = sign(u1*u2)
        #return(vn)
        #return(typ)

    def evtx(self, L, tx, rx):
        """ evtx

        Parameters
        ----------
            L  : Layout
            tx : np.array (2xN)
            rx : np.array (2xM)

        """
        self.pa = tx.reshape(2, 1)
        self.pb = tx.reshape(2, 1)
        self.pc = tx.reshape(2, 1)
        self.typ = np.array([0])
        for k in self.seq:
            if k > 0:
                ta, he = L.Gs.neighbors(k)
                norm1 = L.Gs.node[k]['norm']
                norm = np.array([norm1[0], norm1[1]]).reshape(2, 1)
                pa = np.array(L.Gs.pos[ta]).reshape(2, 1)
                pb = np.array(L.Gs.pos[he]).reshape(2, 1)
                pc = np.array(L.Gs.pos[k]).reshape(2, 1)
                self.pa = np.hstack((self.pa, pa))
                self.pb = np.hstack((self.pb, pb))
                self.pc = np.hstack((self.pc, pc))
                try:
                    self.norm = np.hstack((self.norm, norm))
                except:
                    self.norm = norm
                self.typ = np.hstack((self.typ, np.array([1])))
            else:
                pa = np.array(L.Gs.pos[k]).reshape(2, 1)
                norm = np.array([0, 0]).reshape(2, 1)
                self.pa = np.hstack((self.pa, pa))
                self.pb = np.hstack((self.pb, pa))
                self.pc = np.hstack((self.pc, pa))
                try:
                    self.norm = np.hstack((self.norm, norm))
                except:
                    self.norm = norm
                self.typ = np.hstack((self.typ, np.array([3])))
        self.pa = np.hstack((self.pa, rx.reshape(2, 1)))
        self.pb = np.hstack((self.pb, rx.reshape(2, 1)))
        self.pc = np.hstack((self.pc, rx.reshape(2, 1)))
        self.typ = np.hstack((self.typ, np.array([0])))
        #
        #  vecteur entre deux points adjascents de la signature
        #
        self.v = s.pc[:, 1:] - s.pc[:, :-1]
        self.vn = self.v / np.sqrt(sum(self.v * self.v, axis=0))
        u1 = sum(self.norm * self.vn[:, 0:-1], axis=0)
        u2 = sum(self.norm * self.vn[:, 1:], axis=0)
        self.typ = np.sign(u1 * u2)
        #return(vn)
        #return(typ)


    def image(self, tx):
        """
        Compute the images of tx with respect to the signature segments
        Parameters
        ----------
            tx : numpy.ndarray
        Returns
        -------
            M : numpy.ndarray
        """

        pa = self.pa
        pb = self.pb
        pab = pb - pa
        alpha = np.sum(pab * pab, axis=0)
        zalpha = np.where(alpha == 0.)
        alpha[zalpha] = 1.

        a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
        b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
        c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
                            pa[1, :] * (pa[1, :] - pb[1, :]) *
                            (pb[0, :] - pa[0, :]))
        d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
                            pa[0, :] * (pa[1, :] - pb[1, :]) *
                            (pb[0, :] - pa[0, :]))

        typ = self.typ
        # number of interactions
        N = np.shape(pa)[1]

        S = np.zeros((N, 2, 2))
        S[:, 0, 0] = -a
        S[:, 0, 1] = b
        S[:, 1, 0] = b
        S[:, 1, 1] = a
        blocks = np.zeros((N - 1, 2, 2))
        A = np.eye(N * 2)

        # detect diffraction
        usig = np.nonzero(typ[1:] == 3)[0]
        if len(usig) > 0:
            blocks[usig, :, :] = np.zeros((2, 2))
        # detect transmission
        tsig = np.nonzero(typ[1:] == 2)[0]
        if len(tsig) > 0:
            #blocks[tsig, :, :] = np.zeros((2, 2))
            blocks[tsig, :, :] = -np.eye(2)
        # detect reflexion
        rsig = np.nonzero(typ[1:] == 1)[0]
        if len(rsig) > 0:
            blocks[rsig, :, :] = S[rsig + 1, :, :]
        A = pyu.fill_block_diag(A, blocks, 2, -1)

        y = np.zeros(2 * N)
        if typ[0] == 1:
            vc0 = np.array([c[0], d[0]])
            v0 = np.dot(-S[0, :, :], tx) + vc0
        if typ[0] == 2:
            v0 = tx
        if typ[0] == 3:
            v0 = pa[:, 0]
        y[0:2] = v0
        for i in range(len(typ[1:])):
            if typ[i + 1] == 1:
                y[2 * (i + 1):2 * (i + 1) + 2] = np.array([c[i + 1], d[i + 1]])
            if typ[i + 1] == 2:
                #y[2 * (i + 1):2 * (i + 1) + 2] = y[2*i:2*i+2]
                y[2 * (i + 1):2 * (i + 1) + 2] = np.array([0,0]) 
            if typ[i + 1] == 3:
                y[2 * (i + 1):2 * (i + 1) + 2] = pa[:, i + 1]

        x = la.solve(A, y)
        M = np.vstack((x[0::2], x[1::2]))
        return M
    def backtrace(self, tx, rx, M):
        """
        backtracing step: given the image, tx, and rx, this function
        traces the 2D ray.

        Parameters
        ----------
            tx :  numpy.ndarray
                  transmitter
            rx :  numpy.ndarray
                  receiver
            M  :  numpy.ndarray
                  images obtained using image()

        Returns
        -------
            Y : numpy.ndarray
                2D ray

        Examples
        --------

        .. plot::
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> from pylayers.gis.layout import *
            >>> from pylayers.antprop.signature import *
            >>> L = Layout()
            >>> L.buildGt()
            >>> L.buildGr()
            >>> seq = [1,5,1]
            >>> s = Signature(seq)
            >>> tx = np.array([4,-1])
            >>> rx = np.array([1,1])
            >>> s.ev(L)
            >>> M = s.image(tx)
            >>> Y = s.backtrace(tx,rx,M)
            >>> fig = plt.figure()
            >>> ax = fig.add_subplot(111)
            >>> l1 = ax.plot(tx[0],tx[1],'or')
            >>> l2 = ax.plot(rx[0],rx[1],'og')
            >>> l3 = ax.plot(M[0,:],M[1,:],'ob')
            >>> l4 = ax.plot(Y[0,:],Y[1,:],'xk')
            >>> ray = np.hstack((np.hstack((tx.reshape(2,1),Y)),rx.reshape(2,1)))
            >>> l5 = ax.plot(ray[0,:],ray[1,:],color='#999999',alpha=0.6,linewidth=0.6)
            >>> fig,ax = L.showGs(fig,ax)
            >>> plt.show()

        """

        pa = self.pa
        pb = self.pb
        typ = self.typ

        N = np.shape(pa)[1]
        I2 = np.eye(2)
        z0 = np.zeros((2, 1))

        pkm1 = rx.reshape(2, 1)
        Y = pkm1
        k = 0
        beta = .5
        cpt = 0
        while (((beta <= 1) & (beta >= 0)) & (k < N)):
            if int(typ[k]) != 3:
                # Formula (30) of paper Eucap 2012
                l0 = np.hstack((I2, pkm1 - M[:, N - (k + 1)].reshape(2, 1), z0
                                ))
                l1 = np.hstack((I2, z0,
                                pa[:, N - (k + 1)].reshape(2, 1) -
                                pb[:, N - (k + 1)].reshape(2, 1)
                                ))

                T = np.vstack((l0, l1))
                yk = np.hstack((pkm1[:, 0].T, pa[:, N - (k + 1)].T))
                deT = np.linalg.det(T)
                if abs(deT) < 1e-15:
                    return(None)
                xk = la.solve(T, yk)
                pkm1 = xk[0:2].reshape(2, 1)
                gk = xk[2::]
                alpha = gk[0]
                beta = gk[1]
                Y = np.hstack((Y, pkm1))
            else:
                Y = np.hstack((Y, pa[:, k].reshape((2, 1))))
                pkm1 = pa[:, k].reshape((2, 1))
            k = k + 1
        if ((k == N) & ((beta > 0) & (beta < 1)) & ((alpha > 0) & (alpha < 1))):
            Y = np.hstack((Y, tx.reshape(2, 1)))
            return(Y)
        else:
            return(None)

    def sig2ray(self, L, pTx, pRx):
        """
        convert a signature to a 2D ray
        Parameters
        ----------
            L : Layout
            pTx : ndarray
                2D transmitter position
            pRx : ndarray
                2D receiver position
        Returns
        -------
            Y : numpy.ndarray
        """
        try:
            L.Gr
        except:
            L.build()

        self.ev(L)
        M = self.image(pTx)
        Y = self.backtrace(pTx, pRx, M)
        return Y
 
# def get_sigslist(self, tx, rx):
#        """
#        get signatures (in one list of arrays) between tx and rx
#        Parameters
#        ----------
#            tx : numpy.ndarray
#            rx : numpy.ndarray
#        Returns
#        -------
#            sigslist = numpy.ndarray
#        """
#        try:
#            self.L.Gi
#        except:
#            self.L.build()
#        # all the vnodes >0  from the room
#        #
#        NroomTx = self.L.pt2ro(tx)
#        NroomRx = self.L.pt2ro(rx)
#        print NroomTx,NroomRx
#
#        if not self.L.Gr.has_node(NroomTx) or not self.L.Gr.has_node(NroomRx):
#            raise AttributeError('Tx or Rx is not in Gr')
#
#        #list of interaction 
#        ndt = self.L.Gt.node[self.L.Gr.node[NroomTx]['cycle']]['inter']
#        ndr = self.L.Gt.node[self.L.Gr.node[NroomRx]['cycle']]['inter']
#
#        ndt1 = filter(lambda l: len(eval(l))>2,ndt)
#        ndt2 = filter(lambda l: len(eval(l))<3,ndt)
#        ndr1 = filter(lambda l: len(eval(l))>2,ndr)
#        ndr2 = filter(lambda l: len(eval(l))<3,ndr)
#
#        print ndt1
#        print ndr1
#        ndt1 = filter(lambda l: eval(l)[2]<>NroomTx,ndt1)
#        ndr1 = filter(lambda l: eval(l)[1]<>NroomRx,ndr1)
#
#        ndt = ndt1 + ndt2
#        ndr = ndr1 + ndr2
#
#        ntr = np.intersect1d(ndt, ndr)
#        sigslist = []
#
#        for nt in ndt:
#            print nt
#            for nr in ndr:
#                addpath = False
#                print nr
#                if (nt != nr):
#                    try:
#                        path = nx.dijkstra_path(self.L.Gi, nt, nr)
#                        #paths = nx.all_simple_paths(self.L.Gi,source=nt,target=nr)
#                        addpath = True
#                        showsig(self.L,path,tx,rx)
#                    except:
#                        pass
#                if addpath:
#                    sigarr = np.array([]).reshape(2, 0)
#                    for interaction in path:
#                        it = eval(interaction)
#                        if type(it) == tuple:
#                            if len(it)==2: #reflexion
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]],[1]])))
#                            if len(it)==3: #transmission
#                                sigarr = np.hstack((sigarr,
#                                                np.array([[it[0]], [2]])))
#                        elif it < 0: #diffraction
#                            sigarr = np.hstack((sigarr,
#                                                np.array([[it], [3]])))
#                    sigslist.append(sigarr)
#
#        return sigslist
#
#    def update_sigslist(self):
#        """
#        get signatures taking into account reverberations
#
#        Returns
#        -------
#            sigslist: numpy.ndarry
#
#        Notes
#        -----
#        This is a preliminary function need more investigations
#
#        """
#        pTx = self.pTx
#        pRx = self.pRx
#        NroomTx = self.L.pt2ro(pTx)
#        NroomRx = self.L.pt2ro(pRx)
#        if NroomTx == NroomRx:
#            sigslist = self.get_sigslist(pTx, pRx)
#        else:
#            sigslist = []
#            sigtx = self.get_sigslist(pTx, pTx)
#            sigrx = self.get_sigslist(pRx, pRx)
#            sigtxrx = self.get_sigslist(pTx, pRx)
#            sigslist = sigslist + sigtxrx
#            for sigtr in sigtxrx:
#                for sigt in sigtx:
#                    if (sigt[:, -1] == sigtr[:, 0]).all():
#                        if np.shape(sigtr)[1] == 1 or np.shape(sigt)[1] == 1:
#                            pass
#                        else:
#                            sigslist.append(np.hstack((sigt, sigtr[:, 1:])))
#                for sigr in sigrx:
#                    if (sigr[:, 0] == sigtr[:, -1]).all():
#                        if np.shape(sigtr)[1] == 1 or np.shape(sigr)[1] == 1:
#                            pass
#                        else:
#                            sigslist.append(np.hstack((sigtr, sigr[:, 1:])))
#
#        return sigslist
#
#    def image_ceilfloor(self, tx, pa, pb):
#        """
#        Compute the images of tx with respect to ceil or floor
#        Parameters
#        ----------
#            tx : numpy.ndarray
#            pa : numpy.ndarray
#            pb : numpy.ndarray
#        Returns
#        -------
#            M : numpy.ndarray
#        """
#
#        pab = pb - pa
#        alpha = np.sum(pab * pab, axis=0)
#        zalpha = np.where(alpha == 0.)
#        alpha[zalpha] = 1.
#
#        a = 1 - (2. / alpha) * (pa[1, :] - pb[1, :]) ** 2
#        b = (2. / alpha) * (pb[0, :] - pa[0, :]) * (pa[1, :] - pb[1, :])
#        c = (2. / alpha) * (pa[0, :] * (pa[1, :] - pb[1, :]) ** 2 +
#                            pa[1, :] * (pa[1, :] - pb[1, :]) *
#                            (pb[0, :] - pa[0, :]))
#        d = (2. / alpha) * (pa[1, :] * (pb[0, :] - pa[0, :]) ** 2 +
#                            pa[0, :] * (pa[1, :] - pb[1, :]) *
#                            (pb[0, :] - pa[0, :]))
#
#        S = np.zeros((1, 2, 2))
#        S[:, 0, 0] = -a
#        S[:, 0, 1] = b
#        S[:, 1, 0] = b
#        S[:, 1, 1] = a
#        A = np.eye(2)
#
#        vc0 = np.array([c[0], d[0]])
#        y = np.dot(-S[0, :, :], tx) + vc0
#
#        x = la.solve(A, y)
#        M = np.vstack((x[0::2], x[1::2]))
#        return M
#
#    def backtrace_ceilfloor(self, tx, rx, pa, pb, M):
#        """
#        backtracing step: given the image, tx, and rx, this function
#        traces the 2D ray.
#
#        Parameters
#        ----------
#            tx :  numpy.ndarray
#                  transmitter
#            rx :  numpy.ndarray
#                  receiver
#            M  :  numpy.ndarray
#                  images obtained using image()
#
#        Returns
#        -------
#            Y : numpy.ndarray
#                2D ray
#
#
#        """
#        N = np.shape(pa)[1]
#        I2 = np.eye(2)
#        z0 = np.zeros((2, 1))
#
#        pkm1 = rx.reshape(2, 1)
#        Y = pkm1
#        k = 0
#        beta = .5
#        cpt = 0
#        while (((beta <= 1) & (beta >= 0)) & (k < N)):
#            l0 = np.hstack((I2, pkm1 - M[:, N - (k + 1)].reshape(2, 1), z0
#                            ))
#            l1 = np.hstack((I2, z0,
#                            pa[:, N - (k + 1)].reshape(2, 1) -
#                            pb[:, N - (k + 1)].reshape(2, 1)
#                            ))
#
#            T = np.vstack((l0, l1))
#            yk = np.hstack((pkm1[:, 0].T, pa[:, N - (k + 1)].T))
#            deT = np.linalg.det(T)
#            if abs(deT) < 1e-15:
#                return(None)
#            xk = la.solve(T, yk)
#            pkm1 = xk[0:2].reshape(2, 1)
#            gk = xk[2::]
#            alpha = gk[0]
#            beta = gk[1]
#            Y = np.hstack((Y, pkm1))
#            k += 1
#        if ((k == N) & ((beta > 0) & (beta < 1))):  # & ((alpha > 0) & (alpha < 1))):
#            Y = np.hstack((Y, tx.reshape(2, 1)))
#            return(Y)
#        else:
#            return(None)
#   def sigs2rays(self, sigslist):
#        """ from signatures list to 2D rays
#
#        Parameters
#        ----------
#
#            sigslist : list
#
#        Returns
#        -------
#
#            rays : dict
#
#        """
#        rays = {}
#        for sig in sigslist:
#            s = Signature(sig)
#            Yi = s.sig2ray(self.L, self.pTx[:2], self.pRx[:2])
#            if Yi is not None:
#                #pdb.set_trace()
#                Yi = np.fliplr(Yi)
#                nint = len(sig[0, :])
#                if str(nint) in rays.keys():
#                    Yi3d = np.vstack((Yi[:, 1:-1], np.zeros((1, nint))))
#                    Yi3d = Yi3d.reshape(3, nint, 1)
#                    rays[str(nint)]['pt'] = np.dstack((
#                                                      rays[str(nint)]['pt'], Yi3d))
#                    rays[str(nint)]['sig'] = np.dstack((
#                                                       rays[str(nint)]['sig'],
#                                                       sig.reshape(2, nint, 1)))
#                else:
#                    rays[str(nint)] = {'pt': np.zeros((3, nint, 1)),
#                                       'sig': np.zeros((2, nint, 1))}
#                    rays[str(nint)]['pt'][0:2, :, 0] = Yi[:, 1:-1]
#                    rays[str(nint)]['sig'][:, :, 0] = sig
#        return rays


if __name__ == "__main__":
    doctest.testmod()
