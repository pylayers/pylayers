#-*- coding:Utf-8 -*-
"""
Module : signature

functions

showsig

"""
import doctest
import numpy as np
#import scipy as sp
import scipy.linalg as la
import pdb
import networkx as nx
import pylayers.gis.layout as layout
import pylayers.util.geomutil as geu
import pylayers.util.cone as cone
#import pylayers.util.graphutil as gph
import pylayers.util.pyutil as pyu
import pylayers.util.plotutil as plu
import matplotlib.pyplot as plt
from pylayers.util.project import *
from mpl_toolkits.mplot3d import Axes3D
from pylayers.antprop.rays import Rays
import copy
import pickle
import logging
import time
#from numba import autojit

def showsig(L,s,tx=[],rx=[]):
    """
    
    Parameters
    ----------

    L  : Layout 
    s  : 
    tx :
    rx : 

    """
    L.display['thin']=True
    fig,ax = L.showGs()
    L.display['thin']=False
    L.display['edlabel']=True
    L.showGs(fig=fig,ax=ax,edlist=s,width=4)
    if tx !=[]:
        plt.plot(tx[0],tx[1],'x')
    if rx !=[]:
        plt.plot(rx[0],rx[1],'+')
    plt.title(str(s))
    plt.show()
    L.display['edlabel']=False


def gidl(g):
    """ gi without diffraction

   Returns
   -------

   gr 

   """

    edlist=[]
    pos={}
    for n in g.nodes():
        en = eval(n)
        if type(en) == tuple :
            edlist.append(n)
    gr=g.subgraph(edlist)
    dpos = {k:g.pos[k] for k in edlist}
    gr.pos=dpos
    return(gr)


def frontline(L,nc,v):
    """ determine cycle frontline

    This function calculates the scalar product of the normals of a cycle 
    and returns the indev of segments whith are facing the given direction v.
    scalar product < 0.

    Parameters
    ----------

    L : Layout
    nc : cycle number
    v : direction vector

    Returns
    -------

    nsegf : list 

    Example
    -------

    >>> from pylayers.gis.layout import * 
    >>> L = Layout()
    >>> L.build()
    >>> v = np.array([1,1])
    >>> frontline(L,0,v)
    [3, 4]

    See Also
    --------

    run3

    """
    npt = filter(lambda x: x<0, L.Gt.node[nc]['cycle'].cycle)  # points 
    nseg = filter(lambda x: x>0, L.Gt.node[nc]['cycle'].cycle) # segments
    pt  = map(lambda npt : [L.Gs.pos[npt][0],L.Gs.pos[npt][1]],npt)
    pt1 = np.array(pt)   # convert in ndarray
    n1 = geu.Lr2n(pt1.T) # get the normals of the cycle
    ps = np.sum(n1*v[:,np.newaxis],axis=0) # scalar product with vector v
    u = np.where(ps<0)[0]   # keep segment if scalar product <0
    nsegf = map(lambda n: nseg[n],u)
    return nsegf


def edgeout2(L,g):
    """ filter authorized Gi edges output 

    Parameters
    ----------

    L : Layout
    g : Digraph Gi

    Notes 
    -----

    Let assume a sequence (nstr0,nstr1,{nstr2A,nstr2B,...}) in a signature.
    This function checks that this sequence is feasible
    , whatever the type of nstr0 and nstr1.
    The feasible outputs from nstr0 to nstr1 are stored in an output field of 
    edge (nstr0,nstr1)


    """

    # loop over all edges of Gi
    for e in g.edges():
        # extract  both termination interactions nodes
        i0 = eval(e[0])
        i1 = eval(e[1])
        try:
            nstr0 = i0[0]
        except:
            nstr0 = i0


        try:
            nstr1 = i1[0]
            # Transmission
            if len(i1)>2:
                typ=2
            # Reflexion    
            else :
                typ=1
        # Diffraction        
        except:
            nstr1 = i1
            typ = 3

        # list of authorized outputs, initialized void
        output = []
        # nstr1 : segment number of final interaction
        if nstr1>0:
            pseg1 = L.seg2pts(nstr1).reshape(2,2).T
            cn = cone.Cone()
            if nstr0>0:
                pseg0 = L.seg2pts(nstr0).reshape(2,2).T
                # test if nstr0 and nstr1 are connected segments
                if (len(np.intersect1d(nx.neighbors(L.Gs,nstr0),nx.neighbors(L.Gs,nstr1)))==0):
                    # not connected
                    cn.from2segs(pseg0,pseg1)
                else:
                    # connected 
                    cn.from2csegs(pseg0,pseg1)
            else:
                pt = np.array(L.Gs.pos[nstr0])
                cn.fromptseg(pt,pseg1)
        
            # list all potential successor of interaction i1
            i2 = nx.neighbors(g,str(i1))
            ipoints = filter(lambda x: eval(x)<0 ,i2)
            istup = filter(lambda x : type(eval(x))==tuple,i2)
            isegments = np.unique(map(lambda x : eval(x)[0],istup))
            if len(isegments)>0:
                points = L.seg2pts(isegments)
                pta = points[0:2,:]
                phe = points[2:,:]
                #print points
                #print segments 
                #cn.show()
                if len(i1)==3:
                    bs = cn.belong_seg(pta,phe)
                    #if bs.any():
                    #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                    #if ~bs.any():
                    #    plu.displot(pta[:,~bs],phe[:,~bs],color='k')
                if len(i1)==2:    
                    Mpta = geu.mirror(pta,pseg1[:,0],pseg1[:,1])
                    Mphe = geu.mirror(phe,pseg1[:,0],pseg1[:,1])
                    bs = cn.belong_seg(Mpta,Mphe)
                    #print i0,i1
                    #if ((i0 == (6, 0)) & (i1 == (7, 0))):
                    #    pdb.set_trace()
                    #if bs.any():
                    #    plu.displot(pta[:,bs],phe[:,bs],color='g')
                    #if ~bs.any():
                    #    plu.displot(pta[:,~bs],phe[:,~bs],color='m')
                    #    plt.show()
                    #    pdb.set_trace()
                isegkeep = isegments[bs]     
                output = filter(lambda x : eval(x)[0] in isegkeep ,istup)
                # keep all segment above nstr1 and in Cone if T 
                # keep all segment below nstr1 and in Cone if R 

        g.add_edge(str(i0),str(i1),output=output)

    return(g)
def edgeout(L,g):
    """ filter authorized Gi edges output 

    Parameters
    ----------

    L : Layout
    g : Digraph Gi

    Notes 
    -----

    Let assume a sequence (nstr0,nstr1,{nstr2A,nstr2B,...}) in a signature.
    This function checks that this sequence is feasible
    , whatever the type of nstr0 and nstr1.
    The feasible outputs from nstr0 to nstr1 are stored in an output field of 
    edge (nstr0,nstr1)


    """

    # loop over all edges of Gi
    for e in g.edges():
        # extract  both termination interactions nodes
        i0 = eval(e[0])
        i1 = eval(e[1])
        try:
            nstr0 = i0[0]
        except:
            nstr0 = i0


        try:
            nstr1 = i1[0]
            # Transmission
            if len(i1)>2:
                typ=2
            # Reflexion    
            else :
                typ=1
        # Diffraction        
        except:
            nstr1 = i1
            typ = 3

        # list of authorized outputs, initialized void
        output = []
        # nstr1 : segment number of final interaction
        if nstr1>0:
            #cn = cone.Cone()
            #cn.from2segs(pseg0,pseg1)
            # segment unitary vector
            # l1 : unitary vector along structure segments  
            l1 = L.seguv(np.array([nstr1]))
            #
            # unitary vector along the ray (nstr0,nstr1)
            #
            p0 = np.array(L.Gs.pos[nstr0])
            p1 = np.array(L.Gs.pos[nstr1])
            v01  = p1-p0
            v01m = np.sqrt(np.dot(v01,v01))
            v01n = v01/v01m
            v10n = -v01n
            # next interaction
            # considering all neighbors of i1 in Gi 
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

                # if (reflexion is forward) or (reflexion return to its origin)
                if (d1*d2>=0) or (nstr0 == nstr2) and typ == 1:
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
    """ gathers all signatures from a layout given tx and rx

    Attributes
    ----------

    L : gis.Layout
    pTx : numpy.ndarray
        position of Tx
    pRx : numpy.ndarray
        position of Rx

    """

    def __init__(self,L,source,target,cutoff=3):
        """

        Parameters
        ----------

        L : Layout
        source : int 
            cycle number 
        target : int 
            cycle index

        """
        self.L = L
        self.source = source
        self.target = target
        self.cutoff = cutoff
        self.filename = self.L.filename.split('.')[0] +'_' + str(self.source) +'_' + str(self.target) +'_' + str(self.cutoff) +'.sig'

    def __repr__(self):
        def fun1(x):
            if x==1:
                return('R')
            if x==2:
                return('T')
            if x==3:
                return('D')
        size = {}
        s = self.__class__.__name__ + '\n' + '----------'+'\n'
        #s = s + str(self.__sizeof__())+'\n'
        for k in self:
            size[k] = len(self[k])/2
        s = s + 'from cycle : '+ str(self.source) + ' to cycle ' + str(self.target)+'\n'
        for k in self:
            s = s + str(k) + ' : ' + str(size[k]) + '\n'
            a = np.swapaxes(self[k].reshape(size[k],2,k),0,2)
            # nl x 2 x nsig
            for i in range(k):
                s = s + '   '+ str(a[i,0,:]) + '\n'
                s = s + '   '+ str(a[i,1,:]) + '\n'

        return(s)

    def __len__(self):
        nsig = 0
        for k in self:
            size = len(self[k])/2
            nsig += size
        return(nsig)

    def num(self):
        """ determine the number of signatures
        """
        self.nsig = 0
        self.nint = 0
        for k in self:
            size = len(self[k])/2
            self.nsig += size
            self.nint += size*k

    def info(self):
        """
        """
        # print "Signatures for scenario defined by :"
        # print "Layout"
        # print "======"
        # L = self.L.info()
        # print "================================"
        # print "source : ", self.source
        # print "target : ", self.target
        size = {}
        print self.__class__.__name__ + '\n' + '----------'+'\n'
        #s = s + str(self.__sizeof__())+'\n'
        for k in self:
            size[k] = len(self[k])/2
        print 'from cycle : '+ str(self.source) + ' to cycle ' + str(self.target)+'\n'
        pyu.printout('Reflection',pyu.BLUE)
        print '  '
        pyu.printout('Transmission',pyu.GREEN)
        print '  '
        pyu.printout('Diffraction',pyu.RED)
        print '  \n'
        for k in self:
            print str(k) + ' : ' + str(size[k]) 
            a = np.swapaxes(self[k].reshape(size[k],2,k),0,2)
            # nl x 2 x nsig
            for i in range(k):

                nstr=a[i,0,:]
                typ=a[i,1,:]
                print '[',
                for n,t in zip(nstr,typ):
                    if t==1:
                        pyu.printout(str(n),pyu.BLUE)
                    if t==2:
                        pyu.printout(str(n),pyu.GREEN)
                    if t==3:
                        pyu.printout(str(n),pyu.RED)
                print ']'
            print'\n'
                # s = s + '   '+ str(a[i,0,:]) + '\n'

                # s = s + '   '+ str(a[i,1,:]) + '\n'

    def save(self):
        """ Save signatures
        """
        L=copy.deepcopy(self.L)
        del(self.L)
        filename=pyu.getlong(self.filename,pstruc['DIRSIG'])
        with open(filename, 'wb') as handle:
          pickle.dump(self, handle)
        self.L=L

    def load(self,filename=[]):
        """ Load signatures
        """


        if filename == []:
            _filename = self.filename
        else :
            _filename = filename

        filename=pyu.getlong(_filename,pstruc['DIRSIG'])
        try:
            handle=open(filename, 'rb')
            sitmp = pickle.load(handle)
        except: 
            raise NameError(filename +' does not exist')


        # to load a dictionary, use update 
        self.update(sitmp)
        

        _fileL=pyu.getshort(filename).split('_')[0]+'.ini'
        self.L=layout.Layout(_fileL)
        try:
            self.L.dumpr()
        except:
            self.L.build()
            self.L.dumpw()

    def sp(self,G, source, target, cutoff=None):
        """ 

        """
        if cutoff < 1:
            return
        visited = [source]
        stack = [iter(G[source])]
        while stack:
            children = stack[-1]
            child = next(children, None)
            if child is None:
                stack.pop()
                visited.pop()
            elif len(visited) < cutoff:
                if child == target:
                    for i in range(len(self.ds[source])):
                        s=self.ds[target][i] + visited
                        self.ds[target].append(s)

                    # yield visited +[target]
                elif child not in visited:
                    visited.append(child)
                    stack.append(iter(G[child]))
            else: #len(visited) == cutoff:
                if child == target or target in children:
                    for i in range(len(self.ds[source])):
                        s=self.ds[target][i] + visited
                        self.ds[target].append(s)

                stack.pop()
                visited.pop()



    def procone(self,L, G, source, target, cutoff=1):
        """ seek all simple_path from source to target looking backward

        Parameters
        ----------
        
        L : Layout
        G : networkx Graph Gi
        source : tuple 
            interaction (node of Gi) 
        target : tuple 
            interaction (node of Gi) 
        cutoff : int

        Notes
        -----

        adapted from all_simple_path of networkx 

        1- Determine all nodes connected to Gi 

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
            elif len(visited) < cutoff: # if visited list length is less than cutoff
                if child == target:  # if child is the target point - YIELD A SIGNATURE
                    #print visited + [target]
                    yield visited + [target] # output signature
                else:
                #elif child not in visited: # else visit other node - CONTINUE APPEND CHILD
                    # getting signature until last point
                    diff  = np.where(np.array(visited)<0)[0]
                    if len(diff)==0:
                        brin = visited
                    else:
                        brin = visited[diff[-1]:]
                    # looking backward with a cone
                    if len(brin)>2:
                        # warning visited is also appended visited[-2] is the
                        # last node
                        brin.append(child)
                        s = Signature(brin)
                        s.evf(L)
                        ta,he = s.unfold()
                        cn = cone.Cone()
                        segchild = np.vstack((ta[:,-1],he[:,-1])).T
                        segvm1 = np.vstack((ta[:,-2],he[:,-2])).T
                        cn.from2segs(segchild,segvm1)
                        typ,proba = cn.belong_seg(ta[:,:-2],he[:,:-2])
                        #fig,ax = plu.displot(ta,he)
                        #fig,ax = cn.show(fig=fig,ax=ax)
                        #plt.show()
                        #pdb.set_trace()
                        if (typ==0).any():
                        # child no valid (do nothing)
                            visited.pop()
                        else:
                        # child valid (append child to visited and go forward)
                            stack.append(iter(G[visited[-2]][child]['output']))
                    else:
                        stack.append(iter(G[visited[-1]][child]['output']))
                        visited.append(child)

            else: #len(visited) == cutoff (visited list is too long)
                if child == target or target in children:
                    #print visited + [target]
                    yield visited + [target]
                stack.pop()
                visited.pop()

    def propaths(self,G, source, target, cutoff=1):
        """ seek all simple_path from source to target

        Parameters
        ----------

        G : networkx Graph Gi
        source : tuple
            interaction (node of Gi)
        target : tuple
            interaction (node of Gi)
        cutoff : int

        Notes
        -----

        adapted from all_simple_path of networkx

        1- Determine all nodes connected to Gi 

        """
        #print "source :",source
        #print "target :",target

        if cutoff < 1:
            return

        visited = [source]
        # stack is a list of iterators
        stack = [iter(G[source])]
        # lawp = list of airwall position in visited
        lawp = []

        # while the list of iterators is not void
        # import ipdb
        # ipdb.set_trace()
        while stack: #
            # children is the last iterator of stack

            children = stack[-1]
            # next child
            child = next(children, None)
            # update number of useful segments
            # if there is airwall in visited
            #

            if child is None  : # if no more child
                stack.pop()   # remove last iterator
                visited.pop() # remove from visited list
                try:
                    lawp.pop()
                except:
                    pass

            elif (len(visited) < (cutoff + sum(lawp))):# if visited list length is less than cutoff
                if child == target:  # if child is the target point
                    #print visited + [target]
                    yield visited + [target] # output signature
                elif child not in visited: # else visit other node
                    # only visit output nodes
                    #pdb.set_trace()
                    try:
                        dintpro = G[visited[-1]][child]['output']
                    except:
                        dintpro ={}

                    stack.append(iter(dintpro.keys()))
                    #stack.append(iter(G[visited[-1]][child]['output']))
                    visited.append(child)
                    # check if child (current segment) is an airwall
                    if self.L.di[child][0] in self.L.name['AIR']:
                        lawp.append(1)
                    else:
                        lawp.append(0)



            else: #len(visited) == cutoff (visited list is too long)
                if child == target or target in children:
                    #print visited + [target]
                    yield visited + [target]

                stack.pop()
                visited.pop()
                try:
                    lawp.pop()
                except:
                    pass


    # def propaths(self,G, source, target, cutoff=1, cutprob =0.5):
    #     """ seek all simple_path from source to target

    #     Parameters
    #     ----------

    #     G : networkx Graph Gi
    #     source : tuple 
    #         interaction (node of Gi) 
    #     target : tuple 
    #         interaction (node of Gi) 
    #     cutoff : int

    #     Notes
    #     -----

    #     adapted from all_simple_path of networkx 

    #     1- Determine all nodes connected to Gi 

    #     """
    #     #print "source :",source
    #     #print "target :",target

    #     if cutoff < 1:
    #         return

    #     visited = [source]
    #     # stack is a list of iterators
    #     stack = [iter(G[source])]
    #     ps = [iter([1.0]*len((G[source])))] 
    #     # lawp = list of airwall position in visited
    #     lawp = []

    #     # while the list of iterators is not void
    #     # import ipdb
    #     # ipdb.set_trace()    
    #     while stack: #
    #         # children is the last iterator of stack

    #         children = stack[-1]
    #         pcd = ps[-1]
    #         # next child
    #         child = next(children, None)
    #         pc = next(pcd,None)
    #         # update number of useful segments
    #         # if there is airwall in visited
    #         # 
            
    #         if child is None  : # if no more child
    #             stack.pop()   # remove last iterator
    #             ps.pop()
    #             visited.pop() # remove from visited list
    #             try:
    #                 lawp.pop()
    #             except:
    #                 pass

    #         elif (pc>cutprob): # check proba
    #             if (len(visited) < (cutoff + sum(lawp))):# if visited list length is less than cutoff 
    #                 if child == target:  # if child is the target point
    #                     #print visited + [target]
    #                     yield visited + [target] # output signature
    #                 elif child not in visited: # else visit other node
    #                     # only visit output nodes
    #                     #pdb.set_trace()
    #                     try:
    #                         dintpro = G[visited[-1]][child]['output']
    #                     except:
    #                         dintpro ={}

    #                     # pnc : probability of next children
    #                     # pc : proba of current parent
    #                     # spnc : sum of proba of next children

    #                     # spnc = sum(dintpro.values())
    #                     pnc = [(v*pc) for v in dintpro.values()]

    #                     stack.append(iter(dintpro.keys()))
    #                     ps.append(iter(pnc))
    #                     #stack.append(iter(G[visited[-1]][child]['output']))
    #                     visited.append(child)
    #                     # check if child (current segment) is an airwall
    #                     if self.L.di[child][0] in self.L.name['AIR']:
    #                         lawp.append(1)
    #                     else:
    #                         lawp.append(0)


    #             else :
    #                 stack.pop()
    #                 ps.pop()
    #                 visited.pop()
    #                 lawp.pop()

    #         else: #len(visited) == cutoff (visited list is too long)
    #             if child == target or target in children:
    #                 #print visited + [target]
    #                 yield visited + [target]

    #             stack.pop()
    #             ps.pop()
    #             visited.pop()
    #             lawp.pop()

    def calsig(self,G,dia={},cutoff=None):
        """

        Parameters
        ----------

        G   : graph
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

    def run(self,cutoff=1,dcut=2):
        """ run signature calculation 
        """

        lcil=self.L.cycleinline(self.source,self.target)

        if len(lcil) <= 2:
            print 'run1'
            self.run1(cutoff=cutoff)
        else :
            print 'run2'
            self.run2(cutoff=cutoff,dcut=dcut)


    def run1(self,cutoff=2):
        """ get signatures (in one list of arrays) between tx and rx

        Parameters
        ----------

        cutoff : int 
            limit the exploration of all_simple_path

        Returns
        -------

        sigslist :  numpy.ndarray

        """

        self.cutoff   = cutoff
        self.filename = self.L.filename.split('.')[0] +'_' + str(self.source) +'_' + str(self.target) +'_' + str(self.cutoff) +'.sig'

        try:
            self.L.dGi
        except:
            self.L.buildGi2()
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

        metasig = nx.neighbors(self.L.Gt,self.source)
        metasig = metasig + nx.neighbors(self.L.Gt,self.target)
        metasig = list(np.unique(np.array(metasig)))
        metasig = metasig + [self.source] + [self.target]
        
        #print "metasig",metasig

        # add cycles separated by air walls
        lca=[]
        for cy in metasig:
            try:
                lca.extend(self.L.dca[cy])
            except:
                pass
        metasig = metasig + lca
        metasig = list(np.unique(np.array(metasig)))

        lis = self.L.Gt.node[self.source]['inter']
        lit = self.L.Gt.node[self.target]['inter']

        # source
        #ndt1 = filter(lambda l: len(eval(l))>2,ndt) # Transmission
        #ndt2 = filter(lambda l: len(eval(l))<3,ndt) # Reflexion

        lisT = filter(lambda l: len(eval(l))>2,lis) # Transmission
        lisR = filter(lambda l: len(eval(l))<3,lis) # Reflexion

        # target
        # ndr1 = filter(lambda l: len(eval(l))>2,ndr) # Transmission
        # ndr2 = filter(lambda l: len(eval(l))<3,ndr) # Reflexion

        litT = filter(lambda l: len(eval(l))>2,lit) # Transmission
        litR = filter(lambda l: len(eval(l))<3,lit) # Reflexion

        # tx,rx : attaching rule
        #
        # tx attachs to out transmisision point
        # rx attachs to in transmission point

        #
        # WARNING : room number <> cycle number
        #

        #ncytx = self.L.Gr.node[NroomTx]['cycle']
        #ncyrx = self.L.Gr.node[NroomRx]['cycle']

        #ndt1 = filter(lambda l: eval(l)[2]<>ncytx,ndt1)
        #ndr1 = filter(lambda l: eval(l)[1]<>ncyrx,ndr1)

        lisT = filter(lambda l: eval(l)[2]<>self.source,lisT)
        litT = filter(lambda l: eval(l)[1]<>self.target,litT)

        #ndt = ndt1 + ndt2
        #ndr = ndr1 + ndr2
        lis  = lisT + lisR
        lit  = litT + litR

        #ntr = np.intersect1d(ndt, ndr)
#        li = np.intersect1d(lis, lit)

        li = []
        for ms in metasig:
            li = li + self.L.Gt.node[ms]['inter']
        li = list(np.unique(np.array(li)))

        dpos = {k:self.L.Gi.pos[k] for k in li}

        Gi = nx.subgraph(self.L.Gi,li)
        Gi.pos = dpos
#        for meta in metasig:
#        Gi = nx.DiGraph()
#        for cycle in metasig:
#            Gi = nx.compose(Gi,self.L.dGi[cycle])

#        # facultative update positions
#        Gi.pos = {}
#        for cycle in metasig:
#            Gi.pos.update(self.L.dGi[cycle].pos)
#        pdb.set_trace()
        #
        #
        #
        # remove diffractions from Gi
        Gi = gidl(Gi)
        # add 2nd order output to edges
        #Gi = edgeout(self.L,Gi)
        Gi = edgeout2(self.L,Gi)
        #pdb.set_trace()
        #for interaction source  in list of source interaction 
        for s in lis:
            #for target interaction in list of target interaction
            for t in lit:

                if (s != t):
                    #paths = list(nx.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
                    #paths = list(self.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
                    paths = list(self.propaths(Gi,source=s,target=t,cutoff=cutoff))

                    #paths = [nx.shortest_path(Gi,source=s,target=t)]
                else:
                    #paths = [[nt]]
                    paths = [[s]]
                ### supress the followinfg loops .
                for path in paths:

                    sigarr = np.array([],dtype=int).reshape(2, 0)
                    for interaction in path:

                        it = eval(interaction)
                        if type(it) == tuple:
                            if len(it)==2: #reflexion
                                sigarr = np.hstack((sigarr,
                                                np.array([[it[0]],[1]],dtype=int)))
                            if len(it)==3: #transmission
                                sigarr = np.hstack((sigarr,
                                                np.array([[it[0]],[2]],dtype=int)))
                        elif it < 0: #diffraction
                            sigarr = np.hstack((sigarr,
                                                np.array([[it],[3]],dtype=int)))
                    #print sigarr
                    try:
                        self[len(path)] = np.vstack((self[len(path)],sigarr))
                    except:
                        self[len(path)] = sigarr

    def run4(self,cutoff=2,algo='new',progress=False):
        """ get signatures (in one list of arrays) between tx and rx

        Parameters
        ----------

        cutoff : int
            limit the exploration of all_simple_path

        Returns
        -------

        sigslist :  numpy.ndarray

        """

        self.cutoff   = cutoff
        self.filename = self.L.filename.split('.')[0] +'_' + str(self.source) +'_' + str(self.target) +'_' + str(self.cutoff) +'.sig'

        # Determine meta signature
        # this limits the number of cycles

        #metasig = nx.neighbors(self.L.Gt,self.source)
        #metasig = metasig + nx.neighbors(self.L.Gt,self.target)
        #metasig = list(np.unique(np.array(metasig)))
        #metasig = metasig + [self.source] + [self.target]
        # add cycles separated by air walls
        #lca=[]
        #for cy in metasig:
        #    try:
        #        lca.extend(self.L.dca[cy])
        #    except:
        #        pass
        #metasig = metasig + lca
        #metasig = list(np.unique(np.array(metasig)))

        # list of interaction source
        lis = self.L.Gt.node[self.source]['inter']
        # list of interaction target
        lit = self.L.Gt.node[self.target]['inter']

        # source
        #ndt1 = filter(lambda l: len(eval(l))>2,ndt) # Transmission
        #ndt2 = filter(lambda l: len(eval(l))<3,ndt) # Reflexion

        lisT = filter(lambda l: len(eval(l))>2,lis) # Transmission
        lisR = filter(lambda l: len(eval(l))<3,lis) # Reflexion

        # target
        # ndr1 = filter(lambda l: len(eval(l))>2,ndr) # Transmission
        # ndr2 = filter(lambda l: len(eval(l))<3,ndr) # Reflexion

        litT = filter(lambda l: len(eval(l))>2,lit) # Transmission
        litR = filter(lambda l: len(eval(l))<3,lit) # Reflexion

        # tx,rx : attaching rule
        #
        # tx attachs to out transmisision point
        # rx attachs to in transmission point

        #
        # WARNING : room number <> cycle number
        #

        #ncytx = self.L.Gr.node[NroomTx]['cycle']
        #ncyrx = self.L.Gr.node[NroomRx]['cycle']

        #ndt1 = filter(lambda l: eval(l)[2]<>ncytx,ndt1)
        #ndr1 = filter(lambda l: eval(l)[1]<>ncyrx,ndr1)

        lisT = filter(lambda l: eval(l)[2]<>self.source,lisT)
        litT = filter(lambda l: eval(l)[1]<>self.target,litT)

        #ndt = ndt1 + ndt2
        #ndr = ndr1 + ndr2
        # list of interaction visible from source
        lis  = lisT + lisR
        # list of interaction visible from target
        lit  = litT + litR

        #ntr = np.intersect1d(ndt, ndr)
#        li = np.intersect1d(lis, lit)

        # list of all interactions
        #li = []
        #for ms in metasig:
        #    li = li + self.L.Gt.node[ms]['inter']
        #li = list(np.unique(np.array(li)))
        #
        # dictionnary interaction:position
        #dpos = {k:self.L.Gi.pos[k] for k in li}
        # extracting sub graph of Gi corresponding to metasiganture
        #Gi = nx.subgraph(self.L.Gi,li)
        #Gi.pos = dpos
        Gi = self.L.Gi
        Gi.pos = self.L.Gi.pos
#        for meta in metasig:
#        Gi = nx.DiGraph()
#        for cycle in metasig:
#            Gi = nx.compose(Gi,self.L.dGi[cycle])

#        # facultative update positions
#        Gi.pos = {}
#        for cycle in metasig:
#            Gi.pos.update(self.L.dGi[cycle].pos)
#        pdb.set_trace()
        #
        # TODO : This has to be changed for handling diffraction
        # 
        # remove diffractions from Gi
        Gi = gidl(Gi)
        # add 2nd order output to edges
        #Gi = edgeout(self.L,Gi)
        #Gi = edgeout2(self.L,Gi)
        #pdb.set_trace()
        lmax = len(lis)*len(lit)
        pe = 0
        tic = time.time()
        tic0 = tic
        #for interaction source  in list of source interaction
        for us,s in enumerate(lis):
            #for target interaction in list of target interaction
            for ut,t in enumerate(lit):

                if progress :
                    ratio = np.round((((us)*len(lis)+ut)/(1.*lmax))*10 )
                    if ratio != pe:
                        pe = ratio
                        toc = time.time()
                        print '~%d ' % (ratio*10),
                        print '%',
                        print '%6.3f %6.3f' % (toc-tic, toc-tic0)
                        tic = toc
                if (s != t):
                    #paths = list(nx.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
                    #paths = list(self.all_simple_paths(Gi,source=s,target=t,cutoff=cutoff))
                    if algo=='new':
                        paths = list(self.procone(self.L,Gi,source=s,target=t,cutoff=cutoff))
                    else:
                        paths = list(self.propaths(Gi,source=s,target=t,cutoff=cutoff))

                    #paths = [nx.shortest_path(Gi,source=s,target=t)]
                else:
                    #paths = [[nt]]
                    paths = [[s]]
                ### suppress the following loops .
                for path in paths:

                    sigarr = np.array([],dtype=int).reshape(2, 0)
                    for interaction in path:
                        #print interaction + '->',
                        it = eval(interaction)
                        if type(it) == tuple:
                            if len(it)==2: #reflexion
                                sigarr = np.hstack((sigarr,
                                                np.array([[it[0]],[1]],dtype=int)))
                            if len(it)==3: #transmission
                                sigarr = np.hstack((sigarr,
                                                np.array([[it[0]],[2]],dtype=int)))
                        elif it < 0: #diffraction
                            sigarr = np.hstack((sigarr,
                                                np.array([[it],[3]],dtype=int)))
                    #print sigarr
                    #print ''
                    try:
                        self[len(path)] = np.vstack((self[len(path)],sigarr))
                    except:
                        self[len(path)] = sigarr

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
##
##      1. the sequence of cycle between cycle source cs and
##      cycle target ct are obtained via cycleinline method
##
##      2. all cycles adjscent at least to one of the previously
##      obtained cycles are appended to the list lca (list of cycle around)
##
##      3. It is then required to add all cycles
##      connected to the previous ones via an air wall.
##
##      lca is used to build the sliding graph of interactions
##      it is important that lcil remains ordered this is not the case
##      for lca



        cs = self.source
        ct = self.target
        lcil=self.L.cycleinline(cs,ct)
        lca = [] # list of cycle around
        for cy in lcil:
            ncy = nx.neighbors(self.L.Gt,cy)
            lca = lca+ncy
        lca = list(np.unique(np.array(lca)))
        lca = lcil
        lcair=[]
        for cy in lca:
            try:
                lcair.extend(self.L.dca[cy])
            except:
                pass
        lca = lca + lcair
        lca = list(np.unique(np.array(lca)))
        #
        # extract list of interactions from list of cycles lca
        #
        li = []
        for ms in lca:
            li = li + self.L.Gt.node[ms]['inter']
        # enforce unicity of interactions in list li
        li = list(np.unique(np.array(li)))

        # extract dictionnary of interactions position
        dpos = {k:self.L.Gi.pos[k] for k in li}

        # build the subgraph of L.Gi with only the selected interactions
        Gi = nx.subgraph(self.L.Gi,li)
        Gi.pos = dpos

        Gf = nx.DiGraph()
        Gf.pos = {}
        # remove diffractions points from Gi
        Gi = gidl(Gi)
        # add 2nd order output to edges
        Gi = edgeout(self.L,Gi)
        #for interaction source  in list of source interaction

############################################################
#        filter list of interactions in termination cycles

        # list of interactions belonging to source
        lis = self.L.Gt.node[lcil[0]]['inter']

        # list of interactions belonging to target
        lit = self.L.Gt.node[lcil[-1]]['inter']

        # filter lis remove transmission coming from outside
        lli   = []
        lisR  = filter(lambda l: len(eval(l))==2,lis)
        lisT  = filter(lambda l: len(eval(l))==3,lis)
        lisTo = filter(lambda l: eval(l)[2]<>cs,lisT)
        lli = lisR + lisTo
        # for li in lis:
        #     ei = eval(li)
        #     if len(ei)==2:
        #        lli.append(li)
        #    if len(ei)==3:
        #        if ei[2]<>cs:
        #           lli.append(li)
        # filter lit remove transmission going outside
        llt = []
        litR  = filter(lambda l: len(eval(l))==2,lit)
        litT  = filter(lambda l: len(eval(l))==3,lit)
        litTi = filter(lambda l: eval(l)[2]==ct,litT)
        llt = litR + litTi
        #for li in lit:
        #    ei = eval(li)
        #    if len(ei)==2:
        #        llt.append(li)
        #    if len(ei)==3:
        #        if ei[2]==ct:
        #           llt.append(li)
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

    def run3(self,cutoff=1,dcut=2):
        """ get signatures (in one list of arrays) between tx and rx

        Parameters
        ----------

            cutoff : limit the exploration of all_simple_path

        Returns
        -------

            sigslist = numpy.ndarray

        """

############################################################
##
##      1. the sequence of cycle between cycle source cs and
##      cycle target ct are obtained via cycleinline method
##
##      2. all cycles adjscent at least to one of the previously
##      obtained cycles are appended to the list lca (list of cycle around)
##
##      3. It is then required to add all cycles
##      connected to the previous ones via an air wall.
##
##      lca is used to build the sliding graph of interactions
##      it is important that lcil remains ordered this is not the case
##      for lca

        # cs : cycle source
        cs = self.source
        # ct : cycle target
        ct = self.target
        polys = self.L.Gt.node[cs]['polyg']
        # cps : centroid point source
        cps = polys.centroid.xy
        polyt = self.L.Gt.node[ct]['polyg']
        # cpt : centroid point target
        cpt = polyt.centroid.xy
        ps = np.array([cps[0][0],cps[1][0]])
        pt = np.array([cpt[0][0],cpt[1][0]])
        v = pt-ps
        mv = np.sqrt(np.sum(v*v,axis=0))
        vn = v/mv
        lcil = self.L.cycleinline(cs,ct)

        # dac : dictionary of adjascent cycles
        dac = {}
        # dfl : dictionnary of fronlines
        dfl = {}

        for cy in lcil:
            dfl[cy] = []

            nghb = nx.neighbors(self.L.Gt,cy)
            dac[cy] = nghb
            poly1 = self.L.Gt.node[cy]['polyg']
            cp1 = poly1.centroid.xy
            p1 = np.array([cp1[0][0],cp1[1][0]])

            for cya in nghb:
                if cy == 13:
                    pdb.set_trace()
                poly2 = self.L.Gt.node[cya]['polyg']
                cp2 = poly2.centroid.xy
                p2 = np.array([cp2[0][0],cp2[1][0]])
                vp = p2-p1
                m2 = np.sqrt(np.sum((p2-p1)*(p2-p1),axis=0))
                vpn = vp/m2
                dp = np.dot(vpn,vn)
                # if dot(vn,vpn) >0 cycle cya is ahead
                if dp>0:
                    lsegso = frontline(self.L,cy,vn)
                    for s in lsegso:
                        cyb = filter(lambda n : n <> cy,self.L.Gs.node[s]['ncycles'])
                        if cyb<>[]:
                            dfl[cy].append(str((s,cy,cyb[0])))
            dfl[cy]=np.unique(dfl[cy]).tolist()
        # # list of interactions belonging to source
        # lis = self.L.Gt.node[lcil[0]]['inter']

        # # list of interactions belonging to target
        # lit = self.L.Gt.node[lcil[-1]]['inter']

        # # filter lis remove incoming transmission
        # lli   = []
        # lisR  = filter(lambda l: len(eval(l))==2,lis)
        # lisT  = filter(lambda l: len(eval(l))==3,lis)
        # lisTo = filter(lambda l: eval(l)[2]<>cs,lisT)
        # lis = lisR + lisTo

        # # filter lit remove outgoing transmission
        # llt = []
        # litR  = filter(lambda l: len(eval(l))==2,lit)
        # litT  = filter(lambda l: len(eval(l))==3,lit)
        # litTi = filter(lambda l: eval(l)[2]==ct,litT)
        # lit = litR + litTi
        

        self.ds={}

        for icy in range(len(lcil)-1):
            

            io = dfl[lcil[icy]]
            io_ = dfl[lcil[icy+1]]
            print io
            print io_
            if icy == 0:
                [self.ds.update({k:[[k]]}) for k in io]        

            # remove keys which are not in front line     
            # kds = self.ds.keys()
            # for k in kds :
            #     if k not in io:
            #         self.ds.pop(k)

            for j in io_:
                self.ds[j]=[[]]
                for i in io:
                    self.sp(self.L.Gi,i,j,cutoff=2)
            # [self.ds.pop(k) for k in io]
                    
                    # ds[j]
                    # if len(a) == 1:
                    #     if len(ds[j]) <> 0:
                    #         pdb.set_trace()
                    #         [ds[j][k].extend(a[0][:-1]) for k in range(len(ds[j]))]
                    #     else :    
                    #         ds[j]=a[0][:-1]
                    # elif len(a)> 1:
                    #     if len(ds[j]) <> 0:
                    #         [[ds[j][k].extend(a[l][:-1]) for k in range(len(ds[j]))] for l in range(len(a))]
                    #     else :    
                    #         ds[j]=a[:-1]


            # remove segments which separate two cycles.
            # TODO: See if it worth to implement
            #lsegs = filter(lambda x : x not in interseg,lsegs)
        # add adjascent air cycles
        #lcair=[]
        #for cy in lcil:
        #    try:
        #        lcair.extend(self.L.dca[cy])
        #    except:
        #        pass
        #lca = lca + lcair
        #lca = list(np.unique(np.array(lca)))

        #
        # Reduction of Gi
        #

        #
        # extract list of interactions from list of cycles lca
        #

        li = []
        for cy in dac:
            for cya in dac[cy]:
                li = li + self.L.Gt.node[cya]['inter']
        # enforce unicity of interactions in list li
        li = list(np.unique(np.array(li)))

        # extract dictionnary of interactions position
        dpos = {k:self.L.Gi.pos[k] for k in li}

        # build the subgraph of L.Gi with only the selected interactions
        Gi = nx.subgraph(self.L.Gi,li)
        Gi.pos = dpos

        # remove diffractions points from Gi
        Gi = gidl(Gi)
        # add 2nd order output to edges
        Gi = edgeout(self.L,Gi)
        #for interaction source  in list of source interaction

############################################################
#        filter list of interactions in termination cycles

        # list of interactions belonging to source
        lis = self.L.Gt.node[lcil[0]]['inter']

        # list of interactions belonging to target
        lit = self.L.Gt.node[lcil[-1]]['inter']

        # filter lis remove incoming transmission
        lli   = []
        lisR  = filter(lambda l: len(eval(l))==2,lis)
        lisT  = filter(lambda l: len(eval(l))==3,lis)
        lisTo = filter(lambda l: eval(l)[2]<>cs,lisT)
        lis = lisR + lisTo

        # filter lit remove outgoing transmission
        llt = []
        litR  = filter(lambda l: len(eval(l))==2,lit)
        litT  = filter(lambda l: len(eval(l))==3,lit)
        litTi = filter(lambda l: eval(l)[2]==ct,litT)
        lit = litR + litTi

#################################################
#       propaths (a.k.a. all simple path) per adjacent cycles along cycles in line
#       Obtaining Gf: filtred graph of Gi with Gc ( rename Gt in Gc )

        #
        # Gf : filtered graph
        #
        Gf = nx.DiGraph()
        Gf.pos = {}
        ncycles = len(lcil)

        ltarget = []
        lsource = []
        for ic in np.arange(ncycles-1):

            # determine list of sources and targets
            # The famous so called saute mouton algorithm
            if ic==0:
                lsource = lis
                ltarget = dfl[lcil[0]]
            elif ic==ncycles-2:
                lsource = ltarget
                ltarget = lit
            else:
                lsource = ltarget
                ltarget = dfl[lcil[ic]]

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
        while True:
            culdesac = filter(lambda n: len(nx.neighbors(Gf,n))==0,Gf.nodes())
            culdesac.remove('Rx')
            if len(culdesac)>0:
                Gf.remove_nodes_from(culdesac)
                print culdesac
            else:
                break
        # a =[ 0,  1,  2,  1,  4,  1,  6,  1,  8,  1, 10, 1]
        # aa = np.array(a)
        # X=aa.reshape((2,3,2)) # r x i x 2
        # Y=X.swapaxes(0,2) # 2 x i x r
        self.Gf = Gf
        print 'signatures'
        co = nx.dijkstra_path_length(Gf,'Tx','Rx')
        #pdb.set_trace()
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

    def cones(self,L,i=0,s=0,fig=[],ax=[],figsize=(10,10)):
        """ display cones of an unfolded signature

        Parameters
        ----------

        L : Layout
        i : int
            the interaction block
        s : int
            the signature number in the block
        
        """
        if fig == []:
            fig= plt.figure()
            ax = fig.add_subplot(111)
        elif ax ==[]:
            ax = fig.add_subplot(111)

        
        pta,phe = self.unfold(L,i=i,s=s)
        
        # create a global array or tahe segments

        seg = np.vstack((pta,phe))
        lensi = np.shape(seg)[1]

        for s in range(1,lensi):
            pseg0 = seg[:,s-1].reshape(2,2).T
            pseg1 = seg[:,s].reshape(2,2).T
            #
            # create the cone seg0 seg1 
            #
            cn = cone.Cone()
            cn.from2segs(pseg0,pseg1)
            fig,ax = cn.show(fig = fig, ax = ax, figsize = figsize)

        return (fig,ax)


    def unfold(self,L,i=0,s=0):
        """ unfold a given signature

            return 2 np.ndarray of pta and phe "aligned" (reflexion interaction are mirrored) 

        Parameters
        ----------

        L : Layout
        i : int
            the interaction block
        s : int
            the signature number in the block
        
        See Also
        --------

        Signature.unfold

        """
        
        si = Signature(self[i][(2*s):(2*s)+2])
        si.ev(L)
        pta,phe = si.unfold()

        return pta,phe

    def show(self,L,**kwargs):
        """  plot signatures within the simulated environment


        Parameters
        ----------

        L : Layout
        i : list or -1 (default = all groups)
            list of interaction group numbers 
        s : list or -1 (default = all sig) 
            list of indices of signature in interaction group 
        ctx : cycle of tx (optional)
        crx : cycle of rx (optional)
        graph : type of graph to be displayed


        """
        defaults = {'i':-1,
                   's':-1,
                   'fig':[],
                   'ax':[],
                   'graph':'s',
                    'color':'black',
                    'alphasig':1,
                    'widthsig':0.1,
                    'colsig':'black',
                    'ms':5,
                    'ctx':-1,
                    'crx':-1
                   }
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        #if kwargs['fig'] ==[]:
        #    fig = plt.figure()
        #if kwargs['ax'] ==[]:    
        #    ax = fig.add_subplot(111)

        fig,ax = L.showG(**kwargs)
        if kwargs['ctx']!=-1:
            T=self.L.Gt.node[kwargs['ctx']]['polyg']
            T.coul='r'
            T.plot(fig=kwargs['fig'],ax=kwargs['ax'],color='r')
        if kwargs['crx']!=-1:
            R=self.L.Gt.node[kwargs['crx']]['polyg']
            R.plot(fig=kwargs['fig'],ax=kwargs['ax'],color='g')
        # i=-1 all rays
        # else block of interactions i
        if kwargs['i']==-1:
            lgrint = self.keys()
        else:
            lgrint = [kwargs['i']]
            

        for i in lgrint:
            if kwargs['s']==-1:
                lsig = range(len(self[i])/2)
            else:
                lsig = [kwargs['s']]
            for j in lsig: 
                sig=map(lambda x: self.L.Gs.pos[x],self[i][2*j])
                siga=np.array(sig)
                # sig = np.hstack((self.pTx[0:2].reshape((2, 1)),
                #                  np.hstack((self[i]['pt'][0:2, :, j],
                #                  self.pRx[0:2].reshape((2, 1))))
                #                  ))
                ax.plot(siga[:,0], siga[:,1],
                        alpha=kwargs['alphasig'],color=kwargs['colsig'],linewidth=kwargs['widthsig'])
                ax.axis('off')
        return(fig,ax)  

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


    def rays(self,ptx=0,prx=1):
        """ from signatures dict to 2D rays

        Parameters
        ----------

        tx : numpy.array or int
            Tx coordinates is the center of gravity of the cycle number if
            type(tx)=int
        rx :  numpy.array or int
            Rx coordinates is the center of gravity of the cycle number if
            type(rx)=int

        Returns
        -------

        rays : Rays

        Notes
        -----

        In the same time the signature of the ray is stored in the Rays object

        Todo : Find the best memory implemntation

        See Also
        --------

        Signature.sig2ray

        """

        if type(ptx)==int:
            ptx = np.array(self.L.Gt.pos[ptx])
        if type(prx)==int:
            prx = np.array(self.L.Gt.pos[prx])
        rays = Rays(ptx,prx)
        #
        # detect LOS situation
        #
        lc  = self.L.cycleinline(self.source,self.target)
        # if source  and target in the same cycle
        if len(lc) == 1:
            rays.los=True
#        # if source  and target separated by air walls
#        else :
#            for ic in xrange(len(lc)-1) :
#                # if lc[i] connected to a air wall
#                if lc[ic] in self.L.dca.keys():
#                    # if lc[i] and lc[i+1] are connected by a airwall
#                    if lc[ic+1] in self.L.dca[lc[ic]]:
#                        los = True
#                    else :
#                        los = False
#                        break
#                else :
#                    los = False
#                    break

#        rays[0]['pt']
        for k in self:
            # get signature block with k interactions
            tsig = self[k]
            shsig = np.shape(tsig)
            for l in range(shsig[0]/2):
                sig = tsig[2*l:2*l+2,:]
                s   = Signature(sig)
                isray,Yi  = s.sig2ray(self.L, ptx[:2], prx[:2])
                if isray:
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

        rays.nb_origin_sig = len(self)

        return rays

class Signature(object):
    """ class Signature

    Attributes
    ----------

    seq : list  of interaction point (edges (>0)  or vertices (<0) [int]
    typ : list of interaction type 1-R 2-T 3-D  [int] 
    pa  : tail point of interaction segmenti (2xN) ndarray
    pb  : head point of interaction segment  (2xN) ndarray
    pc  : center point of interaction segment (2xN) ndarray

    """
    def __init__(self, sig):
        """

        Parameters
        ----------

        sig : nd.array or list of interactions

        >>> seq = np.array([[1,5,1],[1,1,1]])
        >>> s = Signature(seq)

        """

        def typinter(l):
            try:
                l = eval(l)
            except:
                pass
            if type(l)==tuple:
                if len(l)==2:
                    return(1)
                if len(l)==3:
                    return(2)
            else:
                return(3)
        
        def seginter(l):
            try:
                l = eval(l)
            except:
                pass
            if type(l)==tuple:
                return l[0]
            else:
                return(l)

        if type(sig)==np.ndarray:
            self.seq = sig[0, :]
            self.typ = sig[1, :]

        if type(sig)==list:    
            self.seq = map(seginter,sig) 
            self.typ = map(typinter,sig) 

    def __repr__(self):
        #s = self.__class__ + ':' + str(self.__sizeof__())+'\n'
        s = ''
        s = s + str(self.seq) + '\n' 
        s = s + str(self.typ) + '\n'
        return s

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


    def ev2(self, L):
        """  evaluation of Signature

        Parameters
        ----------

        L : Layout

        Notes
        -----

        This function converts the sequence of interactions into numpy arrays
        which contains coordinates of segments extremities involved in the 
        signature. At that level the coordinates of extremities (tx and rx) is 
        not known yet.
        
        members data 

        pa  tail of segment  (2xN) 
        pb  head of segment  (2xN)  
        pc  the center of segment (2xN) 

        norm normal to the segment if segment 
        in case the interaction is a point the normal is undefined and then
        set to 0. 

        """
        def seqpointa(k,L=L):
            if k>0:
                ta, he = L.Gs.neighbors(k)
                pa = np.array(L.Gs.pos[ta]).reshape(2,1)
                pb = np.array(L.Gs.pos[he]).reshape(2,1)
                pc = np.array(L.Gs.pos[k]).reshape(2,1)
                nor1 = L.Gs.node[k]['norm']
                norm = np.array([nor1[0], nor1[1]]).reshape(2,1)
            else:    
                pa = np.array(L.Gs.pos[k]).reshape(2,1)
                pb = pa
                pc = pc
                norm = np.array([0, 0]).reshape(2,1)
            return(np.vstack((pa,pb,pc,norm)))

        v = np.array(map(seqpointa,self.seq))

        self.pa = v[:,0:2,:]
        self.pb = v[:,2:4,:]
        self.pc = v[:,4:6,:]
        self.norm = v[:,6:,:] 


    def evf(self, L):
        """  evaluation of Signature (fast version)

        Parameters
        ----------

        L : Layout

        Notes
        -----

        This function converts the sequence of interactions into numpy arrays
        which contains coordinates of segments extremities involved in the 
        signature. 
        
        members data 

        pa  tail of segment  (2xN) 
        pb  head of segment  (2xN)  


        """

        N = len(self.seq)
        self.pa = np.empty((2, N))  # tail
        self.pb = np.empty((2, N))  # head

        for n in range(N):
            k = self.seq[n]
            if k > 0:  # segment
                ta, he = L.Gs.neighbors(k)
                self.pa[:, n] = np.array(L.Gs.pos[ta])
                self.pb[:, n] = np.array(L.Gs.pos[he])
            else:      # node
                pa = np.array(L.Gs.pos[k])
                self.pa[:, n] = pa
                self.pb[:, n] = pa

    def ev(self, L):
        """  evaluation of Signature

        Parameters
        ----------

        L : Layout

        Notes
        -----

        This function converts the sequence of interactions into numpy arrays
        which contains coordinates of segments extremities involved in the 
        signature. At that level the coordinates of extremities (tx and rx) is 
        not known yet.
        
        members data 

        pa  tail of segment  (2xN) 
        pb  head of segment  (2xN)  
        pc  the center of segment (2xN) 

        norm normal to the segment if segment 
        in case the interaction is a point the normal is undefined and then
        set to 0. 

        """
        N = len(self.seq)
        self.pa = np.empty((2, N))  # tail
        self.pb = np.empty((2, N))  # head
        self.pc = np.empty((2, N))  # center
        self.norm = np.empty((2, N))

        for n in range(N):
            k = self.seq[n]
            if k > 0:  # segment
                ta, he = L.Gs.neighbors(k)
                norm1 = np.array(L.Gs.node[k]['norm'])
                norm = np.array([norm1[0], norm1[1]])
                self.pa[:, n] = np.array(L.Gs.pos[ta])
                self.pb[:, n] = np.array(L.Gs.pos[he])
                self.pc[:, n] = np.array(L.Gs.pos[k])
                self.norm[:, n] = norm
            else:      # node
                pa = np.array(L.Gs.pos[k])
                norm = np.array([0, 0])
                self.pa[:, n] = pa
                self.pb[:, n] = pa
                self.pc[:, n] = pa
                self.norm[:, n] = norm
    def unfold(self):
        """ unfold a given signature

            return 2 np.ndarray of pta and phe "aligned" 
            reflexion interactions are mirrored

        """
        
        lensi = len(self.seq)
        pta = np.empty((2,lensi))
        phe = np.empty((2,lensi))

        pta[:,0] = self.pa[:,0]
        phe[:,0] = self.pb[:,0]

        mirror=[]

        for i in range(1,lensi):

            pam = self.pa[:,i].reshape(2,1)
            pbm = self.pb[:,i].reshape(2,1)
                
            if self.typ[i] == 1: # R
                for m in mirror:
                    pam = geu.mirror(pam,pta[:,m],phe[:,m])
                    pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
                pta[:,i] = pam.reshape(2)
                phe[:,i] = pbm.reshape(2)
                mirror.append(i)

            elif self.typ[i] == 2 : # T
                for m in mirror:
                    pam = geu.mirror(pam,pta[:,m],phe[:,m])
                    pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
                pta[:,i] = pam.reshape(2)
                phe[:,i] = pbm.reshape(2)

        return pta,phe

    def evtx(self, L, tx, rx):
        """ evtx ( deprecated ) 

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
        """ Compute the tx's images with respect to the signature segments

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
        """ backtrace given image, tx, and rx
        
        this function
        traces the 2D ray.

        Parameters
        ----------

        tx :  ndarray (2x1)
              transmitter
        rx :  ndarray (2x1) 
              receiver
        M  :  ndarray (2xN) 
              N image points obtained using self.image method

        Returns
        -------

        isvalid : bool
            True if the backtrace ends successfully

        Y : ndarray (2 x (N+2))
            sequence of points corresponding to the seek ray

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
            >>> seq = np.array([[1,5,1],[1,1,1]])
            >>> s = Signature(seq)
            >>> tx = np.array([4,-1])
            >>> rx = np.array([1,1])
            >>> s.ev(L)
            >>> M = s.image(tx)
            >>> isvalid,Y = s.backtrace(tx,rx,M)
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

        Notes
        -----

        For mathematical details see : 

        @INPROCEEDINGS{6546704, 
        author={Laaraiedh, Mohamed and Amiot, Nicolas and Uguen, Bernard}, 
        booktitle={Antennas and Propagation (EuCAP), 2013 7th European Conference on}, 
        title={Efficient ray tracing tool for UWB propagation and
               localization modeling}, 
        year={2013}, 
        pages={2307-2311},}

        """

        pa = self.pa
        pb = self.pb
        typ = self.typ

        N = np.shape(pa)[1]
        I2 = np.eye(2)
        z0 = np.zeros((2, 1))

        pkm1 = rx.reshape(2, 1)
        Y = pkm1
        k = 0     # intercation counter
        beta = .5 # to enter into the loop
        isvalid = True # signature is valid by default 
        epsilon = 1e-2
        # while (((beta <= 1) & (beta >= 0)) & (k < N)):
        while (((beta <= 1-epsilon) & (beta >= epsilon)) & (k < N)):
            if int(typ[k]) != 3: # not a diffraction 
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
                    return(False,(k,None,None))
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
            return isvalid,Y
        else:
            isvalid = False 
            return isvalid,(k,alpha,beta) 


    def sig2beam(self, L, p, mode='incremental'):
        """
        """
        try:
            L.Gr
        except:
            L.build()
        
        # ev transforms a sequence of segment into numpy arrays (points)
        # necessary for image calculation
        self.ev(L)
        # calculates images from pTx
        M = self.image(pTx)
        

    def sig2ray(self, L, pTx, pRx, mode='incremental'):
        """ convert a signature to a 2D ray

        Parameters
        ----------

        L : Layout
        pTx : ndarray
            2D transmitter position
        pRx : ndarray
            2D receiver position
        mod : if mod=='incremental' a set of alternative signatures is return

        Returns
        -------

        Y : ndarray (2x(N+2))

        See Also 
        --------

        Signature.image
        Signature.backtrace
            
        """
        try:
            L.Gr
        except:
            L.build()
        
        # ev transforms a sequence of segment into numpy arrays (points)
        # necessary for image calculation
        self.ev(L)
        # calculates images from pTx
        M = self.image(pTx)
        
        #print self
        #if np.array_equal(self.seq,np.array([5,7,4])):
        #    pdb.set_trace()
        isvalid,Y = self.backtrace(pTx, pRx, M)
        #print isvalid,Y
        # 
        # If incremental mode this function returns an alternative signature
        # in case the signature do not yield a valid ray.
        #
        isray = True
        if mode=='incremental':
            if isvalid:
                return isray,Y
            else:
                isray=False
                # something to do here
                return isray,None
        else:
            if isvalid:
                return isray,Y
            else:
                isray=False
                return isray,None
            
 
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
    plt.ion()
    print "testing pylayers/antprop/signature.py"
    doctest.testmod()
    print "-------------------------------------"
