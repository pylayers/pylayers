#-*- coding:Utf-8 -*-
"""
Class Signatures
================

.. autosummary::
    :toctree: generated/

    Signatures.__init__
    Signatures.__repr__
    Signatures.__len__
    Signatures.num
    Signatures.info
    Signatures.saveh5
    Signatures.loadh5
    Signatures._saveh5
    Signatures._loadh5
    Signatures.load
    Signatures.save
    Signatures.sp
    Signatures.procone
    Signatures.propaths
    Signatures.short_propath
    Signatures.propaths2
    Signatures.propaths3
    Signatures.propaths2015_2
    Signatures.procone2
    Signatures.calsig
    Signatures.exist
    Signatures.run2015
    Signatures.run2015_2
    Signatures.dido
    Signatures.run
    Signatures.run_old
    Signatures.run_exp
    Signatures.run_exp2
    Signatures.meta
    Signatures.lineofcycle
    Signatures.cones
    Signatures.unfold
    Signatures.show
    Signatures.showi
    Signatures.rays
    Signatures.raysv
    Signatures.image
    Signatures.image2

Class Signature
===============

.. autosummary::
    :toctree: generated/

    Signature.__init__
    Signature.__repr__
    Signature.info
    Signature.split
    Signature.ev2
    Signature.evf
    Signature.ev
    Signature.unfold
    Signature.evtx
    Signature.image
    Signature.backtrace
    Signature.sig2beam
    Signature.sig2ray

Utility functions
=================

.. autosummary::
    :toctree: generated/

    showsig
    gidl
    frontline
    edgeout2
    edgeout

"""
import doctest
import numpy as np
#import scipy as sp
import scipy.linalg as la
import pdb
import h5py
import copy
import time
import pickle
import logging
import networkx as nx
import shapely.geometry as shg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylayers.gis.layout as layout
import pylayers.util.geomutil as geu
import pylayers.util.cone as cone
#import pylayers.util.graphutil as gph
import pylayers.util.pyutil as pyu
import pylayers.util.plotutil as plu
from pylayers.antprop.rays import Rays
from pylayers.util.project import *
import heapq
import shapely.geometry as sh
import shapely.ops as sho
#from numba import autojit


def plot_lines(ax, ob, color = []):

    from descartes.patch import PolygonPatch
    for ii,line in enumerate(ob):
        if color == []:
            if ii ==0 : 
                c ='g'
            elif ii == len(ob)-1:
                c ='r'
            else:
                c= 'k'
        else:
            c=color

        x, y = line.xy
        
        ax.plot(x, y, color=c, alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
    return ax
def plot_poly(ax, ob, color = []):
    
    from descartes.patch import PolygonPatch
    for ii,poly in enumerate(ob):
        
        pp = PolygonPatch(poly,alpha=0.3)
        ax.add_patch(pp)

    return ax



def showsig(L,s,tx=[],rx=[]):
    """ show signature

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
        if len(n)>1:
            edlist.append(n)
    gr=g.subgraph(edlist)
    dpos = {k:g.pos[k] for k in edlist}
    gr.pos=dpos
    return(gr)

def shLtmp(L):
    seg_connect = {x:L.Gs.edge[x].keys() for x in L.Gs.nodes() if x >0}

    dpts = {x[0]:(L.Gs.pos[x[1][0]],L.Gs.pos[x[1][1]]) for x in seg_connect.items() }
    L._shseg = {p[0]:sh.LineString(p[1]) for p in dpts.items()}


def valid(lsig,L,tahe=[]):
    """ 
    Check if a signature is valid .
    if a segmetn of a given signature is not in or touvches the polygon
    descibed by the 1st and last segment , the signature is not valid


    Parameters
    ----------

    lsig : list of tuple from run  | signatures
    L : layout
    tahe : 
        lensig , ta|he , x,y

    Returns
    -------

    inside : boolean
        is the signature valid ?


    """


    lensi = len(lsig)
    if lensi<=3:
        return True

    # DEBUG
    # if lensi == 4:
    #     if np.all(lsig == np.array([[ 5,  2, 67, 58],[ 2,  2,  3,  2]]).T):
    #         import ipdb
    #         ipdb.set_trace()

    # ensure compatibility with Signautre.run where
    # lsig is a list of tuple
    if isinstance(lsig,list):
        lsig = np.array([(i[0],len(i)) for i in lsig])

    pta = np.empty((2,lensi))
    phe = np.empty((2,lensi))

    seq = lsig[:,0]
    # upos = np.where(seq>0)[0]
    # uneg = np.where(seq<0)[0]
    
    # tahep = L.seg2pts(seq[upos])
    # tahen = np.array([L.Gs.pos[i] for i in seq[uneg]]).T
    # tahen = np.vstack((tahen,tahen))
    # tahe = np.empty((4,lensi))
    # tahe[:,upos]=tahep
    # try:
    #     tahe[:,uneg]=tahen
    # except:
    #     pass
    # pts = [k for i in seq for k in [L.Gs[i].keys()[0],L.Gs[i].keys()[1]]]
    # if tahe ==[]:
    # print 'run tahe\n',np.array(tahe)
    if tahe == []:
        pts = [L.Gs[i].keys() for i in seq]
        tahe = np.array([[L.Gs.pos[p[0]],L.Gs.pos[p[1]]] for p in pts])

        pta[:,0] = tahe[0,0,:]
        phe[:,0] = tahe[0,1,:]

        typ = lsig[:,1]
        mirror=[]
        # lines = [L._shseg[seq[0]]]
        for i in range(1,lensi):
            # pam = pa[:,i].reshape(2,1)
            # pbm = pb[:,i].reshape(2,1)
            pam = tahe[i,0,:].reshape(2,1)
            pbm = tahe[i,1,:].reshape(2,1)
            if typ[i] == 2: # R
                for m in mirror:
                    pam = geu.mirror(pam,pta[:,m],phe[:,m])
                    pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
                pta[:,i] = pam.reshape(2)
                phe[:,i] = pbm.reshape(2)
                mirror.append(i)

            elif typ[i] == 3 : # T
                for m in mirror:
                    pam = geu.mirror(pam,pta[:,m],phe[:,m])
                    pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
                pta[:,i] = pam.reshape(2)
                phe[:,i] = pbm.reshape(2)
            elif typ[i] == 1 : # D
                pta[:,i] = pam.reshape(2)
                phe[:,i] = pbm.reshape(2)

    else:

        tahe=np.array(tahe)
        pta = tahe[:,0,:]
        phe = tahe[:,1,:]




    # ### ONLY FOR TEST TO BE DELETED
    # pts = [L.Gs[i].keys() for i in seq]
    # tahetest = np.array([[L.Gs.pos[p[0]],L.Gs.pos[p[1]]] for p in pts])
    # ptat = np.empty((2,lensi))
    # phet = np.empty((2,lensi))
    # ptat[:,0] = tahetest[0,0,:]
    # phet[:,0] = tahetest[0,1,:]

    # typ = lsig[:,1]
    # mirror=[]
    # # lines = [L._shseg[seq[0]]]
    # for i in range(1,lensi):
    #     # pam = pa[:,i].reshape(2,1)
    #     # pbm = pb[:,i].reshape(2,1)
    #     pam = tahetest[i,0,:].reshape(2,1)
    #     pbm = tahetest[i,1,:].reshape(2,1)
    #     if typ[i] == 2: # R
    #         for m in mirror:
    #             pam = geu.mirror(pam,ptat[:,m],phet[:,m])
    #             pbm = geu.mirror(pbm,ptat[:,m],phet[:,m])
    #         ptat[:,i] = pam.reshape(2)
    #         phet[:,i] = pbm.reshape(2)
    #         mirror.append(i)

    #     elif typ[i] == 3 : # T
    #         for m in mirror:
    #             pam = geu.mirror(pam,ptat[:,m],phet[:,m])
    #             pbm = geu.mirror(pbm,ptat[:,m],phet[:,m])
    #         ptat[:,i] = pam.reshape(2)
    #         phet[:,i] = pbm.reshape(2)
    #     elif typ[i] == 1 : # D
    #         ptat[:,i] = pam.reshape(2)
    #         phet[:,i] = pbm.reshape(2)

    # tahetest = np.dstack((ptat.T,phet.T)).swapaxes(1,2)
    # if np.sum(tahe-tahetest) != 0:
    #     import ipdb
    #     ipdb.set_trace()
    


    # determine the 2 side of the polygon ( top/bottom = tahe[0]/tahe[-1])
    # vl and vr are 2 director vector lying on the polygon side.
    if not (geu.ccw(pta[:,0],phe[:,0],phe[:,-1]) ^
            geu.ccw(phe[:,0],phe[:,-1],pta[:,-1]) ):
        vl = ( pta[:,0],pta[:,-1])
        vr = ( phe[:,0],phe[:,-1])

        # twisted = True
        # lef = sh.LineString((pta[:,0],pta[:,-1]))
        # rig = sh.LineString((phe[:,0],phe[:,-1]))
    else:    
        vl = ( pta[:,0], phe[:,-1])
        vr = ( phe[:,0],pta[:,-1])
        # twisted = False
        # lef = sh.LineString((pta[:,0],phe[:,-1]))
        # rig = sh.LineString((pta[:,-1],phe[:,0]))
        
       


    # looking situation where Tail and head are not inside the polygon
    # => both tahe are left of vr and vl
    # =>   both tahe are right of vr and vl
    lta = geu.isleft(pta[:,1:-1],vl[0][:,None],vl[1][:,None])
    rta = geu.isleft(pta[:,1:-1],vr[0][:,None],vr[1][:,None])
    lhe =  geu.isleft(phe[:,1:-1],vl[0][:,None],vl[1][:,None])
    rhe = geu.isleft(phe[:,1:-1],vr[0][:,None],vr[1][:,None])

    out = (lta & lhe ) | (~rta & ~rhe)
    inside = ~out

    # debug
    # plt.ion()
    # plt.gcf()
    # plt.title(str(cond))
    # plot_lines(ax=plt.gca(),ob=lines)
    # # plot_lines(ax=plt.gca(),ob=[lef],color='g')
    # # plot_lines(ax=plt.gca(),ob=[rig],color='r')
    # # plt.scatter(pta[0,:],pta[1,:],marker='d',s=70,label='tail')
    # # plt.scatter(phe[0,:],phe[1,:],marker='s',s=70,label='head')
    # plu.displot(vl[0].reshape(2,1),vl[1].reshape(2,1),arrow=True)
    # plu.displot(vr[0].reshape(2,1),vr[1].reshape(2,1),arrow=True)
    # plt.legend()

    return np.all(inside)






class Signatures(PyLayers,dict):
    """ set of Signature given 2 Gt cycle (convex) indices

    Attributes
    ----------

    L : gis.Layout
    source : int
        source convex cycle
    target : int
        target convex cycle

    """

    def __init__(self,L,source,target,cutoff=3):
        """ object constructor

        Parameters
        ----------

        L : Layout
        source : int
            cycle number
        target : int
            cycle index
        cutoff : int
            limiting depth level in graph exploration (default 3)

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

    def check(self):

        OK = Signatures(self.L,self.target,self.source)
        KO = Signatures(self.L,self.target,self.source)
        for i in self:
        
            sigs = self[i]
            for s in range(len(sigs)/2):
                sig = sigs[2*s:2*s+2,:]

                ok = valid(sig.T,self.L)
                if ok :
                    try :
                        OK[i]=np.vstack((OK[i],sig))
                    except:
                        OK[i]=[]
                        OK[i]=sig
                        pass
                else :

                    try :
                        KO[i]=np.vstack((KO[i],sig))
                    except:
                        KO[i]=[]
                        KO[i]=sig
                        pass
        return OK,KO

    def saveh5(self):
        """ save signatures in hdf5 format
        """

        filename=pyu.getlong(self.filename+'.h5',pstruc['DIRSIG'])
        f=h5py.File(filename,'w')

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f.attrs['L']=self.L.filename
            f.attrs['source']=self.source
            f.attrs['target']=self.target
            f.attrs['cutoff']=self.cutoff
            for k in self.keys():
                f.create_dataset(str(k),shape=np.shape(self[k]),data=self[k])
            f.close()
        except:
            f.close()
            raise NameError('Signature: issue when writting h5py file')

    def loadh5(self,filename=[]):
        """ load signatures hdf5 format
        """
        if filename == []:
            _filename = self.filename
        else :
            _filename = filename

        filename=pyu.getlong(_filename+'.h5',pstruc['DIRSIG'])

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f=h5py.File(filename,'r')
            for k in f.keys():
                self.update({eval(k):f[k][:]})
            f.close()
        except:
            f.close()
            raise NameError('Signature: issue when reading h5py file')


        _fileL=pyu.getshort(filename).split('_')[0]+'.ini'
        self.L=layout.Layout(_fileL)
        try:
            self.L.dumpr()
        except:
            self.L.build()
            self.L.dumpw()


    def _saveh5(self,filenameh5,grpname):
        """ Save in hdf5 compliant with Links

        Parameters
        ----------

        filenameh5
        hrpname

        Notes
        -----
        """


        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        # if grpname == '':
        #     grpname = str(self.source) +'_'+str(self.target) +'_'+ str(self.cutoff)
        try:
            # file management
            fh5=h5py.File(filename,'a')
            if not grpname in fh5['sig'].keys():
                fh5['sig'].create_group(grpname)
            else :
                raise NameError('sig/'+grpname +'already exists in '+filenameh5)
            f=fh5['sig/'+grpname]

            # write data
            f.attrs['L']=self.L.filename
            f.attrs['source']=self.source
            f.attrs['target']=self.target
            f.attrs['cutoff']=self.cutoff
            for k in self.keys():
                f.create_dataset(str(k),shape=np.shape(self[k]),data=self[k])
            fh5.close()
        except:
            fh5.close()
            raise NameError('Signature: issue when writting h5py file')


    def _loadh5(self,filenameh5,grpname):
        """ load signatures in hdf5 format compliant with class Links

        Parameters
        ----------

        filenameh5 : string
            filename of the h5py file (from Links Class)
        grpname : string
            groupname of the h5py file (from Links Class)


        See Also
        --------

        pylayers.simul.links

        """

        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        # if grpname =='':
        #     grpname = str(self.source) +'_'+str(self.target) +'_'+ str(self.cutoff)

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            fh5=h5py.File(filename,'r')
            f=fh5['sig/'+grpname]
            for k in f.keys():
                self.update({eval(k):f[k][:]})
            Lname=f.attrs['L']
            fh5.close()
        except:
            fh5.close()
            raise NameError('Signature: issue when reading h5py file')
        self.L=layout.Layout(Lname)
        try:
            self.L.dumpr()
        except:
            self.L.build()
            self.L.dumpw()


    def save(self):
        """ save signatures
        """
        L=copy.deepcopy(self.L)
        del(self.L)
        filename=pyu.getlong(self.filename+'.h5',pstruc['DIRSIG'])
        with open(filename, 'wb') as handle:
          pickle.dump(self, handle)
        self.L=L

    def load(self,filename=[]):
        """ load signatures
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
        """ algorithm for signature determination

        Parameters
        ----------

        G : Graph
        source : tuple or int
        target : tuple or int
        cutoff : int

        See Also
        --------

        pylayers.antprop.signature.run3

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


    def short_propath(self,G,source,target=None,dout={},cutoff=None,weight=False):
        """ updated dijkstra
        """
        if source==target:
            return ({source:0}, {source:[source]})
        dist = {}  # dictionary of final distances
        paths = {source:[source]}  # dictionary of paths
        seen = {source:0}
        fringe=[] # use heapq with (distance,label) tuples
        heapq.heappush(fringe,(0,source))
        firstloop=True
        while fringe:
            if not firstloop:
                oldv = v
            (d,v)=heapq.heappop(fringe)

            if v in dist:
                continue # already searched this node.
            dist[v] = d
            if v == target:
                break
            #for ignore,w,edgedata in G.edges_iter(v,data=True):
            #is about 30% slower than the following
            if firstloop:
                edata = iter(G[v].items())
            else:
                try:
                    edata = iter(G[oldv][v]['output'].items())
                except:
                    break
            for w,edgedata in edata:
                if weight :
                    if not firstloop:
                        vw_dist = dist[v] + edgedata
                    else :
                        vw_dist = dist[v] #+ edgedata.get(weight,1) #<= proba should be add here
                else :
                    vw_dist = dist[v]
                if cutoff is not None:
                    if vw_dist>cutoff:
                        continue
                if w in dist:
                    if vw_dist < dist[w]:
                        raise ValueError('Contradictory paths found:',
                                         'negative weights?')
                elif w not in seen or vw_dist < seen[w]:
                    seen[w] = vw_dist
                    heapq.heappush(fringe,(vw_dist,w))
                    paths[w] = paths[v]+[w]
            firstloop=False


        if paths.has_key(target):
            if dout.has_key(len(paths[target])):
                dout[len(paths[target])].append([[p[0], len(p)] for p in paths[target]])
            else :
                dout[len(paths[target])]=[]
                dout[len(paths[target])].append([[p[0], len(p)] for p in paths[target]])

        return dout


    def propaths(self,G, source, target, cutoff=1,bt=False):
        """ seek all simple_path from source to target

        Parameters
        ----------

        G : networkx Graph Gi
        source : tuple
            interaction (node of Gi)
        target : tuple
            interaction (node of Gi)
        cutoff : int
        bt : bool
            allow backtrace (visite nodes already visited)

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

                elif (child not in visited) or (bt): # else visit other node
                    # only visit output nodes except if bt
                    #pdb.set_trace()
                    try:
                        dintpro = G[visited[-1]][child]['output']
                    except:
                        dintpro ={}

                    stack.append(iter(dintpro.keys()))
                    #stack.append(iter(G[visited[-1]][child]['output']))
                    visited.append(child)
                    # check if child (current segment) is an airwall
                    if child[0] in self.L.name['AIR']:
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

    def propaths3(self,Gi,source,target,cutoff=None):
        """ seek shortest path from source to target

        Parameters
        ----------

        Gi : graph of interactions
        source : source interaction
        target : target interaction
        cutoff : cutoff

        """

        level = 0
        nextlevel={source:Gi[source]}   # list of nodes to check at next level
        paths={source:[source]}         # paths dictionary  (paths to key from source)

        while nextlevel:
            thislevel = nextlevel
            nextlevel = {}
            for v in thislevel:
                for w in thislevel[v]:
                    # reach a node which is not in paths
                    if w not in paths:
                        paths[w]=paths[v]+[w]
                        nextlevel[w]= Gi[v][w]['output'].keys()
                    if w == target:
                        nstr = np.array(map(lambda x: x[0],paths[w]))
                        typ  = np.array(map(lambda x: len(w),paths[w]))
            level=level+1
            if (cutoff is not None and cutoff <= level):  break


    def propaths2(self,G, source, target,dout={}, cutoff=1,bt=False):
        """ seek all simple_path from source to target

        Parameters
        ----------

        G : networkx Graph Gi
        dout : dictionnary
            ouput dictionnary
        source : tuple
            interaction (node of Gi)
        target : tuple
            interaction (node of Gi)
        cutoff : int
        bt : bool
            allow backtrace (visite nodes already visited)

        Returns
        -------

        dout : dictionnary
            key : int
               number of interactions
            values : list

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
                    path = visited + [target]
                    try:
                        dout[len(path)].append([[p[0],len(p)] for p in path])
                    except:
                        dout[len(path)]=[]
                        dout[len(path)].append([[p[0],len(p)] for p in path])
                    #yield visited + [target] # output signature

                elif (child not in visited) or (bt): # else visit other node
                    # only visit output nodes except if bt
                    #pdb.set_trace()
                    try:
                        dintpro = G[visited[-1]][child]['output']
                    except:
                        dintpro ={}

                    stack.append(iter(dintpro.keys()))
                    #stack.append(iter(G[visited[-1]][child]['output']))
                    visited.append(child)
                    # check if child (current segment) is an airwall
                    # warning not efficient if many airwalls
                    if child[0] in self.L.name['AIR']:
                        lawp.append(1)
                    else:
                        lawp.append(0)



            else: #len(visited) == cutoff (visited list is too long)
                if child == target or target in children:
                    path = visited + [target]
                    try:
                        dout[len(path)].append([[p[0],len(p)] for p in path])
                    except:
                        #print "non existing : ",len(path)
                        dout[len(path)]=[]
                        dout[len(path)].append([[p[0],len(p)] for p in path])
                    #print visited + [target]
                    #yield visited + [target]

                stack.pop()
                visited.pop()
                try:
                    lawp.pop()
                except:
                    pass
        return dout



    # def propaths2015(self,G, source, target,dout={}, cutoff=1):
    #     """ seek all simple_path from source to target

    #     Parameters
    #     ----------

    #     G : networkx Graph Gi
    #     dout : dictionnary
    #         ouput dictionnary
    #     source : tuple
    #         interaction (node of Gi)
    #     target : tuple
    #         interaction (node of Gi)
    #     cutoff : int
    #     bt : bool
    #         allow backtrace (visite nodes already visited)

    #     Returns
    #     -------

    #     dout : dictionnary
    #         key : int
    #            number of interactions
    #         values : list of numpy array



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
    #     # lawp = list of airwall position in visited
    #     lawp = []

    #     # while the list of iterators is not void
    #     # import ipdb
    #     # ipdb.set_trace()
    #     while stack: #
    #         # children is the last iterator of stack

    #         children = stack[-1]
    #         # next child

    #         child = next(children, None)

    #         # update number of useful segments
    #         # if there is airwall in visited
    #         if child is None  : # if no more child
    #             stack.pop()   # remove last iterator
    #             visited.pop() # remove from visited list
    #             try:
    #                 lawp.pop()
    #             except:
    #                 pass

    #         elif (len(visited) < (cutoff + sum(lawp))):# if visited list length is less than cutoff
    #             if child == target:  # if child is the target point
    #                 #print visited + [target]
    #                 path = visited + [target]

    #                 try:
    #                     dout[len(path)][0]=np.vstack((dout[len(path)][0],np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='uint16')))
    #                 except:
    #                     dout[len(path)]=[np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='uint16')]

    #                 #yield visited + [target] # output signature

    #             elif (child not in visited): # else visit other node
    #                 # only visit output nodes except if bt
    #                 #pdb.set_trace()
    #                 try:
    #                     dintpro = G[visited[-1]][child]['output']
    #                 except:
    #                     dintpro ={}
    #                 stack.append(iter(dintpro.keys()))
    #                 #stack.append(iter(G[visited[-1]][child]['output']))
    #                 visited.append(child)
    #                 if child[0] in self.L.name['AIR']:
    #                     lawp.append(1)
    #                 else:
    #                     lawp.append(0)




    #         else: #len(visited) == cutoff (visited list is too long)
    #             if child == target or target in children:
    #                 path = visited + [target]

    #                 try:
    #                     dout[len(path)][0]=np.vstack((dout[len(path)][0],np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='uint16')))
    #                 except:
    #                     dout[len(path)]=[np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='uint16')]

    #                 #print visited + [target]
    #                 #yield visited + [target]

    #             stack.pop()
    #             visited.pop()
    #             try:
    #                 lawp.pop()
    #             except:
    #                 pass
    #     return dout


    def propaths2015_2(self,G, source, target,dout={},M={},Mmap=[], cutoff=1):
        """ seek all simple_path from source to target

        Parameters
        ----------

        G : networkx Graph Gi
        dout : dictionnary
            ouput dictionnary
        source : tuple
            interaction (node of Gi)
        target : tuple
            interaction (node of Gi)
        cutoff : int
        bt : bool
            allow backtrace (visite nodes already visited)

        Returns
        -------

        dout : dictionnary
            key : int
               number of interactions
            values : list of numpy array



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
                    path = visited + [target]
                    # M = np.zeros((1,NGs),dtype='bool')
                    out = [i[0] for i in G[visited[-1]][target]['output'].keys()]
                        

                    if Mmap !=[]:
                        M[Mmap[path[-2][0]],Mmap[path[-1][0]],Mmap[out]]=True
                    else: 
                        M[path[-2][0],path[-1][0],out]=True
                    try:
                        dout[len(path)]=np.vstack((dout[len(path)],np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='int16')))
                        # dnvi[len(path)].append([[i[0],len(i)] for i in G[visited[-1]][child]['output'].keys()])
                        # out = [i[0] for i in G[visited[-1]][child]['output'].keys()]
                        # M[path[-2][0],path[-1][0],out]=True
                        # dnvi[len(path)]=np.vstack((dnvi[len(path)],M))

                    except:
                        dout[len(path)]=np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='int16')
                        # dnvi[len(path)]=[[[i[0],len(i)] for i in G[visited[-1]][child]['output'].keys()]]
                        # dnvi[len(path)]=M

                    #yield visited + [target] # output signature

                elif (child not in visited): # else visit other node
                    # only visit output nodes except if bt
                    #pdb.set_trace()
                    try:
                        dintpro = G[visited[-1]][child]['output']
                    except:
                        dintpro ={}
                    stack.append(iter(dintpro.keys()))
                    #stack.append(iter(G[visited[-1]][child]['output']))
                    visited.append(child)
                    if child[0] in self.L.name['AIR']:
                        lawp.append(1)
                    else:
                        lawp.append(0)

            else: #len(visited) == cutoff (visited list is too long)
                # if child == (56, 8, 14):
                #     import ipdb
                #     ipdb.set_trace()
                # print list(children)
                # print ((child == target) or (target in children))

                CC=list(children)
                if ((child == target) or (target in CC)):
                    path = visited + [target]

                    # M = np.zeros((1,NGs),dtype='bool')
                    out = [i[0] for i in G[visited[-1]][target]['output'].keys()]
                    if Mmap != []:
                        M[Mmap[path[-2][0]],Mmap[path[-1][0]],Mmap[out]]=True
                    else :
                        M[path[-2][0],path[-1][0],out]=True
                    try:
                        dout[len(path)]=np.vstack((dout[len(path)],np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='int16')))
                        # dnvi[len(path)].append([[i[0],len(i)] for i in G[visited[-1]][child]['output'].keys()])
                        # dnvi[len(path)].append(np.unique([i[0] for i in G[visited[-1]][child]['output'].keys()]))
                        # M[:,[i[0] for i in G[visited[-1]][child]['output'].keys()]]=True
                        # dnvi[len(path)]=np.vstack((dnvi[len(path)],M))

                    except:
                        dout[len(path)]=np.array([[p[0],len(p)] for p in path],ndmin=3,dtype='int16')
                        # dnvi[len(path)]=[[[i[0],len(i)] for i in G[visited[-1]][child]['output'].keys()]]
                        # dnvi[len(path)]=[np.unique([i[0] for i in G[visited[-1]][child]['output'].keys()])]
                        # M[:,[i[0] for i in G[visited[-1]][child]['output'].keys()]]=True
                        # dnvi[len(path)]=M

                    #print visited + [target]
                    #yield visited + [target]

                stack.pop()
                visited.pop()
                try:
                    lawp.pop()
                except:
                    pass
        return dout

    def procone2(self,L,G, source, target,dout={}, cutoff=1):
        """ seek all simple_path from source to target looking backward

        Parameters
        ----------

        L : Layout
        G : networkx Graph Gi
        dout : dictionnary
            ouput dictionnary
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
                    #yield visited + [target] # output signature
                    path = visited + [target]
                    try:
                        dout[len(path)].append([[p[0],len(p)] for p in path])
                    except:
                        dout[len(path)]=[]
                        dout[len(path)].append([[p[0],len(p)] for p in path])
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
                        cn = cone.Cone()
                        ta,he = s.unfold()
                        segchild = np.vstack((ta[:,-1],he[:,-1])).T
                        segvm1 = np.vstack((ta[:,-2],he[:,-2])).T
                        cn.from2segs(segchild,segvm1)
                        typ,proba = cn.belong_seg(ta[:,:-2],he[:,:-2],proba=False)

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
                    #yield visited + [target]
                    path = visited + [target]
                    try:
                        dout[len(path)].append([[p[0],len(p)] for p in path])
                    except:
                        dout[len(path)]=[]
                        dout[len(path)].append([[p[0],len(p)] for p in path])
                stack.pop()
                visited.pop()
        return dout

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
        """ calculates signature

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

    def exist(self,seq):
        """ verifies if seq exists in signatures

        Parameters
        ----------

        seq : list or np.array()

        Returns
        -------

        boolean

        """
        if type(seq)==list:
            seq = np.array(seq)

        N = len(seq)
        sig = self[N]
        lf = filter(lambda x : (x==seq).all(),sig)
        if len(lf)>0:
            return True,lf
        else:
            return False


    def run_exp(self,source=-1,target=-1,cutoff=1,cutoffbound=1):
        """ EXPERIMENTAL
            Vectorized approach of signature search



        Parameters
        ----------

        source: int (-1)
            source cycle .
            If =-1 => self.source used
        target: int (-1)
            target cycle .
            If =-1 => self.target used
        cutoff= int (1)
            max number of interactions per cycle
            except 1st and last cycle
        cutoffbound= int (1)
            max number of interactions in 1st and last cycle

        Returns
        -------

        Nothing, fill self


        Example
        -------

        >>> from pylayers.simul.link import *
        >>> L=Layout('TA-Office.ini')
        >>> DL=DLink(L=L)
        >>> DL.ca=8
        >>> DL.Cb=13
        >>> DL.eval(force=['sig','ray','Ct','H'],alg=2015,si_reverb=3,cutoff=2,ra_vectorized=True)

        Notes
        -----

        This function makes use of graph Gc. Graph of meged cycles.

        """
        print "Run 2015"
        if source == -1:
            source = self.source
        if target == -1:
            target = self.target
        # list of cycle to reach source -> target. this will be imporve next
        lcil = self.L.cycleinline(source,target)
        llcil=len(lcil)

        # G=nx.Graph(self.L.Gt)

        # G.remove_node(0)
        # lcil = nx.dijkstra_path(G,source,target)
        # llcil=len(lcil)

        #2 determine input signatures for each cycles
        #di key = [input seg, input room, output seg, output room]
        di={}
        # number of points and seg of layout
        NGs = self.L.Ns+self.L.Np
        M = np.zeros((NGs,NGs,NGs),dtype='bool')

        ###
        ### Find interactions per cycles
        ###
        if llcil == 1:
            raise AttributeError('Signatures.run2015 doesn\'t work when source==target')
        for icy,cy in enumerate(lcil):

            vinT=[]
            #valid 'out' interatcion
            voutT=[]
            if self.L.Gt.node[cy].has_key('merged'):
                cym = self.L.Gt.node[cy]['merged']
                lcy = self.L.Gc.node[cym]['merged']
                inter=[]
                [inter.extend(self.L.Gt.node[x]['inter']) for x in lcy]
            else:
                lcy = cy
                inter = self.L.Gt.node[cy]['inter']
            sGi = nx.subgraph(self.L.Gi,inter)

            if icy == 0:
                # the interactions of 1st cycle are kept appart
                # di0 = {}

                outR,outT,outD = self.L.intercyGc2Gt(cy,typ='source')

                for cycle in lcil:
                    fcy = filter(lambda x: cycle == x[2],outT)
                    voutT.extend(fcy)
                vinT = outR + outD

                kdi0 = (0,0,0,voutT[0][0],voutT[0][1],voutT[0][2])

                #for each reverb/diffract interaction,
                # inside 1st cycle, search the output interactions
                
                for o in voutT:
                    io={}
                    for i in vinT:
                        io = self.propaths2015_2(sGi,i,o,dout=io,M=M,cutoff=cutoffbound)


                    di[0,0,0,o[0],o[1],o[2]] = io
                    # add direct signature
                    di[0,0,0,o[0],o[1],o[2]][1]=np.array([o[0],len(o)],ndmin=3)
                # import ipdb
                # ipdb.set_trace()
            elif (icy >=1) and (icy <llcil-1):

                #valid 'in' interatcion

                # select input signatures in regard of selected
                inR,inT,inD = self.L.intercyGc2Gt(cy,typ='target')
                outR,outT,outD = self.L.intercyGc2Gt(cy,typ='source')


                # keep only interactions in identified interesting cycles
                for cycle in lcil:
                    fcy = filter(lambda x: cycle == x[1],inT)
                    vinT.extend(fcy)
                    fcy = filter(lambda x: cycle == x[2],outT)
                    voutT.extend(fcy)

                # for each (identified interesting ) input interactions of the cycle
                #find all path to each (identified interesting) output interactions

                for i in vinT:
                    for o in voutT:
                        io={}
                        if not (i[1],i[2])==(o[2],o[1]):
                            io = self.propaths2015_2(sGi,i,o,dout=io,M=M,cutoff=cutoff)
                            di[i[0],i[1],i[2],o[0],o[1],o[2]] = io
                            # dni[i[0],i[1],i[2],o[0],o[1],o[2]] = ino

            # the interactions of last cycle are kept appart
            elif icy == llcil-1:

                inR,inT,inD = self.L.intercyGc2Gt(cy,typ='target')

                for cycle in lcil:
                    fcy = filter(lambda x: cycle == x[2],inT)
                    vinT.extend( fcy)
                voutT = inR + inD
                kdif = (vinT[0][0],vinT[0][1],vinT[0][2],0,0,0)
                #keep trace of last segments
                sinf = np.array([vinT[i][0] for i in range(len(vinT))])
                # for each (identified interesting ) input interactions,
                #find path to each reverb/diffract interaction of last cycle
                for i in vinT:
                    io={}
                    for o in voutT:
                        io=self.propaths2015_2(sGi,i,o,dout=io,M=M,cutoff=cutoffbound)
                    di[i[0],i[1],i[2],0,0,0] = io
                    # add direct signature
                    di[i[0],i[1],i[2],0,0,0][1]=np.array([i[0],len(i)],ndmin=3)

                # dni[i[0],i[1],i[2],0,0,0] = ino


        #dictionnary of interactions id keys
        # interaction id key are build as tuple Transmission inter in ,
        # Transmission interaction id input ,  Transmission interaction id output
        # e.g. (34, 13, 12, 36, 12, 11).
        #(0,0,0,X,X,X) stands for all interactions from the source
        #(X,X,X,0,0,0) stands for all interactions from the target
        kdi = di.keys()

        # Create 2 arrays with
        # input and output interactions id respectively
        adi0 = np.array(kdi0)
        adif = np.array(kdif)
        adi = np.array(di.keys())
        adii = adi[:,:3]
        adio = adi[:,3:]
        out=[]
        lsig={}

        #initialize loop on the 1st interaction id(0,0,0,X,X,X)

        # uinit = np.unique(np.where(adi[:,:3]==0)[0])
        uinit = np.where(np.sum(adi[:,:3]==0,axis=1)==3)[0]
        oldout=uinit
        stop=False
        dsigiosave={}
        dsigiosave.update({kdi[i][-3:]:di[kdi[i]] for i in uinit})

        def filldinda(d0,d1):
            for kd1 in d1:
                if d0.has_key(kd1):
                    d0[kd1]=np.vstack((d0[kd1],d1[kd1]))
                else:
                    d0[kd1]=d1[kd1]

        firstloop=True
        dsigio={}
        idx = 0

        while not stop:

            # for all detected valid output
            for k in oldout:

                us = np.where(-(adii-adio[k]).T.any(0))[0]
                keep=[]
                for iuus,uus in enumerate(us) :
                    bue = adi[uus][np.array([0,3])]==adi[out][:,np.array([0,3])]
                    ue =np.sum(bue,axis=1)
                    if len(np.where(ue==2)[0]) <=0:
                        keep.append(iuus)
                us = us[keep]

                out.extend(us.tolist())

                #     print ue
                #     if len(ue) == 0:
                #         keep.append(ue)
                #
                #     else:
                #         pass
                # import ipdb
                # ipdb.set_trace()
                # us=us[keep]
                for uus in us:
                    #1st input interactions to all identified a outputs
                    try:
                        lsig=dsigio[kdi[k][3:]]
                    except:
                        lsig = dsigiosave[kdi[k][3:]]
                    # lni=dni[kdi[k]]
                    sigio={}
                    # loop on input interactions
                    for ki in lsig.keys():
                        #loop on output interactions
                        for ko in di[kdi[uus]].keys():

                            # remove impossible signature in terms of cones
                            lin = len(lsig[ki])
                            lout = len(di[kdi[uus]][ko])
                            # manage case 1st interaction with no previous

                            if ki >1 and ko>1:
                                uso = lsig[ki][:,-2:,0]
                                uout = di[kdi[uus]][ko][:,1][:,0]
                                uvi = M[uso[:,0],uso[:,1],:][:,uout]

                                suvi=np.sum(uvi,axis=0)
                            else :
                                uvi = np.ones((lin,lout),dtype='bool')
                                suvi = lin*np.ones(lout)

                            for uv in range(lout):
                                ri = lsig[ki][uvi[:,uv]]
                                ro = np.tile(di[kdi[uus]][ko][uv],(suvi[uv],1,1))

                                # ri = np.repeat(lsig[ki][uvi[:,uv]],suvi[uv],axis=0)
                                # ro = np.tile(di[kdi[uus]][ko],(lin,1,1))

                                # uvi=uvi.reshape(lin*lout,order='F')

                                # ri=ri[uvi]
                                # ro=ro[uvi]

                                asig=np.hstack((ri,ro[:,1:]))

                                try:
                                    sigio[ki+ko-1]=np.vstack((sigio[ki+ko-1],asig))
                                except:
                                    sigio[ki+ko-1]=asig

                            # ri = np.repeat(lsig[ki],lout,axis=0)
                            # ro = np.tile(di[kdi[uus]][ko],(lin,1,1))

                            # uvi=uvi.reshape(lin*lout,order='F')

                            # ri=ri[uvi]
                            # ro=ro[uvi]

                            # asig=np.hstack((ri,ro[:,1:]))

                            # try:
                            #     sigio[ki+ko-1]=np.vstack((sigio[ki+ko-1],asig))
                            # except:
                            #     sigio[ki+ko-1]=asig
                    # key is the output segment
                    if dsigio.has_key(kdi[uus][-3:]):
                        filldinda(dsigio[kdi[uus][-3:]],sigio)
                    else:
                        dsigio[kdi[uus][-3:]]=sigio

            dsigiosave.update(dsigio)
            dsigio={}
            firstloop=False
            if not firstloop:
                if (adi[out][:,3:] == 0).all():
                    stop=True
                    break

            oldout=out
            out=[]
            idx = idx+1



            # # attempt to limit the combinatory
            survive1=adi[oldout][:,2]==lcil[idx]

            survive2 = adi[oldout][:,-1]==lcil[idx+1]
            # survive2 = np.ones((len(oldout)),dtype=bool)
            survive = np.where(survive1&survive2)[0]
            oldout=np.array(oldout)[survive].tolist()




        sig=dsigiosave[(0,0,0)]
        #reshaping to be compliant with signatures format
        sig2= {x:np.swapaxes(sig[x],1,2) for x in sig}
        sig2= {x:sig2[x].reshape(np.prod(sig2[x].shape[:2]),x) for x in sig2}
        self.update(sig2)
        #for debug
        return sig2


    def run_exp2(self,source=-1,target=-1,cutoff=1,cutoffbound=1):
        """ EXPERIMENTAL
            Vectorized approach of signature search



        Parameters
        ----------

        source: int (-1)
            source cycle .
            If =-1 => self.source used
        target: int (-1)
            target cycle .
            If =-1 => self.target used
        cutoff= int (1)
            max number of interactions per cycle
            except 1st and last cycle
        cutoffbound= int (1)
            max number of interactions in 1st and last cycle

        Returns
        -------

        Nothing, fill self


        Example
        -------

        >>> from pylayers.simul.link import *
        >>> L=Layout('TA-Office.ini')
        >>> DL=DLink(L=L)
        >>> DL.ca=8
        >>> DL.Cb=13
        >>> DL.eval(force=['sig','ray','Ct','H'],alg=2015,si_reverb=3,cutoff=2,ra_vectorized=True)

        """

        print "run2015_2"
        if source == -1:
            source = self.source
        if target == -1:
            target = self.target

        if source == target:
            raise AttributeError('Signatures.run2015 doesn\'t work when source==target')
        # # approach 1
        # # list of cycle to reach source -> target. this will be imporve next
        lcil = self.L.cycleinline(source,target)
        llcil=len(lcil)

        # #approach 2
        # #2.1- make a shortest path from source to target cycles
        # G=nx.Graph(self.L.Gt)
        # G.remove_node(0)
        # lcil = nx.dijkstra_path(G,source,target)
        # llcil=len(lcil)
        # nlcil=[]
        # # 2.2- for all cycle in the shortest path, find cycle intersected
        # for ll in range(llcil-1) :
        #     nlcil.extend(self.L.cycleinline(lcil[ll],lcil[ll+1])[:-1])
        # #add target exclued in the loop above
        # nlcil.extend([target])
        # lcil = nlcil
        llcil=len(lcil)

        #2 determine input signatures for each cycles
        #di key = [input seg, input room, output seg, output room]
        di={}
        # number of points and seg of layout
        allpt = np.hstack((self.L.tgs,self.L.ldiffin ))
        #mapping segemnts
        segmapp = self.L.tgs
        # mapping diffraction
        diffmapp = np.empty((-min(self.L.ldiffin)),dtype='int16')
        diffmapp[self.L.ldiffin]=np.arange(len(self.L.ldiffin)) + self.L.Ns+1
        # common mapping diff and segments
        mapp = np.hstack((segmapp,diffmapp))
        lapt = len(allpt)
        M = np.zeros((lapt,lapt,lapt),dtype='bool')

        ###
        ### Find interactions per cycles
        ###

        def dido(cy,lcy):
            """ Difraction In Diffraction Out
                determine, for merged cycles, which diffrxtion
                points get out / in of a cycle

            Parameters
            ----------
                cy : integer
                    cycle to investigate
                lcy : list
                    list of original cycles before merge ( in Gc)

            Return
            ------
                insideD: listr of diffraction inside cy
                inD: list of diffraction points in going in cy
                outD: list of diffraction points out going from cy
                ddin : a dcitionnary for naming ingoing diff poitns (points,cy in , cy out)
                dd : a dcitionnary for naming outgoing diff poitns (points,cy in , cy out)

            """
            if not isinstance(lcy,list):
                lcy=[lcy]
            outR,outT,D = self.L.intercyGc2Gt(cy,typ='source')
            # keep diff points
            D=np.unique(D)
            ddin={}
            ddout={}
            #outgoing diffraction
            outD=[]
            #ingoing diffraction
            inD=[]
            # inside cycle diffraction
            insideD=[]

            for d in D:
                # get cycles involved in diff point d
                dcy = self.L.ptGs2cy(d)
                # keep only current cycle and its merged neighbords in Gc
                dcy = filter(lambda x: x in lcy,dcy)
                #remove the current cycle
                dcy = filter(lambda x: x != cy,dcy)
                if len (dcy) > 0:
                    outD.append((d,))
                    inD.append((d,))
                    for ud in dcy:
                        try:
                            ddout[(d,)].append((d,cy,ud))
                            ddin[(d,)].append((d,ud,cy))
                        except:
                            ddout[(d,)]=[(d,cy,ud)]
                            ddin[(d,)]=[(d,ud,cy)]
                else :
                    insideD.append((d,))
            return insideD,inD,outD,ddin,ddout
        #####END dido

        for icy,cy in enumerate(lcil):

            vinT=[]
            #valid 'out' interatcion
            voutT=[]
            if self.L.Gt.node[cy].has_key('merged'):
                cym = self.L.Gt.node[cy]['merged']
                lcy = self.L.Gc.node[cym]['merged']
                inter=[]
                [inter.extend(self.L.Gt.node[x]['inter']) for x in lcy]
            else:
                lcy = cy
                inter = self.L.Gt.node[cy]['inter']
            sGi = nx.subgraph(self.L.Gi,inter)

            if icy == 0:
                # the interactions of 1st cycle are kept appart
                # di0 = {}

                insideR,outT,D = self.L.intercyGc2Gt(cy,typ='source')
                for cycle in lcil:
                    fcy = filter(lambda x: cycle == x[2],outT)
                    voutT.extend(fcy)
                insideD,inD,outD,ddin,ddout=dido(cy,lcy)


                #outgoing inter
                vout=voutT+outD
                #inside inter
                insideRD = insideR + insideD
                #for each reverb/diffract interaction,
                # inside 1st cycle, search the output interactions

                for o in vout:
                    io={}
                    for i in insideRD:
                        io = self.propaths2015_2(sGi,i,o,dout=io,M=M,Mmap=mapp,cutoff=cutoffbound)
                    if len(o) >1:
                        di[0,0,0,o[0],o[1],o[2]] = io
                        # add direct signature
                        di[0,0,0,o[0],o[1],o[2]][1]=np.array([o[0],len(o)],ndmin=3)
                    else :
                        for oo in ddout[o]:
                            di[0,0,0,oo[0],oo[1],oo[2]]=io
                            di[0,0,0,oo[0],oo[1],oo[2]][1]=np.array([o[0],len(o)],ndmin=3)

            elif (icy >=1) and (icy <llcil-1):

                #valid 'in' interatcion

                # select input signatures in regard of selected
                insideR,inT,D = self.L.intercy(cy,typ='target')
                insideR,outT,D = self.L.intercy(cy,typ='source')

                insideD,inD,outD,ddin,ddout=dido(cy,lcy)


                # keep only interactions in identified interesting cycles
                for cycle in lcil:
                    fcy = filter(lambda x: cycle == x[1],inT)
                    vinT.extend(fcy)
                    fcy = filter(lambda x: cycle == x[2],outT)
                    voutT.extend(fcy)

                #forbiden interactions

                dfi={}

                for ii in inD:
                    try:
                        dfi[ii].append(self.L.Gs[ii[0]].keys())
                    except:
                        dfi[ii]=self.L.Gs[ii[0]].keys()

                #incoming Transmission and diffraction
                vin = vinT + inD
                #outgoing Transmission and diffraction
                vout = voutT + outD
                # for each (identified interesting ) input interactions of the cycle
                #find all path to each (identified interesting) output interactions
                
                for i in vin:
                    for o in vout:
                        io={}
                        # no difraction
                        if len(i)>1 and len(o)>1:
                            # no backward
                            if (i[1],i[2])!=(o[2],o[1]) and (i[1],i[2])!=(o[1],o[2]):
                                io = self.propaths2015_2(sGi,i,o,dout=io,M=M,Mmap=mapp,cutoff=cutoff)
                                di[i[0],i[1],i[2],o[0],o[1],o[2]] = io
                        # input diffraction
                        elif len(i)<2 and len(o)>1:
                            if o[0] not in dfi[i]:
                                io={}
                                io = self.propaths2015_2(sGi,i,o,dout=io,M=M,Mmap=mapp,cutoff=cutoff)
                                for ii in ddin[i]:
                                    di[ii[0],ii[1],ii[2],o[0],o[1],o[2]]=io
                        # output diffraction
                        elif len(i)>1 and len(o)<2:
                            if i[0] not in dfi[o]:
                                io={}
                                io = self.propaths2015_2(sGi,i,o,dout=io,M=M,Mmap=mapp,cutoff=cutoff)
                                for oo in ddout[o]:
                                    di[i[0],i[1],i[2],oo[0],oo[1],oo[2]]=io
                        # input and output diffraction
                        else :
                            if (i)!=(o) :
                                io={}
                                io = self.propaths2015_2(sGi,i,o,dout=io,M=M,Mmap=mapp,cutoff=cutoff)
                                for ii in ddin[i]:
                                    for oo in ddout[o]:
                                        di[ii[0],ii[1],ii[2],oo[0],oo[1],oo[2]]=io
                                # dni[i[0],i[1],i[2],o[0],o[1],o[2]] = ino

            # the interactions of last cycle are kept appart
            elif icy == llcil-1:

                insideR,inT,D = self.L.intercyGc2Gt(cy,typ='target')
                
                insideD,inD,outD,ddin,ddout=dido(cy,lcy)


                for cycle in lcil:
                    # fcy = filter(lambda x: (cycle == x[2]) and (x[1] in lcil),inT)
                    fcy = filter(lambda x: cycle == x[2],inT)
                    vinT.extend( fcy)


                insideRD = insideR + insideD
                vin = vinT + inD

                # for each (identified interesting ) input interactions,
                #find path to each reverb/diffract interaction of last cycle
                for i in vin:
                    io={}
                    if len(i)> 1:
                        for o in insideRD:
                            io=self.propaths2015_2(sGi,i,o,dout=io,M=M,Mmap=mapp,cutoff=cutoffbound)
                        di[i[0],i[1],i[2],0,0,0] = io
                        # add direct signature
                        di[i[0],i[1],i[2],0,0,0][1]=np.array([i[0],len(i)],ndmin=3)
                    else : 
                        for o in insideRD:
                            io=self.propaths2015_2(sGi,i,o,dout=io,M=M,Mmap=mapp,cutoff=cutoffbound)
                        for ii in ddin[i]:
                            di[ii[0],ii[1],ii[2],0,0,0] = io
                            # add direct signature
                            di[ii[0],ii[1],ii[2],0,0,0][1]=np.array([i[0],len(i)],ndmin=3)
                # dni[i[0],i[1],i[2],0,0,0] = ino

        #dictionnary of interactions id keys
        # interaction id key are build as tuple Transmission inter in ,
        # Transmission interaction id input ,  Transmission interaction id output
        # e.g. (34, 13, 12, 36, 12, 11).
        #(0,0,0,X,X,X) stands for all interactions from the source
        #(X,X,X,0,0,0) stands for all interactions from the target
        kdi = di.keys()
        # Create 2 arrays with
        # input and output interactions id respectively
        adi = np.array(di.keys())
        adii = adi[:,:3]
        adio = adi[:,3:]
        out=[]
        lsig={}
        #initialize loop on the 1st interaction id(0,0,0,X,X,X)
        # uinit = np.unique(np.where(adi[:,:3]==0)[0])
        uinit = np.where(np.sum(adi[:,:3]==0,axis=1)==3)[0]
        oldout=uinit
        stop=False
        dsigiosave={}
        dsigiosave.update({kdi[i][-3:]:di[kdi[i]] for i in uinit})

        def filldinda(d0,d1):
            for kd1 in d1:
                if d0.has_key(kd1):
                    d0[kd1]=np.vstack((d0[kd1],d1[kd1]))
                else:
                    d0[kd1]=d1[kd1]

        def filldinda2(d0,d1):
            for kd1 in d1:
                if d0.has_key(kd1):
                    for kkd1 in d1[kd1]:
                        if d0[kd1].has_key(kkd1):
                            d0[kd1][kkd1]=np.vstack((d0[kd1][kkd1],d1[kd1][kkd1]))
                        else:
                            d0[kd1][kkd1]=d1[kd1][kkd1]
                else:
                    d0[kd1]=d1[kd1]

        lastloop=False
        dsigio={}
        idx = 0
        outall=[]

        while not stop:
            # for all detected valid output

            for k in oldout:

                #if not last loop, 
                #find all the output interaction with a given in interaction (adi[k])
                if not lastloop:
                    us = np.where(-(adii-adio[k]).T.any(0))[0]
                    keep=[]
                    # remove output interaction already visited
                    #for lightned computation
                    for iuus,uus in enumerate(us) :
                        if not uus in outall :
                            if adi[uus][-1] != 0 or lastloop:
                                keep.append(iuus)
                                outall.append(uus)

                    us = us[keep]
                    out.extend(us.tolist())

                # else output interaction is the inner interactions
                else :
                    us = [k]

                for uus in us:
                    if not lastloop:
                        lsig = dsigiosave[kdi[k][3:]]
                    else : 
                        lsig = dsigiosave[kdi[k][:3]]
                    sigio={}
                    # loop on input interactions
                    for ki in lsig.keys():
                        #loop on output interactions
                        for ko in di[kdi[uus]].keys():
                          # remove impossible signature in terms of cones
                            lin = len(lsig[ki])
                            lout = len(di[kdi[uus]][ko])
                            # manage case 1st interaction with no previous
                            if ki >1 and ko>1:
                                uso = lsig[ki][:,-2:,0]
                                uout = di[kdi[uus]][ko][:,1][:,0]

                                uvi = M[mapp[uso[:,0]],mapp[uso[:,1]],:][:,mapp[uout]]
                                suvi=np.sum(uvi,axis=0)
                            else :
                                uvi = np.ones((lin,lout),dtype='bool')
                                suvi = lin*np.ones(lout)
                            for uv in range(lout):
                                ri = lsig[ki][uvi[:,uv]]
                                ro = np.tile(di[kdi[uus]][ko][uv],(suvi[uv],1,1))
                                asig=np.hstack((ri,ro[:,1:]))

                                    # ri = lsig[ki][uvi[:,uv]]
                                    # ro = np.tile(di[kdi[uus]][ko][uv],(suvi[uv],1,1))

                                # ri = np.repeat(lsig[ki][uvi[:,uv]],suvi[uv],axis=0)
                                # ro = np.tile(di[kdi[uus]][ko],(lin,1,1))

                                # uvi=uvi.reshape(lin*lout,order='F')

                                # ri=ri[uvi]
                                # ro=ro[uvi]


                                try:
                                    sigio[ki+ko-1]=np.vstack((sigio[ki+ko-1],asig))
                                except:
                                    sigio[ki+ko-1]=asig

                    # key is the output segment
                    if dsigio.has_key(kdi[uus][-3:]):
                        filldinda(dsigio[kdi[uus][-3:]],sigio)
                    else:
                        dsigio[kdi[uus][-3:]]=sigio

            filldinda2(dsigiosave,dsigio)
            dsigio={}


            if lastloop:
                stop=True
                break

            oldout=out
            out=[]
            idx = idx+1

            if lcil[idx] == lcil[-1]:
                lastloop=True

            #filtering for avoiding computing extra 
            #interaction in/out couples
            bo = np.zeros((len(adi)),dtype='bool')
            tmp= np.zeros((len(adi)),dtype='bool')
            # input cycle interaction must be in the already visited list of cycles (lcil)
            for ii in lcil[:idx] :
                bo = bo | (adi[:,1] == ii)
            # output cycle interaction must be in the remaining list of cycles (lcil)
            for ii in lcil[:idx] :
                bo = bo & (adi[:,-1] != ii)
            # output cycle interaction must bethe next cycle in  list of cycles (lcil)
            # or code 0 for inner interaction of last cycle
            for ii in lcil[idx:] :
                if not lastloop:
                    tmp = tmp | (adi[:,-1] == ii)
                else:
                    tmp = tmp | (adi[:,-1] == 0)
            bo = bo & tmp
            #consider only insteraction from the current cycle
            if not lastloop:
                bo = bo & ((adi[:,2]==lcil[idx]) & (adi[:,4]==lcil[idx]))
            oldout=np.where(bo)[0].tolist()



        sig=dsigiosave[(0,0,0)]
        #reshaping to be compliant with signatures format
        sig2= {x:np.swapaxes(sig[x],1,2) for x in sig}
        sig2= {x:sig2[x].reshape(np.prod(sig2[x].shape[:2]),x) for x in sig2 if len(sig2[x])>0 }
        self.update(sig2)
        #for debug
        return sig2


    def run_old(self,cutoff=2,algo='old',bt=False,progress=False,diffraction=True):
        """ get signatures (in one list of arrays) between tx and rx

        Parameters
        ----------

        cutoff : int
            limit the exploration of all_simple_path
        algo: string
            'old' : call propaths2
            'new' : call procone2
        bt : bool
            backtrace (allow to visit already visited nodes in simple path algorithm)
        progress : bool
            display the time passed in the loop


        Returns
        -------

        sigslist :  numpy.ndarray

        See Also
        --------

        pylayers.simul.link.Dlink.eval
        pylayers.antprop.signature.Signatures.propath2
        pylayers.antprop.signature.Signatures.procone2

        """
        print "run old"
        self.cutoff   = cutoff
        self.filename = self.L.filename.split('.')[0] +'_' + str(self.source) +'_' + str(self.target) +'_' + str(self.cutoff) +'.sig'

        # list of interactions visible from source
        lisT,lisR,lisD = self.L.intercy(self.source,typ='source')
        if diffraction:
            lis  = lisT + lisR + lisD
        else:
            lis  = lisT + lisR

        # list of interactions visible from target
        litT,litR,litD = self.L.intercy(self.target,typ='target')

        if diffraction:
           lit  = litT + litR + litD
        else:
           lit  = litT + litR
        #print "source,lis :",self.source,lis
        #print "target,lit :",self.target,lit


        Gi = self.L.Gi
        Gi.pos = self.L.Gi.pos
        #
        # remove diffractions from Gi
        if not diffraction:
            Gi = gidl(Gi)

        # initialize dout dictionnary
        dout = {}

        # progresss stuff...
        lmax = len(lis)*len(lit)
        pe = 0
        tic = time.time()
        tic0 = tic
        # lis=lis+lit
        # lit=lis+lit
        #for interaction source  in list of source interactions
        for us,s in enumerate(lis):
            #for target interaction in list of target interactions
            #print "---> ",s

            for ut,t in enumerate(lit):
                #print "   ---> ",t
                # progress bar
                if progress :

                    ratio = np.round((((us)*len(lit)+ut)/(1.*lmax))*10 )
                    if ratio > pe:
                        pe = ratio
                        toc = time.time()
                        print '~%d ' % (ratio*10),
                        print '%',
                        print '%6.3f %6.3f' % (toc-tic, toc-tic0)
                        tic = toc

                # if source and target interaction are different
                # and R | T
                #if ((type(eval(s))==tuple) & (s != t)):
                if (s != t):
                    if algo=='new':
                        dout = self.procone2(self.L,Gi,dout=dout,source=s,target=t,cutoff=cutoff)
                    elif algo == 'old' :
                        dout = self.propaths2(Gi,source=s,target=t,dout=dout,cutoff=cutoff,bt=bt)
                    elif algo == 'dij':
                        dout = self.short_propath(Gi,source=s,target=t,dout=dout,cutoff=cutoff)
                        # dout = self.short_propath(Gi,source=t,target=s,dout=dout,cutoff=cutoff)
                else:
                    try:
                        if [s[0],len(s)] not in dout[1]:
                            dout[1].append([s[0],len(s)])
                    except:
                        dout[1]=[]
                        dout[1].append([s[0],len(s)])

        for k in dout.keys():
            adout = np.array((dout[k]))
            shad  = np.shape(adout)
            # manage the case of signatures with only 1 interaction
            if k == 1:
                adout = adout.reshape(shad[0],1,shad[1])
                shad = np.shape(adout)
            # rehape (rays * 2 , interaction)
            # the 2 dimensions come from the signature definition :
            # 1st row = segment index
            # 2nd row = type of interaction
            self[k] = adout.swapaxes(1,2).reshape(shad[0]*shad[2],shad[1])




    # @profile                    
    def run(self,cutoff=2,bt=False,progress=False,diffraction=True,threshold=0.1):
        """ get signatures (in one list of arrays) between tx and rx

        Parameters
        ----------

        cutoff : int
            limit the exploration of all_simple_path
        algo: string
            'old' : call propaths2
            'new' : call procone2
        bt : bool
            backtrace (allow to visit already visited nodes in simple path algorithm)
        progress : bool
            display the time passed in the loop


        Returns
        -------

        sigslist :  numpy.ndarray

        See Also
        --------

        pylayers.simul.link.Dlink.eval
        pylayers.antprop.signature.Signatures.propath2
        pylayers.antprop.signature.Signatures.procone2

        """
        print "run"
        self.cutoff   = cutoff
        self.filename = self.L.filename.split('.')[0] +'_' + str(self.source) +'_' + str(self.target) +'_' + str(self.cutoff) +'.sig'

        # list of interactions visible from source
        lisT,lisR,lisD = self.L.intercy(self.source,typ='source')
        if diffraction:
            lis  = lisT + lisR + lisD
        else:
            lis  = lisT + lisR

        # list of interactions visible from target
        litT,litR,litD = self.L.intercy(self.target,typ='target')

        if diffraction:
           lit  = litT + litR + litD
        else:
           lit  = litT + litR
        #print "source,lis :",self.source,lis
        #print "target,lit :",self.target,lit


        Gi = self.L.Gi
        Gi.pos = self.L.Gi.pos
        #
        # remove diffractions from Gi
        if not diffraction:
            Gi = gidl(Gi)

        # initialize dout dictionnary
        dout = {}

        # progresss stuff...
        lmax = len(lis)*len(lit)
        pe = 0
        tic = time.time()
        tic0 = tic
        #for interaction source  in list of source interactions
        for us,s in enumerate(lis):
            #for target interaction in list of target interactions
            #print "---> ",s

            for ut,t in enumerate(lit):
                #print "   ---> ",t
                # progress bar
                if progress :

                    ratio = np.round((((us)*len(lit)+ut)/(1.*lmax))*10 )
                    if ratio > pe:
                        pe = ratio
                        toc = time.time()
                        print '~%d ' % (ratio*10),
                        print '%',
                        print '%6.3f %6.3f' % (toc-tic, toc-tic0)
                        tic = toc

                # if source and target interaction are different
                # and R | T
                #if ((type(eval(s))==tuple) & (s != t)):
                pts = self.L.Gs[s[0]].keys()
                tahe = [np.array([self.L.Gs.pos[pts[0]],self.L.Gs.pos[pts[1]]])]
                # R is a list which contains reflextion matrices(Sn) and translation matrices(vn)
                # for mirror
                # R=[[S0,v0],[S1,v1],...]
                R=[(np.eye(2),np.array([0,0]))]
                if (s != t):

                    visited = [s]
                    # stack is a list of iterators
                    stack = [iter(Gi[s])]
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
                        if child is None  : # if no more child
                            stack.pop()   # remove last iterator
                            visited.pop() # remove from visited list
                            tahe.pop()
                            R.pop()
                            try:
                                lawp.pop()
                            except:
                                pass



                        elif (len(visited) < (cutoff + sum(lawp))) :# if visited list length is less than cutoff
                            if child == t:  # if child is the target point
                                #print visited + [target]
                                path = visited + [t]
                                nstr = np.array(map(lambda x: x[0],path))
                                typ  = np.array(map(lambda x: len(x),path))
                                try:
                                    self[len(typ)]=np.vstack((self[len(typ)],nstr,typ))
                                except:
                                    self[len(typ)]=np.vstack((nstr,typ))
                                #try:
                                #    dout[len(path)].append([[p[0],len(p)] for p in path])
                                #except:
                                #    dout[len(path)]=[]
                                #    dout[len(path)].append([[p[0],len(p)] for p in path])
                                #yield visited + [target] # output signature
                            elif (child not in visited) or (bt): # else visit other node
                                # only visit output nodes except if bt
                                #pdb.set_trace()
                                
                                try:
                                    nexti  = Gi[visited[-1]][child]['output'].keys()
                                    # keyprob  = Gi[visited[-1]][child]['output'].items()
                                    # nexti = map(lambda x:x[0]
                                    #               ,filter(lambda x
                                    #                       :x[1]>threshold,keyprob))

                                except:
                                    nexti = []


                                visited.append(child)
                                seg = visited[-1][0]

                                # diff
                                if len(visited[-2]) == 1:
                                    th = L.Gs.pos[seg]
                                    th = np.array([th,th])
                                    R.append((np.eye(2),np.array([0,0])))

                                # refl
                                if len(visited[-2])==2 and len(visited)> 2:

                                    pts = self.L.Gs[seg].keys()
                                    # th (xy,Npt)
                                    th = np.array([self.L.Gs.pos[pts[0]],self.L.Gs.pos[pts[1]]])
                                    R.append(geu.axmat(tahe[-1][0],tahe[-1][1]))
  
                                # trans
                                else : 
                                    pts = self.L.Gs[seg].keys()
                                    th = np.array([self.L.Gs.pos[pts[0]],self.L.Gs.pos[pts[1]]])
                                    R.append((np.eye(2),np.array([0,0])))

                                # apply symmetry
                                for r in R: 
                                    th = np.einsum('ki,ij->kj',th,r[0])+r[1]
                                tahe.append(th)

                                v = valid(visited,self.L,tahe) 
                                if v:
                                    stack.append(iter(nexti))
                                

                                    # import ipdb
                                    # ipdb.set_trace()


                                    # check if child (current segment) is an airwall
                                    # warning not efficient if many airwalls
                                    if child[0] in self.L.name['AIR']:
                                        lawp.append(1)
                                    else:
                                        lawp.append(0)
                                else:
                                    visited.pop()
                                    tahe.pop()
                                    R.pop()

                        else: #len(visited) == cutoff (visited list is too long)
                            if child == t or t in children:
                                path = visited + [t]
                                nstr = np.array(map(lambda x: x[0],path))
                                typ  = np.array(map(lambda x: len(x),path))
                                try:
                                    self[len(typ)]=np.vstack((self[len(path)],nstr,typ))
                                except:
                                    #print "non existing : ",len(path)
                                    self[len(typ)]=np.vstack((nstr,typ))
                                #print visited + [target]
                                #yield visited + [target]

                            stack.pop()
                            visited.pop()
                            tahe.pop()
                            R.pop()
                            try:
                                lawp.pop()
                            except:
                                pass

                else: # s==t
                    nstr = np.array([s[0]])
                    typ  = np.array([len(s)])
                    try:
                        self[1]=np.vstack((self[1],nstr,typ))
                    except:
                        #print "non existing : ",len(path)
                        self[1]=np.vstack((nstr,typ))




    def plot_cones(self,L,i=0,s=0,fig=[],ax=[],figsize=(10,10)):
        """ display cones of an unfolded signature

        Parameters
        ----------

        L : Layout
        i : int
            the interaction block
        s : int
            the signature number in the block
        fig :
        ax  :
        figsize :

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
            fig,ax = cn.show(fig = fig,ax = ax,figsize = figsize)

        return (fig,ax)


    def unfold(self,L,i=0,s=0):
        """ unfold a given signature

        return 2 np.ndarray of pta and phe "aligned"
        (reflexion interaction are mirrored)

        Parameters
        ----------

        L : Layout
        i : int
            the interaction block
        s : int
            the signature number in the block

        Returns
        -------

        pta,phe

        See Also
        --------

        Signature.unfold

        """

        si = Signature(self[i][(2*s):(2*s)+2])
        si.ev(L)
        pta,phe = si.unfold()
        
        
        return pta,phe

    def pltunfold(self,L,i=0,s=0):
        import shapely.ops as sho
        from descartes.patch import PolygonPatch
        plt.ion()
        plt.gcf()
        plt.clf()
        def plot_lines(ax, ob, color = []):
            for ii,line in enumerate(ob):
                if color == []:
                    if ii ==0 : 
                        c ='g'
                    elif ii == len(ob)-1:
                        c ='r'
                    else:
                        c= 'k'
                else:
                    c=color

                x, y = line.xy
                
                ax.plot(x, y, color=c, alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
            return ax
        def plot_poly(ax, ob, color = []):
            for ii,poly in enumerate(ob):
                
                pp = PolygonPatch(poly,alpha=0.3)
                ax.add_patch(pp)

            return ax
        pta,phe=self.unfold(L=L,i=i,s=s)

        ML =sh.MultiLineString([((pta[0][i],pta[1][i]),(phe[0][i],phe[1][i])) for i in range(pta.shape[1])])
        fig=plt.gcf()
        ax=plt.gca()
        ax = plot_lines(ax,ML)

        s0=sh.LineString([(pta[0,0],pta[1,0]),(phe[0,-1],phe[1,-1])])
        s1=sh.LineString([(phe[0,0],phe[1,0]),(pta[0,-1],pta[1,-1])])
        if s0.crosses(s1):
            s0=sh.LineString([(pta[0,0],pta[1,0]),(pta[0,-1],pta[1,-1])])
            s1=sh.LineString([(phe[0,0],phe[1,0]),(phe[0,-1],phe[1,-1])])
        cross = sh.MultiLineString([s0,s1,ML[0],ML[-1]])

        poly=sho.polygonize(cross)
        # ax = plot_lines(ax,cross,color='b')

        ax = plot_poly(ax,poly)

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
        color : string
        alphasig : float
        widthsig : float
        colsig : string
        ms : int
        ctx  : int
        crx :int



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

        # display layout
        fig,ax = L.showG(**kwargs)


        if kwargs['ctx']!=-1:
            Tpoly = self.L.Gt.node[kwargs['ctx']]['polyg']
            Tpoly.coul='r'
            Tpoly.plot(fig=fig,ax=ax,color='r')

        if kwargs['crx']!=-1:
            Rpoly = self.L.Gt.node[kwargs['crx']]['polyg']
            Rpoly.plot(fig=fig,ax=ax,color='g')

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
                sig = map(lambda x: self.L.Gs.pos[x],self[i][2*j])
                siga = np.array(sig)
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
        fig = plt.figure()

        nit = self.keys()
        ni = nit[uni]
        ust = len(self[ni])/2

        polyS = self.L.Gt.node[self.source]['polyg']
        cp1 = polyS.centroid.xy

        polyT = self.L.Gt.node[self.target]['polyg']
        cp2 = polyT.centroid.xy

        ptx = np.array([cp1[0][0],cp1[1][0]])
        prx = np.array([cp2[0][0],cp2[1][0]])

        st='a'

        while st != 'q':
            inter=[]
            ax = fig.add_subplot(111)
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
                    if ii == 3:
                        inter.append('D')
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

        ptx : numpy.array or int
            Tx coordinates is the center of gravity of the cycle number if
            type(tx)=int
        prx :  numpy.array or int
            Rx coordinates is the center of gravity of the cycle number if
            sigtype(rx)=int

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
        #
        # cycle on a line between 2 cycles
        # lc  = self.L.cycleinline(self.source,self.target)

        #
        # if source and target in the same merged cycle
        # and ptx != prx
        #
        los = shg.LineString(((ptx[0], ptx[1]), (prx[0], prx[1])))

        # convex cycle of each point
        cyptx = self.L.pt2cy(ptx)
        cyprx = self.L.pt2cy(prx)

        # merged cycle of each point
        polyctx = self.L.Gt.node[cyptx]['polyg']
        polycrx = self.L.Gt.node[cyprx]['polyg']

        dtxrx = np.sum((ptx-prx)*(ptx-prx))
        if dtxrx>1e-15:
            if cyptx==cyprx:
                if polyctx.contains(los):
                    rays.los = True
                else:
                    rays.los = False

        # k : Loop on interaction group
        #   l : loop on signature
        # --->
        #  this part should be a generator
        #
        for k in self:
            # print 'block#',k
            # if k ==3:
            #     import ipdb
            #     ipdb.set_trace()
            # get signature block with k interactions
            tsig = self[k]
            shsig = np.shape(tsig)
            for l in range(shsig[0]/2):
                sig = tsig[2*l:2*l+2,:]
                ns0 = sig[0,0]
                nse = sig[0,-1]
                validtx = True
                validrx = True

                if (ns0<0):
                    pD = self.L.Gs.pos[ns0]
                    TxD = shg.LineString(((ptx[0], ptx[1]), (pD[0], pD[1])))
                    seg = polyctx.intersection(TxD)
                    validtx = seg.almost_equals(TxD,decimal=4)
                    if not validtx:
                        pass
                        #print "Signature.rays": ns0

                if (nse<0):
                    pD = self.L.Gs.pos[nse]
                    DRx = shg.LineString(((pD[0], pD[1]), (prx[0], prx[1])))
                    validrx = polyctx.contains(DRx)
                    if not validrx:
                        pass
                        #print nse

                if validtx & validrx:
                    #    print sig
                    #    print pD
                    s  = Signature(sig)
                    #
                    # Transform signature into a ray
                    # --> sig2ray

                    isray,Yi  = s.sig2ray(self.L, ptx[:2], prx[:2])

                    if isray:
                        Yi = np.fliplr(Yi)
                        if k in rays.keys():
                            Yi3d = np.vstack((Yi[:, 1:-1], np.zeros((1, k))))
                            Yi3d = Yi3d.reshape(3, k, 1)
                            rays[k]['pt'] = np.dstack(( rays[k]['pt'], Yi3d))
                            rays[k]['sig'] = np.dstack(( rays[k]['sig'],
                                                        sig.reshape(2, k, 1)))
                        else:

                            rays[k] = {'pt': np.zeros((3, k, 1)),
                                       'sig': np.zeros((2, k, 1),dtype=int)}
                            rays[k]['pt'][0:2, :, 0] = Yi[:, 1:-1]
                            rays[k]['sig'][:, :, 0] = sig

        rays.nb_origin_sig = len(self)
        rays.origin_sig_name = self.filename
        return rays


    def raysv(self,ptx=0,prx=1):

        """ from signatures dict to 2D rays Vectorized version

        Parameters
        ----------

        ptx : numpy.array or int
            Tx coordinates is the center of gravity of the cycle number if
            type(tx)=int
        prx :  numpy.array or int
            Rx coordinates is the center of gravity of the cycle number if
            type(rx)=int

        Returns
        -------

        rays : Rays

        Notes
        -----

        This is a vectorized version of Signatures.rays.
        This implementation take advantage of the np.ndarray
        and calculate images and backtrace for block of signatures.
        A block of signature gather all signatures with the same number of interaction.

        For mathematical details see :

        @phdthesis{amiot:tel-00971809,
          TITLE = {{Design of simulation platform joigning site specific radio propagation and human mobility for localization applications}},
          AUTHOR = {Amiot, Nicolas},
          URL = {https://tel.archives-ouvertes.fr/tel-00971809},
          NUMBER = {2013REN1S125},
          SCHOOL = {{Universit{\'e} Rennes 1}},
          YEAR = {2013},
          MONTH = Dec,
          KEYWORDS = {Electromagnetic wave propagation simulation ; Human mobility simulation ; Wireless localization methods ; Position estimation methods in wireless networks ; Vectorized computation ; Ray-tracing ; Ultra wide band ; Simulateur de propagation {\'e}lectromagn{\'e}tique ; Simulateur de mobilit{\'e} humaine ; M{\'e}thodes de localisation sans fils ; M{\'e}thodes d'estimation de la position dans les r{\'e}seaux sans fils ; Calcul informatique vectoris{\'e} ; Outil de trac{\'e} de rayons ; Ultra large bande},
          TYPE = {Theses},
          HAL_ID = {tel-00971809},
          HAL_VERSION = {v1},
        }

        See Also
        --------

        Signatures.image
        Signatures.backtrace

        """
        if type(ptx)==int:
            ptx = np.array(self.L.Gt.pos[ptx])

        if type(prx)==int:
            prx = np.array(self.L.Gt.pos[prx])


        if len(ptx) == 2:
            ptx= np.r_[ptx,0.5]
        if len(ptx) == 2:
            prx= np.r_[prx,0.5]

        rays = Rays(ptx,prx)
        #
        # detect LOS situation
        #
        #
        # cycle on a line between 2 cycles
        # lc  = self.L.cycleinline(self.source,self.target)

        #
        # if source and target in the same merged cycle
        # and ptx != prx
        #
        los = shg.LineString(((ptx[0], ptx[1]), (prx[0], prx[1])))

        # convex cycle of each point
        cyptx = self.L.pt2cy(ptx)
        cyprx = self.L.pt2cy(prx)


        polyctx = self.L.Gt.node[cyptx]['polyg']
        polycrx = self.L.Gt.node[cyprx]['polyg']

        dtxrx = np.sum((ptx-prx)*(ptx-prx))
        if dtxrx>1e-15:
            if polyctx.contains(los):
                rays.los = True
            else:
                rays.los = False


        M = self.image2(ptx)
        R = self.backtrace(ptx,prx,M)
        rays.update(R)
        rays.nb_origin_sig = len(self)
        rays.origin_sig_name = self.filename
        return rays

    def backtrace(self, tx, rx, M):
        ''' Warning :
            This is an attempt to vectorize the backtrace process.
            Despite it has been tested on few cases with succes,
            this is quite new need to be validated !!!


            Parameters
            ----------

                tx : ndarray
                    position of tx (2,)
                rx : ndarray
                    position of tx (2,)
                M : dict
                    position of intermediate point from self.image()

            Return
            -------

                rayp : dict
                key = number_of_interactions
                value =ndarray positions of interactions for creating rays

            Notes
            -----
            dictionnary of intermediate coordinated :
            key = number_of_interactions
            value = nd array M with shape : (2,nb_signatures,nb_interactions)
            and 2 represent x and y coordinates


        '''


        if len(tx) > 2:
            tx = tx[:2]
        if len(rx) > 2:
            rx = rx[:2]

        rayp={}

        # loop on number of interactions
        for ninter in self.keys():
            signatures = copy.deepcopy(self[ninter])
            #get segment ids of signature with 4 interactions
            #get segment ids of signature with ninter interactions
            seg = self[ninter][::2]
            unegseg=np.where(seg<0)
            uninegseg,idx = np.unique(seg[unegseg],return_inverse=True)
            pneg = np.array([self.L.Gs.pos[x] for x in uninegseg])

            nsig = len(seg)

            # determine positions of points limiting the semgments
            #1 get index in L.tahe
            # 2 get associated position in L.pt


            utahe = self.L.tahe[:,self.L.tgs[seg]]
            # pt : (xycoord (2),pt indexes (2),nb_signatures,nb_interactions)
            pt = self.L.pt[:,utahe]
            

            ####WARNING BIG TRICK HERE :
            #### pa and pb are not set as the same value 
            #### to avoid a singular matrixnext.
            #### set pa =-pb has no incidence but avoid complex and vain code 
            #### modification for handling diffractions
            try:
                pt[:,0,unegseg[0],unegseg[1]]=pneg[idx].T
                pt[:,1,unegseg[0],unegseg[1]]=-pneg[idx].T
            except:
                pass
            # pt shape =
            # 0 : (x,y) coordinates x=0,y=1
            # 1 : 2 points (linking the semgnet) a=0,b=1
            #2 : nb of found signatures/segments
            # 3 : nb interaction
            #shape =
            # 0 : (x,y) coordinates x=0,y=1
            # 1 : 2 points (linking the semgnet) a=0,b=1
            #2 : nb of found signatures/segments
            # 3 : nb interaction
            # how to do this into a while loop
            p=rx

            # creating W matrix required in eq (2.70) thesis Nicolas AMIOT
            #Warning W is rolled after and becomes (nsig,4,4)
            W=np.zeros((4,4,nsig))
            I=np.eye(2)[:,:,np.newaxis]*np.ones((nsig))
            W[:2,:2,...] = I
            W[2:4,:2,...] = I

            # once rolled :
            # W (nsig,4,4)
            W = np.rollaxis(W,-1)


            kinter=ninter-1

            ptr = pt
            Mr = copy.deepcopy(M)

            epsilon = 1e-2
            rayp_i = np.zeros((3,nsig,ninter))
            # rayp_i[:2,:,-1]=rx[:,None]
            #backtrace process
            # if ninter == 6:
            #     print np.where(((signatures[:,0]==42) &(signatures[:,1]==-277) & (signatures[:,2]==135) & (signatures[:,3]==21) & (signatures[:,4]==46) & (signatures[:,5]==319)))
            #     import ipdb
            #     ipdb.set_trace()

            while kinter > -1:

                #Initilization, using the Tx position
                if kinter == ninter-1:
                    p_min_m = p[:,np.newaxis]-Mr[ninter][:,:,kinter]
                else :
                    p_min_m = pvalid[:].T-Mr[ninter][:,:,kinter]

                a_min_b = ptr[:,0,:,kinter]-ptr[:,1,:,kinter]

                # Creating W from  eq (2.71)
                # a_min_b <=> a_{Lh-l}-b_{Lh-l}
                # p_min_m <=> \tilde{p}_{Lh}-\tilde{b}_{Lh-l}
                # W (nsig,4,4)
                # p_min_m (2,nsig)
                # a_min_b (2,nsig)
                W[...,:2,2] = p_min_m.T
                W[...,2:,3] = a_min_b.T

                # create 2nd member from eq (2.72)
                if kinter == ninter-1:
                    y= np.concatenate((p[:,np.newaxis]*np.ones((nsig)),ptr[:,0,:,kinter]))
                else:
                    y= np.concatenate((pvalid.T,ptr[:,0,:,kinter]))

                # y once transposed :
                # y (nsig,4)
                y=y.T


                # search and remove point with singular matrix
                invalid_sig=np.where(abs(np.linalg.det(W))<1e-15)

                W = np.delete(W,invalid_sig,axis=0)
                y = np.delete(y,invalid_sig,axis=0)
                ptr = np.delete(ptr,invalid_sig,axis=2)
                Mr[ninter] = np.delete(Mr[ninter],invalid_sig,axis=1)
                rayp_i = np.delete(rayp_i,invalid_sig,axis=1)

                #remove signatures

                usig = np.repeat(invalid_sig[0],2)
                usig[::2]=usig[::2]*2
                usig[1::2]=usig[1::2]*2+1
                signatures = np.delete(signatures,usig,axis=0)
                # detect diffrac
                uD = signatures[1::2,kinter]==1
                uuD = np.where(signatures[1::2,kinter]==1)[0]


                psolved = np.linalg.solve(W,y)

                #valid ray is : 0 < \alpha < 1 and 0< \beta < 1

                # alpha
                uvalidA= psolved[:,2]>0.
                uvalidB= psolved[:,2]<1.
                #beta
                uvalidC= psolved[:,3] >= epsilon
                uvalidD= psolved[:,3] <=1.-epsilon
                valid = uvalidA & uvalidB & uvalidC & uvalidD
                # consider valid diffraction interactions
                valid = valid | uD
                uvalid = np.where(valid)[0]

                # re-add correct position of diffraction interations
                #indeed diffraction point should not been solved with linalg, 
                # but by setting pa=-pb, no singular matrix appear
                #and diffraction points can be re-add thereafter.
                psolved[uuD,:2] = ptr[:,0,uuD,kinter].T

                pvalid = psolved[uvalid,:2]






                # keep only valid rays for ptr and Mr
                Mr[ninter]=Mr[ninter][:,uvalid,:]
                ptr=ptr[:,:,uvalid,:]
                W = W[uvalid,:,:]


                # remove signatures
                usigv = np.repeat(uvalid,2)
                usigv[::2]=usigv[::2]*2
                usigv[1::2]=usigv[1::2]*2+1
                signatures = signatures[usigv,:]
                rayp_i[:2,uvalid,kinter] = pvalid.T
                rayp_i = rayp_i[:,uvalid,:]
                #if no more rays are valid , then quit block
                # (kinter <0 is the exit while condition)
                if len(uvalid) > 0 :
                    kinter=kinter-1
                else :
                    kinter = -2

            # rayp_i[:2,:,0]=tx[:,None]
            if len(uvalid) !=0:
                sir1=signatures[::2].T.reshape(ninter,len(usigv)/2)
                sir2=signatures[1::2].T.reshape(ninter,len(usigv)/2)
                sig = np.empty((2,ninter,len(usigv)/2))
                sig[0,:,:]=sir1
                sig[1,:,:]=sir2
                rayp_i=np.swapaxes(rayp_i,1,2)
                rayp.update({ninter:{'pt':rayp_i,'sig':sig.astype('int')}})
        return rayp


    def image2(self,tx):
        """ determine rays from images (second implementation)

        Parameters
        ----------

        tx : point

        
        """
        if len(tx) > 2:
            tx = tx[:2]
        dM={}
        # loop on number of interactions
        for ninter in self.keys():

            #get segment ids of signature with ninter interactions
            seg = self[ninter][::2]
            # seek for diffraction
            # negative index points are diffraction points
            unegseg = np.where(seg<0)
            uninegseg,idx = np.unique(seg[unegseg],return_inverse=True)
            pneg = np.array([self.L.Gs.pos[x] for x in uninegseg])
            nsig = len(seg)

            M = np.empty((2,nsig,ninter))
            # determine positions of points limiting the segments
            #1 get index in L.tahe
            # 2 get associated position in L.pt

            utahe = self.L.tahe[:,self.L.tgs[seg]]
            # pt : (xycoord (2),pt indexes (2),nb_signatures,nb_interactions)
            pt = self.L.pt[:,utahe]
            #
            # TODO Upgrading layout for handling slab offsets 
            #
            # uncomment those two lines when the numpy array L.norm and
            # L.offset exist
            #norm    = self.L.normal[:,utahe]
            #offset  = self.L.offset[:,utahe]
            # pt = pt + offset*norm

            try:
                pt[:,0,unegseg[0],unegseg[1]] = pneg[idx].T
                pt[:,1,unegseg[0],unegseg[1]] = pneg[idx].T
            except:
                pass
            # pt shape =
            # 0 : (x,y) coordinates x=0,y=1
            # 1 : 2 points (linking the segment) a=0,b=1
            #2 : nb of found signatures/segments
            # 3 : nb interactions

            ############
            #formula 2.61 -> 2.64 N.AMIOT PH.D thesis
            ############
            sx = pt[0,1,:,:]-pt[0,0,:,:]
            sy = pt[1,1,:,:]-pt[1,0,:,:]
            den = sx**2+sy**2
            # den = ((pt[0,0,:,:]-pt[0,1,:,:])**2+(pt[1,0,:,:]-pt[1,1,:,:])**2)
            # avoiding singularity (should not be possible)
            uz = np.where(den==0)
            den[uz] = 1.

            a = 1 - (2. / den) * (pt[1,0,:, :] - pt[1,1,:, :]) ** 2
            b= (2. / den) * (pt[0,1,:, :] - pt[0,0,:, :]) * (pt[1,0,:, :] - pt[1,1,:, :])
            c = (2. / den) * (pt[0,0,:, :] * (pt[1,0,:, :] - pt[1,1,:, :]) ** 2 +
                              pt[1,0,:, :] * (pt[1,0,:, :] - pt[1,1,:, :]) *
                             (pt[0,1,:, :] - pt[0,0,:, :]))
            d = (2. / den) * (pt[1,0,:, :] * (pt[0,1,:, :] - pt[0,0,:, :]) ** 2 +
                              pt[0,0,:, :] * (pt[1,0,:, :] - pt[1,1,:, :]) *
                             (pt[0,1,:, :] - pt[0,0,:, :]))
            # a = ((pt[0,0,:,:]-pt[0,1,:,:])**2-(pt[1,0,:,:]-pt[1,1,:,:])**2)
            # a=a/(1.*den)

            # b = 2*(pt[0,1,:,:]-pt[0,0,:,:])*(pt[1,1,:,:]-pt[1,0,:,:])
            # b=b/(1.*den)

            # c= 2*(pt[0,0,:,:]*(pt[1,0,:,:]-pt[1,1,:,:])**2+pt[1,0,:,:]*(pt[0,1,:,:]-pt[0,0,:,:])*(pt[1,0,:,:]-pt[1,1,:,:]))
            # c = c/(1.*den)

            # d= 2*(pt[0,0,:,:]*(pt[1,0,:,:]-pt[1,1,:,:])*(pt[0,1,:,:]-pt[0,0,:,:])+pt[1,0,:,:]*(pt[0,1,:,:]-pt[0,0,:,:])**2)
            # d= d/(1.*den)

            # K=np.array([[a,-b],[-b,-a]])
            K = np.array([[a,-b],[-b,-a]])

            # translation vector v (2.60)
            v =np.array(([c,d]))

            ityp = self[ninter][1::2]

            for n in xrange(ninter):
                #get segment ids of signature with ninter interactions
                uT = np.where(ityp[:,n]==3)[0]
                uR = np.where(ityp[:,n]==2)[0]
                uD = np.where(ityp[:,n]==1)[0]
                if n ==0:
                    p=tx[:,None]*np.ones((nsig))
                else :
                    p=M[:,:,n-1]
                #reflexion 0 (2.67)
                M[:,uR,n] = np.einsum('ijk,jk->ik',K[:,:,uR,n],p[:,uR])+v[:,uR,n]
                #transmission 0 (2.67)
                M[:,uT,n] = p[:,uT]
                M[:,uD,n] = pt[:,0,uD,n]

            # if ninter==6:
            #     print np.where(((seg[:,0]==42) & (seg[:,1]==-277) & (seg[:,2]==135) & (seg[:,3]==21)&(seg[:,-1]==319)))
            #     import ipdb
            #     ipdb.set_trace()

            dM.update({ninter:M})
        return dM
    def image(self,tx=np.array([2.7,12.5])):
        ''' Warning :
            This is an attempt to vectorize the image process.
            Despite it has been tested on few cases with succes,
            this is quite new need to be validated !!!


            Parameters
            ----------

                tx : ndarray
                    position of tx (2,)

            Return
            -------

                M : dictionnary

            dictionnary of intermediate coordinated :
            key = number_of_interactions
            value = nd array M with shape : (2,nb_signatures,nb_interactions)
            and 2 represent x and y coordinates


        '''
        if len(tx) > 2:
            tx = tx[:2]

        def nb_split(a):
            nsp = 2
            out=False
            while not out:
                res=a%nsp
                if res!=0:
                    nsp=nsp+1
                else:
                    out=True
            return nsp

        dM={}
        for ninter in self.keys():
            #get segment ids of signature with ninter interactions
            seg = self[ninter][::2]
            nsig = len(seg)
            # determine positions of points limiting the semgments
            #1 get index in L.tahe
            # 2 get associated position in L.pt

            #utahe (2 pt indexes,nb_signatures,nb_interactions)

            utahe = self.L.tahe[:,self.L.tgs[seg]]



            # pt : (xycoord (2),pt indexes (2),nb_signatures,nb_interactions)
            pt = self.L.pt[:,utahe]

            # pt shape =
            # 0 : (x,y) coordinates x=0,y=1
            # 1 : 2 points (linking the semgnet) a=0,b=1
            #2 : nb of found signatures/segments
            # 3 : nb interaction

            ############
            #formula 2.61 -> 2.64 N.AMIOT thesis
            ############
            den = ((pt[0,0,:,:]-pt[0,1,:,:])**2+(pt[1,0,:,:]-pt[1,1,:,:])**2)
            uz = np.where(den ==0)
            den[uz] = 1.

            a = 1 - (2. / den) * (pt[1,0,:, :] - pt[1,1,:, :]) ** 2

            b= (2. / den) * (pt[0,1,:, :] - pt[0,0,:, :]) * (pt[1,0,:, :] - pt[1,1,:, :])

            c = (2. / den) * (pt[0,0,:, :] * (pt[1,0,:, :] - pt[1,1,:, :]) ** 2 +
                                pt[1,0,:, :] * (pt[1,0,:, :] - pt[1,1,:, :]) *
                                (pt[0,1,:, :] - pt[0,0,:, :]))
            d = (2. / den) * (pt[1,0,:, :] * (pt[0,1,:, :] - pt[0,0,:, :]) ** 2 +
                                pt[0,0,:, :] * (pt[1,0,:, :] - pt[1,1,:, :]) *
                                (pt[0,1,:, :] - pt[0,0,:, :]))
            # den = ((pt[0,0,:,:]-pt[0,1,:,:])**2+(pt[1,0,:,:]-pt[1,1,:,:])**2)

            # a = ((pt[0,0,:,:]-pt[0,1,:,:])**2-(pt[1,0,:,:]-pt[1,1,:,:])**2)
            # a=a/(1.*den)

            # b = 2*(pt[0,1,:,:]-pt[0,0,:,:])*(pt[1,1,:,:]-pt[1,0,:,:])
            # b=b/(1.*den)

            # c= 2*(pt[0,0,:,:]*(pt[1,0,:,:]-pt[1,1,:,:])**2+pt[1,0,:,:]*(pt[0,1,:,:]-pt[0,0,:,:])*(pt[1,0,:,:]-pt[1,1,:,:]))
            # c = c/(1.*den)

            # d= 2*(pt[0,0,:,:]*(pt[1,0,:,:]-pt[1,1,:,:])*(pt[0,1,:,:]-pt[0,0,:,:])+pt[1,0,:,:]*(pt[0,1,:,:]-pt[0,0,:,:])**2)
            # d= d/(1.*den)

            #get segment ids of signature with ninter interactions
            ityp = self[ninter][1::2]
            uT = np.where(ityp[:,1:]==3)
            uR = np.where(ityp[:,1:]==2)
            uD=np.where(ityp[:,1:]==1)

            #create matrix AM which is used to create marix A from eq. 2.65
            AM = np.eye(2*ninter)[:,:,np.newaxis]*np.ones(nsig)

            # Reflexion MAtrix K (2.59)
            K=np.array([[a,-b],[-b,-a]])
            # translation vector v (2.60)
            v =np.array(([c,d]))

            ############
            #Create matrix A (2.66) which is fill by blocks
            ############



            blocks=np.zeros((2,2,nsig,ninter-1))

            # Reflexion block
            blocks[:,:,uR[0],uR[1]]=-K[:,:,uR[0],uR[1]+1]
            # Transmission block
            blocks[:,:,uT[0],uT[1]]=-np.eye(2)[:,:,np.newaxis]*np.ones((len(uT[0])))
            # Diff block
            blocks[:,:,uD[0],uD[1]]=0.

            # fill the AM mda on the diagonal below the mda diagonal....
            A=pyu.fill_block_diagMDA(AM,blocks,2,-1)

            # The 2nd member y is firslty completly fill, without taking into account that the 1st line differst from others.
            # 1. find which interaction and signature are R|T|D => create a masked array
            # 2. repeat is created because to each signature/interaction correspond a 2x1 column. Repeat allow to have the correct size to fill y
            # 3. fill the 1st line of y to take into consideration that difference.

            #y is the 2nd memeber from from (2.65) and will be filled following (2.67)
            y = np.zeros((2 * ninter,nsig))

            #######
            # Determine where y has to be filed with R|T|D
            #####
            #find the position where there is T|R|D. non continuous => need mask array
            uTf = np.where(ityp==3)
            uRf = np.where(ityp==2)
            uDf =np.where(ityp==1)

            #postiion in signature <=> 2 lines in y . need to repeat to get the correct size
            uRy2=np.repeat(uRf[0],2)
            uRy1=np.repeat(uRf[1],2)
            uRy1=2*uRy1
            uRy1[1::2]=uRy1[::2]+1

            uDy2=np.repeat(uDf[0],2)
            uDy1=np.repeat(uDf[1],2)
            uDy1=2*uDy1
            uDy1[1::2]=uDy1[::2]+1
            try:
                y[uRy1,uRy2]=v[:,uRf[0],uRf[1]].ravel(order='F')
            except:
                pass #print 'no R'
            try:
                pass
                #uT1mr = np.repeat(uT1m.mask,2,axis=1).T
                #nothing to do. shoould be a zero vector , already initialized by y
            except:
                pass #print 'no T'
            try:
                # NEVER TESTED !!!!!!!!!!!
                y[uDy1,uDy2]=a[uDf]
            except:
                print "signatures.image diffraction line 3672 Not yet tested !"

                pass #print 'no D'

            ######
            #FIRST LINE specific processing of (2.67)
            ######
            uT0 = np.where(ityp[:,0]==3)[0]
            uR0 = np.where(ityp[:,0]==2)[0]
            uD0 =np.where(ityp[:,0]==1)[0]

            #reflexion 0 (2.67)
            r0 = np.einsum('ijk,j->ik',K[:,:,uR0,0],tx)+v[:,uR0,0]
            #trnasmission 0 (2.67)
            t0 = tx[:,np.newaxis]*np.ones(len(uT0))
            #diff 0 (2.67)
            d0 = a[uD0,0]
            #first line
            y[0:2,uR0]=r0
            y[0:2,uT0]=t0
            y[0:2,uD0]=d0

            #reshape for compliant size with linalg
            A=np.rollaxis(A,-1)
            y=np.rollaxis(y,-1)

            leA = len(A)
            res=0
            #trick for memory usage
            if leA > 1e4:
                nsp = nb_split(leA)
                if nsp != leA:
                    lA=np.split(A,nsp)
                    ly=np.split(y,nsp)
                    del A
                    del y
                    print nsp
                    for s in range(nsp):

                        lm=np.linalg.solve(lA[s], ly[s])
                        try:
                            m = np.vstack((m,lm))
                        except:
                            m = lm
                    del lm
                    del lA
                    del ly
                else:
                    m = np.linalg.solve(A, y)
            else :
                m = np.linalg.solve(A, y)
            M=np.array((m[:,0::2],m[:,1::2]))

            dM.update({ninter:M})
        return dM


class Signature(object):
    """ class Signature

    Attributes
    ----------

    seq : list  of interaction point (edges (>0)  or vertices (<0) [int]
    typ : list of interaction type 1-R 2-T 3-D  [int]
    pa  : tail point of interaction segment (2xN) ndarray
    pb  : head point of interaction segment (2xN) ndarray
    pc  : center point of interaction segment (2xN) ndarray

    """
    def __init__(self, sig):
        """ object constructor

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
            return(len(l))

        def seginter(l):
            try:
                l = eval(l)
            except:
                pass
            return l[0]

        if type(sig)==np.ndarray:
            self.seq = sig[0, :]
            self.typ = sig[1, :]

        if type(sig)==list:
            self.seq = map(seginter,sig)
            self.typ = map(typinter,sig)

    def __repr__(self):
        s = ''
        s = s + str(self.seq) + '\n'
        s = s + str(self.typ) + '\n'
        return s

    def info(self):
        for k in self.__dict__.keys():
            print k, ':', self.__dict__[k]

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
        set to 0

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
        signature.

        At that stage coordinates of extremities (tx and rx) is
        not known yet

        members data

        pa  tail of segment  (2xN)
        pb  head of segment  (2xN)
        pc  the center of segment (2xN)

        norm normal to the segment if segment
        in case the interaction is a point the normal is undefined and then
        set to 0.

        """

        # TODO : use map and filter instead of for loop

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

        returns 2 np.ndarray of pta and phe "aligned"
        reflexion interactions are mirrored

        Returns
        -------

        pta : np.array
        phe : np.array

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

            if self.typ[i] == 2: # R
                for m in mirror:
                    pam = geu.mirror(pam,pta[:,m],phe[:,m])
                    pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
                pta[:,i] = pam.reshape(2)
                phe[:,i] = pbm.reshape(2)
                mirror.append(i)

            elif self.typ[i] == 3 : # T
                for m in mirror:
                    pam = geu.mirror(pam,pta[:,m],phe[:,m])
                    pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
                pta[:,i] = pam.reshape(2)
                phe[:,i] = pbm.reshape(2)

            elif self.typ[i] == 1 : # D
                pass
                # TODO not implemented yet

        return pta,phe

    def evtx(self, L, tx, rx):
        """ evaluate transmitter

        Parameters
        ----------

        L  : Layout
        tx : np.array (2xN)
        rx : np.array (2xM)

        DEPRECATED

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
        """ compute the tx's images with respect to the signature segments

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
        usig = np.nonzero(typ[1:] == 1)[0]
        if len(usig) > 0:
            blocks[usig, :, :] = np.zeros((2, 2))
        # detect transmission
        tsig = np.nonzero(typ[1:] == 3)[0]
        if len(tsig) > 0:
            #blocks[tsig, :, :] = np.zeros((2, 2))
            blocks[tsig, :, :] = -np.eye(2)
        # detect reflexion
        rsig = np.nonzero(typ[1:] == 2)[0]
        if len(rsig) > 0:
            blocks[rsig, :, :] = S[rsig + 1, :, :]

        A = pyu.fill_block_diag(A, blocks, 2, -1)

        y = np.zeros(2 * N)
        if typ[0] == 2:
            vc0 = np.array([c[0], d[0]])
            v0 = np.dot(-S[0, :, :], tx) + vc0
        if typ[0] == 3:
            v0 = tx
        if typ[0] == 1:
            v0 = pa[:, 0]

        y[0:2] = v0
        for i in range(len(typ[1:])):
            if typ[i + 1] == 2:
                y[2 * (i + 1):2 * (i + 1) + 2] = np.array([c[i + 1], d[i + 1]])
            if typ[i + 1] == 3:
                #y[2 * (i + 1):2 * (i + 1) + 2] = y[2*i:2*i+2]
                y[2 * (i + 1):2 * (i + 1) + 2] = np.array([0,0])
            if typ[i + 1] == 1:
                y[2 * (i + 1):2 * (i + 1) + 2] = pa[:, i + 1]

        x = la.solve(A, y)
        M = np.vstack((x[0::2], x[1::2]))

        return M

    def backtrace_old(self, tx, rx, M):
        """ backtrace given image, tx, and rx

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
            >>> L.dumpr()
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
            >>> fig,ax = L.showG('s',fig=fig,ax=ax)
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
        #import ipdb
        #pdb.set_trace()
        #import pdb

        pa = self.pa
        pb = self.pb
        typ = self.typ

        N = np.shape(pa)[1]
        I2 = np.eye(2)
        z0 = np.zeros((2, 1))

        pkm1 = rx.reshape(2, 1)
        Y = pkm1
        k = 0          # interaction counter
        beta = .5      # to enter into the loop
        isvalid = True # signature is asumed being valid by default
        epsilon = 1e-2
        # if tuple(self.seq) == ( 42, -277,  135,   21,   46,  319):
        #     import ipdb
        #     ipdb.set_trace()
        # while (((beta <= 1) & (beta >= 0)) & (k < N)):
        while (((beta <= 1-epsilon) & (beta >= epsilon)) & (k < N)):
            #if int(typ[k]) != 1: # not a diffraction (surprisingly it works)
            if int(typ[N-(k+1)]) != 1: # not a diffraction
                # Formula (25) of paper Eucap 2013
                l0 = np.hstack((I2, pkm1 - M[:, N - (k + 1)].reshape(2, 1), z0))
                l1 = np.hstack((I2, z0,
                                pa[:, N - (k + 1)].reshape(2, 1) -
                                pb[:, N - (k + 1)].reshape(2, 1)
                                ))
                # print pkm1
                # import ipdb
                # ipdb.set_trace()
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
                alpha = 0.5 # dummy necessary for the test below
                # fixing #210
                #Y = np.hstack((Y, pa[:, k].reshape((2, 1))))
                #pkm1 = pa[:, k].reshape((2, 1))
                Y = np.hstack((Y, pa[:, N-(k+1)].reshape((2, 1))))
                pkm1 = pa[:, N-(k+1)].reshape((2, 1))
            k = k + 1
        if ((k == N) & ((beta > 0) & (beta < 1)) & ((alpha > 0) & (alpha < 1))):
            Y = np.hstack((Y, tx.reshape(2, 1)))
            return isvalid,Y
        else:
            isvalid = False
            return isvalid,(k,alpha,beta)






if __name__ == "__main__":
    plt.ion()
    print "testing pylayers/antprop/signature.py"
    doctest.testmod()
    print "-------------------------------------"
