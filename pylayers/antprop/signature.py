#-*- coding:Utf-8 -*-
from __future__ import print_function
"""
.. currentmodule:: pylayers.antprop.signature

.. autosummary::
    :members:

"""
import os
import glob
import doctest
import numpy as np
#import scipy as sp
import scipy.linalg as la
import pdb
import h5py
import copy
import time
import pickle
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
from tqdm import tqdm
#from numba import autojit


def plot_lines(ax, ob, color = []):
    """ plot lines with colors

    Parameters
    ----------

    ax : matplotlib axis
    ob : list of lines
    color : list (optional)

    """

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
    """ plot polygon

    Parameters
    ----------

    ax :
    ob :

    """

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

   gr : A graph

   """

    edlist=[]
    pos={}
    for n in g.nodes():
        if len(n)>1:
            edlist.append(n)
    gr = g.subgraph(edlist)
    for k in gr.edges():
        #print(k)
        #di = gr[k[0]][k[1]]
        #ke = di['output'].keys()
        ke  = gr[k[0]][k[1]]['output']
        #ke = di['output'].keys()
        #va = di['output'].values()
        #keva = zip(ke,va)
        #keva_valid = [ x for x in keva if len(x[0])>1]
        ke_valid = [ x for x in ke if len(x)>1 ]
        #gr[k[0]][k[1]]['output'] = dict(keva_valid)
        gr[k[0]][k[1]]['output'] = ke_valid

    dpos = {k:g.pos[k] for k in edlist}
    gr.pos=dpos
    return(gr)

def shLtmp(L):
    seg_connect = {x:L.Gs.edge[x].keys() for x in L.Gs.nodes() if x >0}

    dpts = {x[0]:(L.Gs.pos[x[1][0]],L.Gs.pos[x[1][1]]) for x in seg_connect.items() }
    L._shseg = {p[0]:sh.LineString(p[1]) for p in dpts.items()}

def showsig2(lsig,L,tahe):
    if isinstance(lsig,list):
        lsig = np.array([(i[0],len(i)) for i in lsig])

    for k in lsig:
        k0 = k[0]
        k1 = k[1]
        if k0>0:
            npt = L.Gs[k0].keys()
            pta = np.array(L.Gs.pos[npt[0]])
            phe = np.array(L.Gs.pos[npt[1]])
            if k1==2:
                plu.displot(pta.reshape(2,1),phe.reshape(2,1),color='r',linewidth=2)
            if k1 ==3:
                plu.displot(pta.reshape(2,1),phe.reshape(2,1),color='g',linewidth=2)

    for th in tahe:
        ta = th[0]
        he = th[1]
        plu.displot(ta.reshape(2,1),he.reshape(2,1),color='k',linewidth=1)


    tahe = np.array(tahe) # Nseg x tahe x xy 
    pta = tahe[:,0,:].T  #2 x Nseg
    phe = tahe[:,1,:].T  # 2 x Nseg 
    seq = lsig[:,0]
    if not (geu.ccw(pta[:,0],phe[:,0],phe[:,-1]) ^
            geu.ccw(phe[:,0],phe[:,-1],pta[:,-1]) ):
        vr = ( pta[:,0],phe[:,-1])
        vl = ( phe[:,0],pta[:,-1])
        # twisted = True
        lef = sh.LineString((pta[:,0],phe[:,-1]))
        rig = sh.LineString((phe[:,0],pta[:,-1]))
    else:    
        vr = ( pta[:,0],pta[:,-1])
        vl = ( phe[:,0],phe[:,-1])
        lef = sh.LineString((pta[:,0],pta[:,-1]))
        rig = sh.LineString((phe[:,0],phe[:,-1]))
    plt.ion()
    plt.gcf()
    #L.showG('s',labels=True)
    lines = [L._shseg[seq[0]]]
    plt.title(str(lsig))
    plot_lines(ax=plt.gca(),ob=lines)
    plot_lines(ax=plt.gca(),ob=[lef],color='g')
    plot_lines(ax=plt.gca(),ob=[rig],color='r')
    plt.scatter(pta[0,:],pta[1,:],marker='d',s=70,label='tail')
    plt.scatter(phe[0,:],phe[1,:],marker='s',s=70,label='head')
    #plu.displot(vl[0].reshape(2,1),vl[1].reshape(2,1),arrow=True)
    #plu.displot(vr[0].reshape(2,1),vr[1].reshape(2,1),arrow=True)
    plt.axis('auto')
    plt.legend()
#@profile

def valid(lsig,L,tahe=[]):
    """

    Check if a signature is valid.
    if a segment of a given signature is not in or touches the polygon
    described by the 1st and last segment, the signature is not valid


    Parameters
    ----------

    lsig : list of tuple from run  |signatures
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

    # ensure compatibility with Signature.run where
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
    # if tahe == []:
    #     pts = [L.Gs[i].keys() for i in seq]
    #     tahe = np.array([[L.Gs.pos[p[0]],L.Gs.pos[p[1]]] for p in pts])

    #     pta[:,0] = tahe[0,0,:]
    #     phe[:,0] = tahe[0,1,:]

    #     typ = lsig[:,1]
    #     mirror=[]
    #     # lines = [L._shseg[seq[0]]]
    #     for i in range(1,lensi):
    #         # pam = pa[:,i].reshape(2,1)
    #         # pbm = pb[:,i].reshape(2,1)
    #         pam = tahe[i,0,:].reshape(2,1)
    #         pbm = tahe[i,1,:].reshape(2,1)
    #         if typ[i] == 2: # R
    #             for m in mirror:
    #                 pam = geu.mirror(pam,pta[:,m],phe[:,m])
    #                 pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
    #             pta[:,i] = pam.reshape(2)
    #             phe[:,i] = pbm.reshape(2)
    #             mirror.append(i)

    #         elif typ[i] == 3 : # T
    #             for m in mirror:
    #                 pam = geu.mirror(pam,pta[:,m],phe[:,m])
    #                 pbm = geu.mirror(pbm,pta[:,m],phe[:,m])
    #             pta[:,i] = pam.reshape(2)
    #             phe[:,i] = pbm.reshape(2)
    #         elif typ[i] == 1 : # D
    #             pta[:,i] = pam.reshape(2)
    #             phe[:,i] = pbm.reshape(2)

    # else:

    tahe = np.array(tahe) # Nseg x tahe x xy 
    pta = tahe[:,0,:].T  #2 x Nseg
    phe = tahe[:,1,:].T  # 2 x Nseg 




    # ### ONLY FOR TEST TO BE DELETED
    # pts = [L.Gs[i].keys() for i in seq]
    # tahetest = np.array([[L.Gs.pos[p[0]],L.Gs.pos[p[1]]] for p in pts])
    # ptat = np.empty((2,lensi))
    # phet = np.empty((2,lensi))
    # ptat[:,0] = tahetest[0,0,:]
    # phet[:,0] = tahetest[0,1,:]

    # typ = lsig[:,1]
    # mirror=[]
    #lines = [L._shseg[seq[0]]]
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
    #vl and vr are 2 director vector lying on the polygon side.
    if not (geu.ccw(pta[:,0],phe[:,0],phe[:,-1]) ^
            geu.ccw(phe[:,0],phe[:,-1],pta[:,-1]) ):
        vr = ( pta[:,0],pta[:,-1])
        vl = ( phe[:,0],phe[:,-1])
        # vr = ( pta[:,0],phe[:,-1])
        # vl = ( phe[:,0],pta[:,-1])
        # twisted = True
        #lef = sh.LineString((pta[:,0],pta[:,-1]))
        #rig = sh.LineString((phe[:,0],phe[:,-1]))
    else:
        vr = ( pta[:,0], phe[:,-1])
        vl = ( phe[:,0],pta[:,-1])
        # vr = ( pta[:,0],pta[:,-1])
        # vl = ( phe[:,0],phe[:,-1])
        # twisted = False
        #lef = sh.LineString((pta[:,0],phe[:,-1]))
        #rig = sh.LineString((pta[:,-1],phe[:,0]))




    # looking situation where Tail and head are not inside the polygon
    # => both tahe are left of vr and vl
    #=>   both tahe are right of vr and vl
    lta = geu.isleft(pta[:,1:-1],vl[0][:,None],vl[1][:,None])
    rta = geu.isleft(pta[:,1:-1],vr[0][:,None],vr[1][:,None])
    lhe = geu.isleft(phe[:,1:-1],vl[0][:,None],vl[1][:,None])
    rhe = geu.isleft(phe[:,1:-1],vr[0][:,None],vr[1][:,None])

    out = (lta & lhe ) | (~rta & ~rhe)
    inside = ~out

    # #debug
    # plt.ion()
    # plt.gcf()
    # #plt.title(str(cond))
    # #Ok plot_lines(ax=plt.gca(),ob=lines)
    # plot_lines(ax=plt.gca(),ob=[lef],color='g')
    # plot_lines(ax=plt.gca(),ob=[rig],color='r')
    # plt.scatter(pta[0,:],pta[1,:],marker='d',s=70,label='tail')
    # plt.scatter(phe[0,:],phe[1,:],marker='s',s=70,label='head')
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

    def __init__(self,L,source,target,cutoff=3,threshold = 0.6):
        """ object constructor

        Parameters
        ----------

        L : Layout
        dump : int
        source : int
            cycle number
        target : int
            cycle index
        cutoff : int
            limiting depth level in graph exploration (default 3)

        A signature ia a dict of arrays

        The array is an interleaving between nstr and type of interaction 
        typeInt = 1,2,3 (extremity,diffraction,reflexion,transmission)

        Si[1]
        np.array([5,2,19,2,26,2,72,2]) 

        """
        self.L = L
        self.dump = -1
        self.source = source
        self.target = target
        self.cutoff = cutoff
        self.threshold = threshold
        self.ratio = {}
        self.filename = self.L._filename.split('.')[0] +'_' + str(self.source) +'_' + str(self.target) +'_' + str(self.cutoff) +'.sig'

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
            size[k] = int(len(self[k])/2)
        s = s + 'from cycle : '+ str(self.source) + ' to cycle ' + str(self.target)+'\n'
        if self.dump==-1:
            ldump = self.keys()
        else:
            ldump = self.dump

        for k in ldump:
            s = s + str(k) + ' : ' + str(size[k]) + '\n'
            a = np.swapaxes(self[k].reshape(size[k],2,k),0,2)
            # nl x 2 x nsig

            for l in np.arange(a.shape[2]):
                for i in range(k):
                    if i==k-1:
                        s = s + '('+str(a[i,0,l])+','+str(a[i,1,l])+')'
                    else:
                        s = s + '('+str(a[i,0,l])+','+str(a[i,1,l])+'),'
                s = s+'\n'


        return(s)

    def __len__(self):
        nsig = 0
        for k in self:
            size = int(len(self[k])/2)
            nsig += size
        return(nsig)

    def compl(self,lint,L):
        """ completion from lint

        Parameters
        ----------

        lint : list
            list of interactions

        Examples
        --------

        >>> Si.compl([(6220,3),(6262,3),(6241,3)],DL.L)


        """
        # all group of interactions
        for k in self:
            if k > len(lint):
                Si = self[k]
                Ns,Nb = Si.shape
                # all signatures form a group of interactions
                for l in range(int(Ns/2)):
                    # all interactions
                    b1 = True
                    for i1,it in enumerate(lint):
                        if ((Si[2*l,i1] == it[0]) and
                           (Si[2*l+1,i1] == it[1])):
                            pass 
                        else:
                            b1 = False
                    if b1:
                        sig = Si[2*l:2*l+2,:]
                        sigi = self.sig2inter(L,sig)
                        #print(k,l,' :',sigi)
                        # all 

    def sig2inter(self,L,lsi=[]):
        ''' convert signature to corresponding list of interactions in Gi

            Paramaters:
            ----------

            L : Layout
            lsi : nd.array 
                signature (2xnb_sig,sig_length)

            Examples:
            ---------

            >>> lsi = DL.Si[3]
            >>> DL.Si.sig2inter(DL.L,lsi)

        '''

        assert L.isbuilt,  AttributeError('Layout is not built')
        assert len(lsi)%2==0,   AttributeError('Incorrect signature(s) shape')

        tlinter = []
        for uu in range(0,len(lsi),2):

            si = lsi[uu:uu+2,:]


            lsig = si.shape[1]
            linter = []

            for k in range(lsig):
                # nstr : seg or points
                nstr = si[0,k]
                typ  = si[1,k]
                # cycles connected to seg or point
                seg_cy = copy.deepcopy(L.Gs.node[nstr]['ncycles'])

                if k == 0:
                    cy0 = self.source
                    lcy0 =[cy0]

                if (typ==3) or (typ==2):
                    cy0 = list(set(seg_cy).intersection(set(lcy0)))[0]
                    cy1 = [x for x in seg_cy if x!= cy0 ][0]

                if k == (lsig -1):
                    cy1 = self.target

                if typ == 1:
                    inter = (nstr,)
                    lcy0 = L.Gs.node[nstr]['ncycles']
                elif typ == 2:
                    inter = (nstr,cy0)
                elif typ == 3:
                    inter = (nstr,cy0,cy1)
                    # changing cycle
                    lcy0 = [cy1]
                linter.append(inter)
            tlinter.append(linter)
        if len(lsi) == 2:
            tlinter=tlinter[0]

        return tlinter


    def sig2prob(self,L,lsi):
        """ get signatures probability

        Parameters
        ---------
        L : Layout
        lsi : nd.array 
            signature (2xnb_sig,sig_length)


        Returns
        -------

        tlproba : list  (nb_sig,sig_length-2)
            output proba of each triplet of interaction



        """


        slsi = lsi.shape[1]

        assert L.isbuilt,  AttributeError('Layout is not built')
        assert hasattr(L,'Gi'),  AttributeError('Layout has not Gi Graph')
        assert L.Gi.size != 0,  AttributeError('Gi Graph is empty')
        assert len(lsi)%2==0,   AttributeError('Incorrect signature(s) shape')
        assert slsi>=3, AttributeError('Proba available for signature with at least 3 interacitons')

        linter = self.sig2inter(L,lsi)
        if len(lsi) == 2:
            linter=[linter]
        tlproba = []
        for inter in linter:
            lproba = []
            for k in range(slsi-2):
                proba = L.Gi[inter[k]][inter[k+1]]['output'][inter[k+2]]
                lproba.append(proba)
            tlproba.append(lproba)

        return tlproba


    def num(self):
        """ determine the number of signatures
        """
        self.nsig = 0
        self.nint = 0
        for k in self:
            size = int(len(self[k])/2)
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
        print(self.__class__.__name__ + '\n' + '----------'+'\n')
        #s = s + str(self.__sizeof__())+'\n'
        for k in self:
            size[k] = int(len(self[k])/2)
        print('from cycle : '+ str(self.source) + ' to cycle ' + str(self.target)+'\n')
        pyu.printout('Reflection',pyu.BLUE)
        print('  ')
        pyu.printout('Transmission',pyu.GREEN)
        print('  ')
        pyu.printout('Diffraction',pyu.RED)
        print('  \n')
        for k in self:
            print(str(k) + ' : ' + str(size[k]))
            a = np.swapaxes(self[k].reshape(size[k],2,k),0,2)
            # nl x 2 x nsig
            for i in range(k):

                nstr=a[i,0,:]
                typ=a[i,1,:]
                print('[',)
                for n,t in zip(nstr,typ):
                    if t==1:
                        pyu.printout(str(n),pyu.BLUE)
                    if t==2:
                        pyu.printout(str(n),pyu.GREEN)
                    if t==3:
                        pyu.printout(str(n),pyu.RED)
                print(']')
            print('\n')
                # s = s + '   '+ str(a[i,0,:]) + '\n'

                # s = s + '   '+ str(a[i,1,:]) + '\n'

    def check(self):
        """ check signature

        Returns
        -------

        OK : np.array
        KO : np.array

        """

        OK = Signatures(self.L,self.target,self.source)
        KO = Signatures(self.L,self.target,self.source)
        for i in self:

            sigs = self[i]
            for s in range(int(len(sigs)/2)):
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
            f.attrs['L']=self.L._filename
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
            f.attrs['L']=self.L._filename
            f.attrs['source']=self.source
            f.attrs['target']=self.target
            f.attrs['cutoff']=self.cutoff
            f.attrs['threshold']=self.threshold
            f.create_group('ratio')
            f.create_group('sig')
            for k in self.keys():
                f['sig'].create_dataset(str(k),shape=np.shape(self[k]),data=self[k])
                f['ratio'].create_dataset(str(k),shape=np.shape(self.ratio[k]),data=self.ratio[k])
            fh5.close()
        except:
            fh5.close()
            raise NameError('Signature: issue when writting h5py file')


    def _loadh5(self,filenameh5,grpname,**kwargs):
        """ load signatures in hdf5 format compliant with class Links

        Parameters
        ----------

        filenameh5 : string
            filename of the h5py file (from Links Class)
        grpname : string
            groupname of the h5py file (from Links Class)

        kwargs 
            may contain a L: layout object
            if L =  [] the layout is loaded from the layout name stored
            into the h5 file
            if L = Layout the layout passed in arg is used

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

            # compliant with new h5 format:
            if 'sig' in f.keys():
                for k in f['sig'].keys():
                    self.update({eval(k):f['sig'][k][:]})
                    self.ratio.update({eval(k):f['ratio'][k][:]})
            # old h5 format
            else:
                for k in f.keys():
                    self.update({eval(k):f[k][:]})
            Lname=f.attrs['L']
            self.cutoff = f.attrs['cutoff']

            if 'threshold' in f.attrs.keys():
                self.threshold = f.attrs['threshold']
            # ensure backward compatibility
            else:
                # find threshold
                th = np.min([np.min(self.ratio[x]) 
                                 for x in self.ratio])
                self.threshold = th.round(decimals=2)
            fh5.close()
        except:
            fh5.close()
            raise NameError('Signature: issue when reading h5py file')

        if 'L' in kwargs:
            self.L = kwargs['L']
        else:
            self.L = layout.Layout(Lname)
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

        seq : list of tuple
               [(2,2),(5,3),(7,2)]

        1 : Diffraction
        2 : Reflexion
        3 : Diffraction

        Returns
        -------

        Examples
        --------

        >>> DL=DLink()
        >>> DL.eval()
        >>> seq = [(2,3)] # transmission through segment 2
        >>> DL.Si.exist(seq)

        """
        # Number of interactions
        N = len(seq)
        # signatures with N interaction
        sig = self[N]
        # Number signature with N interaction
        Nsig = int(sig.shape[0]/2)
        nstr = sig[::2,:]
        typ  = sig[1::2,:]
        # List of signat
        lsig = []
        for k in range(Nsig):
            lint = []
            for l in range(N):
                lint.append((nstr[k,l],typ[k,l]))
            lsig.append(lint)

        if seq in lsig:
            return True
        else:
            return False

    #@profile
    def run(self,**kwargs):
        """ evaluate signatures between cycle of tx and cycle of rx

        Parameters
        ----------

        cutoff : int
            limit the exploration of all_simple_path
        bt : boolean
            backtrace (allow to visit already visited nodes in simple path algorithm)
        progress : boolean
            display the time passed in the loop
        diffraction : boolean
            activate diffraction
        threshold : float
            for reducing calculation time
        animations :  boolean
        nD : int
            maximum number of diffraction
        nR : int
            maximum number of reflection
        nT : int
            maximum number of transmission


        See Also
        --------

        pylayers.simul.link.Dlink.eval

        """
        defaults = {'cutoff' : 2,
                    'threshold': 0.1,
                    'delay_excess_max_ns': 400,
                    'nD': 1,
                    'nR': 10,
                    'nT': 10,
                    'bt' : True,
                    'progress': True,
                    'diffraction' : True,
                    'animation' : False
                    }
        self.cpt = 0
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        self.cutoff = kwargs['cutoff']
        if 'threshold' not in kwargs:
            kwargs['threshold'] = self.threshold
        else:
            self.threshold=kwargs['threshold']
        nD = kwargs['nD']
        nT = kwargs['nT']
        nR = kwargs['nR']
        bt = kwargs['bt']
        progress = kwargs['progress']
        diffraction = kwargs['diffraction']
        animation = kwargs['animation']
        delay_excess_max_ns = kwargs['delay_excess_max_ns']
        dist_excess_max = delay_excess_max_ns*0.3


        self.filename = self.L._filename.split('.')[0] +'_' + str(self.source) +'_' + str(self.target) +'_' + str(self.cutoff) +'.sig'

        #
        # AIR : editable AIR separation
        # _AIR : constructed AIR separation
        #

        lair = self.L.name['AIR'] + self.L.name['_AIR']

        # list of interactions visible from source
        lisR, lisT, lisD = self.L.intercy(self.source,typ='source')

        if diffraction:
            lis  = lisT + lisR + lisD
        else:
            lis  = lisT + lisR

        # list of interactions visible from target
        litR, litT, litD = self.L.intercy(self.target,typ='target')

        if diffraction:
           lit  = litT + litR + litD
        else:
           lit  = litT + litR

        pt_source = np.array(self.L.Gt.node[self.source]['polyg'].centroid.coords.xy)
        pt_target = np.array(self.L.Gt.node[self.target]['polyg'].centroid.coords.xy)

        d_source_target = np.linalg.norm(pt_source - pt_target)

        #print("source,lis :",self.source,lis)
        #print("target,lit :",self.target,lit)
        # for u in lit:
        #     print u
        # print "-------------"

        Gi = self.L.Gi
        Gi.pos = self.L.Gi.pos
        #
        # remove diffractions from Gi
        #
        if not diffraction:
            Gi = gidl(Gi)

        # initialize dout dictionnary
        dout = {}

        # progresss stuff...
        lmax = len(lis)*len(lit)
        pe = 0
        tic = time.time()
        tic0 = tic
        #for interaction source in list of source interactions

        bvisu = False
        # signature counter
        cptsig = 0

        if animation:
            fig,ax = self.L.showG('s',aw=1)
            ax.plot(self.L.Gt.pos[self.source][0],self.L.Gt.pos[self.source][1],'ob')
            ax.plot(self.L.Gt.pos[self.target][0],self.L.Gt.pos[self.target][1],'or')
        #
        # Loop over all interactions seen from the source
        #
        # us : loop counter
        # s  : interaction tuple
        # s[0] : point (<0) or segment (>0)a
        # pts : list of neighbour nodes from s[0]
        # tahe : segment extremities or point coordinates (repeated twice)
        lhash = []

        if progress :
            pbar = tqdm(total=100, desc='Signatures')

        for us,s in enumerate(lis):
            if progress:
                pbar.update(100./(1.*len(lis)))

            # start from a segment
            if s[0] > 0:
                pts = list(dict(self.L.Gs[s[0]]).keys())
                tahe = [ np.array([ self.L.Gs.pos[pts[0]], self.L.Gs.pos[pts[1]]]) ]
            # start from a point
            else:
                tahe = [np.array([self.L.Gs.pos[s[0]], self.L.Gs.pos[s[0]]])]

            # R is a list which contains reflexion matrices (Sn) and translation matrices(vn)
            # for interaction mirroring
            # R=[[S0,v0],[S1,v1],...]

            R = [(np.eye(2),np.array([0,0]))]

            # initialize visited list sequence with the first intercation s

            visited = [s]

            # if
            #   s is in target interaction list
            # or
            #   arrival cycle is equal to target cycle
            # then stack a new signature in self[len(typ)]
            #
            # TODO : It concerns self[1] : only one interaction (i.e several single reflection or diffraction)
            #
            if (s in lit) or (s[-1] == self.target):
                anstr = np.array([ x[0] for x in visited ])
                typ = np.array([len(x) for x in visited ])

                assert(len(typ)==1)
                try:
                    self[len(typ)] = np.vstack((self[len(typ)], anstr, typ))
                    self.ratio[len(typ)] = np.append(self.ratio[len(typ)],1.)
                except:
                    self[len(typ)] = np.vstack((anstr, typ))
                    self.ratio[len(typ)] = np.array([1.])
                # update signature counter
                cptsig +=1

            # stack is a list of iterators
            #
            #
            stack = [iter(Gi[s])]
            # air walls do not intervene in the number of transmission (cutoff criteria)
            # lawp is the list of airwall position in visited sequence
            # handle the case of the first segment which can be an airwall
            #
            if len(s)==3:
                nseg = s[0]
                if ((self.L.Gs.node[nseg]['name']=='_AIR') or
                   (self.L.Gs.node[nseg]['name']=='AIR')):
                    lawp = [1]
                else:
                    lawp = [0]
            else:
                lawp = [0]
            # while the stack of iterators is not void
            cpt = 0
            while stack: #
                # iter_on_interactions is the last iterator in the stack
                iter_on_interactions = stack[-1]
                # next interaction child
                interaction = next(iter_on_interactions, None)
                #print visited
                #if ((visited ==[(6236,74,91),(-213,)]) and (interaction==(-1002,))):
                #    print(interaction)
                #    pdb.set_trace()
                #if (visited ==[(6236,74,91),(-213,),(6248,99,111)]):
                #if (visited ==[(6236,74,91),(-213,),(6248,99,111),(6287,111,118)]):
                    #pdb.set_trace()
                #    import ipdb
                # cond1 : there is no more interactions
                # continue if True
                cond1 = not(interaction is None)
                # cond2 : enable reverberation
                #     interaction has not been visited yet
                #     or 
                #     bt : True (allow reentrance) (unconditionnaly)
                # continue if True
                #cond2 = (interaction in visited) and bt (old)
                cond2 = not (interaction in visited) or bt
                # cond3 : test the cutoff condition not get to the limit
                # continue if True
                cond3 = not(len(visited) > (self.cutoff + sum(lawp)))

                uD = [ k for k in range(len(visited)) if len(visited[k])==1 ]
                uR = [ k for k in range(len(visited)) if len(visited[k])==2 ]
                uT = [ k for k in range(len(visited)) if len(visited[k])==3 ]

                if cond1:
                    condD = True
                    condR = True
                    condT = True
                    if ((len(interaction)==1) and (len(uD)==nD)):
                        condD = False
                    if ((len(interaction)==2) and (len(uR)==nR)):
                        condR = False
                    if ((len(interaction)==3) and (len(uT)==nT)):
                        condT = False

                #
                #  animation
                #
                if animation :
                    cpt = cpt+1
                    edge = zip(visited[:-1],visited[1:])
                    N = nx.draw_networkx_nodes(Gi,pos=Gi.pos,
                            nodelist=visited,labels={},
                            node_size=15,ax=ax,fig=fig)
                    E = nx.draw_networkx_edges(Gi,pos=Gi.pos,
                            edgelist=edge,labels={},width=0.1,
                            arrows=False,ax=ax,fig=fig)

                    plt.savefig('./figure/' +str(us) +'_' + str(cpt) +'.png')
                    try:
                        ax.collections.remove(N)
                    except:
                        pass
                    try:
                        ax.collections.remove(E)
                    except:
                        pass

                if (cond1 and cond2 and cond3):
                    if (condD and condR and condT):
                        visited.append(interaction)
                        logger.debug("{}".format(visited))
                        self.cpt+=1
                        #print(visited)
                        # [(44,2,7),(62,7,15),(21,15),(62,15,7),(44,7,2),(16,2)]
                        # if visited ==[(6236,74,91),(141,91)]:
                        #     import ipdb
                        #     ipdb.set_trace()


                        # update list of airwalls
                        if interaction[0] in lair:
                            lawp.append(1)
                        else:
                            lawp.append(0)

                        # update number of useful segments
                        # if there is airwall in visited
                        nstr = interaction[0]
                        # Testing the type of interaction at rank -2
                        # R is a list which contains a rotation matrix
                        # and a translation vector for doing the mirroring
                        # operation

                        # diffraction (retrieve a point)
                        if len(visited[-2]) == 1:
                            #th = self.L.Gs.pos[nstr]
                            R.append((np.eye(2),np.array([0,0])))
                        elif len(visited[-2])==2:
                            #
                            # l'avant dernier point est une reflection
                            #
                            nseg_points = list(dict(self.L.Gs[visited[-2][0]]).keys())
                            ta_seg = np.array(self.L.Gs.pos[nseg_points[0]])
                            he_seg = np.array(self.L.Gs.pos[nseg_points[1]])
                            #
                            # get reflection matrix from segment visited[-2]
                            #
                            R.append(geu.axmat(ta_seg,he_seg))
                            # direct order
                            #R.append(geu.axmat(tahe[-1][0],tahe[-1][1]))
                        # transmission do nothing
                        else :
                            pass
                        # current interaction is of segment type
                        if (nstr>0):
                            nseg_points = list(dict(self.L.Gs[nstr]).keys())
                            th = np.array([self.L.Gs.pos[nseg_points[0]],
                                           self.L.Gs.pos[nseg_points[1]]])
                        else:
                            th = self.L.Gs.pos[nstr]
                            th = np.array([th,th])

                        # current interaction is of point type (diffraction)
                        # apply current chain of symmetries
                        #
                        # th   is the current segment tail-head coordinates
                        # tahe is a list of well mirrored tail-head coordinates

                            #tahe.append(a)
                        #if ((visited[0]==(104,23,17)) and (visited[1]==(1,17))):
                        #    print("th (avant mirror)",th)
                        ik = 1
                        r = R[-ik]
                        #
                        # dtarget :  distance between th and target
                        #
                        pt_th = np.sum(th,axis=0)/2.
                        d_target = np.linalg.norm(pt_target - pt_th)

                        #
                        # mirroring th until the previous point
                        #

                        th_mirror = copy.copy(th)

                        while np.any(r[0] != np.eye(2)):
                            th_mirror = np.einsum('ki,ij->kj',th_mirror,r[0])+r[1]
                            ik = ik + 1
                            r  = R[-ik]

                        pt_mirror = np.sum(th_mirror,axis=0)/2.
                        d_source = np.linalg.norm(pt_source-pt_mirror)
                        d_excess = d_source + d_target - d_source_target

                        # if at least 2 interactions
                        # or previous point is a diffraction

                        if (len(tahe)<2) or (len(visited[-2])==1) or (len(visited[-1])==1):
                            ratio = 1.0
                            ratio2 = 1.0
                        else:
                            # Determine the origin of the cone
                            # either the transmitter (ilast =0)
                            # or the last diffraction point (ilast=udiff[-1] )
                            udiff = [ k for k in range(len(visited)) if len(visited[k])==1 ]
                            if udiff==[]:
                                ilast = 0
                            else:
                                ilast=udiff[-1]

                            #print(tahe)
                            pta0 = tahe[ilast][0]   # tail first segment  (last difraction)
                            phe0 = tahe[ilast][1]   # head first segment

                            #
                            # TODO : it would be better to replace pta_ and phe_ with the intersection
                            # of the previous cone with tahe[-1]
                            #

                            pta_ = tahe[-1][0]  # tail last segment
                            phe_ = tahe[-1][1]  # head last segment

                            #
                            # Calculates the left and right vector of the cone
                            #
                            #  vl left vector
                            #  vr right vector
                            #
                            #
                            # Detect situations of connected segments
                            #
                            # [(60, 2, 8), (61, 8, 11), (15, 11), (61, 11, 8), (60 ,8, 2), (44, 2, 7)]
                            # if visited == [(60, 2, 8), (61, 8, 11), (15, 11), (61, 11, 8), (60 ,8, 2), (44, 2, 7)]:
                            #     print '\n',visited
                            #     import ipdb
                            #     ipdb.set_trace()

                            connected = False
                            if (pta0==pta_).all():
                                apex = pta0
                                connected = True
                                v0 = phe0 - apex
                                v_ = phe_ - apex
                            elif (pta0==phe_).all():
                                apex = pta0
                                connected = True
                                v0 = phe0 - apex
                                v_ = pta_ - apex
                            elif (phe0==pta_).all():
                                apex = phe0
                                connected = True
                                v0 = pta0 - apex
                                v_ = phe_ - apex
                            elif (phe0==phe_).all():
                                apex = phe0
                                connected = True
                                v0 = pta0 - apex
                                v_ = pta_ - apex

                            if connected:
                                if ((np.linalg.norm(v0)==0) or (np.linalg.norm(v_)==0)):
                                    logger.debug("pta0 : %g,%g", pta0[0], pta0[1])
                                    logger.debug("pta_ : %g,%g", pta_[0], pta_[1])
                                    logger.debug("phe0 : %g,%g", phe0[0], phe0[1])
                                    logger.debug("phe_ : %g,%g", phe_[0], phe_[1])
                                    logger.debug("v0 : %g,%g", v0[0], v0[1])
                                    logger.debug("v_ : %g,%g", v_[0], v_[1])
                            #
                            # Does the cone is built from 2 connected segments or
                            # 2 unconnected segments
                            #
                            if not connected:
                                if not (geu.ccw(pta0,phe0,phe_) ^
                                        geu.ccw(phe0,phe_,pta_) ):
                                    vr = (pta0,phe_)
                                    vl = (phe0,pta_)
                                else:  # twisted case
                                    vr = (pta0,pta_)
                                    vl = (phe0,phe_)

                                # cone dot product
                                # print vr
                                # print vl
                                vr_n = (vr[1]-vr[0])/np.linalg.norm(vr[1]-vr[0])
                                vl_n = (vl[1]-vl[0])/np.linalg.norm(vl[1]-vl[0])


                                vrdotvl = np.dot(vr_n,vl_n)
                                # cone angle
                                angle_cone = np.arccos(np.maximum(np.minimum(vrdotvl,1.0),-1.0))
                                #angle_cone = np.arccos(vrdotvl)
                                # prepare lines and seg argument for intersection checking
                                if angle_cone!=0:
                                    linel = (vl[0], vl[1]-vl[0])
                                    liner = (vr[0], vr[1]-vr[0])
                                    # from origin mirrored segment to be tested
                                    seg   = (th_mirror[0], th_mirror[1])

                                    # apex calculation
                                    a0u = np.dot(pta0, vr_n)
                                    a0v = np.dot(pta0, vl_n)
                                    b0u = np.dot(phe0, vr_n)
                                    b0v = np.dot(phe0, vl_n)
                                    #import warnings
                                    #warnings.filterwarnings("error")
                                    try:
                                        kb  = ((b0v-a0v)-vrdotvl*(b0u-a0u))/(vrdotvl*vrdotvl-1)
                                    except:
                                        pdb.set_trace()
                                    apex = phe0 + kb*vl_n



                            else: # cone from connected segments

                                v0n  = v0/np.linalg.norm(v0)
                                try:
                                    v_n  = v_/np.linalg.norm(v_)
                                except:
                                    pdb.set_trace()

                                # import ipdb
                                # ipdb.set_trace()
                                sign = np.sign(np.cross(v_n,v0n))
                                if sign>0:
                                    vr_n = -v0n
                                    vl_n = v_n
                                else:
                                    vr_n = v_n
                                    vl_n = -v0n

                                vrdotvl = np.dot(vr_n,vl_n)
                                # cone angle
                                angle_cone = np.arccos(np.maximum(np.minimum(vrdotvl,1.0),-1.))


                            #
                            # the illuminating cone is defined
                            # the th_mirror to be tested with this cone are known
                            #

                            if ( (not np.isclose(angle_cone,0,atol=1e-6) )
                             and ( not np.isclose(angle_cone,np.pi)) ) :
                                #if self.cpt==16176:
                                #    pdb.set_trace()
                                seg,ratio2 = geu.intersect_cone_seg((apex,vl_n),(apex,vr_n),(th_mirror[0],th_mirror[1]),bvis=False)
                            elif ( not np.isclose(angle_cone,0) ):
                                ratio2 = 1
                            else:
                                ratio2 = 0
                            #print ratio
                            if len(seg)==2:
                                th_mirror = np.vstack((seg[0],seg[1]))
                            else:
                                pass

                            al = np.arctan2(vl_n[1],vl_n[0])
                            ar = np.arctan2(vr_n[1],vr_n[0])
                            if np.allclose(th_mirror[0],apex) or np.allclose(th_mirror[1],apex):
                                ratio2 = 1.

                            # On connecte l'apex du cone courant aux extrémités du segment courant mirroré

                            # Dans certaines circonstances par example un cone emanant d'un point colinéaire 
                            # avec le segment d'arrivé" (-4) (6,4) le point -4 est aligné avec le segment 6
                            # l'ouverture du cone est nul => arret. Cela pourrait être géré dans Gi en interdisant 
                            # la visibilité (-4) (6,4) 

#                            if angle_cone ==0:
#                                ratio = 0
#                            else:
#                                if np.allclose(th_mirror[0],apex) or np.allclose(th_mirror[1],apex):
#                                    ratio = 1.
#                                else:
#                                    wseg0 = th_mirror[0] - apex
#                                    wseg1 = th_mirror[1] - apex
#                                    mod_wseg0 = np.sqrt(np.sum(wseg0*wseg0,axis=0))
#                                    mod_wseg1 = np.sqrt(np.sum(wseg1*wseg1,axis=0))
#
#                                    if np.isclose(mod_wseg0,0):
#                                        #bvisu = True 
#                                        #pdb.set_trace()#
#                                        pass
#                                    if np.isclose(mod_wseg1,0):
#                                        #bvisu = True 
#                                        #pdb.set_trace()#
#                                        pass
#                                    #wseg0_n = wseg0/mod_wseg0
#                                    #wseg1_n = wseg1/mod_wseg1
#                                    wseg0_n = wseg0/np.linalg.norm(wseg0)
#                                    wseg1_n = wseg1/np.linalg.norm(wseg1)
#                                    aseg0 = np.arctan2(wseg0_n[1],wseg0_n[0])
#                                    aseg1 = np.arctan2(wseg1_n[1],wseg1_n[0])
#                                    
#                                    # if al==aseg0 or al==aseg1 or ar==aseg0 or ar==aseg1:
#                                    #     ratio = 1
#                                        #print "toto"
#                                    # else:
#                                    I = geu.angle_intersection2(al,ar,aseg0,aseg1)
#                                    ratio = I/angle_cone
#                                    #if ratio>=1:
#                                    #    pdb.set_trace()
#
#                                # if connected:
#                                #     print "ratio :",ratio
#                                
#
#                            #if visited == [(104, 23, 17), (1, 17), (53, 17)]:
#                            if (bvisu):
#                                fig ,ax = self.L.showG('s',aw=1,labels=0)
#                                #
#                                # magenta : start of the cone
#                                # cyan    :
#                                # yellow  : last interaction
#                                #
#                                ax = geu.linet(ax,pta0,phe0,al=1,color='magenta',linewidth=3)
#                                ax = geu.linet(ax,pta_,phe_,al=1,color='cyan',linewidth=3)
#                                ax = geu.linet(ax,np.array(self.L.Gs.pos[nseg_points[0]]),np.array(self.L.Gs.pos[nseg_points[1]]),al=1,color='yellow',linewidth=4)
#                                # ax = geu.linet(ax,vr[0],vr[1],al=1,color='red',linewidth=3)
#                                # ax = geu.linet(ax,vl[0],vl[1],al=1,color='blue',linewidth=3)
#                                ax = geu.linet(ax,seg[0],seg[1],al=1,color='k',linewidth=3)
#                                ax = geu.linet(ax,th_mirror[0,:],th_mirror[1,:],al=1,color='green',linewidth=3)
#                                nx.draw_networkx_labels(self.L.Gi,
#                                        self.L.Gi.pos,labels={x:str(x) for x in visited},
#                                        ax=ax,fontsize=18)
#                                plt.title(str(visited)+'  '+str(ratio))
#                                ax.plot(apex[0],apex[1],'or')
#                                plt.axis('auto')
#                                pdb.set_trace()
#                            #if visited == [(104, 23, 17), (1, 17), (53, 17), (108, 17, 18)]:
#                            # if visited == [(104, 23, 17), (1, 17), (53, 17)]:
#                            if (1==0):
#                                fig ,ax = self.L.showG('s',aw=1,labels=0)
#                                ax = geu.linet(ax,pta0,phe0,al=1,color='magenta',linewidth=3)
#                                ax = geu.linet(ax,pta_,phe_,al=1,color='cyan',linewidth=3)
#
#                                ax = geu.linet(ax,np.array(self.L.Gs.pos[pts[0]]),np.array(self.L.Gs.pos[pts[1]]),al=1,color='yellow',linewidth=4)
#                                ax = geu.linet(ax,vr[0],vr[1],al=1,color='red',linewidth=3)
#                                ax = geu.linet(ax,vl[0],vl[1],al=1,color='blue',linewidth=3)
#                                #ax = geu.linet(ax,seg[0],seg[1],al=1,color='k',linewidth=3)
#                                ax = geu.linet(ax,th[0,:],th[1,:],al=1,color='green',linewidth=3)
#                                plt.title(str(visited)+'  '+str(ratio))
#                                ax.plot(apex[0],apex[1],'or')
#                                plt.axis('auto')
#                                plt.show()
                    #else:
                    #    th = self.L.Gs.pos[nstr]
                    #    th = np.array([th,th])
                    #    ratio = 1
                        #print self.cpt,ratio,ratio2
                        #if (ratio>0.1) and (ratio2==0):
                        #     pdb.set_trace()
                        #print d_excess,dist_excess_max
                        #if (ratio2 > self.threshold) and (d_excess<dist_excess_max):
                        if (ratio2 > self.threshold) and (d_excess < dist_excess_max):
                        #if (ratio > self.threshold):
                            #
                            # Update sequence of mirrored points
                            #
                            if nstr<0:
                                tahe.append(th)
                            else:
                                tahe.append(th_mirror)
                            #if (tahe[-1][0]==tahe[-1][1]).all():
                            #    pdb.set_trace()
                            # 
                            # Check if the target has been reached
                            # sequence is valid and last interaction is in the list of targets   
                            #if (interaction in lit) or (interaction[-1]==self.target):
                            if (interaction in lit):
                                # idea here is to produce signature without any airwalls
                                # lawp_tmp is a mask where 0 mean no air wall and 1 = airwall
                                # anstr does not contains airwalls
                                # lawp_tmp = [0]+lawp
                                # lll = [x[0] for ix,x in enumerate(visited) if lawp_tmp[ix]==1]
                                # print([self.L.Gs.node[x]['name'] for x in lll])

                                #anstr = np.array([x[0] for ix,x in enumerate(visited) 
                                #                                  if ((lawp[ix]!=1) or (x[0] in self.L.name['AIR']) or (x in (lit+lis)))] )
                                #typ  = np.array([len(x) for ix,x in enumerate(visited) 
                                #                                  if ((lawp[ix]!=1) or (x[0] in self.L.name['AIR']) or (x in (lit+lis)))] )
                                #sig = np.array([anstr,typ])
                                #sighash = hash(str(sig))


                                # if len(anstr) == 2:
                                #     if (anstr == np.array([323,351])).all():
                                #         import ipdb
                                #         ipdb.set_trace()
                                anstr = np.array([x[0] for x in visited ])
                                typ  = np.array([len(x) for x in visited])
                                sig = np.array([anstr,typ])
                                sighash = hash(str(sig))
                                if sighash not in lhash:
                                    lhash.append(sighash)
                                    try:
                                        self[len(typ)] = np.vstack((self[len(typ)],sig))
                                        self.ratio[len(typ)] = np.append(self.ratio[len(typ)],ratio)
                                    except:
                                        self[len(typ)] = np.vstack((sig))
                                        self.ratio[len(typ)] = np.array([ratio])
                                # print ('added',visited)
                                    cptsig +=1

                                if animation:
                                    Nf = nx.draw_networkx_nodes(Gi,pos=Gi.pos,
                                            nodelist=visited,labels={},
                                            node_color='b',
                                            node_size=40,
                                            ax=ax,fig=fig)
                                    Ef = nx.draw_networkx_edges(Gi,pos=Gi.pos,
                                            edgelist=edge,labels={},
                                            width=0.1,arrows=False,
                                            ax=ax,fig=fig)
                                    cpt=cpt+1
                                    plt.savefig('./figure/' +str(us) +'_' + str(cpt) +'.png')
                                    try:
                                        ax.collections.remove(Nf)
                                    except:
                                        pass
                                    try:
                                        ax.collections.remove(Ef)
                                    except:
                                        pass
                            #
                            # output is no longer a dict but a list
                            #
                            outint = Gi[visited[-2]][interaction]['output']
                            #outint = Gi[visited[-2]][interaction]['output'].keys()
                            #
                            # proint not used 
                            #
                            # proint = Gi[visited[-2]][interaction]['output'].values()
                            nexti  = [it for it in outint ]
                            stack.append(iter(nexti))
                        # 1590 ratio <= threshold
                        else:
                            if len(visited)>1:
                                if ((len(visited[-2])==2) or len(visited[-2])==1):
                                    R.pop()
                            last = visited.pop()
                            lawp.pop()
                    # 1389 condR and condT and condD
                    else:
                        pass
                # 1388 cond1 and cond2 and cond3
                else:
                    # if at least 2 interactions
                    # and antepenultiem is a reflexion
                    if len(visited)>1:
                        if ((len(visited[-2])==2) or len(visited[-2])==1):
                            R.pop()
                    last = visited.pop()
                    #
                    # Poping
                    #      tahe
                    #      lawp
                    #      stack
                    #if (tahe[-1][0]==tahe[-1][1]).all():
                    #    pdb.set_trace()
                    tahe.pop()
                    try:
                        lawp.pop()
                    except:
                        pass
                    stack.pop()
                    #stack.pop()

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
        """  plot signatures in the simulated environment

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
                    'crx':-1,
                    'aw':True
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


        if kwargs['s'] == -1:
            for i in lgrint:
                lsig = range(int(len(self[i])/2))
                for j in lsig:
                    sig = [ self.L.Gs.pos[x] for x in self[i][2*j] ]
                    siga = np.array(sig)
                    ax.plot(siga[:,0], siga[:,1],
                            alpha = kwargs['alphasig'],
                            color = kwargs['colsig'],
                            linewidth = kwargs['widthsig'])
                    ax.axis('off')
        else:
            lsig = [kwargs['s']]
            for s1 in lsig:
                sig = [ self.L.Gs.pos[x[0]]  for x in s1]
                siga = np.array(sig)
                ax.plot(siga[:,0], siga[:,1],
                        alpha = kwargs['alphasig'],
                        color = kwargs['colsig'],
                        linewidth = kwargs['widthsig'])
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
                print("incorrect number of interactions")
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
                print("signature index out of bounds of signature")

            line = np.vstack((line,prx))
            ax.plot(line[:,0],line[:,1])
            plt.draw()
            print(inter)
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
                print('press n for next signature')


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

        Todo : Find the best memory implementation

        See Also
        --------

        Signature.sig2ray
        Signature.raysv

        """

        if type(ptx) == int:
            ptx = np.array(self.L.Gt.pos[ptx])

        if type(prx) == int:
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

        #
        # Handling LOS ray
        #  Distance between tx and rx should be sufficient
        #
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


    def raysv(self, ptx=0, prx=1):
        """ transform dict of signatures into 2D rays - default vectorized version

        Parameters
        ----------

        ptx : numpy.array or int
            Tx coordinates is the center of gravity of the cycle ptx if
            type(ptx)=int
        prx :  numpy.array or int
            Rx coordinates is the center of gravity of the cycle prx if
            type(prx)=int

        Returns
        -------

        rays : Rays

        Notes
        -----

        This is a vectorized version of Signatures.rays.
        This implementation takes advantage of the np.ndarray
        and calculates images and backtrace for block of signatures.
        A block of signatures gathers all signatures with the same number of interactions.

        For mathematical details see :

        @phdthesis{amiot:tel-00971809,
          TITLE = {{Design of simulation platform joigning site specific radio propagation and human mobility for localization applications}},
          AUTHOR = {Amiot, Nicolas},
          URL = {https://tel.archives-ouvertes.fr/tel-00971809},
          NUMBER = {2013REN1S125},
          SCHOOL = {{Universit{\'e} Rennes 1}},
          YEAR = {2013},
          MONTH = Dec,
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
            ptx= np.r_[ptx, 0.5]

        if len(ptx) == 2:
            prx= np.r_[prx, 0.5]

        rays = Rays(ptx,prx)

        #
        # detect LOS situation
        #
        #
        # cycle on a line between 2 cycles
        # lc  = self.L.cycleinline(self.source,self.target)

        #
        # if source and target are in the same merged cycle
        # and ptx != prx
        #
        los = shg.LineString(((ptx[0], ptx[1]), (prx[0], prx[1])))

        # convex cycle of each point
        cyptx = self.L.pt2cy(ptx)
        cyprx = self.L.pt2cy(prx)


        polyctx = self.L.Gt.node[cyptx]['polyg']
        polycrx = self.L.Gt.node[cyprx]['polyg']

        # The Line of sight situation is detected here
        # dtxtx : square distance between Tx and Rx
        dtxrx = np.sum((ptx-prx)*(ptx-prx))
        if dtxrx>1e-15:
            if polyctx.contains(los):
                rays.los = True
            else:
                rays.los = False

        M = self.image2(ptx)
        R = self.backtrace(ptx,prx,M)

        #
        # Add LOS ray in ray 2D
        #

        if rays.los:
            R[0]= {'sig':np.zeros(shape=(0,0,1)),'pt': np.zeros(shape=(2,1,0))}

        rays.update(R)
        rays.nb_origin_sig = len(self.keys())
        rays.origin_sig_name = self.filename

        return rays

    def backtrace(self, tx, rx, M):
        ''' backtracing betwen tx and rx

        Parameters
        ----------

            tx : ndarray
                position of tx (2,)
            rx : ndarray
                position of tx (2,)
            M : dict
                position of intermediate points obtained from self.image()

        Returns
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

        See Also
        --------

        pylayers.antprop.signature.image

        '''


        if len(tx) > 2:
            tx = tx[:2]
        if len(rx) > 2:
            rx = rx[:2]

        rayp={}

        # loop on number of interactions
        for ninter in self.keys():
            signatures = copy.deepcopy(self[ninter])
            #get segment ids of signature with ninter interactions
            # seg = self[ninter][::2]
            # unegseg=np.where(seg<0)
            # uninegseg,idx = np.unique(seg[unegseg],return_inverse=True)
            # pneg = np.array([self.L.Gs.pos[x] for x in uninegseg])

            # nsig = len(seg)

            # # determine positions of points limiting the semgments
            # #1 get index in L.tahe
            # # 2 get associated position in L.pt


            # utahe = self.L.tahe[:,self.L.tgs[seg]]
            # # pt : (xycoord (2),pt indexes (2),nb_signatures,nb_interactions)
            # pt = self.L.pt[:,utahe]
            

            # ####WARNING BIG TRICK HERE :
            # #### pa and pb are not set as the same value 
            # #### to avoid a singular matrixnext.
            # #### set pa =-pb has no incidence but avoid complex and vain code 
            # #### modification for handling diffractions
            # try:
            #     pt[:,0,unegseg[0],unegseg[1]]=pneg[idx].T
            #     pt[:,1,unegseg[0],unegseg[1]]=-pneg[idx].T
            # except:
            #     pass
            # pt shape =
            # 0 : (x,y) coordinates x=0,y=1
            # 1 : 2 points (linking the semgnet) a=0,b=1
            #2 : nb of found signatures/segments
            # 3 : nb interaction



            ################################
            ###############################
            ####### This part between hash has been copy/paste from self.image2
            ###### should be considered to become a function

            #get segment ids of signature with ninter interactions
            # nid = node id
            nid = self[ninter][::2]
            nsig = len(nid)


            # pt shape =
            # 0 : (x,y) coordinates x=0,y=1
            # 1 : 2 points (linking the nidment) a=0,b=1
            # 2 : nb of found signatures/nidments
            # 3 : nb interactions
            pt = np.empty((2,2,nsig,ninter))


            # 1 negative points
            # seek for diffraction 
            # negative index points are diffraction points
            upoint = np.where(nid<0)
            unipoint,idx = np.unique(nid[upoint],return_inverse=True)

            #get their coordinates
            #
            # TO BE FIXED 
            #
            #upointcoord = self.L.iupnt[-unipoint]
            #pointcoord = self.L.pt[:,upointcoord]

            pointcoord = np.array([ (self.L.Gs.pos[x][0],self.L.Gs.pos[x][1])  for x in unipoint ]).T
            # #### WARNING BIG TRICK HERE :
            # #### pa and pb are not set as the same value 
            # #### to avoid a singular matrixnext.
            # #### set pa =-pb has no incidence but avoid complex and vain code 
            # #### modification for handling diffractions
            try:
                pt[:,0,upoint[0],upoint[1]] = pointcoord[:,idx]
                pt[:,1,upoint[0],upoint[1]] = -pointcoord[:,idx]
            except:
                pass


            # 2 positive points
            # seek for segments
            useg = np.where(nid>0)
            # removing duplicates ( for increasing speed)
            uniseg,idxp = np.unique(nid[useg],return_inverse=True)

            # determine positions of points limiting the nidments
            #1 get index in L.tahe
            utahe = self.L.tahe[:,self.L.tgs[uniseg]]
            segcoord = self.L.pt[:,utahe]

            pt[:,:,useg[0],useg[1]]=segcoord[:,:,idxp]

            ###################################
            ########################################



            # how to do this into a while loop ?
            p=rx

            # creating W matrix required in eq (2.70) thesis Nicolas AMIOT
            # Warning W is rolled after and becomes (nsig,4,4)
            W = np.zeros((4,4,nsig))
            I = np.eye(2)[:,:,np.newaxis]*np.ones((nsig))

            W[:2,:2,...] = I
            W[2:4,:2,...] = I

            # once rolled :
            # W (nsig,4,4)
            W = np.rollaxis(W,-1)


            kinter=ninter-1

            ptr = pt
            Mr = copy.deepcopy(M)

            epsilon = 1e-12
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
                uvalidA = psolved[:,2]>0.
                uvalidB = psolved[:,2]<1.
                #beta
                uvalidC = psolved[:,3] >= epsilon
                uvalidD = psolved[:,3] <=1.-epsilon
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
                N = int(len(usigv)/2)
                sir1=signatures[::2].T.reshape(ninter,N)
                sir2=signatures[1::2].T.reshape(ninter,N)
                sig = np.empty((2,ninter,N))
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
            # nid = node id
            nid = self[ninter][::2]
            nsig = len(nid)
            M = np.empty((2,nsig,ninter))


            # pt shape =
            # 0 : (x,y) coordinates x=0,y=1
            # 1 : 2 points (linking the nidment) a=0,b=1
            # 2 : nb of found signatures/nidments
            # 3 : nb interactions

            try:
                pt = np.nan*np.zeros((2,2,nsig,ninter))
            except:
                pdb.set_trace()

            #1 negative points
            # seek for diffraction 
            # negative index points are diffraction points
            upoint = np.where(nid<0)
            unipoint,idxpt = np.unique(nid[upoint],return_inverse=True)

            #get their coordinates
            #
            # To be FIXED
            #
            #upointcoord = self.L.iupnt[-unipoint]
            #pointcoord = self.L.pt[:,upointcoord]

            pointcoord = np.array([ (self.L.Gs.pos[x][0],self.L.Gs.pos[x][1])  for x in unipoint ]).T
        
            # try except to handle the case where there is no diffraction point
            try:
                pt[:,0,upoint[0],upoint[1]] = pointcoord[:,idxpt]
                pt[:,1,upoint[0],upoint[1]] = pointcoord[:,idxpt]
            except:
                pass


            #2 positive points
            #seek for segments
            useg = np.where(nid>0)
            # removing duplicates ( for increasing speed)
            uniseg,idxseg = np.unique(nid[useg],return_inverse=True)

            # determine positions of points limiting the nidments
            #1 get index in L.tahe
            utahe = self.L.tahe[:,self.L.tgs[uniseg]]
            segcoord = self.L.pt[:,utahe]

            pt[:,:,useg[0],useg[1]]=segcoord[:,:,idxseg]

            # check every element of pt is filled
            assert not np.isnan(pt).any()
            #
            # TODO Upgrading layout for handling slab offsets 
            #
            # uncomment those two lines when the numpy array L.norm and
            # L.offset exist
            #norm    = self.L.normal[:,utahe]
            #offset  = self.L.offset[:,utahe]
            # pt = pt + offset*norm


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

            for n in np.arange(ninter):
                #get segment ids of signature with ninter interactions
                uT = np.where(ityp[:,n]==3)[0]
                uR = np.where(ityp[:,n]==2)[0]
                uD = np.where(ityp[:,n]==1)[0]
                if n ==0:
                    p = tx[:,None]*np.ones((nsig))
                else :
                    p = M[:,:,n-1]
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

        Parameters
        ----------

        tx : ndarray
            position of tx (2,)

        Returns
        -------

        M : dictionnary

        dictionnary of intermediate coordinates
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
                print("signatures.image diffraction line 3672 Not yet tested !")

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
                    print(nsp)
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


class Signature(PyLayers,object):
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

        if type(sig) == np.ndarray:
            self.seq = sig[0, :]
            self.typ = sig[1, :]

        if type(sig) == list:
            self.seq = map(seginter,sig)
            self.typ = map(typinter,sig)

    def __repr__(self):
        s = ''
        s = s + str(self.seq) + '\n'
        s = s + str(self.typ) + '\n'
        if self.evaluated:
            s = s + str(self.pa)+'\n'
            s = s + str(self.pb)+'\n'
        return s

    def info(self):
        for k in self.__dict__.keys():
            print(k, ':', self.__dict__[k])

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
        self.evaluated = True

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
        self.evaluated = True

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

    def show(self,L,tx,rx,**kwargs):
        """

        Parameters
        ----------
        L : Layout 
        tx : 
        rx : 
        aw

        """
        defaults  = {'aw':True,
                     'axes':True,
                     'labels':False,
                     'fig':[],
                     'ax':[]
                     }
        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        if kwargs['fig']==[]:
            fig = plt.gcf()
        else:
            fig = kwargs['fig']
        if kwargs['ax']==[]:
            fig = plt.gcf()  
        else:
            ax = fig.gca()   

        self.ev(L)
        fig,ax = L.showG('s',labels=kwargs['labels'],
                             aw=kwargs['aw'],
                             axes=kwargs['axes']
                             ,fig=fig,ax=ax)
        M = self.image(tx)
        isvalid,Y,tup = self.backtrace(tx,rx,M)
        l1 = ax.plot(tx[0],tx[1],'or')
        l2 = ax.plot(rx[0],rx[1],'og')
        l3 = ax.plot(M[0,:],M[1,:],'ob')
        l4 = ax.plot(Y[0,:],Y[1,:],'ok')
        ray = np.hstack((np.hstack((rx.reshape(2,1),Y)),tx.reshape(2,1)))
        for k in self.seq:
            ax.annotate(str(k),xy=(L.Gs.pos[k]),xytext=(L.Gs.pos[k]))
        if isvalid:
            l5 = ax.plot(ray[0,:],ray[1,:],color='green',alpha=0.6,linewidth=0.6)
        else:
            l5 = ax.plot(ray[0,:],ray[1,:],color='red',alpha=0.6,linewidth=0.6)

        return fig,ax

    def backtrace(self, tx, rx, M):
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
            >>> L = Layout('defstr.ini')
            >>> s = Signature(seq)
            >>> tx = np.array([760,1113])
            >>> rx = np.array([762,1114])
            >>> s.ev(L)
            >>> M = s.image(tx)
            >>> isvalid,Y = s.backtrace(tx,rx,M)

            >>> fig,ax = L.showG('s',labels=1,aw=1,axes=1)
            >>> l1 = ax.plot(tx[0],tx[1],'or')
            >>> l2 = ax.plot(rx[0],rx[1],'og')
            >>> l3 = ax.plot(M[0,:],M[1,:],'ob')
            >>> l4 = ax.plot(Y[0,:],Y[1,:],'xk')
            >>> ray = np.hstack((np.hstack((tx.reshape(2,1),Y)),rx.reshape(2,1)))
            >>> l5 = ax.plot(ray[0,:],ray[1,:],color='#999999',alpha=0.6,linewidth=0.6)
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
        epsilon = 1e-12
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
                Y = np.hstack((Y, pa[:, N-(k+1)].reshape((2, 1))))
                pkm1 = pa[:, N-(k+1)].reshape((2, 1))
            k = k + 1
        if ((k == N) & ((beta > 0) & (beta < 1)) & ((alpha > 0) & (alpha < 1))):
            Y = np.hstack((Y, tx.reshape(2, 1)))
            return isvalid,Y,(k,alpha,beta)
        else:
            isvalid = False
            return isvalid,Y,(k,alpha,beta)


    def sig2ray(self, L, pTx, pRx):
        """ convert a signature to a 2D ray

        Parameters
        ----------

        L : Layout
        pTx : ndarray
            2D transmitter position
        pRx : ndarray
            2D receiver position

        Returns
        -------

        Y : ndarray (2x(N+2))

        See Also
        --------

        Signature.image
        Signature.backtrace

        """

        # ev transforms a sequence of segment into numpy arrays (points)
        # necessary for image calculation
        self.ev(L)
        # calculates images from pTx
        M = self.image(pTx)

        #print self
        #if np.array_equal(self.seq,np.array([5,7,4])):
        #    pdb.set_trace()
        isvalid,Y,u = self.backtrace(pTx, pRx, M)
        #print isvalid,Y
        #
        # If incremental mode this function returns an alternative signature
        # in case the signature do not yield a valid ray.
        #
        return isvalid,Y,u




if __name__ == "__main__":
    plt.ion()
    print("testing pylayers/antprop/signature.py")
    doctest.testmod()
    print("-------------------------------------")
