#!/usr/bin/python
# -*- coding: latin1 -*-
"""

Class Rays
==========

This modules contains Rays class

.. autosummary::
    :toctree: generated

    Rays.__init__
    Rays.__len__
    Rays.__repr__
    Rays.sort
    Rays.extract
    Rays.mirror
    Rays.to3D
    Rays.locbas
    Rays.fillinter
    Rays.length
    Rays.eval
    Rays.ray
    Rays.typ
    Rays.info
    Rays.signature
    Rays.show
    Rays.show3d
    Rays._show3
    Rays.reciprocal
    Rays.check_reciprocity

"""
import pdb
import os
import copy
import ConfigParser
import glob
import doctest
import networkx as nx
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import struct as stru
import pylayers.util.geomutil as geu
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.antprop.interactions import *
from pylayers.antprop.slab import *
from pylayers.antprop.channel import Ctilde
from pylayers.gis.layout import Layout
import pylayers.signal.bsignal as bs
import h5py
try:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi import mlab
except:
    print 'Layout:Mayavi is not installed'

class Rays(PyLayers,dict):
    """ A set af rays

    Attributes
    ----------

    rays   :
    nbrays :
    rayidx :
    sig    :
    pt     :
    alpha  :


    Methods
    -------

    to3D(H=3,N=1)
        for k in self:   # for all interaction group k
        for k in self:   # for all interaction group k
    locbas(L)
    fillinter(L)
    eval
    show(L)
    mirror(H=3,N=1)
    ray
    typ
    info
    to3D
    signature(L)
    show3d(ray,bdis,bbas,bstruc,col,id,linewidth)
    show3()

    Notes
    -----

    The Rays object is obtained from a signature.
    It is a container for a set of rays between a source
    and a target point defining a radio link.

    Once a Rays object has been obtained in 2D, it is transformed
    in 3D via the **to3D** method. This method has two parameters :
    the height from floor to ceil, and the number N of
    multiple reflections to take into account.

    Once the 3d rays have been calculated, 
    the local basis are evaluated along those rays. This is
    done through the **locbas** method

    Once the local basis have been calculated the different
    interactions along rays can be informed via the **fillinter**
    method.

    Once the interaction are informed the field along rays can
    be evaluated via the **eval** method
    """
    def __init__(self, pTx, pRx):
        """ object constructor

        Parameters
        ----------

        pTx : np.array
        pRx : np.array

        """

        self.pTx = pTx
        self.pRx = pRx
        self.nray = 0
        self.raypt = 0
        self.los = False
        self.is3D = False
        self.isbased = False
        self.filled = False
        self.evaluated = False

    def __len__(self):
        Nray = 0
        for k in self.keys():
            sh = np.shape(self[k]['sig'])
            Nray = Nray+sh[2]
        return Nray


    def __repr__(self):
        s = ''
        ni = 0
        nl = 0

        try:
            if self.is3D:
                s = self.__class__.__name__ + '3D\n' + '----------'+'\n'

                for k in self:
                    r = self[k]['rayidx']
                    nr = len(r)
                    s = s + str(k)+' / '+str(nr)+ ' : '+str(r)+'\n'
                    ni = ni + nr*k
                    nl = nl + nr*(2*k+1)
                s = s + '-----'+'\n'
                s = s+'ni : '+str(ni)+'\n'
                s = s+'nl : '+str(nl)+'\n'
            else:
                s = self.__class__.__name__ + '2D\n' + '----------'+'\n'

                nray = np.sum([np.shape(self[i]['sig'])[2] for i in self.keys()])
                s = 'N2Drays : '+ str(nray) + '\n'
                s = s + 'from '+ str(self.nb_origin_sig) + ' signatures\n'
                s = s + '#Rays/#Sig: '+ str( len(self)/(1.*self.nb_origin_sig) )

                s = s + '\npTx : '+ str(self.pTx) + '\npRx : ' + str(self.pRx)+'\n'

                for k in self:
                    #sk = np.shape(self[k]['sig'])[2]
                    s = s + str(k) + ': '+ str(self[k]['sig'][0,:])+'\n'
                    #s = s + str(sk) + 'rays with' + str(k) + ' interactions'
        except:
            print "problem"
            return(s)

        return(s)


    def saveh5(self,idx=0):
        """ save rays in hdf5 format

        Parameters
        ----------

        idx : int

        See Also
        --------

        loadh5

        """

        filename = self.filename+'_'+str(idx)
        filenameh5=pyu.getlong(filename+'.h5',pstruc['DIRR3D'])



        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f=h5py.File(filenameh5,'w')
            # keys not saved as attribute of h5py file
            notattr = ['I','B','B0','delays','dis']
            for a in self.__dict__.keys():
                if a not in notattr:
                    f.attrs[a]=getattr(self,a)

            for k in self.keys():
                f.create_group(str(k))
                for kk in self[k].keys():
                    if kk == 'sig2d':
                        # Need to find an efficient way to save the signatures
                        # 2d which have created the rays
                        pass
                    elif kk == 'nbrays':
                        f[str(k)].create_dataset(kk,shape=(1,),data=np.array([self[k][kk]]))
                    else:
                        f[str(k)].create_dataset(kk,shape=np.shape(self[k][kk]),data=self[k][kk])
            f.close()
        except:
            f.close()
            raise NameError('Rays: issue when writting h5py file')



    def loadh5(self,filename=[],idx=0):
        """ load rays hdf5 format

        Parameters
        ----------

        """
        if filename == []:
            filenameh5 = self.filename+'_'+str(idx)+'.h5'
        else :
            filenameh5 = filename

        filename=pyu.getlong(filenameh5,pstruc['DIRR3D'])


        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f=h5py.File(filename,'r')
            for k in f.keys():
                self.update({eval(k):{}})
                for kk in f[k].keys():
                    self[eval(k)].update({kk:f[k][str(kk)][:]})

            for a,va in f.attrs.items():
                setattr(self,a,va)
            f.close()

        except:

            f.close()
            raise NameError('Rays: issue when reading h5py file')

        # fill if save was filled

        # temporary solution in order to avoir
        # creating save for Interactions classes
        if self.filled:
            Lname = self.filename.split('_')[0] + '.ini'
            L=Layout(Lname)
            self.fillinter(L)

        if self.evaluated:
            return self.eval(self.fGHz)

    def _saveh5(self,filenameh5,grpname):
        """ Save rays h5py format compliant with Links Class

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

        filenameh5=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:

            fh5=h5py.File(filenameh5,'a')
            if not grpname in fh5['ray'].keys():
                fh5['ray'].create_group(grpname)
            else :
                print 'ray/'+grpname +'already exists in '+filenameh5
            f = fh5['ray/'+grpname]
            # keys not saved as attribute of h5py file
            notattr = ['I','B','B0','delays','dis']
            for a in self.__dict__.keys():
                if a not in notattr:
                    f.attrs[a]=getattr(self,a)

            for k in self.keys():
                f.create_group(str(k))
                for kk in self[k].keys():
                    if kk == 'sig2d':
                        # Need to find an efficient way to save the signatures
                        # 2d which have created the rays
                        pass
                    elif kk == 'nbrays':
                        f[str(k)].create_dataset(kk,shape=(1,),data=np.array([self[k][kk]]))
                    else:
                        f[str(k)].create_dataset(kk,shape=np.shape(self[k][kk]),data=self[k][kk])
            fh5.close()
        except:
            fh5.close()
            raise NameError('Rays: issue when writting h5py file')

    def _loadh5(self,filenameh5,grpname):
        """ load rays  h5py format compliant with Links Class

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
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            fh5=h5py.File(filename,'r')
            f = fh5['ray/'+grpname]
            for k in f.keys():
                self.update({eval(k):{}})
                for kk in f[k].keys():
                    self[eval(k)].update({kk:f[k][str(kk)][:]})

            for a,va in f.attrs.items():
                setattr(self,a,va)
            fh5.close()

        except:

            fh5.close()
            raise NameError('Rays: issue when reading h5py file')

        # fill if save was filled

        # temporary solution in order to avoid
        # creating save for Interactions classes

        if self.filled:
            L=Layout(Lfilename)
            self.fillinter(L)

        if self.evaluated:
            return self.eval(self.fGHz)


    def reciprocal(self):
        """ switch tx and rx

        """


        r = Rays(self.pRx,self.pTx)
        r.is3D = self.is3D
        r.nray = self.nray
        r.origin_sig_name = self.origin_sig_name
        r.nb_origin_sig = self.nb_origin_sig

        for k in self:
            r[k]={}
            r[k]['pt']=self[k]['pt'][:,::-1,:]
            r[k]['sig']=self[k]['sig'][:,::-1,:]
        return(r)


    def check_reciprocity(self,r):
        """ check ray reciprocity

        Parameters
        ----------

        r : Rays
        """
        assert (self.pTx==r.pRx).all()
        assert (self.pRx==r.pTx).all()
        for k in self:
            assert (np.allclose(self[k]['dis'],r[k]['dis']))
            assert (np.allclose(self[k]['pt'],r[k]['pt'][:,::-1,:]))
            assert (np.allclose(self[k]['sig'],r[k]['sig'][:,::-1,:]))
            if (self.isbased) & (r.isbased):
                #assert (np.allclose(self[k]['nstrwall'],r[k]['nstrwall'][:,::-1,:]))
                assert (np.allclose(self[k]['norm'],r[k]['norm'][:,::-1,:]))
                #assert ((np.mod(self[k]['aoa']-r[k]['aod'],2*np.pi)==0).all())
                #assert ((np.mod(self[k]['aod']-r[k]['aoa'],2*np.pi)==0).all())
                assert (np.allclose(self[k]['Bo0'],r[k]['BiN']))
                assert (np.allclose(self[k]['BiN'],r[k]['Bo0']))
                assert (np.allclose(self[k]['vsi'],-r[k]['vsi'][:,::-1,:]))
                assert (np.allclose(abs(self[k]['scpr']),abs(r[k]['scpr'][::-1,:])))
                assert (np.allclose(self[k]['theta'],r[k]['theta'][::-1,:]))
                assert (np.allclose(self[k]['Bi'],r[k]['Bo'][:,:,::-1,:]))
                assert (np.allclose(self[k]['Bo'],r[k]['Bi'][:,:,::-1,:]))
                assert (np.allclose(self[k]['B'],r[k]['B'][:,:,::-1,:].swapaxes(0,1)))

        if self.evaluated :

            for ir in range(self.nray):

                iint1 = self.ray(ir)
                iint2 = r.ray(ir)

                # check Interactions
                A1 = self.I.I[:, iint1, :, :]
                A2 = r.I.I[:, iint2, :, :][:,::-1,:,:]
                assert np.allclose(A1,A2),pdb.set_trace()

                # check bases
                #  ray 1 : B0   | B[0]   | B[1] | B[2] | B[3] | B[4]
                #  ray 2 : B[4] | B[3]  | B[2]  | B[1] | B[0] | B0
                assert np.allclose(self.B0.data[ir,:,:],r.B.data[iint2,:,:][-1,:,:].swapaxes(1,0))
                assert np.allclose(r.B0.data[ir,:,:],self.B.data[iint1,:,:][-1,:,:].swapaxes(1,0))
                assert np.allclose(self.B.data[iint1,:,:][:-1],r.B.data[iint2,:,:][:-1][::-1,:,:].swapaxes(2,1))



    def sort(self):
        """ sort rays
        """
        u = np.argsort(self.dis)

    def extract(self,ni,nr):
        """ Extract a single ray

        Parameters
        ----------

        ni : group of interactions
        nr : ray index in group of interactions

        """


        r = Rays(self.pTx,self.pRx)
        r.is3d = self.is3D
        r[ni] = {}

        for k in self[ni].keys():
            tab  = self[ni][k]
            if type(tab)==np.ndarray:
                r[ni][k] = tab[...,nr][...,np.newaxis]
        r[ni]['nrays']=1
        return(r)

    def show(self,**kwargs):
        """  plot 2D rays within the simulated environment

        Parameters
        ----------

        i : list or -1 (default = all groups)
            list of interaction group numbers
        r : list or -1 (default = all rays)
            list of indices of ray in interaction group
        graph : string t
            type of graph to be displayed
            's','r','t',..
        fig : figure
        ax  : axis
        L   : Layout
        alpharay : float
            1
        widthray : float
            0.1
        colray : string
            'black'
        ms : int
            marker size :  5
        layout : boolean
            True
        points : boolean
            True

        """
        defaults = {'i':-1,
                   'r':-1,
                   'fig':[],
                   'ax':[],
                    'L':[],
                   'graph':'s',
                    'color':'black',
                    'alpharay':1,
                    'widthray':0.1,
                    'colray':'black',
                    'ms':5,
                    'layout':True,
                    'points':True
                   }
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        #if kwargs['fig'] ==[]:
        #    fig = plt.figure()
        #if kwargs['ax'] ==[]:
        #    ax = fig.add_subplot(111)
        if kwargs['layout'] ==True:
            fig,ax = kwargs['L'].showG(**kwargs)
        else:
            fig = kwargs['fig']
            ax = kwargs['ax']
        #
        # display Tx and Rx
        #
        if kwargs['points'] ==True:
            ax.plot(self.pTx[0], self.pTx[1], 'or',ms=kwargs['ms'])
            ax.plot(self.pRx[0], self.pRx[1], 'og',ms=kwargs['ms'])
        # i=-1 all rays
        # else block of interactions i

        if kwargs['i']==-1:
            lgrint = self.keys()
        else:
            lgrint = [kwargs['i']]

        for i in lgrint:
            if kwargs['r']==-1:
                lray = range(len(self[i]['pt'][0, 0, :]))
                if self.filled :
                    ax.set_title('rays index :'+ str(self[i]['rayidx']))
            else:
                lray = [kwargs['r']]
            for j in lray:
                ray = np.hstack((self.pTx[0:2].reshape((2, 1)),
                                 np.hstack((self[i]['pt'][0:2, :, j],
                                 self.pRx[0:2].reshape((2, 1))))
                                 ))
                ax.plot(ray[0, :], ray[1, :],
                        alpha=kwargs['alpharay'],color=kwargs['colray'],linewidth=kwargs['widthray'])
                ax.axis('off')
                if self.filled :
                    ax.set_title('rays index :'+ str(self[i]['rayidx'][lray]))
        return(fig,ax)

    def mirror(self, H=3, N=1):
        """ mirror a ray termination

        Parameters
        ----------

        H : float
            ceil height (default 3m)

        N : int
            handle the number of mirror reflexions

        Returns
        -------

        d : dict
            k : zm  v: alpham
            k : zp  v: alphap

        Examples
        --------

        >>> ptx = np.array([1,1,1.5])
        >>> prx = np.array([2,2,1.2])
        >>> r = Rays(ptx,prx)
        >>> d = r.mirror()
        >>> d[-1.5]
        array([ 0.55555556])

        """
        km = np.arange(-N+1, N+1, 1)
        kp = np.arange(-N, N+1, 1)
        #
        # heights of transmitter and receiver
        #
        ht = self.pTx[2]
        hr = self.pRx[2]

        zkp = 2*kp*H + ht
        zkm = 2*km*H - ht

        d = {}

        for zm in zkm:
            if zm < 0:
                bup = H
                pas = H
                km = int(np.ceil(zm/H))
            else:
                bup = 0
                pas = -H
                km = int(np.floor(zm/H))
            thrm = np.arange(km*H, bup, pas)
            d[zm] = abs(thrm-zm)/abs(hr-zm)

        for zp in zkp:
            if zp < 0:
                bup = H
                pas = H
                kp = int(np.ceil(zp/H))
            else:
                bup = 0
                pas = -H
                kp = int(np.floor(zp/H))
            thrp = np.arange(kp*H, bup, pas)
            d[zp] = abs(thrp-zp)/abs(hr-zp)

            # print "zp",zp
            # print "kp",kp
            # print "thrp",thrp
            # print "alphap",d[zp]

        return(d)

    def to3D(self,L,H=3, N=1):
        """ transform 2D ray to 3D ray

        Parameters
        ----------

        L : Layout object

        H : float
            ceil height (default 3m)

        N : int
            number of mirror reflexions

        returns
        -------

        r3d : Rays

        Notes
        -----

        1.


        """

        tx = self.pTx
        rx = self.pRx

        #
        # Phase 1 : calculate Tx images height and parameterization in the
        # vertical plane
        #

        d = self.mirror(H=H, N=N)

        #
        # Phase 2 : calculate 2D parameterization in the horizontal plane
        #

        # for all group of interactions
        for i in self:

            pts = self[i]['pt'][0:2, :, :]
            sig = self[i]['sig']
            # broadcasting of t and r
            t = self.pTx[0:2].reshape((2, 1, 1)) * \
                np.ones((1, 1, len(pts[0, 0, :])))
            r = self.pRx[0:2].reshape((2, 1, 1)) * \
                np.ones((1, 1, len(pts[0, 0, :])))
            # append t and r to interaction points in 2D
            pts1 = np.hstack((t, np.hstack((pts, r))))
            si1 = pts1[:, 1:, :] - pts1[:, :-1, :]
            # array of all ray segments distances
            si = np.sqrt(np.sum(si1 * si1, axis=0))
            # array of cumulative distance of 2D ray
            al1 = np.cumsum(si, axis=0)

            # initialize parameterization parameter alpha
            self[i]['alpha'] = np.zeros(np.shape(si[:-1, :]))

            for j in range(len(self[i]['alpha'][:, 0])):
                # get alpha
                self[i]['alpha'][j, :] = np.sum(si[0:j+1, :], axis=0) \
                        /np.sum(si, axis=0)
                # get z coordinate
                self[i]['pt'][2, j, :] = tx[2] + self[i]['alpha'][j, :] \
                    * (rx[2] - tx[2])

        #
        #  Phase 3 : Initialize 3D rays dictionnary
        #
        r3d = Rays(tx, rx)
        r3d.los = self.los
        r3d.is3D = True

        #
        # Phase 4 : Fill 3D rays information
        #
        # Two nested loops
        #
        #      for all interaction group
        #          for all type of 3D rays
        #             1) extension
        #             2) sort
        #             3) coordinates as a function of parameter
        #
        for k in self:   # for all interaction group k
            # k = int(k)
            # Number of rays in interaction group k
            Nrayk = np.shape(self[k]['alpha'])[1]

            # get  2D horizontal parameterization
            a1 = self[k]['alpha']

            #if (k==1):
            #    pdb.set_trace()
            # get  2D signature
            sig = self[k]['sig']
            #print "signatures 2D ",sig
            #print "----"
            sigsave = copy.copy(sig)
            # add parameterization of tx and rx (0,1)
            a1 = np.concatenate((np.zeros((1, Nrayk)), a1, np.ones((1, Nrayk))))

            # reshape signature in adding tx and rx
            sig = np.hstack((np.zeros((2, 1, Nrayk), dtype=int),
                             sig,
                             np.zeros((2, 1, Nrayk), dtype=int)))  # add signature of Tx and Rx (0,0))
            # broadcast tx and rx
            Tx = tx.reshape(3, 1, 1)*np.ones((1, 1, Nrayk))
            Rx = rx.reshape(3, 1, 1)*np.ones((1, 1, Nrayk))

            # pte is the sequence of point in 3D ndim =3   ( ndim x k x Nrayk)
            pte = self[k]['pt']

            # ndim x k+2 x Nrayk
            pte = np.hstack((Tx, pte, Rx))

            for l in d:                     # for each vertical pattern (C,F,CF,FC,....)
                #print k,l,d[l]
                Nint = len(d[l])            # number of additional interaction
                #if ((k==1) & (l==5.0)):
                #    pdb.set_trace()
                if Nint > 0:                # if new interaction ==> need extension
                    # a1e : extended horizontal+vertical parameterization
                    a1e = np.concatenate((a1, d[l].reshape(len(d[l]), 1)*
                                          np.ones((1, Nrayk))))
                    # get sorted indices
                    ks = np.argsort(a1e, axis=0)
                    # a1es : extended sorted horizontal + vertical parameterization
                    a1es = np.sort(a1e, axis=0)

                    # #### Check if it exist same parameter value  in horizontal plane
                    # #### and vertical plane. Move parameter is so.
                    da1es = np.diff(a1es,axis=0)
                    pda1es = np.where(da1es<1e-10)
                    a1es[pda1es]=a1es[pda1es]-1e-3


                    # prepare an extended sequence of points ( ndim x  (Nint+k+2) x Nrayk )
                    ptee = np.hstack((pte, np.zeros((3, Nint, Nrayk))))

                    #
                    # Boolean ceil/floor detector
                    #
                    # u is 4 (floor interaction )
                    #      5 (ceil interaction )
                    #  depending on the vertical pattern l.
                    #
                    #  l <0 corresponds to last reflexion on floor
                    #  l >0 corresponds to last reflexion on ceil
                    #
                    # u =0 (floor) or 1 (ceil)
                    if l < 0:
                        u = np.mod(range(Nint), 2)
                    else:
                        u = 1 - np.mod(range(Nint), 2)
                    #
                    u = u + 4
                    #
                    # At that point we introduce the signature of the new
                    # introced points on the ceil and/or floor.
                    #
                    # A signature is compose of two lines
                    # esigs sup line : interaction number
                    # esigi inf line : interaction type
                    #
                    esigs = np.zeros((1, Nint, Nrayk), dtype=int)
                    esigi = u.reshape(1, Nint, 1)* np.ones((1, 1, Nrayk), dtype=int)
                    # esig : extension of the signature
                    esig = np.vstack((esigs, esigi))
                    # sige : signature extended  ( 2 x (Nint+k+2) x Nrayk )
                    sige = np.hstack((sig, esig))

                    #
                    # 2 x (Nint+k+2) x Nrayk
                    #
                    # sort extended sequence of points
                    # and extended sequence of signatures with the sorting
                    # index ks obtained from argsort of merge parametization
                    #
                    # sequence of extended sorted points
                    #
                    ptees = ptee[:, ks, range(Nrayk)]
                    siges = sige[:, ks, range(Nrayk)]
                    # extended and sorted signature
                    iint_f, iray_f = np.where(siges[ 1, :] == 4)  # floor interaction
                    iint_c, iray_c = np.where(siges[ 1, :] == 5)  # ceil interaction
                    #print siges
                    #
                    # find the list of the previous and next point around the
                    # new ceil or floor point. The case of successive ceil or
                    # floor reflexion make
                    #
                    # Tous les points précédents qui ne sont pas des Ceils ou
                    # des floors et tous les points suivants qui ne sont pas
                    # des points de réflexion ceil ou floor
                    #
                    # Afin de tenir compte du rayon et du groupe d'interaction
                    # concerné, il faut passer un tuple qui concatène la valeur
                    # de l'indice d'interaction floor ou ceil et l'indice de
                    # rayons du groupe associé (d'ou le zip)
                    #
                    # Cette séquence d'instruction fixe le bug #133
                    #
                    # Antérieurement il y avait une hypothèse de succession
                    # immediate d'un point 2D renseigné.
                    #
                    try:
                        iintm_f = map(lambda x : np.where( (siges[1,0:x[0],x[1]]<>4) & (siges[1,0:x[0],x[1]]<>5))[0][-1], zip(iint_f,iray_f))
                        iintp_f = map(lambda x : np.where( (siges[1,x[0]:,x[1]]<>4) & (siges[1,x[0]:,x[1]]<>5))[0][0]+x[0], zip(iint_f,iray_f))
                        iintm_c = map(lambda x : np.where( (siges[1,0:x[0],x[1]]<>4) & (siges[1,0:x[0],x[1]]<>5))[0][-1], zip(iint_c,iray_c))
                        iintp_c = map(lambda x : np.where( (siges[1,x[0]:,x[1]]<>4) & (siges[1,x[0]:,x[1]]<>5))[0][0]+x[0], zip(iint_c,iray_c))
                    except:
                        pdb.set_trace()

                    # Update coordinate in the horizontal plane
                    #
                    #
                    # The new interaction ceil or floor has no coordinates in
                    # the horizontal plane.
                    # Those coordinates are evaluated first by finding a sub
                    # parameterization of the point with respect to the two
                    # known adjascent interaction point j-1 and j+1 (Thales)
                    #

                    #iintm_f = iint_f - 1
                    #iintp_f = iint_f + 1

                    #iintm_c = iint_c - 1
                    #iintp_c = iint_c + 1


                    if len(iint_f)>0:
                        a1esm_f = a1es[iintm_f, iray_f]
                        a1esc_f = a1es[iint_f, iray_f]
                        a1esp_f = a1es[iintp_f, iray_f]


                        pteesm_f = ptees[0:2, iintm_f, iray_f]
                        pteesp_f = ptees[0:2, iintp_f, iray_f]

                        coeff_f = (a1esc_f-a1esm_f)/(a1esp_f-a1esm_f)
                        ptees[0:2, iint_f, iray_f] = pteesm_f + coeff_f*(pteesp_f-pteesm_f)

                    if len(iint_c)>0:
                        a1esm_c = a1es[iintm_c, iray_c]
                        a1esc_c = a1es[iint_c, iray_c]
                        a1esp_c = a1es[iintp_c, iray_c]

                        pteesm_c = ptees[0:2, iintm_c, iray_c]
                        pteesp_c = ptees[0:2, iintp_c, iray_c]

                        coeff_c = (a1esc_c-a1esm_c)/(a1esp_c-a1esm_c)
                        ptees[0:2, iint_c, iray_c] = pteesm_c + coeff_c*(pteesp_c-pteesm_c)


                    #    a1es[iint_f+1, iray_f]-a1es[iint_f-1, iray_f])
                    #coeff_f = (a1es[iint_f, iray_f]-a1es[iint_f-1, iray_f])/(
                    #    a1es[iint_f+1, iray_f]-a1es[iint_f-1, iray_f])

                    #coeff_c = (a1es[iint_c, iray_c]-a1es[iint_c-1, iray_c])/(
                    #    a1es[iint_c+1, iray_c]-a1es[iint_c-1, iray_c])
                    #
                    # Update coordinate in the horizontal plane
                    #
                    #
                    #ptees[0:2, iint_f, iray_f] = ptees[0:2, iint_f-1, iray_f] + coeff_f*(
                    #    ptees[0:2, iint_f+1, iray_f]-ptees[0:2, iint_f-1, iray_f])
                    #ptees[0:2, iint_c, iray_c] = ptees[0:2, iint_c-1, iray_c] + coeff_c*(
                    #    ptees[0:2, iint_c+1, iray_c]-ptees[0:2, iint_c-1, iray_c])


                    # vertical plane
                    # WARNING !!
                    #
                    # ptees[2,iint_f,iray_f]   = 0
                    # ptees[2,iint_c,iray_c]   = H
                    #
                    z  = np.mod(l+a1es*(rx[2]-l), 2*H)
                    pz = np.where(z > H)
                    z[pz] = 2*H-z[pz]
                    ptees[2, :] = z
                # recopy old 2D parameterization (no extension)
                else:
                    a1es = a1
                    ks = np.argsort(a1es, axis=0)
                    ptees = pte
                    # fixing bug
                    siges = copy.copy(sig)
                    #print siges

                #---------------------------------
                # handling subsegments (if any)
                #---------------------------------
                #
                #   ptes (3 x i+2 x r )
                if L.Nss>0:
                    # lsss[k] = n means subsegment k belongs to segment n
                    # a same segment n can have several subsegments
                    # (multi-subsegment case) that is the reason of unique
                    lsss = np.unique(np.array(L.lsss))

                    # index of signature which corresponds to subsegment
                    u   = map(lambda x: list(np.where(siges[0,:,:]==x)),lsss)[0]

                    # dimension extension of index u for :
                    #    z coordinate extraction (append line 2 on dimension 0)
                    #    0 signature extraction  (append line 0 on  dimension 0)

                    # v : index 2 is for getting z coordinate
                    # w : index 0 is for getting segment number (first line of
                    # siges)
                    v   = [2*np.ones(len(u[0]),dtype=int)]+u
                    w   = [0*np.ones(len(u[0]),dtype=int)]+u

                    # zss : height of interactions on subsegments
                    zss = ptees[v]
                    #if k==1:
                    #    print "l",l
                    #    print "ptees : ",ptees
                    #    print "u ",u
                    #    print "v ",v
                    #    print "w ",w
                    # structure index of corresponding subsegments

                    # nstrs = [ 1 , 1 , 1 , 1, 1 , 1]
                    nstrs = siges[w]

                    #if k==1:
                    #    print "nstrs: ",nstrs
                    #    print "zss:",zss
                    #
                    # Determine which subsegment has been intersected
                    # k = 0 : no subsegment intersected
                    zinterval = map(lambda x: L.Gs.node[x]['ss_z'],nstrs)


                    # Example
                    #
                    # zinterval = [[(0.0, 2.4), (2.7, 2.8), (2.8, 3)],
                    #              [(0.0, 2.4), (2.7, 2.8), (2.8, 3)],
                    #              [(0.0, 2.4), (2.7, 2.8), (2.8, 3)],
                    #              [(0 .0, 2.4), (2.7, 2.8), (2.8, 3)],
                    #              [(0.0, 2.4), (2.7, 2.8), (2.8, 3)],
                    #              [(0.0, 2.4), (2.7, 2.8), (2.8, 3)]]
                    # zss = array([ 2.62590637,  2.62589727,  1.34152518,
                    # 2.0221785 ,  0.23706671, 0.2378053])
                    #
                    # tab = [[], [], [(0.0, 2.4)], [(0.0, 2.4)], [(0.0, 2.4)], [(0.0, 2.4)], [(0.0, 2.4)]]
                    #
                    tab = map (lambda x: filter(lambda z: ((z[0]<x[1]) &
                                                           (z[1]>x[1])),x[0]),zip(zinterval,zss))
                    #print tab
                    def findindex(x):
                        if len(x[1])>0:
                            k = x[0].index(x[1][0])+1
                            return(k)
                        else:
                            return(0)
                    #
                    # indexss = [0, 0 , 1 , 1 ,1 ,1]
                    # L.stridess[nstrs] = [ 9 , 9 , 9, 9 , 9 , 9 ]
                    # indexnex = [ 9 , 9 , 10 , 10 , 10 , 10 ]
                    #
                    indexss = np.array(map(findindex,zip(zinterval,tab)))
                    uw = np.where(indexss==0)[0]
                    indexnew = L.stridess[nstrs]+indexss
                    indexnew[uw] = nstrs[uw]
                    #ind  = map(lambda x: np.where(L.lsss==x[0])+x[1],zip(nstrs,indexss))
                    #iindexnex = L.isss[ind]
                    #indexnew = map(lambda x: x[0] if x[1]==0 else 1000000+100*x[0]+x[1]-1,zip(nstrs,indexss))
                    #indexnew = map(lambda x: x[0] if x[1]==0 else 1000000+100*x[0]+x[1]-1,zip(nstrs,indexss))
                    # update signature
                    siges[w] = indexnew
                    #if k==3:
                    #    print siges

                    #if k==1:
                    #    print "indexss:",indexss
                    #    print "indexnew:",indexnew
                    #    print "siges",siges
                    #print siges
                    #pdb.set_trace()
                    #pdb.set_trace()
                    # expand dimension add z dimension (2)
                    # tuple concatenation doesn't work with array this is strange!!
                    #
                    # >> a = (1,2,3)
                    # >> b = (3,5,6)
                    # >> a+b
                    # (1,2,3,3,5,6)
                    # but
                    # >> u = (array([1,2]),array([1,2]))
                    # >> v = (array([2,2]))
                    # >> u + v
                    # array([[3,4],[3,4]])  inconsistent !
                    #
                    #   z --> kl subseg level 
                    #   siges[0,:] --> Ms + nstr *Mss + (kl)
                    #
                try:
                    # r3d[k+Nint]['alpha'] = np.hstack((r3d[k+Nint]['alpha'],a1es))
                    # r3d[k+Nint]['ks'] = np.hstack((r3d[k+Nint]['ks'],ks))
                    r3d[k+Nint]['pt']  = np.dstack((r3d[k+Nint]['pt'], ptees))
                    r3d[k+Nint]['sig'] = np.dstack((r3d[k+Nint]['sig'], siges))
                    r3d[k+Nint]['sig2d'].append(sigsave)
                except:
                    r3d[k+Nint] = {}
                    # r3d[k+Nint]['alpha'] = a1es
                    # r3d[k+Nint]['ks'] = ks
                    r3d[k+Nint]['pt'] = ptees
                    r3d[k+Nint]['sig'] = siges
                    r3d[k+Nint]['sig2d'] = [sigsave]

        #
        # Add Line Of Sight ray information
        #   pt =  [tx,rx]
        #   sig = [0,0]
        #

        if (self.los) & (np.sum(tx-rx,axis=0)<>0):
            r3d[0] = {}
            r3d[0]['sig'] = np.zeros((2,2,1))
            r3d[0]['sig2d'] = np.zeros((2,2,1))
            r3d[0]['pt'] = np.zeros((3,2,1))
            r3d[0]['pt'][:,0,:] = tx[:,np.newaxis]
            r3d[0]['pt'][:,1,:] = rx[:,np.newaxis]

        # r3d.nray = reduce(lambda x,y : y + np.shape(r3d[x]['sig'])[2],lnint)
        # count total number of ray
        # evaluate length of ray segment
        #
        # vsi
        # si
        # dis
        #
        val =0
        for k in r3d.keys():
            nrayk = np.shape(r3d[k]['sig'])[2]
            r3d[k]['nbrays'] = nrayk
            r3d[k]['rayidx'] = np.arange(nrayk)+val
            r3d.nray = r3d.nray + nrayk
            val=r3d[k]['rayidx'][-1]+1
            # 3 : x,y,z
            # i : interaction index
            # r : ray index
            #
            # k : group of interactions index
            #
            v = r3d[k]['pt'][:, 1:, :]-r3d[k]['pt'][:, 0:-1, :]
            lsi = np.sqrt(np.sum(v*v, axis=0))
            rlength = np.sum(lsi,axis=0)
            if (lsi.any()==0):
                pdb.set_trace()
            if not (lsi.all()>0):
                pdb.set_trace()
            #assert(lsi.all()>0)

            if (len(np.where(lsi==0.))==0) :
                pdb.set_trace()

            #
            # sort rays w.r.t their length
            #

            u = np.argsort(rlength)
            r3d[k]['pt']  = r3d[k]['pt'][:,:,u]
            r3d[k]['sig'] = r3d[k]['sig'][:,:,u]
            #r3d[k]['sig2d'] = r3d[k]['sig2d'][:,:,u]
            si = v/lsi             # ndim , nint - 1 , nray

            # vsi : 3 x (i+1) x r
            r3d[k]['vsi'] = si[:,:,u]

            # si : (i+1) x r
            r3d[k]['si']  = lsi[:,u]
            r3d[k]['dis'] = rlength[u]

        r3d.origin_sig_name = self.origin_sig_name
        r3d.Lfilename = L.filename
        r3d.filename = L.filename.split('.')[0] + '_' + str(r3d.nray)
        return(r3d)

    def length(self,typ=2):
        """ calculate length of rays

        Parameters
        ----------

        typ : int
            1 : length of all segments
            2 : accumulated length
        """
        dk = {}
        for k in self:   # for all interaction group k
            # 3 x Ni-1 x Nr
            vk = self[k]['pt'][:,1:,:]-self[k]['pt'][:,0:-1,:]
            d1 = np.sqrt(np.sum(vk*vk,axis=0))
            d2 = np.sum(d1,axis=0)
            if typ==1:
                dk[k] = d1
            if typ==2:
                dk[k] = d2
        return(dk)


    def locbas(self, L):
        """ calculate ray local basis

        Parameters
        ----------

        L : Layout

        """

        #
        # extract normal in np.array
        #

        # nsegment x 3
        norm = np.array(nx.get_node_attributes(L.Gs,'norm').values())

        # nsegment x k
        key = np.array(nx.get_node_attributes(L.Gs,'norm').keys())

        # maximum number for refering to segment
        # not to be confused with a segment number

        nsmax = max(L.Gs.node.keys())

        mapping = np.zeros(nsmax+1, dtype=int)
        mapping[key] = np.arange(len(key), dtype=int)

        #
        # Structutre number : nstr
        #   the structure number is < 0 for points
        #                           > 0 for segments
        # A segment can have several subsegments (until 100)
        #  nstrs is the nstr of the segment if subsegment :
        #  nstr  is the glabal which allows to recover the slab values
        #
        idx = np.array(())
        if self.los:
            idxts = 1
            nbrayt = 1
        else:
            idxts = 0
            nbrayt = 0

        for k in self:
            #
            # k is the number of interactions in the block
            #
            if k <> 0:

                # structure number (segment or point)
                # nstr : i x r
                nstr = self[k]['sig'][0, 1:-1, :]

                # ityp : i x r
                ityp = self[k]['sig'][1, 1:-1, :]

                # nstr of underlying segment
                # position of interaction corresponding to a sub segment
                # print nstr
                #
                # uss : index of subsegment
                # subsegments are not nodes of Gs but have positive nst index
                #

                uss   = np.where(nstr>nsmax)

                # print uss

                nstrs = copy.copy(nstr)
                #
                # if subsegments have been found
                #
                if len(uss)>0:
                    ind   = nstr[uss]-nsmax-1
                    nstrs[uss] = np.array(L.lsss)[ind]
                #    print nstr
                #print nstrs
                #pdb.set_trace()
                nray = np.shape(nstr)[1]

                uwall = np.where((ityp == 2) | (ityp == 3))
                udiff = np.where((ityp == 1))
                ufloor= np.where((ityp == 4))
                uceil = np.where((ityp == 5))

                nstrwall  = nstr[uwall[0], uwall[1]]   # nstr of walls
                nstrswall = nstrs[uwall[0], uwall[1]]   # nstrs of walls

                self[k]['nstrwall']  = nstrwall    # store nstr without subsegment
                self[k]['nstrswall'] = nstrswall   # store nstr with subsegment

                self[k]['norm'] = np.zeros((3, k, nray))   # 3 x int x nray

                #
                # Warning : The following commented line assumes that all the segment number are contiguous
                # self[k]['norm'][:,uwall[0],uwall[1]] = norm[nstrwall-1,:].T
                #

                # norm : 3 x i x r
                #
                # norm name is improper norm is in fact the vector associated to the
                # interaction
                # For the diffraction case the normal is replaced by the unit
                # vector along the wedge directed upward.
                #
                self[k]['norm'][:, uwall[0], uwall[1]] = norm[mapping[nstrswall],:].T
                self[k]['norm'][2, ufloor[0], ufloor[1]] = np.ones(len(ufloor[0]))
                self[k]['norm'][2, uceil[0], uceil[1]] = -np.ones(len(uceil[0]))
                self[k]['norm'][2, udiff[0], udiff[1]] = np.ones(len(udiff[0]))

                normcheck = np.sum(self[k]['norm']*self[k]['norm'],axis=0)
                assert normcheck.all()>0.99,pdb.set_trace()


                # 3 : x,y,z
                # i : interaction index
                # r : ray index
                #
                # k : group of interactions index
                #
                #v = self[k]['pt'][:, 1:, :]-self[k]['pt'][:, 0:-1, :]
                #lsi = np.sqrt(np.sum(v*v, axis=0))
                #if (lsi.any()==0):
                #    pdb.set_trace()
                #assert(lsi.all()>0)
                #if (len(np.where(lsi==0.))==0) :
                #    pdb.set_trace()

                #si = v/lsi             # ndim , nint - 1 , nray

                # si : 3 x (i+1) x r
                si = self[k]['vsi']

                # si : (i+1) x r
                #self[k]['si'] = lsi
                #self[k]['dis'] = np.sum(lsi,axis=0)

                # normal : 3 x i x r
                vn = self[k]['norm']
                # s_in : 3 x i x r
                s_in = si[:, 0:-1, :]

                # s_out : 3 x i x r
                s_out = si[:, 1:, :]

                #
                # AOD (rad)
                #

                # th : ,r
                thd = np.arccos(si[2, 0, :])

                # ph : ,r
                phd = np.arctan2(si[1, 0, :], si[0, 0, :])

                # aod : 2 x r  (radians)
                self[k]['aod'] = np.vstack((thd, phd))

                # eth : 3 x r
                eth = np.array([np.cos(thd) * np.cos(phd),
                               np.cos(thd) * np.sin(phd),
                                -np.sin(thd)])
                # eph : 3 x r
                eph = np.array([-np.sin(phd),
                                np.cos(phd),
                                np.zeros(len(phd))])

                # Bo0 : 3 x 2 x r
                Bo0 = np.concatenate((eth[:, np.newaxis, :],
                                      eph[:, np.newaxis, :]), axis=1)

                self[k]['Bo0'] = np.concatenate((si[:, 0, np.newaxis, :],
                                                 eth[:, np.newaxis, :],
                                                 eph[:, np.newaxis, :]), axis=1)

                #
                # scalar product si . norm
                #
                # vn   : 3 x i x r
                # s_in : 3 x i x r

                #
                # scpr : i x r
                #
                #if k ==5:
                #    pdb.set_trace()
                #scpr = np.sum(vn*s_in, axis=0)
                scpr = np.sum(vn*si[:,0:-1,:], axis=0)
                self[k]['scpr'] = scpr
                self[k]['theta'] = np.arccos(abs(scpr))  # *180/np.pi
                #if k ==5:
                    #print "vsi :",self[k]['vsi'][:,:,43]
                    #print "si :",si[:,0:-1,43]
                    #print "vn :",vn[:,:,43]
                    #print "scpr :",self[k]['scpr'][:,43]
                    #print "scpr :",self[k]['scpr'][:,43]
                    #print "theta :",self[k]['scpr'][:,43]

                #
                # Warning need to handle singular case when s_in // vn
                #
                # w : 3 x i x r
                #
                # Handling channel reciprocity s_in --> -s_in
                #
                #w = np.cross(s_in, vn, axisa=0, axisb=0, axisc=0)
                w = np.cross(-s_in, vn, axisa=0, axisb=0, axisc=0)

                # nw : i x r
                #
                #
                # to do fix the colinear bug
                #
                nw = np.sqrt(np.sum(w*w, axis=0))
                if (nw.any()==0):
                    u = np.where(nw==0)
                    # uv = np.array(filter(lambda x : abs(vn[2,x])>0.99,u))
                    # # reshape information for the filter
                    uu = np.array([u[0],u[1]]).T
                    uv = np.array(filter(lambda x : abs(vn[2,x[:,0],x[:,1]])>0.99,[uu]))
                    uh = np.setdiff1d(uu,uv)
                    try:
                        w[:,uv] = np.array(([1,0,0]))[:,np.newaxis,np.newaxis]
                    except:
                        pass
                    try:
                        w[:,uh] = np.array(([0,0,1]))[:,np.newaxis,np.newaxis]
                    except:
                        pass
                #assert(nw.all()>0), pdb.set_trace()
                wn = w/nw
                # Handling channel reciprocity s_in --> -s_in
                #v = np.cross(wn, s_in, axisa=0, axisb=0, axisc=0)
                v = np.cross(wn, -s_in, axisa=0, axisb=0, axisc=0)

                es_in = np.expand_dims(-s_in, axis=1)
                ew = np.expand_dims(wn, axis=1)
                ev = np.expand_dims(v, axis=1)

                #  Bi 3 x 2 x i x r
                Bi = np.concatenate((ew, ev), axis=1)
                #  self[k]['Bi'] 3 x 3 x i x r 
                self[k]['Bi'] = np.concatenate((es_in,ew,ev),axis=1)

                w = np.cross(s_out, vn, axisa=0, axisb=0, axisc=0)
                wn = w/np.sqrt(np.sum(w*w, axis=0))
                v = np.cross(wn, s_out, axisa=0, axisb=0, axisc=0)

                es_out = np.expand_dims(s_out, axis=1)
                ew = np.expand_dims(wn, axis=1)
                ev = np.expand_dims(v, axis=1)

                #  Bi 3 x 2 x i x r
                Bo = np.concatenate((ew, ev), axis=1)
                #  self[k]['Bo'] 3 x 3 x i x r 
                self[k]['Bo'] = np.concatenate((es_out,ew,ev),axis=1)

                #
                # AOA (rad)
                #

                # th : ,r
                # fix doa/dod reciprocity
                #th = np.arccos(si[2, -1, :])
                tha = np.arccos(-si[2, -1, :])

                # th : ,r
                #ph = np.arctan2(si[1, -1, :], si[0, -1, :])
                pha = np.arctan2(-si[1, -1, :], -si[0, -1, :])

                # aoa : 2 x r  (radians)
                self[k]['aoa'] = np.vstack((tha, pha))
                eth = np.array([np.cos(tha) * np.cos(pha),
                               np.cos(tha) * np.sin(pha),
                                -np.sin(tha)])
                eph = np.array([-np.sin(pha),
                                np.cos(pha),
                                np.zeros(len(pha))])
                # Bo0 : 3 x 2 x r
                BiN = np.concatenate((eth[:, np.newaxis, :],
                                      eph[:, np.newaxis, :]), axis=1)

                self[k]['BiN'] = np.concatenate((-si[:,-1,np.newaxis,:],eth[:,np.newaxis,:],
                                                   eph[:,np.newaxis,:]),axis=1)

                #
                # pasting (Bo0,B,BiN)
                #

                # B : 3 x 2 x i x r

                Bo = np.concatenate((Bo0[:, :, np.newaxis, :], Bo), axis=2)
                Bi = np.concatenate((Bi, BiN[:, :, np.newaxis, :]), axis=2)

                # B : 2 x 2 x i x r

                self[k]['B'] = np.einsum('xv...,xw...->vw...', Bi, Bo)

                #BiN = np.array([si[:,-1,:], eth, eph])    # ndim x 3 x Nray
                #self[k]['BiN']=BiN
                # self[k]['B']=np.sum(self[k]['Bi'][:2,:2,np.newaxis]*self[k]['Bo'][np.newaxis,:2,:2],axis=1)

                ## index creation
                ##################
                # create index for retrieving interactions

                # integer offset : total size idx

                idxts = idxts + idx.size

                idx = idxts + np.arange(ityp.size).reshape(np.shape(ityp),order='F')

                nbray = np.shape(idx)[1]

                self[k]['rays'] = idx
                self[k]['nbrays'] = nbray
                self[k]['rayidx'] = nbrayt + np.arange(nbray)

                # create a numpy array to relate the ray index to its corresponding
                # number of interactions

                ray2nbi = np.ones((nbray))


                try:
                    self.ray2nbi=np.hstack((self.ray2nbi,ray2nbi))
                except:
                    self.ray2nbi=ray2nbi

                self.ray2nbi[self[k]['rayidx']]  = k
                nbrayt = nbrayt + nbray
                self.raypt = self.raypt + self[k]['nbrays']
            # if los exists
            else :
                self[k]['nstrwall'] = np.array(())
                self[k]['norm'] = np.array(())
                si = np.sqrt(np.sum((self[0]['pt'][:,0]-self[0]['pt'][:,1])**2,axis=0))
                self[k]['si'] = np.vstack((si,0.))
                self[k]['vsi'] = (self[0]['pt'][:,1]-self[0]['pt'][:,0])/si
                self[k]['dis'] = np.array((si))

                vsi = self[k]['vsi']
                thd = np.arccos(vsi[2])
                phd = np.arctan2(vsi[1], vsi[0])

                self[k]['aod'] = np.vstack((thd, phd))
                self[k]['Bo0'] = np.array(())
                self[k]['scpr'] = np.array(())
                self[k]['theta'] = np.zeros((1,1))

                #
                # The following derivation of the doa is the actual chosen angle convention
                # Those angles are relative to natural spherical coordinates system in the gcs of the scene.
                #
                # for a LOS path :
                #  tha = pi - thd
                #  pha = phd - pi
                #
                self[k]['aoa'] =  np.vstack((np.pi-thd, phd-np.pi))
                E = np.eye(2)[:,:,np.newaxis,np.newaxis]
                self[k]['B'] = np.dstack((E,E))
                ze = np.array([0])
                self[k]['rays'] = np.array(([[0]]))
                self[k]['nbrays'] = 1
                self[k]['rayidx'] = ze
                self.raypt = 1
                self.ray2nbi = ze

            self.isbased=True

    def fillinter(self,L,append=False):
        """  fill ray interactions

        Parameters
        ----------

        L      : Layout
        append : Boolean
            If True append new rays to existing structure


        Returns
        -------

        Update self.I , self.B , self.I0

        """

        # reinitilized ray pointer if not in append mode
        if not append:
            self.raypt = 0

        # stacked interactions
        I = Interactions()

        # rotation basis
        B  = IntB()
        B0 = IntB()

        # LOS Interaction
        Los = IntL()

        # Reflexion
        R = IntR()

        # transmission
        T = IntT()

        # Diffraction
        D = IntD()

        idx = np.array(())

        if self.los:
            idxts = 1
            nbrayt = 1
        else:
            idxts = 0
            nbrayt = 0


        # Transform dictionnary of slab name to array
        # slv = nx.get_node_attributes(L.Gs, "name").values()
        # slk = nx.get_node_attributes(L.Gs, "name").keys()
        # find all material used in simulation
        uslv = np.unique(L.sla[1:])
        #
        # add CEIL and FLOOR
        #
        uslv = np.hstack((uslv, np.array(('CEIL', 'FLOOR'))))

        # create reverse dictionnary with all material as a key
        # and associated point/segment as a value

        #dsla = {}
        #for s in uslv:
        #    dsla[s] = np.where(s == np.array(slv))[0]

        nsmax = max(L.Gs.node.keys())
        #sla = np.zeros((nsmax+1), dtype='S20')

        # array type str with more than 1 character
        # warning use zeros instead of empty because slab zero
        # is virtually used before assigning correct slab to ceil and floor

        #
        # sla is an array of string.
        # each value of Gs node is the index of the corresponding slab
        #

        #sla[slk] = np.array(slv)

        R.dusl = dict.fromkeys(uslv, np.array((), dtype=int))
        T.dusl = dict.fromkeys(uslv, np.array((), dtype=int))

        tsl = np.array(())
        rsl = np.array(())

        # loop on group of interactions
        for k in self:

            if k !=0:

                uR = uT = uD = uRf = uRc = 0.

                # structure number (segment or point)
                # nstr : i x r
                nstr = self[k]['sig'][0, 1:-1, :]

                # ityp : i x r
                ityp = self[k]['sig'][1, 1:-1, :]

                # theta : i x r   ( related to interactions )
                theta = self[k]['theta']

                # (i+1) x r
                si = self[k]['si']

                ## flatten information
                ######################

                # flatten nstr (1 dimension)
                # size1 = i x r
                size1 = nstr.size

                # flatten ityp (method faster than np.ravel() )
                nstrf = np.reshape(nstr,size1,order='F')
                itypf = ityp.reshape(size1,order='F')
                thetaf = theta.reshape(size1,order='F')
                #sif = si[0, :, :].reshape(si[0, :, :].size)

                ## index creation
                ##################
                # create index for retrieving interactions

                # integer offset : total size idx

                idxts = idxts + idx.size

                idx = idxts + np.arange(ityp.size).reshape(np.shape(ityp),order='F')

                nbray = np.shape(idx)[1]

                self[k]['rays'] = idx
                self[k]['nbrays'] = nbray
                self[k]['rayidx'] = nbrayt + np.arange(nbray)

                # create a numpy array to relate the ray index to its corresponding
                # number of interactions

                # ray2nbi = np.ones((nbray))

                #try:
                #    self.ray2nbi=np.hstack((self.ray2nbi,ray2nbi))
                #except:
                #    self.ray2nbi=ray2nbi

                #self.ray2nbi[self[k]['rayidx']]  = k
                nbrayt = nbrayt + nbray
                #self.raypt = self.raypt + self[k]['nbrays']

                idxf = idx.reshape(idx.size,order='F')
                #  (i+1)xr
                size2 = si[:, :].size
                #  ,(i+1)xr
                sif = si[:, :].reshape(size2,order='F')
                # 2x2,(i+1)xr

                #
                # self[k]['B'] 2 x 2 x i x r
                #
                # first unitary matrix (2x2xr)
                b0 = self[k]['B'][:,:,0,:]
                # first unitary matrix 1:
                # dimension i and r are merged
                b  = self[k]['B'][:,:,1:,:].reshape(2, 2, size2-nbray,order='F')

                ## find used slab
                ##################
                # find slab type for the rnstr
                # nstrf is a number of slab
                # this is a problem for handling subsegment
                #

                sl = L.sla[nstrf]

                # seek for interactions position
                ################################

                uD = np.where((itypf == 1))[0]
                uR = np.where((itypf == 2))[0]
                uT = np.where((itypf == 3))[0]
                uRf = np.where((itypf == 4))[0]
                uRc = np.where((itypf == 5))[0]

                # assign floor and ceil slab
                ############################

                # WARNING
                # in future version floor and ceil could be different for each cycle.
                # this information would be directly obtained from L.Gs
                # then the two following lines would have to be  modified

                sl[uRf] = 'FLOOR'
                sl[uRc] = 'CEIL'

                # Fill the used slab
                #####################

                tsl = np.hstack((tsl, sl[uT]))
                rsl = np.hstack((rsl, sl[uR], sl[uRf], sl[uRc]))

    ##            for s in uslv:
    ##
    ##                T.dusl[s]=np.hstack((T.dusl[s],len(T.idx) + np.where(sl[uT]==s)[0]))
    ##                R.dusl[s]=np.hstack((R.dusl[s],len(R.idx) + np.where(sl[uR]==s)[0]))
    ##            R.dusl['FLOOR']=np.hstack((R.dusl['FLOOR'],len(R.idx)+len(uR) + np.where(sl[uRf]=='FLOOR')[0]))
    # R.dusl['CEIL']=np.hstack((R.dusl['CEIL'],len(R.idx)+len(uR)+len(uRf) +
    # np.where(sl[uRc]=='CEIL')[0]))

                # Basis
                # Hugr issue with B index
                # Friedman version Bs was entering in the index
                # maybe B can have the same index that interactions
                # but this must be managed when evaluation of CIR is made

                # BU 10/4/2013
                # .. todo:  This is no longer idxf the good index
                # why the transposition b is first 2x2x(i+1)xr
                #                             idxf is (ixr)
                #
                # need to check how B is used in eval()
                #
                # Warning
                # -------
                # B.idx refers to an interaction index
                # whereas B0.idx refers to a ray number

                B.stack(data=b.T, idx=idxf)
                B0.stack(data=b0.T,idx=self[k]['rayidx'])

                ### Reflexion
                ############
                ### wall reflexion
                R.stack(data=np.array((thetaf[uR], sif[uR], sif[uR+1])).T,
                        idx=idxf[uR])
                # floor reflexion
                R.stack(data=np.array((thetaf[uRf], sif[uRf], sif[uRf+1])).T,
                        idx=idxf[uRf])
                # ceil reflexion
                R.stack(data=np.array((thetaf[uRc], sif[uRc], sif[uRc+1])).T,
                        idx=idxf[uRc])

                ### sl[idxf[uT]]
                # Transmision
                ############
                T.stack(data=np.array((thetaf[uT], sif[uT], sif[uT+1])).T, idx=idxf[uT])

            elif self.los:
                ze = np.array([0])
                #self[k]['rays'] = np.array(([[0]]))
                #self[k]['nbrays'] = 1
                #self[k]['rayidx'] = ze
                #self.raypt = 1
                #self.ray2nbi=ze
                B.stack(data=np.eye(2)[np.newaxis,:,:], idx=ze)
                B0.stack(data=np.eye(2)[np.newaxis,:,:],idx=ze)

        T.create_dusl(tsl)
        R.create_dusl(rsl)

        # create interactions structure
        self.I = I
        self.I.add([T, R])
        # create rotation base B
        self.B = B
        # create rotation base B0
        self.B0 = B0

        self.filled = True

    def eval(self,fGHz=np.array([2.4]),ib=[]):
        """  docstring for eval

        Parameters
        ----------

        fGHz : array
            frequency in GHz array
        ib : list of intercation block

        """

        #print 'Rays evaluation'
        self.fGHz=fGHz
        # evaluation of interaction
        self.I.eval(fGHz)
        # evaluation of base B  (2x2)
        # B and B0 do no depend on frequency
        # just an axis extension (np.newaxis)
        #pdb.set_trace()

        # 1 x i x 2 x 2
        B  = self.B.data[np.newaxis,...]
        # 1 x r x 2 x 2
        B0 = self.B0.data[np.newaxis,...]

        # Ct : f x r x 2 x 2
        Ct = np.zeros((self.I.nf, self.nray, 2, 2), dtype=complex)

        # delays : ,r
        self.delays = np.zeros((self.nray))

        # dis : ,r
        self.dis = np.zeros((self.nray))

        #nf : number of frequency point
        nf = self.I.nf

        aod= np.empty((2,self.nray))
        aoa= np.empty((2,self.nray))
        # loop on interaction blocks
        if ib==[]:
            ib=self.keys()

        for l in ib:
            # ir : ray index
            ir = self[l]['rayidx']
            aoa[:,ir]=self[l]['aoa']
            aod[:,ir]=self[l]['aod']
            if l != 0:
                # l stands for the number of interactions
                r = self[l]['nbrays']
                # reshape in order to have a 1D list of index
                # reshape ray index
                rrl = self[l]['rays'].reshape(r*l,order='F')
                # get the corresponding evaluated interactions
                #
                # reshape error can be tricky to debug.
                #
                # f , r , l , 2 , 2
                A = self.I.I[:, rrl, :, :].reshape(self.I.nf, r, l, 2, 2)
                # get the corresponding unitary matrix B
                # 1 , r , l , 2 , 2
                #Bl = B[:, rrl, :, :].reshape(self.I.nf, r, l, 2, 2,order='F')
                Bl = B[:, rrl, :, :].reshape(1, r, l, 2, 2)
                # get the first uitary matrix B0l
                B0l = B0[:,ir,:, :]
                # get alpha
                # alpha = self.I.alpha[rrl].reshape(r, l,order='F')
                # # get gamma
                # gamma = self.I.gamma[rrl].reshape(r, l,order='F')
                # # get si0
                # si0 = self.I.si0[rrl].reshape(r, l,order='F')
                # # get sout
                # sout = self.I.sout[rrl].reshape(r, l,order='F')

                try:
                    del Z
                except:
                    pass


                #print "\nrays",ir
                #print "-----------------------"
                ## loop on all the interactions of ray with l interactions
                for i in range(0, l):


    ############################################
    ##                # Divergence factor D
    ###                 not yet implementented
    ############################################
    #                if i == 0:
    #                    D0=1./si0[:,1]
    #                    rho1=si0[:,1]*alpha[:,i]
    #                    rho2=si0[:,1]*alpha[:,i]*gamma[:,i]
    #                    D=np.sqrt(
    #                     ( (rho1 ) / (rho1 + sout[:,i]) )
    #                     *( (rho2) / (rho2 + sout[:,i])))
    #                    D=D*D0
    #                    rho1=rho1+(sout[:,i]*alpha[:,i])
    #                    rho2=rho2+(sout[:,i]*alpha[:,i]*gamma[:,i])

    ##                     gerer le loss
    #                    if np.isnan(D).any():
    #                        p=np.nonzero(np.isnan(D))[0]
    #                        D[p]=1./sout[p,1]
    #                else :
    #                    D=np.sqrt(
    #                     ( (rho1 ) / (rho1 + sout[:,i]) )
    #                     *( (rho2) / (rho2 + sout[:,i])))

    #                    rho1=rho1+(sout[:,i]*alpha[:,i])
    #                    rho2=rho2+(sout[:,i]*alpha[:,i]*gamma[:,i])
    ############################################

                    #  A0  (X dot Y)
                    #  |    |     |
                    #  v    v     v
                    ##########################
                    ## B  # I  # B  # I  # B #
                    ##########################
                    #      \_____/   \______/
                    #         |         |
                    #       Atmp(i)   Atmp(i+1)
                    #
                    # Z=Atmp(i) dot Atmp(i+1)

                    #X = A [:, :, i, :, :]
                    #Y = Bl[:, :, i, :, :]
                    # pdb.set_trace()
                    if i == 0:
                    ## First Basis added
                        Atmp = A[:, :, i, :, :]
                        B00 = B0l[:, :,  :, :]
                        Z = np.sum(Atmp[..., :, :, np.newaxis]
                                  *B00[..., np.newaxis, :, :], axis=-2)
                    else:
                        Atmp = A[:, :, i, :, :]
                        BB = Bl[:, :, i-1, :, :]
                        Ztmp = np.sum(Atmp[..., :, :, np.newaxis]
                                  *BB[..., np.newaxis, :, :], axis=-2)


                        Z = np.sum(Ztmp[..., :, :, np.newaxis]
                                  *Z[..., np.newaxis, :, :], axis=-2)

                    if i == l-1:
                        BB = Bl[:, :, i, :, :]
                        Z = np.sum(BB[..., :, :, np.newaxis]
                                  *Z[..., np.newaxis, :, :], axis=-2)


                # fill the C tilde MDA

                Ct[:,ir, :, :] = Z[:, :, :, :]

                # delay computation:
                # sum the distance from antenna to first interaction si0
                # and the sum of all outgoing segments
                #self[l]['dis'] = self.I.si0[self[l]['rays'][0,:]] \
                #        + np.sum(self.I.sout[self[l]['rays']], axis=0)

                # attenuation due to distance
                # will be removed once the divergence factor will be implemented
                Ct[:,ir, :, :] = Ct[:, ir, :, :]*1./(self[l]['dis'][np.newaxis, :, np.newaxis, np.newaxis])
                self.delays[ir] = self[l]['dis']/0.3
                self.dis[ir] = self[l]['dis']
        #
        # true LOS when no interaction
        #
        if self.los:
            Ct[:,0, :, :]= np.eye(2,2)[np.newaxis,np.newaxis,:,:]
            #self[0]['dis'] = self[0]['si'][0]
            # Fris
            Ct[:,0, :, :] = Ct[:,0, :, :]*1./(self[0]['dis'][np.newaxis, :, np.newaxis, np.newaxis])
            self.delays[0] = self[0]['dis']/0.3
            self.dis[0] = self[0]['dis']


        # To be corrected in a future version

        Ct = np.swapaxes(Ct, 1, 0)

        c11 = Ct[:,:,0,0]
        c12 = Ct[:,:,0,1]
        c21 = Ct[:,:,1,0]
        c22 = Ct[:,:,1,1]


        #
        # Construction of the Ctilde channel
        #
        Cn = Ctilde()
        Cn.Cpp = bs.FUsignal(self.I.fGHz, c11)
        Cn.Cpt = bs.FUsignal(self.I.fGHz, c12)
        Cn.Ctp = bs.FUsignal(self.I.fGHz, c21)
        Cn.Ctt = bs.FUsignal(self.I.fGHz, c22)
        Cn.nfreq = self.I.nf
        Cn.nray = self.nray
        Cn.tauk = self.delays
        Cn.fGHz = self.I.fGHz
        # r x 2
        Cn.tang = aod.T
        Cn.tangl = aod.T
        # r x 2
        Cn.rang = aoa.T
        Cn.rangl = aoa.T
        # add aoa and aod

        self.evaluated = True

        return(Cn)


    def ray(self, r):
        """ returns the index of interactions of r

        Parameters
        ----------

        r : integer
            ray index

        Returns
        -------

        ir : index of interactions of r


        """
        raypos = np.nonzero(self[self.ray2nbi[r]]['rayidx'] == r)[0]
        return(self[self.ray2nbi[r]]['rays'][:,raypos][:,0])


    def typ(self, r):
        """ returns interactions list type of a given ray

        Parameters
        ----------

        r : integer
            ray index
        """

        a = self.ray(r)
        return(self.I.typ[a])

    def info(self, r,ifGHz=0):
        """ provides information for a given ray r

        Parameters
        ----------

        r : int
            ray index
        ifGHz : int
            frequency index

        """

        if self.evaluated:
            print '-------------------------'
            print 'Informations of ray #', r
            print '-------------------------\n'

            ray = self.ray(r)
            typ = self.typ(r)
            print '{0:5} , {1:4}, {2:10}, {3:7}, {4:10}, {5:10}'.format('Index',
                                                                        'type',
                                                                        'slab', 'th(rad)', 'alpha', 'gamma2')
            print '{0:5} , {1:4}, {2:10}, {3:7.2}, {4:10.2}, {5:10.2}'.format(r, 'B0', '-', '-', '-', '-')
            for iidx, i in enumerate(typ):
                if i == 'T' or i == 'R':
                    I = getattr(self.I, i)
                    for slab in I.dusl.keys():
    #                    print slab
                        midx = I.dusl[slab]
    #                    print midx
                        Iidx = np.array((I.idx))[midx]
                        th = I.data[I.dusl[slab], 0]
                        gamma = I.gamma[midx]
                        alpha = I.alpha[midx]
                        for ii, Ii in enumerate(Iidx):
                            if Ii == ray[iidx]:
                                print '{0:5} , {1:4}, {2:10}, {3:7.2}, {4:10.2}, {5:10.2}'.format(Ii, i, slab, th[ii], alpha[ii], gamma[ii])

                # else:
                print '{0:5} , {1:4}, {2:10}, {3:7.2}, {4:10.2}, {5:10.2}'.format(ray[iidx], 'B', '-', '-', '-', '-')
                #              print '{0:5} , {1:4}, {2:10}, {3:7}, {4:10}, {5:10}'.format(ray[iidx], i, '-', '-', '-', '-')

            print '\n----------------------------------------'
            print ' Matrix of ray #', r, 'at f=', self.I.fGHz[ifGHz]
            print '----------------------------------------'

            print 'rotation matrix#', 'type: B0'
            print self.B0.data[r,:,:]
            for iidx, i in enumerate(typ):
                print 'interaction #', ray[iidx], 'type:', i
                # f x l x 2 x 2
                print self.I.I[ifGHz, ray[iidx], :, :]
                print 'rotation matrix#',[ray[iidx]], 'type: B'
                print self.B.data[ray[iidx], :, :]
        else:
            print 'Rays have not been evaluated yet'

    def signature(self, ni ,nr):
        """ extract ray signature

        Parameters
        ----------

        ni : int
        nr : int

        Returns
        -------

        sig : ndarray

        Notes
        -----

        Signature of a ray is store as a member

        r[nint]['sig']

        """
        sig = self[ni]['sig'][:,:,nr]
        return(sig)

    def show3d(self,
               ray,
               bdis=True,
               bbas=False,
               bstruc=True,
               col=np.array([1, 0, 1]),
               id=0,
               linewidth=1):
        """ plot a set of 3D rays

        Parameters
        ----------

        ray :
        block : int
            interaction block
        bdis : Boolean
            if False return .vect filename (True)
        bbas : Boolean
            display local basis (False)
        bstruc : Boolean
            display structure (True)
        col  : ndarray() 1x3
            color of the ray  ([1,0,1])
        id   : Integer
            id of the ray (default 0)
        linewidth : Integer
            default 1

        """

        filerac = pyu.getlong("ray" + str(id), pstruc['DIRGEOM'])
        _filerac = pyu.getshort(filerac)
        filename_list = filerac + '.list'
        filename_vect = filerac + '.vect'

        try:
            fo = open(filename_vect, "w")
        except:
            raise NameError(filename)

        fo.write("appearance { linewidth %d }\n" % linewidth)

        fo.write("VECT\n")

        fo.write("1 %d 1\n\n" % len(ray[0, :]))
        fo.write("%d\n" % len(ray[0, :]))
        fo.write("1\n")
        for i in range(len(ray[0, :])):
            fo.write("%g %g %g\n" % (ray[0, i], ray[1, i],
                                     ray[2, i]))
        # fo.write("%d %d %d 0\n" % (col[0],col[1],col[2]))
        fo.write("%g %g %g 0\n" % (col[0], col[1], col[2]))
        fo.close()

        #
        # Ajout des bases locales
        #


        fo = open(filename_list, "w")
        fo.write("LIST\n")
        fo.write("{<" + filename_vect + "}\n")
        if (bstruc):
            # fo.write("{<strucTxRx.off}\n")
            fo.write("{<" + _filestr + ".off}\n")
        filename = filename_list
        fo.close()

        if (bdis):
        #
        # Geomview Visualisation
        #
            chaine = "geomview -nopanel -b 1 1 1 " + filename + \
                " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)

    @mlab.show
    def _show3(self,L=[],ilist=[],rlist=[],newfig=False):
        """ plot 3D rays in environment using Mayavi

        Parameters
        ----------

        L : Layout object
            Layout to be displayed
        ilist : list 
            list of group of interactions
        rlist : list
            list of index rays
        newfig : boolean (default: False)
            if true create a new mayavi figure
            else : use the current

        """
        if newfig:
            mlab.clf()
            f = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
        else :
            f = mlab.gcf()


        if L != []:
            try:
                L.filename
            except:
                raise NameError('L argument must be a layout object')

            L._show3()

        if ilist == []:
            nbi = self.keys()
        else :
            if not isinstance(ilist,list):
                nbi = np.array(([ilist]))
            else :
                nbi = np.array((ilist))

        if not isinstance(rlist,list):
            rlist=[rlist]

        for i in nbi:
            if rlist == []:
                r = range(np.shape(self[i]['pt'])[2])
            else :
                r = np.array((rlist))

            # number of rays
            nbr = len(r) 
            # current number of interactions
            cnbi = i + 2
            pt = self[i]['pt'][:,:,r].reshape(3,cnbi*nbr,order='F')
            # lines = np.arange(cnbi*nbr).reshape(cnbi,nbr)
            lines = np.arange(cnbi*nbr).reshape(nbr,cnbi)
            mesh = tvtk.PolyData(points=pt.T, polys=lines)
            mlab.pipeline.surface(mlab.pipeline.extract_edges(mesh),
                                                 color=(0, 0, 0), )
            f.children[-1].name='Rays with ' + str(i) + 'interactions'

    def show3(self,
              L=[],
              bdis=True,
              bstruc=True,
              bbasi = False,
              bbaso = False,
              id=0,
              ilist=[],
              raylist=[],centered=True):
        """ plot 3D rays within the simulated environment

        Parameters
        ----------

        bdis : boolean
            True
        bstruc : boolean
            True
        bbasi : boolean
            display input basis of each interaction of rays
        bbaso : boolean
            display ouput basis of each interaction of rays
        id : int
        L : Layout object
            Layout to be displayed
        ilist : list of group of interactions
        raylist : list of index rays
        centered : boolean    
            if True center the layout before display
        

        """

        try:
            L.filename
        except:
            raise NameError('L argument must be a layout object')

        if not centered:
            pg=np.array([[0],[0],[0]])

        strucname= L.filename.split('.')[0]
        pg = L.geomfile(centered=centered)
        pg = np.hstack((pg,0.)).reshape(3,1)

        if ilist == []:
            ilist = self.keys()
        pTx = self.pTx.reshape((3, 1))-pg
        pRx = self.pRx.reshape((3, 1))-pg
        filename = pyu.getlong("grRay" + str(id) + ".list", pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        if bstruc:
            fo.write("{<"+strucname+".off}\n")
            if bbasi:
                if not self.isbased:
                    raise NameError('Bases have not been computed (self.locbas(Layout)')
                else:   
                    base_listi = geu.Geomlist('baselisti',clear=True)
                    base_listi.append("LIST\n")
            if bbaso:
                if not self.isbased:
                    raise NameError('Bases have not been computed (self.locbas(Layout)')
                else:   
                    base_listo = geu.Geomlist('baselisto',clear=True)
                    base_listo.append("LIST\n")

            # fo.write("{<strucTxRx.off}\n")

            k = 0
            for i in ilist:
                if raylist == []:
                    rlist = range(np.shape(self[i]['pt'])[2])
                else:
                    rlist = raylist
                for j in rlist:
                    ray = np.hstack((pTx,np.hstack((self[i]['pt'][:, :, j]-pg, pRx))))
                    # ray = rays[i]['pt'][:,:,j]
                    col = np.array([0, 0, 0])
                    # print ray
                    fileray = self.show3d(ray=ray, bdis=False,
                                          bstruc=False, col=col, id=k)
                    k += 1
                    fo.write("{< " + fileray + " }\n")
                    if bbasi:
                        for inter in range(i):
                            filebi = 'bi_' + str(j) + '_' + str(i) + '_' +str(inter)
                            basi = geu.GeomVect(filebi)
                            basi.geomBase(self[i]['Bi'][:,:,inter,j],pt=self[i]['pt'][:,inter+1,j]-pg[:,0])
                            base_listi.append("{<" + filebi +'.vect' "}\n")
                        filebi = 'bi_' + str(j) + '_' + str(i) + '_' +str(inter-1)
                        basi = geu.GeomVect(filebi)
                        basi.geomBase(self[i]['BiN'][:,:,j],pt=self[i]['pt'][:,-1,j]-pg[:,0])
                        base_listi.append("{<" + filebi +'.vect' "}\n")
                    if bbaso:
                        for inter in range(i):
                            filebo = 'bo_' + str(j) + '_' + str(i) + '_' +str(inter)
                            baso = geu.GeomVect(filebo)
                            baso.geomBase(self[i]['Bo'][:,:,inter,j],pt=self[i]['pt'][:,inter+1,j]-pg[:,0])
                            base_listo.append("{<" + filebo +'.vect' "}\n")
                        filebo = 'bo_' + str(j) + '_' + str(i) + '_' +str(inter+1)
                        baso = geu.GeomVect(filebo)
                        baso.geomBase(self[i]['Bo0'][:,:,j],pt=self[i]['pt'][:,0,j]-pg[:,0])
                        base_listo.append("{<" + filebo +'.vect' "}\n")
            if bbasi:
                fo.write("{< " + "baselisti.list}\n")
            if bbaso:  
                fo.write("{< " + "baselisto.list}\n")

        fo.close()
        if (bdis):
            chaine = "geomview " + filename + " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)
   # def show3(self,
   #            bdis=True,
   #            bstruc=True,
   #            id=0,
   #            strucname='defstr',
   #            ilist=[],
   #            raylist=[],pg=np.array([[0],[0],[0]])):
   #      """ plot 3D rays within the simulated environment

   #      Parameters
   #      ----------

   #      bdis : boolean
   #          True
   #      bstruc : boolean
   #          True
   #      id : int
   #      strucname : string
   #          'defstr'
   #      ilist : list of group of interactions
   #      raylist : list of index rays
   #      pg : centroid of the structure
        

   #      """
   #      if ilist == []:
   #          ilist = self.keys()
   #      pTx = self.pTx.reshape((3, 1))-pg
   #      pRx = self.pRx.reshape((3, 1))-pg
   #      filename = pyu.getlong("grRay" + str(id) + ".list", pstruc['DIRGEOM'])
   #      fo = open(filename, "w")
   #      fo.write("LIST\n")
   #      if bstruc:
   #          fo.write("{<"+strucname+".off}\n")
   #          # fo.write("{<strucTxRx.off}\n")
   #          k = 0
   #          for i in ilist:
   #              if raylist == []:
   #                  rlist = range(np.shape(self[i]['pt'])[2])
   #              else:
   #                  rlist = raylist
   #              for j in rlist:
   #                  ray = np.hstack((pTx,np.hstack((self[i]['pt'][:, :, j]-pg, pRx))))
   #                  # ray = rays[i]['pt'][:,:,j]
   #                  col = np.array([2, 0, 1])
   #                  # print ray
   #                  fileray = self.show3d(ray=ray, bdis=False,
   #                                        bstruc=False, col=col, id=k)
   #                  k += 1
   #                  fo.write("{< " + fileray + " }\n")
   #      fo.close()
   #      if (bdis):
   #          chaine = "geomview " + filename + " 2>/dev/null &"
   #          os.system(chaine)
   #      else:
   #          return(filename)


if __name__ == "__main__":
    doctest.testmod()
