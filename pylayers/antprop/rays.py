# -*- coding: latin1 -*-
from __future__ import print_function
"""
.. currentmodule:: pylayers.antprop.rays

.. autosummary::
    :members:

"""
import doctest
import os
import sys
try:
    from mayavi import mlab
except:
    print('Layout:Mayavi is not installed')
import pdb
import copy
if sys.version_info.major==2:
    import ConfigParser
else:
    import configparser
import networkx as nx
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import struct as stru
import pylayers.util.geomutil as geu
import pylayers.util.pyutil as pyu
from pylayers.util.project import *
from pylayers.antprop.interactions import *
from pylayers.antprop.slab import *
from pylayers.antprop.channel import Ctilde
from pylayers.gis.layout import Layout
import pylayers.signal.bsignal as bs
import shapely.geometry as shg
import h5py
import operator


class Rays(PyLayers, dict):
    """ Class handling a set of rays

    Attributes
    ----------

    pTx  : np.array
        transmitter (3,)
    pRx  : np.array
        receiver (3,)
    B    : IntB
    B0   : IntB
    I    : Interactions
    I.I  : np.array
        (f,nI,3,3)
    I.T  : IntT
    I.T.A : np.array
        (f,iT,3,3)
    I.R  : IntR
    I.R.A : np.array
        (f,iR,3,3)
    I.D  : IntD
    I.D.A : np.array
        (f,iD,3,3)
    Lfilename : string
        Layout name
    delays : np.array
        ray delays
    dis : np.array
        ray distance = delays*0.3
    nray : int
        number of rays
    evaluated : boolean
        are rays evaluated ?
    is3D : boolean
        are rays 2d or 3d rays ?
    isbased : boolean
        locbas has been applied ?
    filles : boolean
        filled  has been applied ?
    los : boolean
        Line of sight boolean
    fGHz : np.array
        frequency points for evaluation
    origin_sig_name : string
        signature file which produces the rays


    Notes
    -----

    The Rays object is obtained from a signature.
    It is a container for a set of rays between a source
    and a destination point defining a radio link.

    Once a Rays object has been obtained in 2D, it is transformed
    in 3D via the **to3D** method. This method takes two parameters :
    the height from floor to ceil, and the number N of
    multiple reflections to account for.

    Once the 3d rays have been calculated,
    the local basis are evaluated along those rays. This is
    done through the **locbas** method.

    Once the local basis have been calculated the different
    interactions along rays can be informed via the **fillinter**
    method.

    Once the interactions are informed, the field along rays can
    be evaluated via the **eval** method

    to3D -> fillinter -> eval

    """
    def __init__(self, pTx, pRx):
        """ object constructor

        Parameters
        ----------

        pTx : np.array
            transmitter coordinates
        pRx : np.array
            receiver coordinates

        """

        self.pTx = pTx
        self.pRx = pRx
        self.nray = 0
        self.nray2D = 0
        self.raypt = 0
        self.isdirect = False
        self.los = False
        self.is3D = False
        self.isbased = False
        self.filled = False
        self.evaluated = False

    def __len__(self):
        Nray = 0
        for k in self.keys():
            sh = np.shape(self[k]['sig'])
            Nray = Nray + sh[2]
        return Nray


    # def __add__(self,r):

    #     if (not r.is3D) and (not r.isbased) and (not self.is3D) and (not self.isbased) :
    #         raise AttributeError('both Ray structures must be 3D and based to be added')


    #     for ni in r:
    #         if self.has_key(ni):
    #             import ipdb
    #             ipdb.set_trace()
    #             # check if som rays already exists
    #             # if so, don't add them
    #             lur = np.array([])
    #             for ur in range(self[ni]['pt'].shape[2]):
    #                 udifferent = np.where(np.all(np.all(r[ni]['pt'][...,ur][...,None]!=self[ni]['pt'],axis=0),axis=0))[0]
    #                 lur = np.hstack((lur,udifferent ))
    #                 import ipdb
    #                 ipdb.set_trace()

    #             self[ni]['pt'] = np.concatenate((self[ni]['pt'],r[ni]['pt']),axis=2)
    #             self[ni]['sig'] = np.concatenate((self[ni]['sig'],r[ni]['sig']),axis=2)
    #             self[ni]['si'] = np.concatenate((self[ni]['si'],r[ni]['si']),axis=1)
    #             self[ni]['rayidx'] = np.concatenate((self[ni]['rayidx'],r[ni]['rayidx']),axis=0)
    #             self[ni]['dis'] = np.concatenate((self[ni]['dis'],r[ni]['dis']),axis=0)
    #             self[ni]['vsi'] = np.concatenate((self[ni]['vsi'],r[ni]['vsi']),axis=1)
    #             self[ni]['nbrays'] += 1
    #             if ni != 0:
    #                 self[ni]['BiN'] = np.concatenate((self[ni]['BiN'],r[ni]['BiN']),axis=2)
    #                 self[ni]['Bi'] = np.concatenate((self[ni]['Bi'],r[ni]['Bi']),axis=3)
    #                 self[ni]['Bo'] = np.concatenate((self[ni]['Bo'],r[ni]['Bo']),axis=3)
    #                 self[ni]['Bo0'] = np.concatenate((self[ni]['Bo0'],r[ni]['Bo0']),axis=2)
    #                 self[ni]['scpr'] = np.concatenate((self[ni]['scpr'],r[ni]['scpr']),axis=1)
    #                 self[ni]['norm'] = np.concatenate((self[ni]['norm'],r[ni]['norm']),axis=2)

    #             self[ni]['B'] = np.concatenate((self[ni]['B'],r[ni]['B']),axis=3)
    #             self[ni]['aod'] = np.concatenate((self[ni]['aod'],r[ni]['aod']),axis=1)
    #             self[ni]['aoa'] = np.concatenate((self[ni]['aoa'],r[ni]['aoa']),axis=1)
    #             self[ni]['theta'] = np.concatenate((self[ni]['theta'],r[ni]['theta']),axis=1)

    #             if r[ni].has_key('diffidx'):
    #                 if self[ni].has_key('diffidx'):
    #                     self[ni]['diffidx'] = np.concatenate((self[ni]['diffidx'],r[ni]['diffidx']))
    #                     self[ni]['diffvect'] = np.concatenate((self[ni]['diffvect'],r[ni]['diffvect']),axis=1)
    #                     self[ni]['diffslabs'].append(r[ni]['diffslabs'])
                        
    #                 else:
    #                     self[ni]['diffidx'] = r['diffidx']
    #                     self[ni]['diffvect'] = r['diffvect']
    #                     self[ni]['diffslabs'] = r['diffslabs']

    #         else:
    #             self[ni]=r[ni]


    def __repr__(self):
        s = ''
        ni = 0
        nl = 0
        lgi = list(self.keys())
        lgi.sort()
        if self.is3D:
            s = self.__class__.__name__ + '3D\n' + '----------'+'\n'

            for k in lgi:
                r = self[k]['rayidx']
                nr = len(r)
                s = s + str(k)+' / '+str(nr)+ ' : '+str(r)+'\n'
                ni = ni + nr*k
                nl = nl + nr*(2*k+1)
            nray2D = self.nray2D
        else:
            s = self.__class__.__name__ + '2D\n' + '----------'+'\n'
            nray2D = len(self)

        if self.los:
            s = s + "LOS "
        if self.isbased:
            s = s + "based "
        if self.filled:
            s = s + "filled "

        s = s + '\n'
        s = s + 'N2Drays : '+ str(nray2D) + '\n'

        if hasattr(self,'nb_origin_sig'):
            s = s + 'from '+ str(self.nb_origin_sig) + ' signatures\n'
            s = s + '#Rays/#Sig: '+ str(nray2D/(1.*self.nb_origin_sig) )

        s = s + '\npTx : '+ str(self.pTx) + '\npRx : ' + str(self.pRx)+'\n'

        if not self.is3D:
            ray_cpt = 0
            for k in lgi:
                s = s + str(k) + ':\n'
                if len(self[k]['sig'])>0:
                    sig = self[k]['sig'][0,:]
                    sha0 = sig.shape[0]
                    sha1 = sig.shape[1]
                    for l in np.arange(sha1):
                        s = s + '  ' + str(ray_cpt) + ':'
                        ray_cpt +=1
                        for n in np.arange(sha0):
                            s = s + '      ' + str(sig[n,l])
                        s = s + '\n'

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

        filename = self.filename+'_' + str(idx)
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
        print(filenameh5)


    def loadh5(self,filename=[], idx=0):
        """ load rays hdf5 format

        Parameters
        ----------

        idx : int

        """
        if filename == []:
            filenameh5 = self.filename + '_' + str(idx) + '.h5'
        else :
            filenameh5 = filename

        filename=pyu.getlong(filenameh5,pstruc['DIRR3D'])
        print(filename)

        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            f = h5py.File(filename,'r')
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
        # temporary solution in order to avoid
        # creating save for Interactions classes

        if self.filled:
            #Lname = self.Lfilename
            Lname = '_'.join(self.filename.split('_')[0:-1]) + '.lay'
            #Lname = self.filename.split('_')[0] + '.lay'
            L=Layout(Lname)
            self.fillinter(L)

        if self.evaluated:
            return self.val(self.fGHz)

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
        #try:

        fh5=h5py.File(filenameh5,'a')
        if self.is3D:
            if not grpname in fh5['ray'].keys():
                fh5['ray'].create_group(grpname)
            else :
                print('ray/'+grpname +'already exists in '+filenameh5)
            f = fh5['ray/'+grpname]


        else:
            if not grpname in fh5['ray2'].keys():
                fh5['ray2'].create_group(grpname)
            else :
                print('ray2/'+grpname +'already exists in '+filenameh5)
            f = fh5['ray2/'+grpname]
        # keys not saved as attribute of h5py file
        notattr = ['I','B','B0','dis']
        for a in self.__dict__.keys():
            if a not in notattr:
                if type(a)==str:
                    a.encode('utf-8')
                    if a=='_luw':
                        la = [ x.encode('utf8') for x in getattr(self,a) ] 
                        f.attrs[a] = la
                    else:
                        f.attrs[a] = getattr(self,a)

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
                    if kk=='diffslabs':
                        ldiffslabs = [ x.encode('utf8') for x in self[k][kk] ]
                        f[str(k)].create_dataset(kk,shape=np.shape(self[k][kk]),data=ldiffslabs)
                    else:
                        f[str(k)].create_dataset(kk,shape=np.shape(self[k][kk]),data=self[k][kk])
        fh5.close()
        #except:
        #    fh5.close()
        #    raise NameError('Rays: issue when writting h5py file')

    def _loadh5(self,filenameh5,grpname,**kwargs):
        """ load rays  h5py format compliant with Links Class

        Parameters
        ----------

        filenameh5 : string
            filename of the h5py file (from Links Class)
        grpname : string
            groupname of the h5py file (from Links Class)
        kwargs may contain a L: layout object
            if L =  [] the layout is loaded from the layout name stored
            into the h5 file
            if L = Layout the layout passed in arg is used

        See Also
        --------

        pylayers.simul.links

        """


        filename=pyu.getlong(filenameh5,pstruc['DIRLNK'])
        # try/except to avoid loosing the h5 file if
        # read/write error
        try:
            fh5=h5py.File(filename,'r')

            if self.is3D:
                argfile = 'ray/'+grpname
            else:
                argfile = 'ray2/'+grpname

            f = fh5[argfile]

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
            if 'L' in kwargs:
                self.L=kwargs['L']
            else:
                self.L = Layout(self.Lfilename,bbuild=True)
                try:
                    self.L.dumpr()
                except:
                    self.L.build()
                    self.L.dumpw()
            # L=Layout(self.Lfilename,bbuild=True)
            self.fillinter(self.L)

        # if self.evaluated:
        #     return self.eval(self.fGHz)


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
        """ check ray reciprocity in comparing two reciprocal rays

        Parameters
        ----------

        r : rays reciprocal to self


        """
        # permutation of all termination points
        assert (self.pTx==r.pRx).all()
        assert (self.pRx==r.pTx).all()
        # for all group of interctions
        for k in self:
            # same distances
            assert (np.allclose(self[k]['dis'],r[k]['dis']))
            # same points when reading from right to left
            assert (np.allclose(self[k]['pt'],r[k]['pt'][:,::-1,:]))
            # same signature reading from right to left
            assert (np.allclose(self[k]['sig'],r[k]['sig'][:,::-1,:]))
            # if local basis have been evaluated
            if (self.isbased) & (r.isbased):
                #assert (np.allclose(self[k]['nstrwall'],r[k]['nstrwall'][:,::-1,:]))
                assert (np.allclose(self[k]['norm'],r[k]['norm'][:,::-1,:])), 'interaction block:' + str(k)
                #assert ((np.mod(self[k]['aoa']-r[k]['aod'],2*np.pi)==0).all())
                #assert ((np.mod(self[k]['aod']-r[k]['aoa'],2*np.pi)==0).all())
                # 1st output basis is equal to last input basis of the reciprocal ray
                assert (np.allclose(self[k]['Bo0'],r[k]['BiN'])), 'interaction block:' + str(k)
                # last input basis is equal to 1st output basis of the reciprocal ray
                assert (np.allclose(self[k]['BiN'],r[k]['Bo0'])), 'interaction block:' + str(k)
                # vsi vectors are inversed
                assert (np.allclose(self[k]['vsi'],-r[k]['vsi'][:,::-1,:])), 'interaction block:' + str(k)
                assert (np.allclose(abs(self[k]['scpr']),abs(r[k]['scpr'][::-1,:]))), 'interaction block:' + str(k)
                assert (np.allclose(self[k]['theta'],r[k]['theta'][::-1,:])), 'interaction block:' + str(k)
                assert (np.allclose(self[k]['Bi'],r[k]['Bo'][:,:,::-1,:])), 'interaction block:' + str(k)
                assert (np.allclose(self[k]['Bo'],r[k]['Bi'][:,:,::-1,:])), 'interaction block:' + str(k)
                assert (np.allclose(self[k]['B'],r[k]['B'][:,:,::-1,:].swapaxes(0,1))), 'interaction block:' + str(k)

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

        TODO : not finished

        """
        u = np.argsort(self.dis)


    def rayfromtyp_order(self,nD=[1],nR=[1],nT=[1],llo='&&'):
        """
            Return rays from a given type (R|T|D) to a given order
            ( number of interaction)

            list logic operator : llo ['op0op1'] 

            nD <op0> nR <op1> nT


            Parameters
            ----------

                nD = list|int
                    requested number of Diffraction
                nR = list|int
                    requested number of Reflection
                nT = list|int
                    requested number of Transmission
                llo = list logic operator [op0,op1]
                    nD <op0> nR <op1> nT


            Returns
            -------

                lr : list
                    list of ray index matching the typ & order conditions



        """

        if not isinstance(nD,list):
            nD=[nD]
        if not isinstance(nR,list):
            nR=[nR]
        if not isinstance(nT,list):
            nT=[nT]

        op = {'and' : operator.and_,
              'or': operator.or_,
              '&': operator.and_,
              '|': operator.or_,
              }


        lr=[]
        for ur,r in enumerate(range(self.nray)):
            li = self.ray2ityp(r)
            nRli = li.count('R')
            nTli = li.count('T')
            nDli = li.count('D')



            cD = (nDli in nD)
            cR = (nRli in nR)
            cT = (nTli in nT)

            # if (nDli in nD) and (nRli in nR) and (nTli in nT) :
            if op[llo[1].lower()]( op[llo[0].lower()](cD,cR) , cT):
                lr.append(r)
            elif (self.los) and (1 in nT ) and (0 in nD) and (0 in nR) and (ur == 0):
                lr.append(r)
        return lr


    def extract_typ_order(self,L,nD=[1],nR=[1],nT=[1],llo='&&'):
        """ Extract group of rays from a certain type (R|T|D) 
            at a order ( <=> given number of interaction)

            list logic operator : llo [op0,op1]

            nD <op0> nR <op1> nT


        Parameters
            ----------

        L  : Layout
        nD = list|int
            requested number of Diffraction
        nR = list|int
            requested number of Reflection
        nT = list|int
            requested number of Transmission
        llo = list logic operator [op0,op1]
            nD <op0> nR <op1> nT

        Returns
        -------

        R : Rays object
            New Rays object containing rays matching
            the typ/order conditions


        """

        lr = self.rayfromtyp_order(nD=nD,nR=nR,nT=nT,llo=llo)
        return self.extract(lr,L)


    def extract(self,lnr,L):
        """ Extract a group of rays

        Parameters
        ----------

        lnr : list of rays indexes
        L  : Layout

        """


        if not isinstance(lnr,list):
            lnr=[lnr]

        r = Rays(self.pTx,self.pRx)
        r.is3D = self.is3D

        for unr,nr in enumerate(lnr):

            ni = self.ray2nbi(nr)
            ur = np.where(self[ni]['rayidx']==nr)[0][0]

            if ni == 0:
                los = True
            else:
                los = False

            if 'D' in self.typ(nr):
                diff=True
            else:
                diff=False


            if 'diffvect' in self[ni]:
                # check if the ray has diffraction interaction
                inter = self.ray2iidx(nr)[:,0]
                uD = np.where([i in inter for i in self[ni]['diffidx']])[0]
            else:
                uD=[]

            diffkey = ['diffvect','diffidx','diffslabs']

            cray = {}


            for k in self[ni].keys():

                if ni ==0:

                    cray = self[ni]
                    break

                elif k not in ['nbrays','rayidx','dis','nstrwall','nstrswall']:
                    tab  = self[ni][k]
                    if type(tab)==np.ndarray and k not in diffkey:
                            try:
                                cray[k] = tab[...,ur][...,np.newaxis]
                            except:
                                import ipdb
                                ipdb.set_trace()
                    if diff : 
                        if k in diffkey :
                            if k != 'diffslabs':
                                cray[k]=tab[...,uD][...,np.newaxis]
                            else:
                                if len(uD)>0 :
                                    cray[k]=[tab[uD]]
                                else:
                                    cray[k]=[]


            cray['nbrays'] = unr+1 # keep only one ray
            r.nray = unr+1
            #cray['rayidx']=np.array([self[ni]['rayidx'][nr]]) # ray index in the whole structure
            cray['rayidx'] = np.array([unr])
            cray['dis'] = np.array([self[ni]['dis'][ur]])


            if ni in r:

                # R[ni]['sig2d'].append(self[k]['sig2d'][ur])

                if not los :
                    r[ni]['BiN'] = np.concatenate((r[ni]['BiN'],cray['BiN']),axis=2)
                    r[ni]['Bo'] = np.concatenate((r[ni]['Bo'],cray['Bo']),axis=3)
                    r[ni]['Bi'] = np.concatenate((r[ni]['Bi'],cray['Bi']),axis=3)


                if diff:
                    if 'diffidx' in r[ni]:
                        r[ni]['diffidx'] = np.concatenate((r[ni]['diffidx'],cray['diffidx']))
                        r[ni]['diffvect'] = np.concatenate((r[ni]['diffvect'],cray['diffvect']),axis=1)
                        r[ni]['diffslabs'].append(cray['diffslabs'])

                    else:
                        r[ni]['diffidx'] = cray['diffidx']
                        r[ni]['diffvect'] = cray['diffvect']
                        r[ni]['diffslabs'] = cray['diffslabs']

                r[ni]['nbrays'] += 1
                r[ni]['B'] = np.concatenate((r[ni]['B'], cray['B']), axis=3)

                r[ni]['pt'] = np.concatenate((r[ni]['pt'], cray['pt']), axis=2)
                r[ni]['rayidx'] = np.concatenate((r[ni]['rayidx'], cray['rayidx']), axis=0)
                r[ni]['Bo0'] = np.concatenate((r[ni]['Bo0'],cray['Bo0']), axis=2)

                r[ni]['scpr'] = np.concatenate((r[ni]['scpr'], cray['scpr']), axis=1)
                r[ni]['aod'] = np.concatenate((r[ni]['aod'], cray['aod']), axis=1)
                r[ni]['si'] = np.concatenate((r[ni]['si'], cray['si']), axis=1)
                r[ni]['sig'] = np.concatenate((r[ni]['sig'], cray['sig']), axis=2)
                # r[ni]['sig2d'] = np.concatenate((r[ni]['sig2d'],cray['sig2d']),axis=2)
                r[ni]['aoa'] = np.concatenate((r[ni]['aoa'], cray['aoa']), axis=1)
                r[ni]['vsi'] = np.concatenate((r[ni]['vsi'], cray['vsi']), axis=2)
                r[ni]['theta'] = np.concatenate((r[ni]['theta'], cray['theta']), axis=1)
                r[ni]['norm'] = np.concatenate((r[ni]['norm'], cray['norm']), axis=2)
                r[ni]['dis'] = np.concatenate((r[ni]['dis'], cray['dis']), axis=0)

            else:
                r[ni] = cray

        # r[ni]['rays'] = to be done HERE


        r.locbas(L)
        r.fillinter(L)
        return(r)

    def extract_old(self,nr,L):
        """ Extract a single ray

        Parameters
        ----------

        nr : ray index
        L  : Layout

        """

        r = Rays(self.pTx,self.pRx)
        r.is3D = self.is3D

        r.nray2D = 1
        r.nb_origin_sig = 1

        #ni = self._ray2nbi[nr]
        #ur = np.where(self[ni]['rayidx']==nr)[0][0]

        ni,ur = self.ir2a(nr)

        if 'D' in self.typ(nr):
            diff=True
        else:
            diff=False

        if 'diffvect' in self[ni]:
            # check if the ray has diffraction interaction
            inter = self.ray2iidx(nr)[:,0]
            uD = np.where([i in inter for i in self[ni]['diffidx']])[0]
        else:
            uD=[]

        diffkey = ['diffvect','diffidx','diffslabs']

        r[ni] = {}
        for k in self[ni].keys():
            if k not in ['nbrays','rayidx','dis','nstrwall','nstrswall']:
                tab  = self[ni][k]
                if type(tab)==np.ndarray and k not in diffkey:
                        r[ni][k] = tab[...,ur][...,np.newaxis]
                if diff : 
                    if k in diffkey :
                        if k != 'diffslabs':
                            r[ni][k]=tab[...,uD][...,np.newaxis]
                        else:
                            if len(uD)>0 :
                                r[ni][k]=tab[uD]
                            else:
                                r[ni][k]=[]

        r[ni]['nrays'] = 1 # keep only one ray
        r.nray = 1
        #r[ni]['rayidx']=np.array([self[ni]['rayidx'][nr]]) # ray index in the whole structure
        r[ni]['rayidx'] = np.array([0])
        r[ni]['dis'] = np.array([self[ni]['dis'][ur]])
        r.locbas(L)
        r.fillinter(L)
        return(r)


    def show(self,**kwargs):
        """  plot 2D rays within the simulated environment

        Parameters
        ----------

        rlist : list  (default []= all rays)
            list of indices of ray in interaction group
        graph : string t
            type of graph to be displayed
            's','r','t',..
        fig : figure
        ax  : axis
        L   : Layout
        alpha : float
            1
        linewidth : float
            0.1
        color : string
            'black'
        ms : int
            marker size :  5
        layout : boolean
            True
        points : boolean
            True
        ER : ray energy

        """
        defaults = {'rlist': [],
                    'fig': [],
                    'ax': [],
                    'L': [],
                    'graph': 's',
                    'color': 'black',
                    'alpha': 1,
                    'linewidth': 0.5,
                    'ms': 5,
                    'vmin':0,
                    'vmax':-70,
                    'cmap': plt.cm.hot_r,
                    'layout': True,
                    'points': True,
                    'labels': False,
                    'bcolorbar': False
                   }

        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        if kwargs['fig'] ==[]:
            fig = plt.figure()

        if kwargs['ax'] ==[]:
            ax = fig.add_subplot(111)

        #
        # display the Layout
        #
        if kwargs['layout'] == True:
            if kwargs['L'] != []:
                fig,ax = kwargs['L'].showG(**kwargs)
            else :
                raise AttributeError('Please give a Layout file as argument')
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
        # plot all rays
        if kwargs['rlist'] == []:

            # list of group of interactions
            lgrint = self.keys()

            for i in lgrint:
                # list of rays
                lray = range(len(self[i]['pt'][0, 0, :]))

                #if self.filled :
                #    ax.set_title('rays index :'+ str(self[i]['rayidx']))

                for j in lray:

                    addr_ray = (i,j)
                    index_ray = self.a2ir(addr_ray)

                    ray = np.hstack((self.pTx[0:2].reshape((2, 1)),
                                     np.hstack((self[i]['pt'][0:2, :, j],
                                     self.pRx[0:2].reshape((2, 1))))
                                     ))

                    if 'ER' not in kwargs:
                        ax.plot(ray[0, :], ray[1, :],
                            alpha = kwargs['alpha'],
                            color = kwargs['color'],
                            linewidth = kwargs['linewidth'])
                    else:
                        EdB = 10*np.log10(ER[index_ray])
                        ERdB = 10*np.log10(E)
                        vscale  = 1.-(max(ERdB)-EdB)/(max(ERdB)-min(ERdB))
                        linewidth = 3*vscale
                        alpha = vscale
                        cmap = cm.hot
                        color = cmap(vscale)
                        ax.plot(ray[0, :], ray[1, :],
                                alpha = alpha,
                                color = color,
                                linewidth = linewidth)

                    ax.axis('off')
                    #if self.filled :
                    #    ax.set_title('rays index :'+ str(self[i]['rayidx'][lray]))
        else:
            rlist = kwargs['rlist']
            # 3D ray
            if self.is3D:
                nbi = self._ray2nbi[rlist]
                nr = np.array((nbi,rlist))
                unb = np.unique(nr[0,:])
                unr = {int(i):np.where(nr[0,:]==i)[0] for i in unb}

                for i in unb:
                    raynb = (nr[1,unr[i]]).astype(int)
                    nbr = len(raynb)
                    ptidx = [np.where(self[i]['rayidx']==x)[0][0] for x in raynb]
                    for j in ptidx:

                        ray = np.hstack((self.pTx[0:2].reshape((2, 1)),
                                         np.hstack((self[i]['pt'][0:2, :, j],
                                         self.pRx[0:2].reshape((2, 1))))
                                         ))
                        ax.plot(ray[0, :], ray[1, :],
                                alpha = kwargs['alpha'],
                                color = kwargs['color'],
                                linewidth = kwargs['linewidth'])
                        ax.axis('off')
            # 2D ray
            else:
                for i in rlist:
                    lray = range(len(self[i]['pt'][0, 0, :]))
                    #if self.filled :
                    #    ax.set_title('rays index :'+ str(self[i]['rayidx']))
                    for j in lray:
                        ray = np.hstack((self.pTx[0:2].reshape((2, 1)),
                                         np.hstack((self[i]['pt'][0:2, :, j],
                                         self.pRx[0:2].reshape((2, 1))))
                                         ))
                        ax.plot(ray[0, :], ray[1, :],
                                alpha=kwargs['alpha'],
                                color=kwargs['color'],
                                linewidth=kwargs['linewidth'])
                        ax.axis('off')

        #if kwargs['bcolorbar']:
        #    # axes : left , bottom , width , height
        #    sm = plt.cm.ScalarMappable(cmap = kwargs['cmap'], norm = plt.Normalize(vmin=kwargs['vmin'],vmax=kwargs['vmax']))
        #    sm._A = []  # necessary set_array
        #    cax = fig.add_axes([0.18,0.35, 0.35, 0.025])
        #    #cb = plt.colorbar(sm,cax=cax,orientation='horizontal')
        #    cb = plt.colorbar(sm,cax=cax,orientation='horizontal')
        #    cb.ax.tick_params(labelsize=24)
        #    cb.set_label('Level (dB)', fontsize=24)

        return(fig,ax)

    def mirror(self, H=3, N=1, za = [], zb= []):
        """ mirror a ray termination

        Parameters
        ----------

        H : float
            ceil height (default 3m)
            if H=0 only floor reflection is calculated (outdoor case)
            if H=-1 floor and ceil reflection are inhibited (2D test case)
        N : int
            handle the number of mirror reflexions

        za : float
            height of the point where the parametrization starts ( e.g. pTx[2])

        zb : float
            height of the point where the parametrization ends ( e.g. pRx[2])


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

        Notes
        -----

        d is a dictionnary whose keys are heights along the vertical from where
        are emanating the reflected rays. Values of d are the parameterization
        (0< () <1) along the corresponding to the different reflection
        points.

        See Also
        --------

        to3D

        """



        km = np.arange(-N+1, N+1, 1)
        kp = np.arange(-N, N+1, 1)
        #
        # heights of transmitter and receiver
        #
        if type(za)==list:
            za=self.pTx[2]
        if type(zb)==list:
            zb=self.pRx[2]
        ht = za
        hr = zb

        assert (hr<H or H==0 or H == -1),"mirror : receiver higher than ceil height"
        assert (ht<H or H==0 or H == -1),"mirror : transmitter higher than ceil height"

        zkp = 2*kp*H + ht
        zkm = 2*km*H - ht

        d = {}
        if H>0:
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
        elif H==0:
            d[-ht] = np.array([ht/(ht+hr)])
            d[ht] = np.array([])
        elif H==-1:
            d[ht] = np.array([])

        return d

    def mirror_new(self, N=1, zlevels=[0, 3] , za = None, zb= None):
        """ mirror a ray termination

        Parameters
        ----------

        N : int
            number of mirror reflexions

        zlevels : np.array  (,Nlevels)
            list of levels height

        za : float
            height of the point where the parametrization starts ( e.g. pTx[2])

        zb : float
            height of the point where the parametrization ends ( e.g. pRx[2])


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
        >>> d = r.mirror_new(zlevels=[0,])
        >>> d[-1.5]
        array([ 0.55555556])

        Notes
        -----

        d is a dictionnary whose keys are heights along the vertical from where
        are emanating the reflected rays. Values of d are the parameterization
        (0< () <1) along the corresponding to the different reflection
        points.

        See Also
        --------

        to3D

        """



        km = np.arange(-N+1, N+1, 1)
        kp = np.arange(-N, N+1, 1)
        #
        # heights at both ends
        #
        if za == None:
            za=self.pTx[2]
        if zb == None:
            zb=self.pRx[2]


        #for 
        return d

    def to3D(self, L, H=3, N=1, rmoutceilR=True):
        """ transform 2D ray to 3D ray

        Parameters
        ----------

        L : Layout object

        H : float
            ceil height (default 3m)
            if H = 0 only floor reflection is calculated (outdoor case)
            if H = -1 floor and ceil reflection are deactivated (2D test case)
        N : int
            number of mirror reflexions
        rmoutceilR : bool
            Remove ceil reflexions in cycles (Gt nodes)
            with indoor == False attribute

        Returns
        -------

        r3d : Rays

        See Also
        --------

        mirror

        """

        if H==-1:
            rmoutceilR = False

        tx = self.pTx
        rx = self.pRx

        #
        # Phase 1 : calculate parameterizario of the Tx images height
        #           in the vertical plane
        #

        d = self.mirror(H=H, N=N, za=tx[2], zb=rx[2])

        #
        # Elimination of invalid diffraction point
        # If the diffaction point is a separation between 2 air wall
        # it should be removed.
        #


        #
        # Phase 2 :
        #    calculate 2D parameterization in the horizontal plane
        #

        # for all group of interactions

        for i in self:

            pts = self[i]['pt'][0:2, :, :]
            sig = self[i]['sig']

            if pts.shape[2]!=0:
                # broadcasting of t and r
                t = self.pTx[0:2].reshape((2, 1, 1)) * \
                np.ones((1, 1, len(pts[0, 0, :])))
                r = self.pRx[0:2].reshape((2, 1, 1)) * \
                np.ones((1, 1, len(pts[0, 0, :])))
                pts1 = np.hstack((t, np.hstack((pts, r))))
            else:
                t = self.pTx[0:2].reshape((2, 1, 1))
                r = self.pRx[0:2].reshape((2, 1, 1))
                pts1 = np.hstack((t,r))

            # append t and r to interaction points in 2D
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
        r3d.nray2D = len(self)
        r3d.nb_origin_sig = self.nb_origin_sig

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

            if sig.shape[0]!=0:
                sig = np.hstack((np.zeros((2, 1, Nrayk), dtype=int),
                             sig,
                             np.zeros((2, 1, Nrayk), dtype=int)))  # add signature of Tx and Rx (0,0))
            else:
                sig = np.hstack((np.zeros((2, 1, Nrayk), dtype=int),
                                 np.zeros((2, 1, Nrayk), dtype=int)))
            # broadcast tx and rx
            Tx = tx.reshape(3, 1, 1)*np.ones((1, 1, Nrayk))
            Rx = rx.reshape(3, 1, 1)*np.ones((1, 1, Nrayk))

            if k!=0:
                # pte is the sequence of point in 3D ndim =3   ( ndim x k x Nrayk)
                pte = self[k]['pt']
                # ndim x k+2 x Nrayk
                pte = np.hstack((Tx, pte, Rx))
            else:
                 pte = np.hstack((Tx, Rx))

            # extension
            for l in d:                     # for each vertical pattern (C,F,CF,FC,....)
                #print k,l,d[l]
                Nint = len(d[l])            # number of additional interaction
                #if ((k==1) & (l==5.0)):print
                if Nint > 0:                # if new interaction ==> need extension
                    # a1e : extended horizontal+vertical parameterization
                    a1e = np.concatenate((a1, d[l].reshape(len(d[l]), 1)*
                                          np.ones((1, Nrayk))))
                    # get sorted indices
                    ks = np.argsort(a1e, axis=0)
                    # a1es : extended sorted horizontal + vertical parameterization
                    a1es = np.sort(a1e, axis=0)

                    # #### Check if it exists the same parameter value in the horizontal plane
                    # #### and the vertical plane. Move parameter if so.

                    da1es = np.diff(a1es,axis=0)
                    pda1es = np.where(da1es<1e-10)
                    a1es[pda1es] = a1es[pda1es]-1e-3


                    # prepare an extended sequence of points ( ndim x  (Nint+k+2) x Nrayk )
                    ptee = np.hstack((pte, np.zeros((3, Nint, Nrayk))))

                    #
                    # Boolean ceil/floor detector
                    #
                    # u is 4 (floor interaction )
                    #      5 (ceil interaction )
                    #  depending on the vertical pattern l.
                    #
                    #  l < 0 corresponds to last reflexion on floor
                    #  l > 0 corresponds to last reflexion on ceil
                    #
                    # u == 0 (floor) or 1 (ceil)
                    # if l < 0:
                    #     u = np.mod(range(Nint), 2)
                    # else:
                    #     u = 1 - np.mod(range(Nint), 2)


                    if l < 0 and Nint%2 ==1: # l<0 Nint odd
                        u = np.mod(range(Nint), 2)

                    elif l > 0 and Nint%2 ==1: # l>0 Nint odd
                        u = 1 - np.mod(range(Nint), 2)

                    elif l < 0 and Nint%2 ==0: # l<0 Nint even
                        u = 1 - np.mod(range(Nint), 2)

                    elif l > 0 and Nint%2 ==0: # l>0 Nint even
                        u = np.mod(range(Nint), 2)

                    #
                    u = u + 4
                    #
                    # At that point we introduce the signature of the new
                    # introduced points on the ceil and/or floor.
                    #
                    # A signature is composed of two lines
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
                    # Tous les points prcdents qui ne sont pas des Ceils ou
                    # des floors et tous les points suivants qui ne sont pas
                    # des points de rflexion ceil ou floor
                    #
                    # Afin de tenir compte du rayon et du groupe d'interactions
                    # concerne, il faut passer un tuple qui concatene la valeur
                    # de l'indice d'interaction floor ou ceil et l'indice de
                    # rayons du groupe associe (d'ou le zip)
                    #
                    # Cette sequence d'instruction fixe le bug #133
                    #
                    # Anterieurement il y avait une hypothese de succession
                    # immediate d'un point 2D renseigne.
                    #
                    try:
                        iintm_f = [ np.where( (siges[1,0:x[0],x[1]]!=4) &
                                             (siges[1,0:x[0],x[1]]!=5))[0][-1]
                                   for x in  zip(iint_f,iray_f) ]
                        iintp_f = [ np.where( (siges[1,x[0]:,x[1]]!=4) &
                                             (siges[1,x[0]:,x[1]]!=5))[0][0]+x[0]
                                   for x in  zip(iint_f,iray_f) ]
                        iintm_c = [ np.where( (siges[1,0:x[0],x[1]]!=4) &
                                             (siges[1,0:x[0],x[1]]!=5))[0][-1]
                                   for x in zip(iint_c,iray_c) ]
                        iintp_c = [ np.where( (siges[1,x[0]:,x[1]]!=4) &
                                             (siges[1,x[0]:,x[1]]!=5))[0][0]+x[0]
                                   for x in  zip(iint_c,iray_c) ]
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


                    #
                    # If there are floor points
                    #

                    if len(iint_f)>0:
                        a1esm_f = a1es[iintm_f, iray_f]
                        a1esc_f = a1es[iint_f, iray_f]
                        a1esp_f = a1es[iintp_f, iray_f]


                        pteesm_f = ptees[0:2, iintm_f, iray_f]
                        pteesp_f = ptees[0:2, iintp_f, iray_f]

                        coeff_f = (a1esc_f-a1esm_f)/(a1esp_f-a1esm_f)

                        ptees[0:2, iint_f, iray_f] = pteesm_f + coeff_f*(pteesp_f-pteesm_f)

                    #
                    # If there are ceil points
                    #
                    if len(iint_c)>0:
                        a1esm_c = a1es[iintm_c, iray_c]
                        a1esc_c = a1es[iint_c, iray_c]
                        a1esp_c = a1es[iintp_c, iray_c]

                        pteesm_c = ptees[0:2, iintm_c, iray_c]
                        pteesp_c = ptees[0:2, iintp_c, iray_c]

                        coeff_c = (a1esc_c-a1esm_c)/(a1esp_c-a1esm_c)
                        ptees[0:2, iint_c, iray_c] = pteesm_c + coeff_c*(pteesp_c-pteesm_c)

                    if H != 0:
                        z  = np.mod(l+a1es*(rx[2]-l), 2*H)
                        pz = np.where(z > H)
                        z[pz] = 2*H-z[pz]
                        ptees[2, :] = z
                    # case where ceil reflection are inhibited
                    elif H==0:
                        z  = abs(l+a1es*(rx[2]-l))
                        # pz = np.where(z > H)
                        # z[pz] = 2*H-z[pz]
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
                # handling multi segment (iso segments)
                #    Height of reflexion interaction
                #    Height of diffraction interaction
                #---------------------------------
                #
                #   ptes (3 x i+2 x r )
                if len(L.lsss)>0:
                    #
                    # lsss : list of sub segments ( iso segments siges)
                    # lnss : list of diffaction point involving

                    lsss = np.array(L.lsss)
                    lnss = np.array(L.lnss)

                    # array of structure element (nstr) with TxRx extension  (nstr=0)
                    anstr = siges[0,:,:]
                    # type of interaction
                    typi = siges[1,:,:]

                    # lss : list of subsegments in the current signature
                    #
                    # scalability : avoid a loop over all the subsegments in lsss
                    #
                    lss = [ x for x in lsss if x in anstr.ravel()]

                    ray_to_delete = []
                    for s in lss:
                        u  = np.where(anstr==s)
                        if len(u)>0:
                            zs = ptees[2,u[0],u[1]]
                            zinterval = L.Gs.node[s]['z']
                            unot_in_interval = ~((zs<=zinterval[1]) & (zs>=zinterval[0]))
                            ray_to_delete.extend(u[1][unot_in_interval])

                    # lns : list of diffraction points in the current signature
                    #       involving multi segments (iso)
                    # scalability : avoid a loop over all the points in lnss
                    #
                    lns = [ x for x in lnss if x in anstr.ravel()]

                    #
                    # loop over multi diffraction points
                    #

                    for npt in lns:
                        u  = np.where(anstr==npt)
                        if len(u)>0:
                           # height of the diffraction point
                            zp = ptees[2,u[0],u[1]]

                            #
                            # At which couple of segments belongs this height ?
                            # get_diffslab function answers that question
                            #

                            ltu_seg,ltu_slab = L.get_diffslab(npt,zp)

                            #
                            # delete rays where diffraction point is connected to
                            # 2 AIR segments
                            #
                            [ray_to_delete.append(u[1][i]) for i in range(len(zp))
                            if ((ltu_slab[i][0]=='AIR') & (ltu_slab[i][1]=='AIR'))]
                            # #zinterval = L.Gs.node[s]['z']
                            # # if (zs<=zinterval[1]) & (zs>=zinterval[0]):
                            # if ((tu_slab[0]!='AIR') & (tu_slab[1]!='AIR')):
                            #     #print(npt , zp)
                            #     pass
                            # else:
                            #     ray_to_delete.append(u[1][0])

                    # # nstr : structure number
                    # nstr  = np.delete(nstr,ray_to_delete,axis=1)
                    # typi : type of interaction
                    typi  = np.delete(typi,ray_to_delete,axis=1)
                    # 3d sequence of points
                    ptees = np.delete(ptees,ray_to_delete,axis=2)
                    # extended (floor/ceil) signature
                    siges = np.delete(siges,ray_to_delete,axis=2)

                if rmoutceilR:
                    # 1 determine Ceil reflexion index
                    # uc (inter x ray)
                    uc = np.where(siges[1,:,:]==5)
                    ptc = ptees[:,uc[0],uc[1]]
                    if len(uc[0]) !=0:
                        P = shg.MultiPoint(ptc[:2,:].T)
                        # to determine the cycle where ceil reflexions append
                        # uinter(nb pt x nb cycles)
                        mapnode = list(L.Gt.nodes())
                        uinter = np.array([[L.Gt.node[x]['polyg'].contains(p) for x in mapnode if x>0] for p in P])
                        # import ipdb
                        # ipdb.set_trace()
                        #[plt.scatter(p.xy[0],p.xy[1],c='r') for up,p in enumerate(P) if uinter[0,up]]
                        #[ plt.scatter(p.xy[0],p.xy[1],c='r') for up,p in enumerate(P) if uinter[0,up]]
                        # find points are indoor/outdoor cycles
                        upt,ucy = np.where(uinter)
                        uout = np.where([not L.Gt.node[mapnode[u+1]]['indoor'] for u in ucy])[0] #ucy+1 is to manage cycle 0
                        # 3 remove ceil reflexions of outdoor cycles
                        if len(uout)>0:
                            ptees = np.delete(ptees,uc[1][uout],axis=2)
                            siges = np.delete(siges,uc[1][uout],axis=2)
                            sigsave = np.delete(sigsave,uc[1][uout],axis=2)

                if k+Nint in r3d:
                    r3d[k+Nint]['pt']  = np.dstack((r3d[k+Nint]['pt'], ptees))
                    r3d[k+Nint]['sig'] = np.dstack((r3d[k+Nint]['sig'], siges))
                    r3d[k+Nint]['sig2d'].append(sigsave)
                else:
                    if ptees.shape[2]!=0:
                        r3d[k+Nint] = {}
                        r3d[k+Nint]['pt'] = ptees
                        r3d[k+Nint]['sig'] = siges
                        r3d[k+Nint]['sig2d'] = [sigsave]
                # ax=plt.gca()
                # uu = np.where(ptees[2,...]==3.0)
                # ax.plot(ptees[0,uu[0],uu[1]],ptees[1,uu[0],uu[1]],'ok')
                # import ipdb
                # ipdb.set_trace()
        #
        # Add Line Of Sight ray information
        #   pt =  [tx,rx]
        #   sig = [0,0]
        #
        #pdb.set_trace()
        # if (self.los) & (np.sqrt(np.sum((tx-rx)**2)) !=0) :
        #     r3d[0] = {}
        #     r3d[0]['sig'] = np.zeros((2,2,1))
        #     r3d[0]['sig2d'] = np.zeros((2,2,1))
        #     r3d[0]['pt'] = np.zeros((3,2,1))
        #     r3d[0]['pt'][:,0,:] = tx[:,np.newaxis]
        #     r3d[0]['pt'][:,1,:] = rx[:,np.newaxis]

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

        r3d.delays = np.zeros((r3d.nray))
        for k in r3d.keys():
            ir = r3d[k]['rayidx']
            r3d.delays[ir] = r3d[k]['dis']/0.3


        r3d.origin_sig_name = self.origin_sig_name
        r3d.Lfilename = L._filename
        r3d.filename = L._filename.split('.')[0] + '_' + str(r3d.nray)
        return(r3d)


    def get_rays_slabs(self,L,ir):
        """ return the slabs for a given interaction index 


            Parameters
            ----------

            L : Layout
            ir : interaction block

            Returns
            -------

            numpy array of slabs strings at the shape (ir,r)
            ir : number of interactions ( of the interaction block)
            r : number of rays

        """

        v=np.vectorize( lambda t: L.Gs.node[t]['name'] if (t!=0) and (t>0) else '_')
        return v(self[ir]['sig'][0])


    def remove_aw(self,L):
        """ remove AIR interactions
        """
        # def consecutive(data, stepsize=1):
        #     return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)


        R = Rays(self.pTx,self.pRx)
        R.__dict__.update(self.__dict__)
        # R.is3D=True
        # R.nray = self.nray
        # R.nray2D = self.nray2D
        # R.nray2D = self.nray2D
        # R.nray2D = self.nray2D

        for k in self:
            lr = self[k]['sig'].shape[1]

            inter = self.get_rays_slabs(L,k)

            for ur,r in enumerate(inter.T):

                not_air_mask  = ~((r =='_AIR') | (r == 'AIR' ))
                nb_air = sum(~not_air_mask)
                if nb_air != 0 :
                    new_bi = k-nb_air
                    # +2 : add tx & rx interaciton
                    # -1 : 2 interactions correspond to 1 distance
                    lsi = new_bi + 2 - 1 
                    si = np.zeros(lsi) 
                    si_old = self[k]['si'][:,ur]

                    vsi = np.zeros((3,lsi)) 
                    vsi_old = self[k]['vsi'][...,ur]

                    sig = self[k]['sig'][:,not_air_mask,ur][...,None]
                    # sig2d = self[k]['sig2d'][0][...,ur]
                    pt = self[k]['pt'][:,not_air_mask,ur][...,None]

                    u = 0
                    si_aw = 0

                    # import ipdb
                    # ipdb.set_trace()

                    for uold,b in enumerate(not_air_mask[1:]):
                        if b:
                            # update new si with sum of all 
                            # distance from preceding airwall
                            si[u] = si_old[uold] + si_aw
                            # keep vsi from the last airwall
                            # because vsi don't change on an airwall
                            vsi[:,u] = vsi_old[:,uold] 
                            u += 1
                            si_aw=0
                        else:
                            si_aw += si_old[uold]
                    si = si[...,None]
                    vsi = vsi[...,None]
                    dis = np.array([np.sum(si)])
                    assert np.allclose(dis,np.sum(si_old))



                else:
                    # no air wall case, fill R with self values
                    new_bi = k
                    pt = self[k]['pt'][...,ur][...,None]
                    sig = self[k]['sig'][...,ur][...,None]
                    # sig2d = self[k]['sig2d'][0][...,ur]
                    si = self[k]['si'][:,ur][:,None]
                    vsi = self[k]['vsi'][...,ur][...,None]
                    dis = np.array([self[k]['dis'][ur]])

                if new_bi in R:

                    # R[new_bi]['sig2d'].append(self[k]['sig2d'][ur])
                    R[new_bi]['pt'] = np.concatenate((R[new_bi]['pt'],pt),axis=2)
                    R[new_bi]['sig'] = np.concatenate((R[new_bi]['sig'],sig),axis=2)
                    R[new_bi]['rayidx'] = np.concatenate((R[new_bi]['rayidx'],np.array([self[k]['rayidx'][ur]])))
                    R[new_bi]['si'] = np.concatenate((R[new_bi]['si'],si),axis=1)
                    R[new_bi]['vsi'] = np.concatenate((R[new_bi]['vsi'],vsi),axis=2)
                    R[new_bi]['dis'] = np.concatenate((R[new_bi]['dis'],dis),axis=0)
                else:
                    R[new_bi] = {}
                    # R[new_bi]['sig2d'] = [self[k]['sig2d'][0][...,ur]]
                    R[new_bi]['pt'] = pt
                    R[new_bi]['sig'] = sig
                    R[new_bi]['rayidx'] = np.array([self[k]['rayidx'][ur]])
                    R[new_bi]['si'] = si
                    R[new_bi]['vsi'] = vsi
                    R[new_bi]['dis'] = dis

        if 0 in R:
            R.los=True

        X = [[R[k]['rayidx'][u] for u in range(len(R[k]['rayidx']))] for k in R]
        R._rayidx_aw = sum(X,[])

        return R

    def length(self,typ=2):
        """ calculate length of rays

        Parameters
        ----------

        typ : int
            men1 : length of all segments
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

    def simplify(self):
        if not self.is3D:
            return None

        for ir in self:
            print(self[ik]['si'])

    def locbas(self, L):
        """ calculate ray local bas

        Parameters
        ----------

        L : Layout

        Notes
        -----

        This method adds for each group of interactions the following members

        norm : np.array
            3 x i x r  (interaction vector)
        nstrwall : np.array
            nstr of interactions
        vsi : np.array
            3 x (i+1) x r
        aod : np.array
            2 x r
        aoa : np.array
            2 x r
        BoO : np.array
            3 x 3 x r
        Bi  : np.array
            3 x 3 x r
        Bo  : np.array
            3 x 3 x r
        BiN : np.array
            3 x 3 x r
        scpr : np.array
            i x r
        theta : np.array
            i x r
        rays  : int
        nbrays  : int
        rayidx : np.array
        diffslabs : list
        diffvect :  np.array
            (phi0,phi,beta,NN)


        """

        #
        # extract normal in np.array
        #

        # nsegment x 3
        norm = np.array(list(nx.get_node_attributes(L.Gs,'norm').values()))

        # nsegment x k
        key = np.array(list(dict(nx.get_node_attributes(L.Gs,'norm')).keys()))

        # maximum number for refering to segment
        # not to be confused with a segment number

        nsmax = max(L.Gs.node.keys())

        mapping = np.zeros(nsmax+1, dtype=int)

        mapping[key] = np.arange(len(key), dtype=int)

        #
        # Structure number : nstr
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

        # list of used wedges
        luw=[]

        lgi = list(self.keys())
        lgi.sort()
        for k in lgi:
            #
            # k is the number of interactions in the block
            #
            #print(k,self[11]['rayidx'])
            if k != 0:

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
                # subsegments are not nodes of Gs but have positive nstr index
                #

                uss   = np.where(nstr > nsmax)

                # print uss

                nstrs = copy.copy(nstr)
                #
                # if subsegments have been found
                #
                if len(uss) >0:
                    ind   = nstr[uss]- nsmax-1
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

                # norm : 3 x i x r
                #
                # norm is the vector associated to the interaction
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

                # Bo0 : 3 x 3 x r
                Bo0 = np.concatenate((si[:, 0, None, :],
                                      eth[:, None, :],
                                      eph[:, None, :]), axis=1)

                self[k]['Bo0'] = Bo0

                #
                # scalar product si . norm
                #
                # vn   : 3 x i x r
                # s_in : 3 x i x r

                #
                # scpr : i x r
                #

                scpr = np.sum(vn*si[:,0:-1,:], axis=0)
                self[k]['scpr'] = scpr
                self[k]['theta'] = np.arccos(abs(scpr))  # *180/np.pi

                def fix_colinear(w):
                    """
                    w : vector
                    """
                    nw = np.sqrt(np.sum(w*w, axis=0))
                    u = np.where(nw==0)
                    if len(u[0])!=0:
                        logger.debug('colinear situation detected')
                        if (u[0].any() or u[1].any()) \
                            or (u[0].any()==0 or u[1].any()==0):

                            uu = np.array([u[0],u[1]]).T
                            #determine which interaction and rays
                            #present the colinearity issue
                            uvv = abs(vn[2,uu[:,0],uu[:,1]])>0.99
                            # uv : nbi x nbr colinear index
                            uv = uu[uvv]
                            # uh : nbi x nbr anti-colinear index
                            uh = uu[np.logical_not(uvv)]
                            try:
                                #fix w for colinear index
                                w[:,uv[:,0],uv[:,1]] = np.array(([1,0,0]))[:,None]
                                # update normal
                                nw[uv[:,0],uv[:,1]] = np.sqrt(np.sum(
                                    w[:,uv[:,0],uh[:,1]]*w[:,uv[:,0],uv[:,1]],axis=0))
                            except:
                                pass
                            try:
                                # fix w for anti-colinear index
                                w[:,uh[:,0],uh[:,1]] = np.array(([0,0,1]))[:,None]
                                # update normal
                                nw[uh[:,0],uh[:,1]] = \
                                    np.sqrt(np.sum(w[:,uh[:,0],uh[:,1]]*w[:,uh[:,0],uh[:,1]],axis=0))
                            except:
                                pass
                    return w, nw
                #
                # Warning need to handle singular case when s_in // vn
                #
                # w : 3 x i x r
                #
                w = np.cross(s_in, vn, axisa=0, axisb=0, axisc=0)

                # nw : i x r
                w, nw = fix_colinear(w)

                wn = w/nw
                v = np.cross(wn, s_in, axisa=0, axisb=0, axisc=0)

                es_in = np.expand_dims(s_in, axis=1)

                ew = np.expand_dims(wn, axis=1)
                ev = np.expand_dims(v, axis=1)

                #  Bi 3 x 3 x i x r
                Bi = np.concatenate((es_in,ew,ev),axis=1)
                #  self[k]['Bi'] 3 x 3 x i x r
                self[k]['Bi'] = Bi
                ################################

                w = np.cross(s_out, vn, axisa=0, axisb=0, axisc=0)

                w, nw = fix_colinear(w)
                #wn = w/np.sqrt(np.sum(w*w, axis=0))
                wn = w/nw

                v = np.cross(wn, s_out, axisa=0, axisb=0, axisc=0)

                es_out = np.expand_dims(s_out, axis=1)
                ew = np.expand_dims(wn, axis=1)
                ev = np.expand_dims(v, axis=1)

                #  Bi 3 x 3 x i x r
                Bo = np.concatenate((es_out,ew,ev),axis=1)

                 # self[k]['Bo'] 3 x 3 x i x r
                self[k]['Bo'] = Bo
                #
                # AOA (rad)
                #

                # th : ,r
                # fix doa/dod reciprocity
                #th = np.arccos(si[2, -1, :])
                tha = np.arccos(si[2, -1, :])

                # th : ,r
                #ph = np.arctan2(si[1, -1, :], si[0, -1, :])
                pha = np.arctan2(si[1, -1, :], si[0, -1, :])

                # aoa : 2 x r  (radians)
                self[k]['aoa'] = np.vstack((tha, pha))
                eth = np.array([np.cos(tha) * np.cos(pha),
                               np.cos(tha) * np.sin(pha),
                                -np.sin(tha)])
                eph = np.array([-np.sin(pha),
                                np.cos(pha),
                                np.zeros(len(pha))])
                # Bo0 : 3 x 3 x r
                BiN = np.concatenate((si[:,-1,None,:],
                                      eth[:, None, :],
                                      eph[:, None, :]), axis=1)


                self[k]['BiN'] = BiN
                #self[k]['BiN'] = np.concatenate((-si[:,-1,np.newaxis,:],eth[:,np.newaxis,:],
                #                                   eph[:,np.newaxis,:]),axis=1)

                # Creation of B from Bi and Bo
                # is done after the potential diffraction
                # computation

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
                #pdb.set_trace()
                _ray2nbi = np.ones((nbray), dtype=int)


                try:
                    self._ray2nbi = np.hstack((self._ray2nbi,_ray2nbi))
                except:
                    self._ray2nbi = _ray2nbi


                self._ray2nbi[self[k]['rayidx']]  = k
                nbrayt = nbrayt + nbray
                self.raypt = self.raypt + self[k]['nbrays']

                #################################
                # Start diffraction specific case
                #################################

                if len(udiff[0]) != 0 :
                    Z = np.where(ityp.T==1)
                    udiff=Z[1],Z[0]

                    # diffseg,udiffseg  = np.unique(nstr[udiff],return_inverse=True)
                    diffupt=nstr[udiff]
                    # position of diff seg (- because iupnt accept > 0 reference to points)
                    #
                    # TO BE FIXED
                    #
                    #ptdiff = L.pt[:,L.iupnt[-diffupt]]
                    ptdiff = np.array([ (L.Gs.pos[x][0],L.Gs.pos[x][1])  for x in diffupt ]).T

                    self[k]['diffidx'] = idx[udiff[0],udiff[1]]
                    # get tail head position of seg associated to diff point
                    lair = L.name['AIR'] + L.name['_AIR']
                    #aseg = map(lambda x : filter(lambda y : y not in lair,
                    #                     nx.neighbors(L.Gs,x)),
                    #                     diffupt)

                    aseg = [ [ y for y in nx.neighbors(L.Gs,x) if y not in lair ] for x in diffupt ]

                    #manage flat angle : diffraction by flat segment e.g. door limitation)
                    [aseg[ix].extend(x) for ix,x in enumerate(aseg) if len(x)==1]
                    # get points positions
                    #pdb.set_trace()
                    pts = np.array([ L.seg2pts([x[0],x[1]]) for x in aseg ])

                    #self[k]['diffslabs']=[str(L.sl[L.Gs.node[x[0]]['name']])+'_'
                    #                    + str(L.sl[L.Gs.node[x[1]]['name']]]) for x in aseg]
                    self[k]['diffslabs']=[ L.Gs.node[x[0]]['name']+'@'
                                        +  L.Gs.node[x[1]]['name'] for x in aseg]

                    uwl = np.unique(self[k]['diffslabs']).tolist()
                    luw.extend(uwl)


                    pt1 = pts[:,0:2,0] #tail seg1
                    ph1 = pts[:,2:4,0] #head seg1
                    pt2 = pts[:,0:2,1] #tail seg2
                    ph2 = pts[:,2:4,1] #head seg2


                    #pts is (nb_diffraction_points x 4 x 2)
                    #- The dimension 4 represent the 2x2 points: t1,h1 and t2,h2
                    # tail and head of segment 1 and 2 respectively
                    # a segment
                    #- The dimension 2 is x,y
                    #
                    # The following aims to determine which tails and heads of
                    # segments associated to a given diffraction point
                    # are connected
                    #
                    #

                    # point diff is pt1
                    updpt1 = np.where(np.sum(ptdiff.T==pt1,axis=1)==2)[0]
                    # point diff is ph1
                    updph1 = np.where(np.sum(ptdiff.T==ph1,axis=1)==2)[0]

                    # point diff is pt2
                    updpt2 = np.where(np.sum(ptdiff.T==pt2,axis=1)==2)[0]

                    # point diff is ph2
                    updph2 = np.where(np.sum(ptdiff.T==ph2,axis=1)==2)[0]

                    pa = np.empty((len(diffupt),2))
                    pb = np.empty((len(diffupt),2))

                    ####seg 1 :
                    #if pt1 diff point =>  ph1 is the other point
                    pa[updpt1]= ph1[updpt1]
                    #if ph1 diff point =>  pt1 is the other point
                    pa[updph1]= pt1[updph1]
                    ####seg 2 :
                    #if pt2 diff point =>  ph2 is the other point
                    pb[updpt2]= ph2[updpt2]
                    #if ph2 diff point =>  pt2 is the other point
                    pb[updph2]= pt2[updph2]

                    pt = ptdiff.T

                    # NN : (nb_diffraction_points)
                    # alpha wegde (a.k.a. wedge parameters, a.k.a wedge aperture)

                    NN = (360.-geu.sector(pa.T,pb.T,pt.T))/180.
                    # NN = (2.-NN)*np.pi

                    #angle between face 0, diffraction point and s_in
                    #s_in[:2,udiff[0],udiff[1]]  :
                    # s_in of insteractions udiff (2D) restricted to diffraction points
                    vptpa = pt-pa
                    vptpan = vptpa.T / np.sqrt(np.sum((vptpa)*(vptpa),axis=1))
                    # vpapt= pa-pt # papt : direction vector of face 0
                    # vpaptn = vpapt.T / np.sqrt(np.sum((vpapt)*(vpapt),axis=1))
                    sid = s_in[:,udiff[0],udiff[1]] #s_in restricted to diff
                    sod = s_out[:,udiff[0],udiff[1]] #s_out restricted to diff
                    vnormz = self[k]['norm'][:, udiff[0], udiff[1]]


                    #phi0 = arccos(dot(sid*vpavptn))
                    # phi0 = geu.vecang(sid[:2],vpaptn)
                    uleft = geu.isleft(pa.T,pt.T,pb.T)
                    phi0 = geu.vecang(vptpan,sid[:2])
                    phi0[~uleft] = geu.vecang(sid[:2,~uleft],vptpan[:,~uleft])
                    # phi0 = np.arccos(np.sum(sid[:2]*vpaptn,axis=0))

                    #phi = arccos(dot(sod*vpavptn))
                    # phi = np.arccos(np.sum(-sod[:2]*vpaptn,axis=0))
                    phi = geu.vecang(vptpan,-sod[:2])
                    phi[~uleft] = geu.vecang(-sod[:2,~uleft],vptpan[:,~uleft])
                    # beta
                    #it is important to check if the sid comes from left or right
                    #to this end assume that sid vector is composed
                    #of 2 point : (0,0) and sid
                    # compared to the position of the diffraction point in x
                    # with an elevation=0
                    sidxz = sid[[0,2]]
                    vnormxz = vnormz[[0,2]]
                    zero = np.zeros((2,ptdiff.shape[1]))
                    zdiff = np.vstack((ptdiff[0],zero[0]))
                    left = geu.isleft(zero,sidxz,zdiff)
                    beta = np.arccos(np.sum(vnormz*sid,axis=0))

                    # self[k]['diffvect'] is (4 x Nb_rays )
                    # for axis 0 lenght 4 represent :
                    # 0 => phi0
                    # 1 => phi
                    # 2 => beta
                    # 3 => N (wedge parameter)
                    self[k]['diffvect']=np.array((phi0,phi,beta,NN))

                    ######
                    #Bi diffract
                    #####
                    #w is the \perp \soft in diff
                    w = np.cross(-sid, vnormz, axisa=0, axisb=0, axisc=0)

                    # nw : i x r
                    w, nw = fix_colinear(w)

                    wn = w/nw
                    # Handling channel reciprocity s_in --> -s_in
                    #v = np.cross(wn, s_in, axisa=0, axisb=0, axisc=0)
                    v = np.cross(wn, -sid, axisa=0, axisb=0, axisc=0)

                    e_sid = np.expand_dims(-sid, axis=1)
                    ew = np.expand_dims(wn, axis=1)
                    ev = np.expand_dims(v, axis=1)

                    #  Bid 3 x 3 x (i,r)diff
                    Bid = np.concatenate((e_sid,ev, ew), axis=1)

                    #update Bi for diffracted rays
                    Bi[:,:,udiff[0],udiff[1]] = Bid
                    ######
                    #Bo diffract
                    #####
                    w = np.cross(sod,vnormz, axisa=0, axisb=0, axisc=0)

                    w, nw = fix_colinear(w)
                    wn = w/nw

                    #wn = w/np.sqrt(np.sum(w*w, axis=0))
                    v = np.cross(wn, sod, axisa=0, axisb=0, axisc=0)

                    e_sod = np.expand_dims(sod, axis=1)
                    ew = np.expand_dims(wn, axis=1)
                    ev = np.expand_dims(v, axis=1)
                    #  Bod 3 x 3 x (i,r)diff
                    Bod = np.concatenate((e_sod,ev, ew), axis=1)

                    #update Bo for diffracted rays
                    Bo[:,:,udiff[0],udiff[1]] = Bod
                #################################
                # End of diffraction specific case
                ##################################


                #
                # pasting (Bo0,B,BiN)
                #

                # B : 3 x 3 x i x r

                Bo = np.concatenate((Bo0[:, :, np.newaxis, :], Bo), axis=2)
                Bi = np.concatenate((Bi, BiN[:, :, np.newaxis, :]), axis=2)

                # B : 3 x 3 x i x r

                self[k]['B'] = np.einsum('xv...,xw...->vw...', Bi, Bo)
                #self[k]['B'] = np.einsum('vx...,xw...->vw...', Bi, Bo)

                #BiN = np.array([si[:,-1,:], eth, eph])    # ndim x 3 x Nray
                #self[k]['BiN']=BiN
                # self[k]['B']=np.sum(self[k]['Bi'][:2,:2,np.newaxis]*self[k]['Bo'][np.newaxis,:2,:2],axis=1)


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
                #self[k]['aoa'] =  np.vstack((np.pi-thd, phd-np.pi))
                self[k]['aoa'] =  np.vstack((thd,phd))
                E = np.eye(2)[:,:,np.newaxis,np.newaxis]
                self[k]['B'] = np.dstack((E,E))
                ze = np.array([0])
                self[k]['rays'] = np.array(([[0]]))
                self[k]['nbrays'] = 1
                self[k]['rayidx'] = ze
                self.raypt = 1
                self._ray2nbi = ze
        self._luw = np.unique(luw).tolist()
        self.isbased = True

    def fillinter(self, L, append=False):
        """  fill ray interactions

        Parameters
        ----------

        L      : Layout
        append : Boolean
            If True append new rays to existing structure


        Notes
        -------

        This method adds the following members

        I : Interactions
        B : IntB
        B0 : IntB

        """

        # reinitialized ray pointer if not in append mode
        if not append:
            self.raypt = 0

        # stacked interactions
        I = Interactions(slab=L.sl)

        # rotation basis
        B  = IntB(slab=L.sl)
        B0 = IntB(slab=L.sl)

        # # LOS Interaction
        # Los = IntL()

        # Reflexion
        R = IntR(slab=L.sl)

        # Transmission
        T = IntT(slab=L.sl)

        # Diffraction
        D = IntD(slab=L.sl)

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
        #uslv = np.unique(L.sla[1:])
        uslv = L.sl.keys()
        #
        # add CEIL and FLOOR
        #
        #uslv = np.hstack((uslv, np.array(('CEIL', 'FLOOR'))))

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

        #to be specified and limited to used wedges
        if hasattr(self,'_luw'):
            D.dusl = dict.fromkeys(self._luw, np.array((), dtype=int))

        # transmission/reflection slab array
        tsl = np.array(())
        rsl = np.array(())
        # diffraction wedge list
        dw = np.array(())

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
                # distance in
                s_in = si[0:-1,:]
                # distance in
                s_out = si[1:,:]

                if 'diffvect' in self[k]:
                    dvec = self[k]['diffvect']
                    ldsl = self[k]['diffslabs']
                    dix = self[k]['diffidx']


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

                # ## index creation / already done in rays.locbas
                # ##################
                # # create index for retrieving interactions

                # # integer offset : total size idx

                # idxts = idxts + idx.size

                # idx = idxts + np.arange(ityp.size).reshape(np.shape(ityp),order='F')

                # nbray = np.shape(idx)[1]

                # self[k]['rays'] = idx
                # self[k]['nbrays'] = nbray
                # self[k]['rayidx'] = nbrayt + np.arange(nbray)
                # # create a numpy array to relate the ray index to its corresponding
                # # number of interactions

                # # _ray2nbi = np.ones((nbray))

                # #try:
                # #    self._ray2nbi=np.hstack((self._ray2nbi,_ray2nbi))
                # #except:
                # #    self._ray2nbi=_ray2nbi

                # #self._ray2nbi[self[k]['rayidx']]  = k
                # nbrayt = nbrayt + nbray
                # #self.raypt = self.raypt + self[k]['nbrays']

                idxf = self[k]['rays'].reshape(self[k]['rays'].size,order='F')
                #  (i+1)xr
                # 

                size2 = si[:, :].size
                nbray = self[k]['nbrays']
                # TODO
                # dirty fix
                # nbray is either an int or an array. why ?
                if type(nbray)==np.ndarray:
                    nbray=nbray[0]

                #  ,(i+1)xr
                # sif = si[:, :].reshape(size2,order='F') # TO BE REMOVE
                s_inf = s_in[:, :].reshape(ityp.size,order='F')
                s_outf = s_out[:, :].reshape(ityp.size,order='F')

                # 2x2,(i+1)xr

                #
                # self[k]['B'] 3 x 3 x i x r
                #
                # first unitary matrix (3x3xr)
                b0 = self[k]['B'][:,:,0,:]
                # first unitary matrix 1:
                # dimension i and r are merged
                b = self[k]['B'][:,:,1:,:].reshape(3, 3, size2-nbray,order='F')


                ## find used slab
                ##################
                # find slab type for the rnstr
                # nstrf is a number of slab
                # this is a problem for handling subsegment
                #

                # seek for interactions position
                ################################

                uD  = np.where((itypf == 1))[0]
                uR  = np.where((itypf == 2))[0]
                uT  = np.where((itypf == 3))[0]
                uRf = np.where((itypf == 4))[0]
                uRc = np.where((itypf == 5))[0]

                # assign floor and ceil slab
                ############################

                slT = [ L.Gs.node[x]['name'] for x in nstrf[uT] ]
                slR = [ L.Gs.node[x]['name'] for x in nstrf[uR] ]

                # WARNING
                # in future versions floor and ceil could be different for each cycle.
                # this information would be directly obtained from L.Gs
                # then the two following lines would have to be modified

                slRf = np.array(['FLOOR']*len(uRf))
                slRc = np.array(['CEIL']*len(uRc))


                # Fill the used slab
                #####################

                tsl = np.hstack((tsl, slT))
                rsl = np.hstack((rsl, slR, slRf, slRc))
                if 'diffvect' in self[k]:
                    dw = np.hstack((dw,self[k]['diffslabs']))
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
                # B.stack(data=b.T, idx=idxf)
                # B0.stack(data=b0.T,idx=self[k]['rayidx'])

                B.stack(data=b.T, idx=idxf)
                B0.stack(data=b0.T,idx=self[k]['rayidx'])
                ### Reflexion
                ############
                ### wall reflexion
                #(theta, s_in,s_out)

                R.stack(data=np.array((thetaf[uR], s_inf[uR], s_outf[uR])).T,
                        idx=idxf[uR])
                # floor reflexion
                R.stack(data=np.array((thetaf[uRf], s_inf[uRf], s_outf[uRf])).T,
                        idx=idxf[uRf])
                # ceil reflexion
                R.stack(data=np.array((thetaf[uRc], s_inf[uRc], s_outf[uRc])).T,
                        idx=idxf[uRc])

                # R.stack(data=np.array((thetaf[uR], sif[uR], sif[uR+1])).T,
                #         idx=idxf[uR])
                # # floor reflexion
                # R.stack(data=np.array((thetaf[uRf], sif[uRf], sif[uRf+1])).T,
                #         idx=idxf[uRf])
                # # ceil reflexion
                # R.stack(data=np.array((thetaf[uRc], sif[uRc], sif[uRc+1])).T,
                #         idx=idxf[uRc])

                ### sl[idxf[uT]]
                # Transmision
                ############
                # (theta, s_in,s_out)
                # T.stack(data=np.array((thetaf[uT], sif[uT], sif[uT+1])).T, idx=idxf[uT])
                T.stack(data=np.array((thetaf[uT], s_inf[uT], s_outf[uT])).T, idx=idxf[uT])

                ###
                #Diffraction
                #phi0,phi,si,sd,N,mat0,matN,beta
                #

                if 'diffvect' in self[k]:
                    # self[k]['diffvect'] = ((phi0,phi,beta,N) x (nb_rayxnb_interactions)   )
                    #si and so are stacked at the end of self[k]['diffvect']
                    #as well:
                    #data =  (6 x (nb_rayxnb_interactions) )
                    # ((phi0,phi,beta,N,sin,sout) x (nb_rayxnb_interactions) )
                    data = np.vstack((self[k]['diffvect'],s_inf[uD],s_outf[uD]))
                    D.stack(data=data.T,idx=self[k]['diffidx'])#idxf[uD])

            elif self.los:
                ze = np.array([0])
                #self[k]['rays'] = np.array(([[0]]))
                #self[k]['nbrays'] = 1
                #self[k]['rayidx'] = ze
                #self.raypt = 1
                #self._ray2nbi=ze
                B.stack(data=np.eye(3)[np.newaxis,:,:], idx=ze)
                B0.stack(data=np.eye(3)[np.newaxis,:,:],idx=ze)

        if len(tsl)>0:
            T.create_dusl(tsl)
        if len(rsl)>0:
            R.create_dusl(rsl)
        if len(dw)>0:
            D.create_dusl(dw)
        # create interactions structure
        self.I = I
        self.I.add([T, R, D])
        # create rotation base B
        self.B = B
        # create rotation base B0
        self.B0 = B0

        self.filled = True

    def eval(self,fGHz=np.array([2.4]),bfacdiv=False,ib=[]):
        """  field evaluation of rays

        Parameters
        ----------

        fGHz : array
            frequency in GHz
        ib : list of interactions block

        """

        #print 'Rays evaluation'

        self.fGHz=fGHz

        # evaluation of all interactions
        #
        # core calculation of all interactions is done here
        #

        self.I.eval(fGHz)

        # if np.isnan(self.I.I).any():
        #     pdb.set_trace()
        # evaluation of base B  (2x2)
        # B and B0 do no depend on frequency
        # just an axis extension (np.newaxis)
        #pdb.set_trace()

        # 1 x i x 3 x 3
        B  = self.B.data[np.newaxis,...]
        B  = B.swapaxes(2,3)
        # 1 x r x 3 x 3
        B0 = self.B0.data[np.newaxis,...]
        B0  = B0.swapaxes(2,3)

        # Ct : f x r x 3 x 3
        Ct = np.zeros((self.I.nf, self.nray, 3, 3), dtype=complex)

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

        # loop over group of interactions
        for l in ib:
            # ir : ray index

            ir = self[l]['rayidx']
            aoa[:,ir]=self[l]['aoa']
            aod[:,ir]=self[l]['aod']
            if l != 0:
                # l stands for the number of interactions
                r = self[l]['nbrays']
                # dirty fix should not be an array
                if type(r)==np.ndarray:
                    r = r[0]
                # reshape in order to have a 1D list of index
                # reshape ray index
                rrl = self[l]['rays'].reshape(r*l,order='F')
                # get the corresponding evaluated interactions
                #
                # reshape error can be tricky to debug.
                #
                # f , r , l , 2 , 2
                A = self.I.I[:, rrl, :, :].reshape(self.I.nf, r, l, 3, 3)
                # get the corresponding unitary matrix B
                # 1 , r , l , 2 , 2
                #Bl = B[:, rrl, :, :].reshape(self.I.nf, r, l, 2, 2,order='F')
                Bl = B[:, rrl, :, :].reshape(1, r, l, 3, 3)
                # get the first unitary matrix B0l
                B0l = B0[:,ir,:, :]
                # get alpha
                alpha = self.I.alpha[rrl].reshape(r, l,order='F')
                # # get gamma
                gamma = self.I.gamma[rrl].reshape(r, l,order='F')
                # # get si0
                si0 = self.I.si0[rrl].reshape(r, l,order='F')
                # # get sout
                sout = self.I.sout[rrl].reshape(r, l,order='F')

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
#                    if i == 0:
#                        pdb.set_trace()
#                        D0 = 1./si0[:,1]
#                        rho1 = si0[:,1]*alpha[:,i]
#                        rho2 = si0[:,1]*alpha[:,i]*gamma[:,i]
#                        D =np.sqrt(
#                         ( (rho1 ) / (rho1 + sout[:,i]) )
#                         *( (rho2) / (rho2 + sout[:,i])))
#                        D=D*D0
#                        rho1=rho1+(sout[:,i]*alpha[:,i])
#                        rho2=rho2+(sout[:,i]*alpha[:,i]*gamma[:,i])
#
#   ##                     gerer le loss
#                        if np.isnan(D).any():
#                            p=np.nonzero(np.isnan(D))[0]
#                            D[p]=1./sout[p,1]
#                    else :
#                        D=np.sqrt(
#                         ( (rho1 ) / (rho1 + sout[:,i]) )
#                         *( (rho2) / (rho2 + sout[:,i])))
#
#                        rho1=rho1+(sout[:,i]*alpha[:,i])
#                        rho2=rho2+(sout[:,i]*alpha[:,i]*gamma[:,i])
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

                #
                if bfacdiv:
                    Ct[:,ir, :, :] = Ct[:, ir, :, :]*1./(self[l]['dis'][np.newaxis, :, np.newaxis, np.newaxis])
                else:
                    Ct[:,ir, :, :] = Ct[:, ir, :, :]*1./(self[l]['dis'][np.newaxis, :, np.newaxis, np.newaxis])
                self.delays[ir] = self[l]['dis']/0.3
                self.dis[ir] = self[l]['dis']
        #
        # true LOS when no interaction
        #
        if self.los:
            Ct[:,0, :, :]= np.eye(3,3)[None,None,:,:]
            #self[0]['dis'] = self[0]['si'][0]
            # Fris
            Ct[:,0, :, :] = Ct[:,0, :, :]*1./(self[0]['dis'][None, :, None, None])
            self.delays[0] = self[0]['dis']/0.3
            self.dis[0] = self[0]['dis']


        # To be corrected in a future version
        #
        #  Ct : nf , Nray , theta , phi 
        #
        #  to 
        #
        #  Ct : Nray x nf , theta , phi 
        #
        Ct = np.swapaxes(Ct, 1, 0)

        #c11 = Ct[:,:,0,0]
        #c12 = Ct[:,:,0,1]
        #c21 = Ct[:,:,1,0]
        #c22 = Ct[:,:,1,1]

        c11 = Ct[:,:,1,1]
        c12 = Ct[:,:,1,2]
        c21 = Ct[:,:,2,1]
        c22 = Ct[:,:,2,2]


        #
        # Construction of the Ctilde propagation channel structure
        #
        Cn = Ctilde()

        # Cn.Cpp = bs.FUsignal(self.I.fGHz, c11)
        # Cn.Cpt = bs.FUsignal(self.I.fGHz, c12)
        # Cn.Ctp = bs.FUsignal(self.I.fGHz, c21)
        # Cn.Ctt = bs.FUsignal(self.I.fGHz, c22)
        Cn.Ctt = bs.FUsignal(self.I.fGHz, c11)
        Cn.Ctp = bs.FUsignal(self.I.fGHz, c12)
        Cn.Cpt = bs.FUsignal(self.I.fGHz, c21)
        Cn.Cpp = bs.FUsignal(self.I.fGHz, c22)

        Cn.nfreq = self.I.nf
        Cn.nray = self.nray
        Cn.tauk = self.delays
        Cn.fGHz = self.I.fGHz
        # r x 2
        Cn.tang = aod.T
        Cn.tangl = aod.T
        # r x 2
        #
        # recover angle of arrival convention 
        #
        Cn.rang  = np.hstack([np.pi-aoa.T[:,[0]],aoa.T[:,[1]]-np.pi]) 
        Cn.rangl = np.hstack([np.pi-aoa.T[:,[0]],aoa.T[:,[1]]-np.pi])
        # add aoa and aod

        self.evaluated = True

        return(Cn)


    def rayfromseg(self,ls):
        ''' DEPRECATED 
            use raysfromnstr instead
        '''
        DeprecationWarning('function name update: use raysfromnstr instead')
        return self.rayfromnstr(ls)

    def rayfromnstr(self,ls):
        """ returns the indexes of rays for a given interaction list
        """

        if not isinstance(ls,list):
            ls = [ls]

        lur = []
        for k in self:
            aib = self[k]['sig'][0,...]
            for i in ls :
                # import ipdb
                # ipdb.set_trace()
                ui, ur = np.where(aib == i)
                lur.extend(self[k]['rayidx'][ur].tolist())
        return np.sort(lur)

    def rayfromdelay(self,t0=0,t1=[]):
        """ returns the indexes of rays between 2 timestamps t0 and t1
        """
        if t1 == []:
            t1 = self.delays.max()
        u = np.where((self.delays>t0) & (self.delays<t1))[0]
        return u





    def ray2slab(self,L,ir):
        """ return the slabs for a given interaction index 


            Parameters
            ----------

            L : Layout
            ir : interaction block

            Returns
            -------

            numpy array of slabs strings at the shape (ir,r)
            ir : number of interactions ( of the interaction block)
            r : number of rays

        """

        v=np.vectorize( lambda t: L.Gs.node[t]['name'] if (t!=0) and (t>0) else '_')
        return v(self[ir]['sig'][0])


    def ray(self, r):
        """ returns the index of interactions of r

        Parameters
        ----------

        r : integer
            ray index

        Returns
        -------

        ir : nd.array
            index of interactions of r

        Examples
        --------

        """
        raypos = np.nonzero(self[self._ray2nbi[r]]['rayidx'] == r)[0]
        return(self[self._ray2nbi[r]]['rays'][:,raypos][:,0])

    def ir2a(self,ir):
        """ index ray 2 address ray

        Parameters
        ----------
        ir : integer

        Returns
        -------
        (ni,ux) : tuple address (group of interactions, index)

        """
        assert ir < self.nray, "wrong ray index"
        ni = self._ray2nbi[ir]
        ur = np.where(self[ni]['rayidx']==ir)[0][0]
        return(ni,ur)

    def a2ir(self,t):
        """  address ray 2 index ray

        Parameters
        ----------
        t = (ni,ux) : tuple address (group of interactions, index)
            ray address

        Returns
        -------
        ir : integer
            ray index

        """
        assert t[0] in self.keys(), "wrong number of interactions"
        ir = self[t[0]]['rayidx'][t[1]]
        return(ir)


    def ray2ityp(self,r):
        """ return interaction type for a given ray


        Parameters
        ----------

        r : integer
        ray index


        Returns
        -------
        
        lt : list
            list of type of interactions

        """

        di = {1:'D',2:'R',3:'T',4:'R',5:'R'}
        sig = self.ray2sig(r)
        sig  = sig[1,1:-1]
        return [di[s] for s in sig]


    def ray2nbi(self,r):
        """ Get interaction block/number of interactions of a given ray

        Parameters
        ----------

        r : integer
            ray index

        Returns
        -------

        nbi : int
            interaction block number
        """
        i = self._ray2nbi[r]
        return i 

    def ray2iidx(self,ir):
        """ Get interactions index of a given ray

        Parameters
        ----------

        ir : integer
            ray index

        Returns
        -------

        iidx : array
            interaction index 
        """
        unbi = self.ray2nbi(ir)
        ur = np.where(self[unbi]['rayidx']==ir)[0]
        return self[unbi]['rays'][:,ur]


    def ray2sig(self,ir):
        """ get signature to corresponding ray
        """
        unbi = self.ray2nbi(ir)
        ur = np.where(self[unbi]['rayidx']==ir)[0]
        return self[unbi]['sig'][:,:,ur].squeeze()

    def ray2sig2d(self,ir):
        """ get signature to corresponding ray
        """
        sig = self.ray2sig(ir)
        sig = sig.squeeze()
        sig = sig[:,1:-1] # remove extremal 0
        unfc = np.where(sig[1,:]<4)[0]# index floor cell
        sig2d = sig[:,unfc]
        return sig2d

    def ray2inter(self,ir,L,Si):
        """ get interaction list (Gi style) from a ray

        Parameters
        ----------

        ir : ray index
        L : Layout
        Si : Signatures object

        """
        sig = self.ray2sig2d(ir)
        return Si.sig2inter(L,sig)


    def slab_nb(self, ir):
        """ returns the slab numbers of r

        Parameters
        ----------

        ir : integer
            ray index

        Returns
        -------

        isl : slabs number


        """

        raypos = np.nonzero(self[self._ray2nbi[ir]]['rayidx'] == ir)[0]
        return(self[self._ray2nbi[ir]]['sig'][0,1:-1,raypos[0]])

    def vis(self,ir,L):
        typ = ['Tx'] + self.typ(ir) + ['Rx']
        slab_nb = self.slab_nb(ir)
        slab_nb = np.insert(slab_nb,0,0)
        slab_nb = np.insert(slab_nb,len(slab_nb),0)
        nbi = self._ray2nbi[ir]
        raypos = np.nonzero(self[nbi]['rayidx'] == ir)[0]
        pt = self[nbi]['pt'][:,:,raypos]
        tz  = pt[2].ravel()
        slab = [ L.Gs.node[x]['name'] for x in slab_nb if x > 0]
        st = ''
        for t in typ:
            st = st + t+'      ' 
        print(st)
        st = ''
        for s in slab_nb:
            st = st + str(s)+'     ' 
        print(st)
        st = ''
        for z in tz:
            st = st + str(z)+'     ' 
        print(st)
        print(slab)

    def typ(self, ir,fromR=True):
        """ returns interactions list type of a given ray

        Parameters
        ----------

        ir : integer
            ray index
        fromR : bool
            True : get information from signature in R
            False: get information in R.I

        """
        #
        # In this function we can see that the ceil and floor
        # are hard coded as reflection. This is going to evolve
        # for implementation of multi floor
        #
        if fromR:
            di = {0:'L',1:'D',2:'R',3:'T',4:'R',5:'R'}
            nbi = self._ray2nbi[ir]
            raypos = np.nonzero(self[nbi]['rayidx'] == ir)[0]
            inter = self[nbi]['sig'][1,1:-1,raypos][0]
            return [di[i] for i in inter]
        else:
            a = self.ray(r)
            return(self.I.typ[a])

    def dump(self,ir,L,ifGHz=0,filename='dumpray.ray'):
        """ dump the full information of a ray in a file
        """
        nbi = self._ray2nbi[ir]
        ur = np.where(self[nbi]['rayidx']==ir)[0][0]
        fd=open(filename,'w')
        fd.write('ray #'+str(ir)+'\n')
        fd.write(str(ur)+ ' th ray from the group of ' + str(nbi)+' Interactions' +'\n')
        cy_a = L.pt2cy(self.pTx)
        cy_b = L.pt2cy(self.pRx)

        #fd.write('Tx #'+str(self.pTx)+'\n')
        #fd.write('Rx #'+str(self.pRx)+'\n')
        if self.evaluated:
            ray = self.ray(ir)
            typ = self.typ(ir)
            slabnb = self.slab_nb(ir)
            fd.write('   ray #'+str(ray)+'\n')
            #fd.write('   typ #'+str(typ)+'\n')
            fd.write('   slab #'+str(slabnb)+'\n')
        for k in range(nbi+2):
            if k==0:
                fd.write('Tx  :        ')
            elif k==(nbi+1):
                fd.write('Rx  :        ')
            else:
                six = slabnb[k-1]
                if six==0:
                    slabname='FLOOR'
                    cyc =[-2,-3]
                else:
                    slabname = L.Gs.node[six]['name']
                    cyc = L.Gs.node[six]['ncycles']
                if typ[k-1]=='T':
                    fd.write('T '+slabname +'       ('+str(six)+','+str(cyc[0])+','+str(cyc[1])+')')
                if typ[k-1]=='R':
                    fd.write('R '+slabname +'       ('+str(six)+',)')
                if typ[k-1]=='D':
                    fd.write('D ('+str(six)+') :')

            fd.write(str(self[nbi]['pt'][:,k,ur])+'\n' )
            if k==0:
                fd.write('  '+str(cy_a)+'\n')
            elif k==(nbi+1):
                fd.write('  '+str(cy_b)+'\n')
            if k==0:
                for l in range(3):
                    if l<2:
                        fd.write('\t'+str(self[nbi]['Bo0'][l,:,ur])
                     +'\t'+str(self[nbi]['B'][l,:,0,ur])+'\n')
                    else:
                        fd.write('\t'+str(self[nbi]['Bo0'][l,:,ur]) +'\n')
            elif k==(nbi+1):
                for l in range(3):
                    fd.write('\t'+str(self[nbi]['BiN'][l,:,ur])+'\n')
            else:
                for l in range(3):
                    if l<2:
                        fd.write('\t'+str(self[nbi]['Bi'][l,:,k-1,ur])+'\t'+
                              str(self[nbi]['Bo'][l,:,k-1,ur])
                         +'\t'+str(self[nbi]['B'][l,:,k-1,ur])+'\n')
                    else:
                        fd.write('\t'+str(self[nbi]['Bi'][l,:,k-1,ur])+'\t'+
                              str(self[nbi]['Bo'][l,:,k-1,ur])+'\n')


        fd.close()


    def info(self,ir,ifGHz=0,bB=True,matrix=False):
        """ provides information for a given ray r

        Parameters
        ----------

        ir : int
            ray index
        ifGHz : int
            frequency index
        bB: boolean
            display Basis
        matrix :
            display matrix 
        """

        if self.evaluated:
            print('-------------------------')
            print('Informations of ray #', ir)
            print('-------------------------\n')

            ray = self.ray(ir)
            typ = self.typ(ir)
            slabnb = self.slab_nb(ir)
            # if there is a diffraction, phi0, phi, beta are shown
            if 'D' in typ:
                diff =True
                print('{0:5} , {1:4}, {2:10}, {3:7}, {4:7}, {5:10}, {6:10}, {7:4}, {8:4}, {9:4}'\
                        .format('Index',
                                'type',
                                'slab', 
                                'nstr' ,
                                'th(rad)',
                                'alpha',
                                'gamma2',
                                'phi0',
                                'phi',
                                'beta'))
            else :
                diff =False
                print('{0:5} , {1:4}, {2:10}, {3:7}, {4:7}, {5:10}, {6:10}'\
                     .format('Index',
                        'type',
                        'slab',
                        'nstr',
                        'th(rad)',
                        'alpha',
                        'gamma2'))
            print('{0:5} , {1:4}, {2:10}, {3:7}, {4:7.2}, {5:10.2}, {6:10.2}'\
                  .format(ir, 'B0','-', '-', '-', '-', '-'))

            for iidx, i in enumerate(typ):
                # import ipdb
                # ipdb.set_trace()
                if i == 'T' or i == 'R' or i =='D':
                    I = getattr(self.I, i)
                    for slab in I.dusl.keys():
    #                    print slab
                        midx = I.dusl[slab]
    #                    print midx
                        Iidx = np.array((I.idx))[midx]

                        if i != 'D':
                            th = I.data[I.dusl[slab], 0]
                            gamma = I.gamma[midx]
                            alpha = I.alpha[midx]
                        else : 
                            # from IPython.core.debugger import Tracer
                            # Tracer()()
                            th=['-']*max(max(Iidx),1)
                            gamma = ['NC']*max(max(Iidx),1)
                            alpha = ['NC']*max(max(Iidx),1)
                            udiff = np.where(self.I.D.idx==ray[iidx])[0]
                            phi0 = self.I.D.phi0[udiff][0]
                            phi=self.I.D.phi[udiff][0]
                            beta=self.I.D.beta[udiff][0]
                        for ii, Ii in enumerate(Iidx):
                            if Ii == ray[iidx]:
                                if i=='D': 
                                    print('{0:5} , {1:4}, {2:10}, {3:7}, {4:7.2}, {5:10}, {6:10}, {7:3.4}, {8:3.4}, {9:3.4}'\
                                    .format(Ii, i, slab, slabnb[iidx], th[ii], alpha[ii], gamma[ii],phi0,phi,beta))
                                else:
                                    print('{0:5} , {1:4}, {2:10}, {3:7}, {4:7.2}, {5:10.2}, {6:10.2}'\
                                    .format(Ii, i, slab, slabnb[iidx], th[ii], alpha[ii], gamma[ii]))

                    else:
                        if bB:
                            print('{0:5} , {1:4}, {2:10}, {3:7}, {4:7.2}, {5:10.2}, {6:10.2}'.format(ray[iidx], 'B', '-', '-', '-', '-', '-'))
                #              print '{0:5} , {1:4}, {2:10}, {3:7}, {4:10}, {5:10}'.format(ray[iidx], i, '-', '-', '-', '-')

            if matrix:
                print('\n----------------------------------------')
                print(' Matrix of ray #', ir, 'at f=', self.I.fGHz[ifGHz])
                print('----------------------------------------')
                lmat = []
                ltran = []
                if bB:
                    print('rotation matrix#', 'type: B0')

                    B0 = self.B0.data[ir,:,:]
                    addr = self.ir2a(ir)
                    Bo0 = self[addr[0]]['Bo0'][:,:,addr[1]] 
                    Bi1 = self[addr[0]]['Bi'][:,:,0,addr[1]] 
                    U  = np.dot(Bi1.T,Bo0)
                    assert np.allclose(B0,U) 
                    lmat.append(B0)
                    ltran.append(B0)
                    print(B0)
                for iidx, i in enumerate(typ):
                    print('interaction #', ray[iidx], 'type:', i)
                    # f x l x 2 x 2
                    I = self.I.I[ifGHz, ray[iidx], :, :]
                    print(I)
                    lmat.append(I)

                    if bB:
                        print('rotation matrix#',[ray[iidx]], 'type: B')
                        B = self.B.data[ray[iidx], :, :]
                        print(B) 
                        lmat.append(B)
                        ltran.append(B)
                # evaluate matrix product
                PM0=np.eye(3)
                PM1=np.eye(3)
                for m in lmat[::-1]:
                    PM0=np.dot(PM0,m)
                for m in ltran[::-1]:
                    PM1=np.dot(PM1,m)
                print("matrix product with interactions (dB)")
                print(20*np.log10(np.abs(PM0[1,1])),'  ',20*np.log10(np.abs(PM0[1,2])))
                print(20*np.log10(np.abs(PM0[2,1])),'  ',20*np.log10(np.abs(PM0[2,2])))
                print("matrix product without interactions (dB)")
                print(20*np.log10(np.abs(PM1[1,1])),'  ',20*np.log10(np.abs(PM1[1,2])))
                print(20*np.log10(np.abs(PM1[2,1])),'  ',20*np.log10(np.abs(PM1[2,2])))
                return(PM0)

            else:
                print('\nto display matrix, use matrix=True on call')
        else:
            print('Rays have not been evaluated yet')

    def signature(self, u , typ='full'):
        """ extract ray signature

        Parameters
        ----------

        u : tuple orr int 
            if tuple addr 
            if int index

        Returns
        -------

        sig : ndarray

        Notes
        -----

        Signature of a ray is store as a member

        r[nint]['sig']

        """
        if type(u)==tuple:
            addr = u 
        else:
            addr = self.ir2a(u) 
        if typ=='full':
            sig = self[addr[0]]['sig'][:,:,addr[1]]
        else:
            pass
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

    def _show3(self,L=[],rlist=[],newfig=False,cmap='hot',**kwargs):
        """ plot 3D rays in environment using Mayavi

        Parameters
        ----------

        L : Layout object
            Layout to be displayed
        rlist : list
            list of index rays
        newfig : boolean (default: False)
            if true create a new mayavi figure
            else : use the current
        ER: Ray energy

        """

        if newfig:
            mlab.clf()
            f = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
        else :
            f = mlab.gcf()
            # view=mlab.view()


        if L != []:
            try:
                L._filename
            except:
                raise NameError('L argument must be a layout object')

            L._show3()

        if 'ER' in kwargs:
            ER = kwargs['ER']
            color_range = np.linspace( 0, 1., len(ER))#np.linspace( 0, np.pi, len(ER))
            uER = ER.argsort()[::-1]
            colors= color_range[uER]

        if rlist ==[]:
            nbi = self.keys()
            for i in nbi:
                r = range(np.shape(self[i]['pt'])[2])
                ridx = self[i]['rayidx']
                # number of rays
                nbr = len(r) 
                # current number of interactions
                cnbi = i + 2

                # import ipdb
                # ipdb.set_trace()
                pt = self[i]['pt'][:,:,r].reshape(3,cnbi*nbr,order='F')
                l0 = np.array([np.arange(0,cnbi-1)+k*cnbi for k in range(nbr)]).ravel()
                l1 = l0+1
                connection = np.vstack((l0,l1)).T
                if 'ER' in kwargs:
                    rc = np.repeat(colors[ridx],cnbi)
                    rc[::cnbi]=0
                    src = mlab.pipeline.scalar_scatter(pt[0,:], pt[1,:], pt[2,:],rc,colormap=cmap)
                else: 
                    src = mlab.pipeline.scalar_scatter(pt[0,:], pt[1,:], pt[2,:])

                src.mlab_source.dataset.lines=connection
                src.update()
                lines = mlab.pipeline.stripper(src)
                mlab.pipeline.surface(lines,opacity=0.5,colormap=cmap)
                f.children[-1].name='Rays with ' + str(i) + 'interactions'
        else :

            nbi = self._ray2nbi[rlist]
            nr = np.array((nbi,rlist))
            unb = np.unique(nr[0,:])
            unr = {int(i):np.where(nr[0,:]==i)[0] for i in unb}

            for i in unb:
                raynb = (nr[1,unr[i]]).astype(int)
                nbr=len(raynb)
                ptidx = [np.where(self[i]['rayidx']==x)[0][0] for x in raynb]
                # current number of interactions
                cnbi = i + 2
  
                pt = self[i]['pt'][:,:,ptidx].reshape(3,cnbi*nbr,order='F')

                # lines = np.arange(cnbi*nbr).reshape(cnbi,nbr)
                lines = np.arange(cnbi*nbr).reshape(nbr,cnbi)

                # mesh = tvtk.PolyData(points=pt.T, polys=lines)
                mesh = tvtk.PolyData(points=pt.T, polys=lines)
                mlab.pipeline.surface(mlab.pipeline.extract_edges(mesh),
                                                     color=(0, 0, 0), )
                f.children[-1].name='Rays with ' + str(int(i)) + 'interactions'

        # mlab.view(view[0],view[1],view[2],view[3])
        return(f)

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
            L._filename
        except:
            raise NameError('L argument must be a layout object')

        if not centered:
            pg=np.array([[0],[0],[0]])

        strucname= L._filename.split('.')[0]
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


if __name__ == "__main__":
    doctest.testmod()
