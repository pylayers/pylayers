#!/usr/bin/python
# -*- coding: latin1 -*-
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
import pylayers.signal.bsignal as bs

class Rays(dict):
    """ A set af rays

    Members
    -------

    'rays'   : 
    'nbrays' :
    'rayidx' :
    'sig'    :
    'pt'     :
    'alpha'  :


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

    Rays object is obtained from a signature.
    It is a container for a set of rays between a source 
    and a target point defining a radio link.

    Once a Rays object has been obtained in 2D, it is transform 
    in 3D via the **to3D** method. This method has two parameters : 
    the height from floor to ceil, and the number N of 
    multiple reflections to take into account. 

    Once the 3d rays have been calculated, the  
    the local basis are evaluated along those rays. This is
    done through the **locbas** method

    Once the local basis have been calculated the different
    interactions along rays can be informed via the **fillinter**
    method.

    Once the interaction are informed the field along rays can 
    be evaluated via the **eval** method
    """
    def __init__(self, pTx, pRx):
        self.pTx = pTx
        self.pRx = pRx
        self.nray = 0
        self.los=False

    def __repr__(self):
        s = ''
        ni = 0
        nl = 0
        try:
            for k in self:
                r = self[k]['rayidx']
                nr = len(r)
                s = s + str(k)+' / '+str(nr)+ ' : '+str(r)+'\n'
                ni = ni + nr*k
                nl = nl + nr*(2*k+1)
            s = s + '-----'+'\n'
            s = s+'ni : '+str(ni)+'\n'
            s = s+'nl : '+str(nl)+'\n'
        except:
            return(s)

        return(s)

    def show(self,L,i=-1,r=-1):
        """  plot 2D rays within the simulated environment


        Parameters
        ----------

        L : Layout
        i : list or -1 (default = all groups)
            list of interaction group number 
        r : list or -1 (default = all rays) 
            list of index of ray in interaction group 

        """

        fig = plt.figure()
        ax = fig.add_subplot(111)
        L.showG(fig=fig,ax=ax,graph='s',nodes=False)
        ax.plot(self.pTx[0], self.pTx[1], 'or')
        ax.plot(self.pRx[0], self.pRx[1], 'og')
        if i==-1:
            lgrint = self.keys()
        else:
            lgrint = i
            

        for i in lgrint:
            if r==-1:
                lray = range(len(self[i]['pt'][0, 0, :]))
            else:
                lray = r
            for j in lray: 
                ray = np.hstack((self.pTx[0:2].reshape((2, 1)),
                                 np.hstack((self[i]['pt'][0:2, :, j],
                                 self.pRx[0:2].reshape((2, 1))))
                                 ))
                ax.plot(ray[0, :], ray[1, :], alpha=0.6, linewidth=1.)
                ax.axis('off')

    def mirror(self, H=3, N=1):
        """ mirror

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

        lsss : list od segment with a subsegment    

        returns
        -------

        r3d : Rays 

        Notes
        -----


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
                self[i]['alpha'][j, :] = np.sum(si[0:j+1, :], axis=0) \
                        /np.sum(si, axis=0)
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
            # get  2D signature
            sig = self[k]['sig']            
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
                Nint = len(d[l])            # number of additional interaction
                if Nint > 0:                # if new interaction ==> need extension 
                    # a1e : extended horizontal+vertical parameterization
                    a1e = np.concatenate((a1, d[l].reshape(len(d[l]), 1)*\
                                          np.ones((1, Nrayk)))) 
                    # get sorted indices
                    ks = np.argsort(a1e, axis=0) 
                    # a1es : extended sorted horizontal + vertical parameterization
                    a1es = np.sort(a1e, axis=0)  
                    # prepare an extended sequence of points ( ndim x  (Nint+k+2) x Nrayk )
                    ptee = np.hstack((pte, np.zeros((3, Nint, Nrayk))))     
                    
                    #
                    # Boolean ceil/floor detector 
                    #
                    # u is 4 (floor interaction ) 
                    #      5 (ceil interaction ) 
                    #  depending on the vertical pattern l. The observation is
                    #  that negative l corresponds always to last reflexion on
                    #  the floor and positive l corresponds to last reflexion
                    #  on the ceil 
                    #
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
                    # A signatue is compose of two lines 
                    # esigs sup line : interaction number 
                    # esigi inf line : interaction type 
                    #
                    esigs = np.zeros((1, Nint, Nrayk), dtype=int)
                    esigi = u.reshape(1, Nint, 1)* np.ones((1, 1, Nrayk), dtype=int)
                    # esig : extension of the signature               
                    esig = np.vstack((esigs, esigi))
                    # sige : signature extended  ( 2 x (Nint+k+2) x Nrayk ) 
                    sige = np.hstack((sig, esig))
                    # 2 x (Nint+k+2) x Nrayk
                    #
                    # Now come the time to sort extended sequence of points
                    # and extended sequence of signatures with the sorting
                    # index ks obtained from argsort of merge parametization  
                    #
                    # sequence of extended sorted points
                    ptees = ptee[:, ks, range(Nrayk)]     
                    siges = sige[:, ks, range(Nrayk)]   
                    # extended and sorted signature
                    iint_f, iray_f = np.where(siges[ 1, :] == 4)  # floor interaction
                    iint_c, iray_c = np.where(siges[ 1, :] == 5)  # ceil interaction
                        
                    #
                    # find the list of the previous and next point around the
                    # new ceil or floor point. The case of successive ceil or
                    # floor reflexion make 
                    # Tous les points précédents qui ne sont pas des Ceils ou
                    # des floors et tous les points suivants qui ne sont pas
                    # des points de réflexion ceil ou floor 
                    # Afin de tenir compte du rayon du groupe d'interaction
                    # concerné il faut passer un tuple qui concatène la valeur
                    # de l'indice d'interaction floor ou ceil et l'indice de
                    # rayons du groupe associé (d'ou le zip) 
                    #
                    # Cette séquence d'instruction fixe le bug #133
                    #
                    # Antérieurement il y avait une hypothèse de succesion
                    # immediate d'un point 2D renseigné.
                    #
                    iintm_f = map(lambda x : np.where( (siges[1,0:x[0],x[1]]<>4) & (siges[1,0:x[0],x[1]]<>5))[0][-1], zip(iint_f,iray_f))
                    iintp_f = map(lambda x : np.where( (siges[1,x[0]:,x[1]]<>4) & (siges[1,x[0]:,x[1]]<>5))[0][0]+x[0], zip(iint_f,iray_f))
                    iintm_c = map(lambda x : np.where( (siges[1,0:x[0],x[1]]<>4) & (siges[1,0:x[0],x[1]]<>5))[0][-1], zip(iint_c,iray_c))
                    iintp_c = map(lambda x : np.where( (siges[1,x[0]:,x[1]]<>4) & (siges[1,x[0]:,x[1]]<>5))[0][0]+x[0], zip(iint_c,iray_c))
                    
                   
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
                    siges = sig

                #   ptes (3 x i+2 x r ) 
                u   = map(lambda x: list(np.where(siges[0,:,:]==x)),L.lsss)[0]
                v   = [2*np.ones(len(u[0]),dtype=int)]+u
                w   = [0*np.ones(len(u[0]),dtype=int)]+u
                # zss : height of interaction on subsegments
                zss = ptees[v]
                # structure index of corresponding subsegments 
                nstrs = siges[w]
                #print "nstrs: ",nstrs
                #print "zss:",zss
                #
                # Determine which subsegment has been intersected 
                # k = 0 : no subsegment intersected
                zinterval = map(lambda x: L.Gs.node[x]['ss_z'],nstrs)
                tab = map (lambda x: filter(lambda z: ((z[0]<x[1]) &
                                                       (z[1]>x[1])),x[0]),zip(zinterval,zss))
                #print tab
                def findindex(x):
                    if len(x[1])>0:
                        k = x[0].index(x[1][0])+1
                        return(k)
                    else:
                        return(0)

                indexss = map(findindex,zip(zinterval,tab))
                indexnew = map(lambda x: x[0] if x[1]==0 else 1000000+100*x[0]+x[1]-1,zip(nstrs,indexss))
                # update signature
                siges[w] = indexnew
                #print "indexss:",indexss
                #print "indexnew:",indexnew
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
                    r3d[k+Nint]['pt'] = np.dstack((r3d[k+Nint]['pt'], ptees))
                    r3d[k+Nint]['sig'] = np.dstack((r3d[k+Nint]['sig'], siges))
                except:
                    r3d[k+Nint] = {}
                    # r3d[k+Nint]['alpha'] = a1es
                    # r3d[k+Nint]['ks'] = ks
                    r3d[k+Nint]['pt'] = ptees
                    r3d[k+Nint]['sig'] = siges


        if self.los :
            r3d[0]={}
            r3d[0]['sig']=np.zeros((2,2,1))
            r3d[0]['pt']=np.zeros((3,2,1))
            r3d[0]['pt'][:,0,:]=tx[:,np.newaxis]
            r3d[0]['pt'][:,1,:]=rx[:,np.newaxis]

        return(r3d)

    def length(self,typ=2):
        """ calculate length of rays

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
        """
        Parameters
        ----------

        L : Layout

        """

        #
        # extract normal in np.array
        #

        # nsegment x 3
        norm = np.array(nx.get_node_attributes( L.Gs, 'norm').values())

        # nsegment x k
        key = np.array(nx.get_node_attributes( L.Gs, 'norm').keys())

        nmax = max(L.Gs.node.keys())
        mapping = np.zeros(nmax+1, dtype=int)
        mapping[key] = np.arange(len(key), dtype=int)

        #
        # Structutre number : nstr
        #   the structure number is < 0 for points 
        #                           > 0 for segments
        # A segment can have several subsegments (until 100)
        #  nstrs is the nstr of the segment if sous segment : 1000000 + nstr*100 + k-1
        #  nstr  is the glabal which allows to recover the slab values 
        #
        for k in self:
            if k <> 0:
                nstr = self[k]['sig'][0, 1:-1, :]      # nint x nray
                ityp = self[k]['sig'][1, 1:-1, :]      # nint x nray
                # nstr of underlying segment
                #print nstr
                uss   = np.where(nstr>1000000)
                #print uss
                nstrs = copy.copy(nstr)
                nstrs = nstr
                if len(uss)>0:
                    nstrs[uss] = nstr[uss]/100-10000
                #    print nstr
                #    print nstrs
                #pdb.set_trace()
                nray = np.shape(nstr)[1]

                uwall = np.where((ityp == 1) | (ityp == 2))
                udiff = np.where((ityp == 3))
                ufloor = np.where((ityp == 4))
                uceil = np.where((ityp == 5))

                nstrwall  = nstr[uwall[0], uwall[1]]     # nstr of walls
                nstrswall = nstrs[uwall[0], uwall[1]]   # nstrs of walls

                self[k]['nstrwall']  = nstrwall    # store
                self[k]['nstrswall'] = nstrswall   # store

                self[k]['norm'] = np.zeros((3, k, nray))   # 3 x int x nray

                #
                # Warning : The following commented line assumes that all the segment number are contiguous
                # self[k]['norm'][:,uwall[0],uwall[1]] = norm[nstrwall-1,:].T
                #

                # norm : 3 x i x r
                self[k]['norm'][:, uwall[0], uwall[1]] = norm[mapping[nstrswall],:].T
                self[k]['norm'][2, ufloor[0], ufloor[1]] = np.ones(len(ufloor[0]))
                self[k]['norm'][2, uceil[0], uceil[1]] = -np.ones(len(uceil[0]))

                v = self[k]['pt'][:, 1:, :]-self[k]['pt'][:, 0:-1, :]
                lsi = np.sqrt(np.sum(v*v, axis=0))
                si = v/lsi             # ndim , nint - 1 , nray

                # vsi : 3 x (i+1) x r
                self[k]['vsi'] = si

                # si : (i+1) x r
                self[k]['si'] = lsi
                self[k]['dis'] = np.sum(lsi,axis=0)

                vn = self[k]['norm']

                # s_in : 3 x i x r
                s_in = si[:, 0:-1, :]

                # s_out : 3 x i x r
                s_out = si[:, 1:, :]

                #
                # AOD (rad)
                #

                # th : ,r
                th = np.arccos(si[2, 0, :])

                # ph : ,r
                ph = np.arctan2(si[1, 0, :], si[0, 0, :])

                # aod : 2 x r  (radians)
                self[k]['aod'] = np.vstack((th, ph))

                # eth : 3 x r
                eth = np.array([np.cos(th) * np.cos(ph),
                               np.cos(th) * np.sin(ph),
                                -np.sin(th)])
                # eph : 3 x r
                eph = np.array([-np.sin(ph),
                                np.cos(ph),
                                np.zeros(len(ph))])

                # Bo0 : 3 x 2 x r
                Bo0 = np.concatenate((eth[:, np.newaxis, :],
                                      eph[:, np.newaxis, :]), axis=1)

                self[
                    k]['Bo0'] = np.concatenate((si[:, 0, np.newaxis, :], eth[:, np.newaxis, :],
                                                eph[:, np.newaxis, :]), axis=1)

                #
                # scalar product si . norm
                #

                scpr = np.sum(vn*s_in, axis=0)
                self[k]['scpr'] = scpr
                self[k]['theta'] = np.arccos(abs(scpr))  # *180/np.pi

                #
                # Warning need to handle singular case when s_in//vn
                #
                w = np.cross(s_in, vn, axisa=0, axisb=0, axisc=0)
                wn = w/np.sqrt(np.sum(w*w, axis=0))
                v = np.cross(wn, s_in, axisa=0, axisb=0, axisc=0)

                es_in = np.expand_dims(s_in, axis=1)
                ew = np.expand_dims(wn, axis=1)
                ev = np.expand_dims(v, axis=1)

                #  Bi 3 x 2 x i x r
                Bi = np.concatenate((ew, ev), axis=1)
                # self[k]['Bi'] = np.concatenate((es_in,ew,ev),axis=1)

                w = np.cross(s_out, vn, axisa=0, axisb=0, axisc=0)
                wn = w/np.sqrt(np.sum(w*w, axis=0))
                v = np.cross(wn, s_out, axisa=0, axisb=0, axisc=0)

                es_out = np.expand_dims(s_out, axis=1)
                ew = np.expand_dims(wn, axis=1)
                ev = np.expand_dims(v, axis=1)

                #  Bi 3 x 2 x i x r
                Bo = np.concatenate((ew, ev), axis=1)
                # self[k]['Bo'] = np.concatenate((es_out,ew,ev),axis=1)

                #
                # AOA (rad)
                #

                # th : ,r
                th = np.arccos(si[2, -1, :])

                # th : ,r
                ph = np.arctan2(si[1, -1, :], si[0, -1, :])

                # aoa : 2 x r  (radians)
                self[k]['aoa'] = np.vstack((th, ph))
                eth = np.array([np.cos(th) * np.cos(ph),
                               np.cos(th) * np.sin(ph),
                                -np.sin(th)])
                eph = np.array([-np.sin(ph),
                                np.cos(ph),
                                np.zeros(len(ph))])
                # Bo0 : 3 x 2 x r
                BiN = np.concatenate((eth[:, np.newaxis, :],
                                      eph[:, np.newaxis, :]), axis=1)

                # self[k]['BiN'] = np.concatenate((si[:,-1,np.newaxis,:],eth[:,np.newaxis,:],
                #                                    eph[:,np.newaxis,:]),axis=1)

                #
                # pasting (Bo0,B,BiN)
                #
                # B : 3 x 2 x i x r

                Bo = np.concatenate((Bo0[:, :, np.newaxis, :], Bo), axis=2)
                Bi = np.concatenate((Bi, BiN[:, :, np.newaxis, :]), axis=2)

                # B : 2 x 2 x i x r

                self[k]['B'] = np.einsum('xv...,xw...->vw...', Bo, Bi)

                # BiN = np.array([si[:,-1,:], eth, eph])    # ndim x 3 x Nray
                # self[k]['BiN']=BiN
                # self[k]['B']=np.sum(self[k]['Bi'][:2,:2,np.newaxis]*self[k]['Bo'][np.newaxis,:2,:2],axis=1)

            # if los exists
            else :
                self[k]['nstrwall'] = np.array(())
                self[k]['norm'] = np.array(())
                self[k]['vsi'] = np.array(())
                # si : (i+1) x r
                si = np.sqrt(np.sum((self[0]['pt'][:,0]-self[0]['pt'][:,1])**2,axis=0))
                self[k]['si'] = np.vstack((si,0.))
                self[k]['dis'] = np.array((si))
                self[k]['aod'] = np.array(())
                self[k]['Bo0'] = np.array(())
                self[k]['scpr'] = np.array(())
                self[k]['theta'] = np.zeros((1,1))
                self[k]['aoa'] = np.array(())
                E=np.eye(2)[:,:,np.newaxis,np.newaxis]
                self[k]['B'] = np.dstack((E,E))

    def fillinter(self, L):
        """  docstring for fillinter

        Parameters
        ----------

        L : Layout
        """

        # stacked interactions
        I = Interactions()

        # rotation basis
        B = IntB()
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
        slv = nx.get_node_attributes(L.Gs, "name").values()
        slk = nx.get_node_attributes(L.Gs, "name").keys()

        # find all material used in simulation
        uslv = np.unique(slv)
        uslv = np.hstack((uslv, np.array(('CEIL', 'FLOOR'))))

        # create reverse dictionnary with all material as a key
        # and associated point/segment as a value

        dsla = {}
        for s in uslv:
            dsla[s] = np.where(s == np.array(slv))[0]

        nmax = max(L.Gs.node.keys())
        sla = np.zeros((nmax+1), dtype='S20')

        # array type str with more than 1 character
        # warning use zeros instead of empty because slab zero
        # is virtually used before assigning correct slab to ceil and floor

        #
        # sla is an array of string.
        # each value of Gs node is the index of the corresponding slab
        #

        sla[slk] = np.array(slv)

        R.dusl = dict.fromkeys(uslv, np.array((), dtype=int))
        T.dusl = dict.fromkeys(uslv, np.array((), dtype=int))

        tsl = np.array(())
        rsl = np.array(())
        
        pdb.set_trace()
        # loop on group of interactions 
        for k in self:
            if k !=0:
                
                uR = uT = uD = uRf = uRc = 0.

                # nstr : i x r
                nstr = self[k]['sig'][0, 1:-1, :]

                # ityp : i x r
                ityp = self[k]['sig'][1, 1:-1, :]

                # theta : i x r   ( related to interaction )
                theta = self[k]['theta']

                # (i+1) x r
                si = self[k]['si']

                ## flatten information
                ######################

                # reshape nstr in order to be flat (1 dimension)
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
                ray2nbi=np.ones((nbray))

                
                try:
                    self.ray2nbi=np.hstack((self.ray2nbi,ray2nbi))
                except:
                    self.ray2nbi=ray2nbi

                self.ray2nbi[self[k]['rayidx']]  = k
                nbrayt = nbrayt + nbray
                self.nray = self.nray + self[k]['nbrays']
                idxf = idx.reshape(idx.size,order='F')
                #  (i+1)xr
                size2 = si[:, :].size
                #  ,(i+1)xr
                sif = si[:, :].reshape(size2,order='F')
                # 2x2,(i+1)xr

                b0 = self[k]['B'][:,:,0,:]
                b = self[k]['B'][:,:,1:,:].reshape(2, 2, size2-nbray,order='F')
                ## find used slab
                ##################
                # find slab type for the rnstr

                sl = sla[nstrf]

                # seek for interactions position
                ################################
                uR = np.where((itypf == 1))[0]
                uT = np.where((itypf == 2))[0]
                uD = np.where((itypf == 3))[0]
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
                self[k]['rays'] = np.array(([[0]]))
                self[k]['nbrays'] = 1
                self[k]['rayidx'] = ze
                self.nray=1
                self.ray2nbi=ze
                B.stack(data=np.eye(2)[np.newaxis,:,:], idx=ze)
                B0.stack(data=np.eye(2)[np.newaxis,:,:],idx=ze)

        T.create_dusl(tsl)
        R.create_dusl(rsl)
        self.I = I
        self.I.add([T, R])
        self.B = B
        self.B0 = B0


    def eval(self,fGHz=np.array([2.4])):
        """  docstring for eval

        Parameters
        ----------

        fGHz : array
            
        """

        print 'Rays evaluation'

        self.I.eval(fGHz)
        B=self.B.eval(fGHz)
        B0=self.B0.eval(fGHz)




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
        for l in self:
            if l != 0:
                # l stands for the number of interactions
                r = self[l]['nbrays']

                # reshape in order to have a 1D list of index
                # reshape ray index
                rrl = self[l]['rays'].reshape(r*l,order='F')

                # get the corresponding evaluated interactions
                A = self.I.I[:, rrl, :, :].reshape(self.I.nf, r, l, 2, 2,order='F')
                Bl = B[:, rrl, :, :].reshape(self.I.nf, r, l, 2, 2,order='F')
                B0l = B0[:, self[l]['rayidx'], :, :]
                alpha = self.I.alpha[rrl].reshape(r, l,order='F')
                gamma = self.I.gamma[rrl].reshape(r, l,order='F')
                si0 = self.I.si0[rrl].reshape(r, l,order='F')
                sout = self.I.sout[rrl].reshape(r, l,order='F')


                aoa[:,self[l]['rayidx']]=self[l]['aoa']
                aod[:,self[l]['rayidx']]=self[l]['aod']

                try:
                    del Z
                except:
                    pass



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

                    X = A [:, :, i, :, :]
                    Y = Bl[:, :, i, :, :]
                    ## Dot product interaction X Basis
                    Atmp = np.sum(X[..., :, :, np.newaxis]*Y[
                                  ..., np.newaxis, :, :], axis=-2)   #*D[np.newaxis,:,np.newaxis,np.newaxis]
                    if i == 0:
                    ## First Baspdis added
                        A0 = B0l[:, :,  :, :]
                        Z = np.sum(A0[..., :, :, np.newaxis]*Atmp[
                                   ..., np.newaxis, :, :], axis=-2)
                    else:
                        # dot product previous interaction with latest
                        Z = np.sum(Z[..., :, :, np.newaxis]*Atmp[
                                   ..., np.newaxis, :, :], axis=-2)

                # fill the C tilde
                Ct[:, self[l]['rayidx'], :, :] = Z[:, :, :, :]
                # delay computation:
                # sum the distance from antenna to first interaction si0
                # and the sum of all outgoing segments
                #self[l]['dis'] = self.I.si0[self[l]['rays'][0,:]] \
                #        + np.sum(self.I.sout[self[l]['rays']], axis=0)

                # Power losses due to distances
                # will be removed once the divergence factor will be implemented
                Ct[:, self[l]['rayidx'], :, :] = Ct[:, self[l][
                    'rayidx'], :, :]*1./(self[l]['dis'][np.newaxis, :, np.newaxis, np.newaxis])
                self.delays[self[l]['rayidx']] = self[l]['dis']/0.3
                self.dis[self[l]['rayidx']] = self[l]['dis']
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



        Cn = Ctilde()
        Cn.Cpp = bs.FUsignal(self.I.fGHz, c11)
        Cn.Ctp = bs.FUsignal(self.I.fGHz, c12)
        Cn.Cpt = bs.FUsignal(self.I.fGHz, c21)
        Cn.Ctt = bs.FUsignal(self.I.fGHz, c22)
        Cn.nfreq = self.I.nf
        Cn.nray = self.nray
        Cn.tauk = self.delays
        Cn.fGHz = self.I.fGHz
        # r x 2
        Cn.tang = aod.T
        # r x 2
        Cn.rang = aoa.T
        # add aoa and aod 
        return(Cn)


    def ray(self, r):
        """

        Parameters
        ----------
        r : integer
            ray index
        
        Notes
        -----

            Give the ray number and it returns the index of its interactions
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

    def info(self, r):
        '''
            provides information for a given ray r

        Parameters
        ----------
        r : int
            ray index

        '''
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
        print ' Matrix of ray #', r, 'at f=', self.I.fGHz[0]
        print '----------------------------------------'

        print 'rotation matrix#', 'type: B0'
        print self.B0.data[r,:,:]
        for iidx, i in enumerate(typ):
            print 'interaction #', ray[iidx], 'type:', i
            # f x l x 2 x 2
            print self.I.I[0, ray[iidx], :, :]
            print 'rotation matrix#',[ray[iidx]], 'type: B'
            print self.B.data[ray[iidx], :, :]


    def signature(self, L):
        """
        """
        sig = Signatures(L, self.pTx, self.pRx)
        for k in self:
            sig[k] = self[k]['sig']
        return(sig)

    def show3d(self,
               ray,
               bdis=True,
               bbas=False,
               bstruc=True,
               col=np.array([1, 0, 1]),
               id=0,
               linewidth=1):
        """ plot a 3D ray

        Parameters
        ----------

        bdis :
            display boolean - if False return .vect filename
        bbas :
            display local basis
        bstruc :
            display structure
        col  :
            color of the ray
        id   :
            id of the ray
        linewidth :

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

    def show3(self, bdis=True, bstruc=True, id=0, strucname='defstr', ilist=[], raylist=[],pg=np.array([[0],[0],[0]])):
        """ plot 3D rays within the simulated environment

        Parameters
        ----------
            bdis : boolean
                True
            bstruc : boolean
                True
            id : int 
            strucname : string 
                'defstr'
            ilist : list of group of interactions
            raylist : list of index rays 
            pg : centroid of the structure 
         

        Examples
        -------- 

        """
        if ilist == []:
            ilist = self.keys()
        pTx = self.pTx.reshape((3, 1))-pg
        pRx = self.pRx.reshape((3, 1))-pg
        filename = pyu.getlong("grRay" + str(id) + ".list", pstruc['DIRGEOM'])
        fo = open(filename, "w")
        fo.write("LIST\n")
        if bstruc:
            fo.write("{<"+strucname+".off}\n")
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
                    col = np.array([2, 0, 1])
                    # print ray
                    fileray = self.show3d(ray=ray, bdis=False,
                                          bstruc=False, col=col, id=k)
                    k += 1
                    fo.write("{< " + fileray + " }\n")
        fo.close()
        if (bdis):
            chaine = "geomview " + filename + " 2>/dev/null &"
            os.system(chaine)
        else:
            return(filename)


if __name__ == "__main__":
    doctest.testmod()
