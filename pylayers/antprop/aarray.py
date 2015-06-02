# -*- coding:Utf-8 -*-
import numpy as np
import pylayers.antprop.antenna as ant
import pdb
r"""

.. currentmodule:: pylayers.antprop.array

This module handles antenna arrays

"""
class AntSet(object):
    """ Antenna Set
    """

    def __init__(self):
        defaults = {'Na'   : 8,
                    'A'    : []
                   }

class AntArray(AntSet,Antenna):
    """
    """

    def __init__(self,**kwargs):
        defaults = {'typ': 'UA',
                    'N'   : [8,1,1],
                    'dm'   :[0.075,0,0],
                    'A'    : []
                   }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.typ = kwargs['typ']
        self.N   = np.array(kwargs['N'])
        self.dm  = np.array(kwargs['dm'])

        if kwargs['A']<>[]:
            self.A   = kwargs['A']
        else:
            self.A   = ant.Antenna('Omni')


    def __repr__(self):
         st = "Antenna Array : \n"
         st = st + 'typ : '+self.typ+'\n'
         st = st + 'N : '+str(self.N)+'\n'
         st = st + 'dm : '+str(self.dm)+'\n'
         return(st)

    def F(self,ang,fGHz=2.0):
        """
         Parameters
         ----------

         ang : np.array(Nkx2) [theta,phi]
            array direction angles in radians

         Examples
         --------

         >>> Nk = 180
         >>> theta = np.pi*np.ones(Nk)/2.
         >>> phi = np.linspace(0,np.pi,Nk)
         >>> ang = np.vstack((theta,phi)).T
         >>> A = AntArray()
         >>> A.F(ang)

        """

        lamda = 0.3/fGHz
        k     = 2*np.pi/lamda

        kx = np.sin(ang[:,0])*np.cos(ang[:,1])    # N x 1
        ky = np.sin(ang[:,0])*np.sin(ang[:,1])    # N x 1
        kz = np.cos(ang[:,0])                     # N x 1
        self.k  = np.vstack((kx,ky,kz)).T         # N x 3

        Na = np.prod(self.N)  # number of antennas
        Nx = self.N[0]
        Ny = self.N[1]
        Nz = self.N[2]
        if Nx%2==0:
            Dx = self.dm[0]*np.linspace(-Nx/2,Nx/2,Nx)[None,:,None,None] # 1 x Nx x Ny x Nz
        else:
            Dx = self.dm[0]*np.linspace(-(Nx-1)/2,(Nx-1)/2,Nx)[None,:,None,None] #Nx x Ny x Nz
        if Ny%2==0:
            Dy = self.dm[1]*np.linspace(-Ny/2,Ny/2,Ny)[None,None,:,None] # Nx x Ny x Nz
        else:
            Dy = self.dm[1]*np.linspace(-(Ny-1)/2,(Ny-1)/2,Ny)[None,None,:,None] # Nx x Ny x Nz
        if Nz%2==0:
            Dz = self.dm[2]*np.linspace(-Nz/2,Nz/2,Nz)[None,None,None,:] #  Nx x Ny x Nz
        else:
            Dz = self.dm[2]*np.linspace(-(Nz-1)/2,(Nz-1)/2,Nz)[None,None,None,:] # Nx x Ny x Nz

        self.D = np.zeros((3,Nx,Ny,Nz))
        self.D[0,:,:,:] = Dx
        self.D[1,:,:,:] = Dy
        self.D[2,:,:,:] = Dz

        #
        # F = exp(+jkD)
        #
        kD = np.einsum('ij,jklm->iklm',self.k,self.D)   # k.D
        self.E = np.exp(1j*k*np.einsum('ij,jklm->iklm',self.k,self.D))
        self.F = np.einsum('ijkl->i',self.E)
